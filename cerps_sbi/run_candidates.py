"""
Run the C1-C4 candidate sweeps against the REAL CERPS kernel.

Each candidate is mapped to REAL CERPS parameter axes (no fabricated API).
Verification config (fast, defensible):
    T=1500, n_init=200, seeds=(3, 11), snapshot at 0.95*T.
Production config (use for the paper):
    T=4000, n_init=600, seeds=(3, 11, 21), see analysis/src/mechanism_scan.py.

Honesty flags:
  * C1's P8 (delay) axis is NOT in the code -> only the P1 (niche) axis is run;
    the delay component is reported as "unimplemented".
  * All numbers below are produced by run_cerps(); nothing is hand-written.
"""
import sys, os, time, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cerps_sbi import legacy_config, run_one, NG  # noqa

# ---- verification vs production config -------------------------------------
VERIFY = dict(T=1500.0, n_init=200, seeds=(3, 11))
PROD   = dict(T=4000.0, n_init=600, seeds=(3, 11, 21))

# Selected with the environment variable CERPS_SWEEP_MODE=prod|verify
# (default: verify, i.e. the committed v1 result).  "prod" writes
# results/candidate_sweep_v2.json and is the paper-scale run.
MODE = os.environ.get("CERPS_SWEEP_MODE", "verify").strip().lower()
CUR = PROD if MODE == "prod" else VERIFY
OUT_NAME = "candidate_sweep_v2.json" if MODE == "prod" else "candidate_sweep_v1.json"

SNAP_AT = 0.95  # fraction of T for the attractor snapshot


def run_axis(name, axis_label, levels, seed_list=None, snapshot=False):
    """levels: list of (label, override_dict). Returns list of run dicts."""
    seeds = seed_list or CUR["seeds"]
    snap_t = (CUR["T"] * SNAP_AT,) if snapshot else ()
    rows = []
    for lvl_label, override in levels:
        for sd in seeds:
            t0 = time.time()
            m = run_one(override=override, seed=sd, T=CUR["T"],
                        n_init=CUR["n_init"], snapshot_times=snap_t)
            m["candidate"] = name
            m["axis"] = axis_label
            m["level"] = lvl_label
            m["override"] = override
            m["runtime_s"] = round(time.time() - t0, 1)
            rows.append(m)
            print(f"  [{name}/{lvl_label}/s{sd}] Nfin={m['N_final']:.0f} "
                  f"div={m['div']} termD={m['term_death']} cyc={m['cyc_frac']:.2f} "
                  f"mG={m['mG']:.2f} mP={m['mP']:.2f} ({m['runtime_s']:.0f}s)")
    return rows


def main():
    t0 = time.time()
    out = {"config": dict(CUR, snapshot_at=SNAP_AT), "candidates": {}}

    print("== baseline (legacy anchor) ==")
    base = run_axis("baseline", "none", [("default", {})])
    out["candidates"]["baseline"] = {"axis": "none", "runs": base}

    print("== C1  P1 (niche) x [P8 delay: UNIMPLEMENTED] ==")
    # P1 axis: niche epigenetic reference strength h_ref_S (default 0.40)
    c1 = run_axis("C1", "h_ref_S (niche->epi reference)",
                  [("off", {"h_ref_S": 0.0}),
                   ("default", {"h_ref_S": 0.40}),
                   ("strong", {"h_ref_S": 1.20})])
    out["candidates"]["C1"] = {
        "axis": "P1 only (niche). P8 delay NOT in code -> unverifiable here",
        "runs": c1}

    print("== C2  P5 (epigenetics) x P6 (memory/age) ==")
    # P6 axis: kappa_a (epigenetic-age gating of death; default 0.03)
    c2 = run_axis("C2", "kappa_a (memory-of-age weight)",
                  [("off", {"kappa_a": 0.0}),
                   ("default", {"kappa_a": 0.03}),
                   ("strong", {"kappa_a": 0.15})])
    out["candidates"]["C2"] = {"axis": "P6 via kappa_a", "runs": c2}

    print("== C3  P3 (partitioning) x P5 (epigenetics) ==")
    # P3 axis: s_part (expression partitioning noise at mitosis; default 0.06)
    c3 = run_axis("C3", "s_part (partitioning noise)",
                  [("off", {"s_part": 0.0}),
                   ("default", {"s_part": 0.06}),
                   ("strong", {"s_part": 0.15})])
    out["candidates"]["C3"] = {"axis": "P3 via s_part", "runs": c3}

    print("== C4  GRN mutual-repression (control / negative finding) ==")
    # vary cross-repression bG,bP on the G<->P toggle (default 1.0,1.0)
    c4 = run_axis("C4", "b tuple (toggle cross-repression)",
                  [("weak", {"b": (0.3, 0.3, 1.2, 1.2)}),
                   ("default", {"b": (1.0, 1.0, 1.2, 1.2)}),
                   ("strong", {"b": (3.0, 3.0, 1.2, 1.2)})],
                  snapshot=True)
    out["candidates"]["C4"] = {"axis": "GRN cross-repression", "runs": c4}

    # ---- aggregate means per level -----------------------------------------
    agg = {}
    for cname, cblock in out["candidates"].items():
        levels = {}
        for r in cblock["runs"]:
            lv = r["level"]
            d = levels.setdefault(lv, [])
            d.append(r)
        rows = []
        for lv, rs in levels.items():
            rows.append(dict(
                level=lv,
                N_final=round(float(np_mean(rs, "N_final")), 1),
                div=round(float(np_mean(rs, "div")), 0),
                term_death=round(float(np_mean(rs, "term_death")), 0),
                cyc_frac=round(float(np_mean(rs, "cyc_frac")), 3),
                mG=round(float(np_mean(rs, "mG")), 3),
                mP=round(float(np_mean(rs, "mP")), 3),
                sG=round(float(np_mean(rs, "sG")), 3),
                sP=round(float(np_mean(rs, "sP")), 3),
                frac_both_on=round(float(np_mean(rs, "frac_both_on", default=-1)), 3),
                frac_G_only=round(float(np_mean(rs, "frac_G_only", default=-1)), 3),
                frac_P_only=round(float(np_mean(rs, "frac_P_only", default=-1)), 3),
            ))
        agg[cname] = rows

    out["aggregate"] = agg
    out["elapsed_s"] = round(time.time() - t0, 1)

    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "results", OUT_NAME)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nsaved {path} ({out['elapsed_s']}s)")
    print("\n===== AGGREGATE (mean over seeds) =====")
    for cname, rows in agg.items():
        print(f"\n{cname}: {out['candidates'][cname]['axis']}")
        for r in rows:
            print("  ", r)


def np_mean(rs, key, default=float("nan")):
    vals = [r[key] for r in rs if key in r]
    return float(np.mean(vals)) if vals else default


if __name__ == "__main__":
    main()
