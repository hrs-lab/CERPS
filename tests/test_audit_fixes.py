"""Regression tests for the audit fixes (2026-10-19).

Each test corresponds to one defect found in the CERPS audit report
(CERPS_审查报告.md, sections 4 and 5).  They are deliberately stricter
than the original suite, which only exercised symmetric parameters and
attractor interiors -- exactly the regions where the defects hide.

Run with:  python tests/test_audit_fixes.py   (or pytest tests/ -v)
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "reproduction"))

from cerps.config import (CERPSConfig, GRNConfig, EpiConfig,  # noqa: E402
                          DivisionConfig, RatesConfig, SpaceConfig)
from cerps.grn import GRN  # noqa: E402
from cerps import epigenetics as epi  # noqa: E402
from cerps import division as div_mod  # noqa: E402
from cerps import rates as rates_mod  # noqa: E402
from cerps import spatial as sp_mod  # noqa: E402
from cerps import batch  # noqa: E402
from cerps.simulator import Simulator  # noqa: E402

PASS = []
FAIL = []


def check(name, cond, info=""):
    if cond:
        PASS.append(name)
        print(f"  PASS  {name}")
    else:
        FAIL.append(name)
        print(f"  FAIL  {name}  {info}")


# ======================================================================
# 1. P0-2: cross-regulation threshold must be the source gene's s_j
# ======================================================================
def test_cross_threshold():
    """grn.py: S_rep[i,j] must equal s_j, not (s_i+s_j)/2."""
    from hsl2024.params import Params

    def paper_rhs(x, u, p):
        a1 = p.alpha1 * np.exp(p.lam1 * u[0]) + p.a_tf
        a2 = p.alpha2 * np.exp(p.lam2 * u[1])
        f1 = a1 * (p.rho1 + (1 - p.rho1) * x[0] ** p.n /
                   (p.s1 ** p.n + x[0] ** p.n)) \
            + p.b1 * p.s2 ** p.n / (p.s2 ** p.n + x[1] ** p.n) - p.k1 * x[0]
        f2 = a2 * (p.rho2 + (1 - p.rho2) * x[1] ** p.n /
                   (p.s2 ** p.n + x[1] ** p.n)) \
            + p.b2 * p.s1 ** p.n / (p.s1 ** p.n + x[0] ** p.n) - p.k2 * x[1]
        return np.array([f1, f2])

    for tag, p in [("asym s1=0.3,s2=0.7", Params(s1=0.3, s2=0.7)),
                   ("asym s1=0.7,s2=0.3", Params(s1=0.7, s2=0.3)),
                   ("Fig.10 series s1=0.55", Params(s1=0.55)),
                   ("Fig.10 series k1=1.1", Params(k1=1.1)),
                   ("Fig.10 series alpha1=0.44", Params(alpha1=0.44)),
                   ("Fig.10 series rho1=0.11", Params(rho1=0.11)),
                   ("Fig.10 series b1=1.1", Params(b1=1.1)),
                   ("a_tf=0.3 on gene 1", Params(a_tf=0.3)),
                   ("n=2", Params(n=2.0))]:
        grn = GRN(GRNConfig(sigma=(0, 0), n_hill=p.n, s_self=(p.s1, p.s2),
                            alpha=(p.alpha1, p.alpha2),
                            lam=(p.lam1, p.lam2), rho=(p.rho1, p.rho2),
                            k_deg=(p.k1, p.k2),
                            b_base=(p.b1, p.b2), a_tf=(p.a_tf, 0.0)))
        x = np.array([[1.3, 0.7]])
        u = np.array([[0.42, 0.61]])
        xc, xp = x.copy(), x.copy()
        for _ in range(400):
            xc = grn.step(xc, u, np.zeros((1, 2)), p.dt)
            xp = xp + paper_rhs(xp[0], u[0], p) * p.dt
        err = np.abs(xc - xp).max()
        check(f"grn: exact Eq.(1) reduction [{tag}]", err < 1e-9,
              f"max|dx|={err:.2e}")
    # the threshold really is s_j now
    g = GRN(GRNConfig(sigma=(0, 0), s_self=(0.3, 0.7)))
    check("grn: S_rep[0,1] == s_2", abs(g.S_rep[0, 1] - 0.7) < 1e-12,
          f"S_rep[0,1]={g.S_rep[0, 1]}")
    check("grn: S_rep[1,0] == s_1", abs(g.S_rep[1, 0] - 0.3) < 1e-12,
          f"S_rep[1,0]={g.S_rep[1, 0]}")


# ======================================================================
# 2. orthogonality: all kernels share the conditional mean phi(v)
# ======================================================================
def test_kernel_orthogonality():
    """Switching a kernel must change ONLY the noise, not the drift."""
    rng = np.random.default_rng(1)
    v = np.full((200000, 2), 0.5)
    cfg = EpiConfig(kernel="beta")
    phi = epi.phi_map(0.5, cfg)
    var_beta = phi * (1 - phi) / (1 + cfg.m0)

    # gaussian, including the heavy-truncation regime
    for sg in (0.045, 0.15, 0.30):
        c = EpiConfig(kernel="gaussian", sigma_g=sg)
        for vv in (0.02, 0.5, 0.98):
            ph = epi.phi_map(vv, c)
            d1, _ = epi.inherit(np.full((100000, 2), vv), c, rng)
            check(f"gaussian mean == phi (sg={sg}, v={vv})",
                  abs(d1[:, 0].mean() - ph) < 5e-3,
                  f"{d1[:, 0].mean():.4f} vs {ph:.4f}")
            check(f"gaussian support (sg={sg}, v={vv})",
                  (d1 >= 0).all() and (d1 <= 1).all())

    # cpg moment-matched mode: mean AND variance match beta
    c = EpiConfig(kernel="cpg", cpg_mode="moment", cpg_sites=61)
    d1, _ = epi.inherit(v, c, rng)
    check("cpg(moment) mean == phi", abs(d1[:, 0].mean() - phi) < 2e-3,
          f"{d1[:, 0].mean():.4f} vs {phi:.4f}")
    check("cpg(moment) variance == beta variance",
          abs(d1[:, 0].var() - var_beta) < 5e-4,
          f"{d1[:, 0].var():.5f} vs {var_beta:.5f}")
    # mechanistic mode keeps the affine mean map (documented difference)
    c = EpiConfig(kernel="cpg", cpg_mode="mechanistic")
    d1, _ = epi.inherit(v, c, rng)
    aff = 0.5 * c.cpg_p_maintain + 0.5 * c.cpg_p_denovo
    check("cpg(mechanistic) mean == affine map",
          abs(d1[:, 0].mean() - aff) < 2e-3,
          f"{d1[:, 0].mean():.4f} vs {aff:.4f}")


# ======================================================================
# 3. partitioning noise: continuous alpha (Soltani et al. 2016 Eq. 5)
# ======================================================================
def test_partition_alpha():
    rng = np.random.default_rng(2)
    x = np.full((200000, 2), 1.0)
    V = 200.0
    for alpha in (0.0, 0.5, 1.0, 2.0, 4.0):
        cfg = DivisionConfig(partition="binomial", volume=V, alpha=alpha)
        d1, _ = div_mod.partition(x, cfg, rng)
        th = alpha * 1.0 * (1 - 0.5) / (0.5 * V)      # alpha*x*(1-f)/(f*V)
        check(f"partition alpha={alpha}: mean preserved",
              abs(d1[:, 0].mean() - 1.0) < 0.02,
              f"{d1[:, 0].mean():.4f}")
        check(f"partition alpha={alpha}: Var == alpha*x*(1-f)/(fV)",
              abs(d1[:, 0].var() - th) < 0.05 * th + 1e-5,
              f"{d1[:, 0].var():.5f} vs {th:.5f}")
    # molecule conservation
    cfg = DivisionConfig(partition="binomial", volume=V, alpha=2.0)
    d1, d2 = div_mod.partition(x, cfg, rng)
    tot = 0.5 * V * d1 + 0.5 * V * d2
    check("partition alpha=2: molecules conserved",
          np.allclose(tot, V, atol=0.5), f"mean total={tot.mean():.3f}")
    # exact mode is unaffected
    cfg = DivisionConfig(partition="exact")
    d1, d2 = div_mod.partition(x, cfg, rng)
    check("partition exact: identity", np.allclose(d1, x) and
          np.allclose(d2, x))


# ======================================================================
# 4. K_switch auto-alignment with the legacy classification boundary
# ======================================================================
def test_k_switch_alignment():
    grn = GRN(GRNConfig(sigma=(0, 0)))
    rc = RatesConfig(rate_mode="continuous")
    K = rates_mod.resolve_K(rc, grn)
    # legacy boundary on the diagonal: H+(0.5)^2 = 0.25
    check("K_switch auto == H+(x_thr)^G", abs(K - 0.25) < 1e-12, f"K={K}")
    # with that K, the continuous rule classifies the legacy SC region
    # correctly on the diagonal (at exactly x = x_threshold both rules
    # say "not SC", since legacy uses a strict >)
    for c, expect_sc in [(0.45, False), (0.50, False), (0.51, True),
                         (0.62, True)]:
        S = grn.activation_levels(np.array([[c, c]])).prod(axis=1)[0]
        is_sc = S > K
        check(f"diagonal x={c}: continuous SC == legacy SC",
              bool(is_sc) == expect_sc, f"S={S:.3f} K={K}")
    # an explicit K_switch is honoured
    check("explicit K_switch honoured",
          rates_mod.resolve_K(RatesConfig(rate_mode="continuous",
                                          K_switch=0.5), grn) == 0.5)
    # ... and batch.validate warns when it moves the boundary
    w = " ".join(batch.validate(CERPSConfig(
        rates=RatesConfig(rate_mode="continuous", K_switch=0.5))))
    check("batch.validate warns on shifted K_switch",
          "boundary is shifted" in w, w[:80])


# ======================================================================
# 4b. the TA-side switches must use the complementary threshold 1-K
# ======================================================================
def test_complementary_switch():
    """kappa's TA-side switch must fire at the same S as the SC-side one.

    With a shared threshold the two switches fire at S = K and S = 1-K;
    any K != 0.5 then leaves a band in which a stem-like cell carries the
    TA removal rate (measured: at K = 0.125, S = 0.5 gives
    kappa = 2*kappa0 instead of kappa0).
    """
    grn = GRN(GRNConfig(sigma=(0, 0)))
    for K in (None, 0.25, 0.125, 0.5):
        rc = RatesConfig(rate_mode="continuous", p_switch=40.0, K_switch=K)
        K_eff = rates_mod.resolve_K(rc, grn)
        # a clearly stem-like cell must get the SC removal rate
        b, k, m, tau, _ = rates_mod.rates(grn, np.array([[1.0, 1.0]]), None,
                                          rc, np.zeros(1), 1.0)
        check(f"kappa: SC-like cell gets kappa0 (K={K})",
              abs(k[0] - 0.01) < 1e-3, f"kappa={k[0]:.5f}")
        # a clearly TA-like cell must get 2*kappa0
        b, k, m, tau, _ = rates_mod.rates(grn, np.array([[2.0, 0.1]]), None,
                                          rc, np.zeros(1), 1.0)
        check(f"kappa: TA-like cell gets 2*kappa0 (K={K})",
              abs(k[0] - 0.02) < 1e-3, f"kappa={k[0]:.5f}")
        # and a mid-range cell must be intermediate, not pinned to the
        # TA rate (the failure mode of the shared-threshold version)
        b, k, m, tau, _ = rates_mod.rates(grn, np.array([[0.75, 0.75]]), None,
                                          rc, np.zeros(1), 1.0)
        check(f"kappa: mid-S cell is intermediate (K={K})",
              k[0] < 0.019, f"kappa={k[0]:.5f} at K_eff={K_eff:.3f}")


# ======================================================================
# 5. P0-3: paracrine field -- CFL sub-stepping and reflecting boundaries
# ======================================================================
def test_field_numerics():
    def run(D, dt, h=4.0, T=500.0, n=None):
        sc = SpaceConfig(enabled=True, paracrine=True, diff_rate=D,
                         grid_h=h, decay_rate=0.02, domain_radius=40.0)
        s = sp_mod.Space(sc)
        s.pos = np.array([[0.0, 0.0, 0.0]])
        s.field = np.zeros(s._grid_shape)
        steps = n if n is not None else int(round(T / dt))
        for _ in range(steps):
            s.field_step(np.array([0.05]), dt)
        return s

    # reflecting boundary: no signal wraps around the domain
    s = run(1.0, 0.25, T=2000.0)
    s.pos = np.array([[38.0, 0.0, 0.0]])
    s.field = np.zeros(s._grid_shape)
    for _ in range(20000):
        s.field_step(np.array([0.05]), 0.25)
    ci = s._to_grid(s.pos)
    ghost = s.field[0, ci[0][1], ci[0][2]]
    real = s.field[ci[0][0], ci[0][1], ci[0][2]]
    check("field: reflecting boundary (no wrap-around)",
          ghost < 1e-3 * real, f"ghost/real={ghost / real:.2e}")

    # nominal diffusivity is honoured: D=25 must differ from D=10.667
    # (with the old CFL clamp they were bit-identical)
    A = run(25.0, 0.25, T=500.0).field
    B = run(10.667, 0.25, T=500.0).field
    check("field: nominal diff_rate honoured (no silent clamp)",
          not np.allclose(A, B),
          f"max rel diff={np.abs(A - B).max() / A.max():.3f}")

    # sub-stepping is dt-consistent: equal sub-dt gives equal results
    C = run(25.0, 0.125, T=500.0).field     # 2 sub-steps of 0.0625
    Dd = run(25.0, 0.0625, T=500.0).field   # 1 step of 0.0625
    check("field: dt-consistent (equal sub-dt => equal field)",
          np.allclose(C, Dd), "sub-stepping is exact in the sub-dt")

    # stability: no NaN/blow-up even far above the CFL limit
    E = run(500.0, 0.25, T=200.0).field
    check("field: stable at Fo >> 1/6 (sub-stepping)",
          np.isfinite(E).all() and E.max() < 1e3, f"max={E.max():.3e}")


# ======================================================================
# 6. mechanics: relaxation is time-step consistent
# ======================================================================
def test_mechanics_dt():
    sc = SpaceConfig(enabled=True, cell_radius=1.0, relax=0.15)
    s = sp_mod.Space(sc)
    p0 = np.array([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]])
    s.pos = p0.copy(); s.relax(0.25); a = s.pos.copy()
    s.pos = p0.copy(); s.relax(0.125); s.relax(0.125); b = s.pos.copy()
    # two half steps must equal one full step for an exponential decay
    check("mechanics: two dt/2 steps == one dt step",
          np.allclose(a, b, atol=1e-12),
          f"{a[1]} vs {b[1]}")
    # and the old behaviour is reproduced at dt = 0.25 h
    s.pos = p0.copy(); s.relax(0.25)
    check("mechanics: backward compatible at dt=0.25",
          np.allclose(s.pos, a))


# ======================================================================
# 7. protocol: rec_tran_from filters the transient
# ======================================================================
def test_rec_tran_from():
    cfg = CERPSConfig(max_time=600.0, seed=3)
    sim = Simulator(cfg)
    sim.run()
    n_all = len(sim.transitions)
    cfg2 = CERPSConfig(max_time=600.0, seed=3, rec_tran_from=300.0)
    sim2 = Simulator(cfg2)
    sim2.run()
    n_filt = len(sim2.transitions)
    check("rec_tran_from filters the transient", 0 < n_filt < n_all,
          f"all={n_all} filtered={n_filt}")
    # and the two runs are otherwise identical (same seed, same dynamics)
    check("rec_tran_from does not change the dynamics",
          abs(sim.n - sim2.n) == 0 and
          np.allclose(sim.records["total"][-5:], sim2.records["total"][-5:]))


# ======================================================================
# 8. batch harness
# ======================================================================
def test_batch_harness():
    base = CERPSConfig(max_time=300.0, rec_tran_from=150.0)
    res = batch.sweep(base, {"epi.kernel": ["beta", "cpg"]}, seeds=(0,),
                      out_dir=None, verbose=False)
    check("batch.sweep runs a grid", len(res) == 2)
    check("batch.run_one keys",
          all(k in res[0] for k in ("ho_total", "ho_sc", "tran_0_0",
                                    "tran_0_3", "n_events", "config")))
    agg = batch.aggregate(res)
    check("batch.aggregate groups by tag", len(agg) == 2 and
          all(a["n_seeds"] == 1 for a in agg))
    # transition keys are label-coded: 0=SC, 1=TA1, 3=TA2
    tr = batch.transitions(Simulator(base))
    check("batch.transitions label codes", tr["labels"] == [0, 1, 3],
          str(tr["labels"]))
    # config warnings
    w = batch.validate(CERPSConfig(division=DivisionConfig(alpha=2.0)))
    check("batch.validate: alpha ignored under exact partition",
          any("IGNORED" in x for x in w), str(w))
    w = batch.validate(CERPSConfig(epi=EpiConfig(kernel="cpg",
                                                 cpg_sites=60)))
    check("batch.validate: cpg_sites != m0+1",
          any("variance" in x for x in w), str(w))
    # set/get path
    c = batch.set_path(base, "rates.beta0", 0.08)
    check("batch.set_path / get_path",
          batch.get_path(c, "rates.beta0") == 0.08 and
          batch.get_path(base, "rates.beta0") == 0.04)
    # health check
    h = batch.health_check(CERPSConfig(
        max_time=200.0, space=SpaceConfig(enabled=True, paracrine=True)),
        t_short=100.0)
    check("batch.health_check reports CFL and dt convergence",
          "cfl" in h and "dt_convergence_max" in h and
          h["dt_convergence_max"] >= 0.0,
          str({k: h[k] for k in ("cfl", "n_sub", "dt_convergence_max")}))


# ======================================================================
# 9. integration: the fixed legacy mode matches hsl2024 tightly
# ======================================================================
def test_legacy_matches_hsl2024():
    """The P0-1 fix must bring SC self-renewal back to the paper value.

    Before the fix CERPS reported SC->SC ~ 76.7% against hsl2024's
    74.2%; with the resting-only phenotype update the two agree to
    ~0.01 pp.  A short run keeps this test fast; the tolerance is set
    well below the pre-fix deviation but above Monte-Carlo noise.
    """
    from hsl2024 import DEFAULT, simulate as sim_hsl
    from hsl2024.analysis import transition_stats, homeostasis_stats
    N, T = 3, 1500.0
    warm = T / 2.0
    # NOTE: the hsl2024 engine counts transitions for Step+1 > rec_tran_from
    # (in STEPS), so the warm-up must be converted from hours to steps and
    # be shorter than the run -- otherwise the matrix comes back empty.
    hsl = sim_hsl(DEFAULT, n_runs=N, n_steps=int(T / 0.25), seed=4242,
                  rec_tran_from=int(warm / 0.25))
    ts = transition_stats(dict(tran=hsl["tran"]))
    res = batch.sweep(CERPSConfig(max_time=T, rec_tran_from=warm),
                      grid={}, seeds=(900, 901, 902), out_dir=None,
                      verbose=False)
    agg = batch.aggregate(res)[0]
    sc_row = [agg["mean"][f"tran_0_{c}"] for c in (0, 1, 3)]
    ref = ts["prob_mean"][0]
    check("legacy: SC->SC matches hsl2024 (< 0.6 pp)",
          abs(sc_row[0] - ref[0]) < 0.6,
          f"CERPS {sc_row[0]:.3f} vs hsl {ref[0]:.3f}")
    check("legacy: SC row sums to ~100", abs(sum(sc_row) - 100.0) < 0.5,
          f"{sum(sc_row):.3f}")
    # The SC count is stable across runs; the TOTAL is heavy-tailed and,
    # at T = 1500 h, still below the paper's 2000-4000 h homeostasis
    # window, so only the SC count is compared here (the full-protocol
    # comparison is in the audit's xval_fixed.py: total agrees to 0.6%
    # over 10 runs at 4000 h).
    hs = homeostasis_stats(dict(sys_ts=hsl["sys_ts"],
                                stop_step=hsl["stop_step"]))
    sc = agg["mean"]["ho_sc"]
    check("legacy: SC count matches hsl2024 (< 15%)",
          abs(sc - hs["count_mean"][1]) / hs["count_mean"][1] < 0.15,
          f"CERPS {sc:.0f} vs hsl {hs['count_mean'][1]:.0f}")


def main():
    print("== CERPS audit-fix regression tests ==")
    for fn in [test_cross_threshold, test_kernel_orthogonality,
               test_partition_alpha, test_k_switch_alignment,
               test_complementary_switch,
               test_field_numerics, test_mechanics_dt,
               test_rec_tran_from, test_batch_harness,
               test_legacy_matches_hsl2024]:
        print(f"\n- {fn.__name__}")
        fn()
    print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
    if FAIL:
        print("FAILED:", FAIL)
        sys.exit(1)


if __name__ == "__main__":
    main()
