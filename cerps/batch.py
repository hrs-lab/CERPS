"""Batch experiment harness for the CERPS framework.

The ``cerps/`` package simulates ONE population for ONE configuration.
Mechanism studies need more than that: parameter sweeps, several random
seeds per configuration, checkpoint/resume for long runs, and aggregation
of the statistics that actually carry the conclusions.  This module
supplies that layer.  The interface follows the ``job_*`` pattern of
``reproduction/scripts/common.py``.

Typical use
-----------
>>> from cerps.config import CERPSConfig, RatesConfig
>>> from cerps import batch
>>> base = CERPSConfig(max_time=1000.0)
>>> grid = {"rates.beta0": [0.02, 0.04], "epi.kernel": ["beta", "cpg"]}
>>> res = batch.sweep(base, grid, seeds=(0, 1, 2), out_dir="runs")
>>> agg = batch.aggregate(res)          # mean/std over seeds per config
>>> batch.report(agg)                   # printed table

Every configuration/seed pair is checkpointed to ``out_dir/results.jsonl``
as soon as it finishes, so an interrupted sweep resumes where it stopped.

Statistics
----------
``homeostasis``   time-averaged population counts over a window, exactly as
                  ``hsl2024.analysis.homeostasis_stats`` (per-run mean over
                  the window, then mean/std across runs).
``transitions``   row-normalised mother->daughter transition probabilities,
                  honouring ``CERPSConfig.rec_tran_from``.

Numerical health
----------------
``health_check`` reports the paracrine CFL number (and the number of
sub-steps the solver will take) and a dt-halving convergence estimate, so
that a "mechanism effect" can be told apart from a discretisation effect.
"""
from __future__ import annotations

import itertools
import json
import os
import time
from dataclasses import fields, replace

import numpy as np

from .config import CERPSConfig
from .simulator import Simulator

__all__ = ["run_one", "sweep", "aggregate", "homeostasis", "transitions",
           "health_check", "report", "set_path", "get_path"]


# ----------------------------------------------------------------------
# config path helpers:  "rates.beta0" <-> nested dataclass attribute
# ----------------------------------------------------------------------
def _split(path):
    parts = path.split(".")
    if len(parts) < 2:
        raise ValueError(f"config path {path!r} must be 'section.field'")
    return parts[0], ".".join(parts[1:])


def get_path(cfg: CERPSConfig, path: str):
    sec, fld = _split(path)
    return getattr(getattr(cfg, sec), fld)


def set_path(cfg: CERPSConfig, path: str, value) -> CERPSConfig:
    """Return a copy of ``cfg`` with ``path`` (e.g. 'rates.beta0') set."""
    sec, fld = _split(path)
    sub = getattr(cfg, sec)
    if not hasattr(sub, fld):
        raise ValueError(f"{path!r}: {type(sub).__name__} has no field {fld!r}")
    return replace(cfg, **{sec: replace(sub, **{fld: value})})


def _sections(cfg: CERPSConfig):
    return [f.name for f in fields(cfg)
            if hasattr(getattr(cfg, f.name), "__dataclass_fields__")]


# ----------------------------------------------------------------------
# statistics
# ----------------------------------------------------------------------
def homeostasis(records, t_from=None, t_to=None):
    """Time-averaged counts over [t_from, t_to] of one run.

    Mirrors ``hsl2024.analysis.homeostasis_stats`` (per-run time average
    over the window) but works on the CERPS record format.  Returns
    ``dict(total, sc, ta1, ta2, other, sc_frac)``.
    """
    t = np.asarray(records["t"], float)
    m = np.ones(t.shape, bool)
    if t_from is not None:
        m &= t >= t_from
    if t_to is not None:
        m &= t <= t_to
    if not m.any():
        m = np.ones(t.shape, bool)
    tot = float(np.mean(np.asarray(records["total"])[m]))
    sc = float(np.mean(np.asarray(records["sc"])[m]))
    ta1 = float(np.mean(np.asarray(records["ta1"])[m]))
    ta2 = float(np.mean(np.asarray(records["ta2"])[m]))
    oth = float(np.mean(np.asarray(records["other"])[m]))
    return dict(total=tot, sc=sc, ta1=ta1, ta2=ta2, other=oth,
                sc_frac=(sc / tot if tot > 0 else 0.0))


def transitions(sim):
    """Row-normalised transition probabilities (%) and raw counts.

    Only transitions recorded after ``cfg.rec_tran_from`` are present in
    ``sim.transitions`` (the simulator filters them), so the numbers are
    directly comparable with the paper's Table 3 when
    ``rec_tran_from = 2000``.

    ``prob`` is indexed by the CERPS *label codes*: row/col 0 = SC,
    1 = TA1, 3 = TA2 (and j+1 for single-gene-high states when G > 2).
    Hence ``prob[0, 3]`` is SC -> TA2, matching the paper's Table 3.
    """
    G = sim.grn.G
    labels = [0, 1, 3] if G == 2 else [0] + list(range(1, G + 1))
    idx = {l: i for i, l in enumerate(labels)}
    B = np.zeros((len(labels), len(labels)))
    dropped = 0
    for pl, cl in sim.transitions:
        if pl in idx and cl in idx:
            B[idx[pl], idx[cl]] += 1
        else:
            dropped += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        P = np.where(B.sum(axis=1, keepdims=True) > 0,
                     100.0 * B / B.sum(axis=1, keepdims=True), 0.0)
    return dict(prob=P, counts=B, labels=labels,
                n_events=int(B.sum()), n_dropped=dropped)


# ----------------------------------------------------------------------
# configuration sanity: traps that silently invalidate a comparison
# ----------------------------------------------------------------------
def validate(cfg: CERPSConfig):
    """Return a list of warnings for configurations that are internally
    inconsistent in ways that silently invalidate a mechanism comparison.

    Checks
    ------
    * ``partition="exact"`` with ``alpha != 1``: the partitioning-noise
      strength is ignored (exact mode has no noise at all).
    * ``kernel="cpg"`` with ``cpg_mode="moment"`` but
      ``cpg_sites != m0+1``: the CpG kernel then matches the beta kernel
      in mean but not in variance.
    * ``rate_mode="continuous"`` with an explicit ``K_switch`` far from
      the value that aligns the continuous boundary with the legacy
      ``x_threshold``: the classification boundary moves.
    * ``memory.enabled`` with ``driver != "stemness"``: the option is
      ignored (the integrator always uses the stemness signal).
    * ``memory.enabled`` with ``tau < dt``: the EWMA gain is clamped to 1.
    * ``rec_tran_from == 0``: transition statistics include the transient
      and are not comparable with the paper's Table 3.
    * ``space.paracrine`` with a large Fourier number: reported (the
      solver sub-steps rather than clamping), but it costs runtime.
    """
    from . import rates as rates_mod
    from .grn import GRN
    w = []
    d, e, r = cfg.division, cfg.epi, cfg.rates
    m = cfg.memory
    # P6 memory traps: a dead option and an unresolvable timescale.
    if m.enabled and m.driver != "stemness":
        w.append(f"memory.driver='{m.driver}' is IGNORED: the leaky "
                 f"integrator currently always uses the stemness signal")
    if m.enabled and m.tau < cfg.dt:
        w.append(f"memory.tau={m.tau} < dt={cfg.dt}: the EWMA gain is "
                 f"clamped to 1 (a timescale below one step is "
                 f"unresolvable)")
    if d.partition == "exact" and d.alpha != 1.0:
        w.append(f"division.alpha={d.alpha} is IGNORED because "
                 f"partition='exact' (use partition='binomial')")
    if e.kernel == "cpg" and getattr(e, "cpg_mode", "moment") == "moment" \
            and e.cpg_sites != e.m0 + 1:
        w.append(f"cpg_sites={e.cpg_sites} != m0+1={e.m0 + 1}: the CpG "
                 f"kernel matches the beta kernel in mean but NOT in "
                 f"variance")
    if r.rate_mode == "continuous" and r.K_switch is not None:
        K_auto = rates_mod.resolve_K(
            type(r)(**{f.name: getattr(r, f.name)
                       for f in fields(r) if f.name != "K_switch"}),
            GRN(cfg.grn))
        if abs(r.K_switch - K_auto) > 0.05 * max(K_auto, 1e-9):
            w.append(f"K_switch={r.K_switch} differs from the "
                     f"legacy-aligned value {K_auto:.3f}: the continuous "
                     f"classification boundary is shifted")
    if cfg.rec_tran_from == 0.0 and cfg.max_time > 4 * cfg.dt:
        w.append("rec_tran_from=0: transition statistics include the "
                 "transient (set rec_tran_from=2000 for the paper "
                 "protocol)")
    if cfg.space.paracrine:
        cfl = cfg.space.diff_rate * cfg.dt / cfg.space.grid_h ** 2
        if cfl > 1.0 / 6.0:
            w.append(f"paracrine Fourier number {cfl:.3f} > 1/6: the "
                     f"solver will take "
                     f"{int(np.ceil(cfl / (1/6) - 1e-12))} sub-steps per "
                     f"step (nominal diffusivity preserved)")
    return w


# ----------------------------------------------------------------------
# single run
# ----------------------------------------------------------------------
def run_one(cfg: CERPSConfig, seed=None, tag=None, t_from=None,
            t_to=None, keep_sim=False, warn=True):
    """Run one simulation and return a flat result dictionary."""
    if warn:
        for msg in validate(cfg):
            print(f"  [config warning] {msg}")
    if seed is not None:
        cfg = replace(cfg, seed=seed)
    sim = Simulator(cfg)
    t0 = time.time()
    sim.run()
    runtime = time.time() - t0
    ho = homeostasis(sim.records, t_from, t_to)
    tr = transitions(sim)
    out = dict(
        tag=tag or "",
        seed=int(cfg.seed),
        runtime_s=runtime,
        n_final=int(sim.n),
        max_cells=int(sim.max_cells),
        n_events=tr["n_events"],
        **{f"ho_{k}": v for k, v in ho.items()},
        **{f"tran_{pl}_{cl}": float(tr["prob"][i, j])
           for i, pl in enumerate(tr["labels"])
           for j, cl in enumerate(tr["labels"])},
        config=_cfg_summary(cfg),
    )
    if keep_sim:
        out["sim"] = sim
    return out


def _cfg_summary(cfg: CERPSConfig):
    """JSON-serialisable summary of a configuration."""
    d = {}
    for sec in _sections(cfg):
        sub = getattr(cfg, sec)
        for f in fields(sub):
            v = getattr(sub, f.name)
            if isinstance(v, (int, float, str, bool)) or v is None:
                d[f"{sec}.{f.name}"] = v
            else:
                d[f"{sec}.{f.name}"] = list(np.asarray(v).ravel()) \
                    if np.size(v) else None
    for f in fields(cfg):
        if f.name in _sections(cfg):
            continue
        v = getattr(cfg, f.name)
        if isinstance(v, (int, float, str, bool)) or v is None:
            d[f.name] = v
    return d


# ----------------------------------------------------------------------
# sweep
# ----------------------------------------------------------------------
def sweep(base: CERPSConfig, grid, seeds=(0,), out_dir=None, resume=True,
          t_from=None, t_to=None, verbose=True):
    """Cartesian sweep over ``grid`` x ``seeds``.

    ``grid`` maps config paths (``"rates.beta0"``) to lists of values.
    Checkpoints one JSON object per finished run to
    ``out_dir/results.jsonl``; with ``resume=True`` runs already present
    there are skipped.
    """
    keys = list(grid.keys())
    combos = list(itertools.product(*[grid[k] for k in keys])) \
        if keys else [()]
    jobs = [(combo, seed) for combo in combos for seed in seeds]

    done = set()
    path = None
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        path = os.path.join(out_dir, "results.jsonl")
        if resume and os.path.exists(path):
            # identify finished runs by (tag, seed)
            with open(path) as fh:
                for line in fh:
                    try:
                        r = json.loads(line)
                    except json.JSONDecodeError:
                        continue
                    done.add((r["tag"], r["seed"]))

    results = []
    fh = open(path, "a") if path else None
    try:
        for i, (combo, seed) in enumerate(jobs):
            cfg = base
            tag_parts = []
            for k, v in zip(keys, combo):
                cfg = set_path(cfg, k, v)
                tag_parts.append(f"{k}={v}")
            tag = ",".join(tag_parts)
            if (tag, seed) in done:
                if verbose:
                    print(f"[{i+1}/{len(jobs)}] skip (cached) {tag} seed={seed}")
                continue
            if verbose:
                print(f"[{i+1}/{len(jobs)}] run {tag} seed={seed} ...",
                      flush=True)
            res = run_one(cfg, seed=seed, tag=tag, t_from=t_from, t_to=t_to)
            results.append(res)
            if fh:
                fh.write(json.dumps(res) + "\n")
                fh.flush()
    finally:
        if fh:
            fh.close()
    if out_dir and resume and not results:
        results = load(out_dir)
    return results


def load(out_dir):
    """Read back a checkpointed sweep."""
    path = os.path.join(out_dir, "results.jsonl")
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                out.append(json.loads(line))
    return out


def aggregate(results):
    """Group results by ``tag`` and take mean/std across seeds.

    Returns a list of dicts with ``tag``, ``n_seeds``, ``mean`` and ``std``
    (std is the sample standard deviation across seeds, i.e. the
    Monte-Carlo error of the configuration, not of a single run).
    """
    groups = {}
    for r in results:
        groups.setdefault(r["tag"], []).append(r)
    out = []
    for tag, rs in groups.items():
        numeric = [k for k in rs[0]
                   if isinstance(rs[0][k], (int, float)) and k != "seed"]
        mean = {k: float(np.mean([r[k] for r in rs])) for k in numeric}
        std = {k: float(np.std([r[k] for r in rs], ddof=1))
               if len(rs) > 1 else 0.0 for k in numeric}
        out.append(dict(tag=tag, n_seeds=len(rs), mean=mean, std=std,
                        config=rs[0]["config"]))
    return out


def report(agg, keys=("ho_total", "ho_sc", "ho_ta1", "ho_ta2",
                      "tran_0_0", "tran_0_1", "tran_0_3", "n_events"),
           title="CERPS sweep"):
    """Print a compact table of an aggregated sweep."""
    print(f"\n=== {title} ===")
    head = f"{'configuration':44s} {'n':>3s} " + " ".join(
        f"{k.replace('ho_', '').replace('tran_', 'T'):>9s}" for k in keys)
    print(head)
    print("-" * len(head))
    for a in agg:
        row = f"{a['tag'][:44]:44s} {a['n_seeds']:3d} "
        row += " ".join(f"{a['mean'].get(k, float('nan')):9.4g}" for k in keys)
        print(row)
    print()


# ----------------------------------------------------------------------
# numerical health
# ----------------------------------------------------------------------
def health_check(cfg: CERPSConfig, t_short=None, seed=0):
    """Numerical-health report for a configuration.

    Returns a dict with

    ``cfl``          the paracrine Fourier number D*dt/h^2 (needs <= 1/6;
                     the solver sub-steps instead of clamping)
    ``n_sub``        number of field sub-steps the solver will take
    ``dt_convergence`` relative change of the homeostasis counts when dt
                     is halved (a short run at dt and dt/2).  This is the
                     number that tells a mechanism effect apart from a
                     discretisation effect: if it is comparable to the
                     effect you are measuring, refine dt.
    """
    out = {}
    sc = cfg.space
    if sc.paracrine:
        cfl = sc.diff_rate * cfg.dt / sc.grid_h ** 2
        out["cfl"] = float(cfl)
        out["n_sub"] = int(max(1, int(np.ceil(cfl / (1.0 / 6.0) - 1e-12))))
        out["cfl_ok"] = bool(out["n_sub"] >= cfl / (1.0 / 6.0) - 1e-9)
    else:
        out["cfl"] = None
        out["n_sub"] = 0
        out["cfl_ok"] = True
    # dt halving on a short run
    t_short = cfg.max_time if t_short is None else min(t_short, cfg.max_time)
    t_short = max(t_short, 4 * cfg.dt)
    a = run_one(replace(cfg, max_time=t_short), seed=seed)
    b = run_one(replace(cfg, max_time=t_short, dt=cfg.dt / 2.0), seed=seed)
    out["dt_convergence"] = {
        k: abs(a[f"ho_{k}"] - b[f"ho_{k}"]) / max(abs(a[f"ho_{k}"]), 1e-9)
        for k in ("total", "sc", "ta1", "ta2")
    }
    out["dt_convergence_max"] = max(out["dt_convergence"].values())
    out["runtime_s"] = a["runtime_s"] + b["runtime_s"]
    return out
