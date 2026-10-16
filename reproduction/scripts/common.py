"""Batch job definitions for every figure of the paper.

Each job returns a dict of numpy arrays saved to ``results/data/<name>.npz``
(checkpointing: existing files are skipped unless ``--force``).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPRO = os.path.dirname(HERE)
RESULTS = os.path.join(REPRO, "results")
DATA = os.path.join(RESULTS, "data")
FIGS = os.path.join(RESULTS, "figures")
TABLES = os.path.join(RESULTS, "tables")
for d in (RESULTS, DATA, FIGS, TABLES):
    os.makedirs(d, exist_ok=True)

sys.path.insert(0, REPRO)

from hsl2024 import Params, DEFAULT, simulate          # noqa: E402
from hsl2024.steady_state import sc_steady_state, find_steady_states, \
    a_of_u, epigenetic_fixed_point                     # noqa: E402
from hsl2024.analysis import homeostasis_stats, transition_stats  # noqa: E402

N_RUNS_STD = 30          # standard replicate count (paper)
STEPS_4000H = 16000
STEPS_2000H = 8000
STEPS_1000H = 4000
STEPS_10000H = 40000

# canonical parameter sets ------------------------------------------------
SIGMAS_FIG8 = [0.0, 0.025, 0.05, 0.075, 0.1]
M2S_FIG9 = [0.8, 0.825, 0.85, 0.875, 0.9]
ATF_FIG12 = [0.0, 0.05, 0.1, 0.15, 0.2]
ATF_FIG11 = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
US_FIG3 = [(0.5, 0.5), (0.5, 0.6), (0.6, 0.5), (0.6, 0.6)]
SIGMAS_FIG3 = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
VARIANTS_FIG10 = dict(k1=1.1, alpha1=0.44, rho1=0.11, b1=1.1, s1=0.55)


def _sc_x0(u, p=DEFAULT):
    x = sc_steady_state(u, p)
    if x is None:
        raise RuntimeError(f"no SC attractor for u={u}")
    return x


def _ta_x0(u, p=DEFAULT, which="TA1"):
    a1, a2 = a_of_u(u, p)
    pts, stable = find_steady_states(a1, a2, p)
    want = 1 if which == "TA1" else 3
    from hsl2024.steady_state import classify_state
    for x, st in zip(pts, stable):
        if st and classify_state(x) == want:
            return x
    raise RuntimeError(f"no {which} attractor for u={u}")


# ---------------------------------------------------------------- jobs ----
def job_base(seed=100):
    """sigma=0.05, m2=0.8, alpha=(0.4,0.4): anchors Fig. 4A/5, Table 3,
    Fig. 8 (sigma=0.05), Fig. 9 (m2=0.8), Fig. 12 (a_TF=0), Fig. 13."""
    p = DEFAULT
    res = simulate(p, n_runs=N_RUNS_STD, n_steps=STEPS_4000H, seed=seed,
                   landscape=True, snap_steps=[20, 200, 400, 4000])
    return dict(kind="sweep", xval=0.05,
                sys_ts=res["sys_ts"], tran=res["tran"],
                hist=res["hist"], snap_data=res["snap_data"],
                snap_n=res["snap_n"], stop_step=res["stop_step"])


def job_fig8(sigma, seed=200):
    p = DEFAULT.with_(sigma=sigma)
    res = simulate(p, n_runs=N_RUNS_STD, n_steps=STEPS_4000H, seed=seed,
                   landscape=True)
    return dict(kind="sweep", xval=sigma, sys_ts=res["sys_ts"],
                tran=res["tran"], hist=res["hist"], stop_step=res["stop_step"])


def job_fig9(m2, seed=300):
    p = DEFAULT.with_(m2=m2)
    res = simulate(p, n_runs=N_RUNS_STD, n_steps=STEPS_4000H, seed=seed,
                   landscape=True)
    return dict(kind="sweep", xval=m2, sys_ts=res["sys_ts"],
                tran=res["tran"], hist=res["hist"], stop_step=res["stop_step"])


def job_fig10(variant, m2, seed=400):
    p = DEFAULT.with_(m2=m2, **{variant: VARIANTS_FIG10[variant]})
    res = simulate(p, n_runs=N_RUNS_STD, n_steps=STEPS_4000H, seed=seed)
    return dict(kind="sweep", xval=m2, sys_ts=res["sys_ts"],
                tran=res["tran"], stop_step=res["stop_step"])


def job_fig12(atf, seed=500):
    p = DEFAULT.with_(a_tf=atf)
    res = simulate(p, n_runs=N_RUNS_STD, n_steps=STEPS_4000H, seed=seed,
                   landscape=True)
    return dict(kind="sweep", xval=atf, sys_ts=res["sys_ts"],
                tran=res["tran"], hist=res["hist"], stop_step=res["stop_step"])


def job_fig4(alpha_pair, seed=600):
    p = DEFAULT.with_(alpha1=alpha_pair[0], alpha2=alpha_pair[1])
    res = simulate(p, n_runs=N_RUNS_STD, n_steps=STEPS_2000H, seed=seed)
    return dict(kind="fig4", xval=alpha_pair[0], sys_ts=res["sys_ts"],
                stop_step=res["stop_step"])


def job_fig7(sigma, m2, seed=700):
    p = DEFAULT.with_(sigma=sigma, m2=m2)
    res = simulate(p, n_runs=10, n_steps=STEPS_4000H, seed=seed)
    return dict(kind="fig7", xval=sigma, sys_ts=res["sys_ts"],
                stop_step=res["stop_step"])


def job_fig3(u_pair, sigma, seed=800):
    p = DEFAULT.with_(sigma=sigma)
    x0 = _sc_x0(u_pair, p)
    res = simulate(p, n_runs=100, n_steps=STEPS_1000H, seed=seed,
                   init_mode=1, fixed_u=True, u_init=u_pair, x_init=x0,
                   rec_tran_from=10 ** 9, cap_per_run=20000)
    return dict(kind="fig3", xval=sigma, sys_ts=res["sys_ts"],
                stop_step=res["stop_step"],
                u=u_pair, x0=x0)


def job_fig11(which, atf, seed=900):
    p = DEFAULT.with_(a_tf=atf, n_init=500)
    vstar = epigenetic_fixed_point(p)
    u = (vstar, vstar)
    x0 = _ta_x0(u, p, which)
    res = simulate(p, n_runs=100, n_steps=STEPS_10000H, seed=seed,
                   init_mode=2, u_init=u, x_init=x0,
                   rec_tran_from=10 ** 9, stop_sc_target=100,
                   cap_per_run=20000)
    return dict(kind="fig11", xval=atf, sys_ts=res["sys_ts"],
                stop_step=res["stop_step"], stop_kind=res["stop_kind"],
                u=u, x0=x0)


def job_fig11c(atf=0.2, seed=950):
    """Un-stopped dynamics from 100% TA2 for Fig. 11C."""
    p = DEFAULT.with_(a_tf=atf, n_init=500)
    vstar = epigenetic_fixed_point(p)
    u = (vstar, vstar)
    x0 = _ta_x0(u, p, "TA2")
    res = simulate(p, n_runs=10, n_steps=STEPS_1000H * 2, seed=seed,
                   init_mode=2, u_init=u, x_init=x0,
                   rec_tran_from=10 ** 9, cap_per_run=20000)
    return dict(kind="fig11c", xval=atf, sys_ts=res["sys_ts"],
                stop_step=res["stop_step"])


def build_jobs():
    from functools import partial
    jobs = {}

    def add(name, fn, *a, **kw):
        if a or kw:
            jobs[name] = partial(fn, *a, **kw)
        else:
            jobs[name] = fn

    add("base", job_base)
    for s in SIGMAS_FIG8:
        if abs(s - 0.05) < 1e-12:
            continue                    # base job
        add(f"fig8_s{s:g}", job_fig8, s)
    for m in M2S_FIG9:
        if abs(m - 0.8) < 1e-12:
            continue                    # base job
        add(f"fig9_m{m:g}", job_fig9, m)
    for a in ATF_FIG12:
        if abs(a) < 1e-12:
            continue                    # base job
        add(f"fig12_tf{a:g}", job_fig12, a)
    for key, val in VARIANTS_FIG10.items():
        for m in M2S_FIG9:
            add(f"fig10_{key}_m{m:g}", job_fig10, key, m)
    for ap in [(0.4, 0.4), (0.4, 1.0), (1.0, 0.4), (1.0, 1.0)]:
        add(f"fig4_a{ap[0]:g}_{ap[1]:g}", job_fig4, ap)
    for s, m in [(0.20, 0.80), (0.05, 0.75), (0.05, 1.00)]:
        add(f"fig7_s{s:g}_m{m:g}", job_fig7, s, m)
    for (u1, u2) in US_FIG3:
        for s in SIGMAS_FIG3:
            add(f"fig3_u{u1:g}_{u2:g}_s{s:g}", job_fig3, (u1, u2), s)
    for which in ("TA1", "TA2"):
        for a in ATF_FIG11:
            add(f"fig11_{which}_tf{a:g}", job_fig11, which, a)
    add("fig11c_tf0.2", job_fig11c)
    for r in range(30):
        add(f"ref_run{r:02d}", _job_ref_run, r)
    return jobs


def _job_ref_run(r, seed=2024):
    """Single replicate with the independent numpy reference simulator
    (cross-validation of the numba implementation)."""
    from hsl2024.reference_sim import simulate_reference
    p = DEFAULT
    res = simulate_reference(p, n_runs=1, n_steps=STEPS_4000H, seed=seed + r,
                             rec_tran_from=8000)
    return dict(kind="ref", xval=float(r), sys_ts=res["sys_ts"],
                tran=res["tran"], stop_step=np.full(1, -1, np.int32))


def job_path(name):
    return os.path.join(DATA, f"{name}.npz")
