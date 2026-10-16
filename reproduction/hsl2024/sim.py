"""High-level driver around the numba simulator."""
import numpy as np

from .params import Params
from .numba_sim import run_simulation


def simulate(
    p: Params,
    n_runs: int = 30,
    n_steps: int | None = None,
    seed: int = 0,
    init_mode: int = 0,
    fixed_u: bool = False,
    u_init=None,               # (n_runs, 2) for init_mode != 0
    x_init=None,               # (n_runs, 2)
    rec_tran_from: int = 8000, # record transitions for step+1 > this (paper: 8000)
    stop_sc_target: int = 0,   # >0: stop run when SC count >= target (Fig. 11)
    landscape: bool = False,
    tb_stride: int = 40,
    n_tb: int = 100,
    n_zb: int = 80,
    snap_steps=(),
    cap_per_run: int = 12000,
) -> dict:
    """Run ``n_runs`` replicates of the hybrid model with parameters ``p``.

    init_mode: 0 = paper default (x ~ U(0,3)^2, u ~ U(0,1)^2),
               1 = all cells at (x_init, u_init) (Fig. 3: SC attractor),
               2 = all cells at (x_init, u_init) (Fig. 11: TA attractor).
    """
    if n_steps is None:
        n_steps = int(round(p.max_time / p.dt))
    n_steps = int(n_steps)

    sys_ts = np.zeros((n_steps, n_runs, 4), np.int32)
    tran = np.zeros((n_runs, 3, 3), np.int64)
    stop_step = np.full(n_runs, -1, np.int32)
    stop_kind = np.zeros(n_runs, np.int8)
    max_cells = np.zeros(n_runs, np.int32)
    out_status = np.zeros(1, np.int32)

    if landscape:
        hist = np.zeros((n_runs, n_tb, n_zb), np.int32)
    else:
        hist = np.zeros((0, 0, 0), np.int32)

    snap_steps = np.asarray(snap_steps, np.int32)
    if snap_steps.size:
        snap_data = np.zeros((snap_steps.size, cap_per_run, 4), np.float32)
    else:
        snap_data = np.zeros((0, 0, 0), np.float32)
    snap_n = np.zeros(max(snap_steps.size, 1), np.int32)

    if init_mode == 0:
        u_init = np.zeros((n_runs, 2))
        x_init = np.zeros((n_runs, 2))
    else:
        u_init = np.ascontiguousarray(np.broadcast_to(
            np.asarray(u_init, float), (n_runs, 2)).copy())
        x_init = np.ascontiguousarray(np.broadcast_to(
            np.asarray(x_init, float), (n_runs, 2)).copy())

    cap_total = int(cap_per_run * n_runs)

    n_final = run_simulation(
        p.b1, p.b2, p.k1, p.k2, p.s1, p.s2, p.rho1, p.rho2, p.n,
        p.alpha1, p.alpha2, p.lam1, p.lam2, p.zeta1, p.zeta2, p.sigma, p.a_tf,
        p.mu0, p.kappa0, p.beta0, p.tau0, p.theta0, p.max_div_ta,
        p.m0, p.m1, p.m2, p.m3, p.m4,
        p.dt, n_steps, n_runs, p.n_init, init_mode, int(fixed_u),
        int(rec_tran_from), cap_total, int(stop_sc_target), int(seed),
        sys_ts, tran, stop_step, stop_kind, max_cells,
        int(landscape), hist, int(tb_stride),
        snap_steps, snap_data, snap_n,
        u_init, x_init, out_status,
    )
    if out_status[0] == 1:
        raise RuntimeError(
            f"capacity exceeded (cap_total={cap_total}); raise cap_per_run "
            f"(max cells seen: {max_cells.max()})")
    out = dict(
        params=p, n_runs=n_runs, n_steps=n_steps, dt=p.dt,
        sys_ts=sys_ts, tran=tran, stop_step=stop_step, stop_kind=stop_kind,
        max_cells=max_cells, n_final=n_final,
    )
    if landscape:
        out["hist"] = hist
        out["tb_stride"] = tb_stride
    if snap_steps.size:
        out["snap_data"] = snap_data
        out["snap_n"] = snap_n
        out["snap_steps"] = snap_steps
    return out
