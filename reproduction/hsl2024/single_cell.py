"""Single-cell gene dynamics (no cell cycle): Figs. 2C and 2D.

Vectorised over many independent cells (grid points x replicates).
"""
import numpy as np

from .params import Params


def single_cell_step(x, u, eta, p: Params, dt, rng):
    """One forward-Euler step of Eqs. (1)-(4) for many cells at once.

    x, eta: (M, 2); u: (M, 2) fixed.  Returns updated (x, eta).
    """
    eta = eta - eta * (dt / p.zeta1) + np.sqrt(2.0 * dt / p.zeta1) * rng.standard_normal(size=(x.shape[0], 2))
    a = p.alpha1 * np.exp(p.lam1 * u[:, 0]) * np.exp(p.sigma * eta[:, 0] - 0.5 * p.sigma ** 2)
    a = np.stack([a,
                  p.alpha2 * np.exp(p.lam2 * u[:, 1]) * np.exp(p.sigma * eta[:, 1] - 0.5 * p.sigma ** 2)],
                 axis=1)
    a[:, 0] += p.a_tf
    s1n = p.s1 ** p.n
    s2n = p.s2 ** p.n
    x1p = x[:, 0] ** p.n
    x2p = x[:, 1] ** p.n
    f1 = a[:, 0] * (p.rho1 + (1 - p.rho1) * x1p / (s1n + x1p)) \
        + p.b1 * s2n / (s2n + x2p) - p.k1 * x[:, 0]
    f2 = a[:, 1] * (p.rho2 + (1 - p.rho2) * x2p / (s2n + x2p)) \
        + p.b2 * s1n / (s1n + x1p) - p.k2 * x[:, 1]
    x = x + np.stack([f1, f2], axis=1) * dt
    return x, eta


def classify(x):
    """SC=0, TA1=1, TA2=3, -1 = keep previous (both low)."""
    t = np.full(x.shape[0], -1, np.int8)
    t[(x[:, 0] > 0.5) & (x[:, 1] > 0.5)] = 0
    t[(x[:, 0] > 0.5) & (x[:, 1] < 0.5)] = 1
    t[(x[:, 0] < 0.5) & (x[:, 1] > 0.5)] = 3
    return t


def sc_durations(u_grid, x0_grid, p: Params, t_max=100.0, n_reps=100,
                 seed=0):
    """Average time a cell stays SC, for each (u1, u2) in ``u_grid``.

    u_grid: (G, 2); x0_grid: (G, 2) initial expression (SC steady state,
    or anything else).  Each grid point is replicated ``n_reps`` times with
    independent noise; duration = first exit from SC classification,
    censored at ``t_max``.  Returns durations (G,) in hours.
    """
    rng = np.random.default_rng(seed)
    G = u_grid.shape[0]
    M = G * n_reps
    u = np.repeat(u_grid, n_reps, axis=0)
    x = np.repeat(x0_grid, n_reps, axis=0)
    eta = np.zeros((M, 2))
    cur_type = classify(x)
    exited = cur_type != 0          # already not SC at t=0 -> duration 0
    t_exit = np.zeros(M)
    n_steps = int(round(t_max / p.dt))
    for s in range(n_steps):
        if exited.all():
            break
        x, eta = single_cell_step(x, u, eta, p, p.dt, rng)
        t = (s + 1) * p.dt
        new_type = classify(x)
        newly = (~exited) & (new_type != 0)
        t_exit[newly] = t
        exited = exited | newly
    t_exit[~exited] = t_max
    return t_exit.reshape(G, n_reps).mean(axis=1)


def single_trajectory(u, x0, p: Params, t_max=100.0, seed=0, dt=None):
    """One sample trajectory of (x1(t), x2(t)) for Fig. 2C."""
    rng = np.random.default_rng(seed)
    dt = p.dt if dt is None else dt
    n_steps = int(round(t_max / dt))
    x = np.asarray(x0, float)[None, :].copy()
    u = np.asarray(u, float)[None, :]
    eta = np.zeros((1, 2))
    ts = [0.0]
    xs = [x[0].copy()]
    for s in range(n_steps):
        x, eta = single_cell_step(x, u, eta, p, dt, rng)
        ts.append((s + 1) * dt)
        xs.append(x[0].copy())
    return np.array(ts), np.array(xs)


def fixed_u_end_ratios(u, x0, p: Params, t_max=1000.0, n_runs=100,
                       n_cells=100, seed=0):
    """Fig. 3 protocol: n_runs x n_cells cells start at the SC attractor
    with the epigenetic state fixed at u; pure gene dynamics (no cell
    cycle) for t_max; returns end-state ratios (SC, TA1, TA2) in percent,
    averaged over runs (all cells are simulated in one vectorised batch).

    The phenotype is tracked with the same keep-previous rule as the
    population model (initial label: SC, as all cells start as stem
    cells).
    """
    rng = np.random.default_rng(seed)
    M = n_runs * n_cells
    u_arr = np.tile(np.asarray(u, float), (M, 1))
    x = np.tile(np.asarray(x0, float), (M, 1))
    eta = np.zeros((M, 2))
    cur = np.zeros(M, np.int8)          # all start as SC
    n_steps = int(round(t_max / p.dt))
    for s in range(n_steps):
        x, eta = single_cell_step(x, u_arr, eta, p, p.dt, rng)
        new = classify(x)
        keep = new < 0
        cur = np.where(keep, cur, new)
    # end-state ratios per run
    counts = np.stack([np.mean(cur.reshape(n_runs, n_cells) == t, axis=1)
                       for t in (0, 1, 3)], axis=1)
    return counts * 100.0               # (n_runs, 3)
