"""Deterministic steady-state analysis of the 2-gene GRN (Figs. 2A, 2B)."""
import numpy as np
from scipy.optimize import root

from .params import Params


def rhs(x, a1, a2, p: Params):
    x1, x2 = x
    s1n = p.s1 ** p.n
    s2n = p.s2 ** p.n
    f1 = a1 * (p.rho1 + (1 - p.rho1) * x1 ** p.n / (s1n + x1 ** p.n)) \
        + p.b1 * s2n / (s2n + x2 ** p.n) - p.k1 * x1
    f2 = a2 * (p.rho2 + (1 - p.rho2) * x2 ** p.n / (s2n + x2 ** p.n)) \
        + p.b2 * s1n / (s1n + x1 ** p.n) - p.k2 * x2
    return np.array([f1, f2])


def jacobian(x, a1, a2, p: Params):
    x1, x2 = x
    n = p.n
    s1n = p.s1 ** n
    s2n = p.s2 ** n
    x1n = x1 ** n
    x2n = x2 ** n
    # d/dx1 of a1*(rho + (1-rho)*x1^n/(s1n+x1^n))
    d_self_1 = a1 * (1 - p.rho1) * n * x1 ** (n - 1) * s1n / (s1n + x1n) ** 2
    d_rep_by_x2 = -p.b1 * s2n * n * x2 ** (n - 1) / (s2n + x2n) ** 2
    d_self_2 = a2 * (1 - p.rho2) * n * x2 ** (n - 1) * s2n / (s2n + x2n) ** 2
    d_rep_by_x1 = -p.b2 * s1n * n * x1 ** (n - 1) / (s1n + x1n) ** 2
    return np.array([
        [d_self_1 - p.k1, d_rep_by_x2],
        [d_rep_by_x1, d_self_2 - p.k2],
    ])


def find_steady_states(a1, a2, p: Params, seed_grid=(5, 5), xmax=4.0):
    """All steady states for given rates; returns (points, stable_mask).

    Multistart Newton (scipy hybr) from a coarse grid, then dedupe.
    """
    seeds = []
    for xi in np.linspace(0.02, xmax, seed_grid[0]):
        for xj in np.linspace(0.02, xmax, seed_grid[1]):
            seeds.append([xi, xj])
    sols = []
    for x0 in seeds:
        sol = root(lambda x: rhs(x, a1, a2, p), np.array(x0), tol=1e-11)
        if sol.success and np.all(sol.x > -1e-6):
            x = np.maximum(sol.x, 0.0)
            if not any(np.allclose(x, s, atol=2e-3) for s in sols):
                sols.append(x)
    pts = np.array(sols) if sols else np.zeros((0, 2))
    if pts.shape[0] == 0:
        return pts, np.zeros(0, bool)
    stable = np.zeros(pts.shape[0], bool)
    for i, x in enumerate(pts):
        J = jacobian(x, a1, a2, p)
        ev = np.linalg.eigvals(J)
        stable[i] = np.all(ev.real < 0)
    return pts, stable


def classify_state(x):
    """0=SC, 1=TA1, 3=TA2, -1=unclassified (both low)."""
    if x[0] > 0.5 and x[1] > 0.5:
        return 0
    if x[0] > 0.5 and x[1] < 0.5:
        return 1
    if x[0] < 0.5 and x[1] > 0.5:
        return 3
    return -1


def a_of_u(u, p: Params):
    a1 = p.alpha1 * np.exp(p.lam1 * u[0])
    a2 = p.alpha2 * np.exp(p.lam2 * u[1])
    return a1, a2


def sc_steady_state(u, p: Params):
    """SC-attractor (both x > 0.5) for the given epigenetic state, or None."""
    a1, a2 = a_of_u(u, p)
    pts, stable = find_steady_states(a1, a2, p)
    for x, st in zip(pts, stable):
        if st and classify_state(x) == 0:
            return x
    return None


def epigenetic_fixed_point(p: Params):
    """v* with phi(v*) = v* (stationary mean epigenetic state of TA cells)."""
    v = 0.3
    for _ in range(200):
        g = (p.m3 * v) ** p.m4
        phi = p.m1 + p.m2 * g / (1.0 + g)
        if abs(phi - v) < 1e-12:
            break
        v = phi
    return v
