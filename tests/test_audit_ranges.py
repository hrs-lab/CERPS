"""Range-scan regression tests for audit_2026-10-21 §8 (open item #5).

Pins the two range-sensitive findings of the audit over parameter grids, so a
future refactor cannot silently change them:

  * §8.1  OU stationary variance vs ``zeta`` at dt=0.25: the exact transition
    gives 1 for every zeta, while Euler-Maruyama gives the analytic bias
    ``1/(1 - dt/(2 zeta))`` (+100% at zeta=0.25 ... +3.2% at zeta=4.0).
  * §8.2  explicit-Euler GRN stiffness ``dt*lambda_i`` for the paper's
    parameter sets, via the per-gene linearised upper bound
    ``lambda_i = alpha_i e^{lam_i} (1-rho_i) n/(4 s_i) - k_i``.

Pure NumPy, no simulation, so this file is fast enough for every CI run.
"""
import numpy as np
import pytest

from cerps.config import GRNConfig, CERPSConfig
from cerps.grn import GRN
from cerps import batch

DT = 0.25
ZETA_GRID = [0.25, 0.5, 1.0, 2.0, 4.0]


def _ou_var(exact, zeta, dt=DT, n_ch=3000, n_steps=3000, seed=0):
    cfg = GRNConfig(ou_exact=exact, zeta=(zeta, zeta))
    g = GRN(cfg)
    rng = np.random.default_rng(seed)
    eta = np.zeros((n_ch, 1))
    for _ in range(n_steps):
        eta = g.ou_step(eta, dt, rng)
    return float(eta.var())


@pytest.mark.parametrize("zeta", ZETA_GRID)
def test_ou_exact_variance_is_one_over_zeta_grid(zeta):
    """§8.1: exact OU -> unit stationary variance for every relaxation time."""
    v = _ou_var(True, zeta)
    assert abs(v - 1.0) < 0.06, (zeta, v)


@pytest.mark.parametrize("zeta", ZETA_GRID)
def test_ou_euler_variance_matches_analytic_bias_over_zeta_grid(zeta):
    """§8.1: Euler bias 1/(1-dt/2 zeta) reproduced at every zeta (+100% ... +3.2%)."""
    v = _ou_var(False, zeta)
    bias = 1.0 / (1.0 - DT / (2.0 * zeta))
    assert abs(v / bias - 1.0) < 0.05, (zeta, v, bias)


def _stiff(**kw):
    return batch.grn_stiffness(CERPSConfig(grn=GRNConfig(**kw)))


@pytest.mark.parametrize("kw,expected", [
    (dict(), 0.9535),                            # default (symmetric)
    (dict(alpha=(1.0, 1.0)), 2.759),             # Fig.4 alpha=1.0
    (dict(s_self=(0.3, 0.3)), 1.756),            # small self-threshold
    (dict(s_self=(0.55, 0.5)), 0.9535),          # Fig10 s1=0.55 (gene 2 binds)
    (dict(alpha=(0.44, 0.4)), 1.074),            # Fig10 alpha1=0.44 (crosses 1)
    (dict(k_deg=(1.1, 1.0)), 0.9535),            # Fig10 k1=1.1 (gene 2 binds)
    (dict(rho=(0.11, 0.1)), 0.9535),             # Fig10 rho1=0.11 (gene 2 binds)
])
def test_grn_stiffness_matches_audit_range_scan(kw, expected):
    """§8.2: per-gene stiffness bound reproduces the audit table exactly."""
    assert _stiff(**kw) == pytest.approx(expected, rel=2e-3)


def test_grn_stiffness_crosses_unity_only_where_audit_says():
    """§8.2: only alpha=1.0 / s=0.3 / alpha1=0.44 exceed dt*lambda = 1."""
    over = {name: _stiff(**kw) for name, kw in [
        ("default", dict()), ("alpha1.0", dict(alpha=(1.0, 1.0))),
        ("s0.3", dict(s_self=(0.3, 0.3))), ("s10.55", dict(s_self=(0.55, 0.5))),
        ("alpha10.44", dict(alpha=(0.44, 0.4))),
        ("k11.1", dict(k_deg=(1.1, 1.0))), ("rho10.11", dict(rho=(0.11, 0.1))),
    ] if _stiff(**kw) > 1.0}
    assert set(over) == {"alpha1.0", "s0.3", "alpha10.44"}, over
