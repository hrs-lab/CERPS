"""Regression tests for the 2026-10-21 deep audit.

Covers the three code changes, all of which are OFF by default so the
validated Table-3 baseline is untouched:

  * GRNConfig.ou_exact        -- exact OU transition (stationary var == 1)
  * SpaceConfig.field_interp  -- trilinear sampling + conservative source
  * batch.grn_stiffness       -- explicit-Euler stiffness diagnostic

and asserts the defaults still reproduce the legacy paths.
"""
import numpy as np
import pytest

from cerps.config import GRNConfig, SpaceConfig, CERPSConfig, EpiConfig
from cerps.grn import GRN
from cerps.spatial import Space
from cerps import batch
from cerps import division as div_mod
from cerps import epigenetics as epi


# ---------------------------------------------------------------- OU
def _ou_var(exact, zeta=1.0, dt=0.25, n_ch=4000, n_steps=4000, seed=0):
    cfg = GRNConfig(ou_exact=exact)
    g = GRN(cfg)
    rng = np.random.default_rng(seed)
    eta = np.zeros((n_ch, 1))
    for _ in range(n_steps):
        eta = g.ou_step(eta, dt, rng)
    return float(eta.var())


def test_ou_exact_stationary_variance_is_one():
    v = _ou_var(True)
    assert abs(v - 1.0) < 0.06, v


def test_ou_euler_stationary_variance_is_biased():
    v = _ou_var(False)
    theory = 1.0 / (1.0 - 0.25 / 2.0)          # dt=0.25, zeta=1
    assert abs(v - theory) < 0.10, (v, theory)
    assert v > 1.05                              # clearly biased upward


def test_ou_default_is_euler():
    assert GRNConfig().ou_exact is False


# ---------------------------------------------------------------- field
def _linear_space(interp):
    cfg = SpaceConfig(enabled=True, paracrine=True, domain_radius=10.0,
                      grid_h=1.0, field_interp=interp)
    sp = Space(cfg)
    n = sp._grid_shape[0]
    ii, jj, kk = np.meshgrid(np.arange(n), np.arange(n), np.arange(n),
                             indexing="ij")
    # linear in grid index -> linear in physical coordinates
    sp.field = (0.3 * ii + 0.2 * jj - 0.1 * kk + 0.5).astype(float)
    return sp, n


def test_trilinear_sampling_exact_for_linear_field():
    sp, n = _linear_space("linear")
    rng = np.random.default_rng(1)
    pos = rng.uniform(-9.9, 9.9, size=(2000, 3))
    got = sp.sample(pos)
    f = (pos - sp._grid_origin) / sp.cfg.grid_h
    exact = 0.3 * f[:, 0] + 0.2 * f[:, 1] - 0.1 * f[:, 2] + 0.5
    assert np.max(np.abs(got - exact)) < 1e-9


def test_nearest_sampling_is_piecewise_constant_and_inaccurate():
    sp, n = _linear_space("nearest")
    rng = np.random.default_rng(2)
    pos = rng.uniform(-9.9, 9.9, size=(2000, 3))
    got = sp.sample(pos)
    f = (pos - sp._grid_origin) / sp.cfg.grid_h
    exact = 0.3 * f[:, 0] + 0.2 * f[:, 1] - 0.1 * f[:, 2] + 0.5
    assert np.mean(np.abs(got - exact)) > 0.05


def test_trilinear_weights_sum_to_one():
    sp, n = _linear_space("linear")
    rng = np.random.default_rng(3)
    pos = rng.uniform(-9.9, 9.9, size=(500, 3))
    i0, i1, w0, w1 = sp._tri(pos)
    # per axis w0+w1 == 1 by construction; the 8-point product sums to 1
    total = np.zeros(len(pos))
    for aa in (0, 1):
        wa = w1[:, 0] if aa else w0[:, 0]
        for bb in (0, 1):
            wb = w1[:, 1] if bb else w0[:, 1]
            for cc in (0, 1):
                wc = w1[:, 2] if cc else w0[:, 2]
                total += wa * wb * wc
    assert np.max(np.abs(total - 1.0)) < 1e-12


def test_field_default_is_nearest():
    assert SpaceConfig().field_interp == "nearest"


# ---------------------------------------------------------------- stiffness
def test_grn_stiffness_default_positive():
    st = batch.grn_stiffness(CERPSConfig())
    assert st > 0


def test_grn_stiffness_flags_unstable_dt():
    cfg = CERPSConfig(dt=2.0, grn=GRNConfig(alpha=(3.0, 3.0),
                                            lam=(1.9, 1.9)))
    assert batch.grn_stiffness(cfg) > 1.0
    msgs = batch.validate(cfg)
    assert any("stiffness" in m for m in msgs)


# ---------------------------------------------------------------- asymmetry
def test_ou_exact_asymmetric_zeta_per_gene():
    """Exact OU gives unit stationary variance for EACH gene, even when
    zeta_1 != zeta_2 (Euler-Maruyama would give 1/(1-dt/2 zeta_i))."""
    cfg = GRNConfig(ou_exact=True, zeta=(0.25, 4.0))
    g = GRN(cfg)
    rng = np.random.default_rng(11)
    eta = np.zeros((6000, 2))
    for _ in range(4000):
        eta = g.ou_step(eta, 0.25, rng)
    v = eta.var(axis=0)
    assert np.all(np.abs(v - 1.0) < 0.10), v
    # and the Euler path is strongly biased for the fast gene
    cfg_e = GRNConfig(ou_exact=False, zeta=(0.25, 4.0))
    ge = GRN(cfg_e)
    eta = np.zeros((6000, 2))
    for _ in range(4000):
        eta = ge.ou_step(eta, 0.25, rng)
    ve = eta.var(axis=0)
    assert ve[0] > 1.5 and ve[1] < 1.2, ve


def test_grn_stiffness_is_per_gene():
    """An asymmetric alpha must be reflected by the stiffest gene, not an
    averaged scalar."""
    sym = CERPSConfig(grn=GRNConfig(alpha=(0.4, 0.4)))
    asym = CERPSConfig(grn=GRNConfig(alpha=(0.8, 0.4)))
    assert batch.grn_stiffness(asym) > batch.grn_stiffness(sym)
    # doubling alpha_1 roughly doubles the (a*(1-rho)*n/(4s) - k) part
    assert batch.grn_stiffness(asym) > 1.4 * batch.grn_stiffness(sym) - 0.5


def test_partition_asymmetric_f_matches_theory():
    from cerps.config import DivisionConfig
    cfg = DivisionConfig(partition="binomial", volume=200.0, asymmetry=0.3,
                         alpha=1.0)
    rng = np.random.default_rng(12)
    x = np.ones((200000, 1))
    x1, x2 = div_mod.partition(x, cfg, rng)
    # theory: Var(x1)=x(1-f)/(fV), Var(x2)=xf/((1-f)V)
    th1 = 1.0 * (1 - 0.3) / (0.3 * 200.0)
    th2 = 1.0 * 0.3 / (0.7 * 200.0)
    assert abs(x1.var() - th1) < 0.1 * th1
    assert abs(x2.var() - th2) < 0.1 * th2


def test_product_equals_additive_for_g2_asymmetric():
    """For G=2 the product and additive forms coincide even when s1 != s2."""
    x = np.array([[1.4, 0.2], [0.3, 1.1], [0.8, 0.9]])
    u = np.full_like(x, 0.3)
    eta = np.zeros_like(x)
    gp = GRN(GRNConfig(n_genes=2, regulation="product", s_self=(0.55, 0.5),
                       sigma=(0.0, 0.0)))
    ga = GRN(GRNConfig(n_genes=2, regulation="additive", s_self=(0.55, 0.5),
                       sigma=(0.0, 0.0)))
    dp = gp.step(x, u, eta, 0.25) - x
    da = ga.step(x, u, eta, 0.25) - x
    assert np.max(np.abs(dp - da)) < 1e-12


# ---------------------------------------------------------------- sister corr
def _sister_sample(rho, n=200000, v=0.4):
    cfg = EpiConfig(kernel="beta", sister_corr=rho, m0=60.0)
    rng = np.random.default_rng(21)
    d1, d2 = epi.inherit(np.full((n, 1), v), cfg, rng)
    return d1[:, 0], d2[:, 0], float(epi.phi_map(np.array([v]), cfg)[0])


def test_sister_corr_preserves_beta_marginals():
    d1, d2, phi = _sister_sample(0.6)
    var_th = phi * (1 - phi) / 61.0
    assert abs(d1.mean() - phi) < 0.01
    assert abs(d1.var() - var_th) / var_th < 0.1
    assert abs(d2.mean() - phi) < 0.01


def test_sister_corr_increases_sister_correlation():
    c0 = np.corrcoef(*_sister_sample(0.0)[:2])[0, 1]
    c2 = np.corrcoef(*_sister_sample(0.2)[:2])[0, 1]
    c6 = np.corrcoef(*_sister_sample(0.6)[:2])[0, 1]
    assert abs(c0) < 0.02
    assert 0.0 < c2 < c6 < 0.95


def test_sister_corr_gaussian_raises():
    cfg = EpiConfig(kernel="gaussian", sister_corr=0.5)
    rng = np.random.default_rng(22)
    with pytest.raises(ValueError):
        epi.inherit(np.full((10, 1), 0.4), cfg, rng)


def test_sister_corr_default_zero():
    assert EpiConfig().sister_corr == 0.0


# ---------------------------------------------------------------- sphere mask
def test_sphere_mask_zeros_field_outside_sphere():
    cfg = SpaceConfig(enabled=True, paracrine=True, domain_radius=10.0,
                      grid_h=1.0, sphere_mask=True, diff_rate=1.0,
                      decay_rate=0.1, secretion=1.0)
    sp = Space(cfg)
    sp.init_positions(1, np.random.default_rng(23))
    sp.pos[:] = 0.0                       # one cell at the centre
    sp.field_step(np.array([1.0]), dt=0.1)
    m = sp._sphere_mask
    assert m is not None
    assert np.all(sp.field[~m] == 0.0)
    assert sp.field[m].max() > 0.0


def test_sphere_mask_default_false():
    assert SpaceConfig().sphere_mask is False
