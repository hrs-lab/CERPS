"""CERPS unit and integration tests.

Run with:  python -m pytest tests/ -v   (or python tests/test_cerps.py)
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


# ----------------------------------------------------------------------
def test_grn_legacy_equivalence():
    """R3: the N-gene GRN reduces to Eq. (1) for 2 genes."""
    from hsl2024.single_cell import single_cell_step
    from hsl2024.params import Params
    p = Params()
    grn = GRN(GRNConfig())
    rng = np.random.default_rng(0)
    x_c = np.array([[1.3, 0.7]])
    u_c = np.array([[0.42, 0.61]])
    # same noise sequence
    eta0 = np.zeros((1, 2))
    x1 = x_c.copy()
    for _ in range(200):
        z = rng.standard_normal((1, 2))
        # CERPS
        eta_c = eta0 - eta0 * (p.dt / p.zeta1) + np.sqrt(
            2 * p.dt / p.zeta1) * z
        x1 = grn.step(x1, u_c, eta_c, p.dt)
        eta0 = eta_c
    # reference (hsl2024) with the same noise
    eta0 = np.zeros((1, 2))
    x2 = x_c.copy()
    for _ in range(200):
        z = rng.standard_normal((1, 2))   # different stream; only compare
        x2, eta0 = single_cell_step(x2, u_c, eta0, p, p.dt,
                                    np.random.default_rng(1))
    # deterministic comparison instead (no noise):
    grn_det = GRN(GRNConfig(sigma=(0.0, 0.0)))
    rng = np.random.default_rng(0)
    x1 = x_c.copy()
    e = np.zeros((1, 2))
    for _ in range(400):
        x1 = grn_det.step(x1, u_c, e, p.dt)
    # analytic RHS of the paper's Eq. (1)
    def rhs(x):
        a1 = p.alpha1 * np.exp(p.lam1 * 0.42)
        a2 = p.alpha2 * np.exp(p.lam2 * 0.61)
        f1 = a1 * (p.rho1 + (1 - p.rho1) * x[0, 0] ** 4 / (0.5 ** 4 + x[0, 0] ** 4)) \
            + p.b1 * 0.5 ** 4 / (0.5 ** 4 + x[0, 1] ** 4) - p.k1 * x[0, 0]
        f2 = a2 * (p.rho2 + (1 - p.rho2) * x[0, 1] ** 4 / (0.5 ** 4 + x[0, 1] ** 4)) \
            + p.b2 * 0.5 ** 4 / (0.5 ** 4 + x[0, 0] ** 4) - p.k2 * x[0, 1]
        return np.array([[f1, f2]])
    x2 = x_c.copy()
    for _ in range(400):
        x2 = x2 + rhs(x2) * p.dt
    check("grn: N-gene update == paper Eq.(1) (2-gene)",
          np.allclose(x1, x2, atol=1e-9), f"{x1} vs {x2}")


def test_grn_steady_states():
    """2-gene GRN reproduces the tristable SC/TA1/TA2 attractors
    (u = 0.5 -> a = 0.4*exp(1.9*0.5) = 1.0346)."""
    grn = GRN(GRNConfig(sigma=(0, 0)))
    for seed, expect in [((1.0, 1.0), (1.038, 1.038)),
                         ((2.5, 0.1), (2.029, 0.109)),
                         ((0.1, 2.5), (0.109, 2.029))]:
        x = np.array([seed], float)
        u = np.full((1, 2), 0.5)
        e = np.zeros((1, 2))
        for _ in range(4000):
            x = grn.step(x, u, e, 0.05)
        check(f"grn: attractor from {seed} -> {expect}",
              np.allclose(x[0], expect, atol=5e-3), str(x))


def test_grn_three_genes():
    """R3: a 3-gene ring network (mutual repression) has multiple
    attractors (single-gene-high states)."""
    cfg = GRNConfig(
        n_genes=3, alpha=(0.4, 0.4, 0.4), lam=(1.9, 1.9, 1.9),
        rho=(0.1, 0.1, 0.1), s_self=(0.5, 0.5, 0.5), k_deg=(1, 1, 1),
        b_rep=((0, 1, 1), (1, 0, 1), (1, 1, 0)),
        sigma=(0, 0, 0))
    grn = GRN(cfg)
    finals = []
    for seed in [(3, 0.1, 0.1), (0.1, 3, 0.1), (0.1, 0.1, 3)]:
        x = np.array([seed], float)
        u = np.full((1, 3), 0.5)
        e = np.zeros((1, 3))
        for _ in range(6000):
            x = grn.step(x, u, e, 0.05)
        finals.append(tuple(np.round(x[0], 1)))
    check("grn: 3-gene ring -> 3 distinct committed attractors",
          len(set(finals)) == 3, str(finals))


def test_epi_kernels():
    """R1: moments of every inheritance kernel."""
    rng = np.random.default_rng(0)
    v = np.full((200000, 2), 0.5)
    cfg = EpiConfig(kernel="beta")
    d1, d2 = epi.inherit(v, cfg, rng)
    phi = epi.phi_map(0.5, cfg)
    var = phi * (1 - phi) / (1 + cfg.m0)
    check("epi: beta kernel mean", abs(d1[:, 0].mean() - phi) < 2e-3,
          f"{d1[:, 0].mean():.4f} vs {phi:.4f}")
    check("epi: beta kernel variance", abs(d1[:, 0].var() - var) < 2e-4,
          f"{d1[:, 0].var():.5f} vs {var:.5f}")

    cfg = EpiConfig(kernel="gaussian")
    d1, _ = epi.inherit(v, cfg, rng)
    check("epi: gaussian kernel mean", abs(d1[:, 0].mean() - phi) < 2e-3)
    check("epi: gaussian kernel sd",
          abs(d1[:, 0].std() - cfg.sigma_g) < 2e-3)
    check("epi: gaussian kernel support", (d1 >= 0).all() and (d1 <= 1).all())
    # the mean must be preserved even when the truncation is heavy
    # (sigma_g comparable to the distance of phi from the boundary);
    # the sd is then necessarily smaller than sigma_g
    cfg = EpiConfig(kernel="gaussian", sigma_g=0.30)
    for vv in (0.02, 0.5, 0.98):
        ph = epi.phi_map(vv, cfg)
        d1, _ = epi.inherit(np.full((100000, 2), vv), cfg, rng)
        check(f"epi: gaussian mean preserved (v={vv}, sg=0.3)",
              abs(d1[:, 0].mean() - ph) < 5e-3,
              f"{d1[:, 0].mean():.4f} vs {ph:.4f}")
        check(f"epi: gaussian support (v={vv}, sg=0.3)",
              (d1 >= 0).all() and (d1 <= 1).all())

    # cpg kernel, moment-matched mode (default): mean AND variance agree
    # with the beta kernel, so kernel comparisons are controlled
    cfg = EpiConfig(kernel="cpg", cpg_mode="moment", cpg_sites=61)
    d1, _ = epi.inherit(v, cfg, rng)
    var = phi * (1 - phi) / (1 + cfg.m0)
    check("epi: cpg(moment) mean == phi", abs(d1[:, 0].mean() - phi) < 2e-3,
          f"{d1[:, 0].mean():.4f} vs {phi:.4f}")
    check("epi: cpg(moment) variance == beta variance",
          abs(d1[:, 0].var() - var) < 5e-4,
          f"{d1[:, 0].var():.5f} vs {var:.5f}")
    # cpg kernel, mechanistic mode: classical maintain/de-novo rule, whose
    # mean map is affine and therefore has a DIFFERENT fixed point
    cfg = EpiConfig(kernel="cpg", cpg_mode="mechanistic")
    d1, _ = epi.inherit(v, cfg, rng)
    mu = 0.5 * cfg.cpg_p_maintain + 0.5 * cfg.cpg_p_denovo
    sd = np.sqrt(mu * (1 - mu) / cfg.cpg_sites)
    check("epi: cpg(mechanistic) mean == affine map",
          abs(d1[:, 0].mean() - mu) < 2e-3,
          f"{d1[:, 0].mean():.4f} vs {mu:.4f}")
    check("epi: cpg(mechanistic) sd (site noise)",
          abs(d1[:, 0].std() - sd) < 2e-3)
    check("epi: cpg modes give different mean maps (documented)",
          abs(0.5 * cfg.cpg_p_maintain + 0.5 * cfg.cpg_p_denovo - phi) > 0.02,
          "mechanistic mean map is affine, not the Hill map phi")

    cfg = EpiConfig(kernel="deterministic")
    d1, d2 = epi.inherit(v, cfg, rng)
    check("epi: deterministic kernel", np.allclose(d1, phi) and
          np.allclose(d2, phi))


def test_epi_continuous():
    """R1/R5: continuous dynamics fixed point increases with gene
    activity (GRN -> epigenetics coupling)."""
    cfg = EpiConfig(continuous=True, w0=0.01, w1=0.06, e0=0.01, e1=0.06)
    u_lo = np.full((1, 2), 0.5)
    u_hi = np.full((1, 2), 0.5)
    for _ in range(200000):   # long-time limit
        u_lo = epi.continuous_step(u_lo, np.zeros((1, 2)), cfg, 0.05)
        u_hi = epi.continuous_step(u_hi, np.ones((1, 2)), cfg, 0.05)
    check("epi: active genes -> higher mark level",
          u_hi[0, 0] > u_lo[0, 0] + 0.3,
          f"u_hi={u_hi[0, 0]:.3f} u_lo={u_lo[0, 0]:.3f}")
    # analytic fixed point
    act = 1.0
    w = cfg.w0 + cfg.w1 * act
    e = cfg.e0 + cfg.e1 * (1 - act)
    u_star = w / (w + e)
    check("epi: fixed point matches analytic value",
          abs(u_hi[0, 0] - u_star) < 1e-2, f"{u_hi[0, 0]:.4f} vs {u_star:.4f}")


def test_division_partition():
    """R2: binomial partitioning conserves molecules, exact is identity."""
    rng = np.random.default_rng(0)
    cfg = DivisionConfig(partition="binomial", volume=200.0)
    x = np.full((20000, 2), 1.0)
    x1, x2 = div_mod.partition(x, cfg, rng)
    n_tot = np.round(200 * x)
    check("div: binomial conserves molecules",
          np.allclose(np.round(cfg.volume * 0.5 * x1) +
                      np.round(cfg.volume * 0.5 * x2), n_tot, atol=1))
    check("div: binomial mean preserved",
          abs(x1[:, 0].mean() - 1.0) < 0.01)
    # expected variance of daughter concentration: x / (V f (1-f)) ... for
    # f=1/2: Var = x/V / (f*(1-f)*4) = x/V -> sd = sqrt(x/V)
    sd_th = np.sqrt(1.0 / 200.0)
    check("div: binomial sd ~ sqrt(x/V)",
          abs(x1[:, 0].std() - sd_th) < 0.01,
          f"{x1[:, 0].std():.4f} vs {sd_th:.4f}")
    cfg = DivisionConfig(partition="exact")
    x1, x2 = div_mod.partition(x, cfg, rng)
    check("div: exact partition identity",
          np.allclose(x1, x) and np.allclose(x2, x))


def test_rates_continuous_limit():
    """R4: continuous rates converge to the legacy Table-1 values."""
    grn = GRN(GRNConfig(sigma=(0, 0)))
    x_sc = np.array([[1.0, 1.0]])
    x_ta = np.array([[2.0, 0.1]])
    rc = RatesConfig(rate_mode="continuous", p_switch=40.0, K_switch=0.5)
    Q = 1.0     # one stem-like cell
    b, k, m, tau, sen = rates_mod.rates(grn, x_sc, None, rc, np.zeros(1), Q)
    check("rates: SC beta ~ beta0*theta/(theta+Q)",
          abs(b[0] - 0.04 * 200 / 201) < 5e-4, f"{b[0]:.5f}")
    check("rates: SC kappa ~ kappa0", abs(k[0] - 0.01) < 5e-4)
    check("rates: SC tau ~ tau0", abs(tau[0] - 4.0) < 0.05)
    b, k, m, tau, sen = rates_mod.rates(grn, x_ta, None, rc, np.zeros(1), Q)
    check("rates: TA beta ~ beta0", abs(b[0] - 0.04) < 5e-4, f"{b[0]:.5f}")
    check("rates: TA kappa ~ 2 kappa0", abs(k[0] - 0.02) < 5e-4)
    check("rates: TA tau ~ tau0/2", abs(tau[0] - 2.0) < 0.05)
    # senescence counter
    rc2 = RatesConfig(rate_mode="continuous", p_switch=40.0)
    _, _, _, _, sen = rates_mod.rates(grn, x_ta, None, rc2,
                                       np.array([14.9]), 1.0)
    check("rates: senescence when counter >= 15", bool(sen[0]))


def test_simulator_legacy_vs_hsl2024():
    """Cross-validation: CERPS in legacy mode vs the validated hsl2024
    reproduction (statistical agreement of homeostasis counts)."""
    from hsl2024 import DEFAULT, simulate as sim_hsl
    res = sim_hsl(DEFAULT, n_runs=8, n_steps=8000, seed=101)
    hsl_mean = res["sys_ts"][4000:].mean(axis=(0, 1))   # last 1000 h
    cfg = CERPSConfig(max_time=2000.0, seed=202)
    sim = Simulator(cfg)
    sim.run()
    rec = sim.records
    cerps_mean = np.array([
        np.mean(rec["total"][-250:]), np.mean(rec["sc"][-250:]),
        np.mean(rec["ta1"][-250:]), np.mean(rec["ta2"][-250:])])
    check("sim: CERPS legacy total ~ hsl2024",
          abs(cerps_mean[0] - hsl_mean[0]) / hsl_mean[0] < 0.15,
          f"CERPS {cerps_mean[0]:.0f} vs hsl {hsl_mean[0]:.0f}")
    check("sim: CERPS legacy SC ~ hsl2024",
          abs(cerps_mean[1] - hsl_mean[1]) / max(hsl_mean[1], 1) < 0.35,
          f"CERPS {cerps_mean[1]:.0f} vs hsl {hsl_mean[1]:.0f}")
    check("sim: CERPS legacy TA1 ~ hsl2024",
          abs(cerps_mean[2] - hsl_mean[2]) / hsl_mean[2] < 0.25,
          f"CERPS {cerps_mean[2]:.0f} vs hsl {hsl_mean[2]:.0f}")
    check("sim: CERPS legacy TA2 ~ hsl2024",
          abs(cerps_mean[3] - hsl_mean[3]) / hsl_mean[3] < 0.25,
          f"CERPS {cerps_mean[3]:.0f} vs hsl {hsl_mean[3]:.0f}")


def test_simulator_continuous():
    """R4: continuous-rate mode reaches a homeostasis with the same
    qualitative composition (SC + two TA types coexist)."""
    cfg = CERPSConfig(
        rates=RatesConfig(rate_mode="continuous", p_switch=12.0),
        max_time=1500.0, seed=7)
    sim = Simulator(cfg)
    sim.run()
    rec = sim.records
    sc = np.mean(rec["sc"][-150:])
    ta = np.mean(rec["ta1"][-150:]) + np.mean(rec["ta2"][-150:])
    check("sim: continuous mode coexistence (SC>20, TA>500)",
          sc > 20 and ta > 500, f"SC={sc:.0f} TA={ta:.0f}")


def test_simulator_spatial():
    """R6: 3D spatial growth with contact inhibition stays bounded and
    inside the domain; without interactions it matches the well-mixed
    run statistically."""
    cfg = CERPSConfig(
        rates=RatesConfig(rate_mode="continuous", p_switch=12.0),
        space=SpaceConfig(enabled=True, cell_radius=1.0,
                          domain_radius=40.0, contact_inhibition=True),
        n_init=40, max_time=800.0, seed=11)
    sim = Simulator(cfg)
    sim.run()
    rec = sim.records
    pos = sim.space.pos
    r = np.linalg.norm(pos, axis=1)
    check("space: cells inside the domain",
        len(r) == 0 or r.max() <= 40.0 + 1e-6, f"max r={r.max():.2f}")
    check("space: population bounded by contact inhibition",
          sim.max_cells < 4000, f"max cells {sim.max_cells}")
    check("space: population grew and survived",
          len(pos) > 40, f"n={len(pos)}")


def test_coupling_loop():
    """R5: full bidirectional coupling demonstration:
    epigenetics <-> GRN <-> population (niche feedback)."""
    # 1) epi -> GRN: higher u increases expression
    grn = GRN(GRNConfig(sigma=(0, 0)))
    u_lo = np.full((1, 2), 0.15)
    u_hi = np.full((1, 2), 0.55)
    x = np.array([[1.0, 1.0]])
    e = np.zeros((1, 2))
    xl, xh = x.copy(), x.copy()
    for _ in range(50):
        xl = grn.step(xl, u_lo, e, 0.25)
        xh = grn.step(xh, u_hi, e, 0.25)
    check("coupling: epi -> GRN (higher u -> higher expression)",
          xh[0].min() > xl[0].min(),
          f"x(u=0.55)={xh[0]}, x(u=0.15)={xl[0]}")
    # 2) GRN -> epi: activity raises the epigenetic fixed point
    cfg = EpiConfig(continuous=True)
    u = np.full((1, 2), 0.3)
    for _ in range(400):
        u = epi.continuous_step(u, np.ones((1, 2)), cfg, 0.25)
    check("coupling: GRN -> epi (active gene -> u rises)",
          u[0, 0] > 0.3, f"u={u[0, 0]:.3f}")
    # 3) population -> cell: niche feedback reduces beta when Q grows
    rc = RatesConfig(rate_mode="continuous", p_switch=40.0)
    grn2 = GRN(GRNConfig(sigma=(0, 0)))
    x_sc = np.array([[1.0, 1.0]])
    b1, _, _, _, _ = rates_mod.rates(grn2, x_sc, None, rc, np.zeros(1),
                                      Q_soft=10.0)
    b2, _, _, _, _ = rates_mod.rates(grn2, x_sc, None, rc, np.zeros(1),
                                      Q_soft=400.0)
    check("coupling: population -> cell (niche feedback)",
          b2[0] < b1[0], f"beta(Q=10)={b1[0]:.4f} beta(Q=400)={b2[0]:.4f}")


def test_record_decomposition_and_memory_stability():
    """Regression guard for two defects found 2026-10-20.

    D1 (silent cell loss): ``classify_legacy`` emits labels 0..G (0 = all
        high, j+1 = gene j high), but ``Simulator._record`` hard-coded the
        columns {0, 1, 3} plus an explicit count of label 2.  For G >= 4
        the states above label 3 therefore vanished from every ``ho_*``
        aggregate -- measured 20.3% of cells at G=4 and 30.8% at G=5.
        The invariant sc + ta1 + ta2 + other == total must hold for any G.
    D2 (integrator blow-up): the P6 leaky integrator used gain dt/tau
        unclamped, so tau < dt made ``mem`` diverge (measured |mem| ~ 1e174
        at tau=0.1, dt=0.25).  ``mem`` must stay finite and inside [0, 1]
        because its driver is a stemness fraction.
    """
    from cerps.batch import run_one
    from cerps.config import legacy_config, replace

    base = replace(legacy_config(), max_time=300.0, rec_tran_from=0.0,
                   record_states=False)

    # --- D1: decomposition invariant across gene numbers -------------
    for G in (2, 3, 4, 5):
        cfg = replace(base, grn=replace(base.grn, n_genes=G), seed=1)
        r = run_one(cfg, seed=1, keep_sim=True, warn=False)
        rec = r["sim"].records
        tot = np.asarray(rec["total"], float)
        parts = sum(np.asarray(rec[k], float)
                    for k in ("sc", "ta1", "ta2", "other"))
        gap = np.abs(tot - parts).max()
        check(f"record decomposition G={G} (sc+ta1+ta2+other == total)",
              gap < 0.5, f"max|total-parts|={gap:.1f} at n={tot[-1]:.0f}")

    # --- D2: memory bounded and finite even below the step size ------
    for tau in (10.0, 1.0, 0.25, 0.1, 0.01):
        cfg = replace(base, max_time=200.0, seed=1,
                      memory=replace(base.memory, enabled=True, tau=tau,
                                     mem_couple=0.0))
        r = run_one(cfg, seed=1, keep_sim=True, warn=False)
        m = np.asarray(r["sim"].mem, float)
        check(f"memory bounded tau={tau} (dt/tau={base.dt / tau:.2f})",
              np.isfinite(m).all() and m.min() >= 0.0 and m.max() <= 1.0,
              f"mem in [{m.min():.3g}, {m.max():.3g}]")


def test_memory_driver_dispatch_and_run_reset():
    """Regression guard for the P6 driver option and run() re-entry.

    D3 (2026-10-20): ``MemoryConfig.driver`` was a *dead* option -- declared
        in the config but never read, so any non-default value was silently
        ignored.  It is now dispatched through
        ``rates.memory_driver_signal``: ``"stemness"`` (the validated
        baseline) and ``"u_mean"`` are implemented, anything else RAISES
        rather than falling back silently.
    D4: ``run()`` must reset the record lists as well as the P8 delay
        buffer, so a continuation call cannot concatenate two time origins
        into one series.
    """
    from cerps.batch import run_one, validate
    from cerps.config import legacy_config, replace
    from cerps.simulator import Simulator

    base = replace(legacy_config(), max_time=200.0, record_states=False)

    # --- implemented drivers are bounded and genuinely distinct ---------
    mems = {}
    for drv in ("stemness", "u_mean"):
        cfg = replace(base, memory=replace(base.memory, enabled=True,
                                           driver=drv, mem_couple=0.0))
        r = run_one(cfg, seed=1, keep_sim=True, warn=False)
        m = np.asarray(r["sim"].mem, float)
        mems[drv] = m
        check(f"memory driver '{drv}' stays in [0,1]",
              np.isfinite(m).all() and m.min() >= 0.0 and m.max() <= 1.0,
              f"mem in [{m.min():.3g}, {m.max():.3g}]")
    check("memory drivers 'stemness' and 'u_mean' differ",
          abs(mems["stemness"].mean() - mems["u_mean"].mean()) > 1e-3,
          f"means {mems['stemness'].mean():.4f} vs {mems['u_mean'].mean():.4f}")

    # --- an unimplemented driver fails loudly, and validate() says so ---
    cfg = replace(base, memory=replace(base.memory, enabled=True,
                                       driver="bogus"))
    raised = False
    try:
        run_one(cfg, seed=1, warn=False)
    except ValueError:
        raised = True
    check("unimplemented memory driver raises (no silent fallback)", raised)
    check("validate() pre-warns about an unimplemented driver",
          any("driver" in w for w in validate(cfg)))

    # --- run() resets the record series --------------------------------
    cfg = replace(base, rates=replace(base.rates, delay_steps=20))
    sim = Simulator(cfg)
    sim.run(t_max=100.0)
    n_first = len(sim.records["t"])
    sim.run(t_max=200.0)
    n_second = len(sim.records["t"])
    check("run() resets records (no concatenation across calls)",
          n_first == 101 and n_second == 201,
          f"len after run(100)={n_first}, after run(200)={n_second}")


def main():
    print("== CERPS test suite ==")
    for fn in [test_grn_legacy_equivalence, test_grn_steady_states,
               test_grn_three_genes, test_epi_kernels,
               test_epi_continuous, test_division_partition,
               test_rates_continuous_limit, test_simulator_legacy_vs_hsl2024,
               test_simulator_continuous, test_simulator_spatial,
               test_coupling_loop,
               test_record_decomposition_and_memory_stability,
               test_memory_driver_dispatch_and_run_reset]:
        print(f"- {fn.__name__}")
        fn()
    print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
    if FAIL:
        print("FAILED:", FAIL)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
