"""Numerical verification of the M2 generalisation layer.

Five blocks, all pure numpy/scipy, fixed seed, ~20 s total:

V1  notation alias layer is a relabelling, not a transformation
V2  ProposalGRN with default gates reproduces the native GRN EXACTLY,
    and the proposal's literal frozen form matches the native GRN with
    lam = 0
V3  the epigenetic gates alpha/beta/k evaluate to the proposal's formulas
    and differ from the native behaviour only where the native has no
    counterpart
V4  the H switch: zero freezes, linear relaxes to the basal level, general
    reproduces cerps.epigenetics.continuous_step
V5  the RACIPE enumeration and the sign(alpha-beta) criterion (M1 / R5)

Run:  python tests/test_m2_generalisation.py
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                "..")))

from cerps import h_general as hg                                  # noqa: E402
from cerps import notation as nt                                   # noqa: E402
from cerps.config import CERPSConfig, GRNConfig                     # noqa: E402
from cerps.grn import GRN                                           # noqa: E402
from cerps.simulator import Simulator                               # noqa: E402
from cerps import epigenetics as epi                                # noqa: E402

OUT = {}
HERE = os.path.dirname(os.path.abspath(__file__))


def hr(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def tv(u, v):
    return 0.5 * float(np.abs(u - v).sum())


# ======================================================================
# V1  notation
# ======================================================================
def v1_notation():
    hr("V1  [M2-A] notation alias layer: relabelling, not transformation")
    rng = np.random.default_rng(20261020)
    n, G = 40, 3
    x = rng.uniform(0, 3, (n, G))          # HSL: expression
    u = rng.uniform(0, 1, (n, G))          # HSL: epigenetic mark

    # --- array level -------------------------------------------------
    ue, ve = nt.to_proposal_arrays(x, u)   # (x,u) -> (u_expr, v_epi)
    x2, u2 = nt.to_hsl_arrays(ue, ve)
    print(f"  array swap round-trip: x identical = {x is x2}, "
          f"u identical = {u is u2}")
    assert x is x2 and u is u2

    # --- state view --------------------------------------------------
    sim = Simulator(CERPSConfig(n_init=7, seed=1))
    vh = nt.state_view(sim, nt.HSL)
    vp = nt.state_view(sim, nt.PROPOSAL)
    print(f"  HSL view   : x is sim.x = {vh.x is sim.x}, "
          f"u is sim.u = {vh.u is sim.u}")
    print(f"  proposal   : u is sim.x = {vp.u is sim.x}, "
          f"v is sim.u = {vp.v is sim.u}")
    assert vh.x is sim.x and vh.u is sim.u
    assert vp.u is sim.x and vp.v is sim.u
    # the two views must name OPPOSITE slots
    assert vp.expression is vh.x and vp.epigenetic is vh.u

    # --- gate values -------------------------------------------------
    v = np.array([0.0, 0.5, 1.0])
    for kind, sign in (("alpha", +1.0), ("beta", -1.0), ("k", -1.0)):
        g = nt.epigenetic_gate(kind, 2.0, 1.5, v)
        ref = 2.0 * np.exp(sign * 1.5 * v)
        err = float(np.abs(g - ref).max())
        print(f"  gate {kind:<6}: max|g - ref| = {err:.2e}")
        assert err == 0.0

    # --- ProposalRegulation -> GRNConfig (frozen, exact) --------------
    reg = nt.ProposalRegulation(
        n_genes=2,
        alpha_star=(1.0, 1.0), beta_star=(0.0, 1.0, 1.0, 0.0),
        k_star=(1.0, 1.0), s=(0.5, 0.5), n_hill=(4.0, 4.0),
        lam=(0.7, 0.7), rho=(0.1, 0.1), epi_mode="frozen", v0=(0.4, 0.4))
    cfg, gaps = reg.to_grn_config()
    print(f"  frozen translation gaps = {gaps}")
    assert gaps == []
    a_ref, B_ref, k_ref = reg.gates()
    print(f"  alpha -> GRNConfig.alpha : "
          f"{np.allclose(cfg.alpha, a_ref)}")
    print(f"  beta  -> b_rep           : "
          f"{np.allclose(np.asarray(cfg.b_rep).reshape(2, 2), B_ref)}")
    print(f"  k     -> k_deg           : {np.allclose(cfg.k_deg, k_ref)}")
    print(f"  lam forced to 0          : {np.allclose(cfg.lam, 0.0)}")
    assert np.allclose(cfg.alpha, a_ref)
    assert np.allclose(np.asarray(cfg.b_rep).reshape(2, 2), B_ref)
    assert np.allclose(cfg.k_deg, k_ref)
    assert np.allclose(cfg.lam, 0.0)

    # --- beta gate must depend on v0 (regression for the v == 1 bug) --
    # Independent hand calculation: beta_ij(v) = beta_ij* exp(-lam_i v_i).
    # An earlier version hard-coded v = 1, which was invisible here because
    # both cfg.b_rep and B_ref came from the same wrong formula.
    B0 = np.array([[0.0, 1.0], [1.0, 0.0]])
    lam_v = np.array([0.7, 0.7])
    v0_ = np.array([0.4, 0.4])
    B_hand = B0 * np.exp(-lam_v[:, None] * v0_[:, None])
    print(f"  beta gate at v0=0.4 vs hand calc: "
          f"{np.allclose(B_ref, B_hand)}")
    assert np.allclose(B_ref, B_hand)
    B_wrong = B0 * np.exp(-lam_v[:, None])          # the old v == 1 value
    assert not np.allclose(B_ref, B_wrong)
    print(f"  beta(v0=0.4) differs from the old v==1 value: "
          f"{not np.allclose(B_ref, B_wrong)}  (regression closed)")

    # --- dynamic translation reports the two gaps ---------------------
    reg_d = nt.ProposalRegulation(
        n_genes=2, alpha_star=(1.0, 1.0),
        beta_star=(0.0, 1.0, 1.0, 0.0), k_star=(1.0, 1.0),
        s=(0.5, 0.5), n_hill=(4.0,), lam=(0.7, 0.7), rho=(0.1, 0.1),
        epi_mode="dynamic")
    cfg_d, gaps_d = reg_d.to_grn_config()
    print(f"  dynamic translation gaps = {len(gaps_d)}")
    for g in gaps_d:
        print(f"    - {g[:88]}...")
    assert len(gaps_d) == 2
    assert np.allclose(cfg_d.lam, (0.7, 0.7))      # alpha gate kept
    assert np.allclose(cfg_d.b_rep, (0.0, 1.0, 1.0, 0.0))   # gate dropped
    assert np.allclose(cfg_d.k_deg, (1.0, 1.0))             # gate dropped

    # --- mixed-convention guard --------------------------------------
    try:
        nt.check_no_mixed_convention({"u": 1, "u_expr": 2})
        ok = False
    except nt.NotationError as exc:
        ok = True
        print(f"  mixed-convention guard raises: {str(exc)[:60]}...")
    assert ok
    try:
        nt.state_view(sim, "bogus")
        ok = False
    except nt.NotationError:
        ok = True
    assert ok
    print("  => V1 PASS: the alias layer relabels, never transforms; the two"
          " proposal gates with no native counterpart are REPORTED.")

    OUT["V1_notation"] = dict(
        roundtrip_is_identity=True,
        frozen_gaps=gaps, dynamic_gaps=len(gaps_d))


# ======================================================================
# V2  ProposalGRN vs native GRN
# ======================================================================
def v2_exact_reduction():
    hr("V2  [M2-B] ProposalGRN -> native GRN exact reduction")

    cfg = CERPSConfig(seed=3)
    native = GRN(cfg.grn)
    prop = hg.ProposalGRN.from_native(native)
    print(f"  native : alpha={native.alpha}, lam={native.lam}, "
          f"rho={native.rho}, regulation={native.regulation}")
    print(f"  proposal: alpha={prop.alpha_v}, lam={prop.lam_v}, "
          f"rho={prop.rho_v}, regulation={prop.regulation}")

    rng = np.random.default_rng(7)
    n = 500
    x = rng.uniform(0, 3, (n, native.G))
    u = rng.uniform(0, 1, (n, native.G))
    eta = rng.standard_normal((n, native.G))
    dt = 0.25
    xa = native.step(x, u, eta, dt)
    xb = prop.step(x, u, eta, dt)
    d = float(np.abs(xa - xb).max())
    print(f"  max |native.step - proposal.step| = {d:.3e}   (must be 0)")
    assert d == 0.0, d

    # a longer trajectory: identical cell by cell, step by step
    xs = x.copy(); us = u.copy(); es = eta.copy()
    xn_ = x.copy(); un_ = u.copy(); en_ = eta.copy()
    for _ in range(200):
        xs = prop.step(xs, us, es, dt)
        xn_ = native.step(xn_, un_, en_, dt)
    d2 = float(np.abs(xs - xn_).max())
    print(f"  after 200 coupled steps: max|diff| = {d2:.3e}")
    assert d2 == 0.0, d2

    # --- the proposal's LITERAL form: H == 0, gates baked into alpha ---
    # With H == 0 the mark level v is constant, so alpha_i(v) is a constant
    # and the gate is already evaluated.  Reconstruct that case from the
    # native parameters by freezing v at the initial mean and switching the
    # alpha gate OFF (equivalently lam = 0).
    v0 = u.mean(axis=0)
    frozen = hg.ProposalGRN.from_native(native, epi_gate_alpha=False)
    frozen.alpha_v = native.alpha * np.exp(native.lam * v0)
    xa = native.step(x, u, eta, dt)              # dynamic mark level
    xb = frozen.step(x, u, eta, dt)              # frozen at v0
    d3 = float(np.abs(xa - xb).max())
    print(f"  frozen-mark form vs native (same v0 = {v0.round(4)}): "
          f"max|diff| = {d3:.3e}  (>0 expected: the gate is now constant)")
    assert d3 > 0.0

    # and it must agree when the mark level really is frozen at v0
    u_const = np.broadcast_to(v0, u.shape)
    xa = native.step(x, u_const, eta, dt)
    xb = frozen.step(x, u_const, eta, dt)
    d4 = float(np.abs(xa - xb).max())
    print(f"  frozen-mark form vs native with u held at v0: "
          f"max|diff| = {d4:.3e}  (must be 0)")
    assert d4 < 1e-12, d4

    # --- the whole simulator runs with the swapped GRN ---------------
    sim = Simulator(CERPSConfig(n_init=40, max_time=200.0, seed=3))
    hg.attach(sim, prop)
    rec0 = sim.records["total"][:]
    sim.run()
    print(f"  simulator with attached ProposalGRN: {sim.n} cells after "
          f"200 h (native run for comparison below)")
    sim2 = Simulator(CERPSConfig(n_init=40, max_time=200.0, seed=3))
    sim2.run()
    print(f"  simulator with native GRN         : {sim2.n} cells")
    assert sim.n == sim2.n
    print("  => V2 PASS: default gates reproduce the native GRN exactly;"
          " the frozen-mark form is the proposal's H==0 step.")

    OUT["V2_exact_reduction"] = dict(
        max_abs_diff_one_step=d, max_abs_diff_200_steps=d2,
        frozen_vs_dynamic=d3, frozen_vs_native_const_u=d4,
        simulator_cells_attached=int(sim.n),
        simulator_cells_native=int(sim2.n))


# ======================================================================
# V3  the three gates
# ======================================================================
def v3_gates():
    hr("V3  [M2-B] epigenetic gates alpha_i(v), beta_i(v), k_i(v)")

    G = 2
    base = dict(n_genes=G, alpha=(1.0, 1.2), k_deg=(0.8, 1.1),
                s=(0.5, 0.6), n_hill=4.0, lam=(0.7, 0.9),
                rho=(0.0, 0.0), beta=(0.0, 0.9, 1.1, 0.0),
                regulation="additive")
    off = hg.ProposalGRN(**base, epi_gate_alpha=False,
                         epi_gate_beta=False, epi_gate_k=False)
    on = hg.ProposalGRN(**base, epi_gate_alpha=True,
                        epi_gate_beta=True, epi_gate_k=True)

    v = np.array([0.0, 0.25, 0.5, 1.0])
    a_off = off.alpha_of(v)
    a_on = on.alpha_of(v)
    ref = (np.array([1.0, 1.2])[None, :]
           * np.exp(np.array([0.7, 0.9])[None, :] * v[:, None]))
    print(f"  alpha gate off : {np.round(a_off, 6)}  (constants)")
    print(f"  alpha gate on  : max|g - alpha*exp(+lam v)| = "
          f"{np.abs(a_on - ref).max():.2e}")
    assert np.abs(a_on - ref).max() < 1e-14
    assert np.allclose(a_off, np.array([1.0, 1.2]))

    b_on = on.beta_of(v)                       # (4, G, G)
    gate_b = np.exp(-np.array([0.7, 0.9])[None, :] * v[:, None])  # (4, G)
    ref_b = np.array([[0.0, 0.9], [1.1, 0.0]])[None, :, :] * \
        gate_b[:, :, None]                                  # (4, G, G)
    print(f"  beta  gate on  : max|g - beta*exp(-lam v)| = "
          f"{np.abs(b_on - ref_b).max():.2e}")
    assert np.abs(b_on - ref_b).max() < 1e-14
    assert np.allclose(off.beta_of(v), np.array([[0.0, 0.9], [1.1, 0.0]]))

    k_on = on.k_of(v)
    ref_k = (np.array([0.8, 1.1])[None, :]
             * np.exp(-np.array([0.7, 0.9])[None, :] * v[:, None]))
    print(f"  k     gate on  : max|g - k*exp(-lam v)| = "
          f"{np.abs(k_on - ref_k).max():.2e}")
    assert np.abs(k_on - ref_k).max() < 1e-14

    # --- monotonicity in v, as the proposal's signs require ----------
    # (the diagonal of beta is zero, so it is non-increasing rather than
    #  strictly decreasing; the off-diagonal repression edges are strictly
    #  decreasing)
    print(f"  alpha(v) increasing : "
          f"{bool((np.diff(a_on, axis=0) > 0).all())}")
    print(f"  beta(v)  decreasing : "
          f"{bool((np.diff(b_on, axis=0) <= 0).all())}")
    print(f"  k(v)     decreasing : "
          f"{bool((np.diff(k_on, axis=0) < 0).all())}")
    assert (np.diff(a_on, axis=0) > 0).all()
    assert (np.diff(b_on, axis=0) <= 0).all()
    assert (np.diff(k_on, axis=0) < 0).all()

    # --- dynamics actually differ once the gates are on -------------
    rng = np.random.default_rng(11)
    n = 300
    x = rng.uniform(0, 3, (n, G))
    u = rng.uniform(0, 1, (n, G))
    eta = rng.standard_normal((n, G))
    d = float(np.abs(off.step(x, u, eta, 0.25) -
                     on.step(x, u, eta, 0.25)).max())
    print(f"  gates off vs on, one step: max|diff| = {d:.4f}  (>0 required)")
    assert d > 0.0

    # --- the two-gate model has a different steady state -------------
    # integrate both to steady state and compare the expression profiles
    def steady(g, x0, u0, eta0, dt=0.02, steps=40000):
        xs, us = x0.copy(), u0.copy()
        for _ in range(steps):
            xs = g.step(xs, us, eta0, dt)
        return xs
    xs_off = steady(off, x, u, eta)
    xs_on = steady(on, x, u, eta)
    print(f"  steady-state mean x: gates off = {xs_off.mean():.4f}, "
          f"on = {xs_on.mean():.4f}")
    print(f"  TV between the two steady states = {tv(xs_off, xs_on):.4f}")
    assert tv(xs_off, xs_on) > 0.0

    # --- multiplicative + beta gate is rejected at step time --------
    import numpy as _np
    try:
        g_mul = hg.ProposalGRN(n_genes=G, alpha=(1.0, 1.2), k_deg=(0.8, 1.1),
                               s=(0.5, 0.6), n_hill=4.0, lam=(0.7, 0.9),
                               rho=(0.0, 0.0), beta=(0.0, 0.9, 1.1, 0.0),
                               regulation="multiplicative",
                               epi_gate_beta=True)
        g_mul.step(_np.ones((3, G)), _np.ones((3, G)),
                   _np.zeros((3, G)), 0.25)
        ok = False
    except NotImplementedError:
        ok = True
    assert ok
    print("  multiplicative + epi_gate_beta -> NotImplementedError (not "
          "silently ignored)")
    print("  => V3 PASS: the three gates reproduce the proposal's formulas"
          " exactly and are monotone in v as its sign convention requires.")

    OUT["V3_gates"] = dict(
        alpha_gate_max_err=float(np.abs(a_on - ref).max()),
        beta_gate_max_err=float(np.abs(b_on - ref_b).max()),
        k_gate_max_err=float(np.abs(k_on - ref_k).max()),
        one_step_diff=d, steady_state_TV=tv(xs_off, xs_on))


# ======================================================================
# V4  the H switch
# ======================================================================
def v4_h_switch():
    hr("V4  [M2-B] H(u,v) switch: zero = HSL2024 baseline, then linear/general")

    rng = np.random.default_rng(13)
    n, G = 200, 2
    x = rng.uniform(0, 3, (n, G))
    u = rng.uniform(0, 1, (n, G))
    grn = GRN(GRNConfig())
    act = grn.activation_levels(x)
    dt = 0.25
    epi_cfg = CERPSConfig().epi

    # --- zero == HSL2024 continuous_step, bit for bit ----------------
    # This is the semantic correction recorded in the M2 doc: the proposal
    # writes "H == 0", but HSL2024's own epigenetic dynamics is already an
    # activity-dependent relaxation.  "Step 1" therefore means the HSL2024
    # baseline, and mode="zero" must reproduce it exactly.
    hz = hg.ProposalH.from_epi(epi_cfg, mode="zero")
    vz = hz.step(u, act, dt)
    vref = epi.continuous_step(u, act, epi_cfg, dt)
    print(f"  mode='zero'   : max|H.step - continuous_step| = "
          f"{float(np.abs(vz - vref).max()):.2e}  (must be 0)")
    assert float(np.abs(vz - vref).max()) == 0.0
    # the return value is a fresh array, not the input object
    vz[:] = 0.0
    assert float(np.abs(u).sum()) > 0.0

    # --- linear: baseline + h_gain*(act - 1/2), vanishing at act = 1/2
    hl = hg.ProposalH.from_epi(epi_cfg, mode="linear", h_gain=0.5)
    vl = hl.step(u, act, dt)
    d_lin = float(np.abs(vl - vref).max())
    print(f"  mode='linear' : max|H.step - zero| = {d_lin:.3e}  (>0 required)")
    assert d_lin > 0.0
    act_half = np.full_like(act, 0.5)
    d_half = float(np.abs(hl.step(u, act_half, dt)
                          - hz.step(u, act_half, dt)).max())
    print(f"  mode='linear' : max|H.step - zero| at act = 1/2 = {d_half:.2e}"
          f"  (must be 0)")
    assert d_half == 0.0

    # --- general with h_gain = 0 collapses onto zero -----------------
    hg0 = hg.ProposalH.from_epi(epi_cfg, mode="general", h_gain=0.0)
    d_collapse = float(np.abs(hg0.step(u, act, dt) - vref).max())
    print(f"  mode='general', h_gain=0   : max|H.step - zero| = "
          f"{d_collapse:.2e}  (must be 0)")
    assert d_collapse == 0.0

    # --- general with h_gain > 0 departs from zero -------------------
    hg_ = hg.ProposalH.from_epi(epi_cfg, mode="general", h_gain=0.5)
    vg = hg_.step(u, act, dt)
    d_gen = float(np.abs(vg - vref).max())
    print(f"  mode='general', h_gain=0.5: max|H.step - zero| = {d_gen:.3e}"
          f"  (>0 required)")
    assert d_gen > 0.0

    # --- polarity -1 mirrors the baseline about 1/2 ------------------
    # Exact identity for the baseline:  H_-(v, act) == H_+(v, 1 - act).
    # polarity=-1 sets (w,e) = (w0+w1*(1-a), e0+e1*a), exactly what
    # polarity=+1 produces at activity (1-a).  (The general/linear add-ons
    # are proposal-only and are not covered by this identity.)
    hr_ = hg.ProposalH.from_epi(epi_cfg, mode="zero", polarity=-1.0)
    lhs = hr_.step(u, act, dt)
    rhs = hz.step(u, 1.0 - act, dt)
    print(f"  polarity -1 baseline identity: max|H_-(v,act) - H_+(v,1-act)| = "
          f"{float(np.abs(lhs - rhs).max()):.2e}")
    assert float(np.abs(lhs - rhs).max()) < 1e-14

    # --- the ordering is a genuine refinement ------------------------
    # zero -> general (h_gain = 0.5) must CHANGE the mark distribution,
    # otherwise the generalisation would be vacuous.
    vm = {}
    for mode, gain in (("zero", 0.0), ("linear", 0.5), ("general", 0.5)):
        h = hg.ProposalH.from_epi(epi_cfg, mode=mode, h_gain=gain)
        v = u.copy()
        for _ in range(4000):
            v = h.step(v, act, dt)
        vm[mode] = v
        print(f"  mode={mode:<8}, h_gain={gain:<4}: mean v = {v.mean():.4f},"
              f"  std v = {v.std():.4f}")
    tv_zg = tv(vm["zero"], vm["general"])
    tv_zl = tv(vm["zero"], vm["linear"])
    print(f"  TV(zero, general) = {tv_zg:.4f},  TV(zero, linear) = "
          f"{tv_zl:.4f}  (>0 required)")
    assert tv_zg > 0.0 and tv_zl > 0.0

    try:
        hg.ProposalH(mode="bogus")
        ok = False
    except ValueError:
        ok = True
    assert ok
    print("  => V4 PASS: mode='zero' reproduces the HSL2024 baseline exactly,"
          " and the linear/general levels are non-vacuous refinements.")

    OUT["V4_h_switch"] = dict(
        zero_matches_continuous_step=True,
        linear_departure=d_lin, general_collapse_at_zero_gain=d_collapse,
        general_departure=d_gen,
        TV_zero_vs_general=tv_zg, TV_zero_vs_linear=tv_zl)


# ======================================================================
# V5  RACIPE enumeration + net direction  (M1 / R5)
# ======================================================================
def v5_racipe():
    hr("V5  [M2-A/M1 R5] RACIPE enumeration and sign(alpha-beta)")

    states = hg.racipe_states()
    print("  admissible (alpha, beta) states:")
    for name, a, b in states:
        print(f"    {name:<12} alpha={a:<5} beta={b:<5}")
    assert len(states) == 3
    assert {s[0] for s in states} == {"activation", "repression", "none"}

    # the criterion: sign(alpha - beta)
    cases = [
        # (alpha, beta, expected sign, why)
        (2.0, 1.0, +1, "alpha>beta>0 is still an ACTIVATOR, not repression"),
        (0.0, 1.0, -1, "pure repression needs alpha = 0"),
        (1.0, 0.0, +1, "pure activation"),
        (0.0, 0.0, 0, "no regulation"),
        (1.0, 1.0, 0, "alpha == beta is a net-null edge"),
    ]
    print("\n  net direction  sign(alpha - beta):")
    for a, b, exp, why in cases:
        got = int(np.sign(a - b))
        print(f"    alpha={a:<4} beta={b:<4} -> sign = {got:+d}   ({why})")
        assert got == exp

    # --- the exact HSL2024 toggle -----------------------------------
    # paper Eq. (1): dx1/dt = a1(rho1+(1-rho1)H+(x1)) + b1*H-(x2) - k1 x1
    # => alpha_11 = a1, beta_11 = a1*rho1 (self, activator)
    #    alpha_12 = 0,  beta_12 = b1        (cross, pure repression)
    native = GRN(GRNConfig())                     # default 2-gene mutual rep
    prop = hg.ProposalGRN.from_native(native)
    A = prop.alpha_matrix()
    B = prop.beta_matrix()
    D = prop.net_direction()
    print(f"\n  HSL2024 default 2-gene toggle:")
    print(f"    alpha = \n{np.round(A, 4)}")
    print(f"    beta  = \n{np.round(B, 4)}")
    print(f"    sign(alpha-beta) = \n{D}")
    assert D[0, 0] == +1 and D[1, 1] == +1, "self edges must activate"
    assert D[0, 1] == -1 and D[1, 0] == -1, "cross edges must repress"
    print("    self edges -> +1 (activation), cross edges -> -1 (repression)"
          "  [matches the paper]")

    # --- a naive (alpha>0, beta>0) draw is NOT repression -----------
    naive = hg.ProposalGRN(n_genes=2, alpha=(1.0, 1.0),
                           beta=(0.0, 1.0, 1.0, 0.0), rho=(0.0, 0.0))
    Dn = naive.net_direction()
    print(f"\n  a RACIPE draw with alpha=(1,1), beta_offdiag=1:")
    print(f"    sign(alpha-beta) = \n{Dn}")
    assert Dn[0, 1] == -1, "with alpha_12 = 0 the cross edge IS repression"
    # now the failure mode: BOTH non-zero
    both = hg.ProposalGRN(n_genes=2, alpha=(1.0, 1.0),
                          beta=(0.0, 1.0, 1.0, 0.0), rho=(0.0, 0.0),
                          c_act=(0.0, 1.0, 1.0, 0.0))
    Db = both.net_direction()
    print(f"    with c_act=1 as well (alpha_12 = 1, beta_12 = 1): "
          f"sign = {Db[0, 1]:+d}  <-- NOT repression any more")
    assert Db[0, 1] == 0
    print("    => a sample carrying both terms is a net-null edge, which the"
          "\n       three-state enumeration excludes by construction."
          "\n       [This is the M1/R5 correction: the earlier criterion"
          "\n        'beta > alpha*(1-rho)' was WRONG and is withdrawn.]")

    # --- the bistability window from M1 V3b -------------------------
    print("\n  bistability window check (M1 V3b: a/k in [s, ~3s]):")
    ok_cases = [("in window", 1.0, 1.0, 0.5, 2.0),
                ("too weak", 0.2, 1.0, 0.5, 2.0),
                ("too strong", 5.0, 1.0, 0.5, 2.0)]
    for name, a, k, s, n_ in ok_cases:
        g = hg.ProposalGRN(n_genes=2, alpha=(a, a), beta=(0.0, 1.0, 1.0, 0.0),
                           k_deg=(k, k), s=(s, s), n_hill=n_, rho=(0.0, 0.0),
                           lam=(0.0, 0.0), epi_gate_alpha=False)
        f = g.validate()
        print(f"    {name:<12} a/k={a/k:.2f}: {len(f)} finding(s)"
              + (f"  e.g. {f[0][:60]}..." if f else ""))
    assert len(hg.ProposalGRN(n_genes=2, alpha=(1.0, 1.0),
                              beta=(0.0, 1.0, 1.0, 0.0), k_deg=(1.0, 1.0),
                              s=(0.5, 0.5), n_hill=4.0,
                              rho=(0.0, 0.0), lam=(0.0, 0.0),
                              epi_gate_alpha=False).validate()) == 0
    print("  => V5 PASS: the enumeration and the criterion are the M1"
          " outcome, and the HSL2024 toggle is classified correctly.")

    OUT["V5_racipe"] = dict(
        states=[s[0] for s in states],
        hsl_toggle_signs=D.tolist(),
        mixed_edge_sign=int(Db[0, 1]))


# ======================================================================
# V6  end-to-end: the proposal's model actually runs
# ======================================================================
def v6_end_to_end():
    hr("V6  [M2-B] end-to-end: zero-regression and a non-vacuous H switch")

    # the proposal's configuration: gates baked into alpha, H level as switch
    reg = nt.ProposalRegulation(
        n_genes=2, alpha_star=(1.0, 1.0), beta_star=(0.0, 1.0, 1.0, 0.0),
        k_star=(1.0, 1.0), s=(0.5, 0.5), n_hill=(4.0,), lam=(0.7, 0.7),
        rho=(0.0, 0.0), epi_mode="frozen", v0=(0.5, 0.5))
    cfg_grn, gaps = reg.to_grn_config()
    assert gaps == []

    # ---- V6a  zero-regression: mode="zero" == HSL2024 continuous_step ----
    # Opening the proposal layer must NOT perturb an existing native run,
    # otherwise every CERPS/HSL2024 result obtained so far would be
    # invalidated.  Same seed, same config, native continuous_step vs
    # proposal_h(mode="zero"): the epigenetic trajectories must coincide
    # bit for bit.
    cfgA = CERPSConfig(n_init=60, max_time=200.0, seed=5)
    cfgA.grn = cfg_grn
    cfgA.epi.continuous = True
    simA = Simulator(cfgA)
    simA.run()

    cfgB = CERPSConfig(n_init=60, max_time=200.0, seed=5)
    cfgB.grn = cfg_grn
    cfgB.epi.continuous = True
    simB = Simulator(cfgB)
    hg.attach(simB, None, hg.ProposalH.from_epi(cfgB.epi, mode="zero"))
    simB.run()
    assert simA.n == simB.n, (simA.n, simB.n)
    d0 = float(np.abs(simA.u - simB.u).max())
    print(f"  [V6a zero-regression] max|u(native) - u(H=zero)| = {d0:.3e}"
          f"   (must be 0)")
    assert d0 == 0.0, d0
    u_zero = simB.u.copy()

    # ---- V6b  the H switch is not vacuous: general must move the marks ---
    cfg2 = CERPSConfig(n_init=60, max_time=200.0, seed=5)
    cfg2.grn = cfg_grn
    cfg2.epi.continuous = True
    sim2 = Simulator(cfg2)
    hg.attach(sim2, None,
              hg.ProposalH.from_epi(cfg2.epi, mode="general", h_gain=0.5))
    sim2.run()
    # NOTE: n need not match -- the H generalisation feeds back on the GRN
    # gates and therefore on cell fate, which is exactly what makes the
    # switch non-vacuous.  Compare per-gene mark means, not element-wise.
    d = float(np.abs(sim2.u.mean(axis=0) - u_zero.mean(axis=0)).max())
    print(f"  [V6b H switch] max|mean u(general) - mean u(zero)| = {d:.4f}"
          f"   (>0 required);  n = {simB.n} vs {sim2.n}")
    assert d > 1e-3, d

    # ---- V6c  full proposal form (all three gates on) also runs ---------
    cfg3 = CERPSConfig(n_init=60, max_time=200.0, seed=5)
    sim3 = Simulator(cfg3)
    pgrn3 = hg.ProposalGRN.from_native(sim3.grn, regulation="additive",
                                       epi_gate_beta=True, epi_gate_k=True)
    hg.attach(sim3, pgrn3,
              hg.ProposalH.from_epi(cfg3.epi, mode="general", h_gain=0.5))
    sim3.run()
    print(f"  [V6c full proposal] n = {sim3.n:4d}  "
          f"mean u = {np.round(sim3.u.mean(axis=0), 4)}")
    print(f"       max|mean u(full) - mean u(zero)| = "
          f"{float(np.abs(sim3.u.mean(axis=0) - u_zero.mean(axis=0)).max()):.4f}")

    # ---- V6d  3 genes (the CD4 T direction) -----------------------------
    cfg4 = CERPSConfig(n_init=60, max_time=200.0, seed=5)
    cfg4.grn = GRNConfig(n_genes=3, alpha=(1.0, 1.0, 1.0),
                         b_rep=(0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0),
                         rho=(0.0, 0.0, 0.0), k_deg=(1.0, 1.0, 1.0),
                         s_self=(0.5, 0.5, 0.5), lam=(0.7, 0.7, 0.7))
    cfg4.epi.continuous = True
    sim4 = Simulator(cfg4)
    hg.attach(sim4, hg.ProposalGRN.from_native(sim4.grn, regulation="additive",
                                               epi_gate_beta=True,
                                               epi_gate_k=True),
              hg.ProposalH.from_epi(cfg4.epi, mode="general", h_gain=0.5))
    sim4.run()
    print(f"  [V6d 3-gene triad] n = {sim4.n:4d}  "
          f"labels = {dict(zip(*np.unique(sim4.label, return_counts=True)))}")
    print("  => V6 PASS: mode=zero reproduces the native run exactly,"
          " general moves the marks, and 2/3-gene systems run end to end.")

    OUT["V6_end_to_end"] = dict(
        zero_regression_max_abs_diff=d0,
        h_switch_max_abs_diff=d,
        n_zero=int(simB.n), n_general=int(sim2.n),
        n_full_proposal=int(sim3.n), n_three_gene=int(sim4.n),
        u_zero_mean=u_zero.mean(axis=0).tolist(),
        u_general_mean=sim2.u.mean(axis=0).tolist())


# ======================================================================
# V7  the M2b precondition: b_base != 1 with the AND-gate form
# ======================================================================
def v7_bbase_precondition():
    hr("V7  [M2b] from_native refuses b_base != 1 under 'product'")

    # ---- V7a  product + non-unit b_base -> ValueError (no silent drop) --
    cfg = CERPSConfig(seed=3)
    cfg.grn.b_base = (2.0, 2.0)
    native = GRN(cfg.grn)
    assert native.regulation == "product"
    assert native.has_cross
    raised = None
    try:
        hg.ProposalGRN.from_native(native)
    except ValueError as exc:
        raised = str(exc)
    print(f"  [V7a product+b_base=2] ValueError raised: "
          f"{raised is not None}")
    assert raised is not None, "M2b: b_vec was silently dropped"
    assert "b_base" in raised and "additive" in raised, raised
    # demonstrate the deviation the guard prevents
    rng = np.random.default_rng(7)
    x = rng.uniform(0, 3, (200, native.G))
    u = rng.uniform(0, 1, (200, native.G))
    eta = rng.standard_normal((200, native.G))
    xa = native.step(x, u, eta, 0.25)
    prop = hg.ProposalGRN(n_genes=native.G,
                          alpha=tuple(native.alpha.tolist()),
                          beta=tuple(native.B.ravel().tolist()),
                          k_deg=tuple(native.k_deg.tolist()),
                          s=tuple(native.s_self.tolist()),
                          n_hill=native.n,
                          lam=tuple(native.lam.tolist()),
                          rho=tuple(native.rho.tolist()),
                          c_act=tuple(native.C.ravel().tolist()),
                          regulation="multiplicative")
    d = float(np.abs(xa - prop.step(x, u, eta, 0.25)).max())
    # the dropped term is (b_base-1) * prod_j Hn^B, bounded by 1 in the
    # production rate, so the one-step deviation is bounded by dt = 0.25
    print(f"  [V7a] unguarded drop would deviate by max|diff| = {d:.6f}"
          f"  (<= dt*(b_base-1) = 0.25)")
    assert 0.2 < d <= 0.25 + 1e-9, d

    # ---- V7b  additive + non-unit b_base -> still exact ---------------
    # In the additive composition the native dynamics never uses b_base,
    # so nothing is lost by not carrying it over.
    cfg2 = CERPSConfig(seed=3)
    cfg2.grn.b_base = (2.0, 2.0)
    cfg2.grn.regulation = "additive"
    native2 = GRN(cfg2.grn)
    prop2 = hg.ProposalGRN.from_native(native2)
    d2 = float(np.abs(native2.step(x, u, eta, 0.25)
                      - prop2.step(x, u, eta, 0.25)).max())
    print(f"  [V7b additive+b_base=2] max|diff| = {d2:.3e}   (must be 0)")
    assert d2 == 0.0, d2

    # ---- V7c  product + non-unit b_base but NO cross edges ------------
    # b_vec is inert without any B/C edge, so the guard must not fire.
    cfg3 = CERPSConfig(seed=3)
    cfg3.grn.b_base = (2.0, 2.0)
    cfg3.grn.b_rep = (0.0, 0.0, 0.0, 0.0)
    cfg3.grn.c_act = (0.0, 0.0, 0.0, 0.0)
    native3 = GRN(cfg3.grn)
    assert not native3.has_cross
    prop3 = hg.ProposalGRN.from_native(native3)
    d3 = float(np.abs(native3.step(x, u, eta, 0.25)
                      - prop3.step(x, u, eta, 0.25)).max())
    print(f"  [V7c product+b_base=2, no cross] max|diff| = {d3:.3e}"
          f"   (must be 0)")
    assert d3 == 0.0, d3

    print("  => V7 PASS: the M2b precondition is enforced, and the guard"
          " is inert exactly where b_vec is unused.")
    OUT["V7_bbase_precondition"] = dict(
        v7a_raised=True, v7a_prevented_deviation=d,
        v7b_additive_max_abs_diff=d2, v7c_no_cross_max_abs_diff=d3)


# ======================================================================
# pytest entry points.
#
# The V-checks below are written as plain functions (``v1_..v7_``) so that
# ``python tests/test_m2_generalisation.py`` prints a human-readable report.
# Without these wrappers ``pytest tests/`` collected *nothing* from this
# file -- the whole M2 layer (h_general.py + notation.py) was invisible to
# CI.  The wrappers add no logic: each one asserts that the corresponding
# check runs to completion (the functions already ``assert`` internally).
# ======================================================================
def test_v1_notation_alias_layer():
    v1_notation()


def test_v2_exact_reduction_to_native_grn():
    v2_exact_reduction()


def test_v3_epigenetic_gates():
    v3_gates()


def test_v4_h_switch_levels():
    v4_h_switch()


def test_v5_racipe_enumeration_and_sign_criterion():
    v5_racipe()


def test_v6_end_to_end_simulation():
    v6_end_to_end()


def test_v7_b_base_precondition():
    v7_bbase_precondition()


# ======================================================================
def main():
    v1_notation()
    v2_exact_reduction()
    v3_gates()
    v4_h_switch()
    v5_racipe()
    v6_end_to_end()
    v7_bbase_precondition()

    hr("SUMMARY")
    for k, v in OUT.items():
        print(f"  {k}: {json.dumps(v, default=str)[:150]}")
    path = os.path.join(HERE, "results_m2.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(OUT, fh, indent=2, default=str)
    print(f"\n  results written to {path}")
    print("\n  ALL M2 CHECKS PASSED")


if __name__ == "__main__":
    main()
