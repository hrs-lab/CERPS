"""Cell kinetics: proliferation, removal, apoptosis, cycle duration.

Two modes:

``legacy``      - the discrete Table-1 rules of Huang et al. (2024),
                  with phenotype classification x > 0.5 / x < 0.5.

``continuous``  - all rates are smooth functions of the continuous cell
                  state (no discrete classification).  The core state
                  variable is the *stemness*

                      S = prod_i H+(x_i)      (all-genes-on measure)

                  (configurable per gene).  With a steep Hill switch
                  H_S(v) = v^p/(K^p + v^p), p -> inf, the continuous rules
                  converge to the legacy rules:

                  beta   = beta0 * [1 - H_S(S) * (1 - theta0/(theta0 + Q))]
                           with Q = sum_cells H_S(S) (soft SC count);
                           => SC: beta0*theta/(theta+Q), TA: beta0
                  kappa  = kappa0 * (1 + H_{1-K}(1 - S))
                           => SC: kappa0, TA: 2 kappa0
                  tau    = tau0 * (0.5 + 0.5 * H_S(S))
                           => SC: tau0, TA: tau0/2
                  mu     = mu0
                  senescence: division counter increments by (1 - H_S(S))
                           at each mitosis and decays by factor
                           (1 - H_S(S)) while resting; death at mitosis
                           when the counter >= max_div_ta
                           => TA cells die at the 15th division, SC cells
                              never accumulate (and dedifferentiation
                              resets the counter), as in the C++ code.

Here H_S(v) = v^p/(K^p + v^p) is the steep Hill switch and K is the
stemness threshold (``K_switch``; ``None`` derives it so that the
continuous boundary coincides with the legacy ``x_threshold`` on the
diagonal -- see ``resolve_K``).  The TA-side switches use the
COMPLEMENTARY threshold 1-K, so that both sides fire at the same S = K;
using K on both sides (as an earlier version did) is correct only when
K = 0.5 and leaves stem-like cells with the TA removal rate otherwise.

Contact inhibition (spatial mode) multiplies beta by
theta_rho/(theta_rho + max(0, rho_local - rho0)).
"""
import numpy as np


def hill(v, K, p):
    return v ** p / (K ** p + v ** p)


def resolve_K(cfg_rates, grn):
    """Effective Hill threshold of the stemness switch.

    ``K_switch=None`` (the default) derives the threshold that puts the
    continuous classification boundary on the *legacy* one.  Legacy
    declares a cell SC when every gene exceeds ``x_threshold``, i.e. when
    the stemness S = prod_j H+(x_j) reaches prod_j H+(x_threshold); using
    that value as K makes the two rules agree on the diagonal
    x_1 = ... = x_G = x_threshold.

    The two rules can still disagree off the diagonal: the legacy rule is
    a conjunction of per-gene thresholds (a box in x-space), whereas the
    product rule is a curved region.  Measured over the cube [0,3]^G the
    disagreement covers 5.9% (G=2), 12.6% (G=3) and 19.4% (G=4) -- unavoidable
    for any single-threshold smoothing of a box, and one reason to run the
    continuous mode on the *engine* side (rate_mode="legacy") whenever a
    Table-3-comparable number is required.
    """
    K = getattr(cfg_rates, "K_switch", None)
    if K is not None:
        return float(K)
    G = grn.G
    return float(grn.activation_levels(
        np.full((1, G), cfg_rates.x_threshold)).prod())


def stemness(grn, x):
    """S = product of per-gene activation levels (all-genes-on measure)."""
    return grn.activation_levels(x).prod(axis=1)


def memory_driver_signal(driver, grn, x, u):
    """Per-cell P6 memory driver signal in [0, 1].

    ``driver`` selects what the leaky integrator tracks
    (``MemoryConfig.driver``):

    * ``"stemness"`` -- the stemness fraction S (validated baseline);
    * ``"u_mean"``   -- the mean epigenetic state across genes.

    Both are naturally bounded to [0, 1], which the integrator relies on.
    An unknown name **raises** rather than silently falling back: a silently
    ignored option invalidates a mechanism comparison without any error
    (audit 2026-10-20, defect D3).  This is the single point where the P6
    driver choice is injected, mirroring ``effective_niche_signal`` for P8.
    """
    if driver == "stemness":
        return stemness(grn, x)
    if driver == "u_mean":
        return np.asarray(u, float).mean(axis=1)
    raise ValueError(
        f"memory.driver={driver!r} is not implemented; expected "
        f"'stemness' or 'u_mean'")


def classify_legacy(grn, x, prev, thr):
    """Discrete classification with the keep-previous rule (C++ code)."""
    t = prev.copy()
    hi = x > thr
    lo = x < thr
    if grn.G == 2:
        t[hi[:, 0] & hi[:, 1]] = 0
        t[hi[:, 0] & lo[:, 1]] = 1
        t[lo[:, 0] & hi[:, 1]] = 3
    else:
        # generalisation: label = index of the most-expressed gene if
        # exactly one gene is high, 0 if all high, previous otherwise
        n_hi = hi.sum(axis=1)
        t[(n_hi == grn.G)] = 0
        only = (n_hi == 1)
        for j in range(grn.G):
            m = only & hi[:, j]
            t[m] = j + 1
    return t


def rates(grn, x, u, cfg_rates, div_count, Q_soft, rho_local=None,
          mem=None, mem_couple=0.0):
    """Return (beta, kappa, mu, tau, senesce_at_mitosis_bool).

    beta/kappa/mu/tau: arrays (n_cells,).
    div_count: array (n_cells,) current division counters.
    Q_soft: float, sum over cells of H_S(S) (soft stem-cell count).
    rho_local: optional array of local neighbour counts (spatial).
    """
    r = cfg_rates
    n = x.shape[0]
    if r.rate_mode == "legacy":
        # classification-based (uses previous labels passed via u ignored;
        # caller must pass prev labels through `div_count`'s companion --
        # see simulator: legacy mode uses classify_legacy internally)
        raise NotImplementedError("use legacy_rates() instead")
    S = stemness(grn, x)
    K = resolve_K(r, grn)
    HS = hill(S, K, r.p_switch)
    # NOTE: the "not stem-like" switch must use the COMPLEMENTARY
    # threshold 1-K, not K.  With the same K the two switches fire at
    # S = K and S = 1-K respectively, so any K != 0.5 leaves a band in
    # which a cell is neither properly SC-like nor properly TA-like --
    # e.g. at K = 0.125 a stem-like cell with S = 0.5 receives the TA
    # removal rate kappa = 2*kappa0.  The original hard-coded default
    # K_switch = 0.5 made 1-K == K and hid this.
    HnS = hill(1.0 - S, 1.0 - K, r.p_switch)
    beta = r.beta0 * (1.0 - HS * (1.0 - r.theta0 / (r.theta0 + Q_soft)))
    if rho_local is not None:
        inh = np.maximum(rho_local - r.rho0, 0.0)
        beta = beta * r.theta_rho / (r.theta_rho + inh)
    # P6 coupling: history-dependent memory modulates the proliferation rate.
    # mem in [0,1] (leaky integrator of stemness), centred at 0.5 so the
    # default is neutral; mem_couple == 0 leaves beta untouched.
    if mem is not None and mem_couple != 0.0:
        beta = beta * (1.0 + mem_couple * (mem - 0.5))
    kappa = r.kappa0 * (1.0 + HnS)
    mu = np.full(n, r.mu0)
    tau = r.tau0 * (0.5 + 0.5 * HS)
    senesce = (div_count + 1.0) >= r.max_div_ta
    return beta, kappa, mu, tau, senesce


def legacy_rates(labels, Q_sc, cfg_rates, mem=None, mem_couple=0.0):
    """Table-1 kinetics for the discrete mode (labels: 0 SC, 1 TA1, 3 TA2)."""
    r = cfg_rates
    is_sc = labels == 0
    beta = np.where(is_sc, r.beta0 * r.theta0 / (r.theta0 + Q_sc), r.beta0)
    # P6 coupling (see rates() for semantics)
    if mem is not None and mem_couple != 0.0:
        beta = beta * (1.0 + mem_couple * (mem - 0.5))
    kappa = np.where(is_sc, r.kappa0, 2.0 * r.kappa0)
    mu = np.full(labels.shape, r.mu0)
    tau = np.where(is_sc, r.tau0, 0.5 * r.tau0)
    return beta, kappa, mu, tau


def effective_niche_signal(hist, step, delay_steps, current):
    """Return the niche (stem-cell) signal delayed by ``delay_steps`` steps.

    ``hist`` is the per-step history list, with ``hist[s]`` the signal at
    simulation step ``s``.  When ``step < delay_steps`` the delayed value is
    not yet available, so we fall back to ``current`` (no effective delay
    yet).  ``delay_steps == 0`` reproduces the instantaneous (paper)
    behaviour exactly.  This is the single point where the P8 feedback-time
    delay is injected; it is orthogonal to classification, kernels and the
    spatial modules.
    """
    if delay_steps is None or delay_steps <= 0 or step < delay_steps:
        return current
    return hist[step - delay_steps]
