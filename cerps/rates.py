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
                  kappa  = kappa0 * (1 + H_S(1 - S))
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

Contact inhibition (spatial mode) multiplies beta by
theta_rho/(theta_rho + max(0, rho_local - rho0)).
"""
import numpy as np


def hill(v, K, p):
    return v ** p / (K ** p + v ** p)


def stemness(grn, x):
    """S = product of per-gene activation levels (all-genes-on measure)."""
    return grn.activation_levels(x).prod(axis=1)


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


def rates(grn, x, u, cfg_rates, div_count, Q_soft, rho_local=None):
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
    HS = hill(S, r.K_switch, r.p_switch)
    HnS = hill(1.0 - S, r.K_switch, r.p_switch)
    beta = r.beta0 * (1.0 - HS * (1.0 - r.theta0 / (r.theta0 + Q_soft)))
    if rho_local is not None:
        inh = np.maximum(rho_local - r.rho0, 0.0)
        beta = beta * r.theta_rho / (r.theta_rho + inh)
    kappa = r.kappa0 * (1.0 + HnS)
    mu = np.full(n, r.mu0)
    tau = r.tau0 * (0.5 + 0.5 * HS)
    senesce = (div_count + 1.0) >= r.max_div_ta
    return beta, kappa, mu, tau, senesce


def legacy_rates(labels, Q_sc, cfg_rates):
    """Table-1 kinetics for the discrete mode (labels: 0 SC, 1 TA1, 3 TA2)."""
    r = cfg_rates
    is_sc = labels == 0
    beta = np.where(is_sc, r.beta0 * r.theta0 / (r.theta0 + Q_sc), r.beta0)
    kappa = np.where(is_sc, r.kappa0, 2.0 * r.kappa0)
    mu = np.full(labels.shape, r.mu0)
    tau = np.where(is_sc, r.tau0, 0.5 * r.tau0)
    return beta, kappa, mu, tau
