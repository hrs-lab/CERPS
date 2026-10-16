"""Epigenetic state dynamics: continuous intracellular dynamics and
stochastic inheritance kernels at cell division.

u_i in [0, 1] is the fraction of activating epigenetic marks (methylated
CpG sites or modified nucleosomes) at the regulatory region of gene i.

1. Continuous dynamics (within the cell cycle), GRN -> epigenetics:

       du_i/dt = (w0 + w1 * act_i) * (1 - u_i)
                 - (e0 + e1 * (1 - act_i)) * u_i

   where act_i is the activation level of gene i; active genes recruit
   chromatin writers that lay down activating marks, inactive genes lose
   them (Dodd et al. 2007; Probst et al. 2009; Zhao et al. 2021).

2. Division kernels (epigenetic state of the two daughters given the
   mother's state v), "not limited to the beta distribution":

   beta          - paper Eqs. (6)-(8): u ~ Beta(m0*phi, m0*(1-phi)),
                   phi(v) = m1 + m2*(m3 v)^m4/(1 + (m3 v)^m4)
   gaussian      - u ~ N(phi(v), sigma_g^2) truncated to [0, 1]
   cpg           - mechanistic semi-conservative inheritance: each of
                   M explicit sites is retained with probability
                   p_maintain if marked, gained with probability p_denovo
                   if unmarked (mimicking DNMT1 maintenance and de-novo
                   DNMT3 methylation / passive demethylation)
   deterministic - u = phi(v)
"""
import numpy as np


def phi_map(v, cfg):
    """Conditional mean of the daughter's epigenetic state (Hill map)."""
    g = (cfg.m3 * v) ** cfg.m4
    return cfg.m1 + cfg.m2 * g / (1.0 + g)


def continuous_step(u, act, cfg, dt):
    """One explicit-Euler step of the continuous epigenetic dynamics."""
    w = cfg.w0 + cfg.w1 * act
    e = cfg.e0 + cfg.e1 * (1.0 - act)
    du = w * (1.0 - u) - e * u
    return np.clip(u + du * dt, 0.0, 1.0)


def inherit(mother_u, cfg, rng):
    """Draw daughter epigenetic states.

    mother_u : (n_divisions, G) mother states
    returns (daughter1, daughter2): two (n_divisions, G) arrays
    """
    v = mother_u
    if cfg.kernel == "deterministic":
        d1 = phi_map(v, cfg)
        d2 = d1.copy()
        return d1, d2
    if cfg.kernel == "beta":
        phi = phi_map(v, cfg)
        a = cfg.m0 * phi
        b = cfg.m0 * (1.0 - phi)
        d1 = rng.beta(a, b)
        d2 = rng.beta(a, b)
        return d1, d2
    if cfg.kernel == "gaussian":
        phi = phi_map(v, cfg)
        d1 = _trunc_normal(phi, cfg.sigma_g, rng)
        d2 = _trunc_normal(phi, cfg.sigma_g, rng)
        return d1, d2
    if cfg.kernel == "cpg":
        M = cfg.cpg_sites
        # explicit site-level inheritance: site states ~ Bernoulli(u)
        # in the mother; maintenance / de-novo probabilities per site.
        n, G = v.shape
        out = []
        for d in range(2):
            u_new = np.empty_like(v)
            for g in range(G):
                marked = rng.random((n, M)) < v[:, g:g + 1]
                keep = rng.random((n, M)) < cfg.cpg_p_maintain
                gain = rng.random((n, M)) < cfg.cpg_p_denovo
                daughter_marked = np.where(marked, keep, gain)
                u_new[:, g] = daughter_marked.mean(axis=1)
            out.append(u_new)
        return out[0], out[1]
    raise ValueError(f"unknown kernel {cfg.kernel!r}")


def _trunc_normal(mu, sd, rng):
    """Truncated-to-[0,1] normal draw (reflective clipping)."""
    z = rng.normal(mu, sd)
    # fold back into [0, 1]
    z = np.abs(z) % 2.0
    return np.where(z > 1.0, 2.0 - z, z)
