"""Epigenetic state dynamics: continuous intracellular dynamics and
stochastic inheritance kernels at cell division.

u_i in [0, 1] is the fraction of activating epigenetic marks (methylated
CpG sites or modified nucleosomes) at the regulatory region of gene i.

1. Continuous dynamics (within the cell cycle), GRN -> epigenetics:

       du_i/dt = (w0 + w1 * act_i) * (1 - u_i)
                 - (e0 + e1 * (1 - act_i)) * u_i

   where act_i is the activation level of gene i.  Activity-dependent
   mark writing/erasing: active genes recruit chromatin writers that lay
   down activating marks, inactive genes lose them.  Precedents for the
   two directions of this coupling are Berry, Dean & Howard 2017 (Cell
   Systems 4:445 -- the erasure rate is proportional to the transcription
   rate f, Eq. 3) and Miyamoto, Furusawa & Kaneko 2015 (PLoS Comput Biol
   11:e1004476 -- a slow epigenetic variable relaxes toward a
   transcription-dependent target, tau_epi*dtheta/dt = (Theta - alpha*x)
   - theta).  NOTE: if u denotes a *repressive* mark (H3K27me3-like or
   DNA methylation) rather than an activating one, the sign of the
   coupling must be reversed -- the code does not know which mark u is.

2. Division kernels (epigenetic state of the two daughters given the
   mother's state v).  All kernels share the same conditional mean
   phi(v) so that switching kernels changes ONLY the noise, not the
   deterministic drift -- this is what makes cross-kernel comparisons
   controlled experiments:

   beta          - paper Eqs. (6)-(8): u ~ Beta(m0*phi, m0*(1-phi)),
                   phi(v) = m1 + m2*(m3 v)^m4/(1 + (m3 v)^m4)
   gaussian      - u ~ N(phi(v), sigma_g^2) truncated to [0, 1]
                   (proper truncated normal by rejection)
   cpg           - mechanistic semi-conservative inheritance over M
                   explicit sites.  ``cpg_mode="moment"`` (default) sets
                   the per-site marking probability to phi(v), so that
                   E[u]=phi(v) and Var[u]=phi(1-phi)/M -- with
                   M = m0+1 this matches the beta kernel exactly in both
                   moments.  ``cpg_mode="mechanistic"`` keeps the
                   classical maintain/de-novo rule
                   (p_maintain / p_denovo), whose mean map is the affine
                   function p_m*v + p_d*(1-v) and therefore has a
                   DIFFERENT fixed point (0.143 vs 0.132 for the default
                   parameters): use it only for mechanistic studies, not
                   for controlled comparisons against the other kernels.
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
        mode = getattr(cfg, "cpg_mode", "moment")
        n, G = v.shape
        out = []
        if mode == "moment":
            # per-site marking probability = phi(v): matches the beta
            # kernel in mean AND (with M = m0+1) in variance
            q = np.clip(phi_map(v, cfg), 0.0, 1.0)
            for d in range(2):
                u_new = np.empty_like(v)
                for g in range(G):
                    p = np.repeat(q[:, g:g + 1], M, axis=1)
                    u_new[:, g] = (rng.random((n, M)) < p).mean(axis=1)
                out.append(u_new)
            return out[0], out[1]
        # mechanistic mode: semi-conservative maintenance / de-novo
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


def _trunc_normal(target, sd, rng, a=0.0, b=1.0):
    """Truncated normal on [a, b] whose mean is exactly ``target``.

    Truncating a normal necessarily shifts its mean: for phi = 0.081 and
    sd = 0.30 about 39% of the mass lies below 0, and the conditional
    mean of the truncated law is 0.270, not 0.081.  Neither clipping nor
    folding restores the mean (folding gives 0.248, rejection gives
    0.270 -- both wrong).  The kernel instead *solves* for the
    pre-truncation mean mu* such that E[N(mu*, sd) | a <= X <= b] equals
    the requested target, then draws by rejection.  The result is a true
    truncated normal whose mean matches phi(v) exactly, so the gaussian
    kernel is a controlled noise-only variant of the beta / cpg kernels.
    """
    target = np.asarray(target, float)
    shape = target.shape
    t = target.ravel()
    out = np.empty_like(t)
    if sd <= 0:
        return np.clip(target, a, b)
    # elements whose truncation is negligible need no correction
    far = (t - a > 6.0 * sd) & (b - t > 6.0 * sd)
    mu = t.copy()                           # pre-truncation mean
    need = ~far
    if need.any():
        tn = t[need]
        lo = np.full(tn.shape, a - 12.0 * sd)
        hi = np.full(tn.shape, b + 12.0 * sd)
        # bisection on the monotone map mu -> truncated mean; 1e-7*sd
        # accuracy is far below Monte-Carlo noise, so ~25 iterations
        # suffice and we stop as soon as every element has converged
        m = tn.copy()
        for _ in range(40):
            tm = _trunc_mean(m, sd, a, b)
            too_low = tm < tn
            lo = np.where(too_low, m, lo)
            hi = np.where(too_low, hi, m)
            m = np.where(too_low, 0.5 * (m + hi), 0.5 * (m + lo))
            if np.max(hi - lo) < 1e-7 * max(sd, 1e-12):
                break
        mu[need] = m
    # single rejection pass over all elements
    pending = np.ones(t.shape, bool)
    val = np.empty_like(t)
    for _ in range(200):
        if not pending.any():
            break
        z = rng.normal(mu[pending], sd)
        ok = (z >= a) & (z <= b)
        idx = np.flatnonzero(pending)
        val[idx[ok]] = z[ok]
        pending[idx[ok]] = False
    if pending.any():                    # pathologically small acceptance
        val[pending] = np.clip(t[pending], a, b)
    return val.reshape(shape)


def _trunc_mean(mu, sd, a=0.0, b=1.0):
    """E[N(mu, sd) | a <= X <= b] (vectorised)."""
    from scipy.special import ndtr
    al = (a - mu) / sd
    bl = (b - mu) / sd
    pa = np.exp(-0.5 * al * al) / np.sqrt(2.0 * np.pi)
    pb = np.exp(-0.5 * bl * bl) / np.sqrt(2.0 * np.pi)
    denom = ndtr(bl) - ndtr(al)
    denom = np.where(denom < 1e-300, 1e-300, denom)
    return mu + sd * (pa - pb) / denom
