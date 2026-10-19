"""Gene-expression partitioning at mitosis.

In the paper the two daughters inherit the mother's expression state
exactly ("symmetric division").  CERPS supports:

exact     - daughters inherit x exactly (paper assumption)

binomial  - molecule-level stochastic partitioning.  The strength of the
            partitioning noise is the continuous parameter ``alpha`` of
            Soltani et al. 2016 (PLoS Comput Biol 12:e1004972, Eq. 5):

                E[x_daughter | x]      = x
                Var[x_daughter | x]    = alpha * x * (1-f) / f

            with f the volume fraction of daughter 1 (``asymmetry``).
            The three canonical regimes are

                alpha = 0   deterministic (the discrete analogue of exact)
                alpha = 1   independent binomial partitioning of every
                            molecule (Berg 1978; Huh & Paulsson 2011)
                alpha > 1   over-dispersed: clustered / vesicular
                            segregation (Huh & Paulsson 2011 PNAS)
                0<alpha<1   under-dispersed: active ordered partitioning

            In molecule units (n = V*x molecules, f*n to daughter 1) this
            is Var = alpha*n*f*(1-f); alpha = 1 is implemented as an exact
            binomial draw, alpha > 1 as a beta-binomial with intraclass
            correlation rho = (alpha-1)/(n-1) (which reduces to the
            binomial at alpha = 1), and 0 < alpha < 1 as a variance-
            matched normal approximation.
"""
import numpy as np


def partition(mother_x, cfg, rng):
    """Return (x_daughter1, x_daughter2) for the given mother states."""
    x = mother_x
    if cfg.partition == "exact":
        return x.copy(), x.copy()
    if cfg.partition == "binomial":
        f = cfg.asymmetry
        alpha = getattr(cfg, "alpha", 1.0)
        n = np.maximum(np.round(cfg.volume * x), 0.0).astype(np.int64)
        n1 = _split(n, f, alpha, rng)
        n2 = n - n1
        x1 = n1 / max(f * cfg.volume, 1e-12)
        x2 = n2 / max((1.0 - f) * cfg.volume, 1e-12)
        return x1, x2
    raise ValueError(f"unknown partition mode {cfg.partition!r}")


def _split(n, f, alpha, rng):
    """Molecule counts going to daughter 1.

    alpha = 1 : exact binomial (Var = n f (1-f))
    alpha > 1 : beta-binomial, rho = (alpha-1)/(n-1)  -> same mean, larger
                variance; reduces to the binomial as alpha -> 1
    alpha < 1 : variance-matched normal approximation, clipped to [0, n]
    alpha = 0 : deterministic round(f*n)
    """
    if alpha == 1.0:
        return rng.binomial(n, f)
    if alpha <= 0.0:
        return np.rint(f * n).astype(np.int64)
    if alpha > 1.0:
        out = np.empty_like(n)
        big = n > 1
        if big.any():
            nb = n[big]
            rho = np.clip((alpha - 1.0) / np.maximum(nb - 1, 1), 1e-9,
                          1.0 - 1e-9)
            a = f / rho - f
            b = (1.0 - f) / rho - (1.0 - f)
            p = rng.beta(a, b)
            out[big] = rng.binomial(nb, p)
        small = ~big
        if small.any():
            out[small] = rng.binomial(n[small], f)
        return out
    # 0 < alpha < 1 : normal approximation with the target variance
    sd = np.sqrt(alpha * n * f * (1.0 - f))
    z = rng.standard_normal(n.shape)
    n1 = np.rint(f * n + sd * z)
    return np.clip(n1, 0, n).astype(np.int64)
