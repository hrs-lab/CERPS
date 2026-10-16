"""Gene-expression partitioning at mitosis.

In the paper the two daughters inherit the mother's expression state
exactly ("symmetric division").  CERPS supports:

exact     - daughters inherit x exactly (paper assumption)
binomial  - molecule-level stochastic partitioning: with n_i = V * x_i
            molecules of species i, daughter 1 receives
            Binomial(n_i, f) molecules, daughter 2 the rest; daughter
            concentrations are n/V_i with the daughter volumes f*V and
            (1-f)*V (Huh & Paulsson 2011; Huh et al. 2013).
"""
import numpy as np


def partition(mother_x, cfg, rng):
    """Return (x_daughter1, x_daughter2) for the given mother states."""
    x = mother_x
    if cfg.partition == "exact":
        return x.copy(), x.copy()
    if cfg.partition == "binomial":
        f = cfg.asymmetry
        n = np.maximum(np.round(cfg.volume * x), 0.0).astype(np.int64)
        n1 = rng.binomial(n, f)
        n2 = n - n1
        x1 = n1 / max(f * cfg.volume, 1e-12)
        x2 = n2 / max((1.0 - f) * cfg.volume, 1e-12)
        return x1, x2
    raise ValueError(f"unknown partition mode {cfg.partition!r}")
