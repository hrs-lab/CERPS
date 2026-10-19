"""N-gene regulation network dynamics (vectorised over cells).

Generalises Eq. (1) of Huang et al. (2024) to G genes.  Two cross-term
compositions are available via ``GRNConfig.regulation``:

``product`` (default)  -- b_i * prod_j H-(x_j)^{B_ij}: multiplicative /
    AND-gate combination.  For G = 2 with B_ij in {0,1} this reduces
    exactly to Eq. (1) of the paper, and it coincides with the additive
    form in that case.  Direct precedent: Duddu et al. 2026 (J. R. Soc.
    Interface 23:20251065) Eq. 4.5; thermodynamic tradition: Buchler,
    Gerland & Hwa 2003 (PNAS 100:5136).

``additive`` -- sum_j B_ij * H-(x_j): the additive form of the paper's
    text (the paper notes that its results are consistent with both).

For gene i:

    dx_i/dt = a_i(u_i, eta_i) * (rho_i + (1-rho_i) * H+(x_i))
              + sum_{j != i} B_ij * H-(x_j)          (repression)
              + sum_{j != i} C_ij * H+(x_j)          (activation)
              + a_tf_i + paracrine_i
              - k_i * x_i

with a_i = alpha_i * exp(lam_i * u_i + sigma_i * eta_i - sigma_i^2 / 2)
(exponential-form extrinsic noise, OU-coloured; see Eqs. (2)-(4) of the
paper).  H+(x) = x^n/(s^n + x^n), H-(x) = s^n/(s^n + x^n), and every
cross-term uses the *source* gene's own threshold s_j.

For G = 2 with default matrices this reduces exactly to Eq. (1).
"""
import numpy as np


def _fit(v, G, fill=None):
    """Fit a per-gene parameter tuple to G genes (pad with the last value
    or ``fill``; truncate if longer)."""
    a = list(np.atleast_1d(np.asarray(v, float)))
    if len(a) < G:
        pad = fill if fill is not None else a[-1]
        a = a + [pad] * (G - len(a))
    return np.array(a[:G], float)


class GRN:
    """Precomputed regulation matrices + one-step vectorised update."""

    def __init__(self, cfg):
        self.cfg = cfg
        G = cfg.n_genes
        self.G = G
        self.alpha = _fit(cfg.alpha, G)
        self.lam = _fit(cfg.lam, G)
        self.rho = _fit(cfg.rho, G)
        self.s_self = _fit(cfg.s_self, G)
        self.k_deg = _fit(cfg.k_deg, G)
        self.n = cfg.n_hill
        self.sigma = _fit(cfg.sigma, G, 0.0)
        self.zeta = _fit(cfg.zeta, G, 1.0)
        self.a_tf = _fit(cfg.a_tf, G, 0.0)
        self.B, self.C = cfg.matrices()
        self.b_vec = cfg.b_vector()
        self.regulation = getattr(cfg, "regulation", "product")
        # Regulation thresholds.  The cross-term acting ON gene i and driven
        # BY gene j must use gene j's own half-effective concentration s_j,
        # exactly as in Eq. (1) of Huang et al. (2024):
        #     dx_1/dt = ... + b_1 * s_2^n / (s_2^n + x_2^n) - k_1 x_1
        # (an earlier version averaged s_i and s_j, which broke the exact
        # reduction whenever s_1 != s_2 -- e.g. the s_1 = 0.55 series of the
        # paper's Fig. 10).
        self.S_rep = np.zeros((G, G))
        self.S_act = np.zeros((G, G))
        s = self.s_self
        for i in range(G):
            for j in range(G):
                if i != j:
                    self.S_rep[i, j] = s[j]
                    self.S_act[i, j] = s[j]
                else:
                    # placeholder; diagonal edges carry weight 0 and are
                    # raised to the 0th power, but keep them finite so the
                    # division is well-defined for x = 0
                    self.S_rep[i, j] = 1.0
                    self.S_act[i, j] = 1.0
        self.sself_n = self.s_self ** self.n
        self.S_rep_n = self.S_rep ** self.n
        self.S_act_n = self.S_act ** self.n
        self.has_cross = bool(self.B.any() or self.C.any())

    # ------------------------------------------------------------------
    def activation_levels(self, x):
        """H+(x_i) per gene, shape (n_cells, G)."""
        xn = x ** self.n
        return xn / (self.sself_n[None, :] + xn)

    def step(self, x, u, eta, dt, paracrine=None):
        """One explicit-Euler step; returns new x (array copied).

        x, u, eta: (n_cells, G).  paracrine: (n_cells,) or None.
        """
        a = self.alpha[None, :] * np.exp(
            self.lam[None, :] * u + self.sigma[None, :] * eta
            - 0.5 * self.sigma[None, :] ** 2) + self.a_tf[None, :]
        act = self.activation_levels(x)
        prod = a * (self.rho[None, :] + (1 - self.rho[None, :]) * act)
        if self.has_cross:
            xn = x ** self.n
            Hn = self.S_rep_n[None, :, :] / (self.S_rep_n[None, :, :]
                                             + xn[:, None, :])
            if self.regulation == "product":
                # b_i * prod_j H-(x_j)^B_ij  (exact Eq. (1) for G = 2)
                w = np.where(self.B[None, :, :] > 0, self.B[None, :, :],
                             0.0)
                with np.errstate(divide="ignore"):
                    rep = np.prod(np.power(Hn, w), axis=2)
                prod = prod + self.b_vec[None, :] * rep
            else:
                rep = (self.B[None, :, :] * Hn).sum(axis=2)
                prod = prod + rep
            if self.C.any():
                actx = self.C[None, :, :] * xn[:, None, :] / (
                    self.S_act_n[None, :, :] + xn[:, None, :])
                prod = prod + actx.sum(axis=2)
        if paracrine is not None:
            prod = prod + np.asarray(paracrine)[:, None] * \
                self._paracrine_weight[None, :]
        dx = prod - self.k_deg[None, :] * x
        # safety clamp: production is bounded (<= a+b+c), so x is
        # mathematically bounded by (a+b+c)/k; the clamp only guards the
        # explicit-Euler scheme against rare extreme transients.
        return np.clip(x + dx * dt, 0.0, 1.0e3)

    # paracrine response weights (set by the simulator from SpaceConfig)
    _paracrine_weight = np.zeros(0)

    def set_paracrine_weight(self, w):
        self._paracrine_weight = np.asarray(w, float)

    # ------------------------------------------------------------------
    def ou_step(self, eta, dt, rng):
        """Euler-Maruyama step of the OU noise (Eq. 3)."""
        G = self.G
        zeta = np.maximum(self.zeta, 1e-12)
        return eta - eta * (dt / zeta)[None, :] + \
            np.sqrt(2 * dt / zeta)[None, :] * rng.standard_normal(eta.shape)
