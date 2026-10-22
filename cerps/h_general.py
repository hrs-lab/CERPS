"""The ``F(x) = (G(x), H(x))`` generalisation switch (proposal §3.1(2)①).

The proposal splits the intracellular vector field of ``x = (u, v)`` into
two components and then **turns the second one off as a first step**::

    F(x) = (G(x), H(x)),      assume  H(u, v) == 0,

    G_i = alpha_i(v) u_i^n/(s_i^n + u_i^n)
        + beta_i(v)  s_j^n/(s_j^n + u_j^n) - k_i(v) u_i

    alpha_i(v) = alpha_i* exp(+lam_i v)
    beta_i(v)  = beta_i*  exp(-lam_i v)
    k_i(v)     = k_i*     exp(-lam_i v)

This module implements that split as an explicit, *switchable* object:

:class:`ProposalGRN`
    the ``G`` component.  Additive Hill form with all three epigenetic
    gates.  The ``beta`` and ``k`` gates are OFF by default, in which case
    ``G`` reproduces the native :class:`~cerps.grn.GRN` **exactly** (verified
    in ``tests/test_m2_generalisation.py``, V2); turning them ON is what
    makes the model the proposal's rather than the paper's.

:class:`ProposalH`
    the ``H`` component, with the three-level switch the proposal implies
    (``H≡0`` first, then the general form).  ``mode="zero"`` freezes the
    epigenetic state, which is the proposal's own starting assumption;
    ``mode="general"`` is the validated CERPS continuous dynamics.

Sign conventions and the RACIPE state enumeration follow the M1 outcome
(R5): the net regulation direction of gene ``j`` on gene ``i`` is
``sign(alpha_ij - beta_ij)``, NOT ``beta > alpha*(1-rho)`` -- see
:meth:`ProposalGRN.net_direction` for the derivation and
:func:`racipe_states` for the enumeration.

Everything here is vectorised over cells exactly like the native GRN, so a
:class:`ProposalGRN` can be dropped into a
:class:`~cerps.simulator.Simulator` via :func:`attach`.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "H_MODES", "G_MODES",
    "ProposalGRN", "ProposalH",
    "racipe_states", "net_direction_matrix",
    "attach",
]


# ----------------------------------------------------------------------
# the switches
# ----------------------------------------------------------------------
#: ``H(u, v)`` levels.  ``zero`` is the proposal's own first step; ``linear``
#: is the expression-decoupled relaxation; ``general`` is the CERPS
#: continuous dynamics (w0 + w1*act) writing / (e0 + e1*(1-act)) erasing.
H_MODES = ("zero", "linear", "general")

#: ``G`` composition levels.  ``additive`` is the proposal's form (both
#: terms carry a plus sign); ``multiplicative`` is the thermodynamic
#: / Buchler-Gerland-Hwa AND-gate form, kept for comparison.
G_MODES = ("additive", "multiplicative")


# ----------------------------------------------------------------------
# RACIPE enumeration and the net-direction criterion  (M1 / R5)
# ----------------------------------------------------------------------
def racipe_states():
    """The three regulation states an edge can take (M1, R5).

    ============  =======  =======  =====================================
    state         alpha    beta     meaning
    ============  =======  =======  =====================================
    ``activation``  > 0      0      gene j raises the production of gene i
    ``repression``   0      > 0     gene j lowers it
    ``none``         0       0      no regulation
    ============  =======  =======  =====================================

    The enumeration is ``(alpha, beta) in {(>0,0), (0,>0), (0,0)}`` and it
    is the ONLY admissible set: any sample with ``alpha > 0`` **and**
    ``beta > 0`` is a mixed edge that the proposal's additive form does not
    distinguish from a pure one, because both terms enter with a plus sign.
    Sampling such edges is the bug that puts "repression" draws inside the
    activation region -- see :meth:`ProposalGRN.net_direction`.
    """
    return (("activation", 1.0, 0.0),
            ("repression", 0.0, 1.0),
            ("none", 0.0, 0.0))


def net_direction_matrix(alpha, beta):
    """``sign(alpha - beta)`` for every edge of a ``G x G`` regulation.

    This is the net-regulation criterion established in M1 (R5).  For the
    additive form

        production_i = alpha_ij H(x_j) + beta_ij (1 - H(x_j))

    the two limits are ``production_i(u_j = 0) = beta_ij`` and
    ``production_i(u_j = +inf) = alpha_ij`` (both terms carry a plus sign),
    so the edge activates iff ``alpha_ij > beta_ij``.  In particular a
    **positive** ``beta`` alone is NOT repression unless ``alpha = 0``: with
    ``alpha = 2, beta = 1`` the edge is a weak *activator*.

    Returns an integer array with ``-1`` = repression, ``0`` = none,
    ``+1`` = activation.
    """
    a = np.asarray(alpha, float)
    b = np.asarray(beta, float)
    return np.sign(a - b).astype(int)


# ----------------------------------------------------------------------
# the G component
# ----------------------------------------------------------------------
@dataclass
class ProposalGRN:
    """The ``G`` component of the proposal's ``F(x) = (G(x), H(x))``.

    Parameters
    ----------
    n_genes : int
    alpha : sequence
        ``alpha_i*``, the promotion strength of gene ``i``'s own activation
        term (``u_i^n/(s_i^n+u_i^n)``).
    beta : (G, G) array
        ``beta_ij*``, the strength with which gene ``j`` represses gene
        ``i`` (the term ``s_j^n/(s_j^n+u_j^n)``).
    k_deg : sequence
        ``k_i*``, the degradation rate.
    s : sequence
        per-gene thresholds ``s_i``; a cross term driven by gene ``j`` uses
        ``s_j`` (Eq. (1) of Huang et al. 2024 -- see the note in
        :class:`~cerps.grn.GRN`).
    n_hill : float
        Hill exponent (a scalar, as in the proposal).
    lam : sequence
        epigenetic impact strength ``lam_i``; the sign convention of the
        proposal is ``alpha`` up, ``beta`` and ``k`` down.
    rho : sequence
        basal / max activation ratio of the self term, as in the native
        GRN.  ``rho = 0`` is the proposal's literal form.
    c_act : (G, G) array
        optional activating edges (gene ``j`` activates gene ``i``); the
        proposal does not use them but the native comparison does.
    sigma, zeta, a_tf : sequences
        extrinsic-noise intensity, OU relaxation time, externally induced
        transcription factor -- all as in the native GRN.

    Gates
    -----
    ``epi_gate_alpha``
        Default **True**: the promotion rate carries the exponential
        ``alpha_i(v) = alpha_i* exp(+lam_i v)``, which is the same term as
        Eq. (2) of Huang et al. 2024 and is therefore the validated
        behaviour.  Set it to False for the proposal's literal frozen form
        (``H == 0``, so ``v`` is a constant and the gate has been baked into
        ``alpha``) -- equivalently, use
        :meth:`cerps.notation.ProposalRegulation.to_grn_config` with
        ``epi_mode="frozen"``.
    ``epi_gate_beta`` / ``epi_gate_k``
        Default **False**.  Turning them ON applies the proposal's
        exponential gates ``beta_ij(v) = beta_ij* exp(-lam_i v)`` and
        ``k_i(v) = k_i* exp(-lam_i v)`` at every step -- precisely the two
        gates the native GRN cannot express.  With all gates at their
        defaults the object reproduces the native
        :class:`~cerps.grn.GRN` **exactly** (test V2).
    """

    n_genes: int = 2
    alpha: tuple = (0.4, 0.4)
    beta: tuple = ()
    k_deg: tuple = (1.0, 1.0)
    s: tuple = (0.5, 0.5)
    n_hill: float = 4.0
    lam: tuple = (1.9, 1.9)
    rho: tuple = (0.1, 0.1)
    c_act: tuple = ()
    sigma: tuple = (0.05, 0.05)
    zeta: tuple = (1.0, 1.0)
    a_tf: tuple = (0.0, 0.0)
    regulation: str = "additive"
    epi_gate_beta: bool = False
    epi_gate_k: bool = False
    epi_gate_alpha: bool = True

    # ------------------------------------------------------------------
    def __post_init__(self):
        # the native GRN spells the AND-gate form "product"; accept it as an
        # alias so that ProposalGRN.from_native() never has to rename
        if self.regulation == "product":
            self.regulation = "multiplicative"
        if self.regulation not in G_MODES:
            raise ValueError(
                f"regulation={self.regulation!r} is not implemented; "
                f"expected one of {G_MODES} (or 'product' for "
                f"'multiplicative')")
        self._build()

    # ------------------------------------------------------------------
    def _build(self):
        G = self.n_genes
        self.G = G
        self.alpha_v = _fit(self.alpha, G)
        self.lam_v = _fit(self.lam, G)
        self.rho_v = _fit(self.rho, G, 0.0)
        self.s_v = _fit(self.s, G, 1.0)
        self.k_v = _fit(self.k_deg, G)
        self.sigma_v = _fit(self.sigma, G, 0.0)
        self.zeta_v = _fit(self.zeta, G, 1.0)
        self.a_tf_v = _fit(self.a_tf, G, 0.0)
        self.n = float(self.n_hill)
        B = np.asarray(self.beta, float) if len(self.beta) else \
            np.zeros((G, G))
        if B.ndim == 1:
            B = B.reshape(G, G)
        self.beta_v = B
        C = np.asarray(self.c_act, float) if len(self.c_act) else \
            np.zeros((G, G))
        if C.ndim == 1:
            C = C.reshape(G, G)
        self.c_act_v = C
        # thresholds: a cross term acting on i and driven by j uses s_j
        self.S = np.zeros((G, G))
        for i in range(G):
            for j in range(G):
                self.S[i, j] = self.s_v[j] if i != j else 1.0
        self.sself_n = self.s_v ** self.n
        self.S_n = self.S ** self.n
        self.has_cross = bool(self.beta_v.any() or self.c_act_v.any())
        # effective epigenetic impact: zeroed (not removed) when the alpha
        # gate is off, so the exponential keeps the native term order
        self.lam_eff = self.lam_v if self.epi_gate_alpha else np.zeros(G)
        self._paracrine_weight = np.zeros(0)

    # ------------------------------------------------------------------
    # gate evaluation ---------------------------------------------------
    def _gate(self, lam_coeff, v):
        """Return ``exp(lam_coeff * lam_i * v)`` broadcast against genes.

        Two accepted shapes for ``v``:
          * ``(M,)``  -- a list of *scalar* mark levels; output is
            ``(M, G)`` (the gate evaluated at each level, every gene);
          * ``(n, G)`` -- the per-cell epigenetic state; output is
            ``(n, G)`` (the gate of gene ``j`` uses mark level ``v[:, j]``).
        ``lam_coeff`` is ``+1`` for the alpha gate and ``-1`` otherwise.
        """
        vv = np.asarray(v, float)
        if vv.ndim == 1:
            return np.exp(lam_coeff * self.lam_v[None, :] * vv[:, None])
        return np.exp(lam_coeff * self.lam_v[None, :] * vv)

    def alpha_of(self, v):
        """``alpha_i(v) = alpha_i* exp(+lam_i v)`` if the gate is on,
        else the constant ``alpha_i*``.

        ``v`` is the epigenetic state (the HSL name ``u``).  Shape: ``(M,)``
        (list of mark levels) or ``(n, G)`` (per-cell).
        """
        if self.epi_gate_alpha:
            return self.alpha_v[None, :] * self._gate(+1.0, v)
        return self.alpha_v

    def beta_of(self, v):
        """``beta_ij(v) = beta_ij* exp(-lam_i v)`` if the gate is on.

        Output shape ``(M, G, G)`` for a ``(M,)`` mark-level list, else
        ``(n, G, G)``.
        """
        if self.epi_gate_beta:
            vv = np.asarray(v, float)
            gate = self._gate(-1.0, v)               # (M, G) or (n, G)
            if vv.ndim == 1:
                return self.beta_v[None, :, :] * gate[:, :, None]
            return self.beta_v[None, :, :] * gate[:, :, None]
        return self.beta_v

    def k_of(self, v):
        """``k_i(v) = k_i* exp(-lam_i v)`` if the gate is on.

        Shape ``(M, G)`` for a ``(M,)`` mark-level list, else ``(n, G)``.
        """
        if self.epi_gate_k:
            return self.k_v[None, :] * self._gate(-1.0, v)
        return self.k_v

    # ------------------------------------------------------------------
    # explicit (alpha, beta) pairs -- the RACIPE view -------------------
    def alpha_matrix(self):
        """``G x G`` matrix of ``alpha_ij`` (the ``u_j -> +inf`` limits).

        The diagonal is the self-activation term ``alpha_i*``; the off
        diagonal is filled from ``c_act`` when activating edges exist.
        """
        A = np.diag(self.alpha_v).copy()
        A = A + np.where(self.c_act_v > 0, self.c_act_v, 0.0)
        return A

    def beta_matrix(self):
        """``G x G`` matrix of ``beta_ij`` (the ``u_j -> 0`` limits).

        The diagonal carries the self-term's basal level
        ``alpha_i* rho_i`` (which is why the self edge has
        ``alpha > beta`` and reads as activation); the off diagonal is
        ``beta_ij*``.
        """
        Bm = self.beta_v.copy()
        np.fill_diagonal(Bm, self.alpha_v * self.rho_v)
        return Bm

    def net_direction(self):
        """``sign(alpha - beta)`` per edge (M1 / R5).  See
        :func:`net_direction_matrix` for the derivation."""
        return net_direction_matrix(self.alpha_matrix(), self.beta_matrix())

    # ------------------------------------------------------------------
    # dynamics ----------------------------------------------------------
    def activation_levels(self, x):
        """``H+(x_i)`` per gene, shape ``(n, G)``."""
        xn = x ** self.n
        return xn / (self.sself_n[None, :] + xn)

    def step(self, x, u, eta, dt, paracrine=None):
        """One explicit-Euler step; returns the new ``x``.

        The signature matches :meth:`cerps.grn.GRN.step` exactly so the two
        are interchangeable (V2).  ``u`` is the epigenetic state; under the
        HSL convention it is the same array the native GRN receives.
        """
        # Basal production rate, Eq. (2) of Huang et al. 2024:
        #     a_i = alpha_i * exp(lam_i*u_i + sigma_i*eta_i - sigma_i^2/2)
        #         + a_tf_i
        # The sum inside the exponential is written in the SAME order as the
        # native implementation: floating-point addition is not associative,
        # so regrouping the three terms costs one ulp and would break the
        # exact reduction asserted by test V2.  Turning the alpha gate off
        # is done by zeroing ``lam_eff`` (exact: 0.0 * u == 0.0 for the
        # bounded mark levels), never by dropping the term.
        z = (self.lam_eff[None, :] * u + self.sigma_v[None, :] * eta
             - 0.5 * self.sigma_v[None, :] ** 2)
        a = self.alpha_v[None, :] * np.exp(z) + self.a_tf_v[None, :]
        act = self.activation_levels(x)
        prod = a * (self.rho_v[None, :] + (1 - self.rho_v[None, :]) * act)
        if self.has_cross:
            xn = x ** self.n
            Hn = self.S_n[None, :, :] / (self.S_n[None, :, :] + xn[:, None, :])
            if self.regulation == "additive":
                # proposal: sum_j beta_ij(v) H-(x_j)   (both signs plus)
                B = self.beta_of(u)
                rep = (B * Hn).sum(axis=2)
            else:
                # multiplicative / AND gate: prod_j H-(x_j)^beta_ij.
                # The epigenetic gate is not defined for this composition
                # (a per-cell exponent would break the product identity),
                # so it is rejected rather than silently ignored.
                if self.epi_gate_beta:
                    raise NotImplementedError(
                        "regulation='multiplicative' does not support "
                        "epi_gate_beta: the proposal's exponential gate "
                        "changes the exponent per cell, which is not the "
                        "AND-gate form. Use regulation='additive'.")
                B = self.beta_v
                w = np.where(B[None, :, :] > 0, B[None, :, :], 0.0)
                with np.errstate(divide="ignore"):
                    rep = np.prod(np.power(Hn, w), axis=2)
            prod = prod + rep
            if self.c_act_v.any():
                actx = self.c_act_v[None, :, :] * xn[:, None, :] / (
                    self.S_n[None, :, :] + xn[:, None, :])
                prod = prod + actx.sum(axis=2)
        if paracrine is not None:
            prod = prod + np.asarray(paracrine)[:, None] * \
                self._paracrine_weight[None, :]
        k = self.k_of(u)
        dx = prod - k * x
        return np.clip(x + dx * dt, 0.0, 1.0e3)

    # ------------------------------------------------------------------
    def set_paracrine_weight(self, w):
        self._paracrine_weight = np.asarray(w, float)

    def ou_step(self, eta, dt, rng):
        zeta = np.maximum(self.zeta_v, 1e-12)
        return eta - eta * (dt / zeta)[None, :] + \
            np.sqrt(2 * dt / zeta)[None, :] * rng.standard_normal(eta.shape)

    # ------------------------------------------------------------------
    # validation --------------------------------------------------------
    def validate(self):
        """Check the sampling constraints that make the network capable of
        multistability.  Returns a list of findings (empty = clean).

        Two checks, both from M1 (V3b / V3c of ``verify_m1.py``):

        1. **Empirical bistability window** ``a/k in [s, ~3 s]``.  Below the
           window the self-activation cannot overcome degradation, so no
           stable high state exists; at ``a/k >= 4 s`` both genes sit high
           simultaneously and the toggle disappears (the double-high state
           becomes the only attractor).  [经验规律，数值扫描]
        2. **Positive-feedback-loop requirement** (Thomas 1973 / 1981,
           confirmed numerically in V3c): a necessary condition for
           multistationarity is that every critical non-trivial loop carry
           a positive sign.
        """
        out = []
        A = self.alpha_matrix()
        Bm = self.beta_matrix()
        for i in range(self.G):
            a = A[i, i]
            k = self.k_v[i]
            s = self.s_v[i]
            if k > 0:
                ratio = a / k
                if ratio < s:
                    out.append(
                        f"gene {i}: a/k = {ratio:.3f} < s = {s:.3f} -- the "
                        f"self-activation cannot overcome degradation; no "
                        f"stable high state (V3b lower bound).")
                elif ratio >= 4.0 * s:
                    out.append(
                        f"gene {i}: a/k = {ratio:.3f} >= 4 s = {4*s:.3f} -- "
                        f"both genes can be high simultaneously; the toggle "
                        f"bistability is lost (V3b upper bound).")
        return out

    # ------------------------------------------------------------------
    @classmethod
    def from_native(cls, grn, **overrides):
        """Build a :class:`ProposalGRN` from a native
        :class:`~cerps.grn.GRN`, preserving the current configuration.

        The reverse map is exact: the native ``alpha``/``lam``/``rho``
        become ``alpha``/``lam``/``rho``, ``b_rep`` becomes ``beta`` and
        ``c_act`` becomes ``c_act``.  ``epi_gate_alpha`` defaults to True,
        which is what carries the native ``lam`` term; the two gates the
        native GRN does not have (``beta``, ``k``) stay off, so the result
        reproduces the native dynamics exactly.  Pass
        ``epi_gate_alpha=False`` to obtain the proposal's literal frozen
        form instead.

        Precondition for exactness (M2b, verified by test V7): when the
        inherited composition is the AND-gate form (native
        ``regulation='product'``, the default), the native cross term is
        ``b_i * prod_j H-(x_j)^B_ij`` -- Eq. (1) of Huang et al. (2024)
        -- i.e. it carries the per-target prefactor ``b_i = b_base[i]``.
        The proposal's multiplicative form
        ``prod_j H-(x_j)^beta_ij`` is bounded by 1 and therefore cannot
        express a prefactor ``b_i != 1``.  The reduction is exact only
        for ``b_base`` identically 1 (the default); a non-unit
        ``b_base`` combined with the AND-gate form raises ``ValueError``
        rather than silently dropping ``b_vec``.  In the additive
        composition the native dynamics never uses ``b_base``, so it is
        simply not carried over (test V7b).  Overriding
        ``regulation='additive'`` yields the proposal's literal form --
        a generalisation, not an exact reduction.
        """
        B, C = grn.cfg.matrices()
        # The three epigenetic gates are explicit parameters of the
        # constructor, so extract them from ``overrides`` (if present) to
        # avoid "multiple values for keyword argument" when the caller
        # passes them.  ``epi_gate_alpha`` defaults to True -- that is what
        # carries the native ``lam`` term -- and the two gates the native
        # GRN cannot express (``beta``, ``k``) stay OFF, so the result
        # reproduces the native dynamics bit for bit.
        gate_alpha = overrides.pop("epi_gate_alpha", True)
        gate_beta = overrides.pop("epi_gate_beta", False)
        gate_k = overrides.pop("epi_gate_k", False)
        # ``regulation`` must be overridable: the proposal's G is the
        # ADDITIVE Hill form, and ``epi_gate_beta`` is only defined there,
        # whereas the native GRN defaults to "product".
        regulation = overrides.pop(
            "regulation", getattr(grn, "regulation", "product"))
        # M2b guard: the native AND-gate (product) form multiplies the
        # whole cross-repression term by the per-target prefactor b_i
        # (Eq. (1) of Huang et al. 2024:
        #     dx_i/dt = ... + b_i * prod_j H-(x_j)^B_ij - k_i x_i ),
        # whereas the proposal's multiplicative form is
        #     prod_j H-(x_j)^beta_ij ,
        # which is bounded by 1 and therefore CANNOT carry a prefactor
        # b_i != 1.  Silently dropping b_vec would break the exact
        # reduction asserted above (a single step already deviates by
        # |b_i - 1|, e.g. 0.25 for b_base = 2), so refuse instead.  The
        # guard is inert when the cross term never enters the dynamics
        # (no B and no C edges), because then b_vec is unused.
        if regulation in ("product", "multiplicative"):
            b_vec = np.asarray(
                getattr(grn, "b_vec", np.ones(grn.G)), float)
            if grn.has_cross and not np.allclose(b_vec, 1.0, rtol=0.0,
                                                 atol=1e-12):
                raise ValueError(
                    "from_native: exact reduction unavailable for "
                    "regulation='product' with b_base != 1.  The native "
                    "AND-gate term is b_i * prod_j H-(x_j)^B_ij (Eq. (1) "
                    "of Huang et al. 2024) and carries the per-target "
                    "prefactor b_i = b_base[i]; the proposal's "
                    "multiplicative form prod_j H-(x_j)^beta_ij <= 1 "
                    "cannot express that prefactor.  Either keep "
                    "b_base = 1 (the default) for an exact reduction, or "
                    "pass regulation='additive' for the proposal's "
                    "literal form (a generalisation, not a reduction).")
        return cls(
            n_genes=grn.G,
            alpha=tuple(grn.alpha.tolist()),
            beta=tuple(B.ravel().tolist()),
            k_deg=tuple(grn.k_deg.tolist()),
            s=tuple(grn.s_self.tolist()),
            n_hill=grn.n,
            lam=tuple(grn.lam.tolist()),
            rho=tuple(grn.rho.tolist()),
            c_act=tuple(C.ravel().tolist()),
            sigma=tuple(grn.sigma.tolist()),
            zeta=tuple(grn.zeta.tolist()),
            a_tf=tuple(grn.a_tf.tolist()),
            regulation=regulation,
            epi_gate_alpha=gate_alpha,
            epi_gate_beta=gate_beta,
            epi_gate_k=gate_k,
            **overrides,
        )


# ----------------------------------------------------------------------
# the H component
# ----------------------------------------------------------------------
@dataclass
class ProposalH:
    """The ``H`` component of the proposal's ``F(x) = (G(x), H(x))``.

    ``H`` drives the epigenetic state ``v``.  The proposal only ever states
    its *first* level (``H == 0``); the remaining levels are the natural
    completions.  A critical precision, recorded in the M2 deliverable:

    * The proposal writes ``H(u, v) == 0`` as its starting assumption, but
      HSL2024's own epigenetic dynamics -- which the proposal claims to
      *generalise* -- is :func:`cerps.epigenetics.continuous_step`:
          dv/dt = (w0 + w1*act)(1 - v) - (e0 + e1*(1 - act)) v .
      So ``H == 0`` as "step 1" is imprecise; the mathematically correct
      reading is: **step 1 = the HSL2024 baseline**, and the proposal then
      generalises *away* from it.  We therefore make ``mode="zero"``
      reproduce ``continuous_step`` bit for bit, which guarantees
      zero-regression when the proposal layer is switched on (the existing
      CERPS/HSL2024 regtests must still pass with ``mode="zero"``).

    Levels
    ------
    ``mode="zero"``   -- HSL2024 baseline (activity-dependent relaxation).
                        Identical to ``continuous_step``.
    ``mode="linear"`` -- baseline + linear expression->epigenetics term
                        ``h_gain * polarity * (act - 1/2)``.
    ``mode="general"``-- baseline + nonlinear term
                        ``h_gain * polarity * (act*(1-v) - (1-act)*v)``.

    Mark polarity: the code does not know whether ``v`` is an activating
    or a repressive mark; ``polarity=-1`` reverses the coupling so that an
    active gene LOSES a repressive mark, mirroring the fixed point about
    1/2.  (``polarity=-1`` is a proposal-only extension and is NOT
    expected to equal ``continuous_step``.)
    """
    w0: float = 0.01
    w1: float = 0.06
    e0: float = 0.01
    e1: float = 0.06
    mode: str = "zero"
    polarity: float = 1.0
    h_gain: float = 0.0

    @classmethod
    def from_epi(cls, cfg, mode="zero", **overrides):
        """Build a :class:`ProposalH` from an :class:`EpiConfig`, inheriting
        ``(w0, w1, e0, e1)`` so that ``mode='zero'`` equals ``continuous_step``
        exactly regardless of how the simulator was configured.
        """
        return cls(w0=cfg.w0, w1=cfg.w1, e0=cfg.e0, e1=cfg.e1,
                   mode=mode, **overrides)

    def __post_init__(self):
        if self.mode not in H_MODES:
            raise ValueError(
                f"mode={self.mode!r} is not implemented; expected one of "
                f"{H_MODES}")
        if self.polarity not in (1.0, -1.0):
            raise ValueError(
                f"polarity={self.polarity!r} is not implemented; expected "
                f"1.0 (activating mark) or -1.0 (repressive mark)")

    def step(self, v, act, dt):
        """One explicit-Euler step of ``dv/dt = H(u, v)``.

        Returns the new ``v``.  ``mode="zero"`` is *identical* to
        :func:`cerps.epigenetics.continuous_step`; ``mode="linear"`` and
        ``mode="general"`` add, on top of that baseline, an
        expression->epigenetics generalisation term controlled by
        ``h_gain`` and ``polarity``.
        """
        v = np.asarray(v, float)
        a = np.clip(np.asarray(act, float), 0.0, 1.0)
        # --- HSL2024 baseline (== continuous_step) ---
        if self.polarity > 0:
            w = self.w0 + self.w1 * a
            e = self.e0 + self.e1 * (1.0 - a)
        else:
            w = self.w0 + self.w1 * (1.0 - a)
            e = self.e0 + self.e1 * a
        dv = w * (1.0 - v) - e * v
        # --- proposal generalisation terms (on top of the baseline) ---
        if self.mode == "linear":
            dv = dv + self.h_gain * self.polarity * (a - 0.5)
        elif self.mode == "general":
            dv = dv + self.h_gain * self.polarity * (
                a * (1.0 - v) - (1.0 - a) * v)
        return np.clip(v + dv * dt, 0.0, 1.0)


# ----------------------------------------------------------------------
# attachment
# ----------------------------------------------------------------------
def attach(sim, pgrn=None, ph=None, reclassify=True):
    """Drop a proposal-side ``(G, H)`` pair into a running simulator.

    ``sim`` is a :class:`~cerps.simulator.Simulator`.  ``pgrn`` replaces
    ``sim.grn`` (both expose the same ``step`` / ``activation_levels`` /
    ``set_paracrine_weight`` / ``ou_step`` interface) and ``ph`` is stored
    as ``sim.proposal_h`` for the epigenetic update; ``None`` leaves the
    corresponding component alone.

    After the swap the derived quantities that were computed from the old
    GRN in ``__init__`` are recomputed: the legacy classification labels and
    the resolved stemness-switch threshold.  ``reclassify=False`` skips the
    label recomputation (useful when the caller manages labels itself).

    Returns ``sim``.
    """
    from . import rates as rates_mod
    if pgrn is not None:
        if pgrn.G != sim.grn.G:
            raise ValueError(
                f"gene-count mismatch: simulator has {sim.grn.G} genes, "
                f"proposal GRN has {pgrn.G}")
        sim.grn = pgrn
        if sim.space is not None and sim.space.paracrine:
            pgrn.set_paracrine_weight(sim.cfg.space.response)
    if ph is not None:
        sim.proposal_h = ph
    if reclassify:
        sim.label = rates_mod.classify_legacy(
            sim.grn, sim.x, np.zeros(sim.n, np.int8),
            sim.cfg.rates.x_threshold)
    sim.K_switch = rates_mod.resolve_K(sim.cfg.rates, sim.grn)
    return sim


# ----------------------------------------------------------------------
def _fit(v, G, fill=0.0):
    a = list(np.atleast_1d(np.asarray(v, float)).ravel())
    if len(a) < G:
        a = a + [fill] * (G - len(a))
    return np.array(a[:G], float)
