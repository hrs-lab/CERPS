"""Notation alias layer: isolate the two naming conventions in use.

Two conventions coexist in this project and they are **opposite**:

==============  ===========================  =============================
convention      gene expression              epigenetic modification
==============  ===========================  =============================
``HSL``         ``x``                        ``u``
(CERPS native,  (HSL2024, CERPS core,
 Huang 2024,    ``reproduction/hsl2024``,
 paper1)        A2/direction1)
``PROPOSAL``    ``u``                        ``v``
(NSFC proposal  (proposal §2.1 / §3.1(2)①:
 §2.1, §3.1)    ``x = (u, v)``)
==============  ===========================  =============================

M1 (R6) fixed ONE internal convention -- ``HSL`` -- because the validated
reproduction layer, the batch sweeps and paper1 all speak it.  This module
therefore never renames anything internally; it only

1. **renames on the way in** (a proposal-side parameter block is translated
   into a :class:`~cerps.config.GRNConfig`), and
2. **renames on the way out** (a simulator state is exposed under
   proposal-side names).

The alias layer is *lossless for the HSL convention and lossy for the
proposal convention*: the proposal asks for epigenetic gates on the
repression strength and on the degradation rate, ``beta_i(v)`` and
``k_i(v)`` (proposal §3.1(2)①), which the native GRN does not carry.  Those
two gaps are reported explicitly by :meth:`ProposalRegulation.to_grn_config`
instead of being silently dropped -- a silently ignored option invalidates a
mechanism comparison without any error (audit 2026-10-20, defect D3).
``cerps/h_general.py`` closes them.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "HSL", "PROPOSAL",
    "NotationError",
    "to_hsl_arrays", "to_proposal_arrays",
    "state_view", "StateViewHSL", "StateViewProposal",
    "ProposalRegulation", "epigenetic_gate",
]


# ----------------------------------------------------------------------
# the two conventions
# ----------------------------------------------------------------------
HSL = "hsl"
PROPOSAL = "proposal"


class NotationError(ValueError):
    """Raised when the two naming conventions are mixed in one object."""


# ----------------------------------------------------------------------
# array-level swap (the whole content of the alias layer)
# ----------------------------------------------------------------------
def to_hsl_arrays(u_expr, v_epi):
    """``(u, v)`` (proposal) -> ``(x, u)`` (HSL).

    ``u_expr`` is gene expression and ``v_epi`` the epigenetic mark; the
    return order is ``(expression, epigenetic)``.  Both arguments are
    returned **as the same objects**, so this is a relabelling, not a copy:
    the two conventions differ only in which name goes first.
    """
    return u_expr, v_epi


def to_proposal_arrays(x_expr, u_epi):
    """``(x, u)`` (HSL) -> ``(u, v)`` (proposal)."""
    return x_expr, u_epi


# ----------------------------------------------------------------------
# simulator state views
# ----------------------------------------------------------------------
@dataclass
class StateViewHSL:
    """CERPS-native names for a simulator state."""
    x: object                       # (n, G) gene expression
    u: object                       # (n, G) epigenetic marks
    eta: object = None

    @property
    def expression(self):
        return self.x

    @property
    def epigenetic(self):
        return self.u


@dataclass
class StateViewProposal:
    """Proposal names for the *same* simulator state.

    ``u`` here is gene expression and ``v`` the epigenetic mark -- the exact
    opposite of :class:`StateViewHSL`.  The attributes hold references to
    the simulator arrays, so writing through the view writes to the model.
    """
    u: object                       # (n, G) gene expression
    v: object                       # (n, G) epigenetic marks
    eta: object = None

    @property
    def expression(self):
        return self.u

    @property
    def epigenetic(self):
        return self.v


def state_view(sim, convention=HSL):
    """Expose a simulator state under one convention's names.

    ``sim`` must expose ``x`` (expression) and ``u`` (epigenetic) arrays --
    every CERPS simulator does.  ``convention`` selects the naming only; the
    underlying arrays are shared, never copied.
    """
    if convention == HSL:
        return StateViewHSL(x=sim.x, u=sim.u, eta=getattr(sim, "eta", None))
    if convention == PROPOSAL:
        return StateViewProposal(u=sim.x, v=sim.u,
                                 eta=getattr(sim, "eta", None))
    raise NotationError(
        f"convention={convention!r} is not implemented; expected "
        f"{HSL!r} or {PROPOSAL!r}")


# ----------------------------------------------------------------------
# epigenetic gates  alpha_i(v), beta_i(v), k_i(v)
# ----------------------------------------------------------------------
_GATE_SIGN = {"alpha": +1.0, "beta": -1.0, "k": -1.0}


def epigenetic_gate(kind, base, lam, v):
    """The exponential epigenetic gate of proposal §3.1(2)①.

    ``kind`` selects which of the three gates is evaluated:

    ==========  ==========================================  ==============
    ``kind``    formula                                     HSL analogue
    ==========  ==========================================  ==============
    ``"alpha"`` ``alpha* exp(+lam*v)``                       Eq. (2) of
                                                            Huang 2024
    ``"beta"``  ``beta*  exp(-lam*v)``                       (no analogue)
    ``"k"``     ``k*     exp(-lam*v)``                       (no analogue)
    ==========  ==========================================  ==============

    The sign convention is the proposal's own: a higher mark level ``v``
    strengthens activation and weakens both repression and degradation.  In
    the HSL notation the same statement is "a higher ``u`` strengthens the
    basal production rate ``a_i = alpha_i*exp(lam_i*u_i)``" -- exactly
    Eq. (2) of Huang et al. (2024).  ``base`` may be a scalar or an array
    broadcastable against ``v``.
    """
    if kind not in _GATE_SIGN:
        raise NotationError(
            f"gate kind={kind!r} is not implemented; expected "
            f"'alpha', 'beta' or 'k'")
    return np.asarray(base, float) * np.exp(_GATE_SIGN[kind] * lam * v)


# ----------------------------------------------------------------------
# proposal-side regulation block -> CERPS GRNConfig
# ----------------------------------------------------------------------
@dataclass
class ProposalRegulation:
    """The regulation block of proposal §3.1(2)①, in proposal notation.

    The proposal writes, for gene ``i`` regulated by gene ``j``,

        G_i = alpha_i(v) * u_i^n/(s_i^n + u_i^n)
            + beta_i(v)  * s_j^n/(s_j^n + u_j^n)
            - k_i(v) * u_i

    with the exponential gates

        alpha_i(v) = alpha_i* exp(+lam_i v)
        beta_i(v)  = beta_i*  exp(-lam_i v)
        k_i(v)     = k_i*     exp(-lam_i v)

    (proposal §3.1(2)①, second displayed block).  ``alpha_star`` is the
    promotion strength ``alpha_i*``, ``beta_star`` the ``G x G`` matrix of
    repression strengths ``beta_ij*`` (``beta_star[i, j]`` = strength with
    which gene ``j`` represses gene ``i``), ``k_star`` the degradation
    rates, ``s`` the per-gene thresholds, ``n_hill`` the Hill exponents and
    ``lam`` the epigenetic impact strengths.

    Two ways of *using* this block exist and they are NOT equivalent:

    ``epi_mode="frozen"``
        Treat ``v`` as a constant (the proposal's own first step,
        ``H(u,v) == 0``).  The gates are evaluated once at ``v=v0`` and the
        resulting numbers are handed to the native GRN, whose ``lam`` is set
        to 0.  Fully supported, translation is exact.

    ``epi_mode="dynamic"``
        Keep ``v`` time-dependent.  The ``alpha`` gate is carried by the
        native term ``a_i = alpha_i*exp(lam_i*u_i)`` (Eq. (2) of Huang et
        al. 2024, with ``u`` the HSL name of the proposal's ``v``), but the
        ``beta`` and ``k`` gates have **no native counterpart** and are
        reported in ``gaps``.  Use :mod:`cerps.h_general` to close them.
    """
    n_genes: int = 2
    alpha_star: tuple = ()
    beta_star: tuple = ()            # (G, G) or flat row-major, [i, j]
    k_star: tuple = ()
    s: tuple = ()
    n_hill: tuple = ()
    lam: tuple = ()
    rho: tuple = ()                  # basal / max ratio of the self term
    epi_mode: str = "frozen"
    v0: tuple = ()                   # frozen mark level (frozen mode)

    # ------------------------------------------------------------------
    @staticmethod
    def _fit(v, G, fill=0.0):
        a = list(np.atleast_1d(np.asarray(v, float)).ravel())
        if len(a) < G:
            a = a + [fill] * (G - len(a))
        return np.array(a[:G], float)

    def _beta_matrix(self, G):
        if len(self.beta_star) == 0:
            return np.zeros((G, G))
        b = np.asarray(self.beta_star, float)
        if b.ndim == 1:                       # flat row-major (i, j)
            b = b.reshape(G, G)
        return b

    # ------------------------------------------------------------------
    def gates(self, v0=None):
        """Evaluate the three gates at ``v0``.

        ``v0`` is a length-``G`` array of mark levels (per cell when
        ``epi_mode="dynamic"``).  Returns ``(alpha_v, beta_v, k_v)`` where
        ``beta_v`` is the ``G x G`` matrix ``beta_ij* exp(-lam_i v)``.
        """
        G = self.n_genes
        v = self._fit(self.v0 if v0 is None else v0, G)
        lam = self._fit(self.lam, G)
        a = self._fit(self.alpha_star, G)
        k = self._fit(self.k_star, G)
        B = self._beta_matrix(G)
        # ``beta_ij(v) = beta_ij* exp(-lam_i v_i)``: the gate depends on the
        # mark of the REGULATED gene i (not on j), so v enters as a column
        # vector.  (An earlier version passed ones((1, G)), which silently
        # evaluated the gate at v == 1 -- wrong for every v0 != 1.)
        return (epigenetic_gate("alpha", a, lam, v),
                epigenetic_gate("beta", B, lam[:, None], v[:, None]),
                epigenetic_gate("k", k, lam, v))

    # ------------------------------------------------------------------
    def to_grn_config(self, v0=None, **overrides):
        """Translate into a CERPS :class:`~cerps.config.GRNConfig`.

        Returns ``(cfg, gaps)``.  ``gaps`` lists every proposal-side feature
        the native GRN cannot represent; an empty list means the translation
        is exact.

        Mapping (frozen mode, ``v`` replaced by the constant ``v0``):

        ============================  ====================================
        proposal                      CERPS native
        ============================  ====================================
        ``alpha_i* exp(lam_i v0)``    ``GRNConfig.alpha[i]``
        ``beta_ij* exp(-lam_i v0)``   ``GRNConfig.b_rep[i][j]``
        ``k_i* exp(-lam_i v0)``       ``GRNConfig.k_deg[i]``
        ``s_i``, ``n_i``              ``s_self[i]``, ``n_hill``
        ``rho_i``                     ``rho[i]``
        ============================  ====================================

        With ``epi_mode="dynamic"`` the ``alpha`` gate stays inside the
        native exponential (``lam[i]`` is kept) while the ``beta`` and ``k``
        gates cannot be expressed natively and are listed in ``gaps``.
        """
        from .config import GRNConfig
        G = self.n_genes
        gaps = []
        a, B, k = self.gates(v0)
        lam = self._fit(self.lam, G)
        s = self._fit(self.s, G, fill=1.0)
        n = self._fit(self.n_hill, G, fill=4.0)
        rho = self._fit(self.rho, G, fill=0.0)
        B_raw = self._beta_matrix(G)

        if self.epi_mode == "dynamic":
            # keep the epigenetic impact on the promotion rate; the gates on
            # repression and degradation are NOT representable natively
            lam_out = lam
            a_out = self._fit(self.alpha_star, G)
            B_out = B_raw
            k_out = self._fit(self.k_star, G)
            if B_raw.any():
                gaps.append(
                    "beta_i(v) = beta_i*exp(-lam_i v): the native GRN has a "
                    "CONSTANT repression edge GRNConfig.b_rep[i][j]; the "
                    "epigenetic gate on repression is dropped. Close it "
                    "with cerps.h_general.ProposalGRN(epi_gate_beta=True).")
            if self._fit(self.k_star, G).any():
                gaps.append(
                    "k_i(v) = k_i*exp(-lam_i v): the native GRN has a "
                    "CONSTANT degradation rate GRNConfig.k_deg[i]; the "
                    "epigenetic gate on degradation is dropped. Close it "
                    "with cerps.h_general.ProposalGRN(epi_gate_k=True).")
        else:
            lam_out = np.zeros(G)
            a_out, B_out, k_out = a, B, k

        n_hill = float(n[0]) if np.ptp(n) == 0 else tuple(n.tolist())
        cfg = GRNConfig(
            n_genes=G,
            alpha=tuple(a_out.tolist()),
            lam=tuple(lam_out.tolist()),
            rho=tuple(rho.tolist()),
            s_self=tuple(s.tolist()),
            k_deg=tuple(k_out.tolist()),
            n_hill=n_hill,
            b_rep=tuple(B_out.ravel().tolist()),
            regulation="additive",
            **overrides,
        )
        return cfg, gaps

    # ------------------------------------------------------------------
    def describe(self):
        """The block written out, for the paper text."""
        return (
            "F(x) = (G(x), H(x)) with x = (u, v); the proposal first sets\n"
            "H(u, v) == 0 and\n"
            "    G_i = alpha_i(v) u_i^n/(s_i^n+u_i^n)"
            " + beta_i(v) s_j^n/(s_j^n+u_j^n) - k_i(v) u_i,\n"
            "    alpha_i(v) = alpha_i* exp(+lam_i v),\n"
            "    beta_i(v)  = beta_i*  exp(-lam_i v),\n"
            "    k_i(v)     = k_i*     exp(-lam_i v)."
        )


# ----------------------------------------------------------------------
# mixed-convention guard
# ----------------------------------------------------------------------
def check_no_mixed_convention(arrays):
    """Raise if a single container mixes the two namings.

    ``arrays`` is a mapping ``name -> array``.  The guard rejects a mapping
    that uses ``u`` (HSL: epigenetic mark) together with ``u_expr``
    (proposal: gene expression): the two spellings denote *different* slots,
    so their coexistence signals a convention mix-up.
    """
    present = set(arrays)
    if ({"u", "u_epi"} & present) and ("u_expr" in present):
        raise NotationError(
            "mixed conventions: 'u'/'u_epi' means the epigenetic mark in the "
            "HSL convention, but 'u_expr' means gene expression in the "
            "proposal convention. Pick ONE convention per object.")
    return True
