"""Default parameters of Huang, Situ & Lei (2024), J. Theor. Biol. 577:111664.

All values follow Table 2 of the paper.  Where the released C++ code
(https://github.com/jinzhilei/Cell-type-transition, GPL-3.0) differs from
the text of the paper, the code convention is followed and documented in
``docs/protocol_notes.md``.
"""
from dataclasses import dataclass, field, replace


@dataclass
class Params:
    # ------------------------------------------------------------------
    # Gene regulation network (Eqs. 1-4)
    # ------------------------------------------------------------------
    b1: float = 1.0     # strength of mutual inhibition acting on gene 1
    b2: float = 1.0
    k1: float = 1.0     # degradation rates  [1/h]
    k2: float = 1.0
    s1: float = 0.5     # half-effective concentrations [AU]
    s2: float = 0.5
    rho1: float = 0.1   # basal / maximum expression ratio
    rho2: float = 0.1
    n: float = 4.0      # Hill coefficient
    alpha1: float = 0.4  # mean expression rates [AU/h]
    alpha2: float = 0.4
    lam1: float = 1.9   # impact coefficient of epigenetic modification
    lam2: float = 1.9
    zeta1: float = 1.0  # OU relaxation coefficients [h]
    zeta2: float = 1.0
    sigma: float = 0.05  # extrinsic noise intensity
    a_tf: float = 0.0   # extra TF acting on X1 self-activation (Eq. 9)

    # ------------------------------------------------------------------
    # Cell kinetics (G0 cell cycle model, Table 1)
    # ------------------------------------------------------------------
    mu0: float = 0.0002    # apoptosis rate in proliferating phase [1/h]
    kappa0: float = 0.01   # removal rate from resting phase [1/h]
    beta0: float = 0.04    # maximum proliferation rate [1/h]
    tau0: float = 4.0      # proliferation duration of SC [h]
    theta0: float = 200.0  # half-effective SC number for feedback
    max_div_ta: int = 15   # TA cells die at the 15th division (C++ rule)

    # ------------------------------------------------------------------
    # Epigenetic-state inheritance (Eqs. 6-8)
    # ------------------------------------------------------------------
    m0: float = 60.0
    m1: float = 0.08
    m2: float = 0.8
    m3: float = 2.0
    m4: float = 2.0

    # ------------------------------------------------------------------
    # Numerics / protocol
    # ------------------------------------------------------------------
    dt: float = 0.25      # time step [h]
    n_init: int = 100     # initial cell number
    max_time: float = 4000.0

    # phenotype classification thresholds (C++ UpdateCellType)
    x_hi: float = 0.5

    def with_(self, **kw) -> "Params":
        return replace(self, **kw)

    # codes used by the C++ implementation
    SC = 0
    TA1 = 1
    TA2 = 3   # NOTE: the C++ code uses 3 (not 2) for TA2


DEFAULT = Params()
