"""Configuration of a CERPS model.

The dataclass below defines every component of the framework:

GRN            : N-gene regulation network (self-activation + signed edges)
epigenetics    : continuous intracellular dynamics + division kernel
division       : gene-expression partitioning at mitosis
rates          : continuous phenotype-dependent cell kinetics
               (``rate_mode="legacy"`` reproduces Table 1 of the paper)
space          : 3D positions, mechanics, contact inhibition, paracrine field
"""
from dataclasses import dataclass, field, replace
import numpy as np


# ----------------------------------------------------------------------
# GRN
# ----------------------------------------------------------------------
@dataclass
class GRNConfig:
    n_genes: int = 2
    alpha: tuple = (0.4, 0.4)        # base expression rates [AU/h]
    lam: tuple = (1.9, 1.9)          # epigenetic impact on expression rate
    rho: tuple = (0.1, 0.1)          # basal / max activation ratio
    s_self: tuple = (0.5, 0.5)       # self-activation thresholds
    k_deg: tuple = (1.0, 1.0)        # degradation rates [1/h]
    n_hill: float = 4.0
    # per-gene basal expression without repression (paper b1, b2)
    b_base: tuple = ()
    # repression edges B[i][j] > 0: gene j represses gene i
    b_rep: tuple = ()                # default: full mutual repression
    # cross-term composition: "product" (default; reduces exactly to
    # Eq. (1) for 2 genes) or "additive"
    regulation: str = "product"
    # optional activating edges c[i][j] (gene j activates gene i)
    c_act: tuple = ()
    sigma: tuple = (0.05, 0.05)      # extrinsic noise intensities
    zeta: tuple = (1.0, 1.0)         # OU relaxation times [h]
    a_tf: tuple = (0.0, 0.0)         # externally induced TF (per gene)

    def matrices(self):
        G = self.n_genes
        B = np.zeros((G, G)) if not self.b_rep else np.array(self.b_rep,
                                                             float)
        C = np.zeros((G, G)) if not self.c_act else np.array(self.c_act,
                                                            float)
        if not self.b_rep:                       # full mutual repression
            B[:] = 1.0
            np.fill_diagonal(B, 0.0)
        return B, C

    def b_vector(self):
        G = self.n_genes
        if not self.b_base:
            return np.ones(G)
        a = list(np.atleast_1d(np.asarray(self.b_base, float)))
        if len(a) < G:
            a = a + [a[-1]] * (G - len(a))
        return np.array(a[:G], float)


# ----------------------------------------------------------------------
# epigenetics
# ----------------------------------------------------------------------
@dataclass
class EpiConfig:
    # --- continuous intracellular dynamics (GRN -> epigenetics coupling)
    continuous: bool = False
    w0: float = 0.01                # basal mark-writing rate [1/h]
    w1: float = 0.06                # activity-dependent writing rate
    e0: float = 0.01                # basal mark-erasing rate [1/h]
    e1: float = 0.06                # inactivity-dependent erasing rate
    # --- division inheritance kernel
    kernel: str = "beta"            # beta | gaussian | cpg | deterministic
    # beta kernel (paper Eqs. 6-8)
    m0: float = 60.0
    m1: float = 0.08
    m2: float = 0.8
    m3: float = 2.0
    m4: float = 2.0
    # gaussian kernel
    sigma_g: float = 0.045
    # cpg kernel: M explicit markable sites per gene
    cpg_sites: int = 60
    # cpg_mode="moment"      per-site marking prob = phi(v) -> mean AND
    #                        variance match the beta kernel (exactly, when
    #                        cpg_sites == m0+1); use for comparisons
    # cpg_mode="mechanistic" classical maintain/de-novo rule below; its
    #                        mean map is affine and has a different fixed
    #                        point -- mechanistic studies only
    cpg_mode: str = "moment"
    cpg_p_maintain: float = 0.88    # marked site stays marked
    cpg_p_denovo: float = 0.02      # unmarked site becomes marked
    # deterministic kernel: u_daughter = phi(u_mother)


# ----------------------------------------------------------------------
# division (gene-expression partitioning)
# ----------------------------------------------------------------------
@dataclass
class DivisionConfig:
    partition: str = "exact"        # exact (paper) | binomial
    volume: float = 200.0           # molecules per AU at unit volume
    asymmetry: float = 0.5          # fraction of molecules to daughter 1
    # partitioning-noise strength (Soltani et al. 2016 Eq. 5):
    #   0 = deterministic, 1 = independent binomial, >1 = clustered
    alpha: float = 1.0


# ----------------------------------------------------------------------
# non-Markovian memory (P6)
# ----------------------------------------------------------------------
@dataclass
class MemoryConfig:
    # P6: non-Markovian memory via a per-cell leaky integrator (exponentially
    # weighted moving average) of a driver signal.  enabled=False keeps the
    # validated baseline: mem stays 0 and is never read by the rates.
    # mem_couple=0 tracks memory but does not affect kinetics, so the coupling
    # can be switched on independently of the tracking.  At mitosis the parent's
    # memory is copied to both daughters (state inheritance).
    enabled: bool = False
    tau: float = 10.0                # memory timescale [h] (EWMA window)
    driver: str = "stemness"         # signal integrated: "stemness" -> S in
                                     # [0,1]; "u_mean" -> mean epigenetic
                                     # state in [0,1].  Any other name raises
                                     # ValueError (see rates.memory_driver_signal)
    mem_couple: float = 0.0          # coupling of mem into beta (0 = no effect)


# ----------------------------------------------------------------------
# continuous cell kinetics (rate_mode = "continuous")
# ----------------------------------------------------------------------
@dataclass
class RatesConfig:
    rate_mode: str = "legacy"       # legacy | continuous
    # legacy parameters (Table 2)
    beta0: float = 0.04
    kappa0: float = 0.01
    mu0: float = 0.0002
    tau0: float = 4.0
    theta0: float = 200.0
    max_div_ta: int = 15
    x_threshold: float = 0.5        # legacy classification threshold
    # continuous-rate parameters (sharp limits reproduce legacy values)
    p_switch: float = 12.0          # Hill exponent of the sharp switch
    K_switch: float = None          # Hill threshold of the stemness switch;
                                    # None -> derived so that the continuous
                                    # boundary coincides with the legacy
                                    # x_threshold on the diagonal
                                    # (see rates.resolve_K)
    # contact inhibition (spatial; also usable without space via Q feedback)
    theta_rho: float = 12.0         # half-effective neighbour count
    rho0: float = 6.0
    # --- feedback time delay (P8: non-instantaneous niche sensing) ---
    # If > 0, the niche signal Q (soft SC count in continuous mode, SC count
    # in legacy mode) feeding beta_SC is sampled `delay_steps` steps in the
    # past instead of the current step. 0 = paper (instantaneous) behaviour;
    # default 0 keeps the validated baseline.
    delay_steps: int = 0               # neighbours before inhibition starts


# ----------------------------------------------------------------------
# 3D space
# ----------------------------------------------------------------------
@dataclass
class SpaceConfig:
    enabled: bool = False
    cell_radius: float = 1.0        # [mu-m, arbitrary]
    domain_radius: float = 60.0     # spherical domain, reflecting boundary
    relax: float = 0.15             # overlap relaxation fraction per step
                                    # at a nominal dt = 0.25 h (legacy)
    relax_time: float = None        # mechanical relaxation time [h]; None
                                    # -> derived from `relax` (see
                                    # Space.relax).  Set explicitly for
                                    # dt-independent mechanics.
    motility: float = 0.0           # Brownian motility [mu-m/h^0.5]
    contact_inhibition: bool = True
    inhibition_radius: float = 1.6  # in units of cell diameters
    # paracrine signalling field (population -> intracellular coupling)
    paracrine: bool = False
    diff_rate: float = 25.0         # [mu-m^2/h]
    decay_rate: float = 0.02        # [1/h]
    secretion: float = 0.05         # per stem-like cell per h
    response: tuple = (0.0, 0.0)    # activation added to gene i per unit field
    grid_h: float = 4.0             # field grid spacing [mu-m]
    normalize_source: bool = False  # divide injection by h^3 -> true
                                    # concentration units


# ----------------------------------------------------------------------
# top level
# ----------------------------------------------------------------------
@dataclass
class CERPSConfig:
    grn: GRNConfig = field(default_factory=GRNConfig)
    epi: EpiConfig = field(default_factory=EpiConfig)
    division: DivisionConfig = field(default_factory=DivisionConfig)
    rates: RatesConfig = field(default_factory=RatesConfig)
    space: SpaceConfig = field(default_factory=SpaceConfig)
    memory: MemoryConfig = field(default_factory=MemoryConfig)
    # protocol
    dt: float = 0.25
    n_init: int = 100
    max_time: float = 4000.0
    seed: int = 0
    # discard transitions before this time [h] when filling self.transitions
    # (paper Table-3 protocol: rec_tran_from = 2000 of a 4000 h run)
    rec_tran_from: float = 0.0
    # outputs
    record_stride: int = 4          # record counts every stride steps
    record_states: bool = True      # keep (x, u, pos) snapshots

    def with_(self, **kw):
        return replace(self, **kw)


def legacy_config() -> CERPSConfig:
    """CERPS configuration equivalent to the model of Huang et al. (2024):
    2-gene GRN, beta inheritance, exact expression partitioning, discrete
    Table-1 kinetics, no space."""
    return CERPSConfig()
