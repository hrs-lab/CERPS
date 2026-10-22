"""CERPS: Coupled Epigenetic-Regulatory-Population-Spatial model.

A general cell-differentiation simulation framework extending the hybrid
model of Huang, Situ & Lei (2024).  See docs/model_spec.md for the full
mathematical specification and docs/literature_review.md for the design
rationale.

Modules
-------
config      all parameters of a model (dataclasses)
notation    x=(u,v) / (x,u) naming alias layer isolating the two symbol
            conventions (HSL2024 vs the NSFC proposal)
grn         N-gene regulation network
h_general   the proposal's F(x) = (G(x), H(x)) split as a switchable
            object: G with the epigenetic gates alpha_i(v), beta_i(v),
            k_i(v); H with the zero / linear / general levels
epigenetics continuous mark dynamics + division inheritance kernels
division    gene-expression partitioning at mitosis
rates       cell kinetics (legacy / continuous)
spatial     3-D mechanics, contact inhibition, paracrine field
simulator   individual-based event loop coupling everything
batch       parameter sweeps, multi-seed replicates, checkpointing,
            statistics and numerical-health checks (the experiment layer)
"""
__version__ = "1.2.0"

from . import config, division, epigenetics, grn, rates, spatial  # noqa: F401
from . import notation, h_general                                    # noqa: F401
from . import batch  # noqa: F401
