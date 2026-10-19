"""CERPS: Coupled Epigenetic-Regulatory-Population-Spatial model.

A general cell-differentiation simulation framework extending the hybrid
model of Huang, Situ & Lei (2024).  See docs/model_spec.md for the full
mathematical specification and docs/literature_review.md for the design
rationale.

Modules
-------
config      all parameters of a model (dataclasses)
grn         N-gene regulation network
epigenetics continuous mark dynamics + division inheritance kernels
division    gene-expression partitioning at mitosis
rates       cell kinetics (legacy / continuous)
spatial     3-D mechanics, contact inhibition, paracrine field
simulator   individual-based event loop coupling everything
batch       parameter sweeps, multi-seed replicates, checkpointing,
            statistics and numerical-health checks (the experiment layer)
"""
__version__ = "1.1.0"

from . import config, division, epigenetics, grn, rates, spatial  # noqa: F401
from . import batch  # noqa: F401
