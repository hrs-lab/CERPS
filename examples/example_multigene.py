"""Example: a 3-gene regulatory network (requirement 3).

Genes X1, X2, X3 repress each other (all-to-all); the epigenetic
inheritance kernel is Gaussian and the cell kinetics are continuous.

IMPORTANT - empirically measured attractor distribution
(see docs/multigene_attractors.md, 3-gene all-repression sweep):
the population does NOT settle into a symmetric four-phenotype
homeostasis. Single-gene-high states dominate (~93% across seeds),
a stable double-high tail persists (~7%), and the all-on (stem-like)
state is rare (~0.5-1%). The "label 0 (all-high by label)" line below
tracks `label == 0`, which in continuous mode freezes cells as
all-high at birth and therefore OVER-COUNTS the true all-on fraction;
the right panel's "3 high (true all-on)" bar shows the actual count of
cells with all three genes x > 0.5. Genuine steady state also needs
>= 2000-4000 h, far longer than the 1500 h default.

Gene-to-label map: 0 = all-high (by label), j+1 = gene j high.
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from cerps.config import (CERPSConfig, GRNConfig, EpiConfig,  # noqa: E402
                          RatesConfig)
from cerps.simulator import Simulator  # noqa: E402
from cerps import rates as rates_mod  # noqa: E402

OUT = os.path.join(os.path.dirname(HERE), "reproduction", "results",
                   "figures")


def main(t_max=1500.0):
    cfg = CERPSConfig(
        grn=GRNConfig(n_genes=3, alpha=(0.4, 0.4, 0.4),
                      sigma=(0.05, 0.05, 0.05)),
        epi=EpiConfig(kernel="gaussian"),
        rates=RatesConfig(rate_mode="continuous"),
        n_init=100, max_time=t_max, seed=3)
    sim = Simulator(cfg)
    sim.run()
    rec = sim.records
    t = np.array(rec["t"])

    # phenotype composition from the final snapshot
    snap = sim.snapshots[-1] if sim.snapshots else None
    fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.0))
    ax = axs[0]
    ax.plot(t, rec["sc"], color="#2e9e2e", label="label 0 (all-high by label)")
    ax.plot(t, rec["ta1"], color="#1f4fd8", label="X1 high")
    # NOTE: for G > 2 the recorder keeps the paper's field names, but the
    # label codes are 0 = all-high, j+1 = gene j high.  Hence "ta2" holds
    # label 3 = X3-high and "other" holds label 2 = X2-high.
    # CAUTION: in continuous mode `label == 0` is frozen at birth and
    # over-counts the true all-on state (see docs/multigene_attractors.md).
    ax.plot(t, rec["other"], color="#31c6f2", label="X2 high")
    ax.plot(t, rec["ta2"], color="#d62728", label="X3 high")
    ax.plot(t, rec["total"], "--", color="#555555", label="Total")
    ax.set_xlabel("t (h)")
    ax.set_ylabel("cell number")
    ax.legend(frameon=False, fontsize=6)
    ax.set_title("3-gene GRN population dynamics")
    ax = axs[1]
    if snap is not None:
        x = snap["x"]
        high = x > 0.5
        n_hi = high.sum(axis=1)
        cats = ["all-low", "1 high", "2 high", "3 high (true all-on)"]
        counts = [np.sum(n_hi == j) for j in range(4)]
        ax.bar(range(4), counts, color=["#999", "#1f4fd8", "#888",
                                        "#2e9e2e"])
        ax.set_xticks(range(4))
        ax.set_xticklabels(cats, fontsize=7)
        ax.set_ylabel("cells at t = %.0f h" % snap["t"])
        ax.set_title("final expression states")
    fig.tight_layout()
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, "cerps_example_3gene.png"), dpi=200)
    print("final total:", rec["total"][-1],
          " label-0 (all-high-by-label):", rec["sc"][-1],
          " (NOTE: over-counts true all-on; see docs/multigene_attractors.md)")


if __name__ == "__main__":
    main()
