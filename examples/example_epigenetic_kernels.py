"""Example: comparison of epigenetic inheritance kernels (requirement 1).

Simulates a single TA1 lineage over many divisions with each inheritance
kernel (beta / gaussian / cpg / deterministic) and shows the distribution
of the epigenetic state u across descendants, starting from a stem-like
state u = 0.5.
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from cerps.config import EpiConfig  # noqa: E402
from cerps import epigenetics as epi  # noqa: E402

OUT = os.path.join(os.path.dirname(HERE), "reproduction", "results",
                   "figures")


def lineage(kernel, n_div=60, n_lin=2000, seed=0):
    """Each lineage starts at u=0.5 and divides n_div times (no dynamics
    between divisions); returns the distribution of u in the final
    generation."""
    rng = np.random.default_rng(seed)
    cfg = EpiConfig(kernel=kernel)
    u = np.full((n_lin, 1), 0.5)
    for _ in range(n_div):
        u, _ = epi.inherit(u, cfg, rng)      # track one daughter branch
        u = np.clip(u, 0.0, 1.0)
    return u[:, 0]


def main():
    fig, ax = plt.subplots(figsize=(4.6, 2.8))
    for kernel, col in [("beta", "#1f4fd8"), ("gaussian", "#2e9e2e"),
                        ("cpg", "#b8860b"), ("deterministic", "#d62728")]:
        u = lineage(kernel)
        ax.hist(u, bins=40, range=(0, 0.5), histtype="step", density=True,
                color=col, label=kernel, lw=1.2)
        print(f"{kernel:13s} mean={u.mean():.4f} sd={u.std():.4f}")
    ax.set_xlabel("epigenetic state $u$ after 60 divisions")
    ax.set_ylabel("density")
    ax.set_title("Division inheritance kernels (lineages from $u=0.5$)")
    ax.legend(frameon=False)
    fig.tight_layout()
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, "cerps_example_kernels.png"), dpi=200)


if __name__ == "__main__":
    main()
