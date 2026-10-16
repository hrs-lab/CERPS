"""Example: 3D spatial growth with contact inhibition and paracrine
signalling (requirements 4-6).

A stem-cell spheroid grows from a small cluster.  Proliferation is
continuously regulated by stemness AND by the local density (contact
inhibition); cells secrete a short-range paracrine factor that enhances
gene-1 activation in nearby cells (population -> intracellular coupling),
on top of the bidirectional epigenetics <-> GRN coupling.
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from cerps.config import (CERPSConfig, GRNConfig, EpiConfig,  # noqa: E402
                          RatesConfig, SpaceConfig, DivisionConfig)
from cerps.simulator import Simulator  # noqa: E402
from cerps import rates as rates_mod  # noqa: E402

OUT = os.path.join(os.path.dirname(HERE), "reproduction", "results",
                   "figures")


def main(t_max=600.0):
    cfg = CERPSConfig(
        grn=GRNConfig(),
        epi=EpiConfig(kernel="beta"),
        division=DivisionConfig(partition="binomial", volume=200.0),
        rates=RatesConfig(rate_mode="continuous", p_switch=12.0),
        space=SpaceConfig(
            enabled=True, cell_radius=1.0, domain_radius=45.0,
            relax=0.2, contact_inhibition=True, inhibition_radius=1.6,
            paracrine=True, response=(0.02, 0.0), grid_h=3.0,
            diff_rate=20.0, decay_rate=0.05, secretion=0.05),
        n_init=50, max_time=t_max, seed=5, record_stride=8)
    sim = Simulator(cfg)
    sim.run(t_max=t_max, snapshot_times=(t_max / 2, t_max))
    rec = sim.records
    t = np.array(rec["t"])

    fig = plt.figure(figsize=(9.6, 3.2))
    ax = fig.add_subplot(1, 3, 1)
    ax.plot(t, rec["sc"], color="#2e9e2e", label="SC-like")
    ax.plot(t, np.array(rec["ta1"]) + np.array(rec["ta2"]),
            color="#1f4fd8", label="TA")
    ax.plot(t, rec["total"], "--", color="#555555", label="Total")
    ax.set_xlabel("t (h)")
    ax.set_ylabel("cell number")
    ax.set_title("growth with contact inhibition")
    ax.legend(frameon=False, fontsize=7)

    snap = sim.snapshots[-1]
    pos = snap["pos"]
    x = snap["x"]
    S = rates_mod.stemness(sim.grn, x)
    ax = fig.add_subplot(1, 3, 2, projection="3d")
    sc_ = ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2], c=S, s=4,
                     cmap="viridis", vmin=0, vmax=1)
    ax.set_title(f"spheroid at t={snap['t']:.0f} h (colour: stemness)")
    plt.colorbar(sc_, ax=ax, shrink=0.6, label="stemness S")

    ax = fig.add_subplot(1, 3, 3)
    r = np.linalg.norm(pos, axis=1)
    ax.scatter(r, S, s=4, c="#1f4fd8", alpha=0.4)
    ax.set_xlabel("distance from centre r [$\\mu$m]")
    ax.set_ylabel("stemness S")
    ax.set_title("spatial stemness gradient")

    fig.tight_layout()
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, "cerps_example_3d.png"), dpi=200)
    print(f"final n={len(pos)}, mean S={S.mean():.3f}, max r={r.max():.1f}")
    # spatial correlation of stemness with radius
    if len(r) > 10:
        c = np.corrcoef(r, S)[0, 1]
        print(f"corr(radius, stemness) = {c:+.3f} "
              "(negative: stem cells concentrated in the core)")


if __name__ == "__main__":
    main()
