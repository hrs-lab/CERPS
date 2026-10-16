"""Example 1: legacy mode -- CERPS reproduces the model of Huang et al.

Runs the 2-gene paper configuration (beta inheritance, exact partition,
discrete Table-1 kinetics, no space) and plots the population dynamics
(cf. Fig. 4A of the paper).
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from cerps.config import CERPSConfig  # noqa: E402
from cerps.simulator import Simulator  # noqa: E402

OUT = os.path.join(os.path.dirname(HERE), "reproduction", "results",
                   "figures")


def main(t_max=2000.0, seed=1):
    cfg = CERPSConfig(max_time=t_max, seed=seed)
    sim = Simulator(cfg)
    sim.run()
    rec = sim.records
    t = np.array(rec["t"])

    fig, ax = plt.subplots(figsize=(4.2, 2.6))
    ax.plot(t, rec["sc"], color="#2e9e2e", label="SC")
    ax.plot(t, rec["ta1"], color="#1f4fd8", label="TA1")
    ax.plot(t, rec["ta2"], color="#31c6f2", label="TA2")
    ax.plot(t, rec["total"], "--", color="#555555", label="Total")
    ax.set_xlabel("t (h)")
    ax.set_ylabel("cell number")
    ax.set_title("CERPS legacy mode (= Huang et al. 2024 model)")
    ax.legend(frameon=False)
    fig.tight_layout()
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, "cerps_example_legacy.png"), dpi=200)
    print("final:", rec["total"][-1], "SC:", rec["sc"][-1],
          "TA1:", rec["ta1"][-1], "TA2:", rec["ta2"][-1])
    return rec


if __name__ == "__main__":
    main()
