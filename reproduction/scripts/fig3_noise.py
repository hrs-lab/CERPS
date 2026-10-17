"""Fig. 3: noise-induced differentiation with fixed epigenetic states."""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from common import DATA, FIGS, TABLES, US_FIG3, SIGMAS_FIG3  # noqa: E402
from hsl2024 import style as St  # noqa: E402


def load(u_pair, sigma):
    name = f"fig3_u{u_pair[0]:g}_{u_pair[1]:g}_s{sigma:g}"
    path = os.path.join(DATA, name + ".npz")
    if not os.path.exists(path):
        return None
    return np.load(path)


def main():
    fig, axs = plt.subplots(1, 4, figsize=(9.5, 2.3), sharey=True)
    table = []
    for k, (ax, u_pair) in enumerate(zip(axs, US_FIG3)):
        sc, ta1, ta2 = [], [], []
        for s in SIGMAS_FIG3:
            z = load(u_pair, s)
            if z is None:
                sc.append(np.nan); ta1.append(np.nan); ta2.append(np.nan)
                continue
            r = z["ratios"]                           # (SC, TA1, TA2) %
            sc.append(float(r[0])); ta1.append(float(r[1]))
            ta2.append(float(r[2]))
            table.append((u_pair[0], u_pair[1], s,
                          float(r[0]), float(r[1]), float(r[2])))
        ax.stackplot(SIGMAS_FIG3, sc, ta1, ta2,
                     colors=[St.C_SC, St.C_TA1, St.C_TA2], alpha=0.9)
        ax.set_title(f"$u_1$={u_pair[0]:g}$, $u_2$={u_pair[1]:g}$",
                     loc="left", fontsize=8)
        ax.set_xlabel(r"$\sigma$")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 100)
        St.panel_label(ax, "ABCD"[k])
    axs[0].set_ylabel("Cell ratio (%)")
    axs[0].legend(["SC", "TA1", "TA2"], frameon=False, loc="upper right",
                  ncol=1)
    fig.tight_layout()
    St.save(fig, "fig3_noise_induced", FIGS)
    with open(os.path.join(TABLES, "fig3_ratios.csv"), "w") as f:
        f.write("u1,u2,sigma,SC_pct,TA1_pct,TA2_pct\n")
        for row in table:
            f.write(",".join(f"{v:g}" for v in row) + "\n")


if __name__ == "__main__":
    main()
