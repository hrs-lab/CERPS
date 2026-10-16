"""Fig. 4: population dynamics for four (alpha1, alpha2) parameter sets.
Fig. 5: scatter plots of (x1, x2) and (u1, u2) at t = 5, 50, 100, 1000 h
(taken from the sigma=0.05 base run, same parameters as Fig. 4A)."""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from common import DATA, FIGS, TABLES  # noqa: E402
from hsl2024 import style as St  # noqa: E402

SETS = [(0.4, 0.4), (0.4, 1.0), (1.0, 0.4), (1.0, 1.0)]
TITLES = {0.4: "0.4", 1.0: "1.0"}
SNAP_STEPS = [20, 200, 400, 4000]
SNAP_TIMES = [5, 50, 100, 1000]


def main_fig4():
    fig, axs = plt.subplots(1, 4, figsize=(10.5, 2.4))
    summary = []
    for ax, (a1, a2) in zip(axs, SETS):
        path = os.path.join(DATA, f"fig4_a{a1:g}_{a2:g}.npz")
        if not os.path.exists(path):
            print(f"skip fig4_a{a1:g}_{a2:g} (data missing)")
            continue
        z = np.load(path)
        ts = z["sys_ts"]                      # (steps, runs, 4)
        t = np.arange(ts.shape[0]) * 0.25
        r0 = ts[:, 0, :]                      # representative run
        ax.plot(t, r0[:, 1], color=St.C_SC, label="SC")
        ax.plot(t, r0[:, 2], color=St.C_TA1, label="TA1")
        ax.plot(t, r0[:, 3], color=St.C_TA2, label="TA2")
        ax.plot(t, ts[:, :, 0].mean(axis=1), "--", color=St.C_TOT,
                label="Total (avg)")
        ax.set_xlabel("h")
        ax.set_title(f"$\\alpha_1$={TITLES[a1]}, $\\alpha_2$={TITLES[a2]}",
                     loc="left", fontsize=8)
        ax.set_xlim(0, 2000)
        St.panel_label(ax, "ABCD"[SETS.index((a1, a2))])
        last = ts[-2000:].mean(axis=(0, 1))
        summary.append(((a1, a2), *last))
    axs[0].set_ylabel("Cell number")
    axs[0].legend(frameon=False, ncol=1, loc="upper left")
    fig.tight_layout()
    St.save(fig, "fig4_population_dynamics", FIGS)
    with open(os.path.join(TABLES, "fig4_homeostasis.csv"), "w") as f:
        f.write("alpha1,alpha2,total,SC,TA1,TA2\n")
        for (a1, a2), total, sc, ta1, ta2 in summary:
            f.write(f"{a1:g},{a2:g},{total:.1f},{sc:.1f},{ta1:.1f},"
                    f"{ta2:.1f}\n")


def main_fig5():
    path = os.path.join(DATA, "base.npz")
    if not os.path.exists(path):
        print("skip fig5 (base data missing)")
        return
    z = np.load(path)
    if "snap_data" not in z:
        print("skip fig5 (no snapshots)")
        return
    snap = z["snap_data"]                    # (n_snap, cap, 4)
    snap_n = z["snap_n"]
    fig, axs = plt.subplots(2, 4, figsize=(9.6, 4.6))
    for k, (ax_row, t_lab) in enumerate(zip(axs, ["A", "B"])):
        for j in range(4):
            ax = ax_row[j]
            n = int(snap_n[j])
            d = snap[j, :n]
            if k == 0:
                ax.scatter(d[:, 0], d[:, 1], s=2.5, c="#1f4fd8", alpha=0.35,
                           linewidths=0)
                ax.set_xlabel("$x_1$")
                ax.set_ylabel("$x_2$")
                ax.set_xlim(0, 3)
                ax.set_ylim(0, 3)
            else:
                ax.scatter(d[:, 2], d[:, 3], s=2.5, c="#b8860b", alpha=0.35,
                           linewidths=0)
                ax.set_xlabel("$u_1$")
                ax.set_ylabel("$u_2$")
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)
            ax.set_title(f"t = {SNAP_TIMES[j]} h", fontsize=8)
            if j == 0:
                St.panel_label(ax, t_lab)
    fig.tight_layout()
    St.save(fig, "fig5_scatter_evolution", FIGS)


if __name__ == "__main__":
    main_fig4()
    main_fig5()
