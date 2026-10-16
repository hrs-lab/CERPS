"""Fig. 2: bifurcation analysis (A), phenotype map (B), sample dynamics (C)
and average SC duration (D)."""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from hsl2024 import DEFAULT  # noqa: E402
from hsl2024.steady_state import (  # noqa: E402
    find_steady_states, a_of_u, classify_state)
from hsl2024.single_cell import sc_durations, single_trajectory  # noqa: E402
from hsl2024 import style as St  # noqa: E402
from common import FIGS, TABLES  # noqa: E402

p = DEFAULT


def panel_a(ax):
    a_grid = np.linspace(0.0, 2.2, 221)
    rows = []
    for a in a_grid:
        pts, stable = find_steady_states(a, a, p)
        for x, st in zip(pts, stable):
            rows.append((a, x[0], st))
    rows = np.array(rows)
    for st_val, ls, col in [(True, "-", "tab:grey"), (False, "--", "tab:red")]:
        m = rows[:, 2] == st_val
        ax.plot(rows[m, 0], rows[m, 1], ls, color=col, lw=1.0,
                solid_capstyle="round")
    ax.axvspan(0.75, 1.65, color="0.9", zorder=0)
    ax.axvline(0.75, color="0.6", lw=0.7, ls=":")
    ax.axvline(1.65, color="0.6", lw=0.7, ls=":")
    ax.set_xlabel("a")
    ax.set_ylabel(r"$x^*$")
    ax.set_xlim(0, 2.2)
    ax.set_ylim(0, 3)
    ax.legend(handles=[
        Line2D([], [], color="tab:grey", ls="-", label="stable"),
        Line2D([], [], color="tab:red", ls="--", label="unstable")],
        frameon=False, loc="upper left")
    St.panel_label(ax, "A")


def panel_b(ax):
    n = 41
    us = np.linspace(0.02, 0.98, n)
    count = np.zeros((n, n), int)        # number of stable states
    typeset = [[[] for _ in range(n)] for _ in range(n)]
    for i, u1 in enumerate(us):
        for j, u2 in enumerate(us):
            a1, a2 = a_of_u((u1, u2), p)
            pts, stable = find_steady_states(a1, a2, p)
            ss = [x for x, st in zip(pts, stable) if st]
            count[j, i] = len(ss)
            typeset[j][i] = sorted({classify_state(x) for x in ss} - {-1})
    cmap = ListedColormap(["#d0d0d0", "#7fbf7f", "#2f7f2f", "#1a4a1a"])
    ax.imshow(count, origin="lower", extent=(0, 1, 0, 1), cmap=cmap,
              vmin=0.5, vmax=3.5, interpolation="nearest")
    # region labels
    label_pos = {
        "TA2": (0.2, 0.85), "SC/TA2": (0.45, 0.85), "SC": (0.8, 0.8),
        "SC/TA1": (0.82, 0.45), "TA1": (0.85, 0.2), "TA1/TA2": (0.35, 0.2),
        "SC/TA1/TA2": (0.55, 0.55),
    }
    for lab, (xx, yy) in label_pos.items():
        ax.text(xx, yy, lab, ha="center", va="center", fontsize=7)
    ax.plot(0.34, 0.34, "o", mfc="#ffd700", mec="k", ms=6, zorder=5)
    ax.set_xlabel(r"$u_1$")
    ax.set_ylabel(r"$u_2$")
    St.panel_label(ax, "B")
    return count


def panel_c(ax):
    ts, xs = single_trajectory((0.34, 0.34), (0.8, 0.8), p, t_max=100.0,
                               seed=7)
    ax.plot(ts, xs[:, 0], color=St.C_TA1, label=r"$x_1$")
    ax.plot(ts, xs[:, 1], color=St.C_TA2, label=r"$x_2$")
    ax.set_xlabel("h")
    ax.set_ylabel(r"$x^*$")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 2)
    ax.legend(frameon=False)
    ax.text(0.97, 0.95, r"$u_1=u_2=0.34$" + "\n" + r"$\sigma=0.05$",
            transform=ax.transAxes, ha="right", va="top")
    St.panel_label(ax, "C")


def panel_d(ax):
    n = 21
    us = np.linspace(0.025, 0.975, n)
    U1, U2 = np.meshgrid(us, us, indexing="ij")
    u_grid = np.stack([U1.ravel(), U2.ravel()], axis=1)
    # initial condition: SC steady state where it exists, else the closest
    from hsl2024.steady_state import sc_steady_state
    x0 = np.zeros_like(u_grid)
    for k in range(u_grid.shape[0]):
        xs = sc_steady_state(u_grid[k], p)
        x0[k] = xs if xs is not None else (0.3, 0.3)
    dur = sc_durations(u_grid, x0, p, t_max=100.0, n_reps=100, seed=11)
    D = dur.reshape(n, n).T
    im = ax.imshow(D, origin="lower", extent=(0, 1, 0, 1), cmap="viridis",
                   vmin=10, vmax=100, interpolation="nearest")
    ax.plot(0.34, 0.34, "o", mfc="#ffd700", mec="k", ms=6, zorder=5)
    ax.set_xlabel(r"$u_1$")
    ax.set_ylabel(r"$u_2$")
    cb = plt.colorbar(im, ax=ax, fraction=0.045)
    cb.set_label("Average duration (h)")
    St.panel_label(ax, "D")
    np.save(os.path.join(TABLES, "fig2D_duration_grid.npy"), D)
    return D


def main():
    fig, axs = plt.subplots(1, 4, figsize=(10.5, 2.6))
    panel_a(axs[0])
    panel_b(axs[1])
    panel_c(axs[2])
    panel_d(axs[3])
    fig.tight_layout()
    St.save(fig, "fig2_bifurcation", FIGS)


if __name__ == "__main__":
    main()
