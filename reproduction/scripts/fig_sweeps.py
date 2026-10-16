"""Generic sweep figure (paper Figs. 8, 9, 12 layout, panels A-F).

A: homeostasis cell numbers (lines + error bars);  B: ratios (bars);
C: self-renewal probabilities; D: differentiation (SC->TA1/TA2);
E: dedifferentiation (TA->SC); F: transdifferentiation (TA1<->TA2).
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from common import DATA, FIGS, TABLES  # noqa: E402
from hsl2024.analysis import homeostasis_stats, transition_stats  # noqa: E402
from hsl2024 import style as St  # noqa: E402


def load_sweep(names, xvals, xkey="xval"):
    """Load batches; returns list of dicts with stats, or None entries."""
    out = []
    for name, xv in zip(names, xvals):
        path = os.path.join(DATA, name + ".npz")
        if not os.path.exists(path):
            out.append(None)
            continue
        z = np.load(path)
        res = dict(sys_ts=z["sys_ts"], tran=z["tran"],
                   stop_step=z["stop_step"])
        hs = homeostasis_stats(res)
        tstat = transition_stats(res)
        out.append(dict(x=xv, hs=hs, tr=tstat, z=z))
    return out


def sweep_figure(entries, xlabel, fname, title_prefix, ylims=None):
    """entries: list of dicts from load_sweep (no None)."""
    xs = [e["x"] for e in entries]
    cnt_m = np.array([e["hs"]["count_mean"] for e in entries])
    cnt_s = np.array([e["hs"]["count_std"] for e in entries])
    rat_m = np.array([e["hs"]["ratio_mean"] for e in entries])
    rat_s = np.array([e["hs"]["ratio_std"] for e in entries])
    pm = np.array([e["tr"]["prob_mean"] for e in entries])   # (n, 3, 3)
    ps = np.array([e["tr"]["prob_std"] for e in entries])

    fig, axs = plt.subplots(3, 2, figsize=(6.4, 7.2))
    axA, axB, axC, axD, axE, axF = axs.ravel()

    # A: numbers
    for j, (lab, col) in enumerate(zip(
            ["SC", "TA1", "TA2", "Total"],
            [St.C_SC, St.C_TA1, St.C_TA2, St.C_TOT])):
        col_i = [1, 2, 3, 0][j]
        axA.errorbar(xs, cnt_m[:, col_i], yerr=cnt_s[:, col_i], color=col,
                     marker="o", ms=3, capsize=2, label=lab)
    axA.set_ylabel("Cell number")
    axA.set_xlabel(xlabel)
    axA.legend(frameon=False, ncol=2)
    St.panel_label(axA, "A")

    # B: ratios
    w = 0.8 / 3
    xind = np.arange(len(xs))
    for j, (lab, col) in enumerate(zip(["SC", "TA1", "TA2"],
                                       [St.C_SC, St.C_TA1, St.C_TA2])):
        axB.bar(xind + (j - 1) * w, rat_m[:, j], w, yerr=rat_s[:, j],
                color=col, capsize=2, label=lab)
    axB.set_xticks(xind)
    axB.set_xticklabels([f"{x:g}" for x in xs])
    axB.set_xlabel(xlabel)
    axB.set_ylabel("Cell ratio (%)")
    axB.legend(frameon=False, ncol=3)
    St.panel_label(axB, "B")

    def bars(ax, idx_pairs, labels, title):
        for j, (i, k) in enumerate(idx_pairs):
            ax.bar(xind + (j - 0.5) * (0.8 / len(idx_pairs)),
                   pm[:, i, k], 0.8 / len(idx_pairs),
                   yerr=ps[:, i, k], color=St.COLS[j], capsize=2,
                   label=labels[j])
        ax.set_xticks(xind)
        ax.set_xticklabels([f"{x:g}" for x in xs])
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Transition probability (%)")
        ax.legend(frameon=False, ncol=2)
        St.panel_label(ax, title)

    bars(axC, [(0, 0), (1, 1), (2, 2)], ["SC-SC", "TA1-TA1", "TA2-TA2"], "C")
    bars(axD, [(0, 1), (0, 2)], ["SC-TA1", "SC-TA2"], "D")
    bars(axE, [(1, 0), (2, 0)], ["TA1-SC", "TA2-SC"], "E")
    bars(axF, [(1, 2), (2, 1)], ["TA1-TA2", "TA2-TA1"], "F")

    if ylims:
        for ax, ylim in zip([axA, axB, axC, axD, axE, axF], ylims):
            if ylim:
                ax.set_ylim(*ylim)
    fig.suptitle(title_prefix, y=0.995, fontsize=9)
    fig.tight_layout()
    St.save(fig, fname, FIGS)

    # CSV dump
    csv = os.path.join(TABLES, fname + ".csv")
    with open(csv, "w") as f:
        f.write(f"{xlabel},total,SC,TA1,TA2,SC_pct,TA1_pct,TA2_pct,"
                "SC_SC,SC_TA1,SC_TA2,TA1_SC,TA1_TA1,TA1_TA2,"
                "TA2_SC,TA2_TA1,TA2_TA2\n")
        for i, e in enumerate(entries):
            row = [e["x"], *cnt_m[i], *rat_m[i], *pm[i].ravel()]
            f.write(",".join(f"{v:.5g}" for v in row) + "\n")
    return dict(count_mean=cnt_m, ratio_mean=rat_m, prob_mean=pm)


def main_fig8():
    names = ["fig8_s0", "fig8_s0.025", "base", "fig8_s0.075", "fig8_s0.1"]
    xvals = [0.0, 0.025, 0.05, 0.075, 0.1]
    entries = load_sweep(names, xvals)
    missing = [n for n, e in zip(names, entries) if e is None]
    if missing:
        print("fig8 missing:", missing)
    ok = [e for e in entries if e]
    if len(ok) == len(names):
        sweep_figure(entries, r"$\sigma$", "fig8_noise_sweep",
                     "Impact of extrinsic noise (Fig. 8)",
                     ylims=[(0, None), (0, 85), (0, 160), (0, 20),
                            (0, 0.008), (0, 0.02)])


def main_fig9():
    names = ["base", "fig9_m0.825", "fig9_m0.85", "fig9_m0.875", "fig9_m0.9"]
    xvals = [0.800, 0.825, 0.850, 0.875, 0.900]
    entries = load_sweep(names, xvals)
    missing = [n for n, e in zip(names, entries) if e is None]
    if missing:
        print("fig9 missing:", missing)
    if all(entries):
        sweep_figure(entries, r"$m_2$", "fig9_inheritance_sweep",
                     "Impact of epigenetic inheritance (Fig. 9)",
                     ylims=[(0, None), (0, 85), (0, 160), (0, 20),
                            (0, 2.5), (0, 0.5)])


def main_fig12():
    names = ["base", "fig12_tf0.05", "fig12_tf0.1", "fig12_tf0.15",
             "fig12_tf0.2"]
    xvals = [0.0, 0.05, 0.1, 0.15, 0.2]
    entries = load_sweep(names, xvals)
    missing = [n for n, e in zip(names, entries) if e is None]
    if missing:
        print("fig12 missing:", missing)
    if all(entries):
        sweep_figure(entries, r"$a_{TF}$", "fig12_tf_effect",
                     "Effect of transcription factor (Fig. 12)",
                     ylims=[(0, None), (0, 85), (0, 160), (0, 20),
                            (0, 0.008), (0, 0.02)])


if __name__ == "__main__":
    main_fig8()
    main_fig9()
    main_fig12()
