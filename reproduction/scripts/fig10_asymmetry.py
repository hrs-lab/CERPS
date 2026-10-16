"""Fig. 10: transition probabilities with asymmetric parameters.

Grouped bars over m2 for the default and five single-parameter increases
(k1=1.1, alpha1=0.44, rho1=0.11, b1=1.1, s1=0.55), panels A-D as in the
paper (A: SC-TA1, B: SC-TA2, C: TA1-SC, D: TA2-SC).
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from common import DATA, FIGS, TABLES, M2S_FIG9, VARIANTS_FIG10  # noqa: E402
from hsl2024.analysis import transition_stats  # noqa: E402
from hsl2024 import style as St  # noqa: E402

SERIES = ["default"] + list(VARIANTS_FIG10.keys())
LABELS = {"default": "Default", "k1": r"$k_1=1.1$", "alpha1": r"$\alpha_1=0.44$",
          "rho1": r"$\rho_1=0.11$", "b1": r"$b_1=1.1$", "s1": r"$s_1=0.55$"}
SER_COLORS = ["#888888", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
              "#ff7f00"]


def main():
    # series -> per-m2 probability matrices
    data = {}
    for s in SERIES:
        row = []
        for m in M2S_FIG9:
            if s == "default":
                name = "base" if abs(m - 0.8) < 1e-9 else f"fig9_m{m:g}"
            else:
                name = f"fig10_{s}_m{m:g}"
            path = os.path.join(DATA, name + ".npz")
            if not os.path.exists(path):
                row.append(None)
                continue
            z = np.load(path)
            t = transition_stats(dict(tran=z["tran"]))
            row.append(t["prob_mean"])
        data[s] = row

    n_m = len(M2S_FIG9)
    if any(v is None for s in SERIES for v in data[s]):
        missing = sum(v is None for s in SERIES for v in data[s])
        print(f"fig10: {missing} batches missing, drawing partial")

    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.4))
    panels = [("A", 0, 1, "SC-TA1"), ("B", 0, 2, "SC-TA2"),
              ("C", 1, 0, "TA1-SC"), ("D", 2, 0, "TA2-SC")]
    xind = np.arange(n_m)
    w = 0.8 / len(SERIES)
    csv_rows = []
    for ax, (lab, i, k, name) in zip(axs.ravel(), panels):
        for j, s in enumerate(SERIES):
            vals, errs = [], []
            for v in data[s]:
                if v is None:
                    vals.append(np.nan); errs.append(np.nan)
                else:
                    vals.append(v[i, k])
            ax.bar(xind + (j - len(SERIES) / 2 + 0.5) * w, vals, w,
                   color=SER_COLORS[j], label=LABELS[s])
            if s == "default":
                csv_rows.append((name, *vals))
        ax.set_xticks(xind)
        ax.set_xticklabels([f"{m:.3f}" for m in M2S_FIG9])
        ax.set_xlabel(r"$m_2$")
        ax.set_ylabel("Transition probability (%)")
        ax.set_title(name, fontsize=8, loc="left")
        St.panel_label(ax, lab)
    axs[0, 0].legend(frameon=False, ncol=2, fontsize=6)
    fig.tight_layout()
    St.save(fig, "fig10_asymmetry", FIGS)
    with open(os.path.join(TABLES, "fig10_asymmetry.csv"), "w") as f:
        f.write("series," + ",".join(f"m2={m:g}" for m in M2S_FIG9) + "\n")
        for row in csv_rows:
            f.write(row[0] + "," + ",".join(f"{v:.4g}" for v in row[1:]) + "\n")


if __name__ == "__main__":
    main()
