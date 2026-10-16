"""Fig. 7: extinction / SC-only population dynamics for three (sigma, m2)."""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from common import DATA, FIGS  # noqa: E402
from hsl2024 import style as St  # noqa: E402

CASES = [(0.20, 0.80), (0.05, 0.75), (0.05, 1.00)]


def main():
    fig, axs = plt.subplots(1, 3, figsize=(8.4, 2.4))
    for ax, (s, m), lab in zip(axs, CASES, "ABC"):
        path = os.path.join(DATA, f"fig7_s{s:g}_m{m:g}.npz")
        if not os.path.exists(path):
            print(f"skip fig7_s{s:g}_m{m:g}")
            continue
        z = np.load(path)
        ts = z["sys_ts"][:, 0, :]
        t = np.arange(ts.shape[0]) * 0.25
        ax.plot(t, ts[:, 1], color=St.C_SC, label="SC")
        ax.plot(t, ts[:, 2], color=St.C_TA1, label="TA1")
        ax.plot(t, ts[:, 3], color=St.C_TA2, label="TA2")
        ax.plot(t, ts[:, 0], "--", color=St.C_TOT, label="Total")
        ax.set_xlabel("h")
        ax.set_title(f"$\\sigma$={s:g}, $m_2$={m:.2f}", loc="left",
                     fontsize=8)
        ax.set_xlim(0, 4000)
        St.panel_label(ax, lab)
    axs[0].set_ylabel("Cell number")
    axs[0].legend(frameon=False)
    fig.tight_layout()
    St.save(fig, "fig7_extinction", FIGS)


if __name__ == "__main__":
    main()
