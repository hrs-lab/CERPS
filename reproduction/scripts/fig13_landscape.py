"""Fig. 13: temporal evolution of the Waddington landscape U(t, z),
z = x1 - x2, for varying sigma (row A), m2 (row B) and a_TF (row C)."""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from common import DATA, FIGS, SIGMAS_FIG8, M2S_FIG9, ATF_FIG12  # noqa: E402
from hsl2024.analysis import landscape_from_hist  # noqa: E402
from hsl2024 import style as St  # noqa: E402

ROWS = [
    ("A", r"$\sigma$", [("fig8_s0", 0.0), ("fig8_s0.025", 0.025),
                        ("base", 0.05), ("fig8_s0.075", 0.075),
                        ("fig8_s0.1", 0.1)]),
    ("B", r"$m_2$", [("base", 0.800), ("fig9_m0.825", 0.825),
                     ("fig9_m0.85", 0.85), ("fig9_m0.875", 0.875),
                     ("fig9_m0.9", 0.900)]),
    ("C", r"$a_{TF}$", [("base", 0.0), ("fig12_tf0.05", 0.05),
                        ("fig12_tf0.1", 0.1), ("fig12_tf0.15", 0.15),
                        ("fig12_tf0.2", 0.2)]),
]

TB_STRIDE = 40          # 10 h per time bin
N_TB = 100              # covers the first 1000 h


def main():
    fig, axs = plt.subplots(3, 5, figsize=(10.5, 6.4))
    for r, (rlab, sym, jobs) in enumerate(ROWS):
        for c, (name, xv) in enumerate(jobs):
            ax = axs[r, c]
            path = os.path.join(DATA, name + ".npz")
            if not os.path.exists(path) or "hist" not in np.load(path):
                ax.text(0.5, 0.5, "data missing", ha="center", va="center",
                        transform=ax.transAxes)
                continue
            z = np.load(path)
            hist = z["hist"]
            res = dict(hist=hist, tb_stride=TB_STRIDE)
            U = landscape_from_hist(res)
            im = ax.imshow(U.T, origin="lower", aspect="auto",
                           extent=[0, TB_STRIDE * N_TB * 0.25 / 1000 * 1000,
                                   -4, 4], cmap="viridis_r",
                           vmin=1.5, vmax=4.0, interpolation="nearest")
            ax.set_title(f"${sym.replace('$','')}={xv:g}$", fontsize=8)
            if c == 0:
                ax.set_ylabel("z")
                ax.text(-0.42, 0.5, rlab, transform=ax.transAxes,
                        fontsize=12, fontweight="bold", va="center")
            if r == 2:
                ax.set_xlabel("t (h)")
            ax.set_xlim(0, 1000)
    cb = fig.colorbar(im, ax=axs, fraction=0.02, pad=0.01)
    cb.set_label(r"$U(t,z)=-\log f(t,z)$")
    St.save(fig, "fig13_waddington_landscape", FIGS)


if __name__ == "__main__":
    main()
