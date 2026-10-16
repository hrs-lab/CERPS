"""Plot style helpers (paper-like appearance)."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# paper-like colours
C_SC = "#2e9e2e"      # green
C_TA1 = "#1f4fd8"     # dark blue
C_TA2 = "#31c6f2"     # light blue
C_TOT = "#555555"     # grey/black
COLS = [C_SC, C_TA1, C_TA2, C_TOT]
TYPE_NAMES = ["SC", "TA1", "TA2"]

plt.rcParams.update({
    "font.size": 8,
    "axes.linewidth": 0.7,
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "lines.linewidth": 0.9,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})


def panel_label(ax, s):
    ax.set_title(s, loc="left", fontweight="bold", fontsize=10)


def save(fig, name, figs_dir):
    fig.savefig(f"{figs_dir}/{name}.png")
    fig.savefig(f"{figs_dir}/{name}.pdf")
    plt.close(fig)
    print(f"saved {name}.png/.pdf")
