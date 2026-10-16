"""Fig. 6: dedifferentiation / transdifferentiation trajectories.

Uses the numpy reference simulator with lineage tracking (candidate cells
with inherited u above threshold) to extract (x1, x2) trajectories of
de-/trans-differentiation events, plus insets of x1(t), x2(t).
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from hsl2024 import DEFAULT  # noqa: E402
from hsl2024.reference_sim import simulate_reference  # noqa: E402
from hsl2024 import style as St  # noqa: E402
from common import FIGS, TABLES, DATA  # noqa: E402

PANELS = [("TA1-SC", 1, 0), ("TA2-SC", 3, 0),
          ("TA1-TA2", 1, 3), ("TA2-TA1", 3, 1)]


def collect_events(n_steps=8000, seed=3, t_min=500.0):
    """Run the reference simulator once and return de-/trans-diff events."""
    out_cache = os.path.join(DATA, "fig6_events.npz")
    if os.path.exists(out_cache):
        z = np.load(out_cache, allow_pickle=True)
        return [e for e in z["events"]]

    p = DEFAULT
    res = simulate_reference(p, n_runs=1, n_steps=n_steps, seed=seed,
                             rec_tran_from=10 ** 9, u_track=0.36)
    events = [e for e in res["events"] if e["t"] >= t_min]
    np.savez_compressed(out_cache, events=np.array(events, dtype=object))
    return events


def pick_event(events, old, new):
    """Pick a representative (clean, single-step) event of a given type."""
    cands = [e for e in events if e["old"] == old and e["new"] == new]
    if not cands:
        return None
    # prefer events with moderate duration (50-400 h between birth and the
    # type change) and a long recorded history
    def score(e):
        h = e["history"]
        m = h[:, 0] > 0
        dur = e["t"] - h[m][0, 0] if m.any() else 1e9
        return -abs(np.log(max(dur, 1.0) / 150.0))
    cands.sort(key=score, reverse=True)
    return cands[0]


def main():
    events = collect_events()
    print(f"{len(events)} tracked transition events")
    fig, axs = plt.subplots(1, 4, figsize=(10.5, 2.7))
    got = []
    for ax, (lab, old, new) in zip(axs, PANELS):
        e = pick_event(events, old, new)
        St.panel_label(ax, "ABCD"[PANELS.index((lab, old, new))])
        if e is None:
            ax.text(0.5, 0.5, "no event captured", ha="center", va="center",
                    transform=ax.transAxes)
            got.append(None)
            continue
        h = e["history"]
        h = h[h[:, 0] > 0]
        sc = ax.scatter(h[:, 1], h[:, 2], c=h[:, 0], s=3, cmap="viridis",
                        linewidths=0)
        ax.set_xlabel("$x_1$")
        ax.set_ylabel("$x_2$")
        ax.set_xlim(0, 3)
        ax.set_ylim(0, 3)
        ax.set_title(lab, fontsize=8)
        # inset: x1(t), x2(t)
        ins = ax.inset_axes([0.58, 0.08, 0.38, 0.34])
        ins.plot(h[:, 0], h[:, 1], color=St.C_TA1, lw=0.7)
        ins.plot(h[:, 0], h[:, 2], color=St.C_TA2, lw=0.7)
        ins.axhline(0.5, color="0.7", lw=0.5, ls=":")
        ins.set_xticks([])
        ins.set_yticks([0, 1, 2])
        ins.tick_params(labelsize=5)
        cb = plt.colorbar(sc, ax=ax, fraction=0.045)
        cb.set_label("Time (h)", fontsize=6)
        cb.ax.tick_params(labelsize=6)
        got.append((lab, e["t"], len(h)))
    fig.tight_layout()
    St.save(fig, "fig6_transition_trajectories", FIGS)
    with open(os.path.join(TABLES, "fig6_events_summary.csv"), "w") as f:
        f.write("panel,event_time_h,history_len\n")
        for g in got:
            if g:
                f.write(f"{g[0]},{g[1]:.1f},{g[2]}\n")


if __name__ == "__main__":
    main()
