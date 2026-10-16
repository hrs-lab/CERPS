"""Fig. 11: dynamics of cell reprogramming through TF induction."""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from common import DATA, FIGS, TABLES, ATF_FIG11  # noqa: E402
from hsl2024 import style as St  # noqa: E402


def load(which, a):
    path = os.path.join(DATA, f"fig11_{which}_tf{a:g}.npz")
    if not os.path.exists(path):
        return None
    return np.load(path)


def main():
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(6.8, 2.6))
    probs = {"TA1": [], "TA2": []}
    times = {"TA1": [], "TA2": []}
    terr = {"TA1": [], "TA2": []}
    rows = []
    for which in ("TA1", "TA2"):
        for a in ATF_FIG11:
            z = load(which, a)
            if z is None:
                probs[which].append(np.nan)
                times[which].append(np.nan)
                terr[which].append(np.nan)
                continue
            stop_kind = z["stop_kind"]
            stop_step = z["stop_step"]
            n_runs = len(stop_kind)
            succ = stop_kind == 2
            probs[which].append(100.0 * succ.mean())
            t_succ = stop_step[succ] * 0.25
            times[which].append(t_succ.mean() if succ.any() else np.nan)
            terr[which].append(t_succ.std() if succ.sum() > 1 else 0.0)
            rows.append((which, a, succ.mean() * 100,
                         t_succ.mean() if succ.any() else np.nan,
                         int(succ.sum()), n_runs))
    axA.plot(ATF_FIG11, probs["TA1"], "o-", color=St.C_TA1,
             label="100% TA1")
    axA.plot(ATF_FIG11, probs["TA2"], "s-", color=St.C_TA2,
             label="100% TA2")
    axA.set_xlabel(r"$a_{TF}$")
    axA.set_ylabel("Successful induction probability (%)")
    axA.set_ylim(-5, 105)
    axA.legend(frameon=False)
    St.panel_label(axA, "A")

    axB.errorbar(ATF_FIG11, times["TA1"], yerr=terr["TA1"], fmt="o-",
                 color=St.C_TA1, capsize=2, label="100% TA1")
    axB.errorbar(ATF_FIG11, times["TA2"], yerr=terr["TA2"], fmt="s-",
                 color=St.C_TA2, capsize=2, label="100% TA2")
    axB.set_xlabel(r"$a_{TF}$")
    axB.set_ylabel("Successful induction time (h)")
    axB.legend(frameon=False)
    St.panel_label(axB, "B")
    fig.tight_layout()

    # panel C: time evolution of ratios from 100% TA2, a_TF = 0.2
    path_c = os.path.join(DATA, "fig11c_tf0.2.npz")
    if os.path.exists(path_c):
        z = np.load(path_c)
        ts = z["sys_ts"][:, 0, :]          # representative run
        t = np.arange(ts.shape[0]) * 0.25
        tot = np.maximum(ts[:, 0], 1)
        axC = axB.inset_axes([1.12, 0.0, 1.05, 1.0])
        axC.plot(t, 100 * ts[:, 1] / tot, color=St.C_SC, label="SC")
        axC.plot(t, 100 * ts[:, 2] / tot, color=St.C_TA1, label="TA1")
        axC.plot(t, 100 * ts[:, 3] / tot, color=St.C_TA2, label="TA2")
        axC.set_xlabel("h")
        axC.set_ylabel("Cell ratio (%)")
        axC.set_title(r"$a_{TF}=0.2$, from 100% TA2", fontsize=7)
        axC.set_ylim(0, 100)
        axC.legend(frameon=False, fontsize=6)
        St.panel_label(axC, "C")
        fig.tight_layout()
    St.save(fig, "fig11_reprogramming", FIGS)

    with open(os.path.join(TABLES, "fig11_reprogramming.csv"), "w") as f:
        f.write("init,a_TF,success_pct,mean_time_h,n_success,n_runs\n")
        for r in rows:
            f.write(",".join(str(v) if isinstance(v, str)
                             else (str(v) if isinstance(v, int) else f"{v:g}")
                             for v in r) + "\n")


if __name__ == "__main__":
    main()
