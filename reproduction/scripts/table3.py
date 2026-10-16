"""Table 3: cell-type transition probabilities from mother to daughter
cells at one cell division (statistics window t = 2000-4000 h, 30 runs,
sigma = 0.05, m2 = 0.8, alpha = (0.4, 0.4))."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from common import DATA, TABLES  # noqa: E402
from hsl2024.analysis import transition_stats, transition_table_text  # noqa: E402

PAPER = np.array([
    [73.9402, 12.9889, 13.0709],
    [0.0013, 99.9955, 0.0032],
    [0.0013, 0.0035, 99.9952],
])


def main():
    path = os.path.join(DATA, "base.npz")
    if not os.path.exists(path):
        print("base.npz missing")
        return
    z = np.load(path)
    t = transition_stats(dict(tran=z["tran"]))
    pm = t["prob_mean"]
    txt = ("Table 3: cell type transition probabilities (%) "
           "(mother rows -> daughter columns)\n\n"
           "Reproduction (30 runs, 2000-4000 h):\n"
           + transition_table_text(pm)
           + "\n\nPaper (Table 3):\n" + transition_table_text(PAPER)
           + "\n\nAbsolute difference:\n"
           + transition_table_text(np.abs(pm - PAPER))
           + f"\n\nTotal recorded transitions: {int(t['counts'].sum())}\n")
    with open(os.path.join(TABLES, "table3_transition_probabilities.txt"),
              "w") as f:
        f.write(txt)
    np.savetxt(os.path.join(TABLES, "table3_reproduction.csv"), pm,
               delimiter=",")
    print(txt)


if __name__ == "__main__":
    main()
