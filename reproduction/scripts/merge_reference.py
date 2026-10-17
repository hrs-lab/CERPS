"""Merge the 30 reference-simulator replicates (ref_runXX.npz) and
compare with the numba engine (base.npz): counts + transition matrices."""
import os
import sys
import json

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from common import DATA, TABLES  # noqa: E402
from hsl2024.analysis import homeostasis_stats, transition_stats  # noqa: E402


def main():
    refs = sorted(f for f in os.listdir(DATA) if f.startswith("ref_run"))
    n = len(refs)
    print(f"{n} reference replicates found")
    if n < 30:
        print("(waiting for all 30)")

    sys_all, tran_all = [], np.zeros((1, 3, 3), np.int64)
    for f in refs:
        z = np.load(os.path.join(DATA, f))
        sys_all.append(z["sys_ts"][:, 0, :])       # (steps, 4)
        tran_all += z["tran"]
    L = min(s.shape[0] for s in sys_all)
    sys_ref = np.stack([s[:L] for s in sys_all], axis=1)  # (steps, n, 4)

    hs_ref = homeostasis_stats(dict(sys_ts=sys_ref, stop_step=np.full(n, -1)))
    tr_ref = transition_stats(dict(tran=tran_all))

    z = np.load(os.path.join(DATA, "base.npz"))
    hs_num = homeostasis_stats(dict(sys_ts=z["sys_ts"][:, :n] if n < 30
                                    else z["sys_ts"],
                                    stop_step=np.full(max(n, 1), -1)))
    tr_num = transition_stats(dict(tran=z["tran"]))

    names = ["total", "SC", "TA1", "TA2"]
    print(f"\n{'':8s}{'reference':>12s}{'numba':>12s}{'diff%':>8s}")
    for i, nm in enumerate(names):
        r, m = hs_ref["count_mean"][i], hs_num["count_mean"][i]
        print(f"{nm:8s}{r:12.1f}{m:12.1f}{100*abs(r-m)/max(m,1):8.2f}")
    print("\ntransition matrix (rows mother, cols daughter), %:")
    print("reference:\n", np.round(tr_ref["prob_mean"], 4))
    print("numba:\n", np.round(tr_num["prob_mean"], 4))
    d = np.abs(tr_ref["prob_mean"] - tr_num["prob_mean"])
    print("max abs diff (pp):", d.max(), " (SC row):", d[0].max())

    out = dict(
        n_ref_runs=n,
        ref_counts=hs_ref["count_mean"].tolist(),
        numba_counts=hs_num["count_mean"].tolist(),
        ref_tran=tr_ref["prob_mean"].tolist(),
        numba_tran=tr_num["prob_mean"].tolist(),
        max_tran_diff_pp=float(d.max()),
    )
    json.dump(out, open(os.path.join(TABLES, "reference_vs_numba.json"),
                        "w"), indent=2)
    print(f"\nsaved {TABLES}/reference_vs_numba.json")


if __name__ == "__main__":
    main()
