"""Compare the Python reproduction with the output of the original C++ code.

Reads the .sys/.tran files produced by the authors' C++ program
(main.cpp Test(): 30 runs x 5 m2 values x 4000 h) and compares
homeostasis statistics and transition probabilities with the Python
batches (base + fig9_*).
"""
import os
import sys
import glob
import json

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from common import DATA, TABLES, M2S_FIG9  # noqa: E402
from hsl2024.analysis import homeostasis_stats, transition_stats  # noqa: E402

# Directory holding the .sys/.tran output of the authors' C++ program.
# Set $CELLTYPE_CPP_DIR to override; the default assumes the upstream
# repository (github.com/jinzhilei/Cell-type-transition) is checked out as a
# SIBLING of this repository, i.e. ../Cell-type-transition/output.  No
# absolute machine-specific path is hard-coded.
_REPO_ROOT = os.path.dirname(os.path.dirname(HERE))
CPP_DIR = os.environ.get(
    "CELLTYPE_CPP_DIR",
    os.path.join(os.path.dirname(_REPO_ROOT),
                 "Cell-type-transition", "output"))
A = [0, 1, 3]          # type codes in the C++ .tran files


def read_cpp_sys(path):
    d = np.loadtxt(path)
    return d            # columns: t, total, SC, TA1, TA2


def cpp_homeostasis(path):
    d = read_cpp_sys(path)
    seg = d[8000:, 1:]                    # rows 8000:end (1-based)
    if seg.shape[0] == 0:
        return None
    return seg.mean(axis=0)               # (total, SC, TA1, TA2)


def read_cpp_tran(path):
    d = np.loadtxt(path).astype(int)      # columns: type, parent, step
    B = np.zeros((3, 3))
    for i, pi in enumerate(A):
        m = d[:, 1] == pi
        for j, cj in enumerate(A):
            B[i, j] = ((d[:, 0] == cj) & m).sum()
    return B


def main():
    report = []
    py_names = {0.8: "base", 0.825: "fig9_m0.825", 0.85: "fig9_m0.85",
                0.875: "fig9_m0.875", 0.9: "fig9_m0.9"}
    for j, m2 in enumerate(M2S_FIG9, start=1):
        sys_files = sorted(glob.glob(os.path.join(
            CPP_DIR, "sys", f"*_{j}.sys")))
        tran_files = sorted(glob.glob(os.path.join(
            CPP_DIR, "tran", f"*_{j}.tran")))
        entry = dict(m2=m2, n_cpp_runs=len(sys_files))
        if sys_files:
            H = np.array([cpp_homeostasis(f) for f in sys_files
                          if cpp_homeostasis(f) is not None])
            entry["cpp_counts_mean"] = H.mean(axis=0).tolist()
            entry["cpp_counts_std"] = H.std(axis=0, ddof=1).tolist() \
                if len(H) > 1 else None
        if tran_files:
            Bs = np.array([read_cpp_tran(f) for f in tran_files])
            P = np.zeros_like(Bs, dtype=float)
            for k in range(Bs.shape[0]):
                rs = Bs[k].sum(axis=1, keepdims=True)
                P[k] = np.divide(Bs[k], rs,
                                 out=np.zeros_like(Bs[k], dtype=float),
                                 where=rs > 0) * 100
            entry["cpp_tran_mean"] = P.mean(axis=0).tolist()
        # python side
        py_path = os.path.join(DATA, py_names[m2] + ".npz")
        if os.path.exists(py_path):
            z = np.load(py_path)
            hs = homeostasis_stats(dict(sys_ts=z["sys_ts"],
                                        stop_step=z["stop_step"]))
            entry["py_counts_mean"] = hs["count_mean"].tolist()
            if "cpp_counts_mean" in entry:
                diff = np.array(entry["cpp_counts_mean"]) - hs["count_mean"]
                entry["abs_diff_counts"] = np.abs(diff).tolist()
            ts = transition_stats(dict(tran=z["tran"]))
            entry["py_tran_mean"] = ts["prob_mean"].tolist()
            if "cpp_tran_mean" in entry:
                entry["abs_diff_tran"] = np.abs(
                    np.array(entry["cpp_tran_mean"])
                    - ts["prob_mean"]).tolist()
        report.append(entry)
    out = os.path.join(TABLES, "cpp_vs_python.json")
    json.dump(report, open(out, "w"), indent=2)
    print(json.dumps(report, indent=2)[:3000])
    print(f"\nsaved {out}")


if __name__ == "__main__":
    main()
