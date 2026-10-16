"""Statistics helpers mirroring DataAnalysis.m of the original repository."""
import numpy as np


# Matlab uses rows 8000:end of the .sys file (1-based).  Row k corresponds to
# 0-based step k-1, so the statistics window is step >= 7999 (0-based).
HOMEOSTASIS_FROM_STEP = 8000  # 0-based inclusive lower bound used here


def homeostasis_stats(result, from_step=None):
    """Per-run time-averaged counts/ratios over the statistics window.

    Returns dict with mean/std across runs for (total, SC, TA1, TA2) numbers
    and ratios, exactly as in DataAnalysis.m (per-run mean of each column,
    then mean and std across runs).
    """
    sys_ts = result["sys_ts"]          # (steps, runs, 4)
    n_runs = sys_ts.shape[1]
    if from_step is None:
        from_step = HOMEOSTASIS_FROM_STEP
    stop = result.get("stop_step", None)
    means = np.zeros((n_runs, 4))
    for r in range(n_runs):
        end = n_steps_r = sys_ts.shape[0]
        if stop is not None and stop[r] >= 0:
            end = min(end, stop[r] + 1)
        seg = sys_ts[from_step:end, r, :]
        if seg.shape[0] == 0:          # run ended before the window
            seg = sys_ts[max(end - 1, 0):end, r, :]
        means[r] = seg.mean(axis=0)
    ratios = 100.0 * means[:, 1:] / means[:, [0]]      # SC, TA1, TA2 (%)
    finite = np.isfinite(ratios)
    ratios = np.where(finite, ratios, 0.0)
    return dict(
        count_mean=means.mean(axis=0), count_std=means.std(axis=0, ddof=1),
        ratio_mean=ratios.mean(axis=0), ratio_std=ratios.std(axis=0, ddof=1),
        per_run_counts=means, per_run_ratios=ratios,
    )


def transition_stats(result):
    """Row-normalised mother->daughter transition probabilities (%).

    Mirrors DataAnalysis.m: per run, C = 100*B/rowsum(B); report mean and
    std across runs.  Rows/cols are ordered (SC, TA1, TA2).
    """
    tran = result["tran"]              # (runs, 3, 3)
    n_runs = tran.shape[0]
    perc = np.zeros((n_runs, 3, 3))
    for r in range(n_runs):
        B = tran[r].astype(float)
        rs = B.sum(axis=1)
        for i in range(3):
            if rs[i] > 0:
                perc[r, i] = 100.0 * B[i] / rs[i]
    return dict(prob_mean=perc.mean(axis=0), prob_std=perc.std(axis=0, ddof=1),
                counts=tran.sum(axis=0), per_run=perc)


TYPE_NAMES = ("SC", "TA1", "TA2")


def transition_table_text(prob_mean):
    lines = ["mother\\daughter   SC        TA1       TA2"]
    for i, name in enumerate(TYPE_NAMES):
        row = "  ".join(f"{v:8.4f}" for v in prob_mean[i])
        lines.append(f"{name:<14}{row}")
    return "\n".join(lines)


def landscape_from_hist(result, vmin=1.5, vmax=4.0):
    """Pooled U(t, z) = -log f(t, z) landscape (Fig. 13)."""
    hist = result["hist"]              # (runs, n_tb, n_zb)
    h = hist.sum(axis=0).astype(float)  # (n_tb, n_zb)
    tot = h.sum(axis=1, keepdims=True)
    f = np.divide(h, tot, out=np.zeros_like(h), where=tot > 0)
    with np.errstate(divide="ignore"):
        U = -np.log(f)
    U[~np.isfinite(U)] = vmax
    U = np.clip(U, vmin, vmax)
    return U                            # (n_tb, n_zb)
