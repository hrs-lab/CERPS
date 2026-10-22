"""C1 delay sweep — P8 (feedback time delay) in the class-based CERPS framework.

Rationale
---------
C1 = niche feedback (P1) x feedback delay (P8).  In the `cerps/` framework the
delay primitive is `RatesConfig.delay_steps`: it makes the niche signal Q (the SC
count in legacy mode) that feeds beta_SC be sampled `delay_steps` steps in the
past instead of the current step.  `delay_steps = 0` reproduces the
instantaneous (paper) behaviour, so it is the exact control.

Delayed negative feedback is the textbook generator of oscillations, so the
falsifiable observables are:
  1. amplitude of the stem-cell fraction fluctuation (excess over d = 0),
  2. period of that fluctuation (should grow with the delay),
  3. mean homeostatic level (may shift).

Design
------
  base   : legacy_config()  (2-gene, Table-1 kinetics, no space)
  protocol: rec_tran_from = 2000 h of a 4000 h run  (paper Table-3 protocol)
  grid   : delay_steps in {0, 50, 200, 500, 900, 1200}
            (dt = 0.25 h  ->  0, 12.5, 50, 125, 225, 300 h)
  seeds  : 0, 1, 2  (three replicates => cross-seed sd, as required after the
            [verified-prod] scale lesson: never infer replicate behaviour from
            two seeds)
  output : results/c1_delay_sweep.json

Everything below is measured with `run_cerps`-independent, plain numpy
statistics on the recorded time series.  Nothing is hand-written.
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cerps.config import replace, legacy_config  # noqa: E402
from cerps.batch import run_one  # noqa: E402

# 0 -> 1200 steps.  With dt = 0.25 h this spans 0 -> 300 h (12.5 days), which
# brackets the 9-12 day delayed dispersion effect reported by Dray et al. 2021
# (Cell Stem Cell) - the biological motivation for the P8 axis.  A short-only
# scan (<= 50 h) showed no detectable effect, so the long delays are essential.
DELAYS = (0, 50, 200, 500, 900, 1200)
SEEDS = (0, 1, 2)
T_FROM = 2000.0          # homeostatic window start (paper protocol)
REC_TRAN_FROM = 2000.0
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results", "c1_delay_sweep.json")


def fluctuation(t, sc, tot, t_from):
    """Amplitude / period / regularity of the SC-fraction fluctuation.

    Returns std, dominant period (h) and the autocorrelation value at that lag,
    computed on the de-meaned SC fraction over [t_from, end].
    """
    t = np.asarray(t, float)
    sc = np.asarray(sc, float)
    tot = np.asarray(tot, float)
    m = t >= t_from
    if m.sum() < 20:
        return dict(std=None, period=None, peak_ac=None)
    tt = t[m]
    f = np.where(tot[m] > 0, sc[m] / np.maximum(tot[m], 1e-9), 0.0)
    f = f - f.mean()
    sd = float(f.std())
    if sd < 1e-12:
        return dict(std=0.0, period=None, peak_ac=None)
    # normalised autocorrelation via FFT
    n = len(f)
    nf = int(2 ** np.ceil(np.log2(2 * n)))
    F = np.fft.rfft(f - f.mean(), n=nf)
    ac = np.fft.irfft(F * np.conjugate(F), n=nf)[:n]
    if ac[0] <= 0:
        return dict(std=sd, period=None, peak_ac=None)
    ac /= ac[0]
    # first zero crossing -> first local maximum after it = dominant period
    nz = np.argmax(ac <= 0) if np.any(ac <= 0) else 1
    if nz >= n - 2:
        return dict(std=sd, period=None, peak_ac=None)
    k = nz + int(np.argmax(ac[nz:]))
    return dict(std=sd,
                period=float(tt[k] - tt[0]) if k < len(tt) else None,
                peak_ac=float(ac[k]) if k < len(ac) else None)


def main():
    base = legacy_config()
    # paper protocol: discard transients before 2000 h when filling transitions;
    # record_states=False keeps only the count time series (memory/time).
    base = replace(base, rec_tran_from=REC_TRAN_FROM, record_states=False)

    rows = []
    t0 = time.time()
    for d in DELAYS:
        for sd in SEEDS:
            cfg = replace(base, rates=replace(base.rates, delay_steps=d), seed=sd)
            r = run_one(cfg, seed=sd, tag=f"delay_{d}", t_from=T_FROM,
                        keep_sim=True, warn=False)
            sim = r["sim"]
            fl = fluctuation(sim.records["t"], sim.records["sc"],
                             sim.records["total"], T_FROM)
            row = dict(delay_steps=int(d), delay_h=float(d * base.dt), seed=int(sd),
                       ho_total=r["ho_total"], ho_sc=r["ho_sc"],
                       ho_sc_frac=r["ho_sc_frac"], ho_ta1=r["ho_ta1"],
                       ho_ta2=r["ho_ta2"],
                       n_final=r["n_final"], max_cells=r["max_cells"],
                       sc_frac_std=fl["std"], sc_period_h=fl["period"],
                       sc_peak_ac=fl["peak_ac"],
                       tran_SC_SC=r.get("tran_0_0"), tran_SC_TA1=r.get("tran_0_1"),
                       tran_SC_TA2=r.get("tran_0_3"),
                       runtime_s=round(r["runtime_s"], 1))
            rows.append(row)
            print(f"  d={d:4d} ({d*base.dt:5.2f}h) seed={sd} "
                  f"sc_frac={row['ho_sc_frac']:.3f} std={row['sc_frac_std']:.4f} "
                  f"period={row['sc_period_h']} ac={row['sc_peak_ac']} "
                  f"({row['runtime_s']}s)", flush=True)

    # ---- aggregate over seeds ------------------------------------------------
    agg = []
    for d in DELAYS:
        rs = [r for r in rows if r["delay_steps"] == d]
        def ms(k):
            v = np.array([r[k] for r in rs if r[k] is not None], float)
            return (float(v.mean()), float(v.std(ddof=1))) if len(v) > 1 else \
                   (float(v[0]) if len(v) else None, None)
        sf, sfs = ms("ho_sc_frac")
        st, sts = ms("sc_frac_std")
        pe, pes = ms("sc_period_h")
        to, tos = ms("ho_total")
        agg.append(dict(delay_steps=int(d), delay_h=float(d * base.dt),
                        sc_frac=sf, sc_frac_sd=sfs,
                        fluct_std=st, fluct_std_sd=sts,
                        period_h=pe, period_h_sd=pes,
                        total=to, total_sd=tos))

    out = dict(config=dict(dt=base.dt, max_time=base.max_time, n_init=base.n_init,
                           rec_tran_from=REC_TRAN_FROM, t_from=T_FROM,
                           delays=list(DELAYS), seeds=list(SEEDS),
                           base="legacy_config()"),
               runs=rows, aggregate=agg,
               elapsed_s=round(time.time() - t0, 1))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nsaved {OUT} ({out['elapsed_s']}s)")
    print("\n===== AGGREGATE (mean +- sd over 3 seeds) =====")
    print(f"{'delay_steps':>11} {'delay_h':>8} {'sc_frac':>16} {'fluct_std':>18} {'period_h':>16}")
    for a in agg:
        print(f"{a['delay_steps']:>11} {a['delay_h']:>8.2f} "
              f"{a['sc_frac']:>10.4f}±{a['sc_frac_sd']:.4f} "
              f"{a['fluct_std']:>10.4f}±{a['fluct_std_sd']:.4f} "
              f"{str(round(a['period_h'],1)) if a['period_h'] else 'None':>16}")


if __name__ == "__main__":
    main()
