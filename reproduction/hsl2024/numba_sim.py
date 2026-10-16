"""Numba-JIT individual-based simulator of the hybrid model.

This module re-implements, in Python + numba, the C++ simulation released by
the authors (github.com/jinzhilei/Cell-type-transition, GPL-3.0).  The logic
mirrors ``Header.h`` exactly, including details where the code differs from
the text of the paper (see docs/protocol_notes.md):

* phenotype classification: SC  x1>0.5 & x2>0.5; TA1 x1>0.5 & x2<0.5;
  TA2 x1<0.5 & x2>0.5; both low  ->  previous type is kept (no TC type);
* phenotype is re-evaluated only during the resting phase (stages 0/3);
* the division counter ``div`` increments at mitosis; a TA cell whose
  ``div`` reaches 15 dies at mitosis instead of producing daughters;
* a resting cell classified as SC has ``div`` reset to 0;
* transition events are recorded when a cell enters the proliferating
  phase (stage==1): (current type at entry, parent type at last mitosis);
* beta_SC = beta0 / (1 + Q/theta0) with Q = SC count of the current step
  (after phenotype update, before fate decisions);
* daughters inherit eta (extrinsic noise state) and the gene expression
  of the mother; their epigenetic state is redrawn from
  Beta(m0*phi, m0*(1-phi)) independently per gene (unless ``fixed_u``);
* OU noise is integrated with the Euler-Maruyama scheme of the C++ code,
  and the gene ODE with one forward-Euler step per dt.

The simulator is *batched*: ``n_runs`` replicate runs are advanced
simultaneously (each cell carries its run id).  Statistical behaviour is
identical to running the runs separately.
"""
import numpy as np
from numba import njit


# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------
@njit(cache=True, nogil=True)
def _beta_draw(a, b):
    """Beta(a, b) from two Gamma draws (same construction as Random.h)."""
    x = np.random.gamma(a, 1.0)
    y = np.random.gamma(b, 1.0)
    s = x + y
    if s <= 0.0:
        return 0.5
    return x / s


# ----------------------------------------------------------------------
# main simulator
# ----------------------------------------------------------------------
@njit(cache=True, nogil=True)
def run_simulation(
    # GRN parameters
    b1, b2, k1, k2, s1, s2, rho1, rho2, nhill,
    alpha1, alpha2, lam1, lam2, zeta1, zeta2, sigma, a_tf,
    # cell kinetics
    mu0, kappa0, beta0, tau0, theta0, max_div_ta,
    # epigenetic inheritance
    m0, m1, m2, m3, m4,
    # protocol
    dt, n_steps, n_runs, n_init, init_mode, fixed_u,
    rec_tran_from, cap_total, stop_sc_target, seed,
    # outputs -----------------------------------------------------------
    sys_ts,        # (n_steps, n_runs, 4) int32  rows: total,SC,TA1,TA2
    tran,          # (n_runs, 3, 3) int64  [parent, current] in {SC,TA1,TA2}
    stop_step,     # (n_runs,) int32  step index when run stopped (-1 = n/a)
    stop_kind,     # (n_runs,) int8   0 ran to end, 1 extinct, 2 success
    max_cells_run, # (n_runs,) int32
    # optional landscape histogram of z = x1 - x2 (Fig. 13)
    do_landscape, hist, tb_stride,
    # optional state snapshots of run 0 (Fig. 5)
    snap_steps, snap_data, snap_n,
    # per-run initial (u1,u2) and (x1,x2) for init_mode != 0
    u_init, x_init,
    out_status,    # (1,) int32 : 0 ok, 1 capacity exceeded
):
    np.random.seed(seed)
    out_status[0] = 0
    n_tb = hist.shape[1]
    n_zb = hist.shape[2]

    # precomputed constants ------------------------------------------------
    s1n = s1 ** nhill
    s2n = s2 ** nhill
    hill_pow = (nhill != 4.0) and (nhill != 2.0)  # fast paths below
    ou_a1 = dt / zeta1
    ou_b1 = np.sqrt(2.0 * dt / zeta1)
    ou_a2 = dt / zeta2
    ou_b2 = np.sqrt(2.0 * dt / zeta2)
    half_sig1 = 0.5 * sigma * sigma
    half_sig2 = 0.5 * sigma * sigma
    m3m4_ready = 0.0  # (placeholder; phi computed inline)

    # state arrays ---------------------------------------------------------
    cap = cap_total
    X1 = np.empty(cap, np.float64)
    X2 = np.empty(cap, np.float64)
    U1 = np.empty(cap, np.float64)
    U2 = np.empty(cap, np.float64)
    E1 = np.zeros(cap, np.float64)
    E2 = np.zeros(cap, np.float64)
    STAGE = np.zeros(cap, np.int8)     # 0 rest,1 just entered prolif,2 prolif,3 newborn
    AGE = np.zeros(cap, np.float64)
    DIV = np.zeros(cap, np.int64)
    CTYPE = np.zeros(cap, np.int8)
    PTYPE = np.zeros(cap, np.int8)
    RUN = np.zeros(cap, np.int32)

    div_idx = np.empty(cap, np.int32)
    div_phi1 = np.empty(cap, np.float64)
    div_phi2 = np.empty(cap, np.float64)

    run_active = np.ones(n_runs, np.bool_)
    for r in range(n_runs):
        stop_step[r] = -1
        stop_kind[r] = 0
        max_cells_run[r] = 0

    # initialisation --------------------------------------------------------
    n = 0
    for r in range(n_runs):
        for c in range(n_init):
            if init_mode == 0:
                X1[n] = 3.0 * np.random.random()
                X2[n] = 3.0 * np.random.random()
                U1[n] = np.random.random()
                U2[n] = np.random.random()
            else:
                X1[n] = x_init[r, 0]
                X2[n] = x_init[r, 1]
                U1[n] = u_init[r, 0]
                U2[n] = u_init[r, 1]
            E1[n] = 0.0
            E2[n] = 0.0
            STAGE[n] = 0
            AGE[n] = 0.0
            DIV[n] = 0
            PTYPE[n] = 0
            RUN[n] = r
            t = 0
            if X1[n] > 0.5 and X2[n] > 0.5:
                t = 0
            elif X1[n] > 0.5 and X2[n] < 0.5:
                t = 1
            elif X1[n] < 0.5 and X2[n] > 0.5:
                t = 3
            CTYPE[n] = t
            n += 1

    counts = np.zeros((n_runs, 4), np.int64)

    # main loop -------------------------------------------------------------
    for step in range(n_steps):
        # pre-generated random buffers for this step (one C-level fill each)
        nbuf = np.random.standard_normal(2 * n)
        ubuf = np.random.random(2 * n)
        ucur = 0
        # ---------------- pass A: dynamics + counting ----------------------
        for r in range(n_runs):
            counts[r, 0] = 0
            counts[r, 1] = 0
            counts[r, 2] = 0
            counts[r, 3] = 0
        for i in range(n):
            r = RUN[i]
            if not run_active[r]:
                continue
            # transition record (cells that just entered proliferation)
            if STAGE[i] == 1 and step + 1 > rec_tran_from:
                if PTYPE[i] == 0:
                    pi = 0
                elif PTYPE[i] == 1:
                    pi = 1
                else:
                    pi = 2
                if CTYPE[i] == 0:
                    ci = 0
                elif CTYPE[i] == 1:
                    ci = 1
                else:
                    ci = 2
                tran[r, pi, ci] += 1
            # OU noise (Euler-Maruyama, as in Random-based C++ update)
            e1 = E1[i] - E1[i] * ou_a1 + ou_b1 * nbuf[ucur]
            e2 = E2[i] - E2[i] * ou_a2 + ou_b2 * nbuf[ucur + 1]
            ucur += 2
            E1[i] = e1
            E2[i] = e2
            a1 = alpha1 * np.exp(lam1 * U1[i] + sigma * e1 - half_sig1) + a_tf
            a2 = alpha2 * np.exp(lam2 * U2[i] + sigma * e2 - half_sig2)
            # gene ODE, one forward-Euler step
            x1 = X1[i]
            x2 = X2[i]
            if hill_pow:
                x1p = x1 ** nhill
                x2p = x2 ** nhill
            elif nhill == 4.0:
                x1n = x1 * x1
                x2n = x2 * x2
                x1p = x1n * x1n
                x2p = x2n * x2n
            else:
                x1p = x1 * x1
                x2p = x2 * x2
            f1 = a1 * (rho1 + (1.0 - rho1) * x1p / (s1n + x1p)) \
                + b1 * s2n / (s2n + x2p) - k1 * x1
            f2 = a2 * (rho2 + (1.0 - rho2) * x2p / (s2n + x2p)) \
                + b2 * s1n / (s1n + x1p) - k2 * x2
            x1 = x1 + f1 * dt
            x2 = x2 + f2 * dt
            X1[i] = x1
            X2[i] = x2
            # phenotype update (resting only)
            if STAGE[i] == 0 or STAGE[i] == 3:
                t = CTYPE[i]
                if x1 > 0.5 and x2 > 0.5:
                    t = 0
                elif x1 > 0.5 and x2 < 0.5:
                    t = 1
                elif x1 < 0.5 and x2 > 0.5:
                    t = 3
                CTYPE[i] = t
                if t == 0:
                    DIV[i] = 0
            # counting
            counts[r, 0] += 1
            ct = CTYPE[i]
            if ct == 0:
                counts[r, 1] += 1
            elif ct == 1:
                counts[r, 2] += 1
            elif ct == 3:
                counts[r, 3] += 1
            # landscape histogram (z = x1 - x2)
            if do_landscape and step < tb_stride * n_tb:
                tb = step // tb_stride
                z = x1 - x2
                zb = int((z + 4.0) * (n_zb / 8.0))
                if zb < 0:
                    zb = 0
                elif zb >= n_zb:
                    zb = n_zb - 1
                hist[r, tb, zb] += 1
        # snapshots (run 0)
        for s_i in range(snap_steps.shape[0]):
            if snap_steps[s_i] == step:
                cnt = 0
                for i in range(n):
                    if RUN[i] == 0:
                        snap_data[s_i, cnt, 0] = X1[i]
                        snap_data[s_i, cnt, 1] = X2[i]
                        snap_data[s_i, cnt, 2] = U1[i]
                        snap_data[s_i, cnt, 3] = U2[i]
                        cnt += 1
                snap_n[s_i] = cnt
        # store counts
        for r in range(n_runs):
            sys_ts[step, r, 0] = counts[r, 0]
            sys_ts[step, r, 1] = counts[r, 1]
            sys_ts[step, r, 2] = counts[r, 2]
            sys_ts[step, r, 3] = counts[r, 3]
            if counts[r, 0] > max_cells_run[r]:
                max_cells_run[r] = counts[r, 0]
        # success stopping (Fig. 11)
        if stop_sc_target > 0:
            for r in range(n_runs):
                if run_active[r] and counts[r, 1] >= stop_sc_target:
                    run_active[r] = False
                    stop_step[r] = step
                    stop_kind[r] = 2
        # ---------------- pass B: fate decisions ---------------------------
        w = 0
        ndiv = 0
        alive = np.zeros(n_runs, np.int64)
        for i in range(n):
            r = RUN[i]
            if not run_active[r]:
                continue
            ct = CTYPE[i]
            if ct == 0:
                beta = beta0 / (1.0 + counts[r, 1] / theta0)
                kap = kappa0
                tau = tau0
            else:
                beta = beta0
                kap = 2.0 * kappa0
                tau = 0.5 * tau0
            if STAGE[i] == 0 or STAGE[i] == 3:
                r1 = ubuf[2 * i]
                if r1 < kap * dt:
                    continue                     # removal from resting pool
                if r1 < (kap + beta) * dt:
                    STAGE[i] = 1                 # enter proliferation
                    AGE[i] = 0.0
                else:
                    STAGE[i] = 0                 # remain resting
            else:
                r2 = ubuf[2 * i + 1]
                if r2 < mu0 * dt:
                    continue                     # apoptosis
                if AGE[i] < tau:
                    STAGE[i] = 2
                    AGE[i] = AGE[i] + dt
                else:
                    # mitosis
                    DIV[i] = DIV[i] + 1
                    if DIV[i] == max_div_ta and ct != 0:
                        continue                 # TA senescence death
                    if not fixed_u:
                        # phi from the mother's pre-division epigenetic state
                        v1 = m3 * U1[i]
                        v2 = m3 * U2[i]
                        v1p = v1 ** m4
                        v2p = v2 ** m4
                        phi1 = m1 + m2 * v1p / (1.0 + v1p)
                        phi2 = m1 + m2 * v2p / (1.0 + v2p)
                        U1[i] = _beta_draw(m0 * phi1, m0 * (1.0 - phi1))
                        U2[i] = _beta_draw(m0 * phi2, m0 * (1.0 - phi2))
                        div_idx[ndiv] = w
                        div_phi1[ndiv] = phi1
                        div_phi2[ndiv] = phi2
                        ndiv += 1
                    STAGE[i] = 3
                    PTYPE[i] = ct
            # compact cell i to slot w
            if w != i:
                X1[w] = X1[i]
                X2[w] = X2[i]
                U1[w] = U1[i]
                U2[w] = U2[i]
                E1[w] = E1[i]
                E2[w] = E2[i]
                STAGE[w] = STAGE[i]
                AGE[w] = AGE[i]
                DIV[w] = DIV[i]
                CTYPE[w] = CTYPE[i]
                PTYPE[w] = PTYPE[i]
                RUN[w] = RUN[i]
            w += 1
            alive[r] += 1
        # append daughter cells
        if w + ndiv > cap:
            out_status[0] = 1
            return n
        for d in range(ndiv):
            j = w + d
            mi = div_idx[d]
            X1[j] = X1[mi]
            X2[j] = X2[mi]
            E1[j] = E1[mi]
            E2[j] = E2[mi]
            U1[j] = U1[mi]
            U2[j] = U2[mi]
            if not fixed_u:
                U1[j] = _beta_draw(m0 * div_phi1[d], m0 * (1.0 - div_phi1[d]))
                U2[j] = _beta_draw(m0 * div_phi2[d], m0 * (1.0 - div_phi2[d]))
            STAGE[j] = 3
            AGE[j] = AGE[mi]
            DIV[j] = DIV[mi]
            CTYPE[j] = CTYPE[mi]
            PTYPE[j] = PTYPE[mi]
            RUN[j] = RUN[mi]
            alive[RUN[mi]] += 1
        n = w + ndiv
        # extinction stopping
        if stop_sc_target > 0:
            for r in range(n_runs):
                if run_active[r] and alive[r] == 0:
                    run_active[r] = False
                    stop_step[r] = step
                    stop_kind[r] = 1
        # early exit when every run has stopped
        still = 0
        for r in range(n_runs):
            if run_active[r]:
                still += 1
        if still == 0:
            for s in range(step + 1, n_steps):
                for r in range(n_runs):
                    sys_ts[s, r, 0] = 0
                    sys_ts[s, r, 1] = 0
                    sys_ts[s, r, 2] = 0
                    sys_ts[s, r, 3] = 0
            break
    return n
