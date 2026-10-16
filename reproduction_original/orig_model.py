"""Faithful Python/numba reproduction of:
    Huang, Situ, Lei (2024). Dynamics of cell-type transition mediated by
    epigenetic modifications. J. Theor. Biol. 577:111664.
    Original C++ source: github.com/jinzhilei/Cell-type-transition

This port follows the ORIGINAL C++ SEMANTICS exactly (verified by running the
compiled C++ code side-by-side; see docs/SPEC_MODELS.md section 1):

  * Hill terms: H(x)=x^n/(s^n+x^n), repression R(x)=s^n/(s^n+x^n)  [s^n!]
  * a_i = alpha_i*exp(lam_i*u_i)*exp(sig_i*eta_i - sig_i^2/2);
    TF induction (Eq 9): a_1 -> a_1 + aTF inside the self-activation term.
  * OU noise updated by EXPLICIT EULER: eta += (-eta/zeta)*dt + sqrt(2*dt/zeta)*N(0,1);
    eta initialized to 0 (C++ zeroes all CellInfoDouble at birth).
  * Phenotype codes SC=0, TA1=1, TA2=3 (NO TC code): UpdateCellType only changes
    the label when (x1>0.5 & x2>0.5)->SC, (x1>0.5 & x2<0.5)->TA1,
    (x1<0.5 & x2>0.5)->TA2; in the low-low quadrant the OLD label is retained
    (stale label, keeps its kinetics). Label update only while resting
    (stage 0 = resting, 3 = just-divided daughter). SC label resets division count.
  * Cell stages: -1 dead; 0 resting->resting; 1 resting->proliferating (entry);
    2 proliferating; 3 proliferating->resting (daughter).
  * Rates by CURRENT label: SC: beta=beta0/(1+NumSC/theta), kappa=kappa0,
    mu=mu0, tau=tau0; TA1/TA2: beta=beta0, kappa=2*kappa0, mu=mu0, tau=tau0/2.
    (NumSC = current SC count, computed after phase-1 label updates.)
  * Fate per step (two uniforms r1,r2 drawn per cell): resting: r1<kappa*dt ->
    die; elif r1<(kappa+beta)*dt -> enter proliferating (stage 1, age 0);
    proliferating: r2<mu*dt -> die; elif age<tau -> stay (stage 2, age+=dt);
    else -> DIVISION: daughters stage 3, div = mother_div+1; if daughter div
    >= 15 and label is TA -> both daughters die. Daughters inherit x, eta,
    label, parent label = mother label; u_i resampled per daughter per gene:
    u ~ Beta(m0*phi, m0*(1-phi)), phi = m1+m2*(m3*v)^m4/(1+(m3*v)^m4).
  * Transition logging (Table 3): at PROLIFERATION ENTRY, count
    trans[parent_label, current_label] for events with t > log_after (2000 h),
    i.e. statistics over [2000,4000] h of a 4000 h run, per-run row-normalized
    then averaged over runs (MATLAB DataAnalysis.m protocol).
  * Population hard-capped at MAXCELL=12000 (excess newborns dropped, as in C++).
  * C++ ground truth (1 run, m2=0.8, T=4000h): SC row [74.45, 12.47, 13.08],
    TA rows ~[0.001, 99.995, 0.004]; means over t in [2000,4000]:
    total=6984, SC=177, TA1=3269, TA2=3537.

Usage:
    python orig_model.py [seeds...]      # Table-3 style validation
"""
import numpy as np
from numba import njit

SC, TA1, TA2 = 0, 1, 3
MAXCELL = 12000

PARAM_KEYS = ('b1', 'b2', 'k1', 'k2', 's1', 's2', 'rho1', 'rho2', 'sig1', 'sig2',
              'zeta1', 'zeta2', 'alpha1', 'alpha2', 'lam1', 'lam2', 'n',
              'm0', 'm1', 'm2', 'm3', 'm4',
              'theta', 'beta0', 'kappa0', 'mu0', 'tau0',
              'dt', 'aTF', 'maxdiv', 'thx')


def paper_params(**over):
    P = dict(b1=1., b2=1., k1=1., k2=1., s1=.5, s2=.5, rho1=.1, rho2=.1,
             sig1=.05, sig2=.05, zeta1=1., zeta2=1.,
             alpha1=.4, alpha2=.4, lam1=1.9, lam2=1.9, n=4.,
             m0=60., m1=.08, m2=.8, m3=2., m4=2.,
             theta=200., beta0=.04, kappa0=.01, mu0=.0002, tau0=4.,
             dt=.25, aTF=0., maxdiv=15., thx=.5)
    P.update(over)
    return tuple(float(P[k]) for k in PARAM_KEYS)


@njit(cache=True)
def _new_type(x1, x2, thx, cur):
    """C++ UpdateCellType: stale-label semantics (no TC code)."""
    if x1 > thx and x2 > thx:
        return 0
    if x1 > thx and x2 < thx:
        return 1
    if x1 < thx and x2 > thx:
        return 3
    return cur


@njit(cache=True)
def _beta_u(g, h):
    g1 = np.random.gamma(g, 1.0)
    g2 = np.random.gamma(h, 1.0)
    return g1 / (g1 + g2)


@njit(cache=True)
def _trow(c):
    """Map type code (0,1,3) to matrix row/col (0,1,2)."""
    if c == 0:
        return 0
    if c == 1:
        return 1
    return 2


@njit(cache=True)
def run_sim(seed, P, T, init_mode, iv, fixed_u, stop_sc, snap_times, log_after):
    """One model run with exact C++ semantics.

    init_mode: 0 = 100 cells, x~U(0,3)^2, u~U(0,1)^2, eta=0, labels via rules
               1 = iv[4] cells at state iv[0..3]=(x1,x2,u1,u2)
               3 = SINGLE cell at iv[0..3]; population dynamics disabled;
                   record (t,x1,x2) trajectory (Fig 2C/D)
    fixed_u:   1 = daughters inherit u exactly (Fig 3)
    stop_sc:   >0 -> early stop when SC count reaches this (Fig 11)
    log_after: transition events counted only for t > log_after (default 2000)

    Returns (counts, trans, snaps, snap_n, traj, scal)
      counts : (nH+1, 7) hourly [t, N, SC, TA1, TA2, lowlow, ndiv_cum]
      trans  : (3,3) [parent, current] event counts (t>log_after)
      snaps  : (nS, cap, 4) [x1,x2,u1,u2]; snap_n per-snap counts
      traj   : (nst, 3) [t,x1,x2] for init_mode 3 (else 1 dummy row)
      scal   : [ndiv, t_success, t_ext, n_final, nsc_final, T_used, trunc_events]
    """
    (b1, b2, k1, k2, s1, s2, rho1, rho2, sig1, sig2, zeta1, zeta2,
     alpha1, alpha2, lam1, lam2, nH, m0, m1, m2_, m3, m4,
     theta, beta0, kappa0, mu0, tau0, dt, aTF, maxdiv, thx) = P
    np.random.seed(seed)

    MAX2 = 2 * MAXCELL
    nst = int(round(T / dt))
    nS = snap_times.shape[0]
    s1n = s1 ** nH
    s2n = s2 ** nH
    sq1 = np.sqrt(2.0 * dt / zeta1)
    sq2 = np.sqrt(2.0 * dt / zeta2)

    # cell state (pool)
    x1a = np.zeros(MAXCELL); x2a = np.zeros(MAXCELL)
    u1a = np.zeros(MAXCELL); u2a = np.zeros(MAXCELL)
    e1a = np.zeros(MAXCELL); e2a = np.zeros(MAXCELL)
    cty = np.zeros(MAXCELL, np.int64)   # current type
    pty = np.zeros(MAXCELL, np.int64)   # parent type
    stg = np.zeros(MAXCELL, np.int64)   # stage: 0,1,2,3
    dvn = np.zeros(MAXCELL, np.int64)   # division number
    age = np.zeros(MAXCELL)
    # new pool (temp)
    nx1 = np.zeros(MAX2); nx2 = np.zeros(MAX2)
    nu1 = np.zeros(MAX2); nu2 = np.zeros(MAX2)
    ne1 = np.zeros(MAX2); ne2 = np.zeros(MAX2)
    nct = np.zeros(MAX2, np.int64); npt = np.zeros(MAX2, np.int64)
    nst_ = np.zeros(MAX2, np.int64); ndv = np.zeros(MAX2, np.int64)
    nag = np.zeros(MAX2)

    if init_mode == 0:
        n = 100
        for i in range(n):
            x1a[i] = np.random.rand() * 3.0
            x2a[i] = np.random.rand() * 3.0
            u1a[i] = np.random.rand()
            u2a[i] = np.random.rand()
            cty[i] = _new_type(x1a[i], x2a[i], thx, 0)
            pty[i] = 0
            stg[i] = 0
    elif init_mode == 1:
        n = int(iv[4])
        for i in range(n):
            x1a[i] = iv[0]; x2a[i] = iv[1]; u1a[i] = iv[2]; u2a[i] = iv[3]
            cty[i] = _new_type(x1a[i], x2a[i], thx, 0)
            pty[i] = 0
            stg[i] = 0
    else:
        n = 1
        x1a[0] = iv[0]; x2a[0] = iv[1]; u1a[0] = iv[2]; u2a[0] = iv[3]
        cty[0] = _new_type(x1a[0], x2a[0], thx, 0)
        pty[0] = 0; stg[0] = 0

    nHrec = int(T) + 2
    counts = np.zeros((nHrec, 7))
    trans = np.zeros((3, 3), np.int64)
    SNAP_CAP = 3000
    snaps = np.zeros((nS, SNAP_CAP, 4))
    snap_n = np.zeros(nS, np.int64)
    traj = np.zeros((nst, 3))
    scal = np.zeros(7)

    t = 0.0
    ndiv = 0
    t_success = -1.0
    t_ext = -1.0
    snap_i = 0
    trunc_events = 0
    hour_i = 0

    for step in range(nst):
        # ---------------- phase 1: RunOneStep + label update + type counts
        nsc = 0; nta1 = 0; nta2 = 0; nll = 0
        for i in range(n):
            e1a[i] = e1a[i] + (-e1a[i] / zeta1) * dt + sq1 * np.random.randn()
            e2a[i] = e2a[i] + (-e2a[i] / zeta2) * dt + sq2 * np.random.randn()
            a1 = alpha1 * np.exp(lam1 * u1a[i]) * np.exp(sig1 * e1a[i] - 0.5 * sig1 * sig1)
            a2 = alpha2 * np.exp(lam2 * u2a[i]) * np.exp(sig2 * e2a[i] - 0.5 * sig2 * sig2)
            p1 = x1a[i] ** nH
            p2 = x2a[i] ** nH
            h1 = p1 / (s1n + p1)
            h2 = p2 / (s2n + p2)
            r1 = s2n / (s2n + p2)
            r2 = s1n / (s1n + p1)
            x1a[i] += ((a1 + aTF) * (rho1 + (1.0 - rho1) * h1) + b1 * r1 - k1 * x1a[i]) * dt
            x2a[i] += (a2 * (rho2 + (1.0 - rho2) * h2) + b2 * r2 - k2 * x2a[i]) * dt
            if x1a[i] < 0.0:
                x1a[i] = 0.0
            if x2a[i] < 0.0:
                x2a[i] = 0.0
            if stg[i] == 0 or stg[i] == 3:
                cty[i] = _new_type(x1a[i], x2a[i], thx, cty[i])
                if cty[i] == 0:
                    dvn[i] = 0
            if cty[i] == 0:
                nsc += 1
            elif cty[i] == 1:
                nta1 += 1
            elif cty[i] == 3:
                nta2 += 1
            else:
                nll += 1

        if init_mode == 3:
            traj[step, 0] = t
            traj[step, 1] = x1a[0]
            traj[step, 2] = x2a[0]
            t += dt
            continue

        # ---------------- phase 2: fate decisions (SystemUpdate)
        m = 0  # new pool size
        for i in range(n):
            r1 = np.random.rand()
            r2 = np.random.rand()
            c = cty[i]
            if c == 0:
                kap = kappa0
                bet = beta0 / (1.0 + nsc / theta)
                mu = mu0
                tau = tau0
            else:  # TA1/TA2 (label can only be 0,1,3 here)
                kap = 2.0 * kappa0
                bet = beta0
                mu = mu0
                tau = 0.5 * tau0
            if stg[i] == 0 or stg[i] == 3:
                if r1 < kap * dt:
                    continue  # die
                if r1 < kap * dt + bet * dt:
                    # enter proliferating: stage 1, age 0; log transition
                    if t > log_after:
                        trans[_trow(pty[i]), _trow(c)] += 1
                    nx1[m] = x1a[i]; nx2[m] = x2a[i]
                    nu1[m] = u1a[i]; nu2[m] = u2a[i]
                    ne1[m] = e1a[i]; ne2[m] = e2a[i]
                    nct[m] = cty[i]; npt[m] = pty[i]
                    nst_[m] = 1; ndv[m] = dvn[i]; nag[m] = 0.0
                    m += 1
                else:
                    nx1[m] = x1a[i]; nx2[m] = x2a[i]
                    nu1[m] = u1a[i]; nu2[m] = u2a[i]
                    ne1[m] = e1a[i]; ne2[m] = e2a[i]
                    nct[m] = cty[i]; npt[m] = pty[i]
                    nst_[m] = 0; ndv[m] = dvn[i]; nag[m] = age[i]
                    m += 1
            else:
                if r2 < mu * dt:
                    continue  # die
                if age[i] < tau:
                    nx1[m] = x1a[i]; nx2[m] = x2a[i]
                    nu1[m] = u1a[i]; nu2[m] = u2a[i]
                    ne1[m] = e1a[i]; ne2[m] = e2a[i]
                    nct[m] = cty[i]; npt[m] = pty[i]
                    nst_[m] = 2; ndv[m] = dvn[i]; nag[m] = age[i] + dt
                    m += 1
                else:
                    # ---- division
                    ndiv += 1
                    ddv = dvn[i] + 1
                    if ddv >= maxdiv and (cty[i] == 1 or cty[i] == 3):
                        continue  # both daughters die (TA at max divisions)
                    if fixed_u == 1:
                        d1u1 = u1a[i]; d1u2 = u2a[i]
                        d2u1 = u1a[i]; d2u2 = u2a[i]
                    else:
                        for g in range(2):
                            v = u1a[i] if g == 0 else u2a[i]
                            ph = m1 + m2_ * (m3 * v) ** m4 / (1.0 + (m3 * v) ** m4)
                            ga = m0 * ph
                            gb = m0 * (1.0 - ph)
                            if g == 0:
                                d1u1 = _beta_u(ga, gb)
                                d2u1 = _beta_u(ga, gb)
                            else:
                                d1u2 = _beta_u(ga, gb)
                                d2u2 = _beta_u(ga, gb)
                    # daughter 1
                    nx1[m] = x1a[i]; nx2[m] = x2a[i]
                    nu1[m] = d1u1; nu2[m] = d1u2
                    ne1[m] = e1a[i]; ne2[m] = e2a[i]
                    nct[m] = cty[i]; npt[m] = cty[i]
                    nst_[m] = 3; ndv[m] = ddv; nag[m] = age[i]
                    m += 1
                    # daughter 2
                    nx1[m] = x1a[i]; nx2[m] = x2a[i]
                    nu1[m] = d2u1; nu2[m] = d2u2
                    ne1[m] = e1a[i]; ne2[m] = e2a[i]
                    nct[m] = cty[i]; npt[m] = cty[i]
                    nst_[m] = 3; ndv[m] = ddv; nag[m] = age[i]
                    m += 1

        # swap pools (truncate at MAXCELL like the C++ copy loop)
        if m > MAXCELL:
            trunc_events += 1
            m = MAXCELL
        n = m
        for i in range(n):
            x1a[i] = nx1[i]; x2a[i] = nx2[i]
            u1a[i] = nu1[i]; u2a[i] = nu2[i]
            e1a[i] = ne1[i]; e2a[i] = ne2[i]
            cty[i] = nct[i]; pty[i] = npt[i]
            stg[i] = nst_[i]; dvn[i] = ndv[i]; age[i] = nag[i]

        t += dt

        # hourly records
        if step % 4 == 3 and hour_i < nHrec:
            counts[hour_i, 0] = t
            counts[hour_i, 1] = n
            counts[hour_i, 2] = nsc
            counts[hour_i, 3] = nta1
            counts[hour_i, 4] = nta2
            counts[hour_i, 5] = nll
            counts[hour_i, 6] = ndiv
            hour_i += 1

        # snapshots
        while snap_i < nS and t >= snap_times[snap_i]:
            k = n if n < SNAP_CAP else SNAP_CAP
            for i in range(k):
                snaps[snap_i, i, 0] = x1a[i]
                snaps[snap_i, i, 1] = x2a[i]
                snaps[snap_i, i, 2] = u1a[i]
                snaps[snap_i, i, 3] = u2a[i]
            snap_n[snap_i] = k
            snap_i += 1

        if stop_sc > 0 and nsc >= stop_sc:
            t_success = t
            break
        if n == 0:
            t_ext = t
            break

    scal[0] = ndiv
    scal[1] = t_success
    scal[2] = t_ext
    scal[3] = n
    scal[4] = nsc
    scal[5] = t
    scal[6] = trunc_events
    return counts, trans, snaps, snap_n, traj, scal


# ---------------------------------------------------------------------------
# Python-side helpers
# ---------------------------------------------------------------------------

def dxdt(x1, x2, a1, a2, P, aTF=0.0):
    """Deterministic RHS of Eq (1) (for steady states / bifurcation analysis)."""
    (b1, b2, k1, k2, s1, s2, rho1, rho2, sg1, sg2, z1, z2,
     al1, al2, l1, l2, nH, m0, m1, m2_, m3, m4,
     theta, beta0, kappa0, mu0, tau0, dt, aTF_, maxdiv, thx) = P
    s1n = s1 ** nH; s2n = s2 ** nH
    p1 = x1 ** nH; p2 = x2 ** nH
    h1 = p1 / (s1n + p1); h2 = p2 / (s2n + p2)
    r1 = s2n / (s2n + p2); r2 = s1n / (s1n + p1)
    dx1 = (a1 + aTF) * (rho1 + (1 - rho1) * h1) + b1 * r1 - k1 * x1
    dx2 = a2 * (rho2 + (1 - rho2) * h2) + b2 * r2 - k2 * x2
    return dx1, dx2


def steady_x(a1, a2, P, ic=(2.0, 2.0), T=800.0, dt=0.02, aTF=0.0):
    """Relax deterministic ODE to steady state from ic."""
    x1, x2 = ic
    for s in range(int(T / dt)):
        d1, d2 = dxdt(x1, x2, a1, a2, P, aTF)
        x1 += d1 * dt
        x2 += d2 * dt
        if x1 < 0: x1 = 0.0
        if x2 < 0: x2 = 0.0
    return x1, x2


def a_of_u(u, alpha=0.4, lam=1.9):
    return alpha * np.exp(lam * u)


def phi(v, P):
    """Inheritance expectation map phi(v)."""
    (b1, b2, k1, k2, s1, s2, r1_, r2_, sg1, sg2, z1, z2,
     al1, al2, l1, l2, nH, m0, m1, m2_, m3, m4,
     th, b0, k0, mu0, t0, dt, aTF, md, thx) = P
    return m1 + m2_ * (m3 * v) ** m4 / (1.0 + (m3 * v) ** m4)


def phi_fixed_points(P, lo0=0.15, hi0=0.5):
    u = lo0
    for _ in range(300):
        u = phi(u, P)
    lo = u
    u = hi0
    for _ in range(300):
        u = phi(u, P)
    return lo, u


def run_population(seed, T, params=None, snap_times=None, init='random',
                   init_state=None, fixed_u=False, stop_sc=0, n_init=100,
                   log_after=2000.0):
    """Run one population simulation. init: 'random' | 'given'."""
    P = params if params is not None else paper_params()
    if snap_times is None:
        snap_times = np.zeros(0)
    else:
        snap_times = np.asarray(snap_times, np.float64)
    if init == 'random':
        iv = np.zeros(6)
        mode = 0
    else:
        x1, x2, u1, u2 = init_state
        iv = np.array([x1, x2, u1, u2, n_init, 0.0])
        mode = 1
    return run_sim(seed, P, float(T), mode, iv, 1 if fixed_u else 0,
                   float(stop_sc), snap_times, float(log_after))


def run_single(seed, T, x1, x2, u1, u2, params=None):
    """Single-cell SDE trajectory (no division/death). Fig 2C/D."""
    P = params if params is not None else paper_params()
    iv = np.array([x1, x2, u1, u2, 1.0, 0.0])
    return run_sim(seed, P, float(T), 3, iv, 0, 0.0, np.zeros(0), 0.0)


# ---------------------------------------------------------------------------
# Table 3 validation (MATLAB DataAnalysis.m protocol)
# ---------------------------------------------------------------------------

def validate_table3(seeds, T=4000.0, m2=0.8):
    """Per-run row-normalized transition matrices, averaged across runs."""
    names = ['SC', 'TA1', 'TA2']
    P = paper_params(m2=m2)
    per_run = []
    counts_last = None
    for sd in seeds:
        counts, trans, _, _, _, scal = run_population(
            sd, T, params=P, log_after=2000.0)
        rs = trans.sum(axis=1)
        rs[rs == 0] = 1
        per_run.append(100.0 * trans / rs[:, None])
        counts_last = counts
    M = np.mean(per_run, axis=0)
    print(f"Table 3 validation: {len(seeds)} runs, T={T}h, m2={m2}, "
          f"t>2000h events")
    print("parent\\current  " + "  ".join(f"{n:>9}" for n in names))
    for r in range(3):
        print(f"  {names[r]:>6}       " + "  ".join(f"{M[r, c]:9.4f}" for c in range(3)))
    print("Paper:  SC [73.9402, 12.9889, 13.0709]; TA1 [0.0013, 99.9955, 0.0032];"
          " TA2 [0.0013, 0.0035, 99.9952]")
    w = counts_last[counts_last[:, 0] >= 2000]
    print(f"last seed means t in [2000,4000]: total={w[:, 1].mean():.0f} "
          f"SC={w[:, 2].mean():.0f} TA1={w[:, 3].mean():.0f} TA2={w[:, 4].mean():.0f}"
          f"  (C++ single-run ground truth: 6984/177/3269/3537)")
    return M


if __name__ == '__main__':
    import sys
    seeds = [int(s) for s in sys.argv[1:]] or [1, 2, 3]
    validate_table3(seeds)
