"""Readable NumPy-vectorised reference implementation of the hybrid model.

Independent re-implementation of the same protocol as ``numba_sim.py``
(one pass per step, pure NumPy).  Purposes:

1. cross-validation of the fast numba simulator (statistical agreement),
2. lineage tracking for the transition trajectories of Fig. 6
   (ring-buffer histories of candidate cells; candidate = a daughter whose
   inherited epigenetic state exceeds ``u_track`` in either gene, which is
   a necessary condition for de-/trans-differentiation at sigma = 0.05).

~50x slower than the numba version; use for small jobs and validation.
"""
import numpy as np

from .params import Params

REST, JUST_PROLIF, PROLIF, NEWBORN = 0, 1, 2, 3
_MAP3 = np.array([0, 1, 0, 2])       # type code {0,1,3} -> matrix index


def _phi(v, m1, m2, m3, m4):
    g = (m3 * v) ** m4
    return m1 + m2 * g / (1.0 + g)


def _classify(x, prev):
    t = prev.copy()
    t[(x[:, 0] > 0.5) & (x[:, 1] > 0.5)] = 0
    t[(x[:, 0] > 0.5) & (x[:, 1] < 0.5)] = 1
    t[(x[:, 0] < 0.5) & (x[:, 1] > 0.5)] = 3
    return t


def simulate_reference(p: Params, n_runs=1, n_steps=16000, seed=0,
                       init=None, fixed_u=False, rec_tran_from=8000,
                       u_track=None, ring_len=2400, max_track_slots=1500):
    """Simulate; returns dict(sys_ts, tran, events, track_meta).

    u_track : threshold making a newborn cell a tracking candidate (Fig. 6).
    events  : list of dicts (uid, old, new, t, x, u, history ndarray
              (K,4): t, x1, x2, ctype) for de-/trans-differentiation events.
    """
    rng = np.random.default_rng(seed)
    dt = p.dt
    ring_len = int(ring_len)

    if init is None:
        n = p.n_init * n_runs
        x = rng.uniform(0, 3, size=(n, 2))
        u = rng.uniform(0, 1, size=(n, 2))
        run = np.repeat(np.arange(n_runs), p.n_init).astype(np.int32)
    else:
        x = np.asarray(init["x"], float).copy()
        u = np.asarray(init["u"], float).copy()
        run = np.asarray(init["run"], np.int32).copy()

    n = len(run)
    eta = np.zeros((n, 2))
    stage = np.zeros(n, np.int8)
    age = np.zeros(n)
    div = np.zeros(n, np.int64)
    ctype = _classify(x, np.zeros(n, np.int8))
    ptype = np.zeros(n, np.int8)

    sys_ts = np.zeros((n_steps, n_runs, 4), np.int32)
    tran = np.zeros((n_runs, 3, 3), np.int64)

    # ---- tracking machinery ------------------------------------------------
    track_on = u_track is not None
    slot = np.full(n, -1, np.int32)          # tracking slot per cell
    ring = np.zeros((max_track_slots, ring_len, 4), np.float32)
    ring_pos = np.zeros(max_track_slots, np.int32)
    ring_uid = np.full(max_track_slots, -1, np.int64)
    free_slots = list(range(max_track_slots))
    events = []

    def _try_track(gidx, u_arr):
        """Assign free tracking slots to candidate cells (indices gidx)."""
        if not track_on or not free_slots:
            return
        cand = gidx[(u_arr[:, 0] > u_track) | (u_arr[:, 1] > u_track)]
        for i in cand:
            if not free_slots:
                break
            s = free_slots.pop()
            slot[i] = s
            ring_pos[s] = 0
            ring_uid[s] = uid[i]

    uid = np.arange(n, dtype=np.int64)
    next_uid = n
    if track_on:
        _try_track(np.arange(n), u)

    for step in range(n_steps):
        n = len(run)
        t_now = step * dt
        # ---- gene dynamics -------------------------------------------------
        eta = eta - eta * (dt / p.zeta1) \
            + np.sqrt(2 * dt / p.zeta1) * rng.standard_normal((n, 2))
        a = np.empty((n, 2))
        a[:, 0] = p.alpha1 * np.exp(p.lam1 * u[:, 0]) * \
            np.exp(p.sigma * eta[:, 0] - 0.5 * p.sigma ** 2) + p.a_tf
        a[:, 1] = p.alpha2 * np.exp(p.lam2 * u[:, 1]) * \
            np.exp(p.sigma * eta[:, 1] - 0.5 * p.sigma ** 2)
        s1n = p.s1 ** p.n
        s2n = p.s2 ** p.n
        x1p = x[:, 0] ** p.n
        x2p = x[:, 1] ** p.n
        f1 = a[:, 0] * (p.rho1 + (1 - p.rho1) * x1p / (s1n + x1p)) \
            + p.b1 * s2n / (s2n + x2p) - p.k1 * x[:, 0]
        f2 = a[:, 1] * (p.rho2 + (1 - p.rho2) * x2p / (s2n + x2p)) \
            + p.b2 * s1n / (s1n + x1p) - p.k2 * x[:, 1]
        x = x + np.stack([f1, f2], 1) * dt

        resting = (stage == REST) | (stage == NEWBORN)
        old_type = ctype
        ctype = _classify(x, ctype)
        div[resting & (ctype == 0)] = 0

        # ---- transition records (vectorised) -------------------------------
        jp = (stage == JUST_PROLIF) & (step + 1 > rec_tran_from)
        if jp.any():
            pj = _MAP3[ptype[jp]]
            cj = _MAP3[ctype[jp]]
            flat = run[jp] * 9 + pj * 3 + cj
            tran += np.bincount(flat, minlength=n_runs * 9).reshape(n_runs, 3, 3)
        if track_on:
            # net mother->daughter transitions of tracked candidate cells
            # (Table-3 semantics: type at the birth division vs type at
            # proliferation entry; a transient SC passage in between is
            # invisible to this definition, exactly as in Table 3)
            jtr = (stage == JUST_PROLIF) & (ptype != ctype) & (slot >= 0)
            for i in np.flatnonzero(jtr):
                o, nw = int(ptype[i]), int(ctype[i])
                if o == 0:      # SC-born daughters differentiating are not
                    continue    # de-/trans-differentiation events
                s = slot[i]
                k = ring_pos[s]
                hist = np.concatenate([ring[s, k:], ring[s, :k]])
                events.append(dict(uid=int(uid[i]), old=o, new=nw,
                                   t=t_now, x=x[i].copy(), u=u[i].copy(),
                                   history=hist.copy()))

        # ---- tracking histories --------------------------------------------
        if track_on:
            for i in np.flatnonzero(slot >= 0):
                s = slot[i]
                k = ring_pos[s]
                ring[s, k, 0] = t_now
                ring[s, k, 1] = x[i, 0]
                ring[s, k, 2] = x[i, 1]
                ring[s, k, 3] = ctype[i]
                ring_pos[s] = (k + 1) % ring_len

        # ---- counting (before fate decisions) ------------------------------
        for r in range(n_runs):
            m = run == r
            sys_ts[step, r, 0] = m.sum()
            sys_ts[step, r, 1] = ((ctype == 0) & m).sum()
            sys_ts[step, r, 2] = ((ctype == 1) & m).sum()
            sys_ts[step, r, 3] = ((ctype == 3) & m).sum()
        q_sc = sys_ts[step, :, 1].astype(float)

        # ---- fate decisions --------------------------------------------------
        keep = np.ones(n, bool)
        beta = np.where(ctype == 0, p.beta0 / (1.0 + q_sc[run] / p.theta0),
                        p.beta0)
        kap = np.where(ctype == 0, p.kappa0, 2 * p.kappa0)
        tau = np.where(ctype == 0, p.tau0, 0.5 * p.tau0)
        r1 = rng.random(n)
        r2 = rng.random(n)

        die = resting & (r1 < kap * dt)
        enter = resting & ~die & (r1 < (kap + beta) * dt)
        stage[enter] = JUST_PROLIF
        age[enter] = 0.0
        stage[resting & ~die & ~enter] = REST
        keep &= ~die

        die2 = (~resting) & (r2 < p.mu0 * dt)
        grow = (~resting) & ~die2 & (age < tau)
        stage[grow] = PROLIF
        age[grow] += dt
        mitosis = (~resting) & ~die2 & ~grow
        keep &= ~die2

        senesce = mitosis & ((div + 1) == p.max_div_ta) & (ctype != 0)
        keep &= ~senesce
        mitosis &= ~senesce

        div[mitosis] += 1
        ptype[mitosis] = ctype[mitosis]
        stage[mitosis] = NEWBORN

        n_m = int(mitosis.sum())
        if n_m:
            idx = np.flatnonzero(mitosis)
            if fixed_u:
                u_d1 = u[idx]
                u_d2 = u[idx]
            else:
                phi = _phi(u[idx], p.m1, p.m2, p.m3, p.m4)
                u_d1 = rng.beta(p.m0 * phi, p.m0 * (1 - phi))
                u_d2 = rng.beta(p.m0 * phi, p.m0 * (1 - phi))
                u[idx] = u_d1
            d = dict(
                x=x[idx], u=u_d2, eta=eta[idx], run=run[idx], age=age[idx],
                div=div[idx], ctype=ctype[idx], ptype=ptype[idx],
                uid=np.arange(next_uid, next_uid + n_m, dtype=np.int64))
            next_uid += n_m
            d_stage = np.full(n_m, NEWBORN, np.int8)
            d_slot = np.full(n_m, -1, np.int32)
            x = np.vstack([x, d["x"]])
            u = np.vstack([u, d["u"]])
            eta = np.vstack([eta, d["eta"]])
            run = np.concatenate([run, d["run"]])
            stage = np.concatenate([stage, d_stage])
            age = np.concatenate([age, d["age"]])
            div = np.concatenate([div, d["div"]])
            ctype = np.concatenate([ctype, d["ctype"]])
            ptype = np.concatenate([ptype, d["ptype"]])
            uid = np.concatenate([uid, d["uid"]])
            slot = np.concatenate([slot, d_slot])
            keep = np.concatenate([keep, np.ones(n_m, bool)])
            if track_on:
                _try_track(np.arange(len(run) - n_m, len(run)), u[-n_m:])

        # ---- compaction ------------------------------------------------------
        if not keep.all():
            dead = np.flatnonzero(~keep)
            if track_on:
                for i in dead:
                    if slot[i] >= 0:
                        free_slots.append(int(slot[i]))
                        ring_uid[slot[i]] = -1
            x = x[keep]
            u = u[keep]
            eta = eta[keep]
            run = run[keep]
            stage = stage[keep]
            age = age[keep]
            div = div[keep]
            ctype = ctype[keep]
            ptype = ptype[keep]
            uid = uid[keep]
            slot = slot[keep]
        if len(run) == 0:
            sys_ts[step + 1:] = 0
            break

    return dict(sys_ts=sys_ts, tran=tran, events=events,
                final=dict(x=x, u=u, run=run, stage=stage, age=age,
                           div=div, ctype=ctype, ptype=ptype, eta=eta))
