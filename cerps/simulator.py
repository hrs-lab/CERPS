"""CERPS simulator: couples the GRN, epigenetics, division, continuous
rates and 3D space into one individual-based model.

Event structure follows the G0/G1 cell-cycle scheme of the paper (and its
C++ implementation): resting cells enter proliferation with probability
beta*dt or are removed with probability kappa*dt; proliferating cells
apoptose with probability mu*dt, grow for tau, then divide into two
daughters.  The intracellular state (x: gene expression, u: epigenetic
marks, eta: extrinsic noise) evolves every step for every cell.

Requirement mapping:
R1 intracellular epigenetic dynamics + multiple kernels  -> epi.continuous,
   EpiConfig.kernel in {beta, gaussian, cpg, deterministic}
R2 expression changes at division                        -> DivisionConfig
R3 multi-gene GRN                                        -> GRNConfig
R4 continuous, state-dependent rates (no classification) -> RatesConfig
   rate_mode="continuous"
R5 bidirectional coupling: epigenetics <-> GRN <-> population
   (epi -> GRN: a_i depends on u_i;  GRN -> epi: continuous dynamics;
    population -> cell: niche feedback Q / local density / paracrine;
    cell -> population: rates depend on state)
R6 3D space                                              -> SpaceConfig
"""
import numpy as np

from .config import CERPSConfig
from .grn import GRN
from . import epigenetics as epi
from . import division as div_mod
from . import rates as rates_mod
from .spatial import Space

REST, JUST_PROLIF, PROLIF, NEWBORN = 0, 1, 2, 3


class Simulator:
    def __init__(self, cfg: CERPSConfig, init=None):
        self.cfg = cfg
        self.grn = GRN(cfg.grn)
        self.rng = np.random.default_rng(cfg.seed)
        n = cfg.n_init
        G = cfg.grn.n_genes

        if init is None:
            self.x = self.rng.uniform(0, 3, size=(n, G))
            self.u = self.rng.uniform(0, 1, size=(n, G))
        else:
            self.x = np.asarray(init["x"], float).copy()
            self.u = np.asarray(init["u"], float).copy()
            n = len(self.x)
        self.n = n
        self.eta = np.zeros((n, G))
        self.stage = np.zeros(n, np.int8)
        self.age = np.zeros(n)
        self.div_count = np.zeros(n)      # continuous senescence counter
        self.div_int = np.zeros(n, np.int64)   # integer divisions (record)
        # legacy classification labels (also used for outputs)
        self.label = rates_mod.classify_legacy(
            self.grn, self.x, np.zeros(n, np.int8),
            cfg.rates.x_threshold)
        self.parent_label = np.zeros(n, np.int8)
        # P6 non-Markovian memory: per-cell leaky integrator state (init 0)
        self.mem = np.zeros(n)

        # space
        self.space = None
        if cfg.space.enabled:
            self.space = Space(cfg.space)
            self.space.init_positions(n, self.rng)
            if cfg.space.paracrine:
                self.grn.set_paracrine_weight(cfg.space.response)

        # records
        self.records = dict(t=[], total=[], sc=[], ta1=[], ta2=[],
                            other=[], mean_u=[], mean_S=[])
        self.snapshots = []
        self.transitions = []            # (parent_label, label) at entry
        self.max_cells = 0
        # niche-signal history for the P8 feedback-time delay; hist[s] is the
        # (soft) SC count at step s, used by effective_niche_signal()
        self._q_hist = []
        # resolved stemness-switch threshold (None -> derived from the
        # legacy classification boundary; see rates.resolve_K)
        self.K_switch = rates_mod.resolve_K(cfg.rates, self.grn)

    # ------------------------------------------------------------------
    def _record(self, t, stem_frac=None):
        S = rates_mod.stemness(self.grn, self.x) if stem_frac is None \
            else stem_frac
        r = self.records
        r["t"].append(t)
        r["total"].append(self.n)
        r["sc"].append(int((self.label == 0).sum()))
        r["ta1"].append(int((self.label == 1).sum()))
        r["ta2"].append(int((self.label == 3).sum()))
        r["other"].append(int((self.label == 2).sum()) if
                          (self.label == 2).any() else 0)
        r["mean_u"].append(self.u.mean(axis=0).tolist())
        r["mean_S"].append(float(S.mean()))

    # ------------------------------------------------------------------
    def run(self, t_max=None, snapshot_times=()):
        cfg = self.cfg
        dt = cfg.dt
        t_max = cfg.max_time if t_max is None else t_max
        n_steps = int(round(t_max / dt))
        snap_steps = {min(int(round(ts / dt)), n_steps - 1): ts
                      for ts in snapshot_times}
        legacy = cfg.rates.rate_mode == "legacy"
        contact = cfg.space.enabled and cfg.space.contact_inhibition

        for step in range(n_steps):
            t = step * dt
            n = self.n
            resting = (self.stage == REST) | (self.stage == NEWBORN)

            # ---- paracrine field (population -> intracellular) ---------
            paracrine = None
            if self.space is not None and cfg.space.paracrine:
                S = rates_mod.stemness(self.grn, self.x)
                secretion = cfg.space.secretion * S
                paracrine = self.space.field_step(secretion, dt)

            # ---- intracellular dynamics --------------------------------
            self.eta = self.grn.ou_step(self.eta, dt, self.rng)
            self.x = self.grn.step(self.x, self.u, self.eta, dt, paracrine)
            act = self.grn.activation_levels(self.x)
            if cfg.epi.continuous:
                self.u = epi.continuous_step(self.u, act, cfg.epi, dt)

            # ---- phenotype labels -----------------------------------------
            # Protocol (paper Sec. 2.3 numerical scheme, C++ UpdateCellType,
            # hsl2024.reference_sim): the phenotype is re-evaluated ONLY
            # while the cell is in the resting phase (stages 0/3); it is
            # frozen from proliferation entry to mitosis.  Updating it for
            # proliferating cells as well changes the SC count Q during the
            # cycle, which biases beta_SC = beta0*theta/(theta+Q) upward and
            # inflates SC self-renewal by ~2.4 pp (see audit report 4.5).
            old_label = self.label
            new_label = rates_mod.classify_legacy(
                self.grn, self.x, self.label, cfg.rates.x_threshold)
            self.label = np.where(resting, new_label, self.label)
            if legacy:
                # C++ rule: a RESTING cell classified as SC resets its
                # division counter (dedifferentiation restores potential)
                sc_rest = resting & (self.label == 0)
                self.div_int[sc_rest] = 0
                self.div_count[sc_rest] = 0.0

            # ---- transition records at proliferation entry -------------
            # rec_tran_from discards the transient (paper protocol: only
            # Step > 8000, i.e. t > 2000 h of a 4000 h run, is counted).
            jp = self.stage == JUST_PROLIF
            if jp.any() and t >= cfg.rec_tran_from:
                for pl, cl in zip(self.parent_label[jp], self.label[jp]):
                    self.transitions.append((int(pl), int(cl)))

            # ---- P6 non-Markovian memory: per-cell leaky integrator -------
            # drv in [0,1]; mem converges to the running mean of the driver
            # with timescale tau, giving a genuinely history-dependent state.
            mem_arg = None
            mem_couple = 0.0
            if cfg.memory.enabled:
                drv = rates_mod.stemness(self.grn, self.x)  # per-cell [0,1]
                self.mem = self.mem + (dt / cfg.memory.tau) * (drv - self.mem)
                mem_arg = self.mem
                mem_couple = cfg.memory.mem_couple

            # ---- rates ---------------------------------------------------
            # compute the *current* niche signal for both modes, record it in
            # the history, then feed the DELAYED value to the kinetics (P8).
            if legacy:
                q_now = float((self.label == 0).sum())
            else:
                S = rates_mod.stemness(self.grn, self.x)
                q_now = float(rates_mod.hill(
                    S, self.K_switch, cfg.rates.p_switch).sum())
            self._q_hist.append(q_now)
            q_eff = rates_mod.effective_niche_signal(
                self._q_hist, step, cfg.rates.delay_steps, q_now)

            if legacy:
                beta, kappa, mu, tau = rates_mod.legacy_rates(
                    self.label, q_eff, cfg.rates, mem_arg, mem_couple)
                senesce = (self.div_int + 1 >= cfg.rates.max_div_ta) & \
                    (self.label != 0)
            else:
                rho_local = self.space.local_density() if contact else None
                beta, kappa, mu, tau, senesce = rates_mod.rates(
                    self.grn, self.x, self.u, cfg.rates, self.div_count,
                    q_eff, rho_local, mem_arg, mem_couple)

            # ---- fate decisions ------------------------------------------
            r1 = self.rng.random(n)
            r2 = self.rng.random(n)
            keep = np.ones(n, bool)
            die = resting & (r1 < kappa * dt)
            enter = resting & ~die & (r1 < (kappa + beta) * dt)
            self.stage[enter] = JUST_PROLIF
            self.age[enter] = 0.0
            self.stage[resting & ~die & ~enter] = REST
            keep &= ~die
            die2 = (~resting) & (r2 < mu * dt)
            keep &= ~die2
            grow = (~resting) & ~die2 & (self.age < tau)
            self.stage[grow] = PROLIF
            self.age[grow] += dt
            mitosis = (~resting) & ~die2 & ~grow
            keep &= ~die2

            # senescence at mitosis
            sen_die = mitosis & senesce
            keep &= ~sen_die
            mitosis &= ~sen_die

            # ---- mitosis: partitioning + epigenetic inheritance --------
            idx = np.flatnonzero(mitosis)
            if len(idx):
                self.div_int[idx] += 1
                self.parent_label[idx] = self.label[idx]
                self.stage[idx] = NEWBORN
                if not legacy:
                    HS = rates_mod.hill(
                        rates_mod.stemness(self.grn, self.x[idx]),
                        self.K_switch, cfg.rates.p_switch)
                    self.div_count[idx] += (1.0 - HS)
                    # resting decay toward zero for stem-like cells
                x1, x2 = div_mod.partition(self.x[idx], cfg.division,
                                           self.rng)
                u1, u2 = epi.inherit(self.u[idx], cfg.epi, self.rng)
                self.x[idx] = x1
                self.u[idx] = u1
                new = dict(
                    x=x2, u=u2, eta=self.eta[idx].copy(),
                    stage=np.full(len(idx), NEWBORN, np.int8),
                    age=self.age[idx].copy(),
                    div_int=self.div_int[idx].copy(),
                    div_count=self.div_count[idx].copy(),
                    label=self.label[idx].copy(),
                    parent_label=self.parent_label[idx].copy(),
                    mem=self.mem[idx].copy())
                for key, val in new.items():
                    setattr(self, key, np.concatenate([getattr(self, key),
                                                       val]))
                if self.space is not None:
                    p1, p2 = self.space.division_positions(idx, self.rng)
                    self.space.pos[idx] = p1
                    self.space.add(p2)
                keep = np.concatenate([keep, np.ones(len(idx), bool)])

            # div_count decay while resting (non-legacy): stem-like cells
            # continuously restore their replicative potential
            if not legacy:
                HSr = rates_mod.hill(
                    rates_mod.stemness(self.grn, self.x),
                    self.K_switch, cfg.rates.p_switch)
                resting_now = (self.stage == REST) | \
                    (self.stage == NEWBORN)
                self.div_count[resting_now] *= (1.0 - HSr[resting_now])

            # ---- motility + mechanics ------------------------------------
            if self.space is not None:
                if cfg.space.motility > 0:
                    self.space.pos += cfg.space.motility * np.sqrt(dt) * \
                        self.rng.standard_normal(self.space.pos.shape)
                self.space.relax(dt)

            # ---- compaction ----------------------------------------------
            if not keep.all():
                self.x = self.x[keep]
                self.u = self.u[keep]
                self.eta = self.eta[keep]
                self.stage = self.stage[keep]
                self.age = self.age[keep]
                self.div_count = self.div_count[keep]
                self.div_int = self.div_int[keep]
                self.label = self.label[keep]
                self.parent_label = self.parent_label[keep]
                self.mem = self.mem[keep]
                if self.space is not None:
                    self.space.remove(keep)
            self.n = len(self.x)
            self.max_cells = max(self.max_cells, self.n)

            # ---- records ---------------------------------------------------
            if step % cfg.record_stride == 0 or step == n_steps - 1:
                self._record(t + dt)
            if step in snap_steps:
                self.snapshots.append(dict(
                    t=snap_steps[step], x=self.x.copy(), u=self.u.copy(),
                    pos=None if self.space is None
                    else self.space.pos.copy()))
            if self.n == 0:
                break
        return self.records
