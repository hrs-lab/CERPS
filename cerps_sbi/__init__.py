"""
cerps_sbi — a thin, HONEST adapter layer over the real CERPS kernel.

Why this exists
---------------
The paper-workflow methodology docs (docs/paper-workflow/*.md) describe the
calibration layer with conceptual class names (legacy_config(), Simulator,
RatesConfig, MemoryConfig, EpiConfig, GRNConfig, delay_steps). The *actual*
CERPS code in cerps_model/src/cerps.py exposes only:

    from cerps import run_cerps, CERPS_DEFAULTS
    run_cerps(seed, T, n_init, init, init_x, init_u, init_pos,
              record_after, rec_every, snapshot_times, **params)

This module re-exposes the conceptual API on top of the real kernel so the
methodology pseudocode is executable, and — crucially — records where the
real code diverges from the docs (see HONESTY NOTES below).

HONESTY NOTES (project "honesty boundary" — do NOT fake numerical results)
--------------------------------------------------------------------------------
* There is NO `legacy_config()` class. The calibrated anchor = CERPS_DEFAULTS
  merged with the LOCKED manuscript spatial config BASE (the config used by
  analysis/src/abc.py and analysis/src/mechanism_scan.py).
* P8 (feedback delay) is NOT implemented **in this `cerps_model` kernel** — there
  is no `delay_steps` here, so any delay axis must be simulated elsewhere.
  Scope note (2026-10-20): the class-based `cerps/` framework in the same project
  DOES implement P8 via `RatesConfig.delay_steps` + `rates.effective_niche_signal()`.
  Candidates relying on a delay primitive (C1) are therefore testable in `cerps/`;
  the sweeps produced by THIS adapter run on the `cerps_model` kernel and cannot
  exercise delay. See docs/paper-workflow/api-clarification.md.
* Memory (P6) = the `epa` (epigenetic-age) accumulator, gated by `kappa_a`
  (accelerates death with divisional age) and `div_death_age`.
* Partitioning (P3) = `s_part` (expression partitioning noise at mitosis).
* Epigenetics (P5) = `u` + `k_maint` + `w_u=(w0,wG,wP,wK,wC)` + `D_u`.
* Niche (P1) = the S-block: S0, n_niche, rho_S, delta_S, D_S, uptake, KS, hS,
  h_maint_S, h_ref_S, well_mixed.
* GRN (P4 / mutual repression) = alpha, b (repression weights), rho, k, s,
  n_hill, sigma, zeta, lam.
"""

import os
import sys
import numpy as np

# Robustly locate the real CERPS kernel (cerps_model/src/cerps.py), which is
# NOT part of this adapter. Search order:
#   1. $CERPS_SRC env var (explicit)
#   2. sibling cerps_model/src (when this adapter sits beside cerps_model, as in
#      the upstream hrs-lab/CERPS repo)
#   3. the local CERPS_complete_package checkout used to develop this adapter
HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(HERE)
_CANDIDATES = [
    os.environ.get("CERPS_SRC"),
    os.path.join(REPO_ROOT, "cerps_model", "src"),
    "C:/Users/黄荣生/Downloads/CERPS_complete_package/cerps_model/src",
]
for _c in _CANDIDATES:
    if _c and os.path.isfile(os.path.join(_c, "cerps.py")):
        if _c not in sys.path:
            sys.path.insert(0, _c)
        break

from cerps import run_cerps, CERPS_DEFAULTS  # noqa: E402

# Locked manuscript spatial config (analysis/src/{abc,mechanism_scan}.py).
LOCKED_BASELINE = dict(crowding=1.5, L=8.0, grid=16, n_niche=25, rho_S=2.5)

# ABC-SMC priors (analysis/src/abc.py) — documents the calibrated ranges.
ABC_PRIORS = [
    ("r_eff", 0.45, 0.90), ("r_max", 0.85, 0.99), ("h_ref_S", 0.05, 0.80),
    ("k_maint", 0.01, 0.15), ("s_div", 0.01, 0.08), ("beta0", 0.15, 0.60),
    ("sigma", 0.01, 0.08), ("D_u", 1e-5, 8e-5),
]

NG = 4  # GATA1, SPI1, KLF1, CEBPA


def legacy_config():
    """Return the calibrated anchor config (real: DEFAULTS + locked BASE)."""
    cfg = dict(CERPS_DEFAULTS)
    cfg.update(dict(LOCKED_BASELINE))
    return cfg


def _snapshot_init(n_init, L, seed=3, rng0=None):
    """Build init_x / init_u / init_pos arrays (matches mechanism_scan.py)."""
    if rng0 is None:
        rng0 = np.random.default_rng(0)
    n = n_init
    IX = np.tile(np.array([1.04, 1.04, 0.18, 0.18]), (n, 1)) \
        * (1 + 0.05 * rng0.standard_normal((n, NG)))
    IU = np.tile(np.array([0.50, 0.50, 0.73, 0.73]), (n, 1)) \
        * (1 + 0.05 * rng0.standard_normal((n, NG)))
    IU = np.clip(IU, 0.05, 0.95)
    ip = np.random.default_rng(seed).uniform(0, L, (n, 2))
    return (np.ascontiguousarray(IX, dtype=np.float64),
            np.ascontiguousarray(IU, dtype=np.float64),
            np.ascontiguousarray(ip, dtype=np.float64))


def run_one(override=None, seed=3, T=4000.0, n_init=600,
            snapshot_times=(), record_after=0.7, rec_every=4):
    """Run one CERPS realization on top of the legacy anchor, return metrics.

    `override` is a dict of real CERPS parameter names (e.g. {"kappa_a": 0.5}).
    Returns a flat metrics dict; never fabricates values.
    """
    cfg = legacy_config()
    if override:
        cfg.update(override)
    L = cfg["L"]
    IX, IU, IP = _snapshot_init(n_init, L, seed=seed)
    r = run_cerps(seed=seed, T=float(T), n_init=n_init, init="snapshot",
                  init_x=IX, init_u=IU, init_pos=IP,
                  snapshot_times=snapshot_times,
                  record_after=record_after, rec_every=rec_every, **cfg)
    c = r["counts"]; info = r["info"]; s = r["summary"]; snaps = r["snaps"]
    out = dict(
        seed=int(seed), exit=float(info[0]), t=float(info[1]),
        N_final=float(info[2]), div=int(info[7]), death=int(info[6]),
        term_death=int(info[10]), enter=int(info[8]),
        cyc_frac=float(s[1]),
        mG=float(s[3]), mP=float(s[4]), mK=float(s[5]), mC=float(s[6]),
        sG=float(s[7]), sP=float(s[8]),
        uG=float(s[11]), uP=float(s[12]), uK=float(s[15 - 4 + 2]),
        suG=float(s[15]), suP=float(s[16]),
    )
    # attractor distribution for GRN candidates (needs snapshots)
    if len(snaps):
        xG = snaps[:, 1]; xP = snaps[:, 2]
        on = 1.0
        out["frac_G_on"] = float(np.mean(xG > on))
        out["frac_P_on"] = float(np.mean(xP > on))
        out["frac_both_on"] = float(np.mean((xG > on) & (xP > on)))
        out["frac_both_low"] = float(np.mean((xG < on) & (xP < on)))
        out["frac_G_only"] = float(np.mean((xG > on) & (xP < on)))
        out["frac_P_only"] = float(np.mean((xG < on) & (xP > on)))
    return out


class Simulator:
    """Conceptual wrapper around run_cerps anchored on legacy_config().

    Mirrors the paper-workflow docs' `Simulator` but is backed by the real
    kernel. `config` = real CERPS param overrides; `candidates` sweep axes are
    defined in run_candidates.py.
    """

    def __init__(self, config=None, T=4000.0, n_init=600,
                 seeds=(3, 11, 21), record_after=0.7):
        self.base = legacy_config()
        if config:
            self.base.update(config)
        self.T = T
        self.n_init = n_init
        self.seeds = seeds
        self.record_after = record_after

    def run(self, override=None, snapshot_times=(), seeds=None):
        seeds = seeds or self.seeds
        results = []
        for sd in seeds:
            cfg = dict(self.base)
            if override:
                cfg.update(override)
            L = cfg["L"]
            IX, IU, IP = _snapshot_init(self.n_init, L, seed=sd)
            r = run_cerps(seed=sd, T=float(self.T), n_init=self.n_init,
                          init="snapshot", init_x=IX, init_u=IU, init_pos=IP,
                          snapshot_times=snapshot_times,
                          record_after=self.record_after, rec_every=4, **cfg)
            c = r["counts"]; info = r["info"]; s = r["summary"]; snaps = r["snaps"]
            m = dict(seed=sd, exit=float(info[0]), t=float(info[1]),
                     N_final=float(info[2]), div=int(info[7]),
                     death=int(info[6]), term_death=int(info[10]),
                     cyc_frac=float(s[1]),
                     mG=float(s[3]), mP=float(s[4]), mK=float(s[5]), mC=float(s[6]),
                     sG=float(s[7]), sP=float(s[8]),
                     uG=float(s[11]), uP=float(s[12]),
                     suG=float(s[15]), suP=float(s[16]))
            if len(snaps):
                xG = snaps[:, 1]; xP = snaps[:, 2]; on = 1.0
                m["frac_G_on"] = float(np.mean(xG > on))
                m["frac_P_on"] = float(np.mean(xP > on))
                m["frac_both_on"] = float(np.mean((xG > on) & (xP > on)))
                m["frac_both_low"] = float(np.mean((xG < on) & (xP < on)))
                m["frac_G_only"] = float(np.mean((xG > on) & (xP < on)))
                m["frac_P_only"] = float(np.mean((xG < on) & (xP > on)))
            results.append(m)
        return results
