# cerps_sbi — honest adapter over the real CERPS kernel

This package wraps the **actual** CERPS kernel in `cerps_model/src/cerps.py` so
the `paper-workflow` methodology (docs/paper-workflow/*.md) is **executable**
rather than aspirational. It re-exposes the conceptual API names used in the
methodology docs (`legacy_config()`, `Simulator`, `run_one`) on top of the real
function `run_cerps(**params)`.

## The real API (what actually ships)

```python
from cerps import run_cerps, CERPS_DEFAULTS
run_cerps(seed, T, n_init, init, init_x, init_u, init_pos,
          record_after, rec_every, snapshot_times, **params)
```

There is **no** `CERPSConfig`, `GRNConfig`, `RatesConfig`, `MemoryConfig`,
`EpiConfig`, `DivisionConfig`, `SpaceConfig`, or `delay_steps` in **that
kernel** (`cerps_model/src/cerps.py`). The "config" there is a plain `dict` of
real CERPS parameter names (see `CERPS_DEFAULTS`, 50 keys).

> **Scope note (2026-10-22).** This refers to the `CERPS_complete_package`
> kernel only. The `cerps/` package in **this** repository *is* dataclass-based
> and does define `CERPSConfig`, `GRNConfig`, `EpiConfig`, `DivisionConfig`,
> `SpaceConfig` and P8's `delay_steps` — they are two different codebases, and
> the adapter below does not dispatch to this repository's `cerps/`.

## How to run

```bash
# managed venv already has numpy 2.4.6 + numba 0.68.0
PY="C:/Users/黄荣生/.workbuddy/binaries/python/envs/default/Scripts/python.exe"

# tell the adapter where the kernel lives (or rely on the auto-fallback)
export CERPS_SRC="C:/path/to/CERPS_complete_package/cerps_model/src"

cd C:/Users/黄荣生/CERPS
$PY cerps_sbi/run_candidates.py        # runs C1–C4 sweeps, writes ../results/candidate_sweep_v1.json
```

`legacy_config()` returns the **calibrated anchor** = `CERPS_DEFAULTS` merged with
the locked manuscript spatial config
`BASE = {crowding:1.5, L:8.0, grid:16, n_niche:25, rho_S:2.5}`
(the same anchor used by `analysis/src/abc.py` and `analysis/src/mechanism_scan.py`).

## Honesty flags (do NOT silently ignore)

1. **P8 (feedback delay) is NOT implemented.** There is no `delay_steps` in the
   kernel. C1's delay sub-hypothesis is a *proposed code extension* — flagged,
   never simulated.
2. **Memory (P6)** = the `epa` (epigenetic-age) accumulator, gated by `kappa_a`
   (accelerates death with divisional age) and `div_death_age`.
3. **Partitioning (P3)** = `s_part` (expression partitioning noise at mitosis).
4. **Epigenetics (P5)** = `u` + `k_maint` + `w_u=(w0,wG,wP,wK,wC)` + `D_u`.
5. **Niche (P1)** = the S-block: `S0, n_niche, rho_S, delta_S, D_S, uptake, KS,
   hS, h_maint_S, h_ref_S, well_mixed`.
6. **GRN (P4)** = `alpha, b` (repression weights, a **tuple** `b=(bG,bP,bK,bC)`),
   `rho, k, s, n_hill, sigma, zeta, lam`.

## SBI / calibration

CERPS has **no native parameter inversion**. The real SBI layer is
`analysis/src/abc.py` (ABC-SMC against Paul 2015 mouse-myeloid scRNA;
cross-validation on Buenrostro 2018 scATAC). Run it from the
`CERPS_complete_package` checkout (edit its two hardcoded `/home/z/...` paths
first). `cerps_sbi` only does forward, controlled mechanism sweeps.

## Files

- `__init__.py` — `legacy_config()`, `Simulator`, `run_one`, robust kernel-path
  resolution.
- `run_candidates.py` — C1–C4 sweeps mapped to the REAL parameter axes. Output:
  `../results/candidate_sweep_v1.json`.
