# Two CERPS Implementations — API Clarification (correction of an earlier error)

> **Status.** This document supersedes the "API correction (2026-10-19)" notices previously added to several files in this folder. Those notices were **wrong** and are retracted: they asserted that the class-based `cerps/` API does not exist, whereas it does exist in this repository. The authoritative account of the two implementations is given below.
>
> **Author.** Rongsheng Huang, School of Science, Jimei University, Xiamen 361021, China; hrs@jmu.edu.cn. NSFC 12401641.

---

## 1. What went wrong

An earlier pass inspected **only** the external companion package
(`CERPS_complete_package`, whose kernel is `cerps_model/src/cerps.py` exposing
`run_cerps(**params)` + `CERPS_DEFAULTS`) and then asserted — globally, and inside
*this* repository — that

> "`cerps.config`, `cerps.batch`, `CERPSConfig`, `GRNConfig`, `RatesConfig.delay_steps`,
> `MemoryConfig`, `legacy_config()`, … **None of those exist in the shipped code**."

That assertion is **false for this repository**. Those names exist, are exported,
and are covered by the test suite. The error arose from conflating two distinct
codebases. It has now been corrected everywhere it appeared.

---

## 2. The two implementations (authoritative)

### A. `cerps/` — class-based framework (THIS repository, v1.1.0)

The primary platform. See `cerps/__init__.py`, `docs/model_spec.md`.

| Item | Location | Notes |
|---|---|---|
| `GRNConfig`, `EpiConfig`, `DivisionConfig`, `MemoryConfig`, `RatesConfig`, `SpaceConfig`, `CERPSConfig` | `cerps/config.py` | dataclasses, lines 20 / 66 / 101 / 114 / 131 / 163 / 191 |
| **`legacy_config()`** | `cerps/config.py:214` | returns the 2-gene Table-1 configuration of Huang, Situ & Lei (2024) |
| **`RatesConfig.delay_steps`** | `cerps/config.py:156` | **P8 feedback delay — implemented** (default 0 = instantaneous) |
| `effective_niche_signal()` | `cerps/rates.py:159` | delays the niche (SC-count) signal by `delay_steps` |
| `Simulator._q_hist` | `cerps/simulator.py:81,178–180` | per-step history feeding the delayed signal |
| `cerps.batch` harness | `cerps/batch.py` | `run_one`, `sweep`, `aggregate`, `report`, `homeostasis`, `transitions`, `health_check`, `validate`, `load` |
| Tests | `tests/test_cerps.py`, `tests/test_audit_fixes.py` | comprehensive suite; see §3 |

**P6 (non-Markovian memory) and P8 (feedback delay) are both implemented here.**

### B. `cerps_model/` — numba kernel in the external companion package

A separate artifact (e.g. the local `CERPS_complete_package` checkout):
* `cerps_model/src/cerps.py` — single entry point
  `run_cerps(seed, T, n_init, init, init_x, init_u, init_pos, record_after, rec_every, snapshot_times, **params)`
  plus the 50-parameter `CERPS_DEFAULTS` dict; internally `@njit simulate_cerps`.
* `analysis/src/abc.py` — ABC-SMC against Paul 2015 / Buenrostro 2018.
* **Has no `delay_steps`** (no delay buffer in the kernel).

This kernel is what the `cerps_sbi/` adapter wraps, and it is what produced the
candidate sweeps reported here (`results/candidate_sweep_v1.json`, `v2.json`).

---

## 3. Consequences for the paper workflow

1. **All statements about `h_ref_S`, `kappa_a`, `s_part`, `b`, `epa`, `grid/L/motil/crowding`**
   in `candidate_register.md` §2/§7 refer to **implementation B** (its parameter names),
   because the executed sweeps ran there.
2. **The claim "P8 delay is NOT implemented" applied only to B.** In implementation A
   the primitive is available via `RatesConfig.delay_steps`. Therefore the C1
   **delayed-niche-feedback** sub-hypothesis is **testable now**, using `cerps/`
   rather than requiring a code extension. This converts a previously "impossible"
   item into a concrete next experiment.
3. `cerps_sbi/` remains useful (it is how the numbers were produced and it carries
   the ABC layer), but its "P8 not implemented" note must be read as scoped to
   kernel B only.
4. Nothing about the numerical results changes — they were produced honestly on the
   kernel they claim. Only the *provenance wording* and the *availability of P8*
   were wrong.

---

## 4. Recommended next step for C1 (P8)

Test delayed niche feedback in **implementation A**:

```python
from cerps import batch, legacy_config
from cerps.config import replace

base = legacy_config()                      # Table-1 anchor
for d in (0, 5, 20, 50):
    cfg = replace(base, rates=replace(base.rates, delay_steps=d))
    res = batch.run_one(cfg, seed=3, tag=f"delay_{d}")
```

Expected observable (to be measured, not assumed): a phase lag between the stem-cell
fraction and the niche signal, and a change in clone-size variance — the quantities
already identified as C1's falsifiable observables.
