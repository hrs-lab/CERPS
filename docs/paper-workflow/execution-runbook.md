# CERPS Public-Data Mechanism Mining — Concrete Workflow (Runbook)

> Companion: `research-proposal.md`. This file gives **directly executable** commands and scripts. All CERPS calls are checked against the **real** `CERPS` kernel in `cerps_model/src/cerps.py` (the only model file shipped in this repo).
>
> **API clarification (supersedes the retracted "API correction" notice).** This project contains **two** CERPS implementations. **(A)** `cerps/` in *this repository* — the class-based framework providing `CERPSConfig`, `GRNConfig`, `EpiConfig`, `DivisionConfig`, `MemoryConfig`, `RatesConfig` (including **`delay_steps`** = P8 feedback delay), `SpaceConfig`, `legacy_config()` and the `cerps.batch` harness; covered by `tests/` and reproducing paper Table 3. **(B)** the external companion package's numba kernel `cerps_model/src/cerps.py`, exposing only `run_cerps(**params)` + `CERPS_DEFAULTS` plus the ABC layer `analysis/src/abc.py`. An earlier claim that the names in (A) "do not exist" was **wrong** — it came from inspecting only (B) — and is hereby retracted. **Consequence: the P8 feedback delay IS available in (A) via `RatesConfig.delay_steps`.** Full account in **`api-clarification.md`**. The executed candidate sweeps reported below ran on (B) through `cerps_sbi/`.

---

## Directory structure (concrete, this repo)

```
CERPS_complete_package/                # the model + analysis package
├── cerps_model/src/cerps.py           # THE kernel (run_cerps, CERPS_DEFAULTS)
├── analysis/src/abc.py                # real ABC-SMC SBI layer (Paul2015 target)
├── analysis/src/mechanism_scan.py     # real M1–M4 mechanism scans
├── data/  results/                    # data + outputs
└── cerps_sbi/                         # NEW: honest adapter (this workflow)
    ├── __init__.py                    # legacy_config(), Simulator, run_one
    └── run_candidates.py              # C1–C4 sweeps (this workflow)
```

---

## 00 · Environment (tested on this machine)

```bash
# Use the managed venv that already has numpy + numba.
#   Win: C:/Users/黄荣生/.workbuddy/binaries/python/envs/default/Scripts/python.exe
#   (Python 3.13.14, numpy 2.4.6, numba 0.68.0 — verified working)
PY="C:/Users/黄荣生/.workbuddy/binaries/python/envs/default/Scripts/python.exe"

# The kernel is a single file; no pip install. Just point PYTHONPATH at it:
export CERPS_SRC="C:/Users/黄荣生/Downloads/CERPS_complete_package/cerps_model/src"
# To run the cerps_sbi adapter, also put the package parent on the path:
export CERPS_PKG="C:/Users/黄荣生/Downloads/CERPS_complete_package"

$PY -c "import sys; sys.path.insert(0,'$CERPS_SRC'); from cerps import run_cerps, CERPS_DEFAULTS; print('CERPS OK, params=', len(CERPS_DEFAULTS))"
```

> To build a fresh venv instead: `python -m venv venv && venv/Scripts/activate && pip install "numpy<2.5" numba`. (numba 0.68 needs numpy 2.x; the `default` venv above is already correct.)

---

## 00b · The `cerps_sbi` adapter (honest wrapper)

`cerps_sbi/__init__.py` re-exposes the methodology's conceptual names on top of the real kernel, and **records where the real code diverges**:

- `legacy_config()` → `dict(CERPS_DEFAULTS)` merged with the **locked manuscript spatial config**
  `BASE = {crowding:1.5, L:8.0, grid:16, n_niche:25, rho_S:2.5}` (the anchor used by `analysis/src/abc.py` and `analysis/src/mechanism_scan.py`). This IS the calibrated anchor; there is no other "legacy" object.
- `Simulator(config, T, n_init, seeds, record_after)` → wraps `run_cerps` and returns a list of metric dicts.
- `run_one(override, seed, T, n_init, snapshot_times)` → one realization, returns flat metrics.

**Honesty flags (do not silently ignore):**
1. **P8 (feedback delay) IS implemented** in the `cerps/` framework (`RatesConfig.delay_steps` + `rates.effective_niche_signal()`), but **is absent from the `cerps_model` kernel**. Any candidate relying on a delay axis (C1) can therefore be tested in `cerps/`; the sweeps in this runbook used the `cerps_model` kernel, so delay was not exercised there. See `api-clarification.md`.
2. **Memory (P6)** = the `epa` (epigenetic-age) accumulator, gated by `kappa_a` (accelerates death with divisional age) and `div_death_age`.
3. **Partitioning (P3)** = `s_part` (expression partitioning noise at mitosis).
4. **Epigenetics (P5)** = `u` + `k_maint` + `w_u=(w0,wG,wP,wK,wC)` + `D_u`.
5. **Niche (P1)** = the S-block: `S0, n_niche, rho_S, delta_S, D_S, uptake, KS, hS, h_maint_S, h_ref_S, well_mixed`.
6. **GRN (P4 / mutual repression)** = `alpha, b` (repression weights), `rho, k, s, n_hill, sigma, zeta, lam`.

---

## 01 · Fetch data (public datasets, pin accession before execution)

```bash
# 01_fetch_data.sh  —— GEO series matrix (expression + metadata)
mkdir -p data/raw && cd data/raw
wget -q https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81682/suppl/GSE81682_RAW.tar
wget -q https://ftp.ncbi.nlm.nih.gov/geo/series/GSE89nnn/GSE89433/suppl/GSE89433_RAW.tar
wget -q https://ftp.ncbi.nlm.nih.gov/geo/series/GSE53nnn/GSE53359/suppl/GSE53359_RAW.tar
wget -q https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96772/suppl/GSE96772_RAW.tar
# Unpack; organize per-dataset README into counts matrix + cell annotation.
```

> CARLIN barcodes (Bowling 2020) and RNA-velocity (La Manno 2018) are obtained from the original papers' supplements. CellPhoneDB CLI on scRNA gives ligand–receptor (P4 quorum / P7 spatial evidence).

---

## 02 · Preprocessing (scRNA/scATAC → CERPS init arrays)

```python
# 02_preprocess.py
import numpy as np, scanpy as sc

adata = sc.read_h5ad("data/raw/hspc.h5ad")
sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
genes = ["GATA1", "SPI1", "CEBPA", "GFI1"]          # anchor gene set (4 genes)
X = np.asarray(adata[:, genes].X, float)
Xn = X / np.quantile(X, 0.95, axis=0, keepdims=True) * 1.5   # scale to x≈0..3

labels = adata.obs["cell_type"].values
init_cells = []
for lab in ["HSC", "MPP", "GMP", "MEP"]:
    idx = np.where(labels == lab)[0][:50]
    init_cells.extend(idx.tolist())
init_x = Xn[init_cells]                              # (n, 4)
init_u = np.clip(init_x / 1.5, 0.05, 0.95)           # u initial (open correlate)
init_pos = np.random.default_rng(0).uniform(0, 8.0, (len(init_cells), 2))
np.save("data/init_x.npy", init_x)
np.save("data/init_u.npy", init_u)
np.save("data/init_pos.npy", init_pos)
print("init_x", init_x.shape)
```

> **scATAC → u prior.** Buenrostro 2018 / Lara-Astiaso 2014 accessibility per locus can tighten the `u` initialization and the `w_u` targets (P5). This is a semi-quantitative anchor, not a bit-level reverse-map (honesty boundary).

---

## 03 · Build config from data (real API)

```python
# 03_build_config.py
import sys, numpy as np
sys.path.insert(0, "C:/Users/黄荣生/Downloads/CERPS_complete_package")  # for cerps_sbi
from cerps_sbi import legacy_config

cfg = legacy_config()                 # CERPS_DEFAULTS + locked spatial BASE
# Optionally nudge toward data: e.g. shift epigenetic targets from scATAC
cfg["w_u"] = (0.8, 1.0, 1.0, 1.0, 1.0)     # w0, wG, wP, wK, wC
# Save for later stages
np.save("data/base_cfg.npy", cfg)    # (plain dict; reload via pickle/json)
import json
json.dump(cfg, open("data/base_cfg.json", "w"))
print("anchor params:", {k: cfg[k] for k in ("beta0","kappa0","k_maint","rho_S","n_niche")})
```

> **No `CERPSConfig` class exists.** The "config" is simply a `dict` of real
> CERPS parameter names consumed by `run_cerps(**params)`. `legacy_config()`
> returns the calibrated anchor; every candidate is an `override` dict on top.

---

## 04 · Calibration wrapper `cerps-sbi` (real SBI layer)

CERPS has **no native parameter inversion**. The real SBI layer ships as
`analysis/src/abc.py` (ABC-SMC against the **Paul 2015** mouse-myeloid scRNA
summary + cross-validation against **Buenrostro 2018** scATAC). Its prior ranges
(`ABC_PRIORS` in `cerps_sbi`) are:

```
r_eff 0.45–0.90 | r_max 0.85–0.99 | h_ref_S 0.05–0.80 |
k_maint 0.01–0.15 | s_div 0.01–0.08 | beta0 0.15–0.60 |
sigma 0.01–0.08 | D_u 1e-5–8e-5
```

```python
# 04_calibrate.py  (thin driver over analysis/src/abc.py)
import sys
sys.path.insert(0, "C:/Users/黄荣生/Downloads/CERPS_complete_package/analysis/src")
# abc.py reads /home/z/my-project/CERPS/... hardcoded; edit those two paths
# to the local CERPS_complete_package layout, then:
#   python analysis/src/abc.py
# Output: ABC posterior over the 8 coupling params; posterior-predictive check
# on Paul2015; cross-validation correlation on Buenrostro2018 accessibility.
```

> **Honest boundary.** The ABC target summaries are *population/gene-activity
> layer* (Gata1/Spi1 joint distribution along pseudotime, cycling fraction,
> erythroid:myeloid ratio, locus-accessibility ordering). They are **not** a
> bit-level reverse-map of any single cell. Conclusions stay at the mechanism /
> semi-quantitative level (see `paper-workflow-playbook.md` §honesty).

---

## 05 · Mechanism perturbation sweep (controlled comparison)

```python
# 05_mechanism_sweep.py
import sys, numpy as np, json
sys.path.insert(0, "C:/Users/黄荣生/Downloads/CERPS_complete_package")
from cerps_sbi import run_one

# ---- Example A: P1 niche feedback (C1, real axis) -------------------------
for tag, ov in [("off", {"h_ref_S": 0.0}),
                ("default", {"h_ref_S": 0.40}),
                ("strong", {"h_ref_S": 1.20})]:
    m = run_one(override=ov, seed=3, T=4000.0, n_init=600, snapshot_times=(3800,))
    print(tag, {k: m[k] for k in ("N_final","div","term_death","cyc_frac","mG","mP")})

# ---- Example B: P6 memory-of-age (C2) ------------------------------------
for tag, ov in [("off", {"kappa_a": 0.0}),
                ("default", {"kappa_a": 0.03}),
                ("strong", {"kappa_a": 0.15})]:
    m = run_one(override=ov, seed=3, T=4000.0, n_init=600, snapshot_times=(3800,))
    print(tag, {k: m[k] for k in ("N_final","term_death","cyc_frac")})

# ---- Example C: P3 partitioning noise (C3) -------------------------------
for tag, ov in [("off", {"s_part": 0.0}),
                ("default", {"s_part": 0.06}),
                ("strong", {"s_part": 0.15})]:
    m = run_one(override=ov, seed=3, T=4000.0, n_init=600, snapshot_times=(3800,))
    print(tag, {k: m[k] for k in ("sG","sP","mG","mP")})

# ---- Example D: GRN toggle cross-repression (C4 control) -----------------
for tag, ov in [("weak", {"bG":0.3,"bP":0.3}),
                ("default", {"bG":1.0,"bP":1.0}),
                ("strong", {"bG":3.0,"bP":3.0})]:
    m = run_one(override=ov, seed=3, T=4000.0, n_init=600, snapshot_times=(3800,))
    print(tag, {k: m[k] for k in ("frac_G_only","frac_P_only","frac_both_on")})
```

> **C1's P8 (delay) is UNIMPLEMENTED** — there is no `delay_steps`. To verify
> the delay hypothesis you must first add a delay buffer to `simulate_cerps`
> (store past S or past u and feed `delay_steps` back into the niche/epi target)
> and re-run. Until then C1 is reported as **P1-verified, P8-proposed**.

---

## 06 · Statistical validation

```python
# 06_validate.py
import numpy as np
# (a) multi-seed MC error: run each override over seeds=(3,11,21), report mean±sd;
#     a mechanism effect is real only if it exceeds the sd across seeds.
# (b) permutation / necessity: compare real GRN (b as inferred) vs shuffled-B
#     null on sc_fraction / cycling fraction; one-sided p = mean(real<=null).
# (c) Sobol first-order sensitivity over the 8 ABC priors (reuse analysis/src/abc.py
#     sampling + run_one as the simulator).
# (d) PPC: posterior theta -> simulate summary vs observed Paul2015; PIT coverage ~0.95.
# (e) hold-out: withhold one dataset's time point / perturbation, measure
#     2-Wasserstein(simulated state, observed) using scipy.stats.wasserstein_distance.
```

---

## 07 · Visualization (7 main figures)

```python
# 07_figures.py
import numpy as np, matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
# Fig 2d calibration fit: observed vs simulated summary (from 04).
# Fig 3 Waddington landscape: U = -log KDE(p(x,u)) from snaps[:,1:5] (x) and [5:9] (u).
# Fig 4 parameter-sweep heatmap: read N_final / cyc_frac from 05 jsonl.
# Fig 5 time series: counts[:,0] (t) vs counts[:,1] (N) / [3:7] (mean x).
# Fig 6 PPC coverage / edge bootstrap CI / nested-model Bayes factor.
```

---

## 08 · Orchestration

```bash
# 08_run_all.sh
PY="C:/Users/黄荣生/.workbuddy/binaries/python/envs/default/Scripts/python.exe"
export CERPS_SRC="C:/Users/黄荣生/Downloads/CERPS_complete_package/cerps_model/src"
export CERPS_PKG="C:/Users/黄荣生/Downloads/CERPS_complete_package"
$PY 02_preprocess.py
$PY 03_build_config.py
$PY analysis/src/abc.py                 # 04 calibration (edit hardcoded paths first)
$PY 05_mechanism_sweep.py
$PY 06_validate.py
$PY 07_figures.py
echo "DONE -> results/"
```

---

## Key execution points (against audit conclusions)

1. **Always anchor on `legacy_config()`** (CERPS_DEFAULTS + locked spatial BASE). Mechanism comparison is within-family; do not compare against an unrelated parameterization.
2. **`record_after` is a FRACTION of T** (e.g. 0.7 → record last 30%). Exclude transients before computing steady-state statistics.
3. **`health_check` first** (n/a as a function — instead run a short T=200 pilot and confirm `info[0]==0` (no MAXCELL exit) and N_final ≪ 12000). If `exit==3`, lower `beta0`/`rho_S` or raise `crowding`.
4. **Multi-seed report mean±sd**: conclusions must exceed the MC error across seeds (3,11,21).
5. **SBI is non-skippable**: CERPS has no native inversion; `analysis/src/abc.py` is the layer that gives "data → parameters". Forward `run_one` sweeps only do controlled comparison.
6. **Honest boundaries**: `u` is one scalar per gene; molecular claims need external validation; this pipeline outputs mechanism screening and semi-quantitative conclusions at the gene-activity layer. **P8 delay is not in the code** — flag any delay claim as proposed-until-coded.
