# CERPS Public-Data Mechanism Mining — Manuscript Draft Outline (IMRaD)

> **Purpose.** This file is the "writing skeleton" for `research-proposal.md`, giving a directly expandable IMRaD draft outline and figure-caption points. All numbers come from the completed CERPS audit and mechanism exploration (empirically verified, committed to `hrs-lab/CERPS`); entries marked `[to-be-verified]` need confirmation after fitting real data.
> **Author info (default):** Rongsheng Huang, School of Science, Jimei University, Xiamen 361021, China; hrs@jmu.edu.cn. Acknowledgment: National Natural Science Foundation of China (Grant No. 12401641).
> **LaTeX note:** this machine's MiKTeX lacks `authblk`; use the built-in author block (see memory convention), do not use `\affil`/`\email`.
>
> **⚠ API correction (2026-10-19).** Some module/class names used below (`CERPSConfig`, `GRNConfig`, `RatesConfig`, `MemoryConfig`, `EpiConfig`, `DivisionConfig`, `SpaceConfig`, `delay_steps`, `batch.sweep`) **do not exist in the shipped kernel** (`cerps_model/src/cerps.py` exposes only `run_cerps(**params)` + `CERPS_DEFAULTS`). The real, runnable API + `cerps_sbi` adapter are in `execution-runbook.md` §00b; `legacy_config()` there is the calibrated anchor. **P8 delay is NOT implemented** — mark any delay-based claim as proposed. use the built-in author block (see memory convention), do not use `\affil/\email`.

---

## Title candidates
1. *A controlled in-silico dissection of core mechanisms governing hematopoietic lineage commitment, anchored to single-cell atlases*
2. *From correlation to causation: a mechanistic screen of stem-cell fate decisions using a general differentiation model and public single-cell data*
3. *Waddington landscapes re-engineered: isolating niche, epigenetic, and feedback delays in cell-fate control*

---

## Abstract (skeleton, ~250 words)
Single-cell atlases have mapped the **correlational structure** of lineage differentiation, yet the **causal mechanisms** remain without controlled isolation. We propose an integration framework: using the general cell-differentiation model **CERPS** as a causal sandbox, calibrated with public scRNA/scATAC/lineage data and anchored on the `legacy` mode calibrated against published work, we systematically dissect the core mechanisms of HSC myeloid/lymphoid differentiation. Through controlled in-silico perturbations and simulation-based inference (SBI), we quantify the **independent contribution** of [GRN topology / signaling-pathway delay P8 / epigenetic memory P6] to fate decisions (H1–H4). Key findings: a purely mutual-repression network **cannot** robustly reproduce the "all-on = stemness" metaphor (single-on state dominant ~93%, all-on only ~0.5–1%, opposite in distribution direction to Duddu 2020) [to-be-verified: corresponding conclusion after real-data fitting]; feedback delay P8 and epigenetic memory P6 respectively explain [X%] / [Y%] of fate shift. This work shows that a forward mechanism model + public data can, under the **semi-quantitative-anchoring** paradigm, upgrade descriptive atlases into testable causal claims. [Limitations in Discussion]

---

## 1. Introduction
### 1.1 Background and gap
- Single-cell omics reveal the high-dimensional correlational structure of differentiation, but **observation ≠ mechanism** (correlation cannot separate driver from accompaniment).
- Existing computational methods emphasize "dimensionality reduction / clustering / trajectory inference", lacking **controlled causal-intervention** capability; mechanistic ODE/ABM models are hard to calibrate with real data.
- **CERPS** is positioned as a generalization of HSL2024 (J Theor Biol 577:111664), providing pluggable mechanism modules, but has not previously been systematically coupled with public data.

### 1.2 Falsifiable hypotheses (H1–H4, must pass the adversarial-review gate first)
- **H1** (GRN topology): if repression maintains stemness, CERPS should reproduce "all-on dominant"; **measurement refutes** → instead "repression yields single-on dominant + double-on tail".
- **H2** (signaling pathway / P8 delay): if niche feedback contains delay, lineage-bias phase shifts measurably.
- **H3** (epigenetic dynamics / P6 memory): if epigenetic inheritance contains non-Markovian memory, fate decisions show history dependence.
- **H4** (landscape structure): the reconstructed potential surface should show [multiple attractors + quantifiable barriers].

### 1.3 Research objective
Use CERPS to close the tension between "data phenomenon" and "mechanism explanation", producing statistically supported causal claims.

---

## 2. Methods
### 2.1 Data integration and calibration
- **Data sources:** scRNA-seq (GSE81682/GSE72857/GSE89433/GSE116256), scATAC-seq (GSE53359/GSE96772), spatial + CellPhoneDB, CARLIN lineage, RNA velocity.
- **Mapping:** expression → initial `x` and GRN prior `B/C`; ATAC → epigenetic initial `u` and `s_self`; lineage → fate-transition ground truth; spatial → paracrine source term.
- **Calibration layer `cerps-sbi`:** CERPS has no native inversion; use ABC / neural posterior to estimate `p(θ|s)`; anchor on `legacy_config()` (audit verified SC→SC deviation −0.006 pp, anchor max deviation from published Table 3 = 0.087 pp).

### 2.2 Controlled in-silico perturbations
- Forward activation + reverse ablation + counterfactual triad; `batch.sweep` multi-seed (≥ 8); `rec_tran_from=2000` excludes transients.
- `health_check` first, confirm CFL / `dt` discretization error ≪ mechanism effect.

### 2.3 Statistical validation
- bootstrap edge stability, permutation necessity, Sobol sensitivity, PPC, hold-out 2-Wasserstein.

### 2.4 Honest boundaries (methodological statement)
- semi-quantitative anchoring; `u` single scalar per gene (molecular-resolution gap); continuous mode ≠ legacy (no smooth substitute); non-bit-level reproducibility.

---

## 3. Results
### 3.1 Calibration fidelity (Fig 2)
CERPS anchor reproduces published steady state; hold-out 2-Wasserstein distance [to-be-verified].

### 3.2 Attractor structure of the mutual-repression network (Fig 3 + core finding)
- **Key result:** under 3-gene full mutual repression, single-on dominant ~93%, stable double-on tail ~7%, all-on only ~0.5–1% (from `docs/multigene_attractors.md` measurement).
- **Tension with literature:** Duddu 2020 RACIPE statistics are single-on ~60% / double-on ~36% / triple excluded — **opposite distribution direction**, showing a pure-repression network cannot robustly reproduce the "all-on = stemness" metaphor.
- `label==0` systematically overestimates the true all-on fraction due to continuous-mode label freezing (**must be corrected**).

### 3.3 Mechanism dissection (Fig 4)
P8 delay, P6 memory, GRN rewiring each effect size + Sobol main-effect ranking [to-be-verified: numbers after real fitting].

### 3.4 Signal / epigenetic time series (Fig 5)
After P8 enabled, oscillation phase shifts measurably; after P6 enabled, fate shows history dependence.

### 3.5 Statistical robustness (Fig 6)
All key claims pass ≥2 independent tests; effect size ≫ seed standard deviation.

---

## 4. Discussion
### 4.1 Mechanism interpretation
- Repression yields **diverse phenotypes but biases toward single-on**; stemness maintenance more likely depends on a **mutual-activation core** (Boyer 2005) + niche feedback, not pure repression.
- P8/P6 are gaps "weak in evidence but likely to yield new results", giving testable directions.

### 4.2 Comparison with existing models
- vs purely data-driven (trajectory inference): we provide **controlled causal intervention**;
- vs traditional ODE/ABM: we provide **calibration with public data + SBI posterior**.

### 4.3 Limitations (must write, maps to 2.4)
- semi-quantitative anchoring, molecular-resolution gap, continuous ≠ legacy, no native inversion, non-bit-level reproducibility.
- CERPS conclusions must be compared under same-family classification, not across modes quantitatively.

### 4.4 Outlook
- hand mechanism clues to high-resolution models for molecular validation;
- the `mechanisms.md` mechanism-variable table can serve as a community benchmark.

---

## 5. Figure Captions (points)
- **Fig 1:** study pipeline (data → calibration → perturbation → statistics → paper).
- **Fig 2:** data integration and calibration fit (left: scRNA UMAP anchored genes; right: CERPS vs published steady state).
- **Fig 3:** Waddington landscape (KDE potential) + 3-gene attractor pie (single-on / double-on / all-on).
- **Fig 4:** mechanism-dissection heatmap (rows = mechanism, cols = phenotype fraction) + Sobol bars.
- **Fig 5:** P8/P6 on/off time-series comparison (with phase annotation).
- **Fig 6:** PPC coverage band + permutation necessity + Wasserstein bars.
- **Fig 7 (optional):** spatial niche and paracrine self-organization snapshot.

---

## 6. Supplementary
- S1: data accessions and preprocessing scripts;
- S2: CERPS config schema and `mechanisms.md` variable table;
- S3: `cerps-sbi` calibration-layer description and approximation assumptions;
- S4: adversarial literature review `novelty_register.md`;
- S5: full `results.jsonl` checkpoints and reproduction environment (pinned `requirements`).

---

## 7. To-close list (must close before writing)
- [ ] H1–H4 numerical conclusions after real-data fitting (replace all [to-be-verified])
- [ ] adversarial review `novelty_register.md` run and archived
- [ ] journal selected (Cell Systems / Nat Commun / Cell Stem Cell / PLOS Comp Biol)
- [ ] author block set by LaTeX built-in (not authblk)
