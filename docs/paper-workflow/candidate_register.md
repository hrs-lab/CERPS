# CERPS Candidate Mechanism Register — HSC Differentiation (real-data populated)

> **Companion.** Implements §8 of `candidate-mechanism-selection-protocol.md` with **real, citable public-dataset evidence** (this is the first data-backed candidate list; prior §8 was a placeholder). Read together with `novelty_register.md` (adversarial gate) and `adversarial-literature-review-protocol.md`.
>
> **Discipline reminder.** Mechanisms = combinations of primitives **P1–P8** (closed vocabulary, no free association). Single primitive = ablation control only; two-primitive interaction = main productive zone; ≥3 primitives require a dimensionality-reduction plan first. Every candidate below is falsifiable inside CERPS via a switchable parameter.
>
> **Status tags.** `[lit]` = established by cited literature; `[to-be-verified]` = must be confirmed by CERPS simulation after calibration.
>
> **Author.** Rongsheng Huang, School of Science, Jimei University, Xiamen 361021, China; hrs@jmu.edu.cn. NSFC 12401641.

---

## 1. Data → primitive nomination map (real datasets)

Each row is a *real* dataset whose published findings nominate specific CERPS primitives. Nomination is allowed **only** when a concrete observed phenomenon is cited (protocol §2).

| Dataset (accession) | Publication | Real observed phenomenon | Nominated primitive(s) | Evidence type |
|---|---|---|---|---|
| **GSE81682** | Nestorowa et al., *Blood* 2016 (PMID 27365425) | scRNA-seq of 1,656 mouse HSPCs; pseudotime trajectories for lymphoid / erythroid–megakaryocytic / granulocyte–macrophage branches; GM & EM branches share upregulation of **cell-cycle & oxidative-phosphorylation programs** | **P2** (cycle / temporal driver), **P5** (state), supports continuum | `[lit]` published abstract + GEO summary |
| **GSE89433** | Velten et al., *Nat Cell Biol* 2017 (PMID, doi:10.1038/ncb3493) | scRNA-seq of human HSPCs; **CLOUD-HSPCs** — lineage commitment is a *continuous* process, unilineage-restricted cells emerge from a continuum; **combinatorial gene-expression modules** | **P2**, **P5**, **P3** (division-coupled heterogeneity) | `[lit]` published abstract |
| **GSE53359** | Lara-Astiaso et al., *Science* 2014 (doi:10.1126/science.1256271) | scATAC-seq of mouse HSCs; chromatin-state dynamics; 48,415 enhancer regions, 17,035 lineage-specific *de novo*; enhancers open in undifferentiated states and close with maturation; authors state the process is **"more dynamic and flexible than previously thought"** | **P5** (epigenetic state), **P6** (memory / plasticity window) | `[lit]` published abstract |
| **GSE96772** | Buenrostro et al., *Cell* 2018 (PMID 29706549) | scATAC-seq of human hematopoiesis; **continuous regulatory landscape**; lineage-bias variation in multipotent states; GMP heterogeneity; explicitly frames results on a **Waddington landscape** | **P5**, **P6**, supports Fig 3 landscape framing | `[lit]` published abstract |
| **CARLIN** | Bowling et al., *Cell* 2020 (PMID 32413320) | CRISPR-barcoding of HSC clones; **intrinsic biases in fetal-liver HSC clone activity**; **clonal bottleneck** after injury; HSC clones preferentially localize along long bone (expansion potential **niche-influenced**) | **P3** (partitioning / inheritance), **P5** (state inheritance), **P1** (niche feedback) | `[lit]` published abstract |
| **RNA velocity** | La Manno et al., *Nature* 2018 (doi:10.1038/nature25915) | spliced/unspliced ODE model predicts future cell state; **explicitly assumes a Markovian process** ("cell fate depends only on the current state") | **P6** (non-Markovian memory directly challenges this assumption) | `[lit]` published method + reviews |

> **Active-primitive subset for HSC** = {P1, P2, P3, P5, P6}. P4 (quorum) and P7 (spatial coupling) are *not* directly nominated by the core datasets above and require spatial/CellPhoneDB evidence before promotion (see §5).

---

## 2. Candidate pool (§8 populated with real evidence)

Triple scoring: **Priority = 0.4·Novelty + 0.35·Data + 0.25·Falsifiability** (each 1–5). Novelty reflects residual gap after the adversarial gate in `novelty_register.md`.

### C1 — Niche feedback with delayed response  (P1 × P8)
- **Combination (real API):** niche S-block `h_ref_S` / `rho_S` / `KS` / `hS` (P1) × **P8 delay — implemented in the `cerps/` framework via `RatesConfig.delay_steps`**, though **absent from the `cerps_model` kernel** on which the sweeps below ran (see `api-clarification.md`).
- **Nomination basis (real):** GSE81682 pseudotime shows branch-specific temporal programs; Dray et al. 2021 (*Cell Stem Cell*, PMID 33823144) demonstrates a **9–12 day delayed dispersion effect** among neural stem cells — establishing delayed feedback as a real niche mechanism *in a different stem-cell system*. The HSC-specific delayed-niche-feedback instantiation is our novel angle.
- **Falsifiable H (H2):** *If* the niche-feedback loop drives HSC quiescence, *then* varying the niche epigenetic-reference strength `h_ref_S` (P1) in CERPS anchored on `legacy_config()` changes the cycling fraction and epigenetic-state spread consistent with CARLIN clonal dynamics. The **P8 delay sub-hypothesis is directly testable in `cerps/`** via `RatesConfig.delay_steps`; it has **not** been exercised here because these sweeps ran on the `cerps_model` kernel — it is therefore **unverified, not impossible**. `[verified-prod §7: P1 effect real and large — term_death +145% (124k→305k), N_final −26%; P8 delay available in `cerps/` but NOT yet exercised in these sweeps]`
- **Scores:** Novelty 4 · Data 3 · Falsifiability 5 → **Priority 3.90**.

### C2 — Epigenetic inheritance + non-Markovian memory  (P5 × P6)
- **Combination (real API):** `u`+`k_maint`+`w_u`+`D_u` (P5) × `kappa_a` (gates the `epa` epigenetic-age death term) + `div_death_age` (P6).
- **Nomination basis (real):** GSE53359 & GSE96772 show chromatin is dynamically primed and "more flexible than thought" (P5/P6); 2024 *Trends Cell Biol.* review (PMID-per-CERPS-audit) shows epigenetic landscape differs between HSCs of distinct lineage bias *prior* to transcription (P5); RNA velocity (La Manno 2018) is explicitly Markovian — our P6 is the quantitative challenge to that assumption.
- **Falsifiable H (H3/H4):** *If* epigenetic state carries non-Markovian memory, *then* raising `kappa_a` (P6) makes fate outcome depend on divisional age beyond the current state — measured as elevated `term_death` / altered `N_final` and u-spread — challenging the Markovian assumption of RNA velocity. `[verified-prod §7: kappa_a↑ → N_final −83%, cyc +78%. CORRECTED from v1 — replicates CONVERGE, not polarize (cross-seed sd of mG 0.922→0.149, mP 0.944→0.106); P6 reads as a variance-suppressing divisional-age filter]`
- **Scores:** Novelty 3 · Data 4 · Falsifiability 5 → **Priority 3.85**.

### C3 — Division partitioning + state inheritance  (P3 × P5)
- **Combination (real API):** `s_part` (expression partitioning noise at mitosis, P3) × `u`+`w_u` (epigenetic state, P5).
- **Nomination basis (real):** CARLIN (Bowling 2020) reveals intrinsic clone biases and a post-injury bottleneck (division-coupled history, P3); Wu & Tzanakakis 2012 (*PLoS One* 2012, PMID 23226362) quantify stochastic partitioning as ~17% of hESC NANOG noise (P3 foundation); GSE89433 shows division-coupled heterogeneity (P5).
- **Falsifiable H:** *If* partitioning noise at division is a significant fate-determining input, *then* varying `s_part` (P3) from 0→0.15 shifts the heterogeneity (std of GATA1/SPI1 expression) and lineage-bias spread beyond the mean-field expectation. `[verified-prod §7: s_part 0→0.15 shifts all bulk metrics by only ~1–4% (honest null at bulk level; needs clone-variance metric)]`
- **Scores:** Novelty 3 · Data 3 · Falsifiability 5 → **Priority 3.50**.

### C4 — Mutual-repression network (control / negative finding, not a mechanism candidate)
- **Combination (real API):** 4-gene toggle with cross-repression weights `bG`,`bP` (GATA1↔SPI1 mutual repression), varied from weak→strong (P4 / mutual-repression control).
- **Nomination basis (real):** GSE81682 & GSE89433 scRNA show diverse HSPC states; prior audit (`docs/multigene_attractors.md`) measured **single-on dominant ~93%, double-on tail ~7%, all-on only ~0.5–1%** under full mutual repression — *opposite* in distribution direction to Duddu 2020 (single-on ~60% / double-on ~36% / triple excluded).
- **Role:** **Refutes** the "all-on = stemness" metaphor (Boyer 2005-style mutual-activation core is instead implicated). Becomes a *negative/control finding* in Results, not a positive mechanism candidate.
- **Scores:** Novelty 4 · Data 5 · Falsifiability 5 → **Priority 4.60** (but excluded from the positive-candidate pool by design).

---

## 3. Priority ranking and gate status

| Rank | Candidate | Priority | Adversarial gate (`novelty_register.md`) | Decision |
|---|---|---|---|---|
| 1 | C4 (negative finding) | 4.60 | 🟡 controllable (quantitative refutation) | Results control / Discussion |
| 2 | C1 (P1×P8) | 3.90 | 🟡 adjacent (cite Dray 2021, draw boundary) | **Into Results** |
| 3 | C2 (P5×P6) | 3.85 | 🟡 adjacent (mechanism known, ours = quantitative) | **Into Results** |
| 4 | C3 (P3×P5) | 3.50 | 🟡 adjacent (coupling to HSC bias is the angle) | **Into Results** (provisional) |

No 🔴 (substantive scooping) or un-argued 🟠 entries were found → all four may enter the manuscript pending the simulation verification in §4.

---

## 4. Falsifiability / execution checklist (maps to Runbook)

- [x] Calibrate on `legacy_config()` anchor (CERPS_DEFAULTS + locked spatial BASE) using GSE81682/GSE89433 scRNA initial states + CARLIN clone fractions (Runbook `03_build_config` + `04_calibrate` via `analysis/src/abc.py`).
- [x] C1: verification sweep DONE over niche axis `h_ref_S ∈ {0.0, 0.40, 1.20}` (P1 real). **P8 delay axis NOT RUN** — available in `cerps/` via `RatesConfig.delay_steps`, absent from the `cerps_model` kernel used here.
- [x] C2: verification sweep DONE over `kappa_a ∈ {0.0, 0.03, 0.15}` (P6 memory-of-age); `term_death`, `N_final`, u-spread measured.
- [x] C3: verification sweep DONE over `s_part ∈ {0.0, 0.06, 0.15}` (P3 partitioning); std(GATA1/SPI1) and lineage-bias spread measured.
- [x] C4: verification sweep DONE over GRN cross-repression `b ∈ {(0.3,0.3,1.2,1.2),(1,1,1.2,1.2),(3,3,1.2,1.2)}`; attractor fractions (G-only / P-only / both-on) quantified from snapshots.
- [~] **Promote to production — IN PROGRESS.** Re-run the same axes at **T=4000, N=600, seeds (3,11,21)**; raw output will be `results/candidate_sweep_v2.json` and will supersede §7 once complete. A hold-out 2-Wasserstein / PPC pass is still required before ANY upgrade beyond the sweep-verified level. Numbers currently in §7 use T=1500, N=200, seeds (3,11) and are semi-quantitative anchors, not the publication threshold.
- [x] **Calibration (Runbook §04) attempted → NON-CONVERGENT.** ABC-SMC against Paul 2015 stalls at ε≈4.17 with **acceptance 0**; the saved posterior is **10 edge-collapsed particles**. There is **no identifiable posterior**, so candidates are **not** upgraded to `[derived]` — see `calibration-report.md`. The ABC still contributes *relative* model comparison (`full` < `noGate` < `noNicheRef` in distance), which independently supports the C1 niche and C2 epigenetic-gating mechanisms.

---

## 7. Production-scale verification sweep — real CERPS runs (PROD)

> **What this is.** The same C1–C4 axes re-run at **publication scale** against the **real** CERPS kernel (`cerps_model/src/cerps.py`) via the `cerps_sbi` adapter. It supersedes the earlier verification-scale pass and is the evidence base for the `[verified-prod]` tags in §2 and in `novelty_register.md`. Raw output: **`results/candidate_sweep_v2.json`** (39 runs; the earlier pass is retained as `candidate_sweep_v1.json`).
>
> **Config (PROD).** Spatial anchor (`legacy_config()` = `CERPS_DEFAULTS` + `BASE`); **`T=4000`, `n_init=600`, seeds `(3,11,21)`**, `record_after=0.7`, attractor snapshot at `0.95·T` (t = 3800). Wall clock 2637.6 s. **All 39 runs returned `exit=0`** — no MAXCELL blow-up, no extinction.
>
> **Robustness rule.** Three replicates let us separate signal from seed noise: an effect is called **robust** only when the level-to-level change is large compared with the cross-seed sd. `mG`/`mP` are strongly seed-dependent at most levels, so they are reported **with their sd** and are not over-interpreted.
>
> **Honesty flags.** (i) P8 delay not run — `delay_steps` does not exist in the kernel. (ii) No clone-ID tracking. (iii) Calibration is **non-convergent** (`calibration-report.md`), so every effect is **relative to the legacy anchor**, never an absolute calibrated level.
>
> **⚠ Corrections to the verification-scale pass.** Two v1 claims do **not** survive at PROD scale and are rewritten below: **C2 does not produce seed polarization** (replicates *converge*), and **C4's `frac_both_on` is not 0** (~1.8 % mean).

### 7.1 Aggregate results (mean ± sd over seeds 3 / 11 / 21)

| Candidate / axis | Level | N_final | term_death | cyc_frac | mG | mP | G_only | P_only | both_on |
|---|---|---|---|---|---|---|---|---|---|
| baseline | default | 1973 ± 329 | 143771 ± 9107 | 0.170 ± 0.016 | 1.44 ± 0.45 | 1.03 ± 0.43 | — | — | — |
| **C1** h_ref_S (P1) | off 0.0 | 2052 ± 236 | 124483 ± 10366 | 0.162 ± 0.013 | 0.91 ± 0.96 | 1.33 ± 0.94 | — | — | — |
| | default 0.40 | 1973 ± 329 | 143771 ± 9107 | 0.170 ± 0.016 | 1.44 ± 0.45 | 1.03 ± 0.43 | — | — | — |
| | strong 1.20 | **1524 ± 362** | **305490 ± 7038** | 0.175 ± 0.007 | 0.72 ± 0.06 | **2.33 ± 0.04** | — | — | — |
| **C2** kappa_a (P6) | off 0.0 | 2390 ± 132 | 158474 ± 8478 | 0.144 ± 0.011 | 1.59 ± 0.92 | 0.90 ± 0.94 | — | — | — |
| | default 0.03 | 1973 ± 329 | 143771 ± 9107 | 0.170 ± 0.016 | 1.44 ± 0.45 | 1.03 ± 0.43 | — | — | — |
| | strong 0.15 | **403 ± 196** | **98838 ± 8454** | **0.256 ± 0.014** | 1.63 ± 0.15 | 0.86 ± 0.11 | — | — | — |
| **C3** s_part (P3) | off 0.0 | 2115 ± 178 | 144343 ± 7927 | 0.167 ± 0.014 | 1.46 ± 0.94 | 0.99 ± 0.89 | — | — | — |
| | default 0.06 | 1973 ± 329 | 143771 ± 9107 | 0.170 ± 0.016 | 1.44 ± 0.45 | 1.03 ± 0.43 | — | — | — |
| | strong 0.15 | 2021 ± 42 | 142581 ± 6152 | 0.169 ± 0.002 | 1.45 ± 0.83 | 1.02 ± 0.84 | — | — | — |
| **C4** b tuple (GRN) | weak (0.3) | 2094 ± 103 | **14 ± 5** | 0.129 ± 0.003 | 0.77 ± 0.39 | 0.77 ± 0.44 | 0.43 | 0.55 | 0.017 ± 0.029 |
| | default (1.0) | 1973 ± 329 | 143771 ± 9107 | 0.170 ± 0.016 | 1.44 ± 0.45 | 1.03 ± 0.43 | 0.64 | 0.34 | 0.017 ± 0.017 |
| | strong (3.0) | 1909 ± 186 | **420251 ± 9185** | 0.198 ± 0.004 | 1.71 ± 1.21 | 2.97 ± 1.19 | 0.31 | 0.67 | 0.019 ± 0.019 |

### 7.2 Interpretation (honest, including the two corrections)

- **C1 (P1 niche) — VERIFIED, large effect.** off → strong gives **term_death +145 %** (124 483 → 305 490) and **N_final −26 %** (2052 → 1524); both changes dwarf their sd (10 366 and 7038), so the direction is robust. Notably, at `h_ref_S = 1.20` the *replicate scatter collapses* (mG 0.72 ± 0.06, mP 2.33 ± 0.04, versus ±0.96 / ±0.94 at off) — a strong niche reference drives every replicate to the same SPI1-dominant outcome. The *direction* of the mG/mP shift is non-monotone across levels and must not be quoted as a simple monotone trend. *P8 delay remains proposed-unimplemented.*

- **C2 (P6 memory) — VERIFIED, but the v1 claim is REFUTED and replaced.** `kappa_a` produces the largest effects in the sweep: **N_final −83 %** (2390 → 403), **cycling +78 %** (0.144 → 0.256), term_death −38 %. However — correcting v1 — the replicates do **not** polarize; they **converge**. The cross-seed sd of mG falls **0.922 → 0.451 → 0.149** and of mP **0.944 → 0.431 → 0.106** as `kappa_a` increases. The correct reading is that the epigenetic-age gate acts as a **strong divisional-age selective filter that suppresses stochastic divergence between replicate histories** (it removes old-age cells so aggressively that essentially one composition survives). History-dependence is real — the outcome depends on accumulated divisional age, which is the non-Markovian ingredient — but the specific falsifiable observation written in v1 («opposite fates from identical parameters») is **not** reproduced and has been rewritten here and in `novelty_register.md` N2.

- **C3 (P3 partitioning) — honest null, CONFIRMED at production scale.** `s_part` 0 → 0.15 changes every bulk metric by ≈1–4 % (term_death −1.2 %, cyc +1.3 %, N_final −4.4 %), i.e. within the cross-seed sd and one to two orders of magnitude below the C1/C2 effects. The bulk-level falsifiable H for C3 is **not supported**. C3 requires a clone-level variance metric (not currently tracked) before it can be promoted; this remains a legitimate negative result, not a failure of method.

- **C4 (GRN mutual repression) — negative finding upheld, with a CORRECTED number.** Cross-repression strength `b` swings terminal commitment enormously (term_death 14 ± 5 → 420 251 ± 9185). **Correction to v1:** `frac_both_on` is **not 0**; it averages **0.017 / 0.017 / 0.019** across weak/default/strong (individual runs 0.000–0.050, with seed 21 giving exactly 0 at all three levels). An "all-on" state therefore occurs at the **~2 % level** — still far from a functional "all-on = stemness" attractor, so **the negative finding stands**, but the v1 sentence «frac_both_on = 0.0 at every level» must never be quoted. Its rarity is the point: mutual repression does not produce a substantial all-on population.

### 7.3 What remains open

1. **P8 delay (C1):** available **now** in `cerps/` (`RatesConfig.delay_steps`, `rates.effective_niche_signal()`); run the C1 delay sweep there before making any delayed-feedback claim. It was absent only from the `cerps_model` kernel used for the sweeps above.
2. **C3:** add clone-ID / within-clone variance tracking to test partitioning at clone level; bulk metrics are now demonstrably insensitive (≈1–4 % versus 145 % for C1).
3. **Re-centre C1:** the (degenerate) ABC posterior prefers `h_ref_S ≈ 0.64–0.80`, while the sweep used {0.0, 0.40, 1.20} — add a level near 0.75 to bracket it (`calibration-report.md` §5).
4. **Calibration:** non-convergent, so **no `[derived]` tag may be claimed** (`calibration-report.md`). Hold-out 2-Wasserstein / PPC gating is not satisfiable until a usable posterior exists.
5. **Replicates:** 5–10 seeds would tighten the already-clear C1/C2 separations and would decide the non-monotone mG/mP direction in C1.

---

## 5. Out-of-scope primitives (require additional data before promotion)

- **P4 (quorum sensing)** and **P7 (spatial coupling):** not nominated by the core six datasets. Promote only after adding spatial-transcriptomics (Visium/Slide-seq) + CellPhoneDB/NicheNet ligand–receptor inference; otherwise keep as Discussion outlook, not Results.
- **P2 (cycle driver):** nominated via GSE81682 cell-cycle programs but modeled only as an external `a_tf` periodic term; treat as a covariate, not a standalone candidate.

---

## 6. Honest boundary (carried into the paper's Limitations)

- CERPS `u` is **one scalar per gene** — epigenetic claims stay at the "gene-activity layer"; molecular CpG/TF mechanisms need external validation.
- **Continuous mode ≠ legacy** — all mechanism comparisons use the same-family classification anchored on `legacy_config()`.
- CERPS is **semi-quantitative / mechanism-screening**, not a molecular-resolution predictor.
- **No native parameter inversion** — calibration depends on the external `cerps-sbi` wrapper.
- **No identifiable calibration posterior** — ABC-SMC against Paul 2015 does **not** converge (acceptance 0, 10-particle edge-collapsed posterior; see `calibration-report.md`). All candidate effects must therefore be reported **relative to the legacy anchor**, never as absolute calibrated levels.
- **Non-bit-level reproducibility** — record `numba`/`numpy` versions.
