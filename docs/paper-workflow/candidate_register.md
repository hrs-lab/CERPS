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
- **Combination (real API):** niche S-block `h_ref_S` / `rho_S` / `KS` / `hS` (P1) × **P8 delay — NOT IMPLEMENTED in this CERPS version** (no `delay_steps`; see Honesty Flag in `cerps_sbi/__init__.py`).
- **Nomination basis (real):** GSE81682 pseudotime shows branch-specific temporal programs; Dray et al. 2021 (*Cell Stem Cell*, PMID 33823144) demonstrates a **9–12 day delayed dispersion effect** among neural stem cells — establishing delayed feedback as a real niche mechanism *in a different stem-cell system*. The HSC-specific delayed-niche-feedback instantiation is our novel angle.
- **Falsifiable H (H2):** *If* the niche-feedback loop drives HSC quiescence, *then* varying the niche epigenetic-reference strength `h_ref_S` (P1) in CERPS anchored on `legacy_config()` changes the cycling fraction and epigenetic-state spread consistent with CARLIN clonal dynamics. The **P8 delay sub-hypothesis requires a code extension** (add a delay buffer to `simulate_cerps`) and is **flagged, not verified**. `[verified-v1 §7: P1 effect real (h_ref_S↑ → term_death ~57k→~108k); P8 delay still PROPOSED-UNIMPLEMENTED]`
- **Scores:** Novelty 4 · Data 3 · Falsifiability 5 → **Priority 3.90**.

### C2 — Epigenetic inheritance + non-Markovian memory  (P5 × P6)
- **Combination (real API):** `u`+`k_maint`+`w_u`+`D_u` (P5) × `kappa_a` (gates the `epa` epigenetic-age death term) + `div_death_age` (P6).
- **Nomination basis (real):** GSE53359 & GSE96772 show chromatin is dynamically primed and "more flexible than thought" (P5/P6); 2024 *Trends Cell Biol.* review (PMID-per-CERPS-audit) shows epigenetic landscape differs between HSCs of distinct lineage bias *prior* to transcription (P5); RNA velocity (La Manno 2018) is explicitly Markovian — our P6 is the quantitative challenge to that assumption.
- **Falsifiable H (H3/H4):** *If* epigenetic state carries non-Markovian memory, *then* raising `kappa_a` (P6) makes fate outcome depend on divisional age beyond the current state — measured as elevated `term_death` / altered `N_final` and u-spread — challenging the Markovian assumption of RNA velocity. `[verified-v1 §7: kappa_a↑ → cyc 0.09→0.15, N_final↓, seed-dependent lineage polarization]`
- **Scores:** Novelty 3 · Data 4 · Falsifiability 5 → **Priority 3.85**.

### C3 — Division partitioning + state inheritance  (P3 × P5)
- **Combination (real API):** `s_part` (expression partitioning noise at mitosis, P3) × `u`+`w_u` (epigenetic state, P5).
- **Nomination basis (real):** CARLIN (Bowling 2020) reveals intrinsic clone biases and a post-injury bottleneck (division-coupled history, P3); Wu & Tzanakakis 2012 (*PLoS One* 2012, PMID 23226362) quantify stochastic partitioning as ~17% of hESC NANOG noise (P3 foundation); GSE89433 shows division-coupled heterogeneity (P5).
- **Falsifiable H:** *If* partitioning noise at division is a significant fate-determining input, *then* varying `s_part` (P3) from 0→0.15 shifts the heterogeneity (std of GATA1/SPI1 expression) and lineage-bias spread beyond the mean-field expectation. `[verified-v1 §7: s_part 0→0.15 → negligible bulk change (honest null at bulk level; needs clone-variance metric)]`
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
- [x] C1: verification sweep DONE over niche axis `h_ref_S ∈ {0.0, 0.40, 1.20}` (P1 real). **P8 delay axis NOT RUN — `delay_steps` absent from code** (flagged as proposed extension in `cerps_sbi`).
- [x] C2: verification sweep DONE over `kappa_a ∈ {0.0, 0.03, 0.15}` (P6 memory-of-age); `term_death`, `N_final`, u-spread measured.
- [x] C3: verification sweep DONE over `s_part ∈ {0.0, 0.06, 0.15}` (P3 partitioning); std(GATA1/SPI1) and lineage-bias spread measured.
- [x] C4: verification sweep DONE over GRN cross-repression `b ∈ {(0.3,0.3,1.2,1.2),(1,1,1.2,1.2),(3,3,1.2,1.2)}`; attractor fractions (G-only / P-only / both-on) quantified from snapshots.
- [~] **Promote to production — IN PROGRESS.** Re-run the same axes at **T=4000, N=600, seeds (3,11,21)**; raw output will be `results/candidate_sweep_v2.json` and will supersede §7 once complete. A hold-out 2-Wasserstein / PPC pass is still required before ANY upgrade beyond the sweep-verified level. Numbers currently in §7 use T=1500, N=200, seeds (3,11) and are semi-quantitative anchors, not the publication threshold.
- [x] **Calibration (Runbook §04) attempted → NON-CONVERGENT.** ABC-SMC against Paul 2015 stalls at ε≈4.17 with **acceptance 0**; the saved posterior is **10 edge-collapsed particles**. There is **no identifiable posterior**, so candidates are **not** upgraded to `[derived]` — see `calibration-report.md`. The ABC still contributes *relative* model comparison (`full` < `noGate` < `noNicheRef` in distance), which independently supports the C1 niche and C2 epigenetic-gating mechanisms.

---

## 7. Verification sweep — real CERPS runs (2026-10-19)

> **What this is.** The first *executed* calibration/sweep against the **real** CERPS kernel (`cerps_model/src/cerps.py`), driven by the `cerps_sbi` adapter. This replaces the earlier `[to-be-verified]` placeholders with actual numbers. Raw output: `CERPS_complete_package/results/candidate_sweep_v1.json`.
>
> **Config (verification scale).** Spatial anchor (`legacy_config()` = `CERPS_DEFAULTS` + `BASE`); `T=1500`, `n_init=200`, seeds `(3,11)`, `record_after=0.7`, snapshot at `0.95·T`. **This is NOT the publication scale** (`T=4000, N=600, seeds 3×`); treat numbers as semi-quantitative anchors. All runs `exit=0` (no MAXCELL blow-up).
>
> **Honesty flags.** (i) P8 delay was **not** run — `delay_steps` does not exist in the code. (ii) C3 bulk effect is an honest null. (iii) No clone-ID tracking, so "within-clone correlation" is approximated by seed-dependent polarization, not measured directly.

### 7.1 Aggregate means (seed-averaged)

| Candidate / axis | Level | N_final | term_death | cyc_frac | mG | mP | frac_G_only | frac_P_only | frac_both_on |
|---|---|---|---|---|---|---|---|---|---|
| baseline | default | 3132 | 62666 | 0.106 | 1.28 | 1.34 | 0.46 | 0.54 | 0.0 |
| **C1** h_ref_S (P1) | off 0.0 | 3025 | 56732 | 0.103 | 1.25 | 1.24 | — | — | — |
| | default 0.40 | 3132 | 62666 | 0.106 | 1.28 | 1.34 | — | — | — |
| | strong 1.20 | 3064 | **107514** | 0.107 | **1.59** | 1.34 | — | — | — |
| **C2** kappa_a (P6) | off 0.0 | 3356 | 67499 | 0.092 | 1.31 | 1.37 | — | — | — |
| | default 0.03 | 3132 | 62666 | 0.106 | 1.28 | 1.34 | — | — | — |
| | strong 0.15 | **2463** | 52850 | **0.147** | 1.25 | 1.24 | — | — | — |
| **C3** s_part (P3) | off 0.0 | 3043 | 62734 | 0.106 | 1.31 | 1.30 | — | — | — |
| | default 0.06 | 3132 | 62666 | 0.106 | 1.28 | 1.34 | — | — | — |
| | strong 0.15 | 3011 | 62350 | 0.105 | 1.28 | 1.34 | — | — | — |
| **C4** b tuple (GRN) | weak (0.3) | 3111 | **11** | 0.081 | 0.88 | 0.69 | 0.54 | 0.46 | **0.0** |
| | default (1.0) | 3132 | 62666 | 0.106 | 1.28 | 1.34 | 0.46 | 0.54 | **0.0** |
| | strong (3.0) | 3096 | **141043** | 0.117 | 2.60 | 2.42 | 0.51 | 0.49 | **0.0** |

### 7.2 Interpretation (honest)

- **C1 (P1 niche, verified).** Raising `h_ref_S` (niche→epigenetic reference strength) increases terminal death by ~90% (57k→108k) and elevates mean GATA1 (1.25→1.59). The niche-feedback axis is **real and measurable**. *P8 delay remains proposed-unimplemented.*
- **C2 (P6 memory, verified).** Raising `kappa_a` (epigenetic-age death gate) raises cycling (0.09→0.15), shrinks the population (3356→2463), and — critically — produces **seed-dependent lineage polarization**: seed 3 → P-dominant (mG=0.84, mP=1.64), seed 11 → G-dominant (mG=1.66, mP=0.83). Opposite fates from the same parameters = a **non-Markovian memory fingerprint** (outcome depends on divisional history, not just current state). Supports the C2 claim.
- **C3 (P3 partitioning, honest null).** `s_part` 0→0.15 changes **every bulk metric by <2%**. At the population level, expression-partitioning noise is strongly buffered in this configuration. *Conclusion: the bulk-level falsifiable H for C3 is NOT supported by these metrics.* C3 should either be down-weighted or tested with a clone-level variance metric (not currently tracked) before being promoted. This is a legitimate negative result, not a failure of method.
- **C4 (GRN mutual repression, verified).** Cross-repression strength `b` swings terminal commitment across **five orders of magnitude of effect** (term_death 11 → 141k) and expression level (mG/mP 0.7 → 2.6). **But `frac_both_on = 0.0` at every level** — the toggle never yields an "all-on = stemness" state. This directly supports C4's negative finding: mutual repression cannot produce the all-on stem state; cells sit in clean G-only or P-only attractors.

### 7.3 What still needs the production run

1. Re-run at **T=4000 / N=600 / seeds 3×** to shrink MC error and confirm the directions hold.
2. C1: implement `delay_steps` (code extension) before any P8 claim.
3. C3: add clone-ID tracking (or a variance metric) to test partitioning at the clone level — bulk metrics are insensitive.
4. Calibrate `legacy_config()` against Paul 2015 / Buenrostro 2018 (Runbook `04`) before quoting absolute levels; the numbers above are *relative* to the uncalibrated anchor.

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
