# Research Proposal: Mining Core Mechanisms of Cell Differentiation by Integrating Public Biological Datasets on the CERPS Platform

> Companion documents: this proposal builds on prior work — `CERPS_audit_report` (platform audit and fixes), `docs/mechanisms.md` (mechanism-variable register), `docs/multigene_attractors.md` (multi-gene attractor evidence), `docs/validation_report.md` (reproducibility record). All CERPS module names, config class names, and batch APIs come from the already-committed-and-pushed `hrs-lab/CERPS@main` (`d63c4bb` / `12b4a5d` / `3cfb29c`).

> **API clarification (supersedes the retracted "API correction" notice).** This project contains **two** CERPS implementations. **(A)** `cerps/` in *this repository* — the class-based framework providing `CERPSConfig`, `GRNConfig`, `EpiConfig`, `DivisionConfig`, `MemoryConfig`, `RatesConfig` (including **`delay_steps`** = P8 feedback delay), `SpaceConfig`, `legacy_config()` and the `cerps.batch` harness; covered by `tests/` and reproducing paper Table 3. **(B)** the external companion package's numba kernel `cerps_model/src/cerps.py`, exposing only `run_cerps(**params)` + `CERPS_DEFAULTS` plus the ABC layer `analysis/src/abc.py`. An earlier claim that the names in (A) "do not exist" was **wrong** — it came from inspecting only (B) — and is hereby retracted. **Consequence: the P8 feedback delay IS available in (A) via `RatesConfig.delay_steps`.** Full account in **`api-clarification.md`**. The executed candidate sweeps reported below ran on (B) through `cerps_sbi/`.

---

## 0. Proposal summary (one-sentence positioning)

Using **hematopoietic stem cell (HSC) differentiation toward myeloid / lymphoid** as the demonstration system, we pour public scRNA-seq / scATAC-seq / spatial-transcriptomics / lineage-barcode datasets through a three-step pipeline — **calibration → simulation → validation** — into CERPS, systematically mining three classes of core mechanisms (GRN topology, key signaling pathways, epigenetic dynamics) via **controlled in-silico mechanism perturbations + simulation-based inference (SBI/ABC)**, and outputting publishable evidence in three forms: Waddington landscape, parameter sensitivity, and posterior predictive checks (PPC).

---

## 1. Scientific questions and core mechanism hypotheses

### 1.1 Open scientific questions
1. **GRN topology:** which combinations of mutual repression / mutual activation among genes are **necessary** (not merely sufficient) to maintain the "stemness + multi-lineage commitment" bistable structure?
2. **Signaling pathways:** among niche feedback (Q), contact inhibition, paracrine field, and **feedback delay (P8)**, which dominates HSC quiescence-activation oscillation and lineage bias?
3. **Epigenetic dynamics:** how do epigenetic locking rates (`w0/w1/e0/e1`) and division-inheritance fidelity (cpg / gaussian kernel) determine the **plasticity window** of fate decisions?
4. **Non-Markovian memory (P6):** does a "division-history dependence" beyond the current state exist, and can it be detected in lineage-barcode data?

### 1.2 Falsifiable core hypotheses (all switchable and testable in CERPS)
- **H1 (GRN):** in a G≥3 mutual-repression network, the "all-on = stemness" metaphor is **not robustly reproducible** (pre-verified by `docs/multigene_attractors.md`: single-on dominant ~93%, all-on only ~0.5–1%); a truly multipotent core needs mutual activation (Boyer 2005-style layering).
- **H2 (signaling):** introducing `RatesConfig.delay_steps>0` into the niche-feedback loop produces a measurable phase lag and quiescence oscillation, explaining clone-size fluctuations observed in CARLIN clonal data.
- **H3 (epigenetics):** lowering `EpiConfig` locking rate widens the plasticity window, raising the transitional-cell fraction and making fate decisions more sensitive to initial `u`.
- **H4 (memory):** after enabling `MemoryConfig(enabled=True)`, progeny fate correlation exceeds what is predicted from the current state alone, matching within-clone correlation in lineage barcodes.

---

## 2. Technical basis: CERPS module → mechanism mapping (audited, usable)

| Biological mechanism layer | CERPS module / config | Toggle switch (see `docs/mechanisms.md`) |
|---|---|---|
| Gene-expression dynamics / GRN | `cerps/grn.py` + `GRNConfig(B, C, regulation, n_hill, s_self, k_deg, …)` | `B/C` adjacency matrices, `regulation="product"\|"additive"` |
| Epigenetic dynamics | `cerps/epigenetics.py` + `EpiConfig(kernel, cpg_mode, w0,w1,e0,e1, …)` | kernel type, `u` dynamics rate, division-inheritance fidelity |
| Niche / contact-inhibition signal | `cerps/rates.py` + `RatesConfig(theta0, rho0, theta_rho, …)` | Q feedback, local-density inhibition |
| Paracrine signal field | `cerps/spatial.py` + `SpaceConfig(paracrine, D, …)` | secretion → diffusion-reaction → gene response |
| Feedback delay | `RatesConfig.delay_steps` + `rates.effective_niche_signal` | delay steps (default 0, backward-compatible) |
| Non-Markovian memory | `MemoryConfig(enabled, tau, driver, mem_couple)` | leaky integrator (default off) |
| Fate decision / population statistics | `cerps/simulator.py` + `batch.homeostasis/transitions` | classification threshold, fate-transition matrix |
| Controlled sweep / health diagnostics | `cerps/batch.py` (`sweep/aggregate/health_check/validate`) | Cartesian sweep, CFL convergence estimate |

---

## 3. Required data types and public dataset catalog

> General principle: CERPS is a **forward simulator with no native parameter inversion**; experimental data serve "calibrate initial state + validate output", and parameter inference is done by an external SBI wrapper (see §5.2).

### 3.1 Data layers and formats

| Layer | Data type | Format | How poured into CERPS |
|---|---|---|---|
| Expression state / initial condition | scRNA-seq (sorted HSPC and progeny) | `AnnData(.h5ad)` or 10x `matrix.mtx`+`features/barcodes`; UMI count matrix + cell-type annotation | normalized expression → `Simulator(init={"x":(n,G),"u":(n,G)})`; threshold `s_self` from quantiles |
| GRN-topology prior | scRNA-seq-inferred GRN | edge list `gene_i, gene_j, sign(+/-), weight` | write into `GRNConfig.B/C` matrices |
| Epigenetic accessibility / locking prior | scATAC-seq (or 10x Multiome) | `peak matrix(.h5ad)` + gene-level accessibility score | accessibility → `u` initial / `s_self` correction |
| Signal / niche | spatial transcriptomics (Visium/Slide-seq) + ligand-receptor DB (CellPhoneDB/NicheNet) | `spatial/filtered_feature_bc_matrix` + coordinates `tissue_positions.csv`; LR score matrix | spatial density → `SpaceConfig` contact inhibition; LR expression → paracrine source term |
| Fate / lineage | lineage barcode (CARLIN / CRISPR) or mitochondrial mutation; RNA velocity | `cell × barcode` matrix / VCF + scRNA; spliced/unspliced counts | lineage → `transitions()` ground truth; velocity → pseudotime-direction prior |
| Dynamical time series | live-cell imaging / resampled time series | cell-level time-series table | compare against `records` time-series trajectories |

### 3.2 Recommended public datasets (pin accession at execution)

- **scRNA-seq (anchor HSC differentiation)**
  - Mouse HSPC: `Nestorowa et al. 2016` (*Blood*), **GSE81682**; `Paul et al. 2015` (*Nature*), **GSE72857** (myeloid progenitors).
  - Human HSPC: `Velten et al. 2017` (*Mol Syst Biol*), **GSE89433**; human HSPC barcode `Pellin et al. 2019`, **GSE116256**. Also connect Human Cell Atlas hematopoietic reference.
- **scATAC-seq (epigenetic prior)**
  - Mouse HSC: `Lara-Astiaso et al. 2014` (*Immunity*), **GSE53359**; human HSC: `Buenrostro et al. 2018` (*Cell*), **GSE96772**; or 10x Multiome bone-marrow dataset.
- **Spatial / niche**
  - Bone-marrow spatial transcriptomics (Visium/Slide-seq, recent Nature-level studies); ligand-receptor inference via **CellPhoneDB v5** / **NicheNet** on the scRNA-seq above.
- **Lineage / dynamics**
  - `Bowling et al. 2020` (*Nature*) CARLIN mouse HSC clone tracking; clonal hematopoiesis single-cell `Fabre et al. 2022` (*Science*); RNA velocity `La Manno et al. 2018` (*Nature*).

> Note: accessions must be re-checked for the latest version before execution; this proposal uses hematopoiesis as the demonstration, but the method transfers to other lineages (intestinal epithelium, neurogenesis).

---

## 4. Overall analysis pipeline

```
Public datasets ──► [① preprocessing & pseudotime]
                      │  Scanpy/Monocle: normalize, cluster, UMAP, pseudotime, velocity
                      ▼
                  [② GRN inference] ── pySCENIC/SCENIC ──► B/C matrix (with bootstrap stability)
                  [③ epigenetic prior] ── scATAC ──► u initial / s_self correction
                  [④ signal/niche] ── spatial + LR DB ──► paracrine source / contact-inhibition params
                      │
                      ▼
            [⑤ SBI calibration wrapper]  (external layer, not CERPS-native)
             summary stats s(y_obs) ↔ s(y_sim(θ)), sbi/pyabc estimate p(θ|s)
                      │
                      ▼
            [⑥ mechanism mining: controlled in-silico perturbations]
             KO/ OE genes (B/C) · alter signal (Q/contact/paracrine/delay) · alter epigenetic rate · cut P6/P8
             anchored on legacy_config(), batch.sweep × multi-seed
                      │
            ┌─────────┴──────────┐
            ▼                    ▼
      [⑦ statistical validation]   [⑧ functional interpretation]
      bootstrap/permutation       GO/KEGG/MSigDB enrichment
      Sobol sensitivity           compare vs literature mechanisms
      posterior predictive check (PPC)  generate testable hypotheses
            │                    │
            └─────────┬──────────┘
                      ▼
                [⑨ visualization]  Waddington landscape / parameter-sweep heatmap / PPC coverage
                      ▼
                [⑩ paper integration]
```

---

## 5. Mechanism-mining strategy (three core mechanism classes)

### 5.1 GRN topology (necessary vs sufficient)
- Use pySCENIC on scRNA-seq to infer regulons; take high-confidence edges to build `GRNConfig.B/C`.
- **Necessity test:** in CERPS compare this GRN against an "edge-shuffled null model" (`batch.sweep` with multiple permutations), using a permutation test to judge "is the inferred topology necessary to reproduce the observed proportion distribution" (cf. the four-attractor tension conclusion in `docs/multigene_attractors.md`).
- **Combination effect:** only qualitative two-primitive sweeps (e.g., repression × activation), output directional conclusions; do not claim molecular-level mechanism.

### 5.2 Signaling pathways (incl. P8 delay) — needs SBI wrapper
- CERPS has no native inversion, so **add an external wrapper layer `cerps-sbi`** (lightweight, non-invasive):
  - define summary statistics `s`: cluster proportions, mean-expression trajectory along pseudotime, velocity-field divergence, pseudotime entropy;
  - use `batch.run_one(cfg)` as the simulator, estimate `p(θ|s)` with `sbi` (neural posterior) or `pyabc` (ABC-SMC);
  - **nested model comparison:** compare marginal likelihood / WAIC of nested variants ("no delay" vs "delay_steps>0", "no memory" vs "P6 on") to give Bayes factors.
- P8 validation: after enabling `delay_steps`, the quiescence-oscillation phase should shift measurably; PPC checks whether the simulated oscillation period covers CARLIN clone-size fluctuations.

### 5.3 Epigenetic dynamics
- From scATAC accessibility set `u` initial and `s_self`; sweep `EpiConfig` locking rate, quantify the "plasticity window" (monotonicity of transitional fraction vs locking rate).
- Kernel comparison (cpg `moment` vs gaussian): since already orthogonalized, only controlled comparison, no absolute-truth claim.

---

## 6. Statistical validation plan (mathematically rigorous)

1. **GRN-edge uncertainty:** bootstrap cells → re-infer edges → report edge stability (occurrence frequency); shuffle expression for null distribution → BH controls FDR.
2. **Necessity permutation test:** for "inferred topology vs shuffled topology", simulate with multiple seeds to get the proportion-distribution difference, permutation p-value judges topology necessity.
3. **Global sensitivity:** Sobol first/second-order indices or Morris screening, quantifying which mechanisms (P8 τ, epigenetic `w0/w1`, repression strength) dominate output variance (SC fraction, bistability margin) — directly answers "which mechanism is core".
4. **Convergence / artifact separation:** CFL number and `dt`-halving convergence estimate from `batch.health_check`, separating "mechanism effect" from "discretization artifact"; multi-seed report mean±sd (`batch.aggregate`).
5. **Posterior predictive check (PPC):** simulate from posterior θ, compare simulated vs observed summary-statistic distributions; use PIT (probability integral transform) or coverage to quantify calibration quality.
6. **Hold-out validation:** hold out a time point / perturbation condition, measure distance (e.g., 2-Wasserstein) between simulated and observed cell-state distributions as a prediction-accuracy metric.

---

## 7. Functional-interpretation plan

- **Hub-gene mapping:** map CERPS-identified key GRN hubs (high sensitivity / necessity) to known TFs (e.g., myeloid GATA1–PU.1 antagonism, C/EBPα), compare with literature mechanisms.
- **Pathway enrichment:** GO / KEGG / MSigDB Hallmark enrichment on gene sets related to `u` dynamics and fate switching.
- **Generate testable hypotheses:** e.g., "P8 delay in TGF-β niche feedback produces quiescence oscillation, explaining CARLIN clone-size fluctuations" → propose validation experiment (perturb niche factor, measure clonal output).
- **Cross-validate P6:** if P6 is selected, check lineage-barcode data for "progeny correlation beyond current-state prediction" as empirical evidence of non-Markovian memory.

---

## 8. Visualization plan (figure presentation points)

| Fig | Content | Presentation points |
|---|---|---|
| **Fig 1** | Study overview + pipeline schematic | end-to-end flow: three-layer data → calibration → simulation → validation |
| **Fig 2** | Data integration and calibration | (a) observed scRNA UMAP + cell types; (b) scATAC accessibility trajectory; (c) inferred GRN B/C heatmap (with bootstrap stability); (d) calibration fit (observed vs simulated summary stats) |
| **Fig 3** | Waddington landscape | reconstruct potential U(x) from CERPS simulation (Fokker–Planck / transition-rate estimate), overlay trajectories; juxtapose with data-derived pseudotime landscape |
| **Fig 4** | Mechanism dissection | parameter-sweep heatmap (e.g., repression strength × niche feedback → SC fraction / bistability margin); Sobol sensitivity bar chart |
| **Fig 5** | Signal and epigenetic dynamics | Q(niche)/u(epigenetic)/fate-probability time series; P8 delay phase shift, P6 memory effect |
| **Fig 6** | Statistical validation | hold-out prediction accuracy, edge bootstrap CI, PPC coverage; nested-model Bayes factors |
| **Fig 7** | Spatial niche (optional) | Visium spatial map overlaying paracrine-field prediction |
| **Supplementary** | full GRN, complete parameter table, simulation protocol, full sensitivity table | reproducibility material |

---

## 9. Paper framework (high-impact journal)

### 9.1 Journal targeting (choose by main contribution)
- **Method + discovery combined → `Cell Systems` / `Nature Communications`:** `cerps-sbi` pipeline as method novelty + HSC mechanism discovery as biological result.
- **Pure biological discovery → `Cell Stem Cell` / `Blood`:** if mechanism discovery is strong enough (e.g., P8 delay strongly supported by experimental data).
- **Pure computational method → `PLOS Computational Biology` / `eLife`:** if the SBI wrapper and validation framework are the core contribution.

### 9.2 Paper structure
1. **Title:** mechanism-clear discovery title (include system + core mechanism, e.g., "Feedback-delay-driven HSC quiescence oscillation").
2. **Abstract:** background, data, method (CERPS+SBI), core findings, significance.
3. **Introduction:** gap in differentiation-mechanism research; CERPS as a controlled mechanism-screening platform with its boundaries (cite the "semi-quantitative anchoring" positioning from `CERPS_audit_report`).
4. **Results**
   - R1 data integration and calibration (Fig 2)
   - R2 Waddington landscape reconstruction (Fig 3)
   - R3 GRN-topology necessity (H1, Fig 4a)
   - R4 signaling-pathway dominant factor (H2, P8 delay, Fig 5)
   - R5 epigenetic locking and plasticity window (H3)
   - R6 non-Markovian-memory empirical evidence (H4, P6)
   - R7 statistical validation and model comparison (Fig 6)
5. **Discussion:** biological meaning of mechanism conclusions; dialogue with existing literature (Van Liedekerke 2019 sharp-threshold tension, Duddu 2020 attractor statistics); platform boundaries (`u` single scalar per gene, continuous mode ≠ legacy); testable hypotheses and experimental proposals.
6. **Methods:** data preprocessing, GRN inference, SBI calibration, CERPS config, statistical tests, visualization; with code and data availability.
7. **Supplementary:** parameter table, full sensitivity table, simulation protocol.

### 9.3 Key success criteria
- At least one core mechanism (H1–H4) beats the null model on **hold-out data** with statistical significance;
- SBI-calibrated PPC coverage meets target (empirical coverage ≈ nominal level);
- Sensitivity analysis uniquely identifies the dominant mechanism variable;
- All figures reproducible (`batch` checkpoints + fixed seed + pinned `requirements.txt`).

---

## 10. Risks, mitigations, and limitations

| Risk | Mitigation |
|---|---|
| CERPS has no native inversion → calibration hard | build `cerps-sbi` wrapper (SBI/ABC), §5.2 |
| Resolution limit (`u` single scalar per gene) | molecular claims need external validation; CERPS only descends to "gene-activity layer" |
| continuous mode ≠ legacy boundary | always anchor on `legacy_config()`; mechanism comparison under same-family classification |
| public-data batch effects | cross-dataset harmonization (Harmony/Seurat v4); quantify only near anchor |
| numba non-bit-level reproducibility | statistical reproducibility suffices; record versions (see `validation_report.md`) |
| mechanism conclusion polluted by orthogonality | strictly follow `mechanisms.md` switch register and `batch.validate` trap checks |

---

## 11. Milestones (suggested timeline)

- **M1 (months 1–2):** data acquisition and preprocessing, GRN/epigenetic/signal prior extraction, reproduce `legacy_config()` anchor.
- **M2 (months 2–3):** implement `cerps-sbi` calibration wrapper, complete parameter posterior estimation.
- **M3 (months 3–5):** mechanism mining (H1–H4), `batch.sweep` controlled perturbations + sensitivity analysis.
- **M4 (months 5–6):** statistical validation (PPC, permutation, hold-out), functional interpretation, visualization.
- **M5 (months 6–8):** manuscript writing and submission; release code, data, simulation protocol.

---

*This proposal is a design document; it does not execute actual data download or simulation. At execution, pin dataset accessions per §3.2 and advance per the §11 milestones. All CERPS module names and APIs are empirically usable in `hrs-lab/CERPS@main`.*
