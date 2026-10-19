# CERPS + Public Datasets → High-Impact Paper: A Research Workflow Playbook

> **Scope.** This file is a **strategy-layer (research-process)** guide. It describes how to combine the *CERPS* general cell-differentiation model with real, public biological data to systematically produce publishable, high-impact papers. It does **not** contain shell commands or API-call details — those live in `execution-runbook.md` (the tactics layer); this document is the methodology that sits above it.
>
> **Core thesis.** CERPS is not a black box that "ingests data and emits conclusions"; it is a **controlled causal-inference sandbox**. The contribution of a high-impact paper is not "running a simulation", but "using the simulation to **cut causality, isolate mechanisms, and quantify uncertainty**" — which is CERPS's distinctive selling point over purely data-driven approaches (e.g., scRNA-seq alone yields only correlations).
>
> **⚠ API correction (2026-10-19).** Class/function names used in this document (`CERPSConfig`, `GRNConfig`, `RatesConfig.delay_steps`, `MemoryConfig`, `EpiConfig`, `DivisionConfig`, `SpaceConfig`, `cerps.batch`, `batch.sweep`) **do not exist in the shipped kernel** (`cerps_model/src/cerps.py` exposes only `run_cerps(**params)` + `CERPS_DEFAULTS`). The real, runnable API + `cerps_sbi` adapter are in `execution-runbook.md` §00b; `legacy_config()` there is the calibrated anchor. **P8 delay is NOT implemented.**

---

## 0. Overview: the paper-driving force model

```
Public data (correlational / descriptive)
        │  provides: initial states, GRN priors, phenotypic ground truth,
        │            spatial / lineage constraints
        ▼
CERPS (forward generative model / causal sandbox)
        │  provides: controlled in-silico perturbations, mechanism ablation,
        │            counterfactuals, uncertainty quantification
        ▼
Mechanism claims (falsifiable, statistically supported)
        │
        ▼
Highly-cited paper (Cell Systems / Nat Commun / Cell Stem Cell / PLOS Comp Biol)
```

**The source of a paper's novelty** = the **tension is closed** between "the phenomenon revealed by data" and "the mechanistic explanation supplied by CERPS's counterfactual experiments". Data alone cannot establish causality; model alone cannot establish reality — their coupling is the selling point.

---

## 1. The seven-stage research workflow

Each stage lists: **goal / input / CERPS action / role of public data / key deliverable / the DoD (Definition of Done) gate to the next stage.**

### Stage 1 — Anchor the scientific question and generate falsifiable hypotheses
- **Goal:** converge from "phenomenological curiosity" to "falsifiable mechanistic hypotheses".
- **Input:** domain intuition + preliminary literature.
- **CERPS action:** consult `docs/mechanisms.md` (the six-axis mechanism-variable table) and enumerate candidate mechanisms as **switchable variables** (e.g., P8 `delay_steps`, epigenetics `w1/e1`, P6 `mem_couple`).
- **Role of public data:** use *published conclusions* from public scRNA/scATAC to pin down the "phenomenon to be explained" (e.g., "a certain lineage transition is accompanied by a specific GRN rewiring").
- **Key deliverable:** 3–4 hypotheses H1–Hn, each of the form "**if** mechanism M holds, **then** CERPS should exhibit D (consistent / inconsistent with data)".
- **DoD:** every hypothesis is falsifiable (there exists a CERPS parameter set such that the prediction conflicts with data); an **adversarial literature review** is completed (excluding known results such as Bokes 2012, Jia 2018/2024, Wu–Grima–Jia 2025, to confirm novelty).

### Stage 2 — Inventory data assets and perform integrated calibration
- **Goal:** convert heterogeneous public data into CERPS **calibration anchors** and **validation benchmarks**.
- **Input:** accessions from GEO / Single Cell Portal, etc. (scRNA-seq, scATAC, spatial, lineage barcodes, RNA velocity).
- **CERPS action:** via an external `cerps-sbi` wrapper layer (**CERPS has no native inversion, so this layer is mandatory**) perform simulation-based inference (ABC / neural posterior) to estimate the posterior `p(θ|s)`; use `legacy_config()` as the **calibration zero-point**.
- **Role of public data:**
  - scRNA → initial expression `x` and GRN-topology prior (`B/C`);
  - scATAC → epigenetic initial `u` and threshold `s_self`;
  - lineage barcodes → fate-transition ground truth (validates `transitions()`);
  - spatial / CellPhoneDB → paracrine source terms and contact-inhibition parameters.
- **Key deliverable:** one calibrated `CERPSConfig` + a data→config mapping table + a hold-out partition.
- **DoD:** anchor self-check passes (legacy simulation deviates from published Table 3 by < 0.1 pp); the hold-out set does **not** participate in calibration.

### Stage 3 — Mechanism dissection: controlled in-silico perturbations
- **Goal:** perform **causal cutting** on each hypothesis — isolate "the standalone contribution of mechanism M".
- **Input:** the calibrated config from Stage 2.
- **CERPS action:**
  - **forward activation:** turn on mechanism M (`batch.sweep` over its parameter range, multiple seeds);
  - **reverse ablation:** set M to a null model (e.g., shuffled GRN, P8 `delay_steps=0`, P6 `enabled=False`) as control;
  - **counterfactual:** fix all other parameters, vary only M, observe whether phenotype / landscape changes abruptly.
- **Role of public data:** supply the "real target statistics to be reproduced" (e.g., a certain lineage fraction, transition probability) as the alignment benchmark.
- **Key deliverable:** per-mechanism "effect size + uncertainty" (mean±sd across seeds); necessity/sufficiency matrix.
- **DoD:** conclusions exceed MC error (effect size >> seed standard deviation); `batch.health_check` confirms discretization error is far smaller than the mechanism effect (isolating "numerical artifacts").

### Stage 4 — Statistical validation and functional interpretation
- **Goal:** upgrade "looks correlated" to "statistically supported causal claim".
- **CERPS action:**
  - bootstrap edge stability (are GRN edges robust?);
  - permutation necessity test (does fit degrade significantly after removing M?);
  - Sobol global sensitivity (which parameter drives the conclusion?);
  - posterior predictive check (PPC) (do simulated summaries cover real summaries?);
  - hold-out 2-Wasserstein distance (gap between simulated vs real distribution).
- **Role of public data:** supply empirical null distributions and permutation baselines.
- **Key deliverable:** statistical-validation table (p-values / effect sizes / CIs) + functional interpretation (what mechanism M "means" biologically).
- **DoD:** all key claims pass ≥ 2 independent statistical tests; effect direction is self-consistent with biological interpretation.

### Stage 5 — Visualization and story-line construction
- **Goal:** translate mechanistic conclusions into a causal narrative a reviewer reads at a glance.
- **CERPS action:** export `records` (time series), `snapshots` (high-dimensional states), `transitions()` (fate matrix), `health_check` diagnostics.
- **Key deliverable:** the 7-figure framework (see §3) — Waddington landscape reconstruction, mechanism-dissection heatmap, signal/epigenetic time series, statistical-validation panel.
- **DoD:** every figure answers one clear scientific question, no "decorative" figures; captions include methodological and uncertainty notes.

### Stage 6 — Manuscript writing and journal targeting
- **Goal:** package the work into a manuscript matching the target journal's taste.
- **CERPS action:** none (document production).
- **Key deliverable:** IMRaD-structured manuscript + supplementary material (with reproduction scripts and configs); explicit positioning of "method contribution vs discovery contribution".
- **Three journal options** (by contribution emphasis):
  - method + discovery → *Cell Systems* / *Nature Communications*;
  - pure discovery → *Cell Stem Cell* / *Blood*;
  - pure method → *PLOS Computational Biology*.
- **DoD:** novelty statement is backed by the adversarial review; limitations section is written honestly (see §2).

### Stage 7 — Peer-review response and reproducibility
- **Goal:** use CERPS's "controlled re-run" capability to answer reviewers quickly.
- **CERPS action:** run targeted ablation / sensitivity experiments against reviewer comments (CERPS's advantage peaks here — response experiments can be generated within hours).
- **Key deliverable:** response-experiment bundle + updated reproducible environment (`requirements` pinned, `results.jsonl` checkpoints).
- **DoD:** every reviewer concern has a corresponding simulation evidence; code/data availability statement is complete.

---

## 2. Boundaries that MUST be written honestly into the paper (to avoid a fatal reviewer hit)

These follow from the CERPS audit conclusions and are part of **methodological integrity**; they must appear explicitly in the "Limitations and Risks" section:

1. **Semi-quantitative anchoring positioning.** CERPS is a "mechanism-screening + anchor-based semi-quantitative" platform, not a molecular-resolution or quantitative-prediction engine. Quantitative claims must be benchmarked against the `legacy_config()` anchor; away from the anchor only direction and order of magnitude are guaranteed.
2. **`u` is a single scalar per gene.** Epigenetic state is coarse-grained to one activation ratio per gene; it cannot descend to the molecular mechanism of "which CpG site, which modification, which TF combination". Any molecular-level claim requires external high-resolution validation.
3. **Continuous mode ≠ legacy, cannot serve as smooth substitute.** Their classification-boundary topologies differ (`classify_legacy` is axis-aligned box, `continuous` is a bent hypersurface, `resolve_K` coincides only on the diagonal). Mechanism comparisons must be done under **same-family classification**; it is forbidden to treat continuous mode as a smooth version of legacy and compare quantitatively.
4. **No native parameter inversion.** Calibration depends on the external `cerps-sbi` layer; this wrapper and its approximation assumptions must be declared in the Methods section.
5. **Reproducibility tier.** Statistically reproducible, not bit-level reproducible (numba version causes RNG-stream differences) — the generating `numba/numpy` versions must be recorded.

> **Red line.** Do not claim "CERPS predicted a molecular mechanism" or "continuous-mode results are smooth approximations of legacy". This would invite reviewers to counter with conclusions such as Van Liedekerke 2019.

---

## 3. Paper figure framework (maps to Stage 5)

| Fig | Content | Data source (CERPS / public) | Scientific question answered |
|---|---|---|---|
| **Fig 1** | overall study-pipeline schematic | fusion of both sources | how we did it |
| **Fig 2** | data integration and calibration fit | public data + CERPS anchor | is the model faithfully calibrated to reality |
| **Fig 3** | Waddington landscape reconstruction (potential surface) | CERPS `snapshots` KDE | differentiation barriers and attractor structure |
| **Fig 4** | mechanism-dissection heatmap + Sobol sensitivity | `batch.sweep` + Stage 4 | which mechanism drives the phenomenon |
| **Fig 5** | signal pathway (P8 delay) / epigenetic (P6 memory) time series | `records` | dynamical causal chain |
| **Fig 6** | statistical-validation panel (PPC / permutation / Wasserstein) | Stage 4 | is the conclusion statistically robust |
| **Fig 7** (optional) | spatial niche and paracrine self-organization | `space` module | mechanism in the spatial dimension |

---

## 4. Interface with the tactics layer (Runbook)

| Strategy stage | Tactics segment invoked |
|---|---|
| Stage 2 calibration | Runbook `03_build_config` + `04_calibrate` |
| Stage 3 perturbation | Runbook `05_mechanism_sweep` + `06_validate` (health_check) |
| Stage 4 statistics | Runbook `06_validate` |
| Stage 5 visualization | Runbook `07_figures` |

**Discipline.** Before entering any strategy stage, first confirm its DoD is met; advancing to the next stage while any DoD is unmet means stacking a paper on unvalidated assumptions — the leading cause of rejection at high-impact journals.

---

## 5. One-line summary

> **Feed phenomena with public data, cut mechanisms with CERPS, prove robustness with statistics, and prevent derailment with honest boundaries — this four-step closed loop is the production line of a high-impact paper.**
