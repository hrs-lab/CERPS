# Paper Workflow — CERPS + Public Data → High-Impact Paper

This folder collects the methodology documents for turning the **CERPS** general
cell-differentiation model and public biological datasets into a publishable,
high-impact paper. The documents are split into a **strategy layer** (research
process, causal logic, paper architecture) and a **tactics layer** (executable
shell/Python, API calls).

## Documents

| File | Layer | Content |
|---|---|---|
| `paper-workflow-playbook.md` | Strategy (7 stages) | The research workflow: how to combine CERPS + public data to produce a paper; DoD gates per stage; honest-boundary section. |
| `research-proposal.md` | Proposal | Design of the HSC-differentiation mechanism-mining study: scientific questions, H1–H4, data catalog, pipeline, statistics, milestones. |
| `execution-runbook.md` | Tactics (shell) | Directly executable commands/scripts (`00_env` … `08_run_all`) mapped to the **real** CERPS kernel `cerps_model/src/cerps.py` (`run_cerps(**params)`). |
| `adversarial-literature-review-protocol.md` | Stage-1 DoD gate | Search queries + exclusion checklist (personal red lines + CERPS/HSC domain priors) + novelty decision matrix. |
| `candidate-mechanism-selection-protocol.md` | Mechanism nomination | Closed primitive vocabulary P1–P8, data→primitive mapping, two-primitive composition, triple scoring, adversarial gate, AI-assistant role. |
| `paper-draft-outline.md` | IMRaD skeleton | Title candidates, Abstract, H1–H4, Methods, Results, Discussion, 7-figure captions, Supplementary, to-close list. |
| `candidate_register.md` | §8 populated | **Real-data-backed candidate list** (C1–C4) with GSE81682/GSE89433/GSE53359/GSE96772/CARLIN/RNA-velocity evidence, falsifiable H sentences, triple scoring, priority ranking. |
| `novelty_register.md` | Stage-1 gate output | **Adversarial-review results** (real searches, 2026-10-19): claim inventory, personal red lines (Bokes 2012 / Jia 2018 / Wu–Grima–Jia 2025–26), domain priors, verdicts ⚪🟡🟠🔴, follow-up full-text reads. |
| `calibration-report.md` | Calibration (Runbook §04) | ABC-SMC against Paul 2015 / Buenrostro 2018: **non-convergent** (accept=0, 10-particle edge-collapsed posterior, ε floor ~4.17); relative model comparison still supports niche + gating. Honest identifiability limitation. |

## Core thesis

CERPS is a **controlled causal-inference sandbox**, not a black box. The paper's
novelty comes from **closing the tension** between "the phenomenon revealed by
data" and "the mechanistic explanation given by CERPS counterfactual
experiments". Data alone gives only correlation; model alone gives only
hypothesis — their coupling is the selling point.

## How to use

1. Start with `paper-workflow-playbook.md` for the seven-stage logic.
2. Before any novelty claim, run `adversarial-literature-review-protocol.md`, then record results in `novelty_register.md`.
3. To select which mechanisms to test, follow `candidate-mechanism-selection-protocol.md`; the populated output lives in `candidate_register.md`.
4. For the HSC demonstration design, see `research-proposal.md`; for the writing skeleton, `paper-draft-outline.md`.
5. For execution, follow `execution-runbook.md`.

> **API clarification (supersedes the retracted "API correction" notice).** This project contains **two** CERPS implementations. **(A)** `cerps/` in *this repository* — the class-based framework providing `CERPSConfig`, `GRNConfig`, `EpiConfig`, `DivisionConfig`, `MemoryConfig`, `RatesConfig` (including **`delay_steps`** = P8 feedback delay), `SpaceConfig`, `legacy_config()` and the `cerps.batch` harness; covered by `tests/` and reproducing paper Table 3. **(B)** the external companion package's numba kernel `cerps_model/src/cerps.py`, exposing only `run_cerps(**params)` + `CERPS_DEFAULTS` plus the ABC layer `analysis/src/abc.py`. An earlier claim that the names in (A) "do not exist" was **wrong** — it came from inspecting only (B) — and is hereby retracted. **Consequence: the P8 feedback delay IS available in (A) via `RatesConfig.delay_steps`.** Full account in **`api-clarification.md`**. The executed candidate sweeps reported below ran on (B) through `cerps_sbi/`.

## Honest boundaries (must appear in the paper's Limitations)

- Semi-quantitative anchoring (anchor on `legacy_config()`); `u` is a single scalar per gene; continuous mode ≠ legacy (not a smooth substitute); no native parameter inversion (needs `cerps-sbi`); statistically reproducible, not bit-level reproducible.

---
*Chinese working versions of these documents live in the local WorkBuddy workspace; this folder is the English, repository-versioned copy.*
