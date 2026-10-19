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

> **API correction (2026-10-19).** Earlier drafts referenced `cerps.config` / `cerps.batch` / `cerps.simulator`, `CERPSConfig`, `GRNConfig`, `RatesConfig.delay_steps`, `MemoryConfig`, `legacy_config()`, etc. **None of those exist in the shipped code** (`cerps_model/src/cerps.py` exposes only `run_cerps(seed, T, n_init, init, init_x, init_u, init_pos, record_after, rec_every, snapshot_times, **params)` + `CERPS_DEFAULTS`). `execution-runbook.md`, `candidate_register.md` §2/§4, and `novelty_register.md` N1–N4 have been rewritten to the real API. A thin, honest adapter **`cerps_sbi/`** (repo root) re-exposes `legacy_config()` / `Simulator` / `run_one` on top of the real kernel and **flags the divergences** (notably: **P8 delay is NOT implemented** — `delay_steps` absent). The verification sweep output is committed at **`results/candidate_sweep_v1.json`**; the summary is `candidate_register.md` §7. Use `cerps_sbi` for all candidate sweeps.
>
> **Data-backed status (2026-10-19).** `candidate_register.md` and `novelty_register.md` are the first *real-evidence* deliverables: candidates C1–C4 are tied to specific GEO accessions / CARLIN / RNA velocity, and the adversarial gate found **no 🔴 scooping and no un-argued 🟠** — all claims are 🟡 (adjacent, controllable) or ⚪. A **verification sweep over the real CERPS kernel** is DONE (T=1500, N=200, seeds 2; raw: `CERPS_complete_package/results/candidate_sweep_v1.json`; summary: `candidate_register.md` §7). Key outcomes: **C1(P1) real** (~90% term-death effect), **C2(P6) real** (seed-dependent lineage polarization), **C3(P3) honest bulk-null** (needs clone metric), **C4 GRN** (both-on fraction = 0, supports negative finding). **P8 delay is NOT in the code** — C1's delay sub-hypothesis stays proposed. Production run (T=4000/N=600/seeds3×) + calibration against Paul2015/Buenrostro2018 remain before submission.

## Honest boundaries (must appear in the paper's Limitations)

- Semi-quantitative anchoring (anchor on `legacy_config()`); `u` is a single scalar per gene; continuous mode ≠ legacy (not a smooth substitute); no native parameter inversion (needs `cerps-sbi`); statistically reproducible, not bit-level reproducible.

---
*Chinese working versions of these documents live in the local WorkBuddy workspace; this folder is the English, repository-versioned copy.*
