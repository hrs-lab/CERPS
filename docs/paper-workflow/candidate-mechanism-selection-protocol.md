# CERPS Candidate Biological-Mechanism Selection Protocol

> **Purpose.** Within the CERPS + public-data → paper workflow, this protocol standardizes "how to start from public data and, in a controlled way, select candidate biological mechanisms worth testing with CERPS". It is embedded between Stage 1 and Stage 3 of `paper-workflow-playbook.md` as the "mechanism-nomination sub-process".
> **Core discipline** (from your long-standing mechanism-exploration convention): a mechanism = a combination of primitives P1–P8; **single primitive = low value, two-primitive interaction = the main productive zone, three-or-more primitives require a dimensionality-reduction plan first**; **free association is forbidden**.
> **Companion files:** `adversarial-literature-review-protocol.md` (§5 gate), `docs/mechanisms.md` (switch-variable table), `paper-draft-outline.md` (H1–H4).

---

## §1 Closed mechanism vocabulary (no free association)

Candidate mechanisms can only be composed from the primitives below. **7 are implemented in the shipped kernel; P8 (delay) is NOT — it is a proposed code extension.** Any "new mechanism" that cannot map to an existing switch belongs to the **to-be-implemented components** list, not the candidate pool.

> **API correction.** The "CERPS landing" column uses the **real** parameter names from `cerps_model/src/cerps.py` (consumed by `run_cerps(**params)`). Earlier class names (`CERPSConfig`, `GRNConfig`, `RatesConfig`, `MemoryConfig`, `EpiConfig`, `DivisionConfig`, `SpaceConfig`, `delay_steps`) do **not** exist in the code. `legacy_config()` (in the `cerps_sbi` adapter) returns the calibrated anchor = `CERPS_DEFAULTS` + locked spatial `BASE`.

| Primitive | Biological meaning | CERPS landing (real param) | Status |
|---|---|---|---|
| **P1** niche feedback | niche signal regulates proliferation | `h_ref_S` / `rho_S` / `KS` / `hS` (S-block) | implemented |
| **P2** cycle driver | cell cycle / oscillation | `beta0` / `kappa0` / `tau0` / `h_term` (continuous proliferation) | implemented; periodic forcing = proposed |
| **P3** component dilution | division partitioning noise | `s_part` (expression partitioning at mitosis) | implemented |
| **P4** quorum sensing | density dependence | `crowding` (local occupancy cap) | implemented (coarse) |
| **P5** state inheritance | epigenetic | `u` + `k_maint` + `w_u` + `D_u` | implemented |
| **P6** non-Markovian memory | history dependence | `epa` accumulator gated by `kappa_a` (+ `div_death_age`) | implemented (as age gate) |
| **P7** spatial coupling | paracrine / contact inhibition | `grid` / `L` / `motil` / `crowding` + diffusing `S` | implemented |
| **P8** feedback delay | delayed feedback | — (no `delay_steps` in code) | **NOT implemented — proposed extension** |

---

## §2 Nominate active primitives from public data (data → primitive mapping)

The **source of candidates is data, not imagination.** Reverse-nominate by data type to obtain the "active-primitive subset" of this system — only primitives falling in this subset enter the candidate pool.

| Data type (public repository) | Hinted biological signal | Nominated primitive |
|---|---|---|
| scRNA pseudotime / RNA velocity (GSE81682 etc.) | temporal transition and direction | P2, P8, P5 |
| scATAC-seq (GSE53359 / GSE96772) | chromatin plasticity, memory | P5, P6 |
| spatial transcriptomics + CellPhoneDB / NicheNet | ligand-receptor, neighborhood signal | P1, P7, P4 |
| lineage tracing / CARLIN | division-coupled inheritance and dilution | P3, P5 |
| population oscillation / time series | periodicity | P2, P8 |

> Each nomination must note the **specific dataset + observed phenomenon that supports the primitive** (e.g., "GSE81682 X→Y transition accompanied by delayed up-regulation of Z → nominate P8"); otherwise it is an unsupported nomination and is discarded.

---

## §3 Compose falsifiable candidates (two-primitive is the main productive zone)

Compose per discipline:
- **Single primitive:** only as ablation control, not a standalone candidate;
- **Two-primitive interaction:** main productive zone, prioritize generating candidates. Examples:
  - **P1×P8:** niche feedback with delay
  - **P5×P6:** epigenetic inheritance + non-Markovian memory
  - **P3×P5:** division dilution + state inheritance
  - **P1×P7:** niche signal + spatial coupling
- **Three-or-more primitives:** require a dimensionality-reduction plan first (which two are main effects, how the third is fixed), otherwise not listed.

Each candidate must be written as a falsifiable H (to fill into H1–H4 of `paper-draft-outline.md`):
> **Hk:** if mechanism **{P_i × P_j}** drives phenomenon **X**, then CERPS under {parameters} should observe **D** (consistent / inconsistent with data).

---

## §4 Triple scoring for prioritization (pick "worth doing")

Score each candidate 1–5 on three dimensions; weight to get a priority score.

| Dimension | Meaning | Evidence source |
|---|---|---|
| **Novelty** (w=0.4) | residual gap after adversarial review (⚪ high / 🟡 medium / 🟠 low) | §5 adversarial gate |
| **Data support** (w=0.35) | is there a public dataset to calibrate / validate this primitive | §2 nomination record |
| **Falsifiability** (w=0.25) | can CERPS cleanly isolate it by ablation (switch exists + `health_check` passes) | `mechanisms.md` + batch |

> **Rule of thumb:** high novelty × has data × isolatable → prioritize into Results; novel but no data validation → downgrade one level, only as Discussion outlook.

**Priority score = 0.4·Novelty + 0.35·Data support + 0.25·Falsifiability**, sorted descending.

---

## §5 Adversarial gate (final checkpoint)

Before a candidate enters the paper, run `adversarial-literature-review-protocol.md`:
1. Exclude Bokes 2012 / Jia 2018 / Jia 2024 / Wu–Grima–Jia 2025 (stochastic-expression line);
2. Exclude Duddu 2020/2026, Van Liedekerke 2019, HSL2024 (differentiation / framework line);
3. Output a `novelty_register` entry + ⚪🟡🟠🔴 verdict; **any 🔴/🟠 without a written difference argument must not enter the manuscript**.

---

## §6 Reconcile with CERPS capability (prevent dead ends)

- Check `docs/mechanisms.md` for orthogonality risks, avoid "pseudo-comparisons";
- Confirm candidate primitives are **implemented** (P6/P8 are in; an unimplemented nomination becomes a build task);
- Use `batch.health_check` to pre-judge numerical feasibility (CFL / `dt` convergence).

---

## §7 How the assistant (AI) helps you select candidate mechanisms

The following are actions I can **actually execute** in this process (distinguishing "what I can do" from "what needs your decision"):

### What I can do directly
1. **Data → primitive nomination:** given an accession list, I read / analyze its description and published conclusions, fill the §2 mapping table, and produce the "active-primitive subset".
2. **Code reconciliation:** `grep` / `read` `cerps/config.py`, `docs/mechanisms.md`, confirm candidate primitives are implemented, list switch positions and defaults, flag orthogonality risks.
3. **Combination generation:** per §3 rules, auto-generate the two-primitive candidate list (P1×P8, P5×P6, etc.), removing unimplemented / three-primitive-without-reduction ones.
4. **Scoring matrix:** fill the §4 table with your and my judgments, compute priority scores, output sorted.
5. **Falsifiable-H drafting:** write each candidate in H-form, ready to drop into the draft outline.
6. **Adversarial-review execution:** run the queries of `adversarial-literature-review-protocol.md` (WebSearch / WebFetch), check your personal red lines and domain priors, fill `novelty_register`, give ⚪🟡🟠🔴.
7. **Numerical-feasibility pre-judgment:** use `batch.health_check` thinking to confirm the candidate is isolatable, CFL / `dt` converges.
8. **Checklist & tracking:** build task lists, output the candidate register (MD / CSV), update as iterations proceed.

### What needs your decision (non-delegable)
- **Domain judgment:** the primitive's "expected effect direction", "is this mechanism biologically plausible" — your professional judgment; I only give `[inferred]` and mark for confirmation.
- **Data authorization:** downloading / processing large omics data (network and compute) needs your authorization and environment.
- **Final novelty verdict:** whether 🔴/🟠 constitutes "substantive scooping" is your call.
- **Selection decision:** when candidates are many, which to keep in the paper.

---

## §8 Worked example (HSC-differentiation candidate pool, based on completed audit)

| Candidate | Combination | Nomination basis (data) | Novelty | Data support | Falsifiable | Score | Gate |
|---|---|---|---|---|---|---|---|
| C1 niche feedback w/ delay | P1×P8 | scRNA transition w/ delayed up-regulation | ⚪ | 4 | 5 | 4.45 | to run |
| C2 epigenetic inheritance + memory | P5×P6 | scATAC plasticity + CARLIN | 🟡 | 4 | 5 | 4.35 | to run |
| C3 dilution + inheritance | P3×P5 | lineage tracing | 🟡 | 3 | 5 | 3.95 | to run |
| C4 repression network (control) | single / double P? | 3-gene measured single-on ~93% | ⚪ | 5 | 5 | 4.75 | partial |

> C4 (repression) is novel and data-strong, but the audit shows its "all-on = stemness" metaphor is refuted — so it becomes a **control / negative finding**, not a mechanism candidate.

> **This §8 placeholder is superseded by the real-evidence `candidate_register.md`** (populated 2026-10-19 with GSE81682/GSE89433/GSE53359/GSE96772/CARLIN/RNA-velocity findings and triple-scored priorities). The table above is retained for protocol illustration only.

---

## §9 Outputs
- `candidate_register.md`: the **real-data-populated** candidate-mechanism register (§2 map + §4 matrix + §8 example, with GEO/CARLIN evidence);
- `novelty_register.md`: the adversarial-gate output for every candidate (verdicts ⚪🟡🟠🔴, red-line handling);
- each pooled candidate maps to one H sentence (into the draft outline);
- build-task list for any unimplemented primitive (if present).
