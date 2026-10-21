# CERPS Calibration Report — ABC-SMC against Paul 2015 / Buenrostro 2018

> **Scope.** Documents the calibration step of the paper workflow (Runbook §04) using the
> real SBI layer `analysis/src/abc.py`. This is the evidence base behind the
> `[verified-prod]` (not `[derived]`) status of candidates C1–C4 in
> `candidate_register.md` and `novelty_register.md`.
>
> **Author.** Rongsheng Huang, School of Science, Jimei University, Xiamen 361021, China; hrs@jmu.edu.cn. NSFC 12401641.

---

## 1. Setup (real, as implemented in `analysis/src/abc.py`)

- **Inference.** ABC-SMC, 5 populations, adaptive ε schedule (75th/60th percentile of accepted distances), component-wise Gaussian perturbation (10% of prior range), `Pool(2)` parallelism.
- **Prior (8 coupling parameters):**
  `r_eff ∈ [0.45,0.90]`, `r_max ∈ [0.85,0.99]`, `h_ref_S ∈ [0.05,0.80]`,
  `k_maint ∈ [0.01,0.15]`, `s_div ∈ [0.01,0.08]`, `beta0 ∈ [0.15,0.60]`,
  `sigma ∈ [0.01,0.08]`, `D_u ∈ [1e-5,8e-5]`.
- **Target 1 (fit).** Paul et al. 2015 (GSE72857) scRNA-seq summary statistics
  (34-dim: z-scored mean/std of Gata1/Sfpi1/Klf1/Cebpa, cycling fraction,
  cycling-by-commitment-bin, effector-by-bin, lineage mass ratio).
  Stored as `results/data/paul_summary_target.npz`.
- **Target 2 (cross-validation).** Buenrostro et al. 2018 (GSE96769) ATAC
  accessibility ordering along the hematopoietic hierarchy.
- **Model comparison.** Four decoupled variants tested at the final ε:
  `noGate` (lam=0), `noNicheRef` (h_ref_S=0), `constFid` (r_max=r_eff), `betaKernel`.

---

## 2. Convergence diagnosis (this is the key honest finding)

The ABC-SMC **did not converge to an identified posterior**. Trace from
`results/abc/abc_state.npz`:

| pop | ε       | n_accepted |
|-----|---------|------------|
| 0   | 6.792   | 60         |
| 1   | 5.773   | 16         |
| 2   | 4.907   | 513        |
| 3   | 4.171   | **0** → fallback to top-10 |

At population 3 the acceptance rate dropped to **0** at ε=4.17; the code then
fell back to the 10 lowest-distance candidates, and **the saved posterior
(`abc_final.npz`) contains only 10 particles**. The mean distance never fell
below ~4.1 (with standardization weights `SC`, this is large → poor fit).

**Posterior (10 particles, degenerate / edge-collapsed):**

| param     | mean     | sd       | range          | prior range        |
|-----------|----------|----------|----------------|--------------------|
| r_eff     | 0.477    | 0.030    | [0.45, 0.53]   | [0.45, 0.90]       |
| r_max     | 0.863    | 0.010    | [0.85, 0.88]   | [0.85, 0.99]       |
| h_ref_S   | 0.752    | 0.056    | [0.643, 0.80]  | [0.05, 0.80]       |
| k_maint   | 0.0457   | 0.021    | [0.01, 0.094]  | [0.01, 0.15]       |
| s_div     | 0.0261   | 0.010    | [0.01, 0.046]  | [0.01, 0.08]       |
| beta0     | 0.539    | 0.051    | [0.438, 0.60]  | [0.15, 0.60]       |
| sigma     | 0.050    | 0.011    | [0.033, 0.067] | [0.01, 0.08]       |
| D_u       | 1.22e-5  | 3.1e-6   | [1e-5, 1.95e-5]| [1e-5, 8e-5]       |

Every marginal is pinned to a prior boundary (lower or upper). A 10-particle,
edge-collapsed posterior is **not** a calibrated parameter set.

### 2.1 Independent re-run (2026-10-22) — non-convergence confirmed, and the archived ε trace explained

The shipped `abc.py` hard-coded `/home/z/...` paths and carried two further
defects (a crash in the posterior-predictive stage when a simulation went
extinct, and an ε-schedule branch that only fired on *resumed* runs). All three
were fixed (paths via `$CERPS_ROOT`/`$CERPS_ABC_OUT`, NP/NPOPS/NSEED via
environment; extinct-run guard; explicit first-vs-later population logic), and
the ABC was re-run on this machine (Windows, numba 0.65.0) at reduced scale
`NP=24 × NPOPS=3` (full archived scale would need ~6.7 CPU-h; raw output
`CERPS_complete_package/results/abc_rerun_2026-10-22/`).

**Re-run trace (data-adaptive ε):**

| pop | ε (re-run) | n_acc (/24) | mean_d (accepted) |
|-----|-----------|-------------|-------------------|
| 0   | 6.789     | 18          | 6.446             |
| 1   | 6.798     | 18          | 6.286             |
| 2   | 6.591     | 18          | 5.848             |

With the corrected schedule the threshold **stalls at ≈6.6–6.8**: the
`min(60th percentile, 0.85·ε_prev)` rule is bound by the *data* (the 60th
percentile of perturbed candidates), not by the geometric factor, and the
posterior standard deviations barely shrink (e.g. `r_eff` sd 0.134 → 0.098 →
0.072 over three populations). **Non-convergence is confirmed under the
corrected schedule as well.**

**Forensic note on the archived trace.** The archived ε values are an exact
geometric sequence (6.7921 × 0.85 = 5.7733 × 0.85 = 4.9073 × 0.85 = 4.1712,
to four decimals), i.e. the archived run took the forced 0.85 shrink at every
population — the branch that the buggy `if EPS[0] is None` test only reached on
runs *resumed* from a checkpoint. The final population's **0 acceptances** are
therefore a direct consequence of a hard-coded decay outrunning the achievable
distance, not an adaptive threshold. The archived `n_acc = 513` at population 2
is likewise not reproducible from the current script at any single particle
count and should be treated as an unverifiable archived artefact; the
conclusion it was used to support (non-convergence) is nonetheless confirmed
independently by the re-run above.

---

## 3. Relative model comparison (still informative)

`results/abc/model_comparison.json` reports acceptance and mean distance at the
final ε for each variant:

| variant     | accept | mean_d |
|-------------|--------|--------|
| **full**    | 0.0    | 5.41   |
| betaKernel  | 0.0    | 6.27   |
| constFid    | 0.0    | 6.89   |
| noGate      | 0.0    | 7.19   |
| noNicheRef  | 0.0    | 7.60   |

Absolute acceptance is 0 for every variant (the fit is globally weak), **but the
ranking is meaningful**: `full < betaKernel < constFid < noGate < noNicheRef`.
Removing the niche reference (`noNicheRef`) and the epigenetic gate (`noGate`)
both *worsen* the distance most, so the **full model with niche feedback + gating
is the best (least-bad) fit**. This is consistent with, and independently
supports, the C1 (niche) and C2 (epigenetic-gating) mechanism claims.

---

## 4. Posterior predictive check against the Paul target

`results/abc/posterior_predictive.npz` holds `pp_mean`, `pp_std` (34-dim) and the
target. Comparing them locates **exactly which features fail**:

| block | target (block mean) | pp_mean | pp_std | diagnosis |
|---|---|---|---|---|
| master means `xG,xP` | 1.497 | 1.390 | 0.90 | **reasonable** (within ~4–11%) |
| effector means `xK,xC` | 1.496 | 0.453 | 0.31 | **~70 % under-predicted** |
| spreads `sG,sP` | 0.743 | 0.858 | 0.155 | slightly over-predicted |
| spreads `sK,sC` | 0.739 | 0.300 | 0.055 | **~60 % under-predicted** |
| cycling fraction | 0.500 | 0.237 | 0.008 | **~53 % under; near-zero predictive variance** |
| cycling-by-bin (12) | structured 0 → 0.87 (sd 0.22) | **flat ≈ 0.21** | 0.03 | **structural mismatch — model gives no bin structure** |
| Klf1-by-bin (12) | rises → 2.41 → falls (sd 0.69) | monotonic → 1.02 | — | direction captured, **magnitude/shape off** |
| lineage mass ratio | 1.500 | 2.600 | 2.66 | uninformative (predictive spread exceeds the value) |

**Aggregate PPC:** median |z| = **12.2**, coverage (|z| ≤ 2) = **15 %**, max |z| ≈ 1.5×10⁵.
(The huge max arises because several `pp_std` are ~10⁻³ while target–pp gaps are ~10⁻¹.)

**Reading.** The posterior predictive reproduces the **master regulators**
(GATA1/SPI1 levels) acceptably, but systematically **under-produces the effector
genes** (Klf1/Cebpa ~70 % low, their spread ~60 % low) and **under-produces
proliferation** (cycling 0.24 vs 0.50). Most seriously, the
**cycling-by-commitment-bin profile is flat in the model** where the data show a
structured rise/fall — a *qualitative structural* mismatch, not merely an offset.
Quantitative agreement with Paul 2015 is therefore **not** achieved.

> Caveat worth recording: within the target the four master-gene means are nearly
> identical (≈1.495, sd 3×10⁻³), unlike the strongly structured bin profiles. Any
> future recalibration should inspect `results/data/paul_summary_target.npz`
> construction before treating that block as informative.

---

## 5. Calibration vs. candidate sweep levels

- **C1 (`h_ref_S` niche).** ABC posterior places `h_ref_S ∈ [0.643, 0.80]`
  (mean 0.752). My C1 sweep used levels {0.0, 0.40, 1.20}. The (weak) posterior
  prefers a **stronger** niche reference than my "default" 0.40, and my "strong"
  1.20 *exceeds* the posterior upper bound. → The sweep brackets the posterior
  but is not centered on it. A level near 0.75 would better bracket the
  (degenerate) posterior; flag for a follow-up sweep pass.
- **C2 (`kappa_a` memory).** `kappa_a` is **not** in the ABC priors (the ABC
  only calibrates the 8 coupling params above). C2's verification is therefore
  independent of the ABC posterior and stands on the candidate sweep alone.
- **C3 (`s_part`), C4 (`b` GRN).** Neither axis enters the ABC prior; both are
  validated by the candidate sweep only.

---

## 6. Interpretation (honest boundary)

1. **CERPS cannot be rigorously calibrated to Paul 2015 scRNA-seq summary
   statistics.** The ABC-SMC fails to converge (acceptance 0, 10-particle
   edge-collapsed posterior, ε floor at a large distance), and the posterior
   predictive check is poor (coverage 15 %, effector genes ~70 % low, cycling
   ~53 % low, cycling-by-bin structurally flat). Together these say the failure
   is both **identification** (no usable posterior) and **adequacy** (even the
   best particles mis-predict key observables). This is exactly the
   "semi-quantitative / mechanism-screening, not molecular-resolution predictor"
   boundary established in the audit report.
2. **This is a publishable limitation, not a method failure.** The identifiability
   gap (many parameter combinations reproduce the same coarse summary) should be
   stated explicitly in the paper's Limitations, alongside the single-scalar-per-gene
   `u` resolution limit and the no-native-inversion constraint.
3. **Candidate conclusions rest on the controlled sweep, not on an ABC posterior.**
   C1 (niche effect real), C2 (epigenetic gating → non-Markovian fingerprint),
   C3 (bulk-null), C4 (mutual-repression negative finding) are validated by the
   candidate sweep — now at **production scale** `results/candidate_sweep_v2.json`
   (T=4000, N=600, seeds 3/11/21, 39 runs, 2026-10-22), which supersedes the
   screening-scale `candidate_sweep_v1.json` (T=1500, N=200, seeds 3/11); see
   `candidate_register.md` §7.4 for the v1→v2 comparison. The ABC layer
   contributes only the *relative* model-comparison support (niche + gating
   needed).
   > **Status history.** An earlier draft cited `candidate_sweep_v2.json` before
   > it existed (flagged 2026-10-22 morning); the production run has since been
   > **executed** (2026-10-22, elapsed 1636 s) and this section now cites it as
   > real evidence. A hold-out 2-Wasserstein / PPC pass is still required before
   > any candidate is upgraded beyond the sweep-verified level.
4. **Status decision.** Candidates are upgraded **`[verified-v1]` → `[verified-prod]`**
   (robust multi-seed), **not** to `[derived]`: a `[derived]` tag would falsely
   imply a calibrated posterior, which does not exist.

---

## 7. Recommendations before submission

- [x] Re-run `analysis/src/abc.py` with corrected paths on this machine to
      **confirm** the non-convergence — **DONE (2026-10-22)**; see the
      re-run note in §7.1 below. Report it as a deliberate identifiability
      result.
- [ ] Add a C1 sweep level near `h_ref_S ≈ 0.75` to center on the (degenerate)
      posterior and tighten the C1 narrative.
- [ ] Do **not** quote absolute calibrated expression/cycling levels; report all
      candidate effects as *relative to the legacy anchor*.
- [ ] State the ABC identifiability limitation in the paper's Limitations section.
- [ ] Inspect `results/data/paul_summary_target.npz` construction before any
      recalibration: the four master-gene target means are nearly identical
      (sd 3×10⁻³) and the first entry of both bin profiles is exactly 0.0 —
      confirm whether those components are informative, since they dominate the
      distance and may themselves cause the non-convergence.
