# CERPS Paper — Novelty Register (adversarial literature review, real-search populated)

> **Purpose.** The DoD gate for Stage 1 of `paper-workflow-playbook.md`. Per `adversarial-literature-review-protocol.md`, every novelty/priority claim is stress-tested by *actively seeking refuting literature*. Searches were run on 2026-10-19 via PubMed / Google Scholar / bioRxiv / arXiv across the four source types. Each claim carries a verdict ⚪🟡🟠🔴.
>
> **Principle (from the protocol).** "Prove before assert." Numerical / literature-supported → `[derived]`; numerical-only, not strictly proven → `[to-be-verified]`; conjecture → `[conjecture]`. No `[to-be-verified]`/`[conjecture]` may be written as a proven conclusion.

---

## 0. Claim inventory (mapped to candidate_register.md)

| # | Claim | Candidate |
|---|---|---|
| N1 | Delayed niche feedback (P8) explains HSC quiescence oscillation / clone-size variance | C1 |
| N2 | Epigenetic state carries non-Markovian memory (P6) beyond current state | C2 |
| N3 | Division partitioning noise (P3) is a significant, isolatable fate-determining input in HSC | C3 |
| N4 | A purely mutual-repression G≥3 network **cannot** robustly reproduce "all-on = stemness" (measured single-on ~93%, all-on ~0.5–1%) | C4 |
| N5 (meta) | CERPS as a *controlled causal-inference sandbox* coupling public data to mechanism dissection is a distinctive workflow | framework |

---

## 1. Personal red-line check (universal — Bokes / Jia / Wu–Grima–Jia)

These authors define the **stochastic gene-expression / CME-spectral** line. Our paper is a *differentiation-mechanism causal screen*; collision is therefore expected to be **adjacent, not direct**.

| Author / Year | Verified topic (from real search) | Collision with our claims? | Verdict | Difference argument (our boundary) |
|---|---|---|---|---|
| **Bokes 2012** (J Math Biol, PMID 21656009) | Exact/approx steady-state distributions of mRNA & protein in the **low-copy regime** via the **chemical master equation** (generating-function approach) | No — different subfield | 🟡 | We do **not** analyse the CME spectrum of gene expression; our stochasticity is **extrinsic** (division partitioning + measurement), not intrinsic master-equation noise. No spectral-gap claims. |
| **Jia 2018** (J Chem Phys 148:095102, PMID-per-arXiv 1803.01123) | "Relaxation rates … reveal the feedback signs" via the **spectral gap of the Markov model** of single-cell stochastic gene expression | No — different subfield | 🟡 | Same boundary as above; we use forward simulation + SBI, not CME eigenvalue/spectral-gap analysis. |
| **Wu–Grima–Jia 2025/2026** (J Math Biol, "Spectral theory of stochastic gene expression: a Hilbert space framework", PMID 41733667; bioRxiv 2025.05.19.654994) | Full **CME spectral theory** — eigenvalues, spectral gap, relaxation rate, Hilbert-space framework for constitutive/bursty/autoregulated gene expression | No — different subfield | 🟡 | Explicitly state in Methods/Discussion that our landscape & mechanism claims are **simulation-based and semi-quantitative**, orthogonal to CME spectral exact solutions. Cite them as the *adjacent* gold standard we do not compete with. |

> **Action.** Add a one-paragraph "Relation to the stochastic-gene-expression spectral line" boundary statement to the Discussion, citing Bokes 2012, Jia 2018, Wu–Grima–Jia 2025/2026. No 🔴/🟠 collision → these red-line authors do **not** block any claim.

---

## 2. Domain-prior check (CERPS / HSC context, from the audit)

| Author / Year | Known result | Overlap with our claim | Verdict | Our difference (to write) |
|---|---|---|---|---|
| **Duddu 2020 / 2026** (J R Soc Interface 23:20251065) | Multi-gene GRN multistability: single-on ~60% / double-on ~36% / triple excluded; product combination | N4 (attractor statistics) | 🟡 | Our measurement is **opposite in distribution direction** (single-on ~93%, all-on ~0.5–1%), positively *refuting* the "all-on = stemness" metaphor; this is the point of C4, not a collision. |
| **Van Liedekerke 2019** | Sharp threshold incompatible with data | N5 (model assumption) | 🟡 | We state continuous mode is a modeling assumption, not empirical fact; anchor on `legacy_config()`. |
| **HSL2024** (J Theor Biol 577:111664) | Parent framework of CERPS | N5 (framework novelty) | 🟡 | CERPS is a generalization; declare inheritance explicitly. |
| **Boyer 2005** | Multipotent core is mutually-activating | N4 (stemness mechanism) | 🟡 | We propose a layered (repression + activation) structure; C4 supports the activation-core view. |
| **Berry 2017 / Miyamoto 2015** | Real ODE basis for epigenetic read/write | N2 (epigenetic kernel) | 🟡 | We replace erroneous Dodd/Probst/Zhao citations; use Berry/Miyamoto as the correct basis. |

No 🔴/🟠 without written argument → all domain priors are **controllable**.

---

## 3. Claim-by-claim adversarial results

### N1 — Delayed niche feedback (P8) in HSC  [Candidate C1]
- **Search.** `"delayed feedback stem cell niche"`; `"time delay" OR "delayed feedback" stem cell differentiation`.
- **Real hit.** **Dray et al. 2021**, *Cell Stem Cell* 28(8):1457–1472 (PMID 33823144): in **neural** stem cells, a dividing-NSC dispersion effect is exerted with a **9–12 day delay**; a stochastic lattice model reproduces self-propagating homeostasis.
- **Verdict.** 🟡 adjacent, controllable.
- **Difference argument.** Dray et al. is **neural**, not hematopoietic, and models a *physical dispersion* delay; we model a **transcriptional/niche-feedback delay** in HSC quiescence–activation oscillation. **[Honesty flag: `delay_steps` is NOT in this CERPS version — the P8 delay axis is a proposed code extension, not yet simulated; only the P1 niche axis is verifiable.]` Cite Dray 2021 as the conceptual precedent and explicitly delimit system & feedback modality. HSC-specific delayed-feedback instantiation remains open.
- **Status.** `[verified-prod §7: P1 axis verified at production scale — term_death +145% (124k→305k), N_final −26%; P8 delay still requires the delay_steps code addition]`.

### N2 — Non-Markovian epigenetic memory (P6)  [Candidate C2]
- **Search.** `"epigenetic memory hematopoietic stem cell fate bias"`; `"non-Markovian history-dependent fate"`.
- **Real hits.** 2024 *Trends Cell Biol.* review (epigenetic landscape differs between HSCs of distinct lineage bias, detectable *before* transcription); Nerlov 2023 *Nat Cell Biol.* (lymphoid enhancers as "latent memory"); RNA velocity (La Manno 2018) **explicitly assumes Markovian** fate.
- **Verdict.** 🟡 adjacent, controllable.
- **Difference argument.** The *concept* of epigenetic memory in HSC is established; our novelty is the **quantitative mechanism-screen** — using CERPS `kappa_a` (the epigenetic-age gate, i.e. P6) to test whether a non-Markovian memory term is *necessary and sufficient* to explain CARLIN within-clone correlation, *challenging the Markovian assumption* of RNA velocity. Position as "mechanism quantification of a known bias," not "discovery of epigenetic memory."
- **Status.** `[verified-prod §7: kappa_a↑ → N_final −83%, cyc +78% — P6 has a large, robust effect]`.
- **⚠ CORRECTION (production scale).** The v1 observation of *seed-dependent lineage polarization* is **refuted**: with three replicates the outcomes **converge** as `kappa_a` rises (cross-seed sd of mG 0.922 → 0.149, of mP 0.944 → 0.106). The defensible formulation of N2 is therefore: *the epigenetic-age (P6) term acts as a divisional-age selective filter that **suppresses** stochastic divergence between replicate histories* — the outcome still depends on accumulated divisional history (the non-Markovian ingredient), but it constrains rather than amplifies replicate divergence. **Do not** write "opposite fates from identical parameters".

### N3 — Division partitioning noise (P3)  [Candidate C3]
- **Search.** `"partitioning noise cell division fate decision stem cell"`; `"stochastic partitioning hESC NANOG"`.
- **Real hits.** Wu & Tzanakakis 2012 (*PLoS One* 7:e50715, PMID 23226362): stochastic partitioning ≈ 17% of hESC NANOG noise; Soltani 2016 (audit-cited): analytical binomial partitioning.
- **Verdict.** 🟡 adjacent, controllable.
- **Difference argument.** Partitioning noise theory is established in hESC; our angle is **coupling it to HSC lineage bias** via CERPS `s_part` (partitioning noise at mitosis, P3) and testing isolatability with a permutation necessity test. Cite Wu & Tzanakakis / Soltani as foundation, claim the *HSC-fate-coupling* as ours.
- **Status.** `[verified-prod §7: s_part 0→0.15 shifts all bulk metrics by only ~1–4% (vs 145% for C1) — honest bulk-null confirmed at production scale; needs clone-variance metric]`.

### N4 — Mutual repression cannot give "all-on = stemness"  [Candidate C4, negative finding]
- **Search.** `"mutual repression network stemness all-on state"`; `"mutual antagonism Oct4-Cdx2 Nanog-Gata6 pluripotency"`.
- **Real hits.** Chang 2013 (*PLoS Comput Biol* 9:e1003165): 52-node hESC network whose bistability is driven by mutual antagonism (Oct4–Cdx2, Nanog–Gata6); PMC2566811: trophectoderm/endoderm fates from mutual antagonism. Duddu 2020: single-on dominant ~60%.
- **Verdict.** 🟡 adjacent, controllable (quantitative refutation).
- **Difference argument.** The qualitative "repression alone insufficient; activation core needed" (Boyer 2005) is known; our contribution is a **quantitative, multilocus measurement** (single-on ~93% / all-on ~0.5–1%, opposite to Duddu 2020) showing the "all-on = stemness" metaphor fails under full mutual repression. C4 is a *negative/control finding*, not a positive mechanism claim → low scooping risk.
- **Status.** `[verified-prod §7: GRN b-sweep term_death 14 ± 5 → 420251 ± 9185. CORRECTED: frac_both_on ≈ 0.017–0.019 (individual runs 0.000–0.050, seed 21 always 0) — it is NOT 0. Negative finding upheld at the ~2% level]`.

### N5 — CERPS as controlled causal-inference sandbox (framework)  [meta]
- **Verdict.** ⚪ claimable (no direct competitor found in the searched window combining forward generative model + public-data SBI + mechanism ablation for HSC).
- **Difference argument.** Purely data-driven trajectory inference (e.g., RNA velocity, Palantir) lacks *controlled causal intervention*; traditional ODE/ABM lacks *public-data calibration*. CERPS couples both. Declare HSL2024 inheritance.
- **Status.** `[derived]` (methodological positioning; supported by the audit).

---

## 4. Summary verdict table

| Claim | Verdict | Blocks manuscript? | Required action before submission |
|---|---|---|---|
| N1 P8 delay | 🟡 | No | cite Dray 2021, delimit neural→HSC; **P1 niche axis verified-prod §7 (term_death +145%, N_final −26%); P8 delay still proposed-unimplemented** |
| N2 P6 memory | 🟡 | No | reframe as *variance suppression by divisional-age selection*, not polarization; **verified-prod §7 (N_final −83%, cyc +78%; replicate sd collapses 0.92→0.15)** — v1 polarization claim refuted |
| N3 P3 partitioning | 🟡 | No | cite Wu&Tzanakakis 2012 / Soltani 2016; **verified-prod §7 (bulk-null, ~1–4%; needs clone-variance metric)** |
| N4 mutual-repression refutation | 🟡 | No | **verified-prod §7: both-on ≈ 1.7–1.9% (NOT 0) across GRN strengths; negative finding upheld at the ~2% level**; keep the separate `multigene_attractors.md` single-on ~93% figure with its finer-threshold caveat |
| N5 framework | ⚪ | No | declare HSL2024 inheritance |
| Red lines (Bokes/Jia/Wu–Grima–Jia) | 🟡 | No | add boundary paragraph in Discussion |

**No 🔴 (substantive scooping) and no un-argued 🟠 anywhere.** All claims may proceed to simulation verification.

---

## 5. Follow-up (full-text reads required before final submission)

Per the "read the full text, not just the abstract" rule, the following should be opened in full before the manuscript is finalized (abstracts consulted so far are sufficient for candidate selection, not for final citation):
- Dray et al. 2021 (Cell Stem Cell) — full methods on the 9–12 day delay model.
- Duddu 2020 / 2026 — exact attractor-statistics methodology to confirm the distribution-direction contrast.
- Wu–Grima–Jia 2025/2026 — confirm scope is CME spectral (it is, per abstract) so the boundary paragraph is accurate.
- Nerlov 2023 (Nat Cell Biol, s41556-023-01137-5) — "latent memory" enhancer evidence for C2.
