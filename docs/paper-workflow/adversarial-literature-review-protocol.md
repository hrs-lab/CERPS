# CERPS Paper Adversarial Literature-Review Protocol (search queries + exclusion checklist)

> **Purpose.** Before any "novelty / priority" claim in the paper is committed, this adversarial review MUST be completed. It is the DoD gate for Stage 1 of `paper-workflow-playbook.md`.
> **Principle.** Not "find literature that supports us", but "**actively seek literature that could refute our claim**". The goal is to falsify our own claim; only if we cannot find a refutation does the claim stand.
> **Scope.** This file carries both (a) your long-standing personal red lines (excluding known results Bokes 2012 / Jia 2018 / Jia 2024 / Wu–Grima–Jia 2025) and (b) the domain-prior checklist in the CERPS / HSC-differentiation context.

---

## 1. Operating procedure (5 steps, each logged)

1. **List the claim inventory.** Write each novelty point the paper will assert as a single falsifiable statement (e.g., "We are the first to use feedback delay P8 to explain HSC lineage bias").
2. **Run adversarial searches.** For each claim, run the queries in §2 across the 4 sources (PubMed / Google Scholar / arXiv / bioRxiv).
3. **Compare against the exclusion checklist.** Tick each row of the §3 table, recording "collision? + where is the difference".
4. **Write the difference argument.** For each potential collision, write 1–2 sentences on "our **essential difference** from X" (mechanism / scale / data / conclusion); otherwise downgrade the claim.
5. **Archive.** Write results into `novelty_register.md` (a living document), updated as the manuscript evolves.

---

## 2. Concrete search queries (paste-ready)

> Replace `{TOPIC}` with the current mechanism keyword (e.g., `feedback delay`, `epigenetic memory`, `niche`, `Waddington landscape`, `stemness threshold`).

### 2.1 General adversarial search (English, covers all claims)
```
("cell differentiation" OR "hematopoietic stem cell" OR HSC)
AND ("gene regulatory network" OR GRN OR "epigenetic" OR "niche" OR "signal pathway")
AND ({TOPIC})
```
```
("mechanism" OR "causal" OR "perturbation") AND ("single-cell" OR scRNA) AND ({TOPIC})
```
```
{TOPIC} AND ("agent-based" OR "cellular Potts" OR "off-lattice" OR "hybrid model" OR "ODE model")
```

### 2.2 Targeted search per CERPS mechanism primitive
| Mechanism primitive | Query keywords |
|---|---|
| P1 niche feedback | `"niche feedback" HSC OR "quorum" stem cell` |
| P2 cycle / oscillation driver | `"oscillatory" "stem cell" fate OR "cell cycle" differentiation` |
| P3 component dilution | `"partitioning noise" OR "molecular titration" stem cell` |
| P4 quorum sensing | `"quorum sensing" cell fate OR "density-dependent" differentiation` |
| P5 state inheritance | `"epigenetic memory" OR "transgenerational" cell fate` |
| P6 non-Markovian memory | `"non-Markovian" OR "history-dependent" cell fate OR "memory kernel"` |
| P7 spatial coupling | `"spatial" "contact inhibition" OR "paracrine field" differentiation` |
| P8 feedback delay | `"time delay" OR "delayed feedback" stem cell OR "lag" differentiation` |

### 2.3 Data sources (must cover preprints to avoid "published but unseen")
- **PubMed / Europe PMC:** published journals
- **bioRxiv / medRxiv:** preprints (high priority-risk zone)
- **arXiv (q-bio.MN / q-bio.CB):** methods / theory
- **Google Scholar:** reviews + highly cited (window with `before:2026` and `since:2020`)

---

## 3. Exclusion checklist (tick every row)

### 3.1 Personal long-standing red lines (universal to all papers, from your originality-check protocol)
| Author / Year | Topic | Collision means | Handling |
|---|---|---|---|
| **Bokes 2012** | gene-expression noise / CME spectral analysis | if the paper asserts "stochastic gene-expression spectrum" | must explicitly distinguish our scale and conclusion |
| **Jia 2018** | stochastic gene-expression modeling | same family | same as above |
| **Jia 2024** | stochastic gene-expression modeling (updated) | same family | same as above |
| **Wu–Grima–Jia 2025** | latest stochastic-mechanism results | same family | same as above; if direct collision, must cite and draw boundary |

> ⚠️ The four authors above mainly correspond to the "stochastic gene expression / CME spectrum" line; **even if the CERPS paper does not touch that line, the domain-prior check in §3.2 must still be completed** — do not skip it.

### 3.2 CERPS / HSC-differentiation domain priors (from the CERPS audit, tick each)
| Author / Year | Known result | Potential overlap with CERPS claim | Our difference (to fill) |
|---|---|---|---|
| **Soltani 2016** | analytical binomial partitioning noise at division | if asserting "partitioning-noise implementation" | we use it for calibration, not claiming primacy |
| **Duddu 2020 / 2026 (J R Soc Interface 23:20251065)** | multi-gene GRN multistability statistics (single-on ~60% / double-on ~36% / triple excluded); product combination (4.5) | if asserting "all-repression → all-on = stemness" | our measurement: single-on dominant ~93%, double-on tail ~7%, all-on ~0.5% (opposite distribution direction) → positively refutes that metaphor |
| **Van Liedekerke 2019** | sharp threshold incompatible with data | if asserting "continuous mode sharp limit = Table 1 = reality" | we explicitly state this is a modeling assumption, not empirical fact |
| **Huang & Lei 2017–2019 (HSL2024, J Theor Biol 577:111664)** | the parent framework of CERPS | any "framework novelty" claim | CERPS is a generalization; must declare inheritance |
| **Boyer 2005** | multipotent core is mutually-activating | if asserting "repression maintains stemness" | we propose a layered (repression + activation) structure |
| **Berry 2017 / Miyamoto 2015** | real ODE basis for epigenetic read/write | epigenetic-kernel claim | replaced erroneous Dodd/Probst/Zhao citations |

---

## 4. Novelty-claim decision matrix

| Search result | Verdict | Action |
|---|---|---|
| Found a **directly identical** claim, earlier | 🔴 substantive scooping | delete the claim, reframe as "we independently verify / extend" |
| Found **partial overlap** (same mechanism, different scale / data) | 🟠 partial overlap | write the difference clearly, downgrade to "supporting evidence" |
| Found **adjacent but controllable** (same topic, different question) | 🟡 adjacent, controllable | cite and draw boundary, claim stands |
| **Not found** (within ~5-year window) | ⚪ claimable | claim stands, mark "[to-be-verified]" |

> Follow your "prove-before-assert" labeling: **numerical / literature-supported → [derived]; numerical-only, not strictly proven → [to-be-verified]; conjecture → [conjecture]**. No [to-be-verified] / [conjecture] may be written as a proven conclusion.

---

## 5. Archive template (`novelty_register.md` excerpt)

```markdown
## Claim N: {one sentence}
- Query: {paste}
- Adversarial literature hit: {list}
- Collision verdict: ⚪ / 🟡 / 🟠 / 🔴
- Difference argument: {our essential difference vs X}
- Status: [derived]/[to-be-verified]/[conjecture]
- Last checked: YYYY-MM-DD
```

---

## 6. Red-line reminders (from your long-standing convention)
- **Adversarial review MUST be completed before claiming novelty**; do not "write first, check later".
- Any 🔴/🟠 in the exclusion checklist **without a written difference argument** → that claim must not enter the manuscript.
- Numerical validation ≠ mathematical proof; strict conclusions must be explicitly labeled, and [to-be-verified] must not be treated as proven.
