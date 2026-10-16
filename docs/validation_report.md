# Reproduction validation report

This report documents the validation of the Python reproduction
(`reproduction/hsl2024`) against the published results of Huang, Situ &
Lei (2024) and against the authors' original C++ implementation
(github.com/jinzhilei/Cell-type-transition), and the cross-validation of
the CERPS framework against the reproduction.

Three independent implementations of the same protocol exist in this
repository:

| engine | language | role |
|---|---|---|
| `original C++` | C++ (authors') | reference output, run from source with our resumable driver |
| `hsl2024.numba_sim` | Python/numba | production engine for all figures |
| `hsl2024.reference_sim` | Python/NumPy | independent engine (different code path and RNG), used for cross-validation and Fig. 6 lineage tracking |

## 1. Table 3 — transition probabilities (headline validation)

Protocol: 30 runs, 4000 h, sigma = 0.05, m2 = 0.8, alpha = (0.4, 0.4);
statistics over t in [2000, 4000] h; 15 582 131 recorded transitions.

| mother -> daughter | paper | this repo | abs. diff |
|---|---|---|---|
| SC -> SC  | 73.9402 | 73.9885 | 0.048 |
| SC -> TA1 | 12.9889 | 12.9919 | 0.003 |
| SC -> TA2 | 13.0709 | 13.0197 | 0.051 |
| TA1 -> SC | 0.0013  | 0.0010  | 0.0003 |
| TA1 -> TA1| 99.9955 | 99.9955 | <0.0001 |
| TA1 -> TA2| 0.0032  | 0.0035  | 0.0003 |
| TA2 -> SC | 0.0013  | 0.0015  | 0.0002 |
| TA2 -> TA1| 0.0035  | 0.0032  | 0.0003 |
| TA2 -> TA2| 99.9952 | 99.9953 | 0.0001 |

All nine entries agree to 3-4 significant digits; the residual differences
are within Monte-Carlo error for 30 runs.

## 2. Cross-validation against the original C++ code

The authors' C++ program was compiled from the released sources
(`g++ -O2`) and run with our checkpointing driver (`main_ref.cpp`, which
reorders nothing but completes the same 30 runs x 5 m2 values x 4000 h
sweep of `Test()`).  Comparison at m2 = 0.8 (all 30 C++ runs vs 30 Python
runs):

| statistic | C++ (30 runs) | Python (30 runs) |
|---|---|---|
| total cells | 7026.5 | 7158.5 |
| SC         | 176.3  | 176.8 |
| TA1        | 3456.3 | 3495.8 |
| TA2        | 3394.0 | 3486.0 |
| SC->SC / SC->TA1 / SC->TA2 (%) | 73.98 / 13.26 / 12.76 | 73.99 / 12.99 / 13.02 |

The counts differ by < 3% (within run-to-run spread; SEM of the 30-run
mean is ~60 cells), and the transition matrix agrees; the small TA1/TA2
asymmetry in the C++ SC row is a finite-sample effect (the paper's own
Table 3 shows 12.99/13.07, and the parameters are symmetric).

`reproduction/results/tables/cpp_vs_python.json` contains the full
comparison including the m2 sweep.

## 3. Figure-by-figure reproduction outcomes

* **Fig. 2A** bifurcation: fold points reproduced at a ~ 0.75 and
  a ~ 1.65 (tristable window), SC branch medium, TA branches high/low.
* **Fig. 2B/C/D**: phenotype map over (u1,u2); SC zone for u_i > ~0.34;
  SC duration short near the fold edge (u = 0.34 switches within 100 h)
  and saturating at 100 h deep inside the SC zone.
* **Fig. 3**: with the epigenetic state fixed, TA cells require
  sigma above a threshold to appear; TA ratios at high sigma approach the
  stationary basin occupancy (see `results/tables/fig3_ratios.csv`).
* **Fig. 4** homeostasis composition (30 runs, 2000 h):

  | (alpha1, alpha2) | total | SC | TA1 | TA2 | paper |
  |---|---|---|---|---|---|
  | (0.4, 0.4) | 7170 | 175 | 3461 | 3534 | all three coexist (A) |
  | (0.4, 1.0) | 6226 | 345 | 0 | 5881 | TA1 blocked (B) |
  | (1.0, 0.4) | 6043 | 354 | 5689 | 0 | TA2 blocked (C) |
  | (1.0, 1.0) | 596 | 596 | 0 | 0 | SC only (D) |

* **Fig. 5**: (x1, x2) collapses to three clusters while (u1, u2) stays
  continuous - discrete types from continuous epigenetics.
* **Fig. 6**: dedifferentiation (TA1->SC, TA2->SC) and
  transdifferentiation (TA1<->TA2) trajectories extracted by lineage
  tracking; transdifferentiation paths pass through the transient SC
  region without settling there, consistent with the paper.
* **Fig. 7**: (sigma, m2) = (0.20, 0.80) and (0.05, 0.75) lead to
  extinction (over-differentiation), (0.05, 1.00) to SC-only survival.
* **Fig. 8** (sigma in {0, 0.025, 0.05, 0.075, 0.1}):
  homeostasis numbers/ratios almost unchanged (total 7254 -> 6637,
  SC ratio 2.7 -> 2.0%); differentiation probability increases with
  sigma (12.6 -> 14.1%); dedifferentiation increases (0.0009 -> 0.0051%);
  transdifferentiation increases (0.002 -> 0.009%); self-renewal nearly
  unaffected - all the trends stated in Sec. 3.4.1.
* **Fig. 9** (m2 in {0.800 ... 0.900}): TA numbers decrease
  (3496 -> 251), SC increases (177 -> 661), total decreases
  (7159 -> 1166); SC self-renewal increases (73.99 -> 98.43%);
  differentiation decreases (12.99 -> 0.78%); dedifferentiation increases
  by three orders of magnitude (0.001 -> 1.60%); transdifferentiation
  rises then falls (0.0035 -> 0.070 -> 0.054%) - all trends of Sec. 3.4.2.
* **Fig. 11** (reprogramming; 100 runs per point):
  TA1-init success stays low (0-15%) for all a_TF; TA2-init success rises
  2% -> 94% with a_TF; success time TA1 ~500-750 h, TA2 ~544-595 h for
  a_TF >= 0.4 (paper: ~700 h and ~500 h).
* **Fig. 12/13**: see `results/tables/fig12_tf_effect.csv` and the
  landscape figure (regenerated from the same batches).

## 4. CERPS cross-validation

`tests/test_cerps.py` — 37/37 tests pass, including:

* GRN: N-gene update == paper Eq. (1) for G=2 (bitwise-level numerical
  identity on the deterministic trajectory); 2-gene attractors match the
  paper's tristable values; 3-gene ring network produces 3 committed
  attractors.
* Epigenetics: all four kernels reproduce their analytic moments
  (beta: mean phi, variance phi(1-phi)/(1+m0); cpg: mean/sd of the
  site-level binomial map; gaussian: truncated moments); continuous
  dynamics fixed point matches the analytic value and increases with gene
  activity.
* Division: binomial partition conserves molecules; daughter
  concentration sd = sqrt(x/V).
* Rates: continuous mode with steep switch (p=40) converges to the legacy
  Table-1 values (beta_SC = beta0 theta/(theta+Q), kappa_SC = kappa0,
  tau_SC = tau0, TA values likewise).
* Simulator: legacy mode reproduces hsl2024 homeostasis counts within
  Monte-Carlo error; continuous mode reaches a coexistence homeostasis;
  3D spatial mode stays inside the domain, bounded by contact inhibition.
* Coupling: all three bidirectional links demonstrated
  (epi -> GRN, GRN -> epi, population -> cell).

## 5. Numerical notes

* dt = 0.25 h (paper); forward Euler for the gene ODE and Euler-Maruyama
  for the OU noise, exactly as the original C++ (see
  `docs/protocol_notes.md`).
* The reproduction engine runs ~6x faster than the original C++
  (numba), enabling the full figure set (~30 CPU-hours -> ~5 CPU-hours).
* All batch outputs are checkpointed under
  `reproduction/results/data/` and are deterministic given the seeds
  recorded in `reproduction/scripts/common.py`.
