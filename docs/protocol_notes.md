# Protocol notes: paper text vs released code

The released C++ code (github.com/jinzhilei/Cell-type-transition, GPL-3.0)
is the ground truth used for the Python reproduction.  Where the paper's
text and the code differ, we follow the code; all differences are listed
here.  These details matter: they are required to reproduce Table 3 to
four significant digits.

## Phenotype classification (Eq. 5 vs UpdateCellType)

Paper: SC = both genes medium; TA1 = X1 high, X2 low; TA2 mirrored;
TC = otherwise (with all kinetic rates zero, Table 1).

Code: hard thresholds at x = 0.5.  SC: x1>0.5 & x2>0.5; TA1: x1>0.5 &
x2<0.5; TA2: x1<0.5 & x2>0.5.  When both x < 0.5 the cell **keeps its
previous type** - there is no TC phenotype and no zero-rate state.  The
classification is re-evaluated only during the resting phase (stages 0/3).

## Division counter and senescence (Sec. 2.2 vs CellFateDecision)

Paper: the counter increments when a TA cell enters the proliferating
phase; kappa = infinity once divisions >= 15.

Code: the counter increments **at mitosis**; when it reaches exactly 15
and the cell is TA, the cell **dies at that mitosis instead of dividing**
(there is no kappa = infinity branch).  A resting cell classified as SC
has its counter reset to 0 (dedifferentiation restores replicative
potential).

## Stage machine

0 resting; 1 = entered proliferation in the previous step (transition
records are written for stage-1 cells, i.e. once per cell cycle at
proliferation entry); 2 = proliferating; 3 = newborn daughter (behaves
as resting on the next step; its phenotype is re-evaluated then).

## Transition probabilities (Table 3 / .tran files)

Recorded at every step with `Step > 8000` (i.e. t > 2000 h of a 4000 h
run): the pair (type at proliferation entry, parent type at the mitosis
that produced the cell).  The Matlab analysis row-normalises per run,
then averages percentages across runs (DataAnalysis.m).

## Homeostasis statistics

Per-run time-average of each count over .sys rows 8000:end (t >= 2000 h),
then mean and std across the 30 runs; ratios are computed per run from
the time-averaged counts, then averaged (DataAnalysis.m).

## Noise integration

Euler-Maruyama for the OU processes (not the exact transition), the
expression rates are evaluated with the *new* eta each step, and the gene
ODE is advanced with one forward-Euler step of dt = 0.25 h.

## Daughters

Both daughters inherit the mother's eta, age (reset at the next
proliferation entry), division counter, and expression state; only the
epigenetic state is redrawn (independently per gene and per daughter)
from `Beta(m0 phi, m0 (1-phi))` with phi evaluated at the mother's state.

## Parameters of the sweeps

* Fig. 8: sigma in {0, 0.025, 0.05, 0.075, 0.1}; Fig. 9: m2 in
  {0.800, 0.825, 0.850, 0.875, 0.900} (the released main.cpp implements
  exactly this sweep with m2 = 0.8 + 0.025 (index2 - 1)).
* Fig. 10: the same m2 sweep, six parameter series (default and
  +10% single-parameter increases k1=1.1, alpha1=0.44, rho1=0.11,
  b1=1.1, s1=0.55).
* Fig. 11: 500 initial cells, 100% TA1 or 100% TA2, a_TF in
  {0, 0.1, 0.2, 0.3, 0.4, 0.5}, success = 100 SC cells appearing within
  10 000 h (100 runs per point; caption reports 600 runs = 6 points x 100).
  The initial epigenetic state is the TA fixed point of the inheritance
  map (phi(v*) = v* = 0.1325 at m2 = 0.8) - "aged" TA cells; this choice
  is required to reproduce the low basal reprogramming probability.
* Fig. 12: a_TF in {0, 0.05, 0.1, 0.15, 0.2}, 100 initial SC cells,
  30 runs, statistics over t in [2000, 4000] h.
* Fig. 13: U(t, z) = -log f(t, z), z = x1 - x2, time window 0-1000 h,
  pooled over the 30 runs of the sweeps above.

## Reproducibility of the reference runs

`original_repo` (not redistributed here; see the README) runs 30 runs
x 5 m2 values x 4000 h with `Test()`.  We also provide a resumable driver
(main_ref.cpp, same model code) that checkpoints completed runs.
