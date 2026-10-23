# Multi-axis cross-verification of CERPS — 2026-10-23

**Scope.** Everything in this report was recomputed from the repository state
at `main@02db9d653b` (v1.2.0), i.e. *before* the fixes listed in §3.  The
report is deliberately written as a reviewer's audit: each axis states what was
checked, what method was used, the measured value, and whether the
repository's own claim survived.  Where a claim failed, the failure is
recorded here rather than silently patched.

**Environment.** Python 3.13.12 / 3.13.14, numpy 2.5.4, scipy 1.18.1,
matplotlib 3.11.2, **numba 0.68.0** (the repository pins 0.65.0 — see
§2.7), pytest 9.1.1, Windows/x86-64.  All timings below are for this machine.

---

## 1. What was verified, and how

| # | axis | method | verdict |
|---|---|---|---|
| 1 | test suite actually runs and the claimed counts are real | ran `pytest tests/` and each harness file standalone | **claims confirmed** (69 pytest, 191 standalone) |
| 2 | analytic claims in the code/docs are *numerically* true | 27 independent recomputations, no repo helper reused | **27/27 passed** |
| 3 | doc ↔ code consistency | read every `docs/*.md` claim against the source that implements it | **11 mismatches found** |
| 4 | reproduction data ↔ README | recomputed the C++/Python comparison from the archived JSON | **2 wrong numbers in README** |
| 5 | portability & hygiene | grep for absolute paths, secrets, machine-specific state | **3 hard-coded paths** |
| 6 | packaging / installability | checked for build metadata, imports from a clean tree | **missing entirely** |
| 7 | CI coverage | asked what an automated run would actually execute | **M2 layer invisible** |
| 8 | branch topology / merge safety | compared the two branches tree-by-tree and hash-by-hash | **safe to merge (no-op)** |

---

## 2. Findings

### 2.1 Axis 1 — test counts (claim: 62 pytest, 191 standalone)

Measured:

| command | result | repository claim | match |
|---|---|---|---|
| `python -m pytest tests/ -q` | **62 passed** in 638 s | 62 | yes |
| `python tests/test_cerps.py` | **61 passed, 0 failed** | 61 | yes |
| `python tests/test_audit_fixes.py` | **93 passed, 0 failed** | 93 | yes |
| `python tests/test_audit_2026_10_21.py` | **19 passed** | 19 | yes |
| `python tests/test_audit_ranges.py` | **18 passed** | 18 | yes |
| `python tests/test_m2_generalisation.py` | **ALL M2 CHECKS PASSED** | "V1–V7" | yes |

`test_audit_ranges.py` defines only 4 `test_*` functions but is parametrised
into 18 pytest cases; the two harness files (`test_cerps.py`,
`test_audit_fixes.py`) expose 13/12 of their checks as `test_*`, which is why
62 ≠ 191.  The README explains this correctly.

### 2.2 Axis 2 — independent recomputation of the analytic claims

Written from scratch against `cerps/` (no test helper reused), 27 assertions.
All passed; the two that needed care on the first attempt were *my* estimator
error (a 2 000-chain cross-sectional variance is a 3% estimator, and the
overlap-relaxation measurement saturates at machine precision), not the
code's.

| id | claim | measured | verdict |
|---|---|---|---|
| C1 | E-M OU stationary variance `= 1/(1-dt/(2z))` (+14.3% at z=1, dt=0.25) | 1.1435 vs 1.142857 | pass |
| C1b | `ou_exact=True` gives variance exactly 1 | 1.00014 | pass |
| C2 | beta kernel `E=phi`, `Var=phi(1-phi)/(1+m0)` at v=0.02/0.2/0.5 | all within 1 MC sd | pass |
| C3 | gaussian kernel preserves the mean under truncation (phi=0.081, sd=0.30) | 0.081277 vs 0.081278 (naive truncation: 0.2699, bias +0.189) | pass |
| C4 | cpg(moment) variance `= phi(1-phi)/M`; matches beta only at M=m0+1 | 0.003513 (M=60) / 0.003452 (M=61) | pass |
| C5 | binomial partition: mean preserved, `sd=sqrt(alpha x(1-f)/(fV))`, molecules conserved | 0.4999 / 0.05024 / 0 violations | pass |
| C6 | G=2 product GRN ≡ Eq. (1) written out by hand | max abs err **0.0** (bitwise) | pass |
| C6b | same with asymmetric `s_1 = 0.55` (source-gene threshold) | max abs err **0.0** | pass |
| C7 | sharp limit p=40 converges to Table 1 cell-by-cell | beta/kappa/tau exact to 1e-6 | pass |
| C8 | paracrine: injected mass `= secretion·dt`; FTCS diffusion variance `2DT` | 0.25 exactly; `<x^2>/2DT = 0.98` | pass |
| C9 | overlap decays as `exp(-dt/tau)`, `tau=-0.25/ln(1-relax)` | log-slope −0.162519 vs −0.162519 | pass |
| C10 | `sister_corr` produces a positive sister correlation (rho=0 → +0.0019; rho=0.6 → +0.6007) | as claimed | pass |
| C11 | default `dt·lambda = 0.953` | 0.953461 | pass |

Additional consistency check: the continuous-vs-legacy classification
mismatch quoted in two places as "~8%" (rates.py) and "6.3/12.8/19.1%"
(model_spec.md) was measured as **5.91 / 12.57 / 19.43 %** of `[0,3]^G` for
G=2/3/4.  Both docs were wrong in different directions; the measurement is
now in both.

### 2.3 Axis 3 — doc ↔ code mismatches (11)

`docs/model_spec.md` was the weakest link: it describes **v1.0**, while the
code is v1.2.0.  Item-by-item:

1. **§1** omitted `a_tf` from the GRN equation — `a_tf` is the whole
   reprogramming axis of Figs. 11–12.
2. **§2a** claimed `act_i = rho_i + (1-rho_i) H+(x_i)`; `simulator.py:152`
   feeds `H+(x_i)` only.  A modelling choice, but an *undocumented* one.
3. **§2a** "timescale 1/(w0+w1+e0+e1) ~ 10-100 h" — with the defaults the
   relaxation rate is the constant `w0+w1+e0+e1 = 0.14 /h` → **7.1 h**.
4. **§2b** described `cpg` only as the mechanistic maintain/de-novo rule;
   the code default is `cpg_mode="moment"`, and `sister_corr` was absent.
5. **§3** "Var = x/(V f(1-f)) / 4" — this is symmetric in `f`, but the
   correct daughter-1 variance `alpha x (1-f)/(f V)` is not (the
   smaller-volume daughter gets the larger variance).  `alpha` was missing
   entirely.
6. **§4b** `kappa = kappa0 (1 + H_S(1-S))` — the code uses the
   **complementary** threshold `1-K`; the audit that fixed exactly this bug
   in v1.1.0 never updated the spec, so the spec still documented the bug.
7. **§6** "`dr = relax*overlap` per step" — superseded by the dt-consistent
   exponential relaxation in v1.1.0.
8. **§6** `field_interp` and `sphere_mask` undocumented.
9. **§8** "tests/test_cerps.py 37 unit/integration/cross-validation tests" —
   the file has 61 checks.
10. `docs/validation_report.md` §4 still reported the v1.1.0 suite
    ("11/11 tests pass, full suite 60/60").
11. `cerps/batch.py::__all__` omitted `load` and `dt_convergence_table`;
    `RatesConfig.delay_steps` carried a copy-pasted comment
    ("neighbours before inhibition starts" — it is a step count).

### 2.4 Axis 4 — README numbers vs archived data (2 wrong)

README (validation highlights):

> "homeostasis counts total 7026 (C++) vs 7159 (Python), SC 176.3 vs 176.8"

`reproduction/results/tables/cpp_vs_python.json`, m2 = 0.8, 30 runs:

| | C++ | Python |
|---|---|---|
| total | 7095.98 | 7158.53 |
| SC | 176.88 | 176.76 |

So the Python side was right (7159 / 176.8) and the C++ side was not
(**7096 / 176.9**).  `docs/validation_report.md` §2 had it right (7096), so
the README was the outlier.  The "< 0.3 pp transition-matrix agreement" claim
does hold (max 0.1369 pp at m2 = 0.8).

Everything else in the README's Table-3 block matches
`table3_reproduction.csv` to the last printed digit, and the Fig. 4
homeostasis numbers match `fig4_homeostasis.csv`.

### 2.5 Axis 5 — portability & hygiene (3 hard-coded paths)

| file | offending string | consequence |
|---|---|---|
| `cerps_sbi/__init__.py` | `C:/Users/黄荣生/Downloads/CERPS_complete_package/cerps_model/src` | the adapter imported on **exactly one machine**; everywhere else it failed at `from cerps import run_cerps` |
| `reproduction/scripts/validate_against_cpp.py` | `CPP_DIR = "/home/z/my-project/original_repo/Cell-type transition/output"` | a Linux path that belongs to a different machine, in a Windows-developed repo |
| `docs/paper-workflow/execution-runbook.md`, `cerps_sbi/README.md`, `docs/audit_2026-10-19.md` | 11 further absolute paths (interpreter, kernel, audit working dir) | documentation that cannot be followed by anyone else, plus the author's home directory in a **public** repository |

No secrets, no tokens, no files > 1 MB (repo total 4.3 MB, 46 Python files /
8 414 LOC, 25 markdown documents / 3 492 lines).

### 2.6 Axis 6 & 7 — packaging and CI

* There was **no `pyproject.toml` / `setup.py`**, so `pip install -e .` was
  impossible and `import cerps` only worked from the repository root.
* There was **no CI**, and `tests/test_m2_generalisation.py` defines **zero**
  `test_*` functions.  Consequence: the entire M2 proposal layer
  (`h_general.py` 660 LOC + `notation.py` 364 LOC + its 741-LOC test) would
  never have been executed by any automated run.  It does pass — it was just
  invisible.

### 2.7 One caveat the repository already owns, restated

`requirements.txt` pins `numba==0.65.0`; this machine has 0.68.0 and the
whole suite passes, including the Table-3 comparisons.  `validation_report.md`
§6 already documents that Table 3 is **statistically** but not
**bit-level** reproducible across numba versions (stored SC→SC 73.9885 vs
re-run 73.9016, both within MC error).  Nothing in this pass contradicts
that; the pin is a reproducibility statement, not a correctness requirement,
and it was left untouched.

---

## 3. Fixes applied (v1.2.1)

| axis | fix | file(s) |
|---|---|---|
| 2.4 | 7026→7096, 176.3→176.9 | `README.md` |
| 2.3 | counts 11/11, 60/60 → 61/93/19/18 = 191 + M2 V1–V7 | `docs/validation_report.md` |
| 2.3 | spec rewritten to v1.2.0: `a_tf`, `cpg_mode`/`sister_corr`, `alpha` + corrected partitioning variance, complementary `1-K`, `relax_time`, `field_interp`, `sphere_mask`, P6/P8, corrected timescale, measured boundary mismatch | `docs/model_spec.md` |
| 2.2 | "~8%" → measured 5.9/12.6/19.4% | `cerps/rates.py` |
| 2.5 | env-var/relative kernel resolution + explicit `ImportError`; `CELLTYPE_CPP_DIR` env var; 11 doc paths de-personalised | `cerps_sbi/__init__.py`, `reproduction/scripts/validate_against_cpp.py`, 3 docs |
| 2.6 | packaging metadata | `pyproject.toml` |
| 2.7 | CI running pytest **and** all five standalone harnesses | `.github/workflows/tests.yml` |
| 2.6 | 7 pytest wrappers for the M2 layer (pytest 62 → 69) | `tests/test_m2_generalisation.py` |
| 2.3 | `__all__` completed, `delay_steps` comment corrected | `cerps/batch.py`, `cerps/config.py` |
| 8 | single-branch policy documented and enforced | `CONTRIBUTING.md` + GitHub ruleset |

Nothing in `cerps/` that participates in the validated Table-3 baseline was
changed numerically: the only code edits outside the M2/packaging/portability
surface are a docstring and `__all__`.

---

## 4. Branch topology and the merge

`main@02db9d653b` and `add-paper-workflow@a9a16e1420`:

* compare endpoint: `ahead 0, behind 6, files changed 0`;
* file sets: `main` is a strict superset (5 extra files, all from the
  2026-10-21 hardening line);
* content hash comparison of all 127 shared paths: **zero differences**.

So merging `add-paper-workflow` into `main` is a no-op in content: the
branch is a stale pointer that adds no code, no tests and no data.  It was
therefore **deleted rather than merged** (an empty merge commit would only
misrepresent that the two trees were identical), and the branch that carried
this verification work was committed directly on `main`.

A repository ruleset now enforces the one-branch policy:

| rule | target | effect |
|---|---|---|
| creation | `~ALL` | no new branch of any name |
| deletion | `~DEFAULT_BRANCH` | `main` cannot be deleted |
| non-fast-forward | `~DEFAULT_BRANCH` | `main` history cannot be force-rewritten |

Contributors fork and open pull requests; see `CONTRIBUTING.md`.

---

## 5. What was *not* verified, and why

* **The paper's own Table 2/3 values.**  The "paper" column is transcribed
  from the published article; it cannot be re-derived from this repository.
* **The C++ ground truth.**  `cpp_vs_python.json` is an archived comparison;
  the C++ sources are not redistributed here, so it was not re-run.
* **Figure regeneration.**  `reproduction/scripts/run_batches.py` takes ~4 h
  on 2 cores; the archived figures/tables were read, not regenerated.
* **`cerps_sbi` end-to-end.**  It needs the external `CERPS_complete_package`
  kernel, which is not in this repository.  Its import path was fixed and
  verified to fail loudly; its sweeps were not re-run.
