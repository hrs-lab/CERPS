# CERPS

**C**oupled **E**pigenetic–**R**egulatory–**P**opulation–**S**patial model —
a general framework for simulating cell differentiation, reprogramming and
tissue growth, built as a faithful, fully validated extension of the
hybrid model of

> Rongsheng Huang, Qiaojun Situ, Jinzhi Lei (2024).
> *Dynamics of cell-type transition mediated by epigenetic modifications.*
> Journal of Theoretical Biology **577**, 111664.
> https://doi.org/10.1016/j.jtbi.2023.111664

This repository contains two parts:

1. **`reproduction/`** — a complete, independent **Python re-implementation
   of every result of the paper** (Figs. 2–13 and Table 3), validated
   against both the published numbers and the authors' original C++ code.
2. **`cerps/`** — the **CERPS framework**, which generalises the model
   along six axes (details in `docs/model_spec.md`):

   | # | requirement | CERPS implementation |
   |---|---|---|
   | 1 | intracellular epigenetic dynamics; inheritance not limited to the beta distribution | continuous mark writing/erasing coupled to gene activity + four division kernels (`beta`, `gaussian`, `cpg`, `deterministic`) |
   | 2 | gene-expression changes during cell division | exact (paper) or molecule-level **binomial partitioning** |
   | 3 | multi-gene regulatory networks | N-gene GRN with signed regulation edges (product-form cross-repression, which reduces exactly to Eq. (1) of the paper for 2 genes) |
   | 4 | continuous, state-dependent proliferation/death rates (no discrete classification) | smooth Hill-function kinetics of a continuous *stemness* measure; sharp-limit == Table 1 of the paper |
   | 5 | bidirectional coupling of epigenetics, GRN and population dynamics | epi↔GRN (rates vs. activity-dependent writing), cell↔population (niche feedback, local density, paracrine field) |
   | 6 | 3D spatial position | off-lattice 3D mechanics, contact inhibition, paracrine diffusion field |

---

## Validation highlights (reproduction vs. paper)

**Table 3 — mother→daughter transition probabilities (%)**
(statistics window t ∈ [2000, 4000] h, 30 runs, 15.6 M recorded transitions):

| | SC | TA1 | TA2 |
|---|---|---|---|
| SC (paper) | 73.9402 | 12.9889 | 13.0709 |
| SC (this repo) | **73.9885** | **12.9919** | **13.0197** |
| TA1 (paper) | 0.0013 | 99.9955 | 0.0032 |
| TA1 (this repo) | **0.0010** | **99.9955** | **0.0035** |
| TA2 (paper) | 0.0013 | 0.0035 | 99.9952 |
| TA2 (this repo) | **0.0015** | **0.0032** | **99.9953** |

Cross-validation against the authors' original C++ code (m2 = 0.8 column,
30 runs): homeostasis counts total 7096 (C++) vs 7159 (Python), SC
176.9 vs 176.8; transition matrices agree to < 0.3 percentage points.
Two independent Python engines (numba + NumPy reference) agree within
Monte-Carlo error.  See `docs/validation_report.md` for the full report
and `reproduction/results/` for all figures and tables.

## Quickstart

```bash
pip install -r requirements.txt

# --- CERPS framework -------------------------------------------------
python examples/example_legacy_reproduction.py   # = the paper's model
python examples/example_epigenetic_kernels.py    # 4 inheritance kernels
python examples/example_multigene.py             # 3-gene GRN
python examples/example_3d_spatial.py            # 3D spheroid + paracrine
python tests/test_cerps.py                       # 61 tests (all pass)
python tests/test_audit_fixes.py                 # 93 audit-fix tests
python tests/test_audit_2026_10_21.py            # 19 tests (2026-10-21 audit)
python tests/test_audit_ranges.py                # 18 range-scan tests (audit §8)
python tests/test_m2_generalisation.py           # V1-V7 proposal-layer checks
# every file above runs standalone; python -m pytest tests/ collects the
# 69 pytest-style cases (the two harness files expose 13/12 of their
# 61/93 checks as test_* functions; full 191 via the standalone commands,
# plus the 7 M2 proposal-layer checks)
#
# `pip install -e ".[test,repro]"` then `python -m pytest tests/ -q`
# is the CI-equivalent run (the workflow file is in docs/verification_2026-10-23.md §6).

# --- batch experiments (the experiment layer) -------------------------
python -c "
from cerps.config import CERPSConfig, RatesConfig
from cerps import batch
base = CERPSConfig(max_time=2000.0, rec_tran_from=1000.0)
res = batch.sweep(base, {'rates.beta0': [0.02, 0.04],
                         'epi.kernel': ['beta', 'cpg']},
                  seeds=(0, 1, 2), out_dir='runs')
batch.report(batch.aggregate(res))
print(batch.health_check(base, t_short=500.0))
"

# --- paper reproduction ----------------------------------------------
cd reproduction/scripts
python fig2_bifurcation.py          # Figs. 2A-D
python run_batches.py               # all simulation batches (~4 h, 2 cores)
python table3.py                    # Table 3 + paper comparison
python fig_sweeps.py                # Figs. 8, 9, 12
python fig3_noise.py fig45_population.py fig6_trajectories.py \
        fig7_extinction.py fig10_asymmetry.py fig11_reprogramming.py \
        fig13_landscape.py          # remaining figures
python validate_against_cpp.py      # compare with the original C++ output
```

Each batch checkpoints its results under `reproduction/results/data/`
(skipped on re-run); figures land in `reproduction/results/figures/`,
tables in `reproduction/results/tables/`.

### The batch layer (`cerps/batch.py`)

`Simulator` runs ONE population for ONE configuration.  `cerps.batch`
adds what mechanism studies need:

| function | purpose |
|---|---|
| `sweep(base, grid, seeds, out_dir)` | Cartesian parameter sweep with checkpoint/resume (`results.jsonl`) |
| `run_one(cfg, seed)` | one run + homeostasis counts + transition matrix + runtime |
| `aggregate(results)` | mean/std across seeds per configuration (the MC error) |
| `homeostasis(records, t_from, t_to)` | time-averaged counts over a window |
| `transitions(sim)` | row-normalised transition probabilities, label-coded (0=SC, 1=TA1, 3=TA2) |
| `health_check(cfg)` | paracrine CFL number + dt-halving convergence estimate |
| `validate(cfg)` | warnings for configurations that silently invalidate a comparison |

`health_check` is the one to run before trusting any mechanism effect: it
reports the relative change of the homeostasis counts when `dt` is
halved, so a "mechanism effect" can be told apart from a discretisation
effect.  `validate` catches the traps that make a comparison meaningless
(e.g. sweeping `division.alpha` while `partition="exact"` ignores it).

## Repository layout

```
cerps/                  the CERPS package (config, grn, epigenetics,
                        division, rates, spatial, simulator, batch,
                        h_general, notation)
reproduction/
  hsl2024/              python reproduction of the paper
    numba_sim.py        fast engine (mirrors the C++ protocol exactly)
    reference_sim.py    independent NumPy engine (cross-validation, Fig. 6)
    steady_state.py     bifurcation / attractor analysis (Fig. 2A/B)
    single_cell.py      single-cell dynamics (Fig. 2C/D)
    analysis.py         statistics (mirrors DataAnalysis.m)
  scripts/              one runner per figure + batch orchestration
  results/              generated figures, tables, batch data
tests/                  test suites (test_cerps.py 61, test_audit_fixes.py
                        93, test_audit_2026_10_21.py 19,
                        test_audit_ranges.py 18, test_m2_generalisation.py)
examples/               runnable feature demonstrations
docs/
  model_spec.md         full mathematical specification
  literature_review.md  survey motivating every design choice
  protocol_notes.md     paper-vs-code protocol differences
  validation_report.md  reproduction validation report
  verification_2026-10-23.md  multi-axis cross-verification of this repo
  audit_2026-10-19.md   log of the fixes applied after the code audit
  audit_2026-10-21.md   deep audit: degenerate-limit exactness, new switches
  validity_matrix.md    mechanism x regime x switches x bias lookup table
  paper-workflow/       methodology set for turning CERPS into a paper
cerps_sbi/              adapter over the external CERPS_complete_package
                        kernel (needs $CERPS_SRC / $CERPS_PKG; see its README)
.github/workflows/      CI definition (see the note below -- must be committed
                        once by the repository owner)
CONTRIBUTING.md         branch policy (main only) and contribution rules
pyproject.toml          packaging metadata (`pip install -e ".[test,repro]"`)
```

> **Note on the CI file.**  `.github/workflows/tests.yml` exists in the
> working tree of this change but is **not yet in the repository**: GitHub
> refuses to let a personal access token write under `.github/`
> (workflow-injection protection), so it has to be added once by the
> repository owner (web UI *Add file*, or a normal `git push`).  The full
> content is reproduced in `docs/verification_2026-10-23.md` §6, and the
> equivalent local commands are in `CONTRIBUTING.md`.

## Changelog

### v1.2.1 (2026-10-23) -- multi-axis cross-verification + single-branch policy

An independent verification pass (`docs/verification_2026-10-23.md`) ran
eight axes: test-suite execution, independent numerical recomputation of the
analytic claims, doc/code consistency, portability/secret scan, packaging,
CI coverage, and branch topology.  Every documented number was recomputed;
the code's own claims held, but the *documentation and packaging* had drifted:

| id | finding | fix |
|---|---|---|
| V1 | README quoted C++ totals "7026 / 176.3"; the archived `cpp_vs_python.json` says **7096 / 176.9** (as does `validation_report.md` §2) | README corrected |
| V2 | `docs/validation_report.md` §4 still reported the v1.1.0 suite (11/11, 60/60) | updated to 61/93/19/18 = 191 + the M2 V1-V7 suite |
| V3 | `docs/model_spec.md` described v1.0 only: no `a_tf`, no `cpg_mode="moment"` (the default), no `sister_corr`, no `alpha` in the partitioning variance, `kappa` with the *same* threshold `K` instead of the complementary `1-K`, `dr = relax*overlap` instead of the dt-consistent relaxation, `tests/test_cerps.py 37 tests` | spec brought to v1.2.0, with the corrected complementary threshold and the measured boundary mismatch (5.9/12.6/19.4% for G=2/3/4) |
| V4 | `cerps_sbi/__init__.py` hard-coded `C:/Users/<name>/Downloads/CERPS_complete_package/...`, so the adapter imported on exactly one machine; `validate_against_cpp.py` hard-coded `/home/z/my-project/...` | env-var + relative resolution only, with an explicit `ImportError` explaining how to point at the kernel |
| V5 | `tests/test_m2_generalisation.py` defines **no** `test_*` functions, so `pytest tests/` never ran the M2 layer and CI would have missed it entirely | 7 pytest wrappers added; `pytest tests/` now collects 69 cases |
| V6 | no packaging metadata (`pip install -e .` impossible) and no CI | `pyproject.toml` added; the CI workflow is written out in `docs/verification_2026-10-23.md` §6 but cannot be committed by a PAT (§4b) |
| V7 | `cerps.batch.__all__` omitted `load` and `dt_convergence_table`; `RatesConfig.delay_steps` carried a copy-pasted comment ("neighbours before inhibition starts") | both fixed |
| V8 | the repository carried two branches (`main`, `add-paper-workflow`, 0 ahead / 6 behind, identical tree) | stale branch deleted, `main` is now the only branch; the enforcing ruleset needs admin rights and is documented as a copy-paste block in `CONTRIBUTING.md` |

The verification also **confirmed** the numbers that had already drifted into
legend: 27/27 independent recomputations of the analytic claims passed (OU
variance inflation `1/(1-dt/(2z))`, beta/gaussian/cpg kernel moments, binomial
partitioning mean and variance, the exact G=2 reduction to Eq. (1) for both
symmetric and asymmetric thresholds, the sharp-limit convergence to Table 1,
paracrine mass conservation, the CFL sub-stepping, the exponential overlap
relaxation, and `dt*lambda = 0.953`).  Test counts matched exactly:
61 + 93 + 19 + 18 = 191 standalone, 69 collected by pytest.

### v1.2.0 (2026-10-22) -- hardening fixes merged + M2 proposal layer

Merged the 2026-10-20 hardening line (five measured defects, each with
regression tests) and the M2 generalisation layer into `main`:

| id | defect | measured before | after |
|---|---|---|---|
| D1 | `_record` hard-coded the label columns `{0,1,3}`, so cells with label >= 4 (G >= 4) were silently dropped from every `ho_*` summary | G=4 lost 20.3%, G=5 lost 30.8% | labels >= 4 counted in `other` |
| D2 | P6 memory EWMA gain `dt/tau` unclamped: `tau < dt` diverged | `tau=0.1, dt=0.25`: max\|mem\| = 1.8e42 | clamped, max\|mem\| = 0.967 |
| D3 | `MemoryConfig.driver` was a dead option (never read) | `u_mean` ran bit-identical to `stemness` | driver honoured |
| D4 | `run()` re-entry concatenated records from two time origins | run(100)+run(200) = 1200 records | 800 records, series ends at t=200; P8 delay buffer also reset |
| M2a | `GRNConfig.matrices()` accepted wrong-length flat `b_rep`/`c_act` | crashed later inside `GRN.step` (IndexError) | `ValueError` at `matrices()` |

New **M2 generalisation layer** (`cerps/h_general.py`, `cerps/notation.py`):
the RACIPE-style proposal form `F(x) = (G(x), H(x))` with epigenetic gates
(`alpha`/`beta`/`k`), verified V1-V7 in `tests/test_m2_generalisation.py` --
including bit-exact agreement with the native GRN wherever the reduction is
defined.  `ProposalGRN.from_native` now *refuses* product mode with
`b_base != 1` (which it cannot express) instead of silently dropping
`b_vec`; additive mode stays exact for any `b_base`.

Suite after the merge: every test file is standalone-runnable
(`python tests/<file>.py`) -- `test_cerps.py` 61 + `test_audit_fixes.py`
93 + `test_audit_2026_10_21.py` 19 + `test_audit_ranges.py` 18 = **191
checks passed**, plus the M2 V1-V7 suite.  (`python -m pytest tests/`
collects the 62 pytest-style cases; the two harness files expose only
13/12 of their checks as `test_*` functions, so the standalone commands
are the authoritative full run.)

### v1.1.1 (2026-10-21) -- deep audit: exactness in degenerate limits

An audit of *where the model is exact in degenerate limits but may fail in
general* (`docs/audit_2026-10-21.md`) quantified four biases and **added four
switches, all OFF by default so the validated Table-3 baseline is untouched**.

| switch | default | removes | quantified bias |
|---|---|---|---|
| `GRNConfig.ou_exact` | `False` | Euler-Maruyama OU variance inflation `1/(1-dt/2ζ)` | ζ=1,dt=0.25: +14.3%; ζ=0.25: **+100%** |
| `SpaceConfig.field_interp` | `"nearest"` | O(h) nearest-grid sampling error | mean error 27.3% → **2.0%** (`"linear"`) |
| `SpaceConfig.sphere_mask` | `False` | field leaking into cell-free cube corners | geometric; matters when √(D/k) ≳ 0.5R |
| `EpiConfig.sister_corr` | `0.0` | zero sister-daughter correlation | corr(d1,d2) ≈ −0.006 → ρ-controlled |

Also added: `batch.grn_stiffness` (per-gene `dt·λ_i`, default 0.953),
`batch.dt_convergence_table` (dt/dt/2/dt/4 observed-order table, audit A2),
`batch.validate` warnings, `docs/validity_matrix.md`,
`tests/test_audit_2026_10_21.py` (19 regression tests) and
`tests/test_audit_ranges.py` (18 range-scan tests pinning audit §8). Full
suite at this release: **60 passed** (11 + 12 + 19 + 18; the hardening-line
merge later expanded the first two files to 61 + 93 -- see v1.2.0). See
`docs/audit_2026-10-21.md` §8 for the parameter-range and
symmetry/asymmetry scan.

### v1.1.0 (2026-10-19) -- audit fixes

An independent audit of the framework (literature review + 10 numerical
cross-checks) found four defects that only appear outside the symmetric
parameter regime and the attractor interiors.  All are fixed; the
regression tests live in `tests/test_audit_fixes.py`.

1. **`simulator.py`** -- the phenotype label was re-evaluated for *all*
   cells every step.  The paper's numerical scheme, the C++ code and
   `reference_sim.py` all freeze it outside the resting phase.  The old
   behaviour lowered the SC count `Q` during the cycle and inflated SC
   self-renewal by **+2.4 pp** (76.7% vs 74.2%); with the fix CERPS
   agrees with the `hsl2024` engine to **0.006 pp**.
2. **`grn.py`** -- cross-regulation thresholds used `(s_i+s_j)/2` instead
   of the source gene's `s_j`, so the exact reduction to Eq. (1) of the
   paper failed whenever `s_1 != s_2` (e.g. the `s_1 = 0.55` series of
   Fig. 10).  The reduction is now exact for every tested parameter set.
3. **`spatial.py`** -- the paracrine solver clamped the diffusion number
   to the 3-D CFL limit `1/6`, silently replacing the nominal
   diffusivity (D=25, dt=0.25, h=4 ran at D_eff = 10.67).  It now
   sub-steps instead, and uses reflecting (zero-flux) boundaries that
   match the spherical cell domain instead of `np.roll` periodicity
   (which leaked ~8% of the signal across the domain).
4. **`rates.py`** -- the TA-side Hill switch used the same threshold `K`
   as the SC-side one, so the two fired at `S = K` and `S = 1-K`.  This
   is harmless only at `K = 0.5`; with the (now default) auto-derived
   threshold it gave stem-like cells the TA removal rate.  The TA side
   now uses the complementary threshold `1-K`.

Orthogonality and interface changes:

* **Kernels are now controlled noise-only variants.**  `cpg` defaults to
  `cpg_mode="moment"`, which sets the per-site marking probability to
  `phi(v)` so that the mean *and* (with `cpg_sites = m0+1`) the variance
  match the beta kernel.  The classical maintain/de-novo rule remains
  available as `cpg_mode="mechanistic"` (documented: its mean map is
  affine, so it has a different fixed point).
* **`gaussian` truncation now preserves the mean.**  Truncating a normal
  shifts its mean; the kernel solves for the pre-truncation mean that
  makes the truncated mean exactly `phi(v)` (the old folding gave a
  +0.167 bias at `sigma_g = 0.3, v = 0.02`).
* **`division.alpha`** exposes the continuous partitioning-noise
  strength of Soltani et al. 2016 (Eq. 5): `alpha = 0` deterministic,
  `1` independent binomial, `>1` clustered, `<1` actively ordered.
* **`K_switch=None` (default)** derives the stemness threshold so the
  continuous classification boundary coincides with the legacy
  `x_threshold` on the diagonal.
* **`CERPSConfig.rec_tran_from`** discards the transient from the
  transition statistics (the paper counts only `t > 2000 h`), so
  `sim.transitions` is directly comparable with Table 3.
* **`Space.relax`** is now time-step consistent (exponential overlap
  decay with `relax_time`); `relax` is kept as the legacy per-step
  fraction and is reproduced exactly at `dt = 0.25 h`.
* **New `cerps/batch.py`** experiment layer (see above).

Documentation fixes: `epigenetics.py` cited Dodd 2007 / Probst 2009 /
Zhao 2021 for the activity-dependent writing ODE; none of the three
supports it (Dodd 2007 has no such ODE and no beta kernel, Probst 2009 is
a qualitative review, Zhao 2021 explicitly excludes transcription).  The
correct precedents are Berry, Dean & Howard 2017 and Miyamoto, Furusawa
& Kaneko 2015, now cited.  `params.py` / `protocol_notes.md` now record
that the paper's Table 2 lists `m2 = 1.0` while its own simulations
(Sec. 3.4.2, 3.5, Fig. 12, and the released C++ sweep) use `m2 = 0.8`.

## Upstream code

The authors' original C++/Matlab implementation (GPL-3.0) is at
https://github.com/jinzhilei/Cell-type-transition and is **not**
redistributed here.  `reproduction/docs` documents exactly how the Python
engines mirror its protocol, including the places where the released code
differs from the text of the paper (see `docs/protocol_notes.md`).

## License and credit

GPL-3.0 (this project derives from the model and protocol of the
GPL-3.0-licensed reference implementation).  If you use CERPS or the
reproduction in academic work, please cite the original paper:

```bibtex
@article{huang2024dynamics,
  title={Dynamics of cell-type transition mediated by epigenetic
         modifications},
  author={Huang, Rongsheng and Situ, Qiaojun and Lei, Jinzhi},
  journal={Journal of Theoretical Biology},
  volume={577},
  pages={111664},
  year={2024},
  doi={10.1016/j.jtbi.2023.111664}
}
```

---

## 中文说明（快速上手）

本仓库包含两部分：

1. **`reproduction/`**：用 Python 完整复现论文（Huang, Situ & Lei, 2024）
   的全部结果（图 2–13 与表 3），并与作者发布的 C++ 原始代码进行了
   交叉验证（转移概率矩阵、稳态细胞数均一致；表 3 九项指标全部
   对齐到小数点后 3–4 位）。
2. **`cerps/`**：CERPS 通用框架，在原模型基础上实现六项扩展：
   ① 细胞内表观遗传连续动力学 + beta/高斯/CpG 位点/确定性四种分裂
   遗传内核；② 分裂时基因表达的二项分配；③ 多基因调控网络
   （乘积型交叉抑制，双基因时精确退化为论文式 (1)）；④ 增殖/死亡
   速率连续依赖细胞状态（陡峭 Hill 极限下严格回归论文表 1）；
   ⑤ 表观遗传–基因调控–细胞群体三者的双向耦合（含旁分泌信号场）；
   ⑥ 三维空间（非格子力学、接触抑制、扩散场）。

运行方式见上方 Quickstart；`docs/` 内含完整数学规范、文献综述、
协议差异说明与验证报告。许可证为 GPL-3.0（继承自原始代码）。
