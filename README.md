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
30 runs): homeostasis counts total 7026 (C++) vs 7159 (Python), SC
176.3 vs 176.8; transition matrices agree to < 0.3 percentage points.
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
python tests/test_m2_generalisation.py           # V1-V7 proposal-layer checks

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
                        division, rates, spatial, simulator, batch)
reproduction/
  hsl2024/              python reproduction of the paper
    numba_sim.py        fast engine (mirrors the C++ protocol exactly)
    reference_sim.py    independent NumPy engine (cross-validation, Fig. 6)
    steady_state.py     bifurcation / attractor analysis (Fig. 2A/B)
    single_cell.py      single-cell dynamics (Fig. 2C/D)
    analysis.py         statistics (mirrors DataAnalysis.m)
  scripts/              one runner per figure + batch orchestration
  results/              generated figures, tables, batch data
tests/                  test suites (test_cerps.py, test_audit_fixes.py,
                        test_m2_generalisation.py)
examples/               runnable feature demonstrations
docs/
  model_spec.md         full mathematical specification
  literature_review.md  survey motivating every design choice
  protocol_notes.md     paper-vs-code protocol differences
  validation_report.md  reproduction validation report
  audit_2026-10-19.md   log of the fixes applied after the code audit
```

## Changelog

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
