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
python tests/test_cerps.py                       # 37 tests (all pass)

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

## Repository layout

```
cerps/                  the CERPS package (config, grn, epigenetics,
                        division, rates, spatial, simulator)
reproduction/
  hsl2024/              python reproduction of the paper
    numba_sim.py        fast engine (mirrors the C++ protocol exactly)
    reference_sim.py    independent NumPy engine (cross-validation, Fig. 6)
    steady_state.py     bifurcation / attractor analysis (Fig. 2A/B)
    single_cell.py      single-cell dynamics (Fig. 2C/D)
    analysis.py         statistics (mirrors DataAnalysis.m)
  scripts/              one runner per figure + batch orchestration
  results/              generated figures, tables, batch data
tests/                  CERPS test suite (37 tests)
examples/               runnable feature demonstrations
docs/
  model_spec.md         full mathematical specification
  literature_review.md  survey motivating every design choice
  protocol_notes.md     paper-vs-code protocol differences
  validation_report.md  reproduction validation report
```

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
