# Reproduction of Huang, Situ & Lei (2024)

Python reproduction of *Dynamics of cell-type transition mediated by
epigenetic modifications*, J. Theor. Biol. **577** (2024) 111664,
doi:10.1016/j.jtbi.2023.111664, including every figure (2-13) and
Table 3.  The protocol mirrors the authors' released C++ code
(github.com/jinzhilei/Cell-type-transition, GPL-3.0); the differences
between the paper text and that code are documented in
`../docs/protocol_notes.md`.

## Package `hsl2024`

| module | contents |
|---|---|
| `params.py` | Table-2 parameters (`Params` dataclass) |
| `numba_sim.py` | fast engine; mirrors the C++ stage machine, fate decisions, epigenetic inheritance and recording windows exactly |
| `reference_sim.py` | independent NumPy engine (cross-validation; lineage tracking for Fig. 6) |
| `sim.py` | driver around the numba engine (batched over runs) |
| `steady_state.py` | multistart-Newton steady states, Jacobian stability (Figs. 2A/B) |
| `single_cell.py` | noise-driven single-cell dynamics (Figs. 2C/D) |
| `analysis.py` | statistics mirroring `DataAnalysis.m` (homeostasis windows, transition matrices, Waddington landscape) |
| `style.py` | paper-like plotting style |

## Running

```bash
python scripts/fig2_bifurcation.py       # Fig. 2 (steady-state analysis)
python scripts/run_batches.py            # all population batches (~5 CPU-h)
python scripts/table3.py                 # Table 3 vs. paper
python scripts/fig_sweeps.py             # Figs. 8, 9, 12
python scripts/fig3_noise.py             # Fig. 3
python scripts/fig45_population.py       # Figs. 4, 5
python scripts/fig6_trajectories.py      # Fig. 6 (reference engine)
python scripts/fig7_extinction.py        # Fig. 7
python scripts/fig10_asymmetry.py        # Fig. 10
python scripts/fig11_reprogramming.py    # Fig. 11
python scripts/fig13_landscape.py        # Fig. 13
python scripts/validate_against_cpp.py   # compare with the original C++
```

Batches checkpoint their outputs in `results/data/<job>.npz` (skipped on
re-run, `--force` re-runs).  Figures land in `results/figures/`, tables in
`results/tables/`.  `run_batches.py --budget S` stops submitting new jobs
after S seconds (useful for short sessions).

The C++ reference is expected at `../../../original_repo/...` for the
validation script; see `docs/validation_report.md` for the comparison
results and `scripts` docstrings for the exact protocol of each figure.
