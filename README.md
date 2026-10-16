# CERPS — Coupled Epigenetic, Regulatory, Population and Spatial dynamics of cell-type transitions

This repository accompanies the manuscript:

> Rongsheng Huang (2026). *Mitotic epigenetic erosion and a regeneration–memory
> tradeoff in coupled epigenetic, gene-regulatory and population dynamics of
> stem-cell systems.* (in preparation, PLOS Computational Biology)

It contains:

1. **`reproduction_original/`** — a faithful Python (numba) reproduction of
   Huang, Situ & Lei (2024), *Dynamics of cell-type transition mediated by
   epigenetic modifications*, J. Theor. Biol. 577:111664
   (all main-text results: Figs 2–13 and Table 3).
2. **`cerps/`** — the CERPS model: continuous intracellular epigenetic
   dynamics, a lineage-dependent mitotic inheritance kernel, multi-gene GRNs,
   phenotype-free continuous kinetics, bidirectional coupling and 3-D space.
3. **`experiments/`** — chunked, resumable experiment runners (task queues).
4. **`data_scripts/` + `inference/`** — public-data processing (Paul 2015
   scRNA-seq, Beerman 2013 RRBS) and SMC-ABC posterior inference.
5. **`analysis/`** — publication figures.
6. **`docs/`** — full model specifications and protocols.

See `docs/SPEC_MODELS.md` for exact equations, parameters and validation
targets, and `docs/PLAN.md` for the experiment plan.

## Quick start

```bash
python reproduction_original/orig_model.py 1 2 3   # Table-3 validation runs
python cerps/kernel.py                            # CERPS smoke test
```

Requirements: `numpy`, `numba`, `scipy`, `matplotlib` (and `tectonic` to build
the paper).

## License

MIT (see LICENSE). Simulation code by Rongsheng Huang, Jimei University.
