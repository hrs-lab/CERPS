# Project Plan, Milestones & Session Continuity

## CRITICAL: session-restore procedure (each new session, ~2 min)
1. `cd /home/z/my-project && git clone https://hrs-lab:$(cat .ghtoken)@github.com/hrs-lab/CERPS.git cerps_repo` (or `git -C cerps_repo pull` if dir exists)
2. Read `cerps_repo/docs/SPEC_MODELS.md`, `cerps_repo/docs/PLAN.md`, `cerps_repo/worklog.md`
3. Check `cerps_repo/results/*/status` (task queue JSONs) -> continue pending tasks
4. NEVER leave uncommitted work: commit+push after EVERY completed unit
5. Checkpoint zip of workspace to /home/z/my-project/download/ at session end

## Deliverables (user requirements)
1. Full Python reproduction of Huang/Situ/Lei 2024 (all main results Fig 2-13 + Table 3)
2. CERPS extension (5 modules) + mechanism (mitotic erosion -> regenerative debt)
3. PLOS Comp Biol manuscript (LaTeX + figures + PDF + SI) in ONE zip in download/
4. Code on github.com/hrs-lab/CERPS; commit author:
   hrs-lab <77771845+hrs-lab@users.noreply.github.com>
5. Paper: Bayesian posterior (SMC-ABC), sensitivity (LHS-PRCC), robustness,
   2-dataset cross-validation (Paul 2015 GSE72857 + Beerman 2013 GSE44117)
6. Author: Rongsheng Huang, School of Science, Jimei University, Xiamen 361021,
   hrs@jmu.edu.cn. Acknowledge NSFC Grant No. 12401641. Code availability:
   https://github.com/hrs-lab/CERPS

## Sandbox constraints (CRITICAL)
- ~2h session lifetime; ALL in-sandbox files lost at reset -> push early, push often
- Background processes are reaped at tool-call boundaries -> NO nohup/setsid;
  run everything in foreground chunks: each Bash call <= ~9 min (timeout 565s),
  chunked task-queue pattern: tasks JSON + one .npz per task, resumable
- 2 cores, 3 GB RAM (avoid loading full UMI matrix: single-pass filtered parsing)
- tectonic for LaTeX; numba @njit(cache=True) everywhere

## Milestones (each ~1 session)
- [S1] Infra + docs + GitHub push + orig_model.py rebuild + Table 3 validation
       + CERPS kernel rebuild + smoke test
- [S2] Orig-model experiments: fig2 (bifurcation+duration), fig3, fig4, fig5, table3(30),
       fig7, fig8, fig9, fig10, fig11, fig12, fig13 aggregation -> all orig results
- [S3] CERPS experiments: M1 (12x4), M2 harsh (8x12), M4 dose (6x12), M3, M5,
       SEN (LHS 160), ROB (9x8) -> results npz all pushed
- [S4] Data: download Paul umitab (16MB) + Beerman 12 files (chunked curl -C -),
       process_paul.py (single-pass, 20 target genes), process_beerman.py
       (methylation vs proliferation history); SMC-ABC inference + posterior figs
- [S5] All publication figures (matplotlib, EN labels, 300dpi, PDF+PNG)
- [S6] PLOS LaTeX main + SI + references (verify every citation against paper txt
       + web), tectonic compile, final zip to download/, final GitHub push
- Buffer session for calibration retries / re-runs

## File layout
```
cerps_repo/
  README.md  LICENSE  .gitignore
  docs/SPEC_MODELS.md  docs/PLAN.md  worklog.md
  reproduction_original/orig_model.py (+ analysis scripts)
  cerps/kernel.py  cerps/experiments.py
  experiments/runner.py (chunked task queue)
  data_scripts/process_paul.py  process_beerman.py  download_queue.sh
  inference/smcabc.py
  analysis/figures_*.py
  results/{reproduction,cerps}/jobs/*.npz + status JSONs
  paper/manuscript/ (main.tex, S1_Text.tex, references.bib, figures/)
```

## Paper outline (PLOS Comp Biol)
Title: "Mitotic epigenetic erosion drives a regeneration-versus-memory tradeoff in
coupled gene-regulatory, epigenetic and population dynamics of stem cell systems"
- Abstract (unstructured ~300w), Author Summary
- Intro: cell-type transition, epigenetic inheritance, original model recap, gap:
  inheritance not only at division / discrete classes / no feedback / no space
- Model: CERPS formulation (5 modules, equations above), legacy correspondence
- Results:
  1. Reproduction of original (Table 3 + Fig 4-5 anchor)
  2. Homeostasis in CERPS (steady states, balance score, spatial structure)
  3. Mitotic erosion & regenerative debt (M1 control/acute/serial; M4 dose-response)
  4. Coupling attribution (M2 ablations: q(d) lineage decay necessary)
  5. Aging & spatial modulators (M5, M3)
  6. Bayesian inference from 2 datasets + cross-validation (Paul + Beerman)
  7. Sensitivity (PRCC) & robustness
- Discussion: regenerative bankruptcy, biological predictions (5-FU history,
  transplant seriality), relation to clonal hematopoiesis literature
- S1..Sn SI Text: derivations (u**, S), extra figures, full PRCC tables,
  ABC details, data processing details

## Known pitfalls (from previous session)
- numba: no dynamic tuple indexing (snap_times as ndarray arg); NamedTuples -> flat tuples
- aggregation: per-run row-normalize then average across runs (MATLAB DataAnalysis.m)
- GEO downloads: URL must include &file=; 135070 B = error page; use curl -C -
- CERPS calibration traps: density feedback can mask stress (use ablation protocol);
  senescence must act only on low-S cells; w1 anchor must keep u* safely inside domain
- Hill form x^n/(s^n+x^n) — NOT x^n/(s+x^n)
