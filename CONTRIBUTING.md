# Contributing to CERPS

## Branch policy: `main` only

This repository has **exactly one branch, `main`**, and that is enforced by a
GitHub ruleset, not just by convention:

* creating a new branch is blocked,
* deleting `main` is blocked,
* force-pushing `main` is blocked,
* a pull request must be open (or the commit must be pushed by a
  collaborator) before history is rewritten in any other way.

To contribute, **fork the repository, make your change on your fork, and
open a pull request against `main`**.  Direct pushes to `main` are reserved
for the repository owner.

There is deliberately no `develop`/`dev`/`feature-*` branch.  Work that is
not ready for `main` lives in a fork or in an open pull request, where it is
visible and reviewable.

## Development

```bash
pip install -e ".[test,repro]"
python -m pytest tests/ -q          # 69 pytest cases
python tests/test_cerps.py           # 61 checks (human-readable report)
python tests/test_audit_fixes.py     # 93 checks
python tests/test_audit_2026_10_21.py
python tests/test_audit_ranges.py
python tests/test_m2_generalisation.py   # V1-V7, proposal layer
```

The two hand-written harness files (`test_cerps.py`,
`test_audit_fixes.py`) expose only part of their checks as `test_*`
functions, so **`pytest` and the standalone runs count different things**;
the standalone commands are the authoritative full run.

## Rules for a change

1. **Run the suite.**  A change that alters numerics must keep
   `python -m pytest tests/` green and must not silently move the Table-3
   baseline (`rec_tran_from=2000`, `numba==0.65.0`).
2. **Update the spec.**  `docs/model_spec.md` is the single source of truth
   for the mathematics.  If you add or change a switch, its default must be
   documented there, and the default must leave the validated baseline
   untouched.
3. **Keep the registry honest.**  New mechanism ideas go into
   `docs/paper-workflow/candidate_register.md` with an explicit status; a
   result that cannot be reproduced goes to `GRAVE` in
   `mechanism-lab/hypotheses/_ledger.tsv`, not into the changelog.
4. **No machine-specific paths.**  External kernels are located through
   environment variables (`CERPS_SRC`, `CERPS_PKG`) or relative paths --
   never an absolute path in your home directory.
5. **Record the environment** whenever a number is (re)generated:
   `python -c "import numpy, scipy, numba; print(...)"` plus the OS.  The
   archived Table-3 data are statistically reproducible but not
   bit-reproducible across numba versions (see
   `docs/validation_report.md` §6).

## License

GPL-3.0-or-later.  The model and protocol derive from the GPL-3.0 reference
implementation at <https://github.com/jinzhilei/Cell-type-transition>; please
cite Huang, Situ & Lei (2024), *J. Theor. Biol.* 577:111664.
