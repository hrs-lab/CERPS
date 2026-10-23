# Contributing to CERPS

## Branch policy: `main` only

This repository is meant to have **exactly one branch, `main`**.

**Current state (2026-10-23).** The repository *does* have a single branch:
the stale `add-paper-workflow` branch (0 commits ahead / 6 behind, byte-identical
tree) was deleted, and `main` is the only ref.  The *enforcement* ruleset could
**not** be created by the automation that prepared this change, because
creating a ruleset requires `Administration: write`, which the available token
does not have.  It therefore has to be applied **once, by the repository
owner**:

> Settings → Rules → Rulesets → *New ruleset* → *Import as JSON*, then paste
> the block below (target **Branch**), and leave enforcement **Active**.

```json
{
  "name": "main-only-branch-policy",
  "target": "branch",
  "enforcement": "active",
  "conditions": { "ref_name": { "include": ["~ALL"], "exclude": [] } },
  "rules": [
    { "type": "creation" },
    { "type": "deletion" },
    { "type": "non_fast_forward" }
  ],
  "bypass_actors": [
    { "actor_id": 5, "actor_type": "RepositoryRole", "bypass_mode": "always" }
  ]
}
```

What that buys:

| rule | effect |
|---|---|
| `creation` | no new branch of any name can be created |
| `deletion` | `main` cannot be deleted |
| `non_fast_forward` | `main` history cannot be force-rewritten |
| bypass actor (admin) | the repository owner can still lift it when needed |

Until it is applied the policy is a convention only, so please respect it:
to contribute, **fork the repository, work on your fork, and open a pull
request against `main`**.  There is deliberately no `develop`/`feature-*`
branch; work that is not ready lives in a fork or in an open PR, where it is
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
