# CERPS validity matrix

> Purpose: before running any mechanism, use this table to decide **whether the
> current regime needs a switch turned on** and **how much the default result
> may be off**. Companion documents: `docs/audit_2026-10-21.md` (evidence),
> `batch.validate(cfg)` / `batch.grn_stiffness(cfg)` (automatic warnings).
>
> In one sentence: **the default configuration is ready-to-use only for the
> paper protocol (G=2, legacy, zeta≈1, D≈25, dt=0.25, symmetric, beta kernel,
> independent sisters); outside that cell, the result depends on which switch
> you turned on.**

---

## 1. Default validity domain

| Dimension | Default / safe value | Consequence of leaving it |
|---|---|---|
| Number of genes G | 2 | G≥3 classification mismatch 12.8%/19.1% (§B3) |
| Rate mode | `legacy` | see M5 for `continuous` |
| OU timescale zeta | ≈1 (`ou_exact=False`) | zeta=0.25 noise variance +100% (§A1) |
| Diffusivity D | ≤10.7 (h=4, dt=0.25) | D=25: pre-audit engine ran at 42.7%; the new engine sub-steps (§C2) |
| Time step dt | 0.25 | oscillation when `dt·λ>1` (default already 0.95, §A3) |
| Parameter symmetry | symmetric | asymmetry makes one gene cross first (§A1/§A3) |
| Inheritance kernel | `beta` | cpg/gaussian: see §B1 |
| Sister correlation | independent (0) | lineage correlation suppressed (§B2, M8) |
| Field sampling | `nearest` | mean error 27% (§C2) |

---

## 2. Mechanism × regime × switches × bias

| Mechanism | Default usable? | Key regime | Required/recommended switches | Known bias (quantified) | Evidence |
|---|---|---|---|---|---|
| **Paper baseline** (G=2, legacy) | ✅ ready | zeta≈1, D≤25, dt=0.25 | none | L2 0.0386 pp; L2' 0.1336 pp | §4 |
| **M1 niche** | ✅ | global SC count `Q` | none (P0 fixed) | SC self-renewal +2.4 pp removed | audit #1 |
| **M2 noise** | ⚠️ conditional | zeta<1 or large σ | `ou_exact=True` | zeta=0.25: +100%; zeta=1: +14.3%; zeta=4: +3.2% | §8.1 |
| **M3 epigenetic counter** | ⚠️ | cross-kernel comparison | report tail probabilities; `cpg_sites=m0+1` | cpg vs beta `P(u<0.05)` **2.1×** | §B1 |
| **M4 quorum** | ❌ not by default | spatial paracrine | `field_interp="linear"`; `sphere_mask=True` (when √(D/k)>0.5R); **do not use the pre-audit engine** | nearest 27%→2%; pre-audit CFL D=25→42.7%, D=100→10.7% | §C2/§8.3 |
| **M5 feedback oscillation** | ⚠️ | `rate_mode="continuous"` + feedback | re-check with dt halved | default `dt·λ=0.95`; alpha=1.0→2.76 | §A3 |
| **M6 asymmetric division** | ✅ | `asymmetry f≠0.5` | none | f=0.3/0.7 variance matches theory, n conserved | §8.3 |
| **M7 metabolism/ROS** | ⚠️ | slow variable `r` coupling | dt halved; screen with `grn_stiffness` | same as M5 | §A3 |
| **M8 bet-hedging** | ⚠️ | lineage autocorrelation | `sister_corr>0`; `ou_exact=True` | sister correlation suppressed to ≈0 (should be positive) | §B2 |
| **G≥3 multigene** | ❌ not by default | G≥3 | explicit `regulation` (product=AND / additive=OR), `stemness` measure, `K_switch` | mismatch 12.8%/19.1%; AND vs OR `max|Δ|=1.63` | §B3/B4 |

---

## 3. How to use it

```python
from cerps.config import CERPSConfig, GRNConfig, EpiConfig, SpaceConfig
from cerps import batch

cfg = CERPSConfig(
    grn=GRNConfig(ou_exact=True),                  # fast OU / large sigma
    epi=EpiConfig(kernel="beta", sister_corr=0.3), # M8 lineage correlation
    space=SpaceConfig(paracrine=True,
                      field_interp="linear",        # quantitative M4
                      sphere_mask=True),            # large sqrt(D/k)
)
for msg in batch.validate(cfg):                     # read the warnings first
    print("WARN:", msg)
print("dt*lambda =", batch.grn_stiffness(cfg))      # >1 -> halve dt
```

**Five questions before a run**
1. Is G≥3? → choose `regulation` and `stemness` explicitly.
2. Is OU fast (zeta<1) or sigma large? → `ou_exact=True`.
3. Is spatial paracrine used? → `field_interp="linear"`; `sphere_mask=True`
   when `√(D/k)>0.5R`.
4. Cross-kernel or lineage-correlation work? → report tails; for M8 use
   `sister_corr>0`.
5. Is `grn_stiffness` >1? → if so, halve `dt` and re-run.

---

## 4. Known default traps

| Trap | Trigger | Consequence | Avoidance |
|---|---|---|---|
| OU variance overestimated | `zeta≲1` and `ou_exact=False` | noise amplitude too large; sigma sweep distorted | `ou_exact=True` |
| Nearest field sampling | `field_interp="nearest"` | field error ~27%, grid anisotropy | `"linear"` |
| Pre-audit CFL clamp | `spatial.py` before 2026-10-19 | D silently reduced to 43% or less | use current `main` |
| AND/OR mixed up | G≥3 keeps default `product` | qualitatively different attractors | choose explicitly |
| Independent sisters | M8 uses lineage correlation | signal systematically suppressed | `sister_corr>0` |
| Stiffness over the limit | alpha↑, s↓, lam↑ | explicit-Euler oscillation | halve `dt` |
| Cartesian grid | `√(D/k)≳0.5R` | field leaks into cell-free corners | `sphere_mask=True` |
| Transient included in transitions | `rec_tran_from=0` (**the `CERPSConfig` default**) | transition matrix not comparable with paper Table 3 | `rec_tran_from=2000` (paper protocol); `batch.validate` warns on the default |

---

## 5. Maintenance convention

- Whenever a mechanism or switch is added, **add one row** to this table
  (mechanism, regime, switches, bias, evidence).
- Every warning from `batch.validate` should map to a row here.
- Keep this table in sync with `docs/audit_2026-10-21.md`; cite the audit
  section for evidence.
