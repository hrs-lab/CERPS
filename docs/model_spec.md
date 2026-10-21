# CERPS model specification

CERPS = **C**oupled **E**pigenetic-**R**egulatory-**P**opulation-**S**patial
model.  This document specifies the mathematics; the baseline reduction to
Huang, Situ & Lei (2024) is marked **[legacy]** everywhere.

## 0. State of a cell

Each cell carries

* `x = (x_1..x_G)` gene-expression levels (protein concentrations, AU),
* `u = (u_1..u_G)` epigenetic states: fraction of activating marks
  (methylated CpG sites or modified nucleosomes) at the regulatory region
  of each gene, `u_i in [0,1]`,
* `eta = (eta_1..eta_G)` extrinsic-noise states (OU processes),
* cell-cycle state: resting (G0/G1-like) or proliferating with age `a`,
* a replicative-senescence counter `d` (continuous),
* a 3D position `r` (spatial mode).

## 1. Gene regulation network (R1... of the generalisation: R3)

For gene `i` (Hill functions `H+(x)=x^n/(s^n+x^n)`,
`H-(x)=s^n/(s^n+x^n)`):

```
dx_i/dt = a_i(u_i, eta_i) * (rho_i + (1-rho_i) H+(x_i; s_ii))
          + b_i * PROD_{j != i} H-(x_j; s_ij)^{B_ij}
          + SUM_{j != i} C_ij H+(x_j; s_ij)
          + p_i c(r)                      (paracrine response, optional)
          - k_i x_i
```

with the exponential-form noise-modulated expression rate

```
a_i(u_i, eta_i) = alpha_i * exp(lam_i u_i + sigma_i eta_i - sigma_i^2/2),
d eta_i = -(eta_i/zeta_i) dt + sqrt(2/zeta_i) dW_i.
```

**[legacy]** G=2, `B` full off-diagonal, `C=0`, `b=(1,1)`, `s_ij = s_j`:
the product term reduces to `b_1 H-(x_2)` etc., i.e. exactly Eq. (1)-(4)
of the paper.  Integration: forward Euler with dt (as the original C++).

## 2. Epigenetics

### 2a. Continuous intracellular dynamics (between divisions)

```
du_i/dt = (w0 + w1 act_i) (1 - u_i) - (e0 + e1 (1 - act_i)) u_i
```

where `act_i = rho_i + (1-rho_i) H+(x_i)` is the activation level of gene
`i` (GRN -> epigenetics coupling).  Fixed point increasing in `act_i`;
timescale 1/(w0+w1+e0+e1) ~ 10-100 h (slow vs the GRN timescale ~1 h).
**[legacy]** disabled.

### 2b. Inheritance kernels at mitosis (R1: not limited to beta)

Let `v` be the mother's state; each daughter draws independently:

| kernel | daughter state |
|---|---|
| `beta` | `u ~ Beta(m0*phi, m0*(1-phi))`, `phi(v) = m1 + m2 (m3 v)^m4/(1+(m3 v)^m4)` **[legacy]** |
| `gaussian` | `u ~ N(phi(v), sigma_g^2)` truncated to [0,1] with the pre-truncation mean solved so that `E[u]=phi(v)` exactly (rejection sampling) |
| `cpg` | per-site: marked site retained w.p. `p_maintain`; unmarked site gained w.p. `p_denovo`; `u = #marked / M` with `M` explicit sites |
| `deterministic` | `u = phi(v)` |

The `cpg` kernel is mechanistic (semi-conservative maintenance + de-novo
methylation + passive loss); its single-generation mean map is
`u' = u p_maintain + (1-u) p_denovo` with site-level binomial noise
`sd = sqrt(u'(1-u')/M)`.

## 3. Division (R2)

At mitosis the mother is replaced by two daughters.

* **Expression partitioning**
  - `exact` **[legacy]**: both daughters inherit `x`.
  - `binomial`: with `n_i = round(V x_i)` molecules, daughter 1 receives
    `n_i ~ Binomial(n_i, f)`, daughter 2 the complement; daughter
    concentrations are `n_i / (f V)` and `n_i' / ((1-f) V)`
    (Huh & Paulsson 2011).  Variance of the daughter concentration:
    `Var = x/(V f(1-f)) / 4` (for f=1/2: sd = sqrt(x/V)).
* **Epigenetic state**: drawn from the inheritance kernel (2b).
* **Noise state**: daughters inherit `eta` (as the original C++).
* **Position** (spatial mode): daughters at `r +/- R*axis` with a random
  axis.

## 4. Cell kinetics (R4)

### 4a. Legacy mode (Table 1 of the paper)

Discrete classification (`x>thr` / `x<thr`, keep-previous when both low,
exactly as the released C++ code):

| phenotype | beta | kappa | mu | tau |
|---|---|---|---|---|
| SC | `beta0 theta0/(theta0+Q)` | `kappa0` | `mu0` | `tau0` |
| TA1/TA2 | `beta0` | `2 kappa0` | `mu0` | `tau0/2` |

with `Q` = number of SC; TA cells die at the 15th division (counter
incremented at mitosis, reset to 0 when a resting cell is classified SC).

### 4b. Continuous mode

Let `S = PROD_i H+(x_i; s_ii)` be the **stemness** (all-genes-on measure)
and `H_S(v) = v^p/(K^p + v^p)` a steep Hill switch (default p=12).  With
`K_switch=None` (default) `K = PROD_j H+(x_threshold)` so the continuous
boundary coincides with the legacy box on the diagonal (K=0.25 for G=2,
0.125 for G=3); the two rules still differ off the diagonal, and that
mismatch grows with G (6.3% / 12.8% / 19.1% of x-space for G=2/3/4).

```
beta  = beta0 * [1 - H_S(S) * (1 - theta0/(theta0 + Q_soft))] * CI(rho)
kappa = kappa0 * (1 + H_S(1 - S))
mu    = mu0
tau   = tau0 * (0.5 + 0.5 H_S(S))
Q_soft = SUM_cells H_S(S)                     (soft stem-cell count)
CI     = theta_rho/(theta_rho + max(0, rho_local - rho0))   (spatial mode)
```

Senescence: at each mitosis `d <- d + (1 - H_S(S))`; while resting
`d <- d (1 - H_S(S))`; a cell whose post-increment counter reaches
`max_div_ta` dies at that mitosis (in the sharp limit this is exactly the
legacy 15-division rule with SC reset).

**Sharp-limit consistency** (verified in tests): with p -> inf the
continuous rules converge to 4a cell-by-cell.

## 5. Coupling structure (R5)

```
        epigenetics u  --(a_i = alpha e^{lam u})-->  GRN x
             ^  |                                        |
   continuous |  | inheritance                        rates
   dynamics   |  | kernel                              (beta,kappa,mu,tau)
             |  v                                        |
             GRN x  --(activity-dependent writing)-->    population
                                                        /    ^
              population <--(niche Q, density, paracrine)--/
```

All three pairs are bidirectional:
epigenetics<->GRN (Sec. 1 and 2a),
GRN/epigenetics<->population (Sec. 4; paracrine field Sec. 6).

## 6. Space (R6)

* Off-lattice positions, soft-sphere overlap relaxation
  (`dr = relax * overlap` per step), optional Brownian motility,
  reflecting spherical boundary.
* Contact inhibition: local neighbour count `rho_local` within
  `1.6 * 2R` reduces `beta` (Sec. 4b).
* Paracrine field `c(r,t)`: `dc/dt = D_c laplacian c - lambda c +
  SUM_cells q(S_i) delta(r - r_i)` on a regular grid, solved by FTCS with
  **CFL sub-stepping** (the nominal `D` is honoured; it is never silently
  clamped) and **reflecting (zero-flux) boundaries**.  `q(S) = q0 S`
  (stem-like cells secrete); cells sense `c` which adds `p_i c` to gene-`i`
  production.  Sampling/deposition is nearest-node by default
  (`field_interp="nearest"`, O(h)); `field_interp="linear"` switches to
  trilinear sampling and conservative trilinear source deposition.

## 7. Event scheme

Per time step dt (default 0.25 h, as the paper): update noise, x, u for
all cells; resting cells enter proliferation with probability
`beta dt`, are removed with probability `kappa dt`; proliferating cells
apoptose with probability `mu dt`, mature for `tau` and then divide.
This is the exact scheme of the original C++ implementation.

## 8. Implementation map

```
cerps/config.py      all parameters (dataclasses)
cerps/grn.py         N-gene dynamics, OU noise
cerps/epigenetics.py continuous dynamics + 4 inheritance kernels
cerps/division.py    exact / binomial expression partitioning
cerps/rates.py       legacy Table-1 and continuous rates, stemness
cerps/spatial.py     3D mechanics, contact inhibition, paracrine field
cerps/simulator.py   individual-based event loop, records, snapshots
tests/test_cerps.py  37 unit/integration/cross-validation tests
examples/            runnable demonstrations of each feature
```
