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
          + a_tf_i                        (externally induced TF, Fig. 11/12)
          + p_i c(r)                      (paracrine response, optional)
          - k_i x_i
```

Two cross-term compositions are available (`GRNConfig.regulation`):
`"product"` (default, the form above; identical to the additive one for
`B_ij in {0,1}` and G = 2) and `"additive"` (`SUM_j B_ij H-(x_j)`, the
paper's textual form, in which `b_i` is not applied).  `b_i = 1` unless
`GRNConfig.b_base` is given.

with the exponential-form noise-modulated expression rate

```
a_i(u_i, eta_i) = alpha_i * exp(lam_i u_i + sigma_i eta_i - sigma_i^2/2),
d eta_i = -(eta_i/zeta_i) dt + sqrt(2/zeta_i) dW_i.
```

The OU is integrated by Euler-Maruyama (the C++ scheme, default) or, with
`GRNConfig.ou_exact=True`, by the exact one-step transition
`eta <- eta e^{-dt/zeta} + sqrt(1-e^{-2dt/zeta}) N(0,1)` whose stationary
variance is exactly 1 for every `dt`/`zeta` (Euler-Maruyama inflates it to
`1/(1-dt/(2 zeta))`: +14.3% at `zeta = 1`, `dt = 0.25`, +100% at
`zeta = 0.25`).

**[legacy]** G=2, `B` full off-diagonal, `C=0`, `b=(1,1)`, `s_ij = s_j`:
the product term reduces to `b_1 H-(x_2)` etc., i.e. exactly Eq. (1)-(4)
of the paper.  Integration: forward Euler with dt (as the original C++).

## 2. Epigenetics

### 2a. Continuous intracellular dynamics (between divisions)

```
du_i/dt = (w0 + w1 act_i) (1 - u_i) - (e0 + e1 (1 - act_i)) u_i
```

`act_i` is the activation level driving the mark dynamics.  **What the code
actually feeds in** (`simulator.py`) is the pure Hill activation
`act_i = H+(x_i)`, *not* the `rho_i + (1-rho_i) H+(x_i)` combination used in
the production term of Sec. 1.  This is a modelling choice (a sharper,
threshold-like mark response), not a typo: the `rho` mixture was introduced
for the expression ODE only.  Fixed point `u* = (w0+w1 act)/(w0+w1 act+e0+e1(1-act))`
increases with `act_i`; relaxation rate `w0+w1 act + e0+e1(1-act)` (with the
defaults `w1 = e1` this is the constant `w0+w1+e0+e1 = 0.14 /h`, i.e. a
~7.1 h timescale, slower than the GRN timescale `k_i ~ 1 /h`).  **[legacy]**
disabled (`EpiConfig.continuous=False`); the paper freezes `u` between
divisions.

### 2b. Inheritance kernels at mitosis (R1: not limited to beta)

Let `v` be the mother's state; each daughter draws independently:

| kernel | daughter state |
|---|---|
| `beta` | `u ~ Beta(m0*phi, m0*(1-phi))`, `phi(v) = m1 + m2 (m3 v)^m4/(1+(m3 v)^m4)` **[legacy]** |
| `gaussian` | `u ~ N(phi(v), sigma_g^2)` truncated to [0,1] with the pre-truncation mean solved so that `E[u]=phi(v)` exactly (rejection sampling) |
| `cpg` | `M` explicit markable sites per gene; two modes (see below) |
| `deterministic` | `u = phi(v)` |

`cpg` has two modes (`EpiConfig.cpg_mode`):

* `"moment"` (**default**) — the per-site marking probability is set to
  `phi(v)`, so `E[u] = phi(v)` and `Var[u] = phi(1-phi)/M`.  With
  `M = m0+1` this matches the beta kernel in **both** moments, which makes
  it a controlled noise-only variant; with `M != m0+1` only the mean matches
  (`batch.validate` warns).
* `"mechanistic"` — the classical semi-conservative rule: a marked site is
  retained w.p. `p_maintain`, an unmarked one gains w.p. `p_denovo`.  Its
  single-generation mean map is the affine `u' = u p_maintain + (1-u) p_denovo`
  and therefore has a **different** fixed point from the beta kernel
  (0.143 vs 0.132 for the default parameters); use it for mechanistic
  studies, not for controlled comparisons.

`EpiConfig.sister_corr` (default 0 = independent sisters, the paper's
kernel) induces a controlled correlation between the two daughters via a
Gaussian copula on the beta kernel (latent correlation `rho`, exact beta
marginals) and correlated per-site uniforms for `cpg(moment)`.  It is not
implemented for `gaussian`, where `inherit()` raises instead of silently
returning independent sisters.

## 3. Division (R2)

At mitosis the mother is replaced by two daughters.

* **Expression partitioning**
  - `exact` **[legacy]**: both daughters inherit `x`.
  - `binomial`: with `n_i = round(V x_i)` molecules, daughter 1 receives
    `n_i` molecules out of `n_i` total split with fraction `f`
    (`DivisionConfig.asymmetry`), daughter 2 the complement; daughter
    concentrations are `n_1/(f V)` and `n_2/((1-f) V)` (Huh & Paulsson 2011).
    With the partitioning-noise strength `alpha` (Soltani et al. 2016 Eq. 5)

    ```
    E[x_daughter | x]   = x
    Var[x_daughter | x] = alpha * x * (1-f) / (f * V)
    ```

    so `alpha = 0` is deterministic, `alpha = 1` the exact binomial draw,
    `alpha > 1` a beta-binomial (over-dispersed, clustered segregation) and
    `0 < alpha < 1` a variance-matched normal approximation (under-dispersed,
    active ordering).  For `f = 1/2`, `alpha = 1`: `Var = x/V`, i.e.
    `sd = sqrt(x/V)`.  Note that daughter 1's variance is *not* symmetric in
    `f` -- the smaller-volume daughter inherits the larger variance.
    `batch.validate` warns when `alpha != 1` is combined with
    `partition="exact"`, which silently ignores it.
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
0.125 for G=3, 0.0625 for G=4); the two rules still differ off the
diagonal -- the legacy rule is a box, the product rule a curved region --
and the mismatch grows with G: **5.9% / 12.6% / 19.4%** of the cube
`[0,3]^G` for G = 2 / 3 / 4 (measured, unavoidable for any
single-threshold smoothing of a box).

```
beta  = beta0 * [1 - H_S(S) * (1 - theta0/(theta0 + Q_soft))] * CI(rho)
kappa = kappa0 * (1 + H_{1-K}(1 - S))          <-- complementary threshold
mu    = mu0
tau   = tau0 * (0.5 + 0.5 H_S(S))
Q_soft = SUM_cells H_S(S)                     (soft stem-cell count)
CI     = theta_rho/(theta_rho + max(0, rho_local - rho0))   (spatial mode)
```

Senescence: at each mitosis `d <- d + (1 - H_S(S))`; while resting
`d <- d (1 - H_S(S))`; a cell whose post-increment counter reaches
`max_div_ta` dies at that mitosis (in the sharp limit this is exactly the
legacy 15-division rule with SC reset).

The TA-side switches use the **complementary** threshold `1-K`, not `K`:
with the same `K` the two sides fire at `S = K` and `S = 1-K`, leaving a
band in which a cell is neither properly SC-like nor properly TA-like (at
`K = 0.125` a stem-like cell with `S = 0.5` would receive the TA removal
rate `kappa = 2 kappa0`).  The original hard-coded `K_switch = 0.5` hid
this because there `1-K == K`.

**Sharp-limit consistency** (verified in tests): with p -> inf the
continuous rules converge to 4a cell-by-cell.

Two history-dependent extensions are orthogonal to everything above:

* **P6 non-Markovian memory** (`MemoryConfig`): a per-cell leaky integrator
  `mem <- mem + min(1, dt/tau) (drv - mem)` of a driver signal
  (`"stemness"` or `"u_mean"`; any other name raises `ValueError`),
  optionally modulating `beta` by `1 + mem_couple (mem - 0.5)`.  Both
  daughters inherit the parent's `mem`.  `enabled=False` keeps the
  validated baseline untouched.
* **P8 feedback-time delay** (`RatesConfig.delay_steps`, in steps):
  the niche signal `Q` fed to `beta` is sampled `delay_steps` steps in the
  past (`rates.effective_niche_signal`).  `0` (default) is exactly the
  paper's instantaneous feedback.

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

* Off-lattice positions, soft-sphere overlap relaxation.  The overlap decays
  exponentially with the mechanical relaxation time,
  `overlap <- overlap * exp(-dt/tau)` with `tau = -0.25/ln(1-relax)`
  (so the legacy per-step fraction `relax` is reproduced exactly at the
  nominal `dt = 0.25 h` and the mechanics become dt-consistent); set
  `SpaceConfig.relax_time` explicitly for dt-independent mechanics.
  Optional Brownian motility (`pos += motility sqrt(dt) N`), reflecting
  spherical boundary at `r = domain_radius`.
* Contact inhibition: local neighbour count `rho_local` within
  `inhibition_radius * 2R` (default `1.6 * 2R`) reduces `beta` (Sec. 4b).
* Paracrine field `c(r,t)`: `dc/dt = D_c laplacian c - lambda c +
  SUM_cells q(S_i) delta(r - r_i)` on a regular grid, solved by FTCS with
  **CFL sub-stepping** (the nominal `D` is honoured; it is never silently
  clamped) and **reflecting (zero-flux) boundaries**.  `q(S) = q0 S`
  (stem-like cells secrete); cells sense `c` which adds `p_i c` to gene-`i`
  production.  Sampling/deposition is nearest-node by default
  (`field_interp="nearest"`, O(h) error, mean ~27% at `h = 4`);
  `field_interp="linear"` switches to trilinear sampling and conservative
  trilinear source deposition (error ~2%).  `sphere_mask=True` forces the
  field to zero outside the spherical cell domain, removing the
  cell-free-cube-corner reservoir that appears once the decay length
  `sqrt(D/k)` approaches the domain radius.

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
cerps/rates.py       legacy Table-1 and continuous rates, stemness, P6/P8
cerps/spatial.py     3D mechanics, contact inhibition, paracrine field
cerps/simulator.py   individual-based event loop, records, snapshots
cerps/batch.py       sweeps, statistics, numerical-health checks
cerps/h_general.py   proposal layer F = (G, H) with epigenetic gates
cerps/notation.py    (x,u) vs (u,v) naming alias layer
tests/               191 standalone checks / 69 pytest cases
examples/            runnable demonstrations of each feature
```
