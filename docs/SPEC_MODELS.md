# Model Specifications — Original Model & CERPS
(Verified against paper text `upload/paper_full.txt` on 2026-10-17. This file is the
single source of truth for all implementations. Keep in sync.)

## 1. Original model (Huang, Situ & Lei 2024, JTB 577:111664)

### 1.1 GRN dynamics (per cell, 2 genes)
```
dx1/dt = a1*(rho1 + (1-rho1)*H(x1)) + b1*R(x2) - k1*x1
dx2/dt = a2*(rho2 + (1-rho2)*H(x2)) + b2*R(x1) - k2*x2
H(x) = x^n/(s^n + x^n)          <-- CRITICAL: s is RAISED to n (50% eff conc: H(s)=0.5)
R(x) = s^n/(s^n + x^n)          <-- repression (decreasing in other gene's x)
a_i = alpha_i * exp(lam_i*u_i) * exp(sig_i*eta_i - sig_i^2/2)
d eta_i = -(eta_i/zeta_i) dt + sqrt(2/zeta_i) dW_i     (OU, stationary var = 1)
```
With TF induction (Eq 9, Fig 11/12): a1 -> (a1 + aTF) in the SELF-ACTIVATION term only.

### 1.2 Table 2 parameters (paper values)
b1=b2=1, k1=k2=1, s1=s2=0.5, rho1=rho2=0.1, sig1=sig2=0.05, zeta1=zeta2=1,
alpha1=alpha2=0.4, lam1=lam2=1.9, n=4, m0=60, m1=0.08, m2=0.8, m3=2.0, m4=2.0,
theta=200, beta0=0.04 /h, kappa0=0.01 /h, mu0=0.0002 /h, tau0=4 h, dt=0.25 h.

### 1.3 Phenotypes (Eq 5). Type codes: SC=0, TA1=1, TC=2, TA2=3
```
SC : x1>0.5 and x2>0.5      TA1: x1>0.5 and x2<=0.5
TA2: x1<=0.5 and x2>0.5     TC : x1<=0.5 and x2<=0.5
```
Label updated ONLY during resting phase (each dt step, after x-update).
When label becomes SC -> division counter div resets to 0.

### 1.4 Kinetic rates (Table 1)
- SC : beta = beta0*theta/(theta+Q), kappa=kappa0, mu=mu0, tau=tau0   (Q = #SC, computed after label updates, before fates)
- TA1/TA2: beta=beta0, mu=mu0, tau=tau0/2, kappa = 2*kappa0 if div<15 else INF (immediate removal)
- TC: beta=kappa=mu=tau=0 (frozen until label changes; TC is transient)

### 1.5 Cell cycle algorithm (per dt step, per cell)
1. Integrate x by Euler dt; update OU eta exactly: eta <- eta*exp(-dt/zeta) + sqrt(1-exp(-2*dt/zeta))*N(0,1)
2. If resting: update label from x (reset div if SC)
3. Compute rates (needs Q)
4. Resting: with p=kappa*dt -> remove; elif p<kappa*dt+beta*dt -> enter proliferating
   (age=0; if label TA: div += 1)   [single uniform draw]
5. Proliferating: if age>=tau -> MITOSIS; elif rand<mu*dt -> remove; else age+=dt
6. Mitosis: daughters inherit x (both, identical), eta, div, label; birth_type = mother label;
   each daughter u_i ~ Beta(g_i,h_i) INDEPENDENTLY per gene per daughter:
   phi_i = m1 + m2*(m3*v_i)^m4/(1+(m3*v_i)^m4)   (v = mother u)
   g_i = m0*phi_i, h_i = m0*(1-phi_i)   [sample via 2 gammas: u = G1/(G1+G2)]
   Daughters start resting, age=0. Transition log at mitosis: trans[birth_type, label] += 1.
7. fixed_u mode (Fig 3): daughters inherit u exactly (no resampling).

### 1.6 Verified steady states / consistency checks (a = 0.4*exp(1.9*u))
- SC: u=(0.45,0.45) -> a=0.98 -> x≈(0.98,0.98) (both>0.5 ✓)
- TA1: u=(0.36,0.24) -> x≈(1.75,0.07); TA2 mirror
- SC attractor exists iff a>0.75 ⟺ u>0.331 (per gene); TA states exist for a<1.65
- SC->TA per-division probability ≈ P(u_daughter < 0.331) ≈ 13%/gene ✓ matches Table 3
- phi fixed points: hi≈0.405, lo≈0.175 (iterate u <- 0.08+0.8*(2u)^2/(1+(2u)^2))

### 1.7 Experiment protocols (Fig 2-13, Table 3)
- Table 3: 30 runs, init 100 cells x~U(0,3), u~U(0,1), resting; T=2000h(or 4000h);
  transitions logged at mitosis (birth_type -> settled label of dividing cell),
  row-normalize per run, average across runs. TARGETS: SC row [73.94, 12.99, -, 13.07],
  TA rows diag >99.99%. (Previous session port hit 74.5% on seed 12.)
- Fig 2A: deterministic bifurcation over a in (0.3,2.2), a1=a2=a; find all FPs
  (2D Newton from multi ICs), classify stability by Jacobian. Stable=solid, unstable=dashed.
- Fig 2B: grid (u1,u2), a_i=0.4*exp(1.9*u_i); color = # stable steady states (1/2/3).
- Fig 2C: single cell, u=(0.34,0.34), sigma=0.05, T=100h, record x traj. No division.
- Fig 2D: grid (u1,u2), start at SC attractor x(a(u)), fixed u, sigma=0.05,
  run until label != SC (censor 2000h), avg over ~20 reps.
- Fig 3: fixed_u populations, u in {(0.5,0.5),(0.5,0.6),(0.6,0.5),(0.6,0.6)},
  sigma in {0,0.1,...,1.0} (11 pts), 100 runs, T=1000h, init 100 SC at x=SS(a(u)).
  Final SC/TA ratios. TA emerges only above noise threshold (~0.15 or ~0.25).
- Fig 4: 4 param sets (alpha1,alpha2)={(0.4,0.4),(0.4,1.0),(1.0,0.4),(1.0,1.0)},
  30 runs, T=2000h, count dynamics.
- Fig 5: snapshots (x1,x2),(u1,u2) at t=5,50,100,1000 from Fig4A runs.
- Fig 6: transition trajectories: 2-pass same-seed runs; pass1 detect label-change
  events (t_event), pass2 record all cells (id,t,x1,x2,label) in window; extract cell trace.
- Fig 7: (sigma,m2) = (0.20,0.80),(0.05,0.75),(0.05,1.00) single illustrative runs.
- Fig 8: sigma in {0,0.02,...,0.1} (6 pts) x 30 runs, T=4000h, stats 2000-4000h.
- Fig 9: m2 in {0.80,...,0.90} (6 pts) x 30 runs, same protocol.
- Fig 10: asymmetry: vary k1,s1,alpha1,rho1 (x0.5/x1.5 etc), ~10 variants x 30 runs.
- Fig 11: aTF in {0,0.05,...,0.5} (11 pts) x {TA1-init, TA2-init} x ~27 runs, T=10000h.
  Init 500 cells at TA steady state. Success = SC count >= 100 (early stop). Fig11C from
  TA2 + aTF=0.4-0.5. TA1-init: P(success) stays low ~ basal; TA2-init: P->1 by aTF>0.4,
  t_success ~ 700h (TA1) / ~500h (TA2, aTF>0.4).
- Fig 12: init 100 SC (x=SS(u=0.45)), aTF in {0,0.05,...,0.2} (5 pts) x 30 runs,
  T=4000h, stats 2000-4000h.
- Fig 13: Waddington U(t,z) = -log f(t,z), z=x1-x2, from x-snapshots (every 50h):
  (A) sigma=0.02 vs 0.1; (B) m2=0.8 vs 0.9; (C) aTF=0 vs 0.2.

## 2. CERPS model (Coupled Epigenetic-Regulatory-Population-Spatial)

### 2.1 Five extensions over original
E1 continuous intracellular epigenetic dynamics (not only at division, not beta-only)
E2 division kernel: lineage-dependent maintenance fidelity + TF-guided de novo + binomial
   mark allocation (non-beta); multi-gene GRN (G genes, repression matrix T)
E3 NO discrete cell classes: rates depend continuously on cell state (balance score S)
E4 bidirectional coupling epi <-> GRN <-> population (local density feedback)
E5 3D off-lattice space (sphere domain, niche, motility, contact inhibition)

### 2.2 State per cell
x[G] expression, u[G] epigenetic (fraction of activating marks), eta[G] OU noise,
div (lineage division count, daughter = mother+1 at birth), age, phase, pos[3].

### 2.3 Interphase epigenetic dynamics (birth-death on mark fraction)
```
sig_i = sigmoid(g_grn*(x_i/2 - theta_grn))            # TF activity proxy
du_i/dt = (w0 + w1*sig_i)*(1-u_i) - e0*u_i
u*_i = (w0+w1*sig_i)/(w0+w1*sig_i+e0); relaxation time tau_epi = 1/(w0+w1*sig+e0) ~ 39h
```
Check: SC-like (x=1, sig=0.29): u*≈0.42 (observed 0.47 with mild erosion);
ON gene (x~1.8, sig~0.93): u*≈0.61; OFF gene: u*≈0.26. Documented steady states:
SC x=(1,1) u=(0.47,0.47); TA1 x=(1.74,0.07) u=(0.33,0.23); TA2 mirror.

### 2.4 Division kernel
```
q(d) = q_inf + (q_keep-q_inf)*exp(-d/d_erode)        # maintenance fidelity, lineage depth d
p_dn_i = p_dn0 + p_dn1*sig_i                          # TF-guided de novo marking prob
per daughter, per gene: u_d = Binomial(c_epi, q*u_m + p_dn_i*(1-u_m))/c_epi
E[u_d] = q*u_m + p(1-u_m); erosion fixed point u** = p/(1-q+p) << u* (interphase target)
```
Example: SC-like sig=0.29 -> p≈0.00416, q=0.92 -> u**≈0.05 vs u*≈0.42.

### 2.5 Continuous rates (no discrete types)
```
S = (geo(x)/mean(x))^2    (G=2: S = 4*x1*x2/(x1+x2)^2)   # balance/plasticity score, 1=balanced SC-like
beta = beta0*(eps_b + (1-eps_b)*S^p_sens)/(1+(rho_loc/rho50)^h_inhib)   # contact inhibition
tau  = tau0*(tau_min_frac + (1-tau_min_frac)*S)
kappa = kappa0*(1 + c_kappa*(1-S) + c_sen*(1-S)*max(0,div-d_half)/d_half) * niche_prot
   # senescence term acts ONLY on low-S (polarized) cells with deep division history
   # (replaces the hard 15-division TA cap of the original model)
mu = mu0 (const).  niche_prot = 1 - 0.5*(r_pos < r_niche)   # protective niche
```

### 2.6 Calibrated parameters (previous session; validation targets in 2.8)
GRN: G=2, rho=.1, b=1, k=1, s=.5, n_hill=4, alpha=.4, lam=1.9, zeta=1, sigma=.05
epi: w0=.005, w1=.0198, e0=.015, g_grn=8, theta_grn=.6 (x normalized x/2), c_epi=60, d_epi=1.0
div: q_keep=.92, q_inf=.80, d_erode=45, p_dn0=.003, p_dn1=.004
rates: mu0=.0002, kappa0=.01, beta0=.04, tau0=4, tau_min_frac=.5, eps_b=1.0, p_sens=2,
       c_kappa=1.0, c_sen=.5, d_half=60
space: r_dom=30, r_niche=10, rho50=.003, h_inhib=3, diff_mot=.5, sigma_div=2.5,
       r_loc=10 (neighbor radius), grid cell 10

### 2.7 Spatial module
Sphere radius r_dom; hash grid (cell size 10); local density rho_loc =
#neighbors within r_loc / (4/3*pi*r_loc^3), recomputed every 4 steps (hourly);
motility: resting cells pos += sqrt(2*diff_mot*dt)*randn per axis;
daughters: mother pos + randn*sigma_div per axis; reflect into sphere at boundary.
Well-mixed variant: rho_loc = N/V_dom. Niche: kappa halved inside r_niche.

### 2.8 Injury/protocol module + validation targets (M1, 12 seeds each, T=13000h)
- control: no injury. TARGET: stable, bal≈49.3% at end
- acute: single ablation (kill 90%) at t=6000. TARGET: recovers, bal≈49.8%
- serial4: 90% ablations at 6000,7500,9000,10500. TARGET: bal≈45.2%
- serial6: 6 ablations every 1500h from 6000. TARGET: bal≈23.9%
- harsh (M2/M4 protocol): rounds_period=1500, frac=0.9, continuing -> 3/6 arms collapse,
  bal@13k≈25.9%, collapse ~8630h
- bal = fraction of cells with x1>0.5 and x2>0.5 (SC-analog). Also record nhighu
  (u1>0.5 & u2>0.5), meanS, meandiv, ndiv.
- stress variant: beta *= beta_stress in windows (t_stress_on/off/dur/period)

### 2.9 Records
counts hourly (7 cols): [t, N, nbal, nhighu, meanS, meandiv, ndiv_cum]
div_rec per division (9 cols): [t, d, q, u1_m, u2_m, S_m, u1_d, u2_d, S_d] (daughter 1)
snaps at snap_times: [x1,x2,u1,u2,px,py,pz] per cell (cap 3000/snap)

### 2.10 Mechanism (mitotic erosion -> regenerative debt -> bankruptcy)
Interphase pulls u up toward u* (TF-guided writing); every division pulls u toward
u** = p/(1-q+p) ≈ 0.05 (maintenance failure). Quiescent tissue: rare divisions,
u ≈ u*, high S -> no senescence. Regeneration: bursts of divisions ->
(i) each division erodes u by (1-q)*(u-u**); (ii) lineage depth d grows -> q(d)->q_inf
(0.80) -> stronger erosion per division. u drifts below ~0.34 -> a=0.4e^{1.9u} < 0.75
-> SC attractor vanishes -> x polarizes -> S drops -> senescence term (low-S, deep-d)
removes cells faster than replenishment -> "replication debt" accumulates ->
regenerative bankruptcy (collapse) under serial but not acute injury.
Counter-intuitive claim: regeneration itself, not the injury, erodes epigenetic memory;
q(d) lineage decay (Beerman 2013 proliferation-history methylation erosion) is the
key nonlinearity (removing q(d) decay -> serial injury tolerated).

### 2.11 Legacy cross-validation mode
CERPS kernel flags: legacy_rates (orig Table-1 discrete rates via classify + maxdiv 15 +
beta_SC=beta0*theta/(theta+Q)), legacy_inherit (Beta(m0*phi,m0*(1-phi)) u-sampling),
no interphase epi dynamics, well_mixed. Must reproduce orig Table 3 within noise.

## 3. Experiment suites (CERPS)
- M1 mechanism (above, 12 seeds x 4 arms)
- M2 coupling ablations, harsh serial protocol, 8 variants x 12 seeds:
  {full, no_epi_to_grn (lam=0), no_grn_to_epi (w1=0), no_div_kernel (q=1 const, p_dn=0),
   no_lineage_decay (q=q_keep const), no_density_fb (rho50=inf), no_senescence
   (c_kappa=c_sen=0), well_mixed (no space)}
  -> which coupling link is necessary/sufficient for debt accumulation
- M3 spatial: 8 variants (full vs well-mixed x niche on/off x motility x density-fb)
- M4 dose-response: ablation frac f in {0.5,...,0.95} x 12 seeds, period 1500
- M5 aging: young (d_half=60) vs aged (d_half=30 or c_sen=1.0) x 10 seeds
- SEN: LHS 160 points, 12 params +-30% (w1,e0,q_keep,q_inf,d_erode,p_dn1,c_kappa,c_sen,
  d_half,rho50,beta0,kappa0) -> PRCC on outputs (bal@end, collapse time, meanS)
- ROB: 9 perturbation variants x 8 seeds (+-20% key params)

## 4. Data & inference
- Paul et al. 2015 (GSE72857) scRNA-seq UMI: GSE72857_umitab.txt.gz (16413252 bytes),
  URL: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1833856&format=file&file=GSE72857_umitab.txt.gz
  (also experimental_design file). Summary stats (previous session, recomputed next
  session): Sfpi1/Gata1 z-quantiles, corr=-0.27, GMM 3-comp weights [0.42,0.36,0.22],
  pro_frac=0.234. Map: gene1=Sfpi1-like, gene2=Gata1-like.
- Beerman et al. 2013 (GSE44117) RRBS methylation, HSC proliferation history
  (5-FU 1x/2x/4x, serial transplant, aged): 12 bed.gz files, file list saved in
  data/GSE44117_filelist.txt; URL pattern
  https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSMxxxxxx&format=file&file=NAME
  NOTE: 135070-byte files are error HTML -> delete & retry. NCBI https ~73KB/s.
  Use: methylation drift vs proliferation history -> constrain (q_keep,q_inf,d_erode,p_dn).
- SMC-ABC (Toni 2009): params (w1,e0,q_keep,d_erode,p_dn1,c_epi), summary stats from
  short CERPS runs (3000h, from-100 init): Paul stats (x-quantiles, corr, GMM weights,
  pro frac) + Beerman stat (u vs div slope). Synthetic-recovery test first; then fit
  Paul; cross-predict Beerman stat; posterior 1D/2D plots. Budget: N=50-100 particles,
  5-8 generations, 2 cores parallel.
