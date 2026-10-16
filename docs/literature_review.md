# Literature review underlying the CERPS design

This review surveys the modelling literature that motivates each component
of CERPS (Coupled Epigenetic-Regulatory-Population-Spatial model).  The
baseline model is the hybrid framework of Huang, Situ & Lei (2024)
[J. Theor. Biol. 577:111664], which couples a 2-gene regulation network,
stochastic beta-distributed inheritance of epigenetic states at cell
division, and a G0/G1 cell-cycle model of regeneration.  CERPS generalises
that model along six axes (intracellular epigenetic dynamics, expression
partitioning at mitosis, multi-gene networks, continuous kinetics,
bidirectional coupling, 3D space); each generalisation is grounded in the
literature below.

## 1. Epigenetic inheritance and intracellular epigenetic dynamics

* **Dodd, Micheelsen, Sneppen, Thon (2007)**, "Theoretical analysis of
  epigenetic cell memory by nucleosome modifications", *Cell* 129:813-822.
  The canonical model of cooperative histone-modification memory: positive
  feedback between neighbouring nucleosomes produces bistable
  (active/repressed) chromatin domains that are robust to noise and
  restored after replication.  Motivates CERPS's *continuous* intracellular
  epigenetic dynamics (mark writing/erasing with cooperativity through the
  gene-activity readout) in addition to division-time inheritance.
* **Probst, Dunleavy, Almouzni (2009)**, "Epigenetic inheritance during the
  cell cycle", *Nat. Rev. Mol. Cell Biol.* 10:192-206.  Semi-conservative
  inheritance: parental marks are diluted at replication and re-established
  by reader-writer complexes; the source of the random transition of
  epigenetic states modelled at division.
* **Huang & Lei (2018, 2019)** (+ the present paper): computational models
  of histone-modification inheritance showing that the daughter-cell
  modification level is captured by a **beta distribution** whose mean
  depends on the mother's state - the kernel used in the baseline model.
* **Zhao et al. (2021)**, "Mathematical modeling of histone modifications
  reveals dynamic roles of H3K4me3 and H3K27me3" (bivalent domains);
  **Ramakrishnan et al. (2022)**, threshold-filling of histone marks
  post-replication (information-theoretic).  Confirm that (i) antagonistic
  marks, (ii) activity-dependent writing and (iii) finite-site stochasticity
  are the key ingredients - CERPS's continuous dynamics and its explicit
  CpG-site kernel mirror these ingredients.
* **Jones & Liang (2009)**, "Rethinking how DNA methylation patterns are
  maintained" (*Nat. Rev. Genet.*); **Moore, Le & Fan (2013)**, "DNA
  methylation and its basic function" (*Neuropsychopharmacology*):
  semi-conservative maintenance methylation (DNMT1/UHRF1) plus de-novo
  methylation (DNMT3) and passive dilution - the mechanistic template of
  CERPS's `cpg` inheritance kernel (per-site maintenance / de-novo /
  loss probabilities).
* **Ashwin & Sasai (2015)**, effects of collective histone-state dynamics
  on reprogramming kinetics; **Sediyara et al.** and follow-ups on
  epigenetic noise and lineage fidelity.

## 2. Gene-expression partitioning at mitosis

* **Huh & Paulsson (2011)**, "Non-genetic heterogeneity from stochastic
  partitioning at cell division", *Nat. Genet.* 43:95-100 (see also
  Huh et al. 2010, *PNAS* 107:6824-6830): random segregation of molecules
  creates daughter-cell variance scaling as ~1/n; symmetric division does
  not imply equal daughters.  CERPS's `binomial` partition option
  implements exactly this (molecule counts n_i = V x_i, daughter 1 receives
  Binomial(n_i, f)), replacing the baseline's exact-symmetry assumption.
* **Landgraf et al. (2012)**, segregation of proteins at division
  (*PLoS Biol.*): quantitative confirmation of largely binomial
  partitioning for cytosolic proteins.
* **Cadart et al. (2018)**, size control in mammalian cells
  (*Nat. Commun.* 9:3275) and **Cadart et al. (2022)** (*eLife* 11:e70816);
  **Bjorklund (2019)** review; **Zatulovskiy & Skotheim (2020)** (paper's
  own reference): cell-volume growth and its coupling to the cycle -
  the reason CERPS keeps an optional explicit molecule/volume bookkeeping
  rather than assuming constant concentrations.

## 3. Multi-gene regulatory networks and cell fate

* The baseline 2-gene toggle (PU.1-GATA1 / SOX2-OCT4-type motifs) and the
  additive vs multiplicative regulation formalisms are discussed in
  **Jolly et al. (2015)** and **Xu et al. (2014)** (both cited by the
  paper).  For G > 2 the *multiplicative* composition used by CERPS
  (b_i * prod_j H-(x_j)) is the natural extension - it reduces exactly to
  the paper's Eq. (1) for G = 2 and avoids the artefact of additive
  basals from non-expressed repressors.
* **Huang, S. (2012)**, "The molecular and mathematical basis of
  Waddington's epigenetic landscape", *Phys. Biol.* 9:045002 (with
  **Huang, Kauffman & Shmulevich 2007**): GRNs define the landscape;
  attractors = cell types.  CERPS's multi-gene module is attractor-based in
  the same sense.
* **Wu et al. (2017)**, synthetic quadrastable gene network
  (*eLife* 6:e23702): engineered >2-gene fate networks exist and show the
  expected multistability - validating multi-gene generalisations.
* **Moris, Pina & Arias (2016)**, transition states in epigenetic
  landscapes (*Nat. Rev. Genet.*); **Mojtahedi et al. (2016)**
  (*PLoS Biol.*): high-dimensional critical transitions in fate decisions.
* **Rand et al. (2021)**, geometry of gene regulatory dynamics
  (*PNAS* 118:e2109729118): the formal footing for reading CERPS's
  N-dimensional expression space.
* **Hoppe et al. (2016)** (*Nature* 535:299-302): PU.1/GATA1 protein
  ratios are NOT the trigger of myeloid fate choice; multiple generations
  are required - the key empirical motivation of the baseline paper (fate
  transitions driven by division-coupled epigenetic inheritance, not by
  fast expression fluctuations) which CERPS preserves.

## 4. Continuous, state-dependent cell kinetics

The baseline model assigns discrete kinetic rates per classified phenotype
(SC/TA1/TA2; Table 1 of the paper).  CERPS replaces the classification by
continuous functions of the cell state, following:

* **Lei (2020a)**, "A general mathematical framework for understanding
  the behaviour of heterogeneous stem cell regeneration",
  *J. Theor. Biol.* 492:110196; **Lei, Levin & Nie (2014)**
  (*PNAS* 111:E880-E887): phenotype-structured population models in which
  birth/death rates are functions of a continuous cellular state - the
  direct theoretical ancestor of CERPS's continuous-rate mode.
* **Lei (2020b)**, evolutionary dynamics of cancer from epigenetic
  regulation to cell population dynamics (*Sci. China Math.* 63:411-424):
  explicitly couples epigenetic state, rates and population - the same
  coupling structure CERPS adopts.
* **Bernard, Belair & Mackey (2003)** (cyclical neutropenia model;
  cited by the paper): Hill-function rate modulation derived from
  signalling/receptor arguments - justification for the Hill forms of
  beta(kappa, S) etc.
* **Stiehl & Marciniak-Czochra (2011, 2012)**; **Pazdziorek et al. (2014)**
  (*BMC Syst. Biol.*; mathematical models of stem-cell driven tissue
  organisation): continuous maturity-structured dynamics, and
  **Pedersen et al. (2023)** (*Math. Biosci. Eng.* / *Springer* review)
  for HSC-specific formulations.
* Sharp-limit consistency: CERPS's continuous rates are constructed so
  that steep Hill switches (p -> inf) recover Table 1 exactly; this is
  verified numerically in `tests/test_cerps.py`.

## 5. Spatial and mechanical population models

* **Osborne et al. (2017)**, "Comparing individual-based approaches to
  modelling the self-organisation of multicellular tissues"
  (*PLoS Comput. Biol.* 13:e1005387): centre-based (off-lattice) models -
  the class CERPS's space module implements.
* **Van Liedekerke et al. (2018)**, off-lattice agent-based models for
  cell and tumour growth (review); **Drasdo & Hoehme** and
  **Sadhukhan et al. (2021)** multiscale lattice-free tumour ABMs.
* **Germano et al. (2023)**, free and interfacial boundaries in
  individual-based models with contact-inhibited cycles (*PMC10560655*);
  **Zimmermann et al. (2016)** (*PNAS* 113:E61) contact inhibition of
  locomotion/proliferation in MDCK colonies - basis of CERPS's
  density-dependent beta factor.
* **Ghaffarizadeh et al. (2018)** PhysiCell and **Letort et al. (2019)**
  PhysiBoSS (*Bioinformatics* 35:1188): agent-based cells + substrate
  diffusion fields + intracellular (Boolean) regulation - the reference
  architecture for CERPS's optional paracrine field (FTCS
  diffusion-decay on a grid, cell-based secretion, state-dependent
  response).
* **Yu et al. (2020)** ARCADE (*Front. Bioeng. Biotechnol.* 8:249):
  multiscale ABM of heterogeneous populations with paracrine signalling.

## 6. Noise modelling

* The exponential-form extrinsic noise (multiplicative log-normal rate
  perturbations driven by OU processes) follows the baseline paper
  (Eqs. (2)-(4)), with justifications from **Austin et al. (2006)**
  (*Nature* 439:608) for log-normal expression rates and **Lei (2011)**
  for the chemical-potential argument.
* Stochastic partitioning noise (Sec. 2) provides an additional,
  division-coupled noise source absent from the baseline model.

## 7. Mapping to CERPS design decisions

| CERPS feature | Literature basis |
|---|---|
| continuous epigenetic dynamics (activity-dependent writing/erasing) | Dodd 2007; Zhao 2021; Probst 2009 |
| beta inheritance kernel (baseline) | Huang & Lei 2018/2019 |
| gaussian / deterministic kernels | generic phenomenological alternatives matched to the same conditional moments |
| explicit CpG-site kernel | Jones 2009; Moore 2013 |
| binomial expression partitioning | Huh & Paulsson 2011; Landgraf 2012 |
| N-gene GRN, multiplicative cross-repression | Jolly 2015; Wu 2017; Huang 2012 |
| continuous rates (stemness S, soft counts, senescence counter) | Lei 2014/2020a/2020b; Bernard 2003 |
| 3D off-lattice mechanics + contact inhibition | Osborne 2017; Zimmermann 2016; Germano 2023 |
| paracrine diffusion field | Ghaffarizadeh 2018; Letort 2019; Yu 2020 |
| exponential-form OU extrinsic noise | baseline paper; Austin 2006; Lei 2011 |

## References (baseline paper and upstream code)

* Huang, R., Situ, Q., Lei, J. (2024). Dynamics of cell-type transition
  mediated by epigenetic modifications. *J. Theor. Biol.* 577:111664.
  https://doi.org/10.1016/j.jtbi.2023.111664
* Original C++/Matlab source code (GPL-3.0):
  https://github.com/jinzhilei/Cell-type-transition
