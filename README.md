
# **A repository for the multi-state diel occupancy model (MSDOM) and associated manuscripts:** 

### River, K, Fidino, M, Farris, ZJ, Murphy, A, Magle, SB, and Gerber, BD. Rethinking habitat occupancy modeling and the role of diel activity in an anthropogenic world. In Review.
---


**This repository includes 6 subfolders: Chicago coyote, Data Processing, JAGS, Makira Fosa2, RNP Fosa, and Simulation Files**
1) The **Chicago coyote** folder includes data, R scripts, and plots specific to case study on coyotes that uses the dynamic MSDOM.
2) The **Data Procesing** folder includes R scripts and example data on how to prepare data for the MSDOM.
3) The **JAGS** folder includes JAGS models for the static and dynamic MSDOM, including the full, reduced, and null parameterizations.
4) The **Makira Fosa2** folder includes fosa data from Makira Natural Park and R scripts for fitting the static MSDOM, including model comparison using CPO.
5) The **RNP Fosa** folder includes fosa data from Ranomafana National Park and R scripts for fitting the MSDOM, including model comparison using CPO.
6) The **Simulation Files** folder includes scripts for simulating data under different versions of the MSDOM (full, reduced, null) and fitting these models using JAGS.

---

<div align="center"><img width="150" height="auto" src="raccoon.jpg" alt="A silhouette of a raccoon." /></div>

<div align="center"> <h2>JAGS folder</h2> </div>
<div align="left"> <h3>Dynamic MSDOM File</h3> </div>
**jags.dynamic.fake.multistate.R**

**jags.dynamic.multistate.covars.lasso.R**

**jags.dynamic.multistate.covars.R**

**jags.dynamic.multistate.null.R**

**temporal_multi_varying_covars.R**

**null_temporal_multi_varying_covars.R**

<div align="left"> <h3>Static MSDOM</h3> </div>
<div align="center"> <h4>FULL Model</h3> </div>
**jags.multistate.occ.full.R** - Full MSDOM with probabilities estimated directly (no logit transformation) without site covariates.

**jags.multistate.occ.full.alt.R**- Full MSDOM with probabilities estimated on the logit scale without site covariates. 

**jags.multistate.occ.full.alt.RE.R** - Full MSDOM with probabilities estimated on the logit scale; detection and occupancy parameters are realizations from a higher order Normal distribution (i.e. random effect). No covariates.

**jags.multistate.occ.full.site.covs.by.state.R** - Full MSDOM with probabilities estimated on the logit scale, in which occupancy parameters are modeled using any set or size of site covariates, separately for each state.

**jags.multistate.occ.full.site.covs.R** - Full MSDOM with probabilities estimated on the logit scale, in which occupancy parameters are modeled using the same single site covariate, separately for each state.

**jags.multistate.occ.full.site.covs.RE.R** - Full MSDOM with probabilities estimated on the logit scale, in which occupancy parameters are modeled using separate covariates for each state and detection and occupancy parameters are realizations from a higher order Normal distribution (i.e. random effect).   

<div align="center"> <h4>REDUCED MSDOM</h3> </div>
**jags.multistate.occ.reduced.R**

**jags.multistate.occ.reduced.alt.R**

**jags.multistate.occ.red.det.null.R**

**jags.multistate.occ.reduced.site.covs.R**     

**jags.multistate.occ.reduced.alt.RE.R**

**jags.multistate.occ.reduced.site.covs.RE.R**

**jags.multistate.occ.reduced.site.covs.by.state.R**

**jags.multistate.occ.null.R**

**jags.multistate.occ.null.alt.R**

**jags.multistate.occ.null.det.null.R**

**jags.multistate.occ.null.alt.RE.R**

**jags.multistate.occ.null.site.covs.R**

**jags.multistate.occ.null.site.covs.by.state.R**

**jags.multistate.occ.null.site.covs.RE.R**




---


**city_mean_occupancy.R:** This is a single-species hierarchical occupancy model with among-city effects.

---

<div align="center"><img width="150" height="auto" src="coyote.jpg" alt="A silhouette of a coyote." /></div>

<div align="center"> <h3>Data</h3> </div>
