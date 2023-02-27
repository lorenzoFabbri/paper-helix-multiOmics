# Childhood exposure to non-persistent endocrine disrupting chemicals and multi-omic profiles: a panel study

Repository containing the code to reproduce the results of the first paper of my PhD.

DOI: [https://doi.org/10.1016/j.envint.2023.107856](https://doi.org/10.1016/j.envint.2023.107856).

## Abstract
### Background
Individuals are exposed to environmental pollutants with endocrine disrupting activity (endocrine disruptors, EDCs) and the early stages of life are particularly susceptible to these exposures. Previous studies have focused on identifying molecular signatures associated with EDCs, but none have used repeated sampling strategy and integrated multiple omics. We aimed to identify multi-omic signatures associated with childhood exposure to non-persistent EDCs.

### Methods
We used data from the HELIX Child Panel Study, which included 156 children aged 6 to 11. Children were followed for one week, in two time periods. Twenty-two non-persistent EDCs (10 phthalate, 7 phenol, and 5 organophosphate pesticide metabolites) were measured in two weekly pools of 15 urine samples each. Multi-omic profiles (methylome, serum and urinary metabolome, proteome) were measured in blood and in a pool urine samples. We developed visit-specific Gaussian Graphical Models based on pairwise partial correlations. The visit-specific networks were then merged to identify reproducible associations. Independent biological evidence was systematically sought to confirm some of these associations and assess their potential health implications.

### Results
950 reproducible associations were found among which 23 were direct associations between EDCs and omics. For 9 of them, we were able to find corroborating evidence from previous literature: DEP - serotonin, OXBE - cg27466129, OXBE - dimethylamine, triclosan - leptin, triclosan - serotonin, MBzP - Neu5AC, MEHP - cg20080548, oh-MiNP - kynurenine, oxo-MiNP - 5-oxoproline. We used these associations to explore possible mechanisms between EDCs and health outcomes, and found links to health outcomes for 3 analytes: serotonin and kynurenine in relation to neuro-behavioural development, and leptin in relation to obesity and insulin resistance.

### Conclusions
This multi-omics network analysis at two time points identified biologically relevant molecular signatures related to non-persistent EDC exposure in childhood, suggesting pathways related to neurological and metabolic outcomes.

## Structure of the repository
* `code/`
  * `ggm/`
    * `ggm.R`: main functions to fit GGMs to data
      * `main.pipeline.ggm()`: sets parameters (e.g., scaling method) and calls fitting functions (`perform.analysis()` and `fit.ggm()`)
      * `perform.analysis()`: loads data and covariates taking care of right visit and calls main fitting function (`fit.ggm()`)
      * `fit.ggm()`: adjusts for covariates, transforms data, computes correlations based on chosen method (e.g., `corpcor`)
      * `process.ggms()`: computes p-values using `GeneNet` and filters networks based on chosen method (e.g., probability)
      * `net.properties()`: computes simple properties of the fitted GGMs
      * `merge.networks()`: merges visit-specific networks to obtain the *merged network*
  * `hpc/`
    * `scripts_hpc.R`: driver script to fit GGMs to data and perform bootstrapping
  * `methylome/`: functions to pre-process methylation data and perform EWAS
    * `main.R`: script to load results EWAS and find common CpG sites across visits
    * `methylome.R`: main functions to run EWAS using `ewaff` (the main function being `methylome.ewaff()`)
    * `methylomeRaw.R`: script to pre-process methylation data with `minfi` and write to disk
    * `run_methylome1.R`, `run_methylome2.R`: driver scripts to run EWAS by visit
  * `multivariate_analysis/`
    * `dictionaries.R`: helper functions and dictionaries (e.g., mapping to chemical classes to color for plots, list of confounders)
    * `model.R`: helper functions to perform 2-stage residual-outcome regression to adjust omics for confounders
    * `preprocess.R`: helper functions to transform the data (scaling of exposures and omics)
  * `paper/`: functions to reproduce figures (`images.R`) and tables (`tables.R`) from the manuscript
    * `main.R` is the driver script to generate them
  * scripts
    * `diagnose.Rmd`: initial pre-processing of exposure and omics data by visit (tidying column names, imputation of missing exposures); creation of metadata (i.e., covariates); saving of datasets to disk
* `results/`
  * `final_material_paper_v2/`: material (i.e., figures and tables) that appear in the manuscript.
