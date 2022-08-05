# Childhood exposure to non-persistent endocrine disrupting chemicals and multi-omic markers: a repeated-measures panel study

Repository containing the code to reproduce the results of the first paper of my PhD.
DOI: TODO.

## Abstract
*Background*: Individuals are exposed to multiple environmental pollutants with endocrine disrupting activity (endocrine disruptors, EDCs) and the early stages of life are recognized to be particularly susceptible to these stressors. Previous studies have focused on identifying molecular signatures associated with EDCs to elucidate potential mechanisms, but few have integrated multiple omics layers. We aimed to identify multi-omic signatures associated with childhood exposure to non-persistent EDCs using an integrative network approach.

*Methods*: We used data from the HELIX (Human Early-Life Exposome) Child Panel Study, which included 152 children aged 6 to 11. Children were followed for 2 weeks 6 months apart. Twenty-two non-persistent EDCs (10 phthalate metabolites, 7 phenols, and 5 organophosphate pesticide metabolites) were measured in the 2 weekly pools of 15 urine samples each, to account for their intra-individual variability. Multi-omic profiles (methylome, serum and urinary metabolome, proteome) were measured in blood and in a pool of 2 urine samples collected at the end of each week. We developed Gaussian Graphical Models (GGMs) based on shrinkage estimates of pairwise partial correlations between EDCs and covariate-adjusted molecular features in each week. The obtained networks for the two visits were merged in order to identify reproducible associations.

*Results*: The visit-specific networks included 4,083 and 4,908 edges of comparable strength and statistical significance. Most of the reproducible associations across visits included features of the same layer. Among the mixed associations, those between exposures and either proteins or urinary metabolites showed the highest reproducibility scores. The reproducible inter-layer associations were corroborated with data from the literature. Among the most significant, we found associations between the organophosphate pesticide diethyl phosphate and serotonin, triclosan and serotonin, mono‑4‑methyl‑7‑hydroxyoctyl phthalate and kynurenine, and triclosan and leptin.

*Conclusions*: Through partial correlation network analyses, repeat time points and weekly EDC exposure assessment, we identified several biologically relevant molecular signatures related to non-persistent EDC exposure in childhood.

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
