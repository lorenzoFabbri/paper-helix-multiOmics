- comparison_methods.R
	Functions to fit different linear and non-linear models and compare them. Abandoned project
- comparison_methods.Rmd
	Notebook to run functions from comparison_methods.R. Abandoned project
- diagnose.Rmd
	Notebook to analyze variables in Exposome and -omics datasets and eventually identify "problematic" variables (e.g., non-normality)
- exploratory_dataAnalysis.R
	Functions to perform EDA on both Exposome and -omics datasets (e.g., distribution variables, heatmaps, correlations, PCA)
- exploratory_dataAnalysis.Rmd
	Notebook to run functions from exploratory_dataAnalysis.R. It also serves as notebook to create the datasets for both time points
- hierarchical_communityNetwork.R
	Functions to perform enrichment analysis on communities of Exposures and -omics variables. Abandoned project
- hierarchical_communityNetwork.Rmd
	Notebook to run functions from hierarchical_communityNetwork.R. Abandoned project
- nultivariate_analysis.Rmd
	Main notebook to fit (s)PLS models. Not necessary anymore but it contains some notes on the interpretation of the outputs from mixOmics
- prepare_datasets.R
	Functions to load and prepare all datasets used in downstream analyses
- regularized_regression.Rmd
	Notebook to fit regularized regression models (e.g., GLMNET). Abandoned project
- deprecated/
	- multivariate_analysis.R
		Functions to fit (s)PLS models from mixOmics. Deprecated
	- run_main.R
		Functions to actually fit (s)PLS models from mixOmics. Deprecated
- ggm/
	- ggm.R
		Functions to generate partial correlation networks (i.e., GGMs) for Exposures and -omics datasets. This is the main script containing all the functions related to the GGMs part of the paper
- health/
	- health.R
		Functions to interrogate databases of interest to automatically retrieve the health effects of the considered chemicals (i.e., exposures)
- multivariate_analysis/
	- compare_time.R
		Functions to compare results (s)PLS for the 2 time points
	- dictionaries.R
		Functions to retrieve all the variables in all the datasets stratified by e.g., chemical group or -omic layer. Functions to map different names to different codes/IDs
	- main.R
		Main script containing functions to actually perform all the analyses related to the PLS part of the paper. Functions to fit models for all the combinations of parameters of interest
	- model.R
		Functions that essentially call mixOmics functions to fit (s)PLS models. Helper functions to perform adjustment for covariates
	- networks.R
		Functions to generate Relevance Networks from PLS models and to post-process and compare them
	- plot.R
		Functions for custom plots from mixOmics fitted (s)PLS models using ggplot2
	- preprocess.R
		Functions to scale, transform, pre-process datasets prior to fitting (s)PLS models
	- sgpls_preliminary.R
		Functions to fit sgPLS models. Abandoned project since the sgPLS R package does not work anymore
	- tables.R
		Functions to generate publication-ready, interactive tables of results (e.g., scores of PLS models) and to arrange the plots to compare the time points
	- validate.R
		Functions to perform validation of fitted (s)PLS models. Abandoned project or at least not in development as of 11/Oct/21

