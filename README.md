# Power PheWAS 
This pipeline was used in the following publication:

Moysés-Oliveira et al 2025. Pleiotropic effects of APOE variants on a sleep-based adult epidemiological cohort. Sleep Medicine. https://doi.org/10.1016/j.sleep.2025.106490

# Usage

This Github contains the workflow pipeline that was used to run power analysis with simulation to select traits for PheWAS. The workflow is organized as a [snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline. 

To run the pipeline with snakemake:
```console
cd  PowerPheWAS/
mamba activate snakemake
snakemake --cores 7 --resources mem_mb=28672
# number of cores can be modified accordingly
```

# Input data
- N per trait: 
.txt file with three columns phenotype name (Pheno), sample size (N) and PheWAS method (method)

- Unique N: 
.txt file with two columns: sample size (N) and PheWAS method (method). This file can be generated from the previous, retaining unique combinations of N and method

For sample size, total N for continuous traits and N per group split by ; for categorical data (e.g 797;23;21;40). Method can be one of: LinearRegression, LogisticRegression, OrdinalRegression

# Config file
*config_power.yaml* with input data and simulation parameters:
- alpha, number of simulations, MAF, effect size and power for filtering
- Number of runs per regression type: maf * alpha * effect * trait

# How is power calculated?
The pipeline is designed to calculate power by simulating data under the alternative hypothesis (true effect of SNV on trait). Linear, logistic, or ordinal logistic regression models are applied to test the simulated SNV effect. This process is repeated n_sim times. For each simulation, the p-value is retained, and power is estimated as the proportion of simulations yielding a p-value below the specified alpha. 

An additive model is tested for the SNV. 

Unordered categorical with more than 2 categories can be tested with logistic regression with groups being considered in a pairwise manner.

For simulation with ordinal model, simulations that resulted in less than 5 individuals with a genotype or that failed proportional odds assumptions were discarted and a new simulations was done until the total of n_sim is achieved.

# Tools
- R 
