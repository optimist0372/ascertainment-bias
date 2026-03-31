# Ascertainment-project

Framework for detecting and quantifying ascertainment bias in large biobanks using summary statistics.

---

## 📦 Data

The tutorial dataset used for illustration is available on Zenodo:

👉 https://doi.org/10.5281/zenodo.19340785

After downloading the dataset from Zenodo, extract it and copy all files into the `data/` folder of this repository, ensuring that the directory structure remains unchanged.

This dataset is provided for reproducibility and demonstration purposes.  
All underlying data sources are publicly available (see manuscript for details).

---

## 🚀 Quick start
```r
source("scripts/get_ancestry_loadings.R")
source("scripts/estimate_asc_with_1KGP.R")
source("scripts/estimate_asc_with_census_rep.R")
---

## 📁 Repository structure
.
├── README.md
├── LICENSE
├── scripts/
│   ├── get_ancestry_loadings.R
│   ├── estimate_asc_with_census_rep.R
│   └── estimate_asc_with_1KGP.R
├── simulations/
│   ├── Deconvolution_drift_model.R
│   ├── Deconvolution_extreme_model.R
│   ├── Impact_of_attenuation_bias_drift_model.R
│   ├── Impact_of_attenuation_bias_extreme_model.R
│   ├── Impact_of_missing_ref_drift_model.R
│   ├── Impact_of_missing_ref_extreme_model.R
│   ├── Population_stratification_drift_model.R
│   ├── Population_stratification_extreme_model.R
│   └── Power_theta.R
├── data/ ##get this folder using the Zenodo link
│   ├── beta/
│   ├── independent_snps/
│   ├── sample_with_census_rep/
│   └── sample_with_1KGP/
└── results/
    ├── result_anc_loadings/
    ├── result_cen_rep/
    └── result_1kgp_rep/

```
## 📊 Overview

This framework provides two complementary approaches for detecting ascertainment bias:

### 1. Census-representative reference (preferred)

- Uses observed reference allele frequencies (`pr`)
- Minimizes bias from population stratification
- Recommended when representative reference data are available

### 2. Ancestry-deconvolution approach

- Estimates ancestry loadings from the 1000 Genomes Project (1KGP)
- Constructs an ancestry-matched synthetic reference
- Enables analysis when no census reference is available

---

## 🧬 Data conventions

All input files follow consistent formatting.

### Allele-frequency files

Must contain:

- `SNP`: SNP identifier  
- `CHR`: chromosome  
- `REF`: reference allele  
- `ALT`: alternate allele  

### Effect-size files

- Effect sizes (`BETA`) are aligned to `ALT`

### Alignment assumption

- `ALT` corresponds to the trait-increasing allele  
- Allele frequencies are aligned to the same allele as effect sizes  

### Pre-processing (already applied in tutorial data)

- MAF filtering  
- SNP matching  
- Harmonization of allele coding  

---

## ⚠️ Minimum SNP requirement

To produce valid estimates:

- At least **10 SNPs** must remain after filtering  
- Reported as `M` in output  

---

## 1️⃣ Census-based approach

- Use this when you have a census-representative reference.

```r
source("scripts/estimate_asc_with_census_rep.R")
```

### Required inputs:
#### Allele-frequency columns
Each cohort must have:

- [cohort]_ps: study sample frequency
- [cohort]_pr: reference frequency

Example:

- SweGen_ps, SweGen_pr
- TWB_ps, TWB_pr

#### Cohort metadata
```r

cohort_info <- data.frame(
  cohort = c("SweGen", "TWB"),
  Ns = c(1000, 92615),
  Nr = c(1000, 100000),
  stringsAsFactors = FALSE
)
```

- Ns: study sample size
- Nr: reference sample size

#### LD reference panel
```r

ld_ref_info <- data.frame(
  cohort = c("SweGen", "TWB"),
  ld_ref = c("EUR", "EAS"),
  stringsAsFactors = FALSE
)
```

#### Output

For each trait–cohort pair:

- theta1: standardized PGS difference
- theta2: regression-based estimator
- I2: intercept (diagnostic)
- standard errors and p-values
- M, M_prior



## 2️⃣ Ancestry-deconvolution approach

Use when no census reference is available.

### Step 1: Estimate ancestry loadings
```r
source("scripts/get_ancestry_loadings.R")
```

#### Input
- cohort allele frequencies
- 26 ancestral population allele frequencies from 1KGP
  
#### Outputs (in result_anc_loadings/)
- ancestry_population_weights.tsv
- ancestry_superpop_weights.tsv
- ancestry_dominant.tsv
- ancestry_ld_panel_info.tsv

#### Plots
- population heatmap
- super-population bar plot
- dominant ancestry plot

#### Interpretation
- Dominant ancestry defined as ≥90% contribution
- Mixed ancestry cohorts use combined LD panels

### Step 2: Estimate ascertainment (1KGP reference)
```r
source("scripts/estimate_asc_with_1KGP.R")
```
#### Required inputs
```r
cohort_info <- data.frame(
  cohort = c("SweGen", "TWB"),
  Ns = c(1000, 92615),
  stringsAsFactors = FALSE
)
```
#### Additional inputs
- ancestry weights (Step 1 output: ancestry_population_weights.tsv)
- LD panel mapping (Step 1 output: ancestry_ld_panel_info.tsv)
- beta files
- allele-frequency files

#### Output
Same structure as census-based approach:

- theta1, theta2, I2
- standard errors
- p-values
- SNP counts

### 📊 Interpretation of estimators
- theta1: mean PGS difference estimator
- theta2: regression-based estimator (more robust)
- I2: captures systematic deviation not explained by model

  
## 🧪 Simulations

Simulation scripts used in the manuscript are provided in:
```r
simulations/
```
These reproduce:
- statistical power
- drift models
- extreme scenarios
- missing reference effects
- attenuation bias


## 🛠 Using your own data

To apply this framework:

- match input file structure
- align alleles consistently
- ensure effect sizes correspond to trait-increasing allele frequency 
- provide correct sample sizes
- use consistent cohort naming

## 📄 Citation

If you use this repository or dataset, please cite:

Using summary data to detect and quantify ascertainment in biobanks

Zenodo dataset:
https://doi.org/10.5281/zenodo.19340785
