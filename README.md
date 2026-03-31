# Ascertainment 

Framework for detecting and quantifying ascertainment bias in large biobanks using summary statistics.

---

## 📦 Data

The tutorial dataset used for illustration is available on Zenodo:

👉 https://doi.org/10.5281/zenodo.19340785

After downloading and extracting the dataset, ensure that the `data/` directory is placed directly in the root of this repository without modifying its internal structure.

The expected layout is:
```
repo/
├── data/
│   ├── beta/
│   ├── independent_snps/
│   ├── sample_with_census_rep/
│   └── sample_with_1KGP/
```

Scripts will fail if the directory structure is altered.

This dataset is provided for reproducibility and demonstration purposes.  
All underlying data sources are publicly available (see manuscript for details).

---

## 📦 Required R packages

Before running the scripts, install and load the required R packages.

### Install packages (run once)

```r
install.packages(c(
  "data.table",
  "ggplot2",
  "pheatmap",
  "reshape2",
  "nnls"
))
```

### Notes
These packages are required for data processing, visualization, and ancestry deconvolution.
Make sure all packages are installed without errors before running the scripts.

## 💻 Execution environment

All scripts are designed to run in a standard R environment on a local machine.

- No high-performance computing (HPC) resources are required for the tutorial dataset  
- The provided example can be run on a typical laptop or desktop  
- For larger datasets, users may optionally run the scripts on an HPC system  

The code does not rely on any HPC-specific configuration.

## 🚀 Quick start
```r
source("scripts/estimate_asc_with_census_rep.R")
source("scripts/get_ancestry_loadings.R")
source("scripts/estimate_asc_with_1KGP.R")

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

> ⚠️ **Note on reference data**
>
> The reference allele frequencies used for SweGen and TWB in this example are not true census-representative estimates.
>
> They are included for illustration and practical demonstration only.
>
> For real analyses, users should use appropriately matched, census-representative reference data where available.

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

- $\theta_1$: mean PGS difference estimator  
- $\theta_2$: regression-based estimator (more robust)  
- $I_2$: intercept capturing systematic deviation (similar to LD score intercept)
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

- $\theta_1$:, $\theta_2$, $I_2$
- standard errors
- p-values
- SNP counts

### 📊 Interpretation of estimators
- $\theta_1$: mean PGS difference estimator  
- $\theta_2$: regression-based estimator (more robust)  
- $I_2$: intercept capturing systematic deviation (similar to LD score intercept)
  
## 🧪 Simulations

Simulation scripts used in the manuscript are provided in:
```r
simulations/
```
These reproduce:
- statistical power
- impact of population stratification with drift and extreme model
- deconvolution approach
- impact of missing reference during deconvolution
- impact of attenuation bias in the reference during deconvolution


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
