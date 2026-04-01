# Ascertainment 

Framework for detecting and quantifying ascertainment bias in large biobanks using summary statistics.

---

## 📦 Data

The tutorial dataset used for illustration is available on Zenodo:

👉 https://doi.org/10.5281/zenodo.19340785

After downloading and extracting the dataset, place the data/ directory in the root of this repository without modifying its internal structure.

The expected layout is:
```
repo/
├── data/
│   ├── beta/
│   ├── independent_snps/
│   ├── sample_with_census_rep/
│   └── sample_with_1KGP/
```
Scripts assume this directory structure.

This dataset is provided for reproducibility and demonstration purposes.  
All underlying data sources are publicly available (see manuscript for details).

---

## 📦 Requirements

Install required R packages:

```r
install.packages(c(
  "data.table",
  "ggplot2",
  "pheatmap",
  "reshape2",
  "nnls"
))
```

## 💻 Execution environment

All scripts are designed to run in a standard R environment.

- No high-performance computing (HPC) resources are required for the tutorial dataset  
- The example runs on a typical laptop or desktop 
- Larger datasets may benefit from HPC environments  

## 🚀 Quick start

```r
source("scripts/estimate_asc_with_census_rep.R")
source("scripts/get_ancestry_loadings.R")
source("scripts/estimate_asc_with_1KGP.R")
```
---

## 📁 Repository structure
```
.
├── README.md
├── LICENSE
├── scripts/
│   ├── get_ancestry_loadings.R
│   ├── estimate_asc_with_census_rep.R
│   └── estimate_asc_with_1KGP.R
├── simulations/
│   ├── deconvolution_drift_model.R
│   ├── deconvolution_extreme_model.R
│   ├── impact_of_attenuation_bias_drift_model.R
│   ├── impact_of_attenuation_bias_extreme_model.R
│   ├── impact_of_missing_ref_drift_model.R
│   ├── impact_of_missing_ref_extreme_model.R
│   ├── population_stratification_drift_model.R
│   ├── population_stratification_extreme_model.R
│   └── power_theta.R
├── data/
│   ├── beta/
│   ├── independent_snps/
│   ├── sample_with_census_rep/
│   └── sample_with_1KGP/
├── result_cen_rep/
├── result_anc_loadings/
└── result_1kgp_rep/
```
> ⚠️ Output directories (`result_*`) are generated automatically at the repository root when scripts are executed.


## 📊 Overview

This framework provides two complementary approaches for detecting ascertainment bias:

### 1. Census-representative reference (preferred)

- Uses observed reference allele frequencies (`pr`)
- Minimizes bias from population stratification
- Recommended when representative reference data are available

> ⚠️ **Note on reference data**
>
> The reference allele frequencies used in this tutorial are not true census-representative estimates.
>
> They are included for illustration and practical demonstration only.
>
> For real analyses, appropriately matched census-representative reference data should be used.

### 2. Ancestry-deconvolution approach

- Estimates ancestry loadings from the 1000 Genomes Project (1KGP)
- Constructs an ancestry-matched synthetic reference
- Enables analysis when no census reference is available

---

## 🧬 Data conventions

All input files follow a consistent and harmonized format.

### Allele-frequency files

Required columns:

- `SNP`: SNP identifier  
- `CHR`: chromosome  
- `REF`: reference allele  
- `ALT`: alternate allele  
- Frequency columns: allele frequency of `ALT` for each cohort 

### Effect-size files

Required columns:

- `SNP`: SNP identifier  
- `A1`: effect allele  
- `A2`: non-effect allele  
- `BETA`: effect size corresponding to `A1`

#### File naming

- Beta file names must follow:

  `<Trait>_<Study>_<Year>`

Example:

  `Height_Yengo_et_al_2022`

- The trait name is inferred from the file name as the substring preceding the first underscore.
- Trait names should not contain underscores.

  For example:
  - ❌ `Adult_Height_Yengo...`
  - ✅ `AdultHeight_Yengo...`

- File names must match those in `data/independent_snps/` for correct SNP matching.

### Independent SNP sets

The files in `data/independent_snps/` contain approximately independent variants selected from GWAS summary statistics using LD clumping (PLINK):

- LD threshold: $r^2 < 0.01$  
- Window size: 1 Mb  
- Association threshold: $P < 5 \times 10^{-3}$  

This follows the specification described in the main manuscript.

- File names in `independent_snps/` must match those in `data/beta/`.

#### File format

Each file in `data/independent_snps/` contains a single column:

- `SNP`: list of independent SNP identifiers

SNP identifiers must match those in the corresponding beta and frequency files.
No additional columns are expected.
  
### Allele alignment

- Effect sizes (`BETA`) are aligned to the effect allele (`A1`)  
- Allele frequencies are expressed with respect to the `ALT` allele  
- The pipeline enforces:  
  **`ALT` ≡ `A1` (effect allele / trait-increasing allele)**  

This ensures that:

- `BETA` and allele frequencies refer to the same allele  


### Pre-processing (already applied in tutorial data)

The provided dataset has been pre-processed to ensure compatibility:

- Minor allele frequency (MAF) filtering  
- SNP matching across datasets  
- Harmonization of allele coding between frequency and effect-size files  

---

## ⚠️ Minimum SNP requirement

To produce valid estimates:

- At least **10 SNPs** must remain after filtering for independent snps
- Reported as `M` in output  

---

## 1️⃣ Census-based approach

- Use this when you have a census-representative reference.

```r
source("scripts/estimate_asc_with_census_rep.R")
```

### Required inputs:
- Allele-frequency columns:
 
Each cohort must provide:

  - `[cohort]_ps`: study allele frequency  
  - `[cohort]_pr`: reference allele frequency  

   **Example (SweGen):**
   
    `SNP`, `CHR`, `REF`, `ALT`, `SweGen_ps`, `SweGen_pr`
    
- Cohort metadata
```r

cohort_info <- data.frame(
  cohort = c("SweGen", "TWB"),
  Ns = c(1000, 92615),
  Nr = c(1000, 100000),
  stringsAsFactors = FALSE
)
```
 - `Ns`: study sample size
 - `Nr`: reference sample size


- LD reference panel
```r

ld_ref_info <- data.frame(
  cohort = c("SweGen", "TWB"),
  ld_ref = c("EUR", "EAS"),
  stringsAsFactors = FALSE
)

```
> ⚠️ Note: LD reference panels used here are predefined.
> Guidance on selecting appropriate LD panels is provided in the ancestry-deconvolution approach.

#### Output

For each trait–cohort pair:

- $\theta_1$: mean PGS difference estimator  
- $\theta_2$: regression-based estimator (more robust)  
- $I_2$: intercept capturing systematic deviation (similar to LD score intercept)
- standard errors and p-values
- SNP counts(`M`, `M_prior`)



## 2️⃣ Ancestry-deconvolution approach

Use when no census reference is available.

### Step 1: Estimate ancestry loadings
```r
source("scripts/get_ancestry_loadings.R")
```

#### Input
- Allele-frequency must contain:
  - Cohort allele frequency (e.g. `SweGen`)
  - 26 populations allele frequencies from 1KGP (`ACB` , `ASW` , ..., `YRI`)

     **Example (SweGen):**

    `SNP`, `CHR`, `REF`, `ALT`, `SweGen`, `ACB` , `ASW` , ...,`YRI` 

- Cohort 
```r
cohort_names <- c("SweGen", "TWB")
```
> ⚠️ Note: ancestry weights are estimated using ~100,000 randomly sampled SNPs.

#### Outputs (in result_anc_loadings/)
Files:
  - `ancestry_population_weights.tsv`: population-level weights (1KGP)
  - `ancestry_superpop_weights.tsv` : superpopulation weights (`AFR`, `AMR` , `EAS` , `EUR` , `SAS`)
  - ancestry_dominant.tsv: dominant ancestry per cohort
  - ancestry_ld_panel_info.tsv: selected LD reference panel
    
> ⚠️ Note: Dominant ancestry defined as ≥90% contribution.
> Mixed-ancestry cohorts that do not meet the ≥90% dominance threshold may require combined LD reference panels.

Plots:
  - population heatmap
  - super-population bar plot
  - dominant ancestry plot

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
- allele frequency files (same as step 1)
- beta files
- independent snps files

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

To apply this framework to external datasets:

- Ensure input files follow the required structure (see *Data conventions*)
- Align alleles consistently across all files:
  - Effect sizes (`BETA`) must correspond to the same allele as allele frequencies (ALT)
- Verify that allele frequencies and effect sizes are harmonized (no strand or allele mismatches)
- Provide accurate sample sizes:
  - `Ns`: study sample size  
  - `Nr`: reference sample size (if applicable)
- Use consistent cohort naming across all inputs (frequency files, metadata, LD panel mapping)
- Ensure SNP overlap between:
  - effect-size files  
  - allele-frequency files  
  - independent SNP sets


## 📄 Citation

If you use this repository or dataset, please cite:

Using summary data to detect and quantify ascertainment in biobanks


