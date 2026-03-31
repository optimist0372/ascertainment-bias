# Ascertainment-project

Framework for detecting and quantifying ascertainment bias in large biobanks using summary statistics.

---

## 📦 Data

The tutorial dataset used for illustration is available on Zenodo:

👉 https://doi.org/10.5281/zenodo.19340785

This dataset is provided for reproducibility and demonstration purposes.  
All underlying data sources are publicly available (see manuscript for details).

---

## 🚀 Quick start

```r
source("scripts/get_ancestry_loadings.R")
source("scripts/estimate_asc_with_1KGP.R")
source("scripts/estimate_asc_with_census_rep.R")

📁 Repository structure

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
├── data/
│   ├── beta/
│   ├── independent_snps/
│   ├── sample_with_census_rep/
│   └── sample_with_1KGP/
└── results/
    ├── result_anc_loadings/
    ├── result_cen_rep/
    └── result_1kgp/


## 📊 Overview

This framework provides two complementary approaches for detecting ascertainment bias:

1. Census-representative reference (preferred)
  * Uses observed census represenatative reference allele frequencies (pr)
  - Minimizes bias from population stratification
  + Recommended when representative reference data are available

2. Ancestry-deconvolution approach
  - Estimates ancestry loadings from 1000 Genomes Project (1KGP)
  - Constructs an ancestry-matched synthetic reference
  - Enables analysis when no census reference is available


🧬 Data conventions

All input files follow consistent formatting.

Allele-frequency files

Must contain:

- `SNP`: SNP identifier
- `CHR`: chromosome
- `REF`: reference allele
- `ALT`: alternate allele
- `Effect-size files
  Effect sizes (BETA) are aligned to ALT
  Alignment assumption
  ALT corresponds to the trait-increasing allele
  Allele frequencies are aligned to the same allele as effect sizes


Pre-processing (already applied in tutorial data):
MAF filtering
SNP matching
Harmonization of allele coding

⚠️ Minimum SNP requirement

To produce valid estimates:

At least 10 SNPs must remain after filtering
Reported as M in output


1️⃣ Census-based approach

Use this when you have a census-representative reference.

source("scripts/estimate_asc_with_census_rep.R")

Required inputs
Allele-frequency columns

Each cohort must have:

[cohort]_ps: study sample frequency
[cohort]_pr: reference frequency

Example:

SweGen_ps, SweGen_pr
TWB_ps, TWB_pr
