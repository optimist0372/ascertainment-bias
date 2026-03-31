## set working directory ##

#setwd("/Users/uqbolase/Documents/github/")

## load required packages ###
pkgs <- c("data.table", "dplyr", "nnls", "quadprog")
invisible(lapply(pkgs, library, character.only = TRUE))

# ============================================================
# User inputs
# ============================================================

beta_dir <- "data/beta"
freq_dir <- "data/sample_with_census_rep"
ind_snp_dir <- "data/independent_snps"
output_dir <- "result_cen_rep"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_file <- file.path(output_dir, "ascertainment_results.tsv")

cohort_info <- data.frame(
  cohort = c("SweGen", "TWB"),
  Ns = c(1000, 92615),
  Nr = c(1000, 100000),
  stringsAsFactors = FALSE
)

ld_ref_info <- data.frame(
  cohort = c("SweGen", "TWB"),
  ld_ref = c("EUR", "EAS"),
  stringsAsFactors = FALSE
)

min_snp_count <- 10L

# ============================================================
# Helpers
# ============================================================

get_trait <- function(path) {
  strsplit(basename(path), "_")[[1]][1]
}

get_cohort_from_freq <- function(dt) {
  ps_cols <- grep("_ps$", names(dt), value = TRUE)
  
  if (length(ps_cols) != 1L) {
    stop(
      "Expected exactly one *_ps column in allele-frequency file; found: ",
      paste(ps_cols, collapse = ", ")
    )
  }
  
  cohort <- sub("_ps$", "", ps_cols)
  
  if (!paste0(cohort, "_pr") %in% names(dt)) {
    stop("Missing matching _pr column for cohort: ", cohort)
  }
  
  cohort
}

get_sample_sizes <- function(cohort, cohort_info) {
  idx <- match(cohort, cohort_info$cohort)
  
  if (is.na(idx)) {
    stop("Missing Ns/Nr for cohort: ", cohort)
  }
  
  list(
    Ns = cohort_info$Ns[idx],
    Nr = cohort_info$Nr[idx]
  )
}

get_ld_ref <- function(cohort, ld_ref_info) {
  idx <- match(cohort, ld_ref_info$cohort)
  
  if (is.na(idx)) {
    stop("Missing LD reference mapping for cohort: ", cohort)
  }
  
  ld_ref_info$ld_ref[idx]
}

get_ld_snp_files <- function(ind_snp_dir, ld_ref) {
  ld_path <- file.path(ind_snp_dir, ld_ref)
  
  if (!dir.exists(ld_path)) {
    stop("LD reference directory not found: ", ld_path)
  }
  
  list.files(ld_path, full.names = TRUE)
}

prepare_analysis_data <- function(freq_dt, beta_dt, ind_snp_dt) {
  if (!"SNP" %in% names(beta_dt)) {
    stop("beta file is missing SNP column")
  }
  
  if (!"SNP" %in% names(ind_snp_dt)) {
    stop("independent SNP file is missing SNP column")
  }
  
  beta_sub <- beta_dt[SNP %in% ind_snp_dt$SNP]
  
  analysis_dt <- merge(freq_dt, beta_sub, by = "SNP", all = FALSE)
  
  if (!all(c("ALT", "A1", "BETA") %in% names(analysis_dt))) {
    stop("Merged data are missing one or more of: ALT, A1, BETA")
  }
  
  analysis_dt[, BETA := as.numeric(BETA)]
  analysis_dt[, beta_aligned := fifelse(ALT == A1, BETA, -BETA)]
  
  analysis_dt
}

estimate_ascertainment_with_census_rep <- function(
    dt,
    cohort,
    Ns,
    Nr,
    min_snp_count = 10L
) {
  ps_col <- paste0(cohort, "_ps")
  pr_col <- paste0(cohort, "_pr")
  
  keep <- !is.na(dt$beta_aligned) &
    !is.na(dt[[ps_col]]) &
    !is.na(dt[[pr_col]])
  
  dt <- dt[keep]
  
  M <- nrow(dt)
  
  if (M < min_snp_count) {
    return(list(
      M = M,
      theta1 = NA_real_,
      se_theta1 = NA_real_,
      pval_theta1 = NA_real_,
      theta2 = NA_real_,
      se_theta2 = NA_real_,
      pval_theta2 = NA_real_,
      I2 = NA_real_,
      se_I2 = NA_real_,
      pval_I2 = NA_real_
    ))
  }
  
  b <- dt$beta_aligned
  p_s <- dt[[ps_col]]
  p_r <- dt[[pr_col]]
  
  u <- sum(2 * b * (p_s - p_r))
  v <- sum(2 * b^2 * p_r * (1 - p_r))
  
  if (!is.finite(v) || v <= 0) {
    return(list(
      M = M,
      theta1 = NA_real_,
      se_theta1 = NA_real_,
      pval_theta1 = NA_real_,
      theta2 = NA_real_,
      se_theta2 = NA_real_,
      pval_theta2 = NA_real_,
      I2 = NA_real_,
      se_I2 = NA_real_,
      pval_I2 = NA_real_
    ))
  }
  
  theta1 <- u / sqrt(v)
  
  var_ps <- p_s * (1 - p_s) / (2 * Ns)
  var_pr <- p_r * (1 - p_r) / (2 * Nr)
  
  var_u <- sum(4 * b^2 * (var_ps + var_pr))
  var_v <- sum(4 * b^4 * (var_pr * (1 - 2 * p_r)^2 + 2 * var_pr^2))
  cov_uv <- sum(-4 * b^3 * (1 - 2 * p_r) * var_pr)
  
  var_theta1 <- (1 / v) * (
    var_u +
      (u^2 * var_v) / (4 * v^2) -
      (u * cov_uv / v)
  )
  
  se_theta1 <- if (is.finite(var_theta1) && var_theta1 >= 0) {
    sqrt(var_theta1)
  } else {
    NA_real_
  }
  
  pval_theta1 <- if (!is.na(se_theta1) && se_theta1 > 0) {
    pchisq((theta1 / se_theta1)^2, df = 1, lower.tail = FALSE)
  } else {
    NA_real_
  }
  
  lhs <- p_s - p_r
  rhs <- b * p_r * (1 - p_r) / sqrt(v)
  
  ok <- is.finite(lhs) & is.finite(rhs) & is.finite(var_pr) & var_pr > 0
  
  theta2 <- se_theta2 <- pval_theta2 <- NA_real_
  I2 <- se_I2 <- pval_I2 <- NA_real_
  
  if (sum(ok) >= 3L) {
    fit <- lm(lhs[ok] ~ rhs[ok], weights = 1 / var_pr[ok])
    coef_tab <- summary(fit)$coefficients
    
    if (nrow(coef_tab) >= 2L) {
      theta2 <- coef_tab[2, 1]
      se_theta2 <- coef_tab[2, 2]
      pval_theta2 <- coef_tab[2, 4]
    }
    
    I2 <- coef_tab[1, 1]
    se_I2 <- coef_tab[1, 2]
    pval_I2 <- coef_tab[1, 4]
  }
  
  list(
    M = M,
    theta1 = theta1,
    se_theta1 = se_theta1,
    pval_theta1 = pval_theta1,
    theta2 = theta2,
    se_theta2 = se_theta2,
    pval_theta2 = pval_theta2,
    I2 = I2,
    se_I2 = se_I2,
    pval_I2 = pval_I2
  )
}

# ============================================================
# File discovery
# ============================================================

beta_files <- list.files(beta_dir, full.names = TRUE)
freq_files <- list.files(freq_dir, full.names = TRUE)

if (length(beta_files) == 0L) {
  stop("No beta files found in: ", beta_dir)
}

if (length(freq_files) == 0L) {
  stop("No allele-frequency files found in: ", freq_dir)
}

results <- vector("list", 0L)

# ============================================================
# Main loop
# ============================================================

for (freq_file in freq_files) {
  message("Reading frequency file: ", basename(freq_file))
  
  freq_dt <- fread(freq_file)
  cohort <- get_cohort_from_freq(freq_dt)
  
  sample_sizes <- get_sample_sizes(cohort, cohort_info)
  ld_ref <- get_ld_ref(cohort, ld_ref_info)
  
  message("Cohort: ", cohort, " | LD reference: ", ld_ref)
  
  ld_snp_files <- get_ld_snp_files(ind_snp_dir, ld_ref)
  ld_snp_lookup <- setNames(ld_snp_files, basename(ld_snp_files))
  
  for (beta_file in beta_files) {
    file_id <- basename(beta_file)
    
    if (!file_id %in% names(ld_snp_lookup)) {
      message("  Skipping (no LD SNP match): ", file_id)
      next
    }
    
    snp_file <- ld_snp_lookup[[file_id]]
    trait <- get_trait(beta_file)
    
    message("  Processing trait: ", trait)
    
    beta_dt <- fread(beta_file)
    ind_snp_dt <- fread(snp_file)
    
    analysis_dt <- prepare_analysis_data(
      freq_dt = freq_dt,
      beta_dt = beta_dt,
      ind_snp_dt = ind_snp_dt
    )
    
    est <- estimate_ascertainment_with_census_rep(
      dt = analysis_dt,
      cohort = cohort,
      Ns = sample_sizes$Ns,
      Nr = sample_sizes$Nr,
      min_snp_count = min_snp_count
    )
    
    results[[paste(trait, cohort, sep = "_")]] <- data.frame(
      Trait = trait,
      cohort = cohort,
      Trait_name_cohort = paste(trait, cohort, sep = "_"),
      M = est$M,
      M_prior = nrow(ind_snp_dt),
      theta1 = est$theta1,
      se_theta1 = est$se_theta1,
      pval_theta1 = est$pval_theta1,
      theta2 = est$theta2,
      se_theta2 = est$se_theta2,
      pval_theta2 = est$pval_theta2,
      I2 = est$I2,
      se_I2 = est$se_I2,
      pval_I2 = est$pval_I2,
      stringsAsFactors = FALSE
    )
  }
}

# ============================================================
# Save
# ============================================================

if (length(results) == 0L) {
  warning("No results generated.")
} else {
  results_dt <- rbindlist(results, fill = TRUE)
  
  fwrite(
    results_dt,
    output_file,
    sep = "\t",
    append = FALSE,
    row.names = FALSE
  )
  
  message("Results written to: ", output_file)
}
