suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# 1. Metadata required for 1KGP-weighted reference
# ============================================================

onekgp_info <- data.frame(
  population = c(
    "ACB","ASW","BEB","CDX","CEU","CHB","CHS","CLM","ESN","FIN",
    "GBR","GIH","GWD","IBS","ITU","JPT","KHV","LWK","MSL","MXL",
    "PEL","PJL","PUR","STU","TSI","YRI"
  ),
  Nr = c(
    96,61,86,93,99,103,105,94,99,99,
    91,103,113,107,102,104,99,99,85,64,
    85,96,104,102,107,108
  ),
  stringsAsFactors = FALSE
)

# ============================================================
# 2. Helpers
# ============================================================

validate_onekgp_info <- function(onekgp_info) {
  required_cols <- c("population", "Nr")
  missing_cols <- setdiff(required_cols, names(onekgp_info))
  
  if (length(missing_cols) > 0L) {
    stop(
      "onekgp_info is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (anyDuplicated(onekgp_info$population)) {
    stop("Duplicated population names found in onekgp_info.")
  }
  
  if (any(is.na(onekgp_info$Nr)) || any(onekgp_info$Nr <= 0)) {
    stop("All Nr values in onekgp_info must be positive and non-missing.")
  }
  
  invisible(TRUE)
}

validate_cohort_info <- function(cohort_info) {
  required_cols <- c("cohort", "Ns")
  missing_cols <- setdiff(required_cols, names(cohort_info))
  
  if (length(missing_cols) > 0L) {
    stop(
      "cohort_info is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (anyDuplicated(cohort_info$cohort)) {
    stop("Duplicated cohort names found in cohort_info.")
  }
  
  invisible(TRUE)
}

get_trait_name <- function(file_path) {
  strsplit(basename(file_path), "_")[[1]][1]
}

get_ld_panel_files <- function(ind_snp_dir, ld_ref) {
  ld_path <- file.path(ind_snp_dir, ld_ref)
  
  if (!dir.exists(ld_path)) {
    stop("LD panel directory not found: ", ld_path)
  }
  
  list.files(ld_path, full.names = TRUE)
}

get_cohort_sample_size <- function(cohort_name, cohort_info) {
  idx <- match(cohort_name, cohort_info$cohort)
  
  if (is.na(idx)) {
    stop("Sample size Ns not found for cohort: ", cohort_name)
  }
  
  cohort_info$Ns[idx]
}

get_cohort_ld_ref <- function(cohort_name, ld_panel_info) {
  idx <- match(cohort_name, ld_panel_info$cohort)
  
  if (is.na(idx)) {
    stop("LD reference not found for cohort: ", cohort_name)
  }
  
  ld_panel_info$ld_ref[idx]
}

get_cohort_weights <- function(cohort_name, weights_dt, population_order) {
  if (!"cohort" %in% names(weights_dt)) {
    stop("weights table must contain a 'cohort' column")
  }
  
  missing_pops <- setdiff(population_order, names(weights_dt))
  if (length(missing_pops) > 0L) {
    stop(
      "weights table is missing population columns: ",
      paste(missing_pops, collapse = ", ")
    )
  }
  
  idx <- match(cohort_name, weights_dt$cohort)
  
  if (is.na(idx)) {
    stop("No ancestry weights found for cohort: ", cohort_name)
  }
  
  w <- as.numeric(weights_dt[idx, ..population_order])
  
  if (length(w) != length(population_order)) {
    stop("Weight vector length does not match population order.")
  }
  
  if (!is.finite(sum(w)) || sum(w) <= 0) {
    stop("Invalid ancestry weights for cohort: ", cohort_name)
  }
  
  w / sum(w)
}

compute_weighted_reference_variance <- function(ref_matrix, Nr_vec, weights_vec) {
  ref_var <- sweep(ref_matrix * (1 - ref_matrix), 2, 2 * Nr_vec, FUN = "/")
  as.vector(ref_var %*% (weights_vec^2))
}

prepare_analysis_data_1kgp <- function(
    freq_dt,
    beta_dt,
    ind_snp_dt,
    cohort_name,
    population_order,
    weights_vec,
    outlier_sd_threshold = 4
) {
  required_beta <- c("SNP", "A1", "BETA")
  missing_beta <- setdiff(required_beta, names(beta_dt))
  
  if (length(missing_beta) > 0L) {
    stop(
      "beta file is missing required columns: ",
      paste(missing_beta, collapse = ", ")
    )
  }
  
  if (!"SNP" %in% names(ind_snp_dt)) {
    stop("independent SNP file is missing SNP column")
  }
  
  required_freq <- c("SNP", "REF", "ALT", cohort_name, population_order)
  missing_freq <- setdiff(required_freq, names(freq_dt))
  
  if (length(missing_freq) > 0L) {
    stop(
      "frequency file is missing required columns: ",
      paste(missing_freq, collapse = ", ")
    )
  }
  
  beta_sub <- beta_dt[SNP %in% ind_snp_dt$SNP]
  
  freq_sub <- freq_dt[, c("SNP", "REF", "ALT", cohort_name, population_order), with = FALSE]
  
  analysis_dt <- merge(freq_sub, beta_sub, by = "SNP", all = FALSE)
  
  if (nrow(analysis_dt) == 0L) {
    return(analysis_dt)
  }
  
  analysis_dt[, BETA := as.numeric(BETA)]
  analysis_dt[, beta_aligned := fifelse(ALT == A1, BETA, -BETA)]
  
  ref_matrix <- as.matrix(analysis_dt[, ..population_order])
  storage.mode(ref_matrix) <- "numeric"
  
  analysis_dt[, weighted_pr := as.numeric(ref_matrix %*% weights_vec)]
  
  diff_vec <- analysis_dt[[cohort_name]] - analysis_dt$weighted_pr
  diff_sd <- sd(diff_vec, na.rm = TRUE)
  
  if (is.finite(diff_sd) && diff_sd > 0) {
    keep_idx <- which(abs(diff_vec) <= outlier_sd_threshold * diff_sd)
    analysis_dt <- analysis_dt[keep_idx]
  }
  
  analysis_dt
}

estimate_ascertainment_1kgp <- function(
    analysis_dt,
    cohort_name,
    Ns,
    Nr_vec,
    weights_vec,
    population_order,
    min_snp_count = 10L
) {
  keep <- !is.na(analysis_dt$beta_aligned) &
    !is.na(analysis_dt[[cohort_name]]) &
    !is.na(analysis_dt$weighted_pr)
  
  analysis_dt <- analysis_dt[keep]
  
  M <- nrow(analysis_dt)
  
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
  
  b <- analysis_dt$beta_aligned
  p_s <- analysis_dt[[cohort_name]]
  p_r_mat <- as.matrix(analysis_dt[, ..population_order])
  storage.mode(p_r_mat) <- "numeric"
  w_pr <- analysis_dt$weighted_pr
  
  u <- sum(2 * b * (p_s - w_pr))
  v <- sum(2 * b^2 * w_pr * (1 - w_pr))
  
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
  var_pr <- compute_weighted_reference_variance(
    ref_matrix = p_r_mat,
    Nr_vec = Nr_vec,
    weights_vec = weights_vec
  )
  
  var_u <- sum(4 * b^2 * (var_ps + var_pr))
  var_v <- sum(4 * b^4 * (var_pr * (1 - 2 * w_pr)^2 + 2 * var_pr^2))
  cov_uv <- sum(-4 * b^3 * (1 - 2 * w_pr) * var_pr)
  
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
  
  lhs <- p_s - w_pr
  rhs <- b * w_pr * (1 - w_pr) / sqrt(v)
  
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
# 3. Main public function
# ============================================================

run_ascertainment_with_1kgp_reference <- function(
    beta_dir,
    freq_dir,
    ind_snp_dir,
    weights_path,
    ld_panel_path,
    cohort_info,
    onekgp_info,
    output_dir = "result_1kgp_rep",
    output_prefix = "ascertainment_1kgp",
    min_snp_count = 10L,
    outlier_sd_threshold = 4
) {
  validate_onekgp_info(onekgp_info)
  validate_cohort_info(cohort_info)
  
  if (!dir.exists(beta_dir)) {
    stop("beta_dir does not exist: ", beta_dir)
  }
  
  if (!dir.exists(freq_dir)) {
    stop("freq_dir does not exist: ", freq_dir)
  }
  
  if (!dir.exists(ind_snp_dir)) {
    stop("ind_snp_dir does not exist: ", ind_snp_dir)
  }
  
  if (!file.exists(weights_path)) {
    stop("weights_path does not exist: ", weights_path)
  }
  
  if (!file.exists(ld_panel_path)) {
    stop("ld_panel_path does not exist: ", ld_panel_path)
  }
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  population_order <- onekgp_info$population
  Nr_vec <- onekgp_info$Nr
  
  weights_dt <- fread(weights_path)
  ld_panel_info <- fread(ld_panel_path)
  
  beta_files <- list.files(beta_dir, full.names = TRUE)
  freq_files <- list.files(freq_dir, full.names = TRUE)
  
  if (length(beta_files) == 0L) {
    stop("No beta files found in: ", beta_dir)
  }
  
  if (length(freq_files) == 0L) {
    stop("No frequency files found in: ", freq_dir)
  }
  
  results <- vector("list", 0L)
  
  for (freq_file in freq_files) {
    message("Reading frequency file: ", basename(freq_file))
    
    freq_dt <- as.data.table(fread(freq_file))
    
    fixed_cols <- c("SNP", "CHR", "REF", "ALT")
    freq_cols <- setdiff(names(freq_dt), fixed_cols)
    cohort_candidates <- setdiff(freq_cols, population_order)
    
    if (length(cohort_candidates) != 1L) {
      stop(
        "Expected exactly one cohort column in 1KGP frequency file; found: ",
        paste(cohort_candidates, collapse = ", "),
        " in file ", basename(freq_file)
      )
    }
    
    cohort_name <- cohort_candidates[1]
    message("Cohort identified as: ", cohort_name)
    
    Ns <- get_cohort_sample_size(cohort_name, cohort_info)
    ld_ref <- get_cohort_ld_ref(cohort_name, ld_panel_info)
    weights_vec <- get_cohort_weights(cohort_name, weights_dt, population_order)
    
    message("Using LD panel ", ld_ref, " for cohort ", cohort_name)
    
    ld_panel_files <- get_ld_panel_files(ind_snp_dir, ld_ref)
    ld_panel_lookup <- setNames(ld_panel_files, basename(ld_panel_files))
    
    for (beta_file in beta_files) {
      beta_basename <- basename(beta_file)
      
      if (!beta_basename %in% names(ld_panel_lookup)) {
        message("  Skipping (no independent SNP file match): ", beta_basename)
        next
      }
      
      snp_file <- ld_panel_lookup[[beta_basename]]
      trait_name <- get_trait_name(beta_file)
      
      message("  Processing trait: ", trait_name)
      
      beta_dt <- as.data.table(fread(beta_file))
      ind_snp_dt <- as.data.table(fread(snp_file))
      
      analysis_dt <- prepare_analysis_data_1kgp(
        freq_dt = freq_dt,
        beta_dt = beta_dt,
        ind_snp_dt = ind_snp_dt,
        cohort_name = cohort_name,
        population_order = population_order,
        weights_vec = weights_vec,
        outlier_sd_threshold = outlier_sd_threshold
      )
      
      est <- estimate_ascertainment_1kgp(
        analysis_dt = analysis_dt,
        cohort_name = cohort_name,
        Ns = Ns,
        Nr_vec = Nr_vec,
        weights_vec = weights_vec,
        population_order = population_order,
        min_snp_count = min_snp_count
      )
      
      result_name <- paste(trait_name, cohort_name, sep = "_")
      
      results[[result_name]] <- data.frame(
        Trait = trait_name,
        cohort = cohort_name,
        Trait_name_cohort = result_name,
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
  
  if (length(results) == 0L) {
    warning("No results were generated.")
    return(invisible(NULL))
  }
  
  results_dt <- rbindlist(results, fill = TRUE)
  
  output_file <- file.path(output_dir, paste0(output_prefix, "_results.tsv"))
  fwrite(results_dt, output_file, sep = "\t")
  
  message("Results written to: ", output_file)
  
  invisible(results_dt)
}

# ============================================================
# 4. Example usage
# ============================================================

cohort_info <- data.frame(
  cohort = c("SweGen", "TWB"),
  Ns = c(1000, 10000),
  stringsAsFactors = FALSE
)

res <- run_ascertainment_with_1kgp_reference(
  beta_dir = "beta",
  freq_dir = "sample_with_1KGP",
  ind_snp_dir = "independent_snps",
  weights_path = "result_anc_loadings/ancestry_population_weights.tsv",
  ld_panel_path = "result_anc_loadings/ancestry_ld_panel_info.tsv",
  cohort_info = cohort_info,
  onekgp_info = onekgp_info,
  output_dir = "result_1kgp_rep",
  output_prefix = "ascertainment_1kgp",
  min_snp_count = 10L,
  outlier_sd_threshold = 4
)
