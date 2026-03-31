## set working directory ##

#setwd("/Users/uqbolase/Documents/github/")

## load required packages ###
pkgs <- c("data.table", "dplyr", "nnls", "quadprog", "pheatmap", "ggplot2")
invisible(lapply(pkgs, library, character.only = TRUE))

# ============================================================
# 1. Metadata
# ============================================================

onekgp_info <- data.frame(
  population = c(
    "ACB","ASW","BEB","CDX","CEU","CHB","CHS","CLM","ESN","FIN",
    "GBR","GIH","GWD","IBS","ITU","JPT","KHV","LWK","MSL","MXL",
    "PEL","PJL","PUR","STU","TSI","YRI"
  ),
  super_pop = c(
    "AFR","AFR","SAS","EAS","EUR","EAS","EAS","AMR","AFR","EUR",
    "EUR","SAS","AFR","EUR","SAS","EAS","EAS","AFR","AFR","AMR",
    "AMR","SAS","AMR","SAS","EUR","AFR"
  ),
  Nr = c(
    96,61,86,93,99,103,105,94,99,99,
    91,103,113,107,102,104,99,99,85,64,
    85,96,104,102,107,108
  ),
  stringsAsFactors = FALSE
)

ann_colors <- list(
  super_pop = c(
    AFR = "black",
    AMR = "pink",
    EAS = "brown",
    EUR = "forestgreen",
    SAS = "purple"
  )
)

# ============================================================
# 2. Helpers
# ============================================================

validate_onekgp_info <- function(onekgp_info) {
  required_cols <- c("population", "super_pop", "Nr")
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

read_frequency_files <- function(freq_dir, cohort_names = NULL) {
  if (!dir.exists(freq_dir)) {
    stop("Directory does not exist: ", freq_dir)
  }
  
  freq_files <- list.files(freq_dir, full.names = TRUE)
  
  if (length(freq_files) == 0L) {
    stop("No files found in freq_dir: ", freq_dir)
  }
  
  if (is.null(cohort_names)) {
    cohort_names <- tools::file_path_sans_ext(basename(freq_files))
  }
  
  if (length(cohort_names) != length(freq_files)) {
    stop(
      "cohort_names must have the same length as the number of files in freq_dir."
    )
  }
  
  data.frame(
    cohort = cohort_names,
    file = freq_files,
    stringsAsFactors = FALSE
  )
}

check_frequency_columns <- function(dt, cohort_name, onekgp_info) {
  fixed_cols <- c("SNP", "CHR", "REF", "ALT")
  missing_fixed <- setdiff(fixed_cols, names(dt))
  
  if (length(missing_fixed) > 0L) {
    stop(
      "Missing required columns: ",
      paste(missing_fixed, collapse = ", ")
    )
  }
  
  freq_cols <- setdiff(names(dt), fixed_cols)
  
  if (!cohort_name %in% freq_cols) {
    stop("Cohort column not found in file: ", cohort_name)
  }
  
  predictor_pops <- setdiff(freq_cols, cohort_name)
  
  missing_pops <- setdiff(predictor_pops, onekgp_info$population)
  
  if (length(missing_pops) > 0L) {
    stop(
      "These predictor populations are not found in onekgp_info: ",
      paste(missing_pops, collapse = ", ")
    )
  }
  
  predictor_pops
}

estimate_corrected_weights <- function(X, y, Nr_vec) {
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  
  y <- as.matrix(y)
  
  if (ncol(y) != 1L) {
    stop("y must be a single-column matrix or vector.")
  }
  
  if (length(Nr_vec) != ncol(X)) {
    stop("Length of Nr_vec must equal the number of columns in X.")
  }
  
  Cx <- cov(X)
  Cy <- cov(X, y)
  
  s <- colMeans(sweep(X * (1 - X), 2, Nr_vec, "/"))
  A <- Cx - 0.5 * diag(s)
  b <- Cy[, 1]
  
  fit <- nnls(A, b)
  w <- coef(fit)
  
  if (!is.finite(sum(w)) || sum(w) <= 0) {
    stop("Estimated ancestry weights are invalid or sum to zero.")
  }
  
  w / sum(w)
}

collapse_to_superpop <- function(population_weights_long, onekgp_info) {
  dt <- merge(
    as.data.table(population_weights_long),
    as.data.table(onekgp_info)[, .(population, super_pop)],
    by = "population",
    all.x = TRUE
  )
  
  if (dt[is.na(super_pop), .N] > 0L) {
    stop("Some populations could not be mapped to super_pop.")
  }
  
  dt[
    ,
    .(weight = sum(weight)),
    by = .(cohort, super_pop)
  ][order(cohort, -weight)]
}

choose_ld_panel <- function(superpop_dt,
                            dominant_threshold = 0.90,
                            min_component = 0.10) {
  superpop_dt <- as.data.table(copy(superpop_dt))[order(-weight)]
  
  if (nrow(superpop_dt) == 0L) {
    return(NA_character_)
  }
  
  if (superpop_dt$weight[1] >= dominant_threshold) {
    return(superpop_dt$super_pop[1])
  }
  
  keep <- superpop_dt[weight >= min_component, super_pop]
  
  if (length(keep) == 0L) {
    return(superpop_dt$super_pop[1])
  }
  
  paste(keep, collapse = "_")
}

make_dominant_ancestry_table <- function(superpop_weights_long) {
  as.data.table(superpop_weights_long)[
    order(cohort, -weight),
    .SD[1],
    by = cohort
  ][
    ,
    .(
      cohort = cohort,
      dominant_super_pop = super_pop,
      dominant_weight = weight
    )
  ]
}

# ============================================================
# 3. Main public function
# ============================================================

estimate_ancestry_loadings <- function(
    freq_dir,
    cohort_names = NULL,
    onekgp_info,
    sample_n_snp = 1e5,
    dominant_threshold = 0.90,
    min_component = 0.10
) {
  validate_onekgp_info(onekgp_info)
  
  file_map <- read_frequency_files(freq_dir = freq_dir, cohort_names = cohort_names)
  
  population_weights_list <- vector("list", nrow(file_map))
  names(population_weights_list) <- file_map$cohort
  
  for (i in seq_len(nrow(file_map))) {
    cohort_name <- file_map$cohort[i]
    freq_file <- file_map$file[i]
    
    message("Processing cohort: ", cohort_name)
    
    dt <- as.data.table(fread(freq_file))
    predictor_pops <- check_frequency_columns(
      dt = dt,
      cohort_name = cohort_name,
      onekgp_info = onekgp_info
    )
    
    dt <- dt[, c("SNP", "CHR", "REF", "ALT", cohort_name, predictor_pops), with = FALSE]
    
    if (nrow(dt) > sample_n_snp) {
      dt <- dt[sample(.N, sample_n_snp)]
    }
    
    y <- as.matrix(dt[[cohort_name]])
    X <- as.matrix(dt[, ..predictor_pops])
    
    Nr_vec <- onekgp_info$Nr[match(colnames(X), onekgp_info$population)]
    
    if (anyNA(Nr_vec)) {
      stop("Could not match predictor populations to onekgp_info$population.")
    }
    
    weights_cor <- estimate_corrected_weights(X = X, y = y, Nr_vec = Nr_vec)
    
    population_weights_list[[cohort_name]] <- data.table(
      cohort = cohort_name,
      population = colnames(X),
      weight = weights_cor
    )
  }
  
  population_weights_long <- rbindlist(
    population_weights_list,
    use.names = TRUE,
    fill = TRUE
  )
  
  population_weights_long <- population_weights_long[order(cohort, population)]
  
  population_weights_wide <- dcast(
    population_weights_long,
    cohort ~ population,
    value.var = "weight",
    fill = 0
  )
  
  population_weights_wide <- population_weights_wide[
    ,
    c("cohort", onekgp_info$population),
    with = FALSE
  ]
  
  superpop_weights_long <- collapse_to_superpop(
    population_weights_long = population_weights_long,
    onekgp_info = onekgp_info
  )
  
  superpop_weights_wide <- dcast(
    superpop_weights_long,
    cohort ~ super_pop,
    value.var = "weight",
    fill = 0
  )
  
  ld_panel_info <- superpop_weights_long[
    ,
    .(
      ld_ref = choose_ld_panel(
        .SD,
        dominant_threshold = dominant_threshold,
        min_component = min_component
      )
    ),
    by = cohort
  ]
  
  dominant_ancestry <- make_dominant_ancestry_table(superpop_weights_long)
  
  list(
    population_weights_long = population_weights_long,
    population_weights_wide = population_weights_wide,
    superpop_weights_long = superpop_weights_long,
    superpop_weights_wide = superpop_weights_wide,
    dominant_ancestry = dominant_ancestry,
    ld_panel_info = ld_panel_info
  )
}

# ============================================================
# 4. Plotting helpers
# ============================================================

plot_population_heatmap <- function(
    population_weights_wide,
    onekgp_info,
    main = "Corrected ancestry loadings",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_row = 10,
    fontsize_col = 10,
    angle_col = "270"
) {
  
  superpop_levels <- c("AFR", "AMR", "EAS", "EUR", "SAS")
  
  onekgp_ordered <- onekgp_info[
    order(match(onekgp_info$super_pop, superpop_levels),
          onekgp_info$population),
  ]
  
  pop_order <- onekgp_ordered$population
  
  if (!all(pop_order %in% names(population_weights_wide))) {
    stop("population_weights_wide is missing one or more required population columns.")
  }
  
  mat <- as.matrix(population_weights_wide[, ..pop_order])
  rownames(mat) <- population_weights_wide$cohort
  
  col_annot <- data.frame(
    super_pop = onekgp_ordered$super_pop,
    row.names = onekgp_ordered$population
  )
  
  pheatmap(
    mat = mat,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    annotation_col = col_annot,
    annotation_colors = ann_colors,
    angle_col = angle_col,
    main = main,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col
  )
}


plot_superpop_bar <- function(
    superpop_weights_long,
    dominant_threshold = 0.90
) {
  superpop_levels <- c("AFR", "AMR", "EAS", "EUR", "SAS")
  superpop_weights_long <- copy(as.data.table(superpop_weights_long))
  superpop_weights_long[, super_pop := factor(super_pop, levels = superpop_levels)]
  
  ggplot(
    superpop_weights_long,
    aes(x = cohort, y = weight, fill = super_pop)
  ) +
    geom_col(width = 0.8) +
    geom_hline(yintercept = dominant_threshold, linetype = 2) +
    scale_fill_manual(values = ann_colors$super_pop, drop = FALSE) +
    labs(
      x = NULL,
      y = "Ancestry loading",
      fill = "Super-population"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

plot_dominant_ancestry <- function(dominant_ancestry) {
  superpop_levels <- c("AFR", "AMR", "EAS", "EUR", "SAS")
  dominant_ancestry <- copy(as.data.table(dominant_ancestry))
  dominant_ancestry[
    ,
    dominant_super_pop := factor(dominant_super_pop, levels = superpop_levels)
  ]
  
  ggplot(
    dominant_ancestry,
    aes(x = cohort, y = dominant_weight, fill = dominant_super_pop)
  ) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = ann_colors$super_pop, drop = FALSE) +
    labs(
      x = NULL,
      y = "Dominant ancestry loading",
      fill = "Dominant super-population"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

# ============================================================
# 5. Convenience wrapper
# ============================================================

run_ancestry_pipeline <- function(
    freq_dir,
    cohort_names = NULL,
    onekgp_info,
    output_prefix = "ancestry",
    output_dir = "result_anc_loadings",
    sample_n_snp = 1e5,
    dominant_threshold = 0.90,
    min_component = 0.10
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  res <- estimate_ancestry_loadings(
    freq_dir = freq_dir,
    cohort_names = cohort_names,
    onekgp_info = onekgp_info,
    sample_n_snp = sample_n_snp,
    dominant_threshold = dominant_threshold,
    min_component = min_component
  )
  
  fwrite(
    res$population_weights_wide,
    file.path(output_dir, paste0(output_prefix, "_population_weights.tsv")),
    sep = "\t"
  )
  
  fwrite(
    res$superpop_weights_wide,
    file.path(output_dir, paste0(output_prefix, "_superpop_weights.tsv")),
    sep = "\t"
  )
  
  fwrite(
    res$dominant_ancestry,
    file.path(output_dir, paste0(output_prefix, "_dominant.tsv")),
    sep = "\t"
  )
  
  fwrite(
    res$ld_panel_info,
    file.path(output_dir, paste0(output_prefix, "_ld_panel_info.tsv")),
    sep = "\t"
  )
  
  png(
    filename = file.path(output_dir, paste0(output_prefix, "_population_heatmap.png")),
    width = 12,
    height = 6,
    units = "in",
    res = 300
  )
  plot_population_heatmap(
    population_weights_wide = res$population_weights_wide,
    onekgp_info = onekgp_info
  )
  dev.off()
  
  ggsave(
    filename = file.path(output_dir, paste0(output_prefix, "_superpop_barplot.png")),
    plot = plot_superpop_bar(
      superpop_weights_long = res$superpop_weights_long,
      dominant_threshold = dominant_threshold
    ),
    width = 10,
    height = 5,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(output_dir, paste0(output_prefix, "_dominant_ancestry.png")),
    plot = plot_dominant_ancestry(res$dominant_ancestry),
    width = 8,
    height = 5,
    dpi = 300
  )
  
  invisible(res)
}

# ============================================================
# 6. Example usage
# ============================================================

# Expected input directory structure:
# sample_with_1KGP/
#   SweGen
#   TWB
#
# Each file should contain:
# SNP, CHR, REF, ALT, <cohort_name>, and 1KGP population columns.

cohort_names <- c("SweGen", "TWB")
freq_dir <- "data/sample_with_1KGP"

res <- run_ancestry_pipeline(
  freq_dir = freq_dir,
  cohort_names = cohort_names,
  onekgp_info = onekgp_info,
  output_prefix = "ancestry",
  output_dir = "result_anc_loadings",
  sample_n_snp = 1e5,
  dominant_threshold = 0.90,
  min_component = 0.10
)

res$population_weights_wide
res$superpop_weights_wide
res$dominant_ancestry
res$ld_panel_info
