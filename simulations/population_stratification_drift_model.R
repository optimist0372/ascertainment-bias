## simulation to investigate the impact of population stratification ##
## on theta1 and theta2 using drift model ##

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

Sim_PopStrat_drift <- function(M,Pi_s,Pi_r,  Fst, N_s, N_r){
  
  #────────────────────────────────────────────
  # 1. Sample size decomposition by ancestry
  #────────────────────────────────────────────
  N1_s <- N_s * Pi_s
  N2_s <- N_s * (1 - Pi_s)
  
  N1_r <- N_r * Pi_r
  N2_r <- N_r * (1 - Pi_r)
  
  #────────────────────────────────────────────
  # 2. Ancestral allele frequency generation
  #────────────────────────────────────────────
  Ep0  <- 0.05 + 0.9 * rbeta(n=M,shape1=0.8,shape2=0.8)
  
  # population-specific frequencies under Fst
  Ep1  <- rbeta(M,shape1=Ep0*(1-Fst)/Fst,shape2=(1-Ep0)*(1-Fst)/Fst) # Pop 1
  Ep2  <- rbeta(M,shape1=Ep0*(1-Fst)/Fst,shape2=(1-Ep0)*(1-Fst)/Fst) # Pop 2 
  
  #────────────────────────────────────────────
  # 3. Generate trait-increasing alllele effect sizes 
  #────────────────────────────────────────────
  vb   <- ( Ep0 * (1-Ep0) )**(-0.5)
  b    <- abs( rnorm(n=M,sd=sqrt(vb)) )
  
  #────────────────────────────────────────────
  # 4. Sample allele frequencies (study population)
  #────────────────────────────────────────────
  p1_s <- rbeta(M,shape1=Ep1*(2*N1_s-1),shape2=(1-Ep1)*(2*N1_s-1))
  p2_s <- rbeta(M,shape1=Ep2*(2*N2_s-1),shape2=(1-Ep2)*(2*N2_s-1))
  p_s  <- Pi_s * p1_s + (1-Pi_s) * p2_s
  
  #────────────────────────────────────────────
  # 5. Sample allele frequencies (reference population)
  #────────────────────────────────────────────
  p1_r  <- rbeta(M,shape1=Ep1*(2*N1_r-1),shape2=(1-Ep1)*(2*N1_r-1))
  p2_r  <- rbeta(M,shape1=Ep2*(2*N2_r-1),shape2=(1-Ep2)*(2*N2_r-1))
  p_r   <- Pi_r * p1_r + (1-Pi_r) * p2_r
  
  #────────────────────────────────────────────
  # 6. Sampling variance of allele frequencies
  #────────────────────────────────────────────
  # var_p_s   <- p_s*(1-p_s)/(2 * N_s)
  # var_p_r   <- p_r*(1-p_r)/(2 * N_r)
  # 
  #────────────────────────────────────────────
  # 6. Theta1 estimation (score-based statistic)
  #────────────────────────────────────────────
  u   <- sum( 2*b * (p_s - p_r) )
  v   <- sum(b^2 * 2 * p_r * (1 - p_r))
  theta1 <- u / sqrt(v)
  
  #────────────────────────────────────────────
  # 7. Variance of theta1 (delta-method expansion)
  #────────────────────────────────────────────
  var_ps <- p_s * (1 - p_s) / (2 * N_s)
  var_pr <- p_r * (1 - p_r) / (2 * N_r)
  
  var_u <- sum(4 * b^2 * (var_ps + var_pr))
  var_v <- sum(4 * b^4 * (var_pr * (1 - 2 * p_r)^2 + 2 * var_pr^2))
  cov_uv <- sum(-4 * b^3 * (1 - 2 * p_r) * var_pr)
  
  var_theta1 <- (1 / v) * (var_u + (u^2 * var_v) / (4 * v^2) - (u * cov_uv / v))
  sd_theta1 <- sqrt(var_theta1)
  
  #────────────────────────────────────────────
  # 8. Hypothesis test for theta1 (Wald test)
  #────────────────────────────────────────────
  ts <- (theta1 / sd_theta1)^2
  pval_theta1 <- pchisq(ts, df = 1, lower.tail = FALSE)
  
  #────────────────────────────────────────────
  # 9. Theta2 estimation (weighted regression)
  #────────────────────────────────────────────
  rhs <- b * p_r * (1 - p_r) / sqrt(v)
  lhs <- p_s - p_r
  iv_weight <- 1 / var_pr
  
  fit <- lm(lhs ~ rhs, weights = iv_weight)
  coef_tab <- summary(fit)$coefficients
  
  theta2 <- coef_tab[2, 1]
  pval_theta2 <- coef_tab[2, 4]
  
  #────────────────────────────────────────────
  # 10. Intercept (I2) from regression
  #────────────────────────────────────────────
  I2 <- coef_tab[1, 1]
  pval_I2 <- coef_tab[1, 4]
  
  #────────────────────────────────────────────
  # 11. Output summary statistics
  #────────────────────────────────────────────
  stats <- c(theta1     = theta1,
             p_theta1   = pval_theta1,
             
             theta2     = theta2,
             p_theta2   = pval_theta2,
             
             I2 = I2,
             p_I2 = pval_I2)
}


# ------------------------------------------------------------
# 1. Simulation settings
# ------------------------------------------------------------
pars <- list(
  M    = 1000,
  N_s  = 5000,
  N_r  = 2500,
  Fst  = c(0.001, 0.01),
  Pi_s = 0.5,
  Pi_r = seq(0.1, 0.9, by = 0.2),
  nrep = 10000,
  alpha = 0.05
)

# ------------------------------------------------------------
# 2. Design grid
# ------------------------------------------------------------
make_design <- function(pars) {
  design <- expand.grid(
    M    = pars$M,
    Pi_s = pars$Pi_s,
    Pi_r = pars$Pi_r,
    Fst  = pars$Fst,
    N_s  = pars$N_s,
    N_r  = pars$N_r
  )
  
  do.call(
    "rbind",
    lapply(seq_len(pars$nrep), function(k) cbind(Rep = k, design))
  )
}

Design <- make_design(pars)

# ------------------------------------------------------------
# 3. Simulation runner
# ------------------------------------------------------------
run_simulations <- function(Design, sim_fun) {
  pb <- txtProgressBar(min = 0, max = nrow(Design), style = 3)
  on.exit(close(pb), add = TRUE)
  
  res <- lapply(seq_len(nrow(Design)), function(i) {
    setTxtProgressBar(pb, i)
    
    out <- sim_fun(
      M    = Design[i, "M"],
      Pi_s = Design[i, "Pi_s"],
      Pi_r = Design[i, "Pi_r"],
      Fst  = Design[i, "Fst"],
      N_s  = Design[i, "N_s"],
      N_r  = Design[i, "N_r"]
    )
    
    as.data.frame(as.list(out))
  })
  
  out <- rbindlist(res, fill = TRUE)
  
  # keep names explicit and stable
  expected <- c("theta1", "p_theta1", "theta2", "p_theta2", "I2", "p_I2")
  missing_cols <- setdiff(expected, names(out))
  if (length(missing_cols) > 0) {
    stop("Simulation output is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  out[, ..expected]
}

simulations <- run_simulations(Design, Sim_PopStrat_drift)

# one master table for everything downstream
results <- cbind(as.data.table(Design), simulations)

###summarize result across replicates and plot

# # ------------------------------------------------------------
# # 4. Generic summary helpers
# # ------------------------------------------------------------
# summarise_estimates <- function(results, value_cols, group_cols = c("Pi_s", "Pi_r", "Fst")) {
#   mean_dt <- results[, lapply(.SD, mean), by = group_cols, .SDcols = value_cols]
#   sd_dt   <- results[, lapply(.SD, sd),   by = group_cols, .SDcols = value_cols]
#   
#   nrep <- uniqueN(results$Rep)
#   se_dt <- copy(sd_dt)
#   se_names <- paste0(value_cols, "_se")
#   
#   for (j in seq_along(value_cols)) {
#     set(se_dt, j = value_cols[j], value = se_dt[[value_cols[j]]] / sqrt(nrep))
#   }
#   setnames(se_dt, old = value_cols, new = se_names)
#   
#   merge(mean_dt, se_dt, by = group_cols)
# }
# 
# summarise_type1 <- function(results, p_cols, alpha = 0.05, group_cols = c("Pi_s", "Pi_r", "Fst")) {
#   nrep <- uniqueN(results$Rep)
#   
#   out_list <- lapply(p_cols, function(pcol) {
#     nm <- sub("^p_", "", pcol)
#     err_nm <- paste0("type1_error_", nm)
#     se_nm  <- paste0("se_", nm)
#     
#     tmp <- results[, .(
#       stat = nm,
#       type1 = mean(get(pcol) < alpha)
#     ), by = group_cols]
#     
#     tmp[, (se_nm) := sqrt(type1 * (1 - type1) / nrep)]
#     setnames(tmp, "type1", err_nm)
#     tmp
#   })
#   
#   out_list
# }
# 
# # ------------------------------------------------------------
# # 5. Estimate summaries
# # ------------------------------------------------------------
# theta_sum <- summarise_estimates(results, value_cols = c("theta1", "theta2"))
# I2_sum    <- summarise_estimates(results, value_cols = c("I2"))
# 
# type1_list <- summarise_type1(results, p_cols = c("p_theta1", "p_theta2"), alpha = pars$alpha)
# type1_theta1 <- type1_list[[1]]
# type1_theta2 <- type1_list[[2]]
# 
# type1_theta <- merge(
#   type1_theta1[, stat := NULL],
#   type1_theta2[, stat := NULL],
#   by = c("Pi_s", "Pi_r", "Fst"),
#   all = TRUE
# )
# 
# type1_I2 <- summarise_type1(results, p_cols = c("p_I2"), alpha = pars$alpha)[[1]]
# 
# # ------------------------------------------------------------
# # 6. Tidy/long-format helpers
# # ------------------------------------------------------------
# make_theta_long <- function(df) {
#   bind_rows(
#     transmute(df, Pi_s, Pi_r, Fst, stat = "theta1", theta = theta1, se = theta1_se),
#     transmute(df, Pi_s, Pi_r, Fst, stat = "theta2", theta = theta2, se = theta2_se)
#   ) %>%
#     mutate(
#       Pi_s = as.numeric(Pi_s),
#       Fst  = as.numeric(Fst),
#       stat = factor(stat, levels = c("theta1", "theta2"))
#     )
# }
# 
# make_I2_long <- function(df) {
#   transmute(df, Pi_s, Pi_r, Fst, stat = "I2", theta = I2, se = I2_se) %>%
#     mutate(
#       Pi_s = as.numeric(Pi_s),
#       Fst  = as.numeric(Fst),
#       stat = factor(stat, levels = "I2")
#     )
# }
# 
# make_type1_long <- function(df) {
#   bind_rows(
#     transmute(df, Pi_s, Pi_r, Fst, stat = "theta1",
#               theta = type1_error_theta1, se = se_theta1),
#     transmute(df, Pi_s, Pi_r, Fst, stat = "theta2",
#               theta = type1_error_theta2, se = se_theta2)
#   ) %>%
#     mutate(
#       Pi_s = as.numeric(Pi_s),
#       Fst  = as.numeric(Fst),
#       stat = factor(stat, levels = c("theta1", "theta2"))
#     )
# }
# 
# make_type1_I2_long <- function(df) {
#   transmute(df, Pi_s, Pi_r, Fst, stat = "I2",
#             theta = type1_error_I2, se = se_I2) %>%
#     mutate(
#       Pi_s = as.numeric(Pi_s),
#       Fst  = as.numeric(Fst),
#       stat = factor(stat, levels = "I2")
#     )
# }
# 
# add_parsed_labels <- function(df) {
#   df %>%
#     mutate(
#       Pi_s_lab = sprintf("pi[s]==%.1f", Pi_s),
#       Fst_lab  = sprintf("F[ST]==%s", Fst)
#     )
# }
# 
# theta_long   <- make_theta_long(theta_sum)   %>% add_parsed_labels()
# I2_long      <- make_I2_long(I2_sum)         %>% add_parsed_labels()
# type1_long   <- make_type1_long(type1_theta) %>% add_parsed_labels()
# type1_long2  <- make_type1_I2_long(type1_I2) %>% add_parsed_labels()
# 
# # ------------------------------------------------------------
# # 7. Plot helper
# # ------------------------------------------------------------
# plot_theta <- function(df,
#                        x_var = "Pi_r",
#                        y_var = "theta",
#                        se_var = "se",
#                        group_var = "stat",
#                        stat_labels,
#                        colors,
#                        linetypes,
#                        shapes = NULL,
#                        yintercept = 0,
#                        y_lab = NULL,
#                        x_lab = NULL,
#                        y_limits = NULL,
#                        y_breaks = NULL,
#                        zero_floor = FALSE,
#                        dodge = FALSE,
#                        ci = FALSE) {
#   
#   pos <- if (dodge) position_dodge(width = 0.15) else position_identity()
#   ci_mult <- if (ci) 1.96 else 1
#   
#   levels_stat <- unique(df[[group_var]])
#   df1 <- df[df[[group_var]] == levels_stat[1], ]
#   df2 <- if (length(levels_stat) > 1) df[df[[group_var]] == levels_stat[2], ] else NULL
#   
#   p <- ggplot(df, aes_string(
#     x = x_var,
#     y = y_var,
#     colour = group_var,
#     linetype = group_var,
#     group = group_var
#   )) +
#     geom_hline(yintercept = yintercept, linetype = "solid", colour = "grey60")
#   
#   p <- p +
#     geom_errorbar(
#       data = df1,
#       aes(
#         x = .data[[x_var]],
#         ymin = .data[[y_var]] - ci_mult * .data[[se_var]],
#         ymax = .data[[y_var]] + ci_mult * .data[[se_var]]
#       ),
#       width = 0.03,
#       linewidth = 0.8,
#       colour = colors[1],
#       linetype = "solid",
#       inherit.aes = FALSE,
#       position = pos
#     )
#   
#   if (!is.null(df2)) {
#     p <- p +
#       geom_errorbar(
#         data = df2,
#         aes(
#           x = .data[[x_var]],
#           ymin = .data[[y_var]] - ci_mult * .data[[se_var]],
#           ymax = .data[[y_var]] + ci_mult * .data[[se_var]]
#         ),
#         width = 0.03,
#         linewidth = 0.8,
#         colour = colors[2],
#         linetype = "solid",
#         inherit.aes = FALSE,
#         position = pos
#       )
#   }
#   
#   p <- p +
#     geom_line(position = pos, linewidth = 0.9) +
#     geom_point(aes_string(shape = group_var), size = 2, position = pos) +
#     facet_grid(rows = vars(Fst_lab), cols = vars(Pi_s_lab), labeller = label_parsed) +
#     scale_x_continuous(breaks = sort(unique(df[[x_var]]))) +
#     scale_colour_manual(values = colors, labels = stat_labels, name = NULL) +
#     scale_linetype_manual(values = linetypes, labels = stat_labels, name = NULL)
#   
#   if (!is.null(shapes)) {
#     p <- p + scale_shape_manual(values = shapes, labels = stat_labels, name = NULL)
#   }
#   
#   if (!is.null(y_limits) || !is.null(y_breaks)) {
#     p <- p + scale_y_continuous(limits = y_limits, breaks = y_breaks)
#   } else {
#     p <- p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.10)))
#   }
#   
#   if (isTRUE(zero_floor)) {
#     p <- p + coord_cartesian(ylim = c(0, NA))
#   }
#   
#   p +
#     labs(x = x_lab, y = y_lab) +
#     theme_bw(base_size = 12) +
#     theme(
#       strip.background = element_rect(fill = "grey90"),
#       strip.text = element_text(face = "bold"),
#       axis.text.x = element_text(angle = 90, vjust = 0.5),
#       legend.position = "right",
#       panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2)
#     )
# }
# 
# # ------------------------------------------------------------
# # 8. Build panels
# # ------------------------------------------------------------
# stat_labels12 <- c(expression(hat(theta)[1]), expression(hat(theta)[2]))
# cols12 <- c("#1b9e77", "#d95f02")
# lt12   <- c("solid", "dashed")
# shp12  <- c(16, 17)
# 
# pA1 <- plot_theta(
#   theta_long,
#   stat_labels = stat_labels12,
#   colors = cols12,
#   linetypes = lt12,
#   shapes = shp12,
#   y_lab = expression(hat(theta) ~ "\u00B1" ~ "95% CI"),
#   x_lab = expression(pi[r]),
#   ci = TRUE
# )
# 
# pC1 <- plot_theta(
#   type1_long,
#   stat_labels = stat_labels12,
#   colors = cols12,
#   linetypes = lt12,
#   shapes = shp12,
#   y_lab = expression("Type I error (" * hat(theta) ~ "\u00B1" ~ SE * ")"),
#   x_lab = expression(pi[r]),
#   yintercept = 0.05,
#   zero_floor = TRUE,
#   ci = FALSE
# )
# 
# stat_labels2 <- c(expression(I[2]))
# cols2 <- "#7570b3"
# lt2   <- "solid"
# shp2  <- 16
# 
# pB1 <- plot_theta(
#   I2_long,
#   stat_labels = stat_labels2,
#   colors = cols2,
#   linetypes = lt2,
#   shapes = shp2,
#   y_lab = expression(I[2] ~ "\u00B1" ~ "95% CI"),
#   x_lab = expression(pi[r]),
#   ci = TRUE
# )
# 
# pD1 <- plot_theta(
#   type1_long2,
#   stat_labels = stat_labels2,
#   colors = cols2,
#   linetypes = lt2,
#   shapes = shp2,
#   y_lab = expression("Type I error (" * "I"[2] ~ "\u00B1" ~ SE * ")"),
#   x_lab = expression(pi[r]),
#   yintercept = 0.05,
#   y_limits = c(0, 0.2),
#   y_breaks = seq(0, 0.2, by = 0.1),
#   ci = FALSE
# )
# 
# # ------------------------------------------------------------
# # 9. Final figure
# # ------------------------------------------------------------
# final_plot <- ((pA1 | pC1) / (pB1 | pD1)) +
#   plot_layout(guides = "collect") +
#   plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
#   theme(legend.position = "right")
# 
# print(final_plot)