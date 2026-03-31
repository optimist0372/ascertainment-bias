###simulation to assess how deconvolution approach mitigate bias due ##
# to population stratification using extreme model ####

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})


Sim_deconv_extreme <- function(M, Pi_s, Fst, N_s, N1_r, N2_r, C=0.01){
  
  #────────────────────────────────────────────
  # 1. Sample size decomposition by ancestry
  #────────────────────────────────────────────
  N1_s <- N_s * Pi_s
  N2_s <- N_s * (1 - Pi_s)
  
  #────────────────────────────────────────────
  # 2. Ancestral allele frequency generation
  #────────────────────────────────────────────
  
  Ep0  <- 0.05 + 0.9 * rbeta(n=M,shape1=0.8,shape2=0.8)
  
  #enforce systemic difference between Ep1 and Ep2
  p01 <-   Ep0 -C/2
  p02 <-   Ep0 +C/2
  
  Ep1  <- rbeta(M,shape1=p01*(1-Fst)/Fst,shape2=(1-p01)*(1-Fst)/Fst) # Pop 1
  Ep2  <- rbeta(M,shape1=p02*(1-Fst)/Fst,shape2=(1-p02)*(1-Fst)/Fst)# Pop 2
  
  
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
  
  
  # 6. Re-estimate sample ancestral proportion from independent SNPs (Mr) using nnls
  
  Mr <- 10* M
  
  p1_s2 <- rbeta(Mr,shape1=Ep1*(2*N1_s-1),shape2=(1-Ep1)*(2*N1_s-1))
  p2_s2 <- rbeta(Mr,shape1=Ep2*(2*N2_s-1),shape2=(1-Ep2)*(2*N2_s-1))
  p_s2  <- Pi_s * p1_s2 + (1-Pi_s) * p2_s2
  
  p1_r2 <- rbeta(Mr,shape1=Ep1*(2*N1_r-1),shape2=(1-Ep1)*(2*N1_r-1))
  p2_r2 <- rbeta(Mr,shape1=Ep2*(2*N2_r-1),shape2=(1-Ep2)*(2*N2_r-1))
  
  Pi_s_hat <- coef(nnls(cbind(p1_r2, p2_r2), p_s2))[1]
  
  #weighted reference using estimated ancestral proportions
  p_r   <- Pi_s_hat * p1_r + (1-Pi_s_hat) * p2_r
  
  #────────────────────────────────────────────
  # 7. Theta1 estimation 
  #────────────────────────────────────────────
  u   <- sum( 2*b * (p_s - p_r) )
  v   <- sum(b^2 * 2 * p_r * (1 - p_r))
  theta1 <- u / sqrt(v)
  
  #────────────────────────────────────────────
  # 8. Variance of theta1 (delta-method expansion)
  #────────────────────────────────────────────
  var_ps <- p_s * (1 - p_s) / (2 * N_s)
  var_pr <- (Pi_s_hat^2) * p1_r * (1 - p1_r)/(2*N1_r) + ( (1-Pi_s_hat)^2) * p2_r * (1 - p2_r)/(2*N2_r) 
  
  var_u <- sum(4 * b^2 * (var_ps + var_pr))
  var_v <- sum(4 * b^4 * (var_pr * (1 - 2 * p_r)^2 + 2 * var_pr^2))
  cov_uv <- sum(-4 * b^3 * (1 - 2 * p_r) * var_pr)
  
  var_theta1 <- (1 / v) * (var_u + (u^2 * var_v) / (4 * v^2) - (u * cov_uv / v))
  sd_theta1 <- sqrt(var_theta1)
  
  #────────────────────────────────────────────
  # 9. Hypothesis test for theta1 (Wald test)
  #────────────────────────────────────────────
  ts <- (theta1 / sd_theta1)^2
  pval_theta1 <- pchisq(ts, df = 1, lower.tail = FALSE)
  
  #────────────────────────────────────────────
  # 10. Theta2 estimation (weighted regression)
  #────────────────────────────────────────────
  rhs <- b * p_r * (1 - p_r) / sqrt(v)
  lhs <- p_s - p_r
  iv_weight <- 1 / var_pr
  
  fit <- lm(lhs ~ rhs, weights = iv_weight)
  coef_tab <- summary(fit)$coefficients
  
  theta2 <- coef_tab[2, 1]
  pval_theta2 <- coef_tab[2, 4]
  
  #────────────────────────────────────────────
  # 11. Intercept (I2) from regression
  #────────────────────────────────────────────
  I2 <- coef_tab[1, 1]
  pval_I2 <- coef_tab[1, 4]
  
  #────────────────────────────────────────────
  # 12. Output summary statistics
  #────────────────────────────────────────────
  stats <- c(theta1     = theta1,
             p_theta1   = pval_theta1,
             
             theta2     = theta2,
             p_theta2   = pval_theta2,
             
             I2 = I2,
             p_I2 = pval_I2)
}

# ============================================================
# 1. Simulation settings
# ============================================================

M    <- 1000
N_s  <- 5000
N1_r <- N2_r <- 1000
Fst  <- c(0.001, 0.01)
Pi_s <- seq(0.1, 0.9, by = 0.2)

nreps <- 100
alpha <- 0.05

# ============================================================
# 2. Simulation design
# ============================================================

design <- expand.grid(
  M    = M,
  Pi_s = Pi_s,
  Fst  = Fst,
  N_s  = N_s,
  N1_r = N1_r,
  N2_r = N2_r
)

Design <- do.call(
  "rbind",
  lapply(seq_len(nreps), function(k) cbind(Rep = k, design))
)

# ============================================================
# 3. Run simulations
# ============================================================

run_simulations <- function(Design, sim_fun) {
  pb <- txtProgressBar(min = 0, max = nrow(Design), style = 3)
  on.exit(close(pb), add = TRUE)
  
  out <- lapply(seq_len(nrow(Design)), function(i) {
    setTxtProgressBar(pb, i)
    
    sim_fun(
      M    = Design[i, "M"],
      Pi_s = Design[i, "Pi_s"],
      Fst  = Design[i, "Fst"],
      N_s  = Design[i, "N_s"],
      N1_r = Design[i, "N1_r"],
      N2_r = Design[i, "N2_r"]
    )
  })
  
  out <- as.data.table(do.call("rbind", out))
  
  setnames(
    out,
    old = names(out),
    new = c("theta1", "p_theta1", "theta2", "p_theta2", "theta2_int", "p_theta2_int")
  )
  
  out
}

simulations <- run_simulations(Design, Sim_deconv_extreme)

results <- cbind(as.data.table(Design), simulations)


###summarize result across replicates and plot
# 
# # ============================================================
# # 4. Summary helpers
# # ============================================================
# {
#   summarise_estimates <- function(dt, value_cols, by_cols = c("Pi_s", "N2_r", "Fst")) {
#     mean_dt <- dt[, lapply(.SD, mean), by = by_cols, .SDcols = value_cols]
#     sd_dt   <- dt[, lapply(.SD, sd),   by = by_cols, .SDcols = value_cols]
#     
#     se_dt <- copy(sd_dt)
#     for (nm in value_cols) {
#       set(se_dt, j = nm, value = se_dt[[nm]] / sqrt(uniqueN(dt$Rep)))
#     }
#     setnames(se_dt, old = value_cols, new = paste0(value_cols, "_se"))
#     
#     merge(mean_dt, se_dt, by = by_cols)
#   }
#   
#   summarise_type1_error <- function(dt, p_cols, alpha = 0.05, by_cols = c("Pi_s", "N2_r", "Fst")) {
#     n_rep <- uniqueN(dt$Rep)
#     
#     out <- lapply(p_cols, function(p_col) {
#       stat_nm <- sub("^p_", "", p_col)
#       
#       tmp <- dt[, .(
#         value = mean(get(p_col) < alpha)
#       ), by = by_cols]
#       
#       tmp[, se := sqrt(value * (1 - value) / n_rep)]
#       tmp[, stat := stat_nm]
#       tmp
#     })
#     
#     rbindlist(out)
#   }
#   
#   # ============================================================
#   # 5. Summaries
#   # ============================================================
#   
#   theta_summary <- summarise_estimates(results, c("theta1", "theta2"))
#   theta2int_summary <- summarise_estimates(results, "theta2_int")
#   
#   theta_type1_summary <- summarise_type1_error(
#     results,
#     p_cols = c("p_theta1", "p_theta2"),
#     alpha = alpha
#   )
#   
#   theta2int_type1_summary <- summarise_type1_error(
#     results,
#     p_cols = "p_theta2_int",
#     alpha = alpha
#   )
#   
#   # ============================================================
#   # 6. Plot-ready data
#   # ============================================================
#   
#   add_facet_labels <- function(df) {
#     df %>%
#       mutate(
#         Pi_s = as.numeric(Pi_s),
#         Fst  = as.numeric(Fst),
#         Fst_lab = sprintf("F[ST]==%s", Fst)
#       )
#   }
#   
#   theta_est_plot <- bind_rows(
#     transmute(theta_summary, Pi_s, N2_r, Fst, stat = "theta1", value = theta1, se = theta1_se),
#     transmute(theta_summary, Pi_s, N2_r, Fst, stat = "theta2", value = theta2, se = theta2_se)
#   ) %>%
#     mutate(stat = factor(stat, levels = c("theta1", "theta2"))) %>%
#     add_facet_labels()
#   
#   theta_type1_plot <- theta_type1_summary %>%
#     mutate(stat = factor(stat, levels = c("theta1", "theta2"))) %>%
#     add_facet_labels()
#   
#   theta2int_est_plot <- theta2int_summary %>%
#     transmute(Pi_s, N2_r, Fst, stat = "theta2_int", value = theta2_int, se = theta2_int_se) %>%
#     mutate(stat = factor(stat, levels = "theta2_int")) %>%
#     add_facet_labels()
#   
#   theta2int_type1_plot <- theta2int_type1_summary %>%
#     mutate(stat = factor(stat, levels = "theta2_int")) %>%
#     add_facet_labels()
#   
#   # ============================================================
#   # 7. Plot helper
#   # ============================================================
#   
#   plot_sim_panel <- function(df,
#                              stat_labels,
#                              colors,
#                              linetypes,
#                              shapes,
#                              y_lab,
#                              x_lab = expression(pi[s]),
#                              yintercept = 0,
#                              interval = c("ci95", "se"),
#                              zero_floor = FALSE) {
#     
#     interval <- match.arg(interval)
#     mult <- if (interval == "ci95") 1.96 else 1
#     
#     p <- ggplot(
#       df,
#       aes(x = Pi_s, y = value,
#           colour = stat, linetype = stat, shape = stat, group = stat)
#     ) +
#       geom_hline(yintercept = yintercept, linetype = "dashed", colour = "grey60") +
#       geom_errorbar(
#         aes(ymin = value - mult * se,
#             ymax = value + mult * se),
#         width = 0.03,
#         position = position_dodge(0.015)
#       ) +
#       geom_line(position = position_dodge(0.015), linewidth = 0.8) +
#       geom_point(position = position_dodge(0.015), size = 2) +
#       facet_grid(cols = vars(Fst_lab), labeller = label_parsed) +
#       scale_x_continuous(breaks = sort(unique(df$Pi_s))) +
#       scale_colour_manual(values = colors, labels = stat_labels, name = NULL) +
#       scale_linetype_manual(values = linetypes, labels = stat_labels, name = NULL) +
#       scale_shape_manual(values = shapes, labels = stat_labels, name = NULL) +
#       labs(x = x_lab, y = y_lab) +
#       theme_bw(base_size = 12) +
#       theme(
#         strip.background = element_rect(fill = "grey90"),
#         strip.text = element_text(face = "bold"),
#         axis.text.x = element_text(angle = 90, vjust = 0.5),
#         legend.position = "right",
#         panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2)
#       ) +
#       scale_y_continuous(expand = expansion(mult = c(0.05, 0.10)))
#     
#     # enforce non-negative scale when needed
#     if (isTRUE(zero_floor)) {
#       p <- p + coord_cartesian(ylim = c(0, NA))
#     }
#     
#     p
#   }
#   # ============================================================
#   # 8. Panels
#   # ============================================================
#   
#   theta_labels <- c(expression(hat(theta)[1]), expression(hat(theta)[2]))
#   theta_cols   <- c("#1b9e77", "#d95f02")
#   theta_types  <- c("solid", "dashed")
#   theta_shapes <- c(16, 17)
#   
#   pA <- plot_sim_panel(
#     theta_est_plot,
#     stat_labels = theta_labels,
#     colors = theta_cols,
#     linetypes = theta_types,
#     shapes = theta_shapes,
#     y_lab = expression(hat(theta) ~ "\u00B1" ~ "95% CI"),
#     interval = "ci95"
#   )
#   
#   pC <- plot_sim_panel(
#     theta_type1_plot,
#     stat_labels = theta_labels,
#     colors = theta_cols,
#     linetypes = theta_types,
#     shapes = theta_shapes,
#     y_lab = expression("Type I error (" * hat(theta) ~ "\u00B1" ~ SE * ")"),
#     yintercept = 0.05,
#     interval = "se",
#     zero_floor = TRUE
#   )
#   
#   theta2int_labels <- c(expression(hat(theta)[2]^int))
#   theta2int_cols   <- "#7570b3"
#   theta2int_types  <- "solid"
#   theta2int_shapes <- 16
#   
#   pB <- plot_sim_panel(
#     theta2int_est_plot,
#     stat_labels = theta2int_labels,
#     colors = theta2int_cols,
#     linetypes = theta2int_types,
#     shapes = theta2int_shapes,
#     y_lab = expression(hat(theta)[2]^int ~ "\u00B1" ~ "95% CI"),
#     interval = "ci95"
#   )
#   
#   pD <- plot_sim_panel(
#     theta2int_type1_plot,
#     stat_labels = theta2int_labels,
#     colors = theta2int_cols,
#     linetypes = theta2int_types,
#     shapes = theta2int_shapes,
#     y_lab = expression("Type I error (" * hat(theta)[2]^int ~ "\u00B1" ~ SE * ")"),
#     yintercept = 0.05,
#     interval = "se",
#     zero_floor = TRUE
#   )
#   
#   # ============================================================
#   # 9. Final figure
#   # ============================================================
#   
#   final_plot <- ((pA | pC) / (pB | pD)) +
#     patchwork::plot_layout(guides = "collect") +
#     patchwork::plot_annotation(
#       tag_levels = "a",
#       tag_prefix = "(",
#       tag_suffix = ")"
#     ) &
#     theme(legend.position = "right")
#   
# }
# 
# print(final_plot)
