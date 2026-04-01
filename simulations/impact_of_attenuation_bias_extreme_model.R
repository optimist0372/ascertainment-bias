##Simulation to assess the robustness of theta1 and theta2 to attenuation bias
## in the reference following deconvolution (extreme model) ##


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(gridExtra)
})



Sim_deconc_attenuation_bias_extreme <- function(M, Pi_s, Fst, N_s, N1_r, N2_r ,C=0.01){
  
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
  
  # 7. Sample ancestral proportion uncorrected for atteniation bias
  Pi_s_hat_unc <- coef(nnls(cbind(p1_r2, p2_r2), p_s2))[1]
  
  
  
  
  # 8. Sample ancestral proportion *corrected* for atteniation bias
  Cx <- cov( cbind ( p1_r2 , p2_r2 ) )
  Cy <- cov( cbind ( p1_r2 , p2_r2 ), p_s2 )
  s  <- c( mean(p1_r2 * ( 1 - p1_r2 )) / N1_r, mean( p2_r2* ( 1- p2_r2 )) / N2_r)
  
  A <- Cx - 0.5 * diag(s)
  Y <- Cy[, 1]
  
  fit_cor <- nnls(A, Y)
  weights_cor <- coef(fit_cor)
  
  Pi_s_hat_cor <- weights_cor[1]
  
  # 9. Weighted reference using estimated uncorrected and corrected ancestral proportions
  p_r_unc   <- Pi_s_hat_unc * p1_r + (1-Pi_s_hat_unc) * p2_r
  p_r_cor   <- Pi_s_hat_cor * p1_r + (1-Pi_s_hat_cor) * p2_r
  
  
  # 10. Theta1 estimation: 
  #uncorrected
  u_unc   <- sum( 2*b * (p_s - p_r_unc) )
  v_unc   <- sum(b^2 * 2 * p_r_unc * (1 - p_r_unc))
  theta1_unc <- u_unc / sqrt(v_unc)
  
  #corrected
  u_cor   <- sum( 2*b * (p_s - p_r_cor) )
  v_cor   <- sum(b^2 * 2 * p_r_cor * (1 - p_r_cor))
  theta1_cor <- u_cor / sqrt(v_cor)
  
  #────────────────────────────────────────────
  # 11. Variance of theta1 
  #────────────────────────────────────────────
  #uncorrected
  var_ps <- p_s * (1 - p_s) / (2 * N_s)
  
  var_pr_unc <- (Pi_s_hat_unc^2) * p1_r * (1 - p1_r)/(2*N1_r) + ( (1-Pi_s_hat_unc)^2) * p2_r * (1 - p2_r)/(2*N2_r) 
  
  var_u_unc <- sum(4 * b^2 * (var_ps + var_pr_unc))
  var_v_unc <- sum(4 * b^4 * (var_pr_unc * (1 - 2 * p_r_unc)^2 + 2 * var_pr_unc^2))
  cov_uv_unc <- sum(-4 * b^3 * (1 - 2 * p_r_unc) * var_pr_unc)
  
  var_theta1_unc <- (1 / v_unc) * (var_u_unc + (u_unc^2 * var_v_unc) / (4 * v_unc^2) - (u_unc * cov_uv_unc / v_unc))
  sd_theta1_unc <- sqrt(var_theta1_unc)
  
  
  #corrected
  var_pr_cor <- (Pi_s_hat_cor^2) * p1_r * (1 - p1_r)/(2*N1_r) + ( (1-Pi_s_hat_cor)^2) * p2_r * (1 - p2_r)/(2*N2_r) 
  
  var_u_cor <- sum(4 * b^2 * (var_ps + var_pr_cor))
  var_v_cor <- sum(4 * b^4 * (var_pr_cor * (1 - 2 * p_r_cor)^2 + 2 * var_pr_cor^2))
  cov_uv_cor <- sum(-4 * b^3 * (1 - 2 * p_r_cor) * var_pr_cor)
  
  var_theta1_cor <- (1 / v_cor) * (var_u_cor + (u_cor^2 * var_v_cor) / (4 * v_cor^2) - (u_cor * cov_uv_cor / v_cor))
  sd_theta1_cor <- sqrt(var_theta1_cor)
  
  
  
  # 12. Hypothesis test for theta1 (Wald test):
  #uncorrected
  ts_unc <- (theta1_unc / sd_theta1_unc)^2
  pval_theta1_unc <- pchisq(ts_unc, df = 1, lower.tail = FALSE)
  
  #corrected
  ts_cor <- (theta1_cor / sd_theta1_cor)^2
  pval_theta1_cor <- pchisq(ts_cor, df = 1, lower.tail = FALSE)
  
  #
  # 10. Theta2 and I2 estimation (weighted regression)
  #uncorrected
  rhs_unc <- b * p_r_unc * (1 - p_r_unc) / sqrt(v_unc)
  lhs_unc <- p_s - p_r_unc
  iv_weight_unc <- 1 / var_pr_unc
  
  
  
  fit_unc <- lm(lhs_unc ~ rhs_unc, weights = iv_weight_unc)
  coef_tab_unc <- summary(fit_unc)$coefficients
  
  theta2_unc <- coef_tab_unc[2, 1]
  pval_theta2_unc <- coef_tab_unc[2, 4]
  
  I2_unc <- coef_tab_unc[1, 1]
  pval_I2_unc <- coef_tab_unc[1, 4]
  
  #corrected
  rhs_cor <- b * p_r_cor * (1 - p_r_cor) / sqrt(v_cor)
  lhs_cor <- p_s - p_r_cor
  iv_weight_cor <- 1 / var_pr_cor
  
  fit_cor <- lm(lhs_cor ~ rhs_cor, weights = iv_weight_cor)
  coef_tab_cor <- summary(fit_cor)$coefficients
  
  theta2_cor <- coef_tab_cor[2, 1]
  pval_theta2_cor <- coef_tab_cor[2, 4]
  
  I2_cor <- coef_tab_cor[1, 1]
  pval_I2_cor <- coef_tab_cor[1, 4]
  
  #────────────────────────────────────────────
  # 12. Output summary statistics
  #────────────────────────────────────────────
  
  stats <- c(theta1_unc     = theta1_unc,
             p_theta1_unc   = pval_theta1_unc,
             
             theta1_cor     = theta1_cor,
             p_theta1_cor   = pval_theta2_cor,
             
             theta2_unc     = theta2_unc,
             p_theta2_unc   = pval_theta2_unc,
             I2_unc = I2_unc,
             p_I2_unc = pval_I2_unc,
             
             theta2_cor     = theta2_cor,
             p_theta2_cor   = pval_theta2_cor,
             I2_cor = I2_cor,
             p_I2_cor = pval_I2_cor,
             Pi_s_unc      = Pi_s_hat_unc,
             Pi_s_cor      = Pi_s_hat_cor)
}




M=1000
N_s <- 5000
N1_r <- 1000
N2_r <- c(1000, 500, 100)
Fst=c(0.001,0.01)
Pi_s <- seq(0.1,0.9 ,by=0.2)

design=expand.grid(M=M,Pi_s=Pi_s, Fst=Fst,N_s=N_s,N1_r=N1_r,N2_r=N2_r )

nreps  <-10000

Design <- do.call("rbind",lapply(1:nreps,function(k) cbind(Rep=k,design)))


pb <- txtProgressBar(min = 0, max = nrow(Design), style = 3)


simulations <- do.call("rbind", lapply(1:nrow(Design), function(i) {
  setTxtProgressBar(pb, i)  
  Sim_deconc_attenuation_bias_extreme(M=Design[i, "M"],Fst=Design[i, "Fst"], Pi_s=Design[i, "Pi_s"],  
                                      N_s=Design[i, "N_s"], N1_r=Design[i, "N1_r"], N2_r=Design[i, "N2_r"])
  
}))

close(pb)

###summarize result across replicates and plot

# {
#   ## ============================================================
#   sim_dt <- as.data.table(simulations)
#   
#   # setnames(
#   #   sim_dt,
#   #   c(
#   #     "theta1_unc",       "p_theta1_unc",
#   #     "theta1_cor",       "p_theta1_cor",
#   #     "theta2_unc",       "p_theta2_unc",
#   #     "theta2_int_unc",   "p_theta2_int_unc",
#   #     "theta2_cor",       "p_theta2_cor",
#   #     "theta2_int_cor",   "p_theta2_int_cor",
#   #     "pi_s_unc",         "pi_s_cor"
#   #   )
#   # )
#   
#   design_dt <- as.data.table(Design)
#   
#   ## ============================================================
#   ## 1. Summary for Pi_s
#   ## ============================================================
#   Pi_dt <- cbind(
#     design_dt[, .(Pi_s, N2_r, Fst)],
#     sim_dt[, .(Pi_s_unc, Pi_s_cor)]
#   )
#   
#   Pi_summary <- Pi_dt[
#     ,
#     .(
#       Pi_s_unc_mean = mean(Pi_s_unc),
#       Pi_s_cor_mean = mean(Pi_s_cor),
#       Pi_s_unc_sd   = sd(Pi_s_unc),
#       Pi_s_cor_sd   = sd(Pi_s_cor)
#     ),
#     by = .(Pi_s, N2_r, Fst)
#   ]
#   
#   Pi_summary[, N2_r := factor(N2_r, levels = sort(unique(N2_r)))]
#   Pi_summary[, Fst_lab := factor(
#     ifelse(Fst == 0.001, "F[ST]==0.001", "F[ST]==0.01"),
#     levels = c("F[ST]==0.001", "F[ST]==0.01")
#   )]
#   
#   ## ============================================================
#   ## 2. Helper to summarise Type I error for one p-value column
#   ## ============================================================
#   make_type1_summary <- function(design_dt, sim_dt, p_col) {
#     
#     tmp <- cbind(
#       design_dt[, .(Pi_s, N2_r, Fst)],
#       pval = sim_dt[[p_col]]
#     )
#     
#     out <- tmp[
#       ,
#       .(value = mean(pval < 0.05)),
#       by = .(Pi_s, N2_r, Fst)
#     ]
#     
#     out[, se := sqrt(value * (1 - value) / nreps)]
#     
#     n2_levels <- sort(unique(out$N2_r), decreasing = TRUE)
#     out[, n2_lab := factor(
#       sprintf("N[2]^('r')==%s", N2_r),
#       levels = sprintf("N[2]^('r')==%s", n2_levels)
#     )]
#     
#     fst_vals <- sort(unique(out$Fst))
#     fst_cols <- setNames(
#       c("#1b9e77", "#d95f02")[seq_along(fst_vals)],
#       as.character(fst_vals)
#     )
#     
#     list(data = out, fst_cols = fst_cols)
#   }
#   
#   ## ============================================================
#   ## 3. Build Type I error summaries
#   ## ============================================================
#   theta1_unc_obj <- make_type1_summary(
#     design_dt = design_dt,
#     sim_dt    = sim_dt,
#     p_col     = "p_theta1_unc"
#   )
#   
#   theta1_cor_obj <- make_type1_summary(
#     design_dt = design_dt,
#     sim_dt    = sim_dt,
#     p_col     = "p_theta1_cor"
#   )
#   
#   theta2_unc_obj <- make_type1_summary(
#     design_dt = design_dt,
#     sim_dt    = sim_dt,
#     p_col     = "p_theta2_unc"
#   )
#   
#   theta2_cor_obj <- make_type1_summary(
#     design_dt = design_dt,
#     sim_dt    = sim_dt,
#     p_col     = "p_theta2_cor"
#   )
#   
#   ## ============================================================
#   ## 4. Plot helper for Type I error
#   ## ============================================================
#   make_type1_plot <- function(obj, y_lab_expr, ylim_max, correction_label) {
#     
#     ggplot(
#       obj$data,
#       aes(
#         x = Pi_s,
#         y = value,
#         colour = factor(Fst),
#         group = Fst
#       )
#     ) +
#       geom_hline(yintercept = 0.05, linetype = "dashed", colour = "red") +
#       geom_errorbar(
#         aes(
#           ymin = value - se,
#           ymax = value + se
#         ),
#         width = 0.02,
#         linewidth = 0.4
#       ) +
#       geom_line(linewidth = 0.7) +
#       geom_point(size = 1.5) +
#       facet_grid(
#         . ~ n2_lab,
#         labeller = labeller(n2_lab = label_parsed)
#       ) +
#       scale_colour_manual(
#         name = expression(F[ST]),
#         values = obj$fst_cols
#       ) +
#       scale_x_continuous(
#         limits = c(0, 1),
#         breaks = seq(0.1, 0.9, by = 0.2)
#       ) +
#       coord_cartesian(ylim = c(0, ylim_max)) +
#       labs(
#         x = expression(pi[s]),
#         y = y_lab_expr,
#         title = correction_label
#       ) +
#       theme_bw(base_size = 10) +
#       theme(
#         legend.position = "right",
#         strip.background = element_rect(fill = "grey90"),
#         plot.title = element_text(face = "bold", hjust = 0.5, size = 11)
#       )
#   }
#   
#   ## ============================================================
#   ## 5. Top row plots: Pi_s
#   ## ============================================================
#   p_Pi_unc <- ggplot(
#     Pi_summary,
#     aes(
#       x = Pi_s,
#       y = Pi_s_unc_mean,
#       ymin = Pi_s_unc_mean - Pi_s_unc_sd,
#       ymax = Pi_s_unc_mean + Pi_s_unc_sd,
#       colour = N2_r,
#       group = N2_r
#     )
#   ) +
#     geom_abline(slope = 1, intercept = 0, colour = "grey40") +
#     geom_errorbar(width = 0.03, linewidth = 0.5) +
#     geom_line(linewidth = 0.8) +
#     geom_point(size = 1.8) +
#     facet_wrap(~ Fst_lab, nrow = 1, labeller = label_parsed) +
#     scale_x_continuous(
#       limits = c(0, 1),
#       breaks = seq(0.1, 0.9, by = 0.2)
#     ) +
#     scale_y_continuous(
#       limits = c(0, 1),
#       breaks = seq(0.1, 0.9, by = 0.2)
#     ) +
#     scale_colour_brewer(
#       palette = "Dark2",
#       name = expression(N[2]^"(r)")
#     ) +
#     labs(
#       x = expression(True ~ pi[s]),
#       y = expression(hat(pi)[s] %+-% sd),
#       title = "No correction"
#     ) +
#     theme_bw(base_size = 11) +
#     theme(
#       legend.position = "right",
#       plot.title = element_text(face = "bold", hjust = 0.5, size = 11)
#     )
#   
#   p_Pi_cor <- ggplot(
#     Pi_summary,
#     aes(
#       x = Pi_s,
#       y = Pi_s_cor_mean,
#       ymin = Pi_s_cor_mean - Pi_s_cor_sd,
#       ymax = Pi_s_cor_mean + Pi_s_cor_sd,
#       colour = N2_r,
#       group = N2_r
#     )
#   ) +
#     geom_abline(slope = 1, intercept = 0, colour = "grey40") +
#     geom_errorbar(width = 0.03, linewidth = 0.5) +
#     geom_line(linewidth = 0.8) +
#     geom_point(size = 1.8) +
#     facet_wrap(~ Fst_lab, nrow = 1, labeller = label_parsed) +
#     scale_x_continuous(
#       limits = c(0, 1),
#       breaks = seq(0.1, 0.9, by = 0.2)
#     ) +
#     scale_y_continuous(
#       limits = c(0, 1),
#       breaks = seq(0.1, 0.9, by = 0.2)
#     ) +
#     scale_colour_brewer(
#       palette = "Dark2",
#       name = expression(N[2]^"(r)")
#     ) +
#     labs(
#       x = expression(True ~ pi[s]),
#       y = expression(hat(pi)[s] %+-% sd),
#       title = "Corrected"
#     ) +
#     theme_bw(base_size = 11) +
#     theme(
#       legend.position = "right",
#       plot.title = element_text(face = "bold", hjust = 0.5, size = 11)
#     )
#   
#   ## ============================================================
#   ## 6. Middle row plots: theta1 Type I error
#   ## ============================================================
#   p_theta1_unc <- make_type1_plot(
#     theta1_unc_obj,
#     expression("Type I error (" * hat(theta)[1] %+-% SE * ")"),
#     ylim_max = 1,
#     correction_label = "No correction"
#   )
#   
#   p_theta1_cor <- make_type1_plot(
#     theta1_cor_obj,
#     expression("Type I error (" * hat(theta)[1] %+-% SE * ")"),
#     ylim_max = 1,
#     correction_label = "Corrected"
#   )
#   
#   ## ============================================================
#   ## 7. Bottom row plots: theta2 Type I error
#   ## ============================================================
#   p_theta2_unc <- make_type1_plot(
#     theta2_unc_obj,
#     expression("Type I error (" * hat(theta)[2] %+-% SE * ")"),
#     ylim_max = 1,
#     correction_label = "No correction"
#   )
#   
#   p_theta2_cor <- make_type1_plot(
#     theta2_cor_obj,
#     expression("Type I error (" * hat(theta)[2] %+-% SE * ")"),
#     ylim_max = 1,
#     correction_label = "Corrected"
#   )
#   
#   ## ============================================================
#   ## 8. Add panel tags
#   ## ============================================================
#   add_tag <- function(p, tag_label) {
#     p + labs(tag = tag_label) +
#       theme(
#         plot.tag = element_text(face = "bold", size = 14),
#         plot.tag.position = c(0.02, 0.98)
#       )
#   }
#   
#   p_Pi_unc     <- add_tag(p_Pi_unc, "(a)")
#   p_Pi_cor     <- add_tag(p_Pi_cor, "(b)")
#   p_theta1_unc <- add_tag(p_theta1_unc, "(c)")
#   p_theta1_cor <- add_tag(p_theta1_cor, "(d)")
#   p_theta2_unc <- add_tag(p_theta2_unc, "(e)")
#   p_theta2_cor <- add_tag(p_theta2_cor, "(f)")
#   
#   ## ============================================================
#   ## 9. Arrange final figure
#   ## ============================================================
#   row1 <- grid.arrange(p_Pi_unc,     p_Pi_cor,     ncol = 2)
#   row2 <- grid.arrange(p_theta1_unc, p_theta1_cor, ncol = 2)
#   row3 <- grid.arrange(p_theta2_unc, p_theta2_cor, ncol = 2)
#   
#   grid.arrange(row1, row2, row3, nrow = 3)
# }
