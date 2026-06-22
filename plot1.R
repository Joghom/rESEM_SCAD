# ============================================================
# LOAD ALL OUTPUT FILES AND CREATE COMBINED FACET BOXPLOTS
# ============================================================

rm(list = ls())

library(ggplot2)
library(dplyr)

# ------------------------------------------------------------
# SET WORKING DIRECTORY
# ------------------------------------------------------------

current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

# ------------------------------------------------------------
# LOAD SIMULATION INFO
# ------------------------------------------------------------

load("DATA/Info_simulation.RData")

Info_matrix <- Infor_simulation$design_matrix_replication
Ndatasets   <- Infor_simulation$n_data_sets

# ------------------------------------------------------------
# FUNCTION TO LOAD RESULTS
# ------------------------------------------------------------

load_results <- function(folder, prefix, method_name){
  
  results_all <- list()
  
  for(i in 1:Ndatasets){
    
    # load result
    load(paste0(folder, "/", prefix, i, ".RData"))
    
    # --------------------------------------------------------
    # SELECT OBJECT DEPENDING ON METHOD
    # --------------------------------------------------------
    
    if(method_name == "SCAD"){
      measures <- as.data.frame(output_SCAD$measure_summary)
    }
    
    if(method_name == "rESEMCC"){
      measures <- as.data.frame(output$measure_summary)
    }
    
    if (method_name == 'Lasso') {
      measures <- as.data.frame(output_l1$measure_summary)
    }
    
    # --------------------------------------------------------
    # ADD CONDITIONS
    # --------------------------------------------------------
    
    measures$dataset     <- i
    measures$N           <- Info_matrix$sample[i]
    measures$reliability <- Info_matrix$reliability[i]
    measures$nfactors    <- Info_matrix$nfactors[i]
    measures$factor_corr <- Info_matrix$factor_corr[i]
    measures$CL          <- Info_matrix$CL[i]
    
    measures$method <- method_name
    
    results_all[[i]] <- measures
  }
  
  bind_rows(results_all)
}

# ------------------------------------------------------------
# LOAD BOTH METHODS
# ------------------------------------------------------------

results_SCAD <- load_results(
  folder = "SCAD_out1",
  prefix = "scad",
  method_name = "SCAD"
)

results_rESEMCC <- load_results(
  folder = "Complete_rESEMCC_out",
  prefix = "rESEM",
  method_name = "rESEMCC"
)

results_lasso <- load_results(
  folder = "Complete_rESEMl1_out",
  prefix = 'Lasso',
  method_name = 'Lasso'
)

# ------------------------------------------------------------
# COMBINE
# ------------------------------------------------------------

results_df <- bind_rows(results_SCAD, results_rESEMCC, results_lasso)

# ------------------------------------------------------------
# FACTORS
# ------------------------------------------------------------

results_df$method <- factor(
  results_df$method,
  levels = c("SCAD", "rESEMCC", 'Lasso')
)

results_df$N <- factor(results_df$N)

results_df$reliability <- factor(results_df$reliability)

results_df$nfactors <- factor(results_df$nfactors)

results_df$factor_corr <- factor(results_df$factor_corr)

results_df$CL <- factor(results_df$CL)

results_uncorrelated <- subset(results_df, factor_corr == 0)
results_correlated <- subset(results_df, factor_corr == 0.3)

# ============================================================
# PLOT 1: ZERO/NONZERO RECOVERY
# ============================================================

p_zero_corr <- ggplot(
  results_correlated,
  aes(
    x = method,
    y = zero_rec,
    fill = CL
  )
) +
  
  ggtitle('Correlated Zero/Nonzero recovery rate') +
  
  geom_boxplot(
    width = 0.6,
    outlier.size = 0.5
  ) +
  
  geom_hline(
    yintercept = 0.85,
    linetype = "dashed"
  ) +
  
  facet_grid(
    rows = vars(N),
    cols = vars(reliability, nfactors),
    labeller = label_both
  ) +
  
  coord_cartesian(ylim = c(0, 1.05)) +
  
  labs(
    x = "Method",
    y = "Zero/nonzero loadings recovery rate"
  ) +
  
  theme_bw(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    
    strip.background = element_rect(
      fill = "grey85",
      colour = "black"
    ),
    
    panel.grid.major = element_line(
      colour = "grey90"
    ),
    
    panel.grid.minor = element_blank(),
    
    legend.position = "bottom"
  )

print(p_zero_corr)
ggsave('p_zero_corr.PNG', plot = p_zero_corr, path = 'plots/', width=6.76, height = 6.76)


p_zero_uncorr <- ggplot(
  results_uncorrelated,
  aes(
    x = method,
    y = zero_rec,
    fill = CL
  )
) +
  
  ggtitle('Uncorrelated Zero/Nonzero recovery rate') +
  
  geom_boxplot(
    width = 0.6,
    outlier.size = 0.5
  ) +
  
  geom_hline(
    yintercept = 0.85,
    linetype = "dashed"
  ) +
  
  facet_grid(
    rows = vars(N),
    cols = vars(reliability, nfactors),
    labeller = label_both
  ) +
  
  coord_cartesian(ylim = c(0, 1.05)) +
  
  labs(
    x = "Method",
    y = "Zero/nonzero loadings recovery rate"
  ) +
  
  theme_bw(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    
    strip.background = element_rect(
      fill = "grey85",
      colour = "black"
    ),
    
    panel.grid.major = element_line(
      colour = "grey90"
    ),
    
    panel.grid.minor = element_blank(),
    
    legend.position = "bottom"
  )

print(p_zero_uncorr)
ggsave('p_zero_uncorr.PNG', plot = p_zero_uncorr, path = 'plots/', width=6.76, height = 6.76)

# ============================================================
# PLOT 2: TUCKER CONGRUENCE
# ============================================================

p_tucker_corr <- ggplot(
  results_correlated,
  aes(
    x = method,
    y = tucker,
    fill = CL
  )
) +
  
  ggtitle('Correlated Tucker Congruence') +
  
  geom_boxplot(
    width = 0.6,
    outlier.size = 0.5
  ) +
  
  geom_hline(
    yintercept = 0.95,
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = 0.85,
    linetype = "dotted"
  ) +
  
  facet_grid(
    rows = vars(N),
    cols = vars(reliability, nfactors),
    labeller = label_both
  ) +
  
  coord_cartesian(ylim = c(0, 1.05)) +
  
  labs(
    x = "Method",
    y = "Tucker congruence"
  ) +
  
  theme_bw(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    
    strip.background = element_rect(
      fill = "grey85",
      colour = "black"
    ),
    
    panel.grid.major = element_line(
      colour = "grey90"
    ),
    
    panel.grid.minor = element_blank(),
    
    legend.position = "bottom"
  )

print(p_tucker_corr)
ggsave('p_tucker_corr.PNG', plot = p_tucker_corr, path = 'plots/', width=6.76, height = 6.76)


p_tucker_uncorr <- ggplot(
  results_uncorrelated,
  aes(
    x = method,
    y = tucker,
    fill = CL
  )
) +
  
  ggtitle('Uncorrelated Tucker Congruence') +
  
  geom_boxplot(
    width = 0.6,
    outlier.size = 0.5
  ) +
  
  geom_hline(
    yintercept = 0.95,
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = 0.85,
    linetype = "dotted"
  ) +
  
  facet_grid(
    rows = vars(N),
    cols = vars(reliability, nfactors),
    labeller = label_both
  ) +
  
  coord_cartesian(ylim = c(0, 1.05)) +
  
  labs(
    x = "Method",
    y = "Tucker congruence"
  ) +
  
  theme_bw(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    
    strip.background = element_rect(
      fill = "grey85",
      colour = "black"
    ),
    
    panel.grid.major = element_line(
      colour = "grey90"
    ),
    
    panel.grid.minor = element_blank(),
    
    legend.position = "bottom"
  )

print(p_tucker_uncorr)
ggsave('p_tucker_uncorr.PNG', plot = p_tucker_uncorr, path = 'plots/', width=6.76, height = 6.76)