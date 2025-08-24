library(figsimR)
library(dplyr)
library(forcats)
library(ggplot2)
library(scales)
library(DiagrammeR)
library(vegan)
library(patchwork)
library(purrr)
library(stringr)
library(tibble)
library(readr)
library(tidyr)

# Set custom font (Arial) for figure output
# This section ensures consistent font rendering across figures
# NOTE: This step may be platform-specific (Windows). Users on macOS/Linux may need to adjust font path or use a different font.
library(extrafont)
library(showtext)

list.files("C:/Windows/Fonts", pattern = "\\.ttf$", full.names = TRUE)  # optional: inspect available fonts

font_import(paths = "C:/Windows/Fonts", pattern = "arial", prompt = FALSE)  # import Arial font (Windows only)
loadfonts(device = "win")


font_add("Arial", regular = "C:/Windows/Fonts/arial.ttf")  # register Arial font for showtext
showtext_auto()  # enable showtext for all plots

####################################### Figure 1 #########################################################


dot_graph <- "
digraph figsimR_architecture {

  // ---- Global ----
  graph [layout=dot, rankdir=TB, splines=ortho, fontname=Helvetica, bgcolor=\"#FFFFFF\"];
  node  [style=filled, fontname=Helvetica, shape=box, fontsize=11, color=gray60];
  edge  [color=gray40, arrowsize=0.8, fontname=Helvetica, fontsize=10];

  // ---- Inputs ----
  subgraph cluster_inputs {
    label = \"1. Inputs\"; labelloc = t; fontsize = 16; fontcolor = \"#37474F\";
    style = \"rounded,dashed\"; penwidth = 1.2; bgcolor = \"#FAFAFA\"; color = \"#90A4AE\";
    params  [label = \"Biological Parameters\",        fillcolor = \"#ECEFF1\", color = \"#546E7A\"];
    species [label = \"Species Roles & Interactions\", fillcolor = \"#ECEFF1\", color = \"#546E7A\"];
    obs     [label = \"Empirical Observed Data\",      fillcolor = \"#ECEFF1\", color = \"#546E7A\"];
  }

  // ---- Core Simulation (parallel modules acting on a shared state) ----
  subgraph cluster_process {
    label = \"Core Simulation Engine\"; labelloc = t; fontsize = 16; fontcolor = \"#37474F\";
    style = \"rounded,dashed\"; penwidth = 1.2; bgcolor = \"#FAFAFA\"; color = \"#90A4AE\";
    state  [label = \"Per-fig\\nSimulation State\", shape=doublecircle,  fillcolor = \"#E3F2FD\", color = \"#1E88E5\",
             fontsize = 13, fixedsize = true, width = 1.3];
    init   [label = \"Initialize\\nFig Environment\", shape=ellipse, fillcolor = \"#CFD8DC\"];
    entry  [label = \"Arrival &\\nEntry\",            shape=ellipse, fillcolor = \"#CFD8DC\"];
    ovip   [label = \"Oviposition &\\nCompetition\",  shape=ellipse, fillcolor = \"#CFD8DC\"];
    drop   [label = \"Host\\nSanctions\",             shape=ellipse, fillcolor = \"#CFD8DC\"];
    layer  [label = \"Spatial\\nPartitioning\",       shape=ellipse, fillcolor = \"#CFD8DC\"];
    para   [label = \"Parasitism\",                   shape=ellipse, fillcolor = \"#CFD8DC\"];
    emerge [label = \"Wasp Emergence\",               shape=cds,     fillcolor = \"#B3E5FC\"];
  }

  // ---- Output & Diagnosis ----
  subgraph cluster_output {
    label = \"3. Output & Diagnosis\"; labelloc = t; fontsize = 16; fontcolor = \"#37474F\";
    style = \"rounded,dashed\"; penwidth = 1.2; bgcolor = \"#FAFAFA\"; color = \"#90A4AE\";
    sim_data [label = \"Simulated Community Data\",                 shape=rounded, fillcolor = \"#E8F5E9\", color = \"#43A047\"];
    metrics  [label = \"Calculate Metrics &\\nQuantify Variation Gap\", shape=rounded, fillcolor = \"#E8F5E9\", color = \"#43A047\"];
  }

  // ---- Key Modular Controls (knockout targets) ----
  subgraph cluster_controls {
    label = \"Key Modular Controls (Knockout Targets)\"; labelloc = t; fontsize = 14; fontcolor = \"#37474F\";
    style = \"rounded,dashed\"; penwidth = 1.2; bgcolor = \"#FFFDE7\"; color = \"#FBC02D\";
    c_prio  [label = \"Priority Effects\",        shape=box, fillcolor = \"#FFF9C4\", color = \"#FBC02D\", fontsize=12];
    c_phase [label = \"Phase-specific Success\",  shape=box, fillcolor = \"#FFF9C4\", color = \"#FBC02D\", fontsize=12];
    c_drop  [label = \"Host Sanctions\",          shape=box, fillcolor = \"#FFF9C4\", color = \"#FBC02D\", fontsize=12];
    c_layer [label = \"Spatial Partitioning\",    shape=box, fillcolor = \"#FFF9C4\", color = \"#FBC02D\", fontsize=12];
    c_para  [label = \"Parasitism Probability\",  shape=box, fillcolor = \"#FFF9C4\", color = \"#FBC02D\", fontsize=12];
  }

  // ---- Edges ----
  edge [style=bold, arrowhead=normal];
  params  -> init;
  init    -> state;
  state   -> emerge;
  emerge  -> sim_data;

  species -> entry [style=bold, arrowhead=normal];

  edge [style=bold, arrowhead=none];
  entry -> state  [style=dashed];
  ovip  -> state  [style=dashed];
  drop  -> state  [style=dashed];
  layer -> state  [style=dashed];
  para  -> state  [style=dashed];

  edge [style=solid, arrowhead=normal, color=gray40];
  sim_data -> metrics;
  obs      -> metrics [style=dotted];

  metrics -> params [
    style=dashed, arrowhead=normal, color=\"#0D47A1\",
    label=\"  Calibration Loop\\n(LHS / Optimization)  \",
    fontcolor=\"#0D47A1\", fontsize=12, constraint=false
  ];

  c_prio  -> entry [style=dashed, arrowhead=odot, color=\"#F57F17\", penwidth=1.2, constraint=false];
  c_phase -> ovip  [style=dashed, arrowhead=odot, color=\"#F57F17\", penwidth=1.2, constraint=false];
  c_drop  -> drop  [style=dashed, arrowhead=odot, color=\"#F57F17\", penwidth=1.2, constraint=false];
  c_layer -> layer [style=dashed, arrowhead=odot, color=\"#F57F17\", penwidth=1.2, constraint=false];
  c_para  -> para  [style=dashed, arrowhead=odot, color=\"#F57F17\", penwidth=1.2, constraint=false];

  { rank=same; params; species; obs }
  c_prio  -> entry [style=invis, weight=10];
  c_phase -> ovip  [style=invis, weight=10];
  c_drop  -> drop  [style=invis, weight=10];
  c_layer -> layer [style=invis, weight=10];
  c_para  -> para  [style=invis, weight=10];
  { rank=same; entry; ovip; drop; layer; para }
  { rank=same; sim_data; metrics }
}
"


grViz(dot_graph)

g <- DiagrammeR::grViz(dot_graph)

svg <- DiagrammeRsvg::export_svg(g)
rsvg::rsvg_pdf(charToRaw(svg), "fig1_architecture.pdf")


###################### 
set.seed(42)

data("observed_data")
data("species_list")
data("parameter_list_default")

file_path <- system.file("extdata", "lhs_optimize_result.rds", package = "figsimR")
lhs_optimize_result <- readRDS(file_path)

file_path <- system.file("extdata", "simulated_metric_distribution.rds", package = "figsimR")
resample_simulated_metrics_df_results <- readRDS(file_path)

# Rebuild the full parameter list object from the best result row
lhs_optimize_parameter <- rebuild_parameter_list_from_row(lhs_optimize_result[1, ],
                                                                 parameter_list_default)

sim_df <- simulate_figwasp_community(
  num_figs                = 1000,
  fecundity_mean          = lhs_optimize_parameter$fecundity_mean,
  fecundity_dispersion    = lhs_optimize_parameter$fecundity_dispersion,
  entry_mu                = lhs_optimize_parameter$entry_mu,
  entry_size              = lhs_optimize_parameter$entry_size,
  entry_priority          = lhs_optimize_parameter$entry_priority,
  species_roles           = lhs_optimize_parameter$species_roles,
  max_entry_table         = lhs_optimize_parameter$max_entry_table,
  enable_drop             = TRUE,
  record_individual       = FALSE,
  drop_cancels_emergence  = FALSE,
  entry_distribution      = "lognormal",
  interaction_matrix      = lhs_optimize_parameter$interaction_matrix,
  interaction_weight      = 0,
  egg_success_prob        = lhs_optimize_parameter$egg_success_prob,
  egg_success_prob_by_phase = lhs_optimize_parameter$egg_success_prob_by_phase,
  layer_preference        = lhs_optimize_parameter$layer_preference,
  use_layering            = FALSE,
  seed                    = 42
)

names(sim_df)

sim_df_comm <- sim_df$summary

# observed data distribution
obs_df_comm<-observed_data

obs_df <- obs_df_comm
obs_df$total <- rowSums(obs_df)
obs_df <- obs_df[obs_df$total != 0, ]



# resample observed data
resample_observed_metrics_df_results  <- resample_observed_all_metrics(
  observed_df = observed_data,       # empirical data
  seed = 42,                         # set seed for reproducibility
  wasp_cols = species_list,          # columns corresponding to species counts
  n_draws = 500,                     # number of bootstrap replicates
  sample_n = 200,                    # number of figs sampled per replicate
  calc_func = calc_all_metrics       # metric calculation function
)


# Summary statistics across all bootstrap replicates
summary(resample_observed_metrics_df_results )

head(resample_simulated_metrics_df_results)



metrics_keep <- c("mean_richness", "mean_evenness", "mean_bray_curtis",
                  "connectance", "nestedness", "modularity")

metric_labels <- c(
  mean_richness    = "Richness",
  mean_evenness    = "Pielou’s evenness",
  mean_bray_curtis = "Mean Bray–Curtis",
  connectance      = "Connectance",
  nestedness       = "Nestedness",
  modularity       = "Modularity"
)

# --- Sim bootstrap summary
sim_long <- resample_simulated_metrics_df_results %>%
  select(all_of(metrics_keep)) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "val")

sim_sum <- sim_long %>%
  group_by(metric) %>%
  summarise(
    sim_mu = mean(val),
    sim_sd = sd(val),
    q80_lo = quantile(val, 0.10),
    q80_hi = quantile(val, 0.90),
    q95_lo = quantile(val, 0.025),
    q95_hi = quantile(val, 0.975),
    .groups = "drop"
  ) %>%
  mutate(sim_sd = pmax(sim_sd, .Machine$double.eps))  # 保护

# --- Observed mean (from its bootstrap-of-means)
obs_mu <- resample_observed_metrics_df_results %>%
  select(all_of(metrics_keep)) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "val") %>%
  group_by(metric) %>%
  summarise(obs_mu = mean(val), .groups = "drop")

# --- Combine & standardize
dash <- sim_sum %>%
  inner_join(obs_mu, by = "metric") %>%
  mutate(
    z_obs  = (obs_mu - sim_mu) / sim_sd,
    z80_lo = (q80_lo - sim_mu) / sim_sd,
    z80_hi = (q80_hi - sim_mu) / sim_sd,
    z95_lo = (q95_lo - sim_mu) / sim_sd,
    z95_hi = (q95_hi - sim_mu) / sim_sd,
    band = case_when(
      z_obs < z95_lo | z_obs > z95_hi ~ ">95%",
      z_obs < z80_lo | z_obs > z80_hi ~ "80–95%",
      TRUE                            ~ "Inside 80%"
    ),
    
    metric_lab = fct_reorder(metric_labels[metric], abs(z_obs), .desc = TRUE),
    
    z_cap = scales::squish(z_obs, range = c(-3, 3)),
    capped = if_else(z_obs != z_cap, TRUE, FALSE),
    diff_raw = obs_mu - sim_mu  
  )

# SUGGESTION: Use a color-blind friendly palette and slightly more descriptive labels.
# A dark gray for the "good" points can make the colored outliers pop more.
cols <- c("Within 80% PI" = "gray20", "80-95% PI" = "#E69F00", "Outside 95% PI" = "#D55E00")
band_labels <- c("Within 80% PI", "80-95% PI", "Outside 95% PI")
dash$band_labeled <- factor(case_when(
  dash$band == "Inside 80%" ~ "Within 80% PI",
  dash$band == "80–95%" ~ "80-95% PI",
  TRUE ~ "Outside 95% PI"
), levels = band_labels)


figure_2 <- ggplot(dash, aes(y = metric_lab)) +
  # Predictive Intervals
  # 95% band (outer, lighter)
  geom_segment(aes(x = z95_lo, xend = z95_hi, yend = metric_lab),
               linewidth = 2, color = "#AED6F1", lineend = "round") +
  # 80% band (inner, darker)
  geom_segment(aes(x = z80_lo, xend = z80_hi, yend = metric_lab),
               linewidth = 4, color = "#5DADE2", lineend = "round") +
  
  # FIM Mean (reference line)
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  
  # Observed Mean points (using the capped values for plotting)
  geom_point(aes(x = z_cap, color = band_labeled), size = 4) +
  
  # Text annotation for capped outliers
  geom_text(
    data = subset(dash, capped),
    aes(
      x     = ifelse(z_obs > 0, 3.2, -3.2),
      label = sprintf("%.1f", z_obs),
      hjust = ifelse(z_obs > 0, 0, 1)
    ),
    vjust = 0.5, size = 3.5, fontface = "bold",
    color = "#D55E00" # Match the outlier point color
  ) +
  
  # Scales and Labels
  scale_color_manual(values = cols, name = "Observed Mean Position:", breaks = band_labels) +
  scale_x_continuous(limits = c(-3.5, 4.5), breaks = seq(-2, 2, 2)) +
  labs(
    #title = "Figure 2. The FIM Accurately Recovers Empirical Central Tendencies",
    #subtitle = "Most observed community metric means fall within the FIM's 80% predictive interval.",
    y = NULL,
    x = "Standardized Deviation from FIM Mean"
  ) +
  
  # Theme
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(
    legend.position = "bottom",
    plot.title.position = "plot",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30", margin = margin(b = 10)),
    axis.title.x = element_text(margin = margin(t = 5))
  )

figure_2
ggsave("Figure 2.pdf", figure_2, width = 183, height = 100, units = "mm", dpi = 500)


## ======================
## Figure 3 (A–C) — code
## ======================


## ---------- Inputs you already have ----------
## 1) Community matrices (rows = figs, cols = species)
obs_comm_mat <- obs_df_comm
sim_comm_mat <- sim_df_comm[,c(18:23)]
colnames(sim_comm_mat) <- gsub("^emergence_", "", colnames(sim_comm_mat))
sim_comm_mat <- sim_comm_mat[, match(colnames(obs_comm_mat), colnames(sim_comm_mat))]

## 2) Bootstrap tables used also for Fig.2
## resample_observed_metrics_df_results
## resample_simulated_metrics_df_results

## 3) Metric list and pretty labels

pal_obs <- "#D55E00"  # orange
pal_sim <- "#2E86DE"  # blue

## =========================================================
## Panel A — Multivariate Separation (Hellinger–PCoA)
## =========================================================
# Hellinger transform
hell_obs <- decostand(as.matrix(obs_comm_mat), method = "hellinger")
hell_sim <- decostand(as.matrix(sim_comm_mat), method = "hellinger")

# Combine for a single Euclidean geometry
X        <- rbind(hell_obs, hell_sim)
grp      <- factor(c(rep("Observed", nrow(hell_obs)),
                     rep("FIM",      nrow(hell_sim))),
                   levels = c("FIM", "Observed"))


# obs_comm_mat / sim_comm_mat: rows = figs, cols = species (counts)
dispersion_test <- function(obs_comm_mat, sim_comm_mat,
                            transform = c("hellinger", "none"),
                            permutations = 999, seed = 1) {
  transform <- match.arg(transform)
  
  # 1) Transform (Hellinger) and bind
  if (transform == "hellinger") {
    Xobs <- decostand(as.matrix(obs_comm_mat), "hellinger")
    Xsim <- decostand(as.matrix(sim_comm_mat), "hellinger")
  } else {
    Xobs <- as.matrix(obs_comm_mat)
    Xsim <- as.matrix(sim_comm_mat)
  }
  grp <- factor(c(rep("Observed", nrow(Xobs)),
                  rep("FIM",      nrow(Xsim))),
                levels = c("FIM", "Observed"))
  X <- rbind(Xobs, Xsim)
  
  # 2) Distances & betadisper
  D  <- dist(X, method = "euclidean")
  bd <- betadisper(D, group = grp)
  
  # 3) Permutation test (Pr(>F))
  set.seed(seed)
  prm <- permutest(bd, permutations = permutations)
  
  # 4) Means & ratio
  mu  <- tapply(bd$distances, grp, mean)  # mean distance to centroid per group
  rho <- unname(mu["Observed"] / mu["FIM"])
  p   <- prm$tab[1, "Pr(>F)"]
  
  list(
    rho = rho, p = p,
    means = mu,
    betadisper = bd, perm = prm,
    distances = bd$distances, group = grp   # for plotting Panel B
  )
}

# 
res <- dispersion_test(obs_comm_mat, sim_comm_mat, permutations = 999)
res$rho  # dispersion ratio: 1.331619
res$p    # permutation p-value: 0.001
res$means
# FIM  Observed 
# 0.3835673 0.5107653



#########################################
# PCoA (cmdscale on Euclidean)
D        <- dist(X, method = "euclidean")
pc        <- cmdscale(D, k = 2, eig = TRUE)
scores    <- as.data.frame(pc$points)
colnames(scores) <- c("PC1", "PC2")
scores$group     <- grp

# Variance explained (pc$eig can be negative small due to rounding; zero-bound them)
eig_pos  <- pmax(pc$eig, 0)
var_exp1 <- round(100 * eig_pos[1] / sum(eig_pos), 1)
var_exp2 <- round(100 * eig_pos[2] / sum(eig_pos), 1)

fg3_panela <- ggplot(scores, aes(PC1, PC2, color = group)) +
  geom_point(alpha = 0.55, size = 1) +
  stat_ellipse(level = 0.95, linewidth = 0.2) +
  scale_color_manual(values = c(FIM = pal_sim, Observed = pal_obs), name = NULL) +
  labs(#title = "A   Multivariate separation (PCoA on Hellinger–Euclidean)",
       x = paste0("PC1 (", var_exp1, "%)"),
       y = paste0("PC2 (", var_exp2, "%)")) +
  theme_classic(base_size = 9, base_family = "Arial") +
  theme(legend.position = "bottom")

## =========================================================
## Panel B — Multivariate dispersion (betadisper)
## =========================================================
bd   <- betadisper(D, grp)          # distances to within-group centroids
perm <- permutest(bd, permutations = 999)

# Distances for density plot
dist_df <- tibble(distance = bd$distances,
                  group    = grp)

# Effect size: ratio of mean distances
mean_dist <- dist_df %>% group_by(group) %>% summarise(mu = mean(distance), .groups = "drop")
rho <- with(mean_dist, mu[group=="Observed"] / mu[group=="FIM"])
rho_lab <- sprintf("Dispersion ratio ρ = %.2f;  perm p = %.3f",
                   rho, perm$tab[1, "Pr(>F)"])

fg3_panelb <- ggplot(dist_df, aes(x = distance, fill = group, color = group)) +
  geom_density(alpha = 0.25, adjust = 1) +
  geom_vline(data = mean_dist, aes(xintercept = mu, color = group),
             linetype = "dashed") +
  scale_fill_manual(values = c(FIM = pal_sim, Observed = pal_obs), name = NULL) +
  scale_color_manual(values = c(FIM = pal_sim, Observed = pal_obs), name = NULL) +
  labs(#title = "B   Multivariate dispersion (distance to group centroid)",
       #subtitle = rho_lab,
       x = "Distance to centroid", y = "Density") +
  theme_classic(base_size = 9, base_family = "Arial") +
  theme(legend.position = "bottom")

## =========================================================
## Panel C — Univariate coverage & quantile residuals
## =========================================================
# Helper: bootstrap a quantile from a vector v
boot_quantile <- function(v, q, B = 1000) {
  replicate(B, { quantile(sample(v, replace = TRUE), q, type = 8) })
}

# Gather replicate vectors for each metric
sim_long <- resample_simulated_metrics_df_results %>%
  select(all_of(metrics_keep)) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "val_sim")

obs_long <- resample_observed_metrics_df_results %>%
  select(all_of(metrics_keep)) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "val_obs")

sim_nst <- sim_long %>%
  group_by(metric) %>%
  summarise(v_sim = list(val_sim), .groups = "drop")

obs_nst <- obs_long %>%
  group_by(metric) %>%
  summarise(v_obs = list(val_obs), .groups = "drop")

# 
repd <- left_join(sim_nst, obs_nst, by = "metric")

# 
q_levels <- c(0.05, 0.25, 0.5, 0.75, 0.95)

# 
boot_quantile <- function(v, q, B = 1000) {
  replicate(B, quantile(sample(v, replace = TRUE), q, type = 8))
}

# 
qc <- repd %>%
  mutate(panel = map2(v_sim, v_obs, ~{
    map_dfr(q_levels, function(q) {
      bq  <- boot_quantile(.x, q, B = 1000)
      tibble(
        q      = q,
        lo95   = quantile(bq, 0.025, type = 8),
        hi95   = quantile(bq, 0.975, type = 8),
        lo80   = quantile(bq, 0.10,  type = 8),
        hi80   = quantile(bq, 0.90,  type = 8),
        simmid = median(bq),
        obs    = quantile(.y, q, type = 8)
      )
    })
  })) %>%
  select(metric, panel) %>%
  unnest(panel) %>%
  mutate(
    metric_lab = factor(metric_labels[metric], levels = metric_labels[metrics_keep]),
    status95   = if_else(obs < lo95 | obs > hi95, "Outside 95%", "Inside 95%"),
    q_lab      = factor(paste0(q*100, "%"), levels = paste0(q_levels*100, "%"))
  )

# 
qc_order <- qc %>%
  group_by(metric_lab) %>%
  summarise(max_tail = max(abs(obs - if_else(q < 0.5, lo95, hi95))), .groups = "drop") %>%
  arrange(desc(max_tail)) %>% pull(metric_lab)
qc <- qc %>% mutate(metric_lab = factor(metric_lab, levels = qc_order))

# 
fg3_panelc <- ggplot(qc, aes(x = q)) +
  # 95% 
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "grey90", alpha = 0.8) +
  # 80% 
  geom_ribbon(aes(ymin = lo80, ymax = hi80), fill = "grey75", alpha = 0.8) +
  
  # 
  geom_line(aes(y = simmid), color = "#0072B2", linetype = "dashed", linewidth = 0.6) +
  
  # 
  geom_line(aes(y = obs, group = 1), color = "#D55E00", linewidth = 1) +
  geom_point(aes(y = obs, fill = status95), shape = 21, size = 1.5, stroke = 0.6) +
  
  # 
  scale_fill_manual(
    values = c("Inside 95%" = "white", "Outside 95%" = "#D55E00"),
    name = "Observed Quantile:"
  ) +
  scale_x_continuous(
    breaks = q_levels,
    labels = scales::percent_format(accuracy = 1)
  ) +
  facet_wrap(~ metric_lab, ncol = 2, scales = "free_y") + # 改为两列，让每个面板更宽
  
  # 
  labs(
    #title = "C   Scalar Coverage and Quantile Diagnostics",
    #subtitle = "Observed quantiles (orange line) compared to the FIM's predictive intervals (grey bands)",
    x = "Quantile of the Distribution",
    y = "Metric Value"
  ) +
  theme_bw(base_size = 9, base_family = "Arial") +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30", margin = margin(b = 10))
  )



## =========================================================
## Assemble Figure 3
## =========================================================
figure_3 <- fg3_panela + fg3_panelb + fg3_panelc + plot_layout(widths = c(1, 1, 1.2)) +
  plot_annotation(tag_levels = "A", tag_prefix = "", tag_suffix = "", 
                  theme = theme(plot.tag = element_text(face = "bold", size = 10)))
figure_3

ggsave("Figure 3.pdf", figure_3, width = 240, height = 100, units = "mm", dpi = 500)

############ Figure 4 ###################
######################## Experiment III: Dissecting the Importance of Intrinsic Mechanisms via Perturbation Analysis ########################

# Run a set of reduced model configurations to test the necessity of specific intrinsic mechanisms.
# Each entry modifies one key mechanism in the Full Intrinsic Model (FIM) while keeping other parameters fixed.

fim_reduced_results <- list(
  run_fim_reduced_experiment("FIM", lhs_optimize_parameter, species_list),  # Remove spatial layering
  run_fim_reduced_experiment("FIM_No_space", lhs_optimize_parameter, species_list, use_layering = FALSE, use_layer_preference = FALSE),  # Remove spatial layering
  run_fim_reduced_experiment("FIM_No_sanction", lhs_optimize_parameter, species_list, host_sanction = 1),                                # Disable host sanctioning
  run_fim_reduced_experiment("FIM_No_temporal", lhs_optimize_parameter, species_list, use_egg_success_by_phase = FALSE)                # Remove phenological effects
  )

# Combine the outputs from each model variant into a single data frame
metrics_reduced_models <- bind_rows(fim_reduced_results)

# observed data matrics
bootstrap_metrics_observed <- resample_observed_metrics_df_results %>%
  select(all_of(metrics_keep))

bootstrap_metrics_observed$Source <- "Observed"

colnames(bootstrap_metrics_observed) <- c("Richness", "Pielou's Evenness", "Bray-Curtis Dissimilarity",
                                             "Nestedness", "Connectance", "Modularity","Source")

# delta coverage & delta loss


# ---------- 0) clean data  ----------
metric_cols <- c("Richness","Pielou's Evenness","Bray-Curtis Dissimilarity",
                 "Nestedness","Connectance","Modularity")

#  KO + FIM
sim_boot <- metrics_reduced_models %>%
  select(all_of(metric_cols), Source)

# 
obs_boot <- bootstrap_metrics_observed %>%
  select(all_of(metric_cols), Source)

# 
unique(sim_boot$Source)
# "FIM", "FIM_No_space", "FIM_No_sanction", "FIM_No_temporal"

# ---------- 1) ΔLoss （FIM）----------
# Loss(M) = sum_i ( mean_sim_i(M) - mean_obs_i )^2

obs_means <- obs_boot %>%
  select(all_of(metric_cols)) %>%
  pivot_longer(everything(), names_to="metric", values_to="val") %>%
  group_by(metric) %>% summarise(obs_mu = mean(val), .groups="drop")

sim_means <- sim_boot %>%
  select(all_of(metric_cols), Source) %>%
  pivot_longer(cols = all_of(metric_cols), names_to="metric", values_to="val") %>%
  group_by(Source, metric) %>%
  summarise(sim_mu = mean(val), .groups="drop")

loss_by_model <- sim_means %>%
  inner_join(obs_means, by="metric") %>%
  group_by(Source) %>%
  summarise(Loss = sum((sim_mu - obs_mu)^2), .groups="drop")

Loss_FIM <- loss_by_model %>% filter(Source == "FIM") %>% pull(Loss)

delta_loss <- loss_by_model %>%
  filter(Source != "FIM") %>%
  mutate(dLoss = Loss - Loss_FIM) %>%
  select(Source, dLoss)


# ---------- 2) ΔCoverage（pp， FIM）----
q_levels <- c(0.05, 0.25, 0.50, 0.75, 0.95)

boot_quantile <- function(v, q, B=1000){
  replicate(B, quantile(sample(v, replace=TRUE), q, type=8))
}

predictive_envelope_full <- function(v, q_levels=q_levels, B=1000){
  purrr::map_dfr(q_levels, function(q){
    bq <- boot_quantile(v, q, B=B)
    tibble::tibble(
      q    = q,
      mean = mean(bq),
      sd   = sd(bq),
      lo95 = quantile(bq, 0.025, type=8),
      hi95 = quantile(bq, 0.975, type=8)
    )
  })
}

envelope_for_model <- function(df_one_model, metric_cols, B=1000){
  df_one_model |>
    dplyr::select(dplyr::all_of(metric_cols)) |>
    tidyr::pivot_longer(everything(), names_to="metric", values_to="val") |>
    dplyr::group_by(metric) |>
    dplyr::summarise(env = list(predictive_envelope_full(val, q_levels, B)),
                     .groups="drop") |>
    tidyr::unnest(env)
}

# 
obs_q <- obs_boot |>
  dplyr::select(dplyr::all_of(metric_cols)) |>
  tidyr::pivot_longer(everything(), names_to="metric", values_to="val") |>
  dplyr::group_by(metric) |>
  dplyr::summarise(
    `5%`  = quantile(val, 0.05, type=8),
    `25%` = quantile(val, 0.25, type=8),
    `50%` = quantile(val, 0.50, type=8),
    `75%` = quantile(val, 0.75, type=8),
    `95%` = quantile(val, 0.95, type=8),
    .groups="drop"
  ) |>
  tidyr::pivot_longer(-metric, names_to="q_lab", values_to="obs_qval") |>
  dplyr::mutate(q = readr::parse_number(q_lab)/100)

# 
mean_exceed_norm <- function(env_df, obs_q_df){
  env_df |>
    dplyr::inner_join(obs_q_df, by=c("metric","q")) |>
    dplyr::mutate(
      half = (hi95 - lo95)/2,
      ctr  = (hi95 + lo95)/2,
      
      exceed = dplyr::case_when(
        obs_qval > hi95 ~ (obs_qval - hi95)/half,
        obs_qval < lo95 ~ (lo95 - obs_qval)/half,
        TRUE            ~ 0
      )
    ) |>
    dplyr::summarise(m = mean(exceed), .groups="drop") |>
    dplyr::pull(m)
}

# 
env_FIM   <- envelope_for_model(sim_boot |> dplyr::filter(Source=="FIM"), metric_cols, B=1000)
marg_FIM  <- mean_exceed_norm(env_FIM, obs_q)

# 
sources_KO <- setdiff(unique(sim_boot$Source), "FIM")
delta_cov_cont <- purrr::map_dfr(sources_KO, function(src){
  env_KO  <- envelope_for_model(sim_boot |> dplyr::filter(Source==src), metric_cols, B=1000)
  marg_KO <- mean_exceed_norm(env_KO, obs_q)
  tibble::tibble(Source = src, dCov_cont = marg_KO - marg_FIM)
})


######### community data


# 
sim_FIM <- simulate_figwasp_community(
  num_figs                 = 1000,
  fecundity_mean           = lhs_optimize_parameter$fecundity_mean,
  fecundity_dispersion     = lhs_optimize_parameter$fecundity_dispersion,
  entry_mu                 = lhs_optimize_parameter$entry_mu,
  entry_size               = lhs_optimize_parameter$entry_size,
  entry_priority           = lhs_optimize_parameter$entry_priority,
  species_roles            = lhs_optimize_parameter$species_roles,
  max_entry_table          = lhs_optimize_parameter$max_entry_table,
  enable_drop              = TRUE,
  drop_cancels_emergence   = FALSE,
  entry_distribution       = "lognormal",
  interaction_matrix       = lhs_optimize_parameter$interaction_matrix,
  interaction_weight       = 0,
  egg_success_prob_by_phase= lhs_optimize_parameter$egg_success_prob_by_phase,
  layer_preference         = lhs_optimize_parameter$layer_preference,
  seed                     = 42
)

sim_FIM <- sim_FIM$summary
sim_FIM <- sim_FIM[,c(18:23)]
colnames(sim_FIM) <- gsub("^emergence_", "", colnames(sim_FIM))
sim_FIM <- sim_FIM[, match(colnames(obs_comm_mat), colnames(sim_FIM))]


# KO: sim_KO_space
sim_KO_space <- simulate_figwasp_community(
  num_figs                 = 1000,
  fecundity_mean           = lhs_optimize_parameter$fecundity_mean,
  fecundity_dispersion     = lhs_optimize_parameter$fecundity_dispersion,
  entry_mu                 = lhs_optimize_parameter$entry_mu,
  entry_size               = lhs_optimize_parameter$entry_size,
  entry_priority           = lhs_optimize_parameter$entry_priority,
  species_roles            = lhs_optimize_parameter$species_roles,
  max_entry_table          = lhs_optimize_parameter$max_entry_table,
  enable_drop              = TRUE,
  drop_cancels_emergence   = FALSE,
  entry_distribution       = "lognormal",
  interaction_matrix       = lhs_optimize_parameter$interaction_matrix,
  interaction_weight       = 0,
  egg_success_prob_by_phase= lhs_optimize_parameter$egg_success_prob_by_phase,
  layer_preference         = lhs_optimize_parameter$layer_preference,
  use_layering             = FALSE,  
  use_layer_preference     = FALSE,
  seed                     = 42
)

sim_KO_space <- sim_KO_space$summary
sim_KO_space <- sim_KO_space[,c(18:23)]
colnames(sim_KO_space) <- gsub("^emergence_", "", colnames(sim_KO_space))
sim_KO_space <- sim_KO_space[, match(colnames(obs_comm_mat), colnames(sim_KO_space))]

# KO: sim_KO_sanction
sim_KO_sanction <- simulate_figwasp_community(
  num_figs                 = 1000,
  fecundity_mean           = lhs_optimize_parameter$fecundity_mean,
  fecundity_dispersion     = lhs_optimize_parameter$fecundity_dispersion,
  entry_mu                 = lhs_optimize_parameter$entry_mu,
  entry_size               = lhs_optimize_parameter$entry_size,
  entry_priority           = lhs_optimize_parameter$entry_priority,
  species_roles            = lhs_optimize_parameter$species_roles,
  max_entry_table          = lhs_optimize_parameter$max_entry_table,
  enable_drop              = TRUE,
  drop_cancels_emergence   = FALSE,
  entry_distribution       = "lognormal",
  interaction_matrix       = lhs_optimize_parameter$interaction_matrix,
  interaction_weight       = 0,
  egg_success_prob_by_phase= lhs_optimize_parameter$egg_success_prob_by_phase,
  layer_preference         = lhs_optimize_parameter$layer_preference,
  host_sanction            = 1,      
  seed                     = 42
)

sim_KO_sanction <- sim_KO_sanction$summary
sim_KO_sanction <- sim_KO_sanction[,c(18:23)]
colnames(sim_KO_sanction) <- gsub("^emergence_", "", colnames(sim_KO_sanction))
sim_KO_sanction <- sim_KO_sanction[, match(colnames(obs_comm_mat), colnames(sim_KO_sanction))]



# KO: sim_KO_temporal
sim_KO_temporal <- simulate_figwasp_community(
  num_figs                 = 1000,
  fecundity_mean           = lhs_optimize_parameter$fecundity_mean,
  fecundity_dispersion     = lhs_optimize_parameter$fecundity_dispersion,
  entry_mu                 = lhs_optimize_parameter$entry_mu,
  entry_size               = lhs_optimize_parameter$entry_size,
  entry_priority           = lhs_optimize_parameter$entry_priority,
  species_roles            = lhs_optimize_parameter$species_roles,
  max_entry_table          = lhs_optimize_parameter$max_entry_table,
  enable_drop              = TRUE,
  drop_cancels_emergence   = FALSE,
  entry_distribution       = "lognormal",
  interaction_matrix       = lhs_optimize_parameter$interaction_matrix,
  interaction_weight       = 0,
  egg_success_prob_by_phase= lhs_optimize_parameter$egg_success_prob_by_phase,
  layer_preference         = lhs_optimize_parameter$layer_preference,
  use_egg_success_by_phase = FALSE,  
  seed                     = 42
)

sim_KO_temporal <- sim_KO_temporal$summary
sim_KO_temporal <- sim_KO_temporal[,c(18:23)]
colnames(sim_KO_temporal) <- gsub("^emergence_", "", colnames(sim_KO_temporal))
sim_KO_temporal <- sim_KO_temporal[, match(colnames(obs_comm_mat), colnames(sim_KO_temporal))]


compute_rho_matched_fast <- function(obs_comm, sim_comm, reps = 200, seed = 42) {
  set.seed(seed)
  # 1) Hellinger
  Xobs <- decostand(as.matrix(obs_comm), "hellinger")
  Xsim <- decostand(as.matrix(sim_comm), "hellinger")
  
  # 2) 
  n_obs <- nrow(Xobs)
  c_obs <- colMeans(Xobs)
  m_obs <- mean(sqrt(rowSums((Xobs - c_obs)^2)))
  
  # 3) 
  rhos <- replicate(reps, {
    idx   <- sample.int(nrow(Xsim), n_obs, replace = nrow(Xsim) < n_obs)
    Ymod  <- Xsim[idx, , drop = FALSE]
    c_mod <- colMeans(Ymod)
    m_mod <- mean(sqrt(rowSums((Ymod - c_mod)^2)))
    m_obs / m_mod
  })
  
  c(mean = mean(rhos),
    lo   = quantile(rhos, 0.025),
    hi   = quantile(rhos, 0.975))
}


rho_FIM   <- compute_rho_matched_fast(obs_comm_mat, sim_FIM, reps = 200)
rho_space <- compute_rho_matched_fast(obs_comm_mat, sim_KO_space, reps = 200)
rho_sanc  <- compute_rho_matched_fast(obs_comm_mat, sim_KO_sanction, reps = 200)
rho_temp  <- compute_rho_matched_fast(obs_comm_mat, sim_KO_temporal, reps = 200)



dRho_tbl <- tibble(
  Source   = c("FIM_No_space","FIM_No_sanction","FIM_No_temporal"),
  dRho_pct = c(100*(rho_space["mean"]/rho_FIM["mean"]-1),
               100*(rho_sanc["mean"]/rho_FIM["mean"]-1),
               100*(rho_temp["mean"]/rho_FIM["mean"]-1))
)



# ---- 5) ）----
scale_safe <- function(x) if (sd(x, na.rm=TRUE)==0) rep(0, length(x)) else as.numeric(scale(x))

fig4_summary <- delta_loss |>
  dplyr::inner_join(delta_cov_cont, by="Source") |>        
  dplyr::inner_join(dRho_tbl,     by="Source") |>
  dplyr::mutate(
    label = dplyr::recode(Source,
                          FIM_No_space    = "KO: Spatial partitioning",
                          FIM_No_sanction = "KO: Host sanctions",
                          FIM_No_temporal = "KO: Phase-specific success"
    ),
    zLoss = scale_safe(dLoss),
    zCov  = scale_safe(dCov_cont),
    rank_score = abs(zLoss) + abs(zCov)
  ) |>
  dplyr::arrange(dplyr::desc(rank_score)) |>
  dplyr::mutate(label = factor(label, levels = label))




# 
# Panel A: ΔLoss
fg4_panela <- ggplot(fig4_summary, aes(x = dLoss, y = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_col(fill = "#5DADE2", width = 0.4) +
  #scale_x_continuous(labels = label_number(accuracy = 1)) +
  scale_x_continuous(labels = label_number(big.mark=",")) +
  labs(
    #title = "A. Impact on Mean Fit (ΔLoss)",
    x = "Increase in Total Loss (relative to FIM)",
    y = NULL
  ) +
  theme_classic(base_size = 9, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(margin = margin(t = 6))
  )

# Panel B: ΔCoverage
fg4_panelb <- ggplot(fig4_summary, aes(x = dCov_cont, y = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_col(fill = "#F5B041", width = 0.4) +
  geom_text(aes(label = sprintf("Δρ = %.1f%%", dRho_pct),
                x = ifelse(dCov_cont >= 0, dCov_cont + 0.02, dCov_cont - 0.02)),
            hjust = ifelse(fig4_summary$dCov_cont >= 0, 0, 1),
            vjust = 0.5, size = 2) +
  labs(
    #title = "B. Impact on Variance Adequacy (continuous gap)",
       x = "Increase in normalized out-of-band gap (relative to FIM)",
       y = NULL) +
  theme_classic(base_size = 9, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(margin = margin(t = 6)),
    plot.margin = margin(t = 5, r = 20, b = 5, l = 5) 
  )


figure_4 <- fg4_panela + fg4_panelb +  plot_layout(ncol = 2, nrow = 2, widths = c(1, 1.1)) +
            plot_annotation(tag_levels = 'A')
figure_4

ggsave("Figure 4.pdf", figure_4, width = 183, height = 100, units = "mm", dpi = 500)






