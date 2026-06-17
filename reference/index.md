# Package index

## Core simulation

- [`simulate_figwasp_community()`](https://dongyiyi.github.io/figsimR/reference/simulate_figwasp_community.md)
  :

  Simulate Fig Wasp Community Assembly (Observed on *Ficus racemosa*)

## Built-in data

- [`observed_data`](https://dongyiyi.github.io/figsimR/reference/observed_data.md)
  : Fig wasp community data from Aung et al. 2022

- [`parameter_list_default`](https://dongyiyi.github.io/figsimR/reference/parameter_list_default.md)
  : Default Parameter List for Simulating Fig Wasp Communities

- [`species_list`](https://dongyiyi.github.io/figsimR/reference/species_list.md)
  :

  List of All Fig Wasp Species Observed on *Ficus racemosa*

## Metrics and diagnostics

- [`calc_all_metrics()`](https://dongyiyi.github.io/figsimR/reference/calc_all_metrics.md)
  : Calculate Alpha, Beta, and Network Metrics from Community Data
- [`calc_primary_metrics()`](https://dongyiyi.github.io/figsimR/reference/calc_primary_metrics.md)
  : Calculate primary diagnostic metrics from a community matrix
- [`calc_secondary_metrics()`](https://dongyiyi.github.io/figsimR/reference/calc_secondary_metrics.md)
  : Calculate emergent secondary community metrics
- [`calculate_metric_loss()`](https://dongyiyi.github.io/figsimR/reference/calculate_metric_loss.md)
  : Calculate Total Loss Between Simulated and Observed Metrics
- [`make_secondary_metric_interval_df()`](https://dongyiyi.github.io/figsimR/reference/make_secondary_metric_interval_df.md)
  : Create interval summaries for secondary metrics
- [`plot_distribution_metrics_grid_pvalue()`](https://dongyiyi.github.io/figsimR/reference/plot_distribution_metrics_grid_pvalue.md)
  : Plot Metric Distributions with Wilcoxon Test Significance
- [`plot_distribution_metrics_grid_pvalue_legend()`](https://dongyiyi.github.io/figsimR/reference/plot_distribution_metrics_grid_pvalue_legend.md)
  : Plot Metric Distributions with KS Test and Shared Legend
- [`resample_metrics_from_simulation()`](https://dongyiyi.github.io/figsimR/reference/resample_metrics_from_simulation.md)
  : Resample and Compute Community Metrics from Simulated Data
- [`resample_observed_all_metrics()`](https://dongyiyi.github.io/figsimR/reference/resample_observed_all_metrics.md)
  : Bootstrap Resampling of Observed Community Metrics
- [`resample_simulated_metrics()`](https://dongyiyi.github.io/figsimR/reference/resample_simulated_metrics.md)
  : Bootstrap Resampling of Simulated Community Metrics
- [`summarize_simulated_metrics()`](https://dongyiyi.github.io/figsimR/reference/summarize_simulated_metrics.md)
  : Summarize Simulated Metrics from Ficus Wasp Community Output

## Primary and emergent metrics

- [`calc_richness_per_fig()`](https://dongyiyi.github.io/figsimR/reference/calc_richness_per_fig.md)
  : Calculate species richness per fig
- [`calc_per_fig_richness_evenness()`](https://dongyiyi.github.io/figsimR/reference/calc_per_fig_richness_evenness.md)
  : Calculate fig-level richness and Pielou's evenness
- [`calc_bray_distribution()`](https://dongyiyi.github.io/figsimR/reference/calc_bray_distribution.md)
  : Calculate pairwise Bray-Curtis dissimilarity distribution
- [`calc_dispersion_distribution()`](https://dongyiyi.github.io/figsimR/reference/calc_dispersion_distribution.md)
  : Calculate distance-to-centroid distribution
- [`calc_primary_metrics()`](https://dongyiyi.github.io/figsimR/reference/calc_primary_metrics.md)
  : Calculate primary diagnostic metrics from a community matrix
- [`calc_secondary_metrics()`](https://dongyiyi.github.io/figsimR/reference/calc_secondary_metrics.md)
  : Calculate emergent secondary community metrics
- [`collect_simulated_secondary_distributions()`](https://dongyiyi.github.io/figsimR/reference/collect_simulated_secondary_distributions.md)
  : Collect simulated distributions of secondary metrics
- [`compute_alpha_diversity()`](https://dongyiyi.github.io/figsimR/reference/compute_alpha_diversity.md)
  : Compute Alpha Diversity at the Fig Level
- [`compute_rho_matched_fast()`](https://dongyiyi.github.io/figsimR/reference/compute_rho_matched_fast.md)
  : Estimate observed-to-simulated dispersion ratio by matched
  resampling
- [`dispersion_test()`](https://dongyiyi.github.io/figsimR/reference/dispersion_test.md)
  : Compare multivariate dispersion between observed and simulated
  communities

## Parameter optimization and resampling

- [`convert_parameter_list_to_param_ranges()`](https://dongyiyi.github.io/figsimR/reference/convert_parameter_list_to_param_ranges.md)
  : Convert a Parameter List into a Param Range List for Sampling
- [`rebuild_parameter_list_from_row()`](https://dongyiyi.github.io/figsimR/reference/rebuild_parameter_list_from_row.md)
  : Rebuild Parameter List from a Row of Flat Parameter Values
- [`explore_parameter_space_lhs()`](https://dongyiyi.github.io/figsimR/reference/explore_parameter_space_lhs.md)
  : Explore Parameter Space Using Latin Hypercube Sampling
- [`resample_simulated_metrics()`](https://dongyiyi.github.io/figsimR/reference/resample_simulated_metrics.md)
  : Bootstrap Resampling of Simulated Community Metrics
- [`resample_observed_all_metrics()`](https://dongyiyi.github.io/figsimR/reference/resample_observed_all_metrics.md)
  : Bootstrap Resampling of Observed Community Metrics

## Other functions

- [`boot_quantile()`](https://dongyiyi.github.io/figsimR/reference/boot_quantile.md)
  : Bootstrap a quantile
- [`envelope_for_model()`](https://dongyiyi.github.io/figsimR/reference/envelope_for_model.md)
  : Build predictive envelopes for model metrics
- [`mean_exceed_norm()`](https://dongyiyi.github.io/figsimR/reference/mean_exceed_norm.md)
  : Calculate mean normalized exceedance outside a predictive envelope
- [`predictive_envelope_full()`](https://dongyiyi.github.io/figsimR/reference/predictive_envelope_full.md)
  : Build a bootstrap predictive envelope
- [`run_all_pairwise_tests_tidyverse()`](https://dongyiyi.github.io/figsimR/reference/run_all_pairwise_tests_tidyverse.md)
  : Perform Pairwise Dunn's Tests for Multiple Metrics (Tidyverse
  Version)
- [`run_fim_reduced_experiment()`](https://dongyiyi.github.io/figsimR/reference/run_fim_reduced_experiment.md)
  : Run a Reduced FIM Simulation and Bootstrap Key Metrics
- [`run_kruskal_dunn()`](https://dongyiyi.github.io/figsimR/reference/run_kruskal_dunn.md)
  : Perform Kruskal-Wallis and Dunn's Post Hoc Test for Multiple Metrics
- [`run_sim_and_bootstrap_multiple_weights()`](https://dongyiyi.github.io/figsimR/reference/run_sim_and_bootstrap_multiple_weights.md)
  : Run Simulations and Bootstrap Community Metrics Across Interaction
  Weights
- [`safe_cor()`](https://dongyiyi.github.io/figsimR/reference/safe_cor.md)
  : Calculate a correlation with safety checks
- [`safe_mean()`](https://dongyiyi.github.io/figsimR/reference/safe_mean.md)
  : Calculate a mean with empty-vector protection
- [`scale_safe()`](https://dongyiyi.github.io/figsimR/reference/scale_safe.md)
  : Scale a numeric vector with zero-variance protection
