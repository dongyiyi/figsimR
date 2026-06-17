# Run a Reduced FIM Simulation and Bootstrap Key Metrics

This function runs a reduced version of the Full Intrinsic Model (FIM)
simulation using configurable ecological mechanisms, and computes
bootstrapped summaries of key diversity and network metrics.

## Usage

``` r
run_fim_reduced_experiment(
  version_label,
  lhs_param_list,
  species_list,
  num_figs = 1000,
  host_sanction = 0.8,
  use_layering = TRUE,
  use_layer_preference = TRUE,
  use_supplemental_parasitism = FALSE,
  enable_drop = TRUE,
  use_egg_success_by_phase = TRUE,
  drop_cancels_emergence = FALSE,
  use_sink_strength = FALSE,
  sink_min_prop = 0.2,
  sink_max_prop = 0.95,
  sink_controls_drop = FALSE,
  sink_drop_condition = c("outside_range", "below_min", "above_max"),
  seed = 42,
  resample_n = 500,
  resample_size = 200
)
```

## Arguments

- version_label:

  Character. Label to assign to the output (e.g., "NoLayering").

- lhs_param_list:

  List of input parameters (entry, fecundity, etc.).

- species_list:

  List of species to include in network summary.

- num_figs:

  Integer. Number of figs to simulate. Default is 1000.

- host_sanction:

  Numeric. Host sanction strength (0–1).

- use_layering:

  Logical. Whether to enable spatial layering. Default TRUE.

- use_layer_preference:

  Logical. Whether to apply layer preference. Default TRUE.

- use_supplemental_parasitism:

  Logical. Whether to activate supplemental parasitism. Default FALSE.

- enable_drop:

  Logical. Enable legacy drop (host sanction) mechanism. If flower use
  proportion exceeds host_sanction, fig may abort. In the output
  "drop_prob" column: Ranges from 0 (never drop) to 1 (always drop);
  "is_dropped" column: no drop (0) and drop (1). When sink strength is
  used (use_sink_strength = TRUE) and enable_drop = FALSE (host sanction
  is disabled); the columns "is_dropped" and "drop_prob" in the outputs
  are not meaningful.

- use_egg_success_by_phase:

  Logical. Use egg success probabilities by phase. Default TRUE.

- drop_cancels_emergence:

  Logical. Whether fruit drop cancels emergence. Default FALSE.

- use_sink_strength:

  Logical. If TRUE, compute sink-strength metrics but do NOT alter
  legacy drop decision: use_sink_strength = TRUE and enable_drop = FALSE
  and the columns "is_dropped" and "drop_prob" are not meaningful. The
  fig is dropped when the value of "sink_prop" falls outside the
  thresholds defined by sink_min_prop and sink_max_prop; The columns
  "sink_below_min" and "sink_above_max" are logical: TRUE = drop; FALSE
  = no drop. Adds columns: sink_strength, sink_prop, sink_below_min,
  sink_above_max. Sink strength threshold calculation: sink_strength =
  sink_linear_coef \* (sink_w_gall \* galled_total + sink_w_seed \*
  seed_count).

- sink_min_prop:

  Numeric in \[0,1\]. Reference lower bound of sink proportion (default
  0.20).

- sink_max_prop:

  Numeric in \[0,1\]. Reference upper bound of sink proportion (default
  0.95).

- sink_controls_drop:

  Logical. If `TRUE`, sink-strength thresholds are used to mark figs as
  dropped after simulation.

- sink_drop_condition:

  Character. Sink-based drop rule; one of `"outside_range"`,
  `"below_min"`, or `"above_max"`.

- seed:

  Integer. Random seed for reproducibility.

- resample_n:

  Integer. Number of bootstrap iterations. Default 500.

- resample_size:

  Integer. Sample size per bootstrap iteration. Default 200.

## Value

A data frame of 6 metrics (plus source label), including:

- Richness

- Pielou's Evenness

- Bray-Curtis Dissimilarity

- Nestedness

- Connectance

- Modularity

## See also

[`simulate_figwasp_community`](https://dongyiyi.github.io/figsimR/reference/simulate_figwasp_community.md),
[`resample_metrics_from_simulation`](https://dongyiyi.github.io/figsimR/reference/resample_metrics_from_simulation.md),
[`summarize_simulated_metrics`](https://dongyiyi.github.io/figsimR/reference/summarize_simulated_metrics.md)
