# Summarize Simulated Metrics from Ficus Wasp Community Output

This function computes summary statistics from the output of the
\`simulate_figwasp_community()\` function, including per-species entry,
egg, and emergence means, dropped fruit summary, fig-level richness, and
species presence proportions by richness class.

## Usage

``` r
summarize_simulated_metrics(all_metrics, species_list, version_label = "v")
```

## Arguments

- all_metrics:

  Data frame. The \`\$summary\` output from
  \`simulate_figwasp_community()\`.

- species_list:

  Character vector of species names (used to extract relevant columns).

- version_label:

  Optional string for labeling plots. Default is "v".

## Value

A named list with the following elements:

- stats_table:

  Data frame with Entry, Egg, and Emergence means per species

- drop_summary:

  Table summarizing dropped vs non-dropped figs

- richness_vector:

  Integer vector of species richness per fig

- species_by_richness:

  Data frame with species presence proportions by richness

- heatmap:

  The richness–species presence ggplot object (printed automatically)

## Details

A heatmap is plotted to visualize how each species' presence frequency
varies across figs with different species richness levels.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- simulate_figwasp_community(...)
summary <- summarize_simulated_metrics(result$summary, species_list)
} # }
```
