# Perform Pairwise Dunn's Tests for Multiple Metrics (Tidyverse Version)

This function applies Dunn's test (with Benjamini-Hochberg adjustment)
to all pairwise comparisons of groups for each unique metric. It uses a
tidyverse pipeline to nest and map over metrics.

## Usage

``` r
run_all_pairwise_tests_tidyverse(
  long_df,
  value_col = "Value",
  group_col = "Source",
  metric_col = "Metric"
)
```

## Arguments

- long_df:

  A long-format data frame containing one row per observation.

- value_col:

  Column name for the numeric values to compare. Default is `"Value"`.

- group_col:

  Column name representing the groupings (e.g., "Source" or
  "Treatment"). Default is `"Source"`.

- metric_col:

  Column name indicating the metric type. Default is `"Metric"`.

## Value

A tidy data frame with the following columns:

- `Metric`:

  Metric name.

- `Comparison`:

  Pairwise group comparison (e.g., "Simulated - Observed").

- `Z`:

  Z statistic from Dunn's test.

- `P.unadj`:

  Unadjusted p-value.

- `P.adj`:

  Benjamini-Hochberg adjusted p-value.

## Details

This function is a shortcut for performing multiple group comparisons
across several ecological or network metrics, assuming Kruskal-Wallis
significance. Dunn’s test is performed using the
[`FSA::dunnTest()`](https://fishr-core-team.github.io/FSA/reference/dunnTest.html)
function.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- data.frame(
  Metric = rep(c("Richness", "Shannon"), each = 30),
  Source = rep(c("Observed", "Sim1", "Sim2"), times = 20),
  Value = rnorm(60)
)
result <- run_all_pairwise_tests_tidyverse(df)
print(result)
} # }
```
