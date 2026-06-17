# Perform Kruskal-Wallis and Dunn's Post Hoc Test for Multiple Metrics

This function applies the Kruskal-Wallis test to each unique metric in
the input data and conducts a post hoc Dunn's test with
Benjamini-Hochberg adjustment if the Kruskal-Wallis test is significant
(p \< 0.05).

## Usage

``` r
run_kruskal_dunn(
  df,
  metric_col = "Metric",
  group_col = "Source",
  value_col = "Value"
)
```

## Arguments

- df:

  A long-format data frame containing one row per metric observation.

- metric_col:

  Name of the column that identifies the metric (default: "Metric").

- group_col:

  Name of the column that contains grouping information (e.g., model
  type).

- value_col:

  Name of the column containing the numeric values to compare.

## Value

A data frame with one row per comparison, including:

- `Metric`:

  The metric name.

- `KW_p`:

  The Kruskal-Wallis test p-value.

- `Comparison`:

  Group comparison name (e.g., "A - B").

- `Z`:

  The z-score from Dunn's test.

- `P.unadj`:

  Unadjusted p-value from Dunn's test.

- `P.adj`:

  BH-adjusted p-value.

If the Kruskal-Wallis test is not significant for a metric, only the KW
p-value will be shown, and other fields will be `NA`.

## Examples

``` r
if (FALSE) { # \dontrun{
long_df <- data.frame(
  Metric = rep(c("Richness", "Shannon"), each = 30),
  Source = rep(c("Observed", "Simulated1", "Simulated2"), 20),
  Value = rnorm(60)
)
result <- run_kruskal_dunn(long_df)
print(result)
} # }
```
