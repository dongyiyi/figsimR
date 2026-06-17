# Compute Alpha Diversity at the Fig Level

Calculates sample-wise (e.g., per-fig) alpha diversity metrics,
including species richness, Shannon diversity, Simpson diversity, and
Pielou's evenness.

## Usage

``` r
compute_alpha_diversity(data, sample_id_col = NULL)
```

## Arguments

- data:

  A data frame where rows represent individual samples (e.g., figs) and
  columns represent species abundances.

- sample_id_col:

  Optional. Name of the column containing sample (fig) IDs. If provided,
  it is preserved in the output. If `NULL`, row names are used.

## Value

A data frame with alpha diversity metrics for each sample (row),
including:

- `Sample_ID`: The sample identifier.

- `Richness`: Number of species (non-zero abundances).

- `Shannon`: Shannon diversity index.

- `Simpson`: Simpson diversity index.

- `Evenness`: Pielou's evenness, computed as `Shannon / log(Richness)`.

All values are rounded to 3 decimal places.

## Details

This function is useful for fig-level community analysis, allowing
quantification of diversity across replicates or treatments. All input
species columns must be numeric. Non-numeric columns (except
`sample_id_col`) are ignored.

## Examples

``` r
if (FALSE) { # \dontrun{
data <- data.frame(
  FigID = paste0("F", 1:3),
  SpeciesA = c(3, 0, 2),
  SpeciesB = c(1, 4, 0),
  SpeciesC = c(0, 2, 1)
)
compute_alpha_diversity(data, sample_id_col = "FigID")
} # }
```
