# Calculate Alpha, Beta, and Network Metrics from Community Data

This function computes key community-level summary metrics from a fig
wasp emergence matrix. It includes alpha diversity indices (richness,
Shannon, Simpson, evenness), beta diversity (Bray-Curtis, Jaccard), and
bipartite network structure metrics (nestedness, connectance, linkage
density, modularity).

## Usage

``` r
calc_all_metrics(df)
```

## Arguments

- df:

  A data.frame where each row represents a fig (or sampling unit) and
  each column corresponds to a wasp species. Values are species
  abundances (counts). Column names may include or omit the prefix
  "emergence\_".

## Value

A single-row data.frame containing the following metrics:

- `mean_richness`: Mean number of species per fig

- `mean_shannon`: Mean Shannon diversity

- `mean_simpson`: Mean Simpson diversity

- `mean_evenness`: Mean Shannon evenness (Shannon/log(richness))

- `mean_bray_curtis`: Mean Bray-Curtis dissimilarity (beta diversity)

- `mean_jaccard`: Mean Jaccard dissimilarity (beta diversity,
  presence-absence)

- `nestedness`: Nestedness based on NODF metric (from bipartite)

- `connectance`: Proportion of realized links in the species-fig network

- `links_per_species`: Average number of links per species (linkage
  density)

- `modularity`: Modularity score from community detection (igraph)

## Details

Input is typically a matrix of species emergence abundances across
multiple figs, with optional "emergence\_" prefixes in column names.

If the input matrix contains zero-only rows, they will be excluded from
beta diversity and network structure calculations. Alpha diversity is
computed using vegan, beta diversity metrics are derived from
Bray-Curtis and Jaccard distances, and network metrics are calculated
using bipartite and igraph.

Column names with the "emergence\_" prefix will be automatically
stripped. If fewer than two rows remain after removing zero-only rows,
beta and network metrics return NA.

## Examples

``` r
# Simulated emergence data
set.seed(1)
df <- data.frame(
  emergence_A = sample(0:10, 100, replace = TRUE),
  emergence_B = sample(0:5, 100, replace = TRUE),
  emergence_C = sample(0:3, 100, replace = TRUE)
)
metrics <- calc_all_metrics(df)
print(metrics)
#>   mean_richness mean_shannon mean_simpson mean_evenness mean_bray_curtis
#> 1          2.43    0.7033545    0.4440592     0.8378314        0.3962443
#>   mean_jaccard nestedness connectance links_per_species modularity
#> 1    0.2963994    27.0643   0.8181818          42.27984 0.07406561
```
