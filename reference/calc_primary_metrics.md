# Calculate primary diagnostic metrics from a community matrix

Calculates process-proximal diagnostic metrics from a fruit-level
community matrix, including total abundance, species-specific abundance,
occurrence frequency, guild-level contribution, host-parasitoid
conditional relationships, mean-variance relationships, and
extreme-value figs.

## Usage

``` r
calc_primary_metrics(
  comm_mat,
  dataset_name,
  guild_lookup = NULL,
  host_pairs = NULL,
  low_total_cutoff = NULL,
  high_total_cutoff = NULL,
  dominance_threshold = 0.8
)
```

## Arguments

- comm_mat:

  A matrix or data frame with figs as rows and wasp species as columns.

- dataset_name:

  Character. Label for the dataset, such as \`"Observed"\` or
  \`"Simulated"\`.

- guild_lookup:

  Optional data frame with columns \`species\` and \`guild\`.

- host_pairs:

  Optional data frame with columns \`parasitoid\` and \`host\`.

- low_total_cutoff:

  Optional numeric cutoff defining very low total abundance.

- high_total_cutoff:

  Optional numeric cutoff defining very high total abundance.

- dominance_threshold:

  Numeric. Relative-abundance threshold used to define single-species
  dominance.

## Value

A list of data frames containing primary diagnostic summaries.
