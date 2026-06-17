# Estimate observed-to-simulated dispersion ratio by matched resampling

Estimates the ratio of observed to simulated multivariate dispersion by
repeatedly sampling simulated figs to match the number of observed figs.

## Usage

``` r
compute_rho_matched_fast(obs_comm, sim_comm, reps = 200, seed = 42)
```

## Arguments

- obs_comm:

  Observed community matrix with figs as rows and species as columns.

- sim_comm:

  Simulated community matrix with figs as rows and species as columns.

- reps:

  Integer. Number of matched resampling replicates.

- seed:

  Integer. Random seed.

## Value

A named numeric vector with the mean ratio and 95 bounds: `mean`, `lo`,
and `hi`.
