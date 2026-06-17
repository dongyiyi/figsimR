# Calculate mean normalized exceedance outside a predictive envelope

Computes the mean normalized exceedance of observed quantiles relative
to simulated predictive-envelope bounds.

## Usage

``` r
mean_exceed_norm(env_df, obs_q_df)
```

## Arguments

- env_df:

  A data frame containing predictive-envelope columns `metric`, `q`,
  `lo95`, and `hi95`.

- obs_q_df:

  A data frame containing observed quantiles with columns `metric`, `q`,
  and `obs_qval`.

## Value

A numeric value giving the mean normalized exceedance. Values of zero
indicate that observed quantiles fall within the envelope.
