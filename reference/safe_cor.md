# Calculate a correlation with safety checks

Returns the Pearson correlation between two numeric vectors. If either
vector has fewer than two values or zero variance, the function returns
`NA_real_`.

## Usage

``` r
safe_cor(x, y)
```

## Arguments

- x:

  A numeric vector.

- y:

  A numeric vector.

## Value

A numeric value.
