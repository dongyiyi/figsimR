# Scale a numeric vector with zero-variance protection

Standardizes a numeric vector. If the vector has zero variance, the
function returns a vector of zeros with the same length.

## Usage

``` r
scale_safe(x)
```

## Arguments

- x:

  A numeric vector.

## Value

A numeric vector.
