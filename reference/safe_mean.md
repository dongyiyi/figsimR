# Calculate a mean with empty-vector protection

Returns the mean of a numeric vector while ignoring missing values. If
the input vector has length zero, the function returns `NA_real_`.

## Usage

``` r
safe_mean(x)
```

## Arguments

- x:

  A numeric vector.

## Value

A numeric value.
