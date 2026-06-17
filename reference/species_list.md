# List of All Fig Wasp Species Observed on *Ficus racemosa*

This character vector contains the standardized names of all fig wasp
species observed on *Ficus racemosa*, used in the fig wasp community
simulator. Each species name is formatted using an underscore to join
the genus and species, e.g., \`Genus_species\`. This vector defines the
species order used in all parameter lists (e.g., \`entry_mu\`,
\`fecundity_mean\`, etc.), and it serves as the default species pool for
simulations.

## Usage

``` r
data(species_list)
```

## Format

A character vector of length 6.

## See also

[`parameter_list_default`](https://dongyiyi.github.io/figsimR/reference/parameter_list_default.md)
for parameters corresponding to each species.

## Examples

``` r
data(species_list)
print(species_list)
#> [1] "Sycophaga_testacea"  "Apocrypta_sp"        "Sycophaga_mayri"    
#> [4] "Ceratosolen_sp"      "Sycophaga_agraensis" "Apocrypta_westwoodi"
# [1] "Sycophaga_testacea"   "Apocrypta_sp"
#     "Sycophaga_mayri"     "Ceratosolen_sp"
#     "Sycophaga_agraensis" "Apocrypta_westwoodi"
```
