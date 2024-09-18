
<!-- README.md is generated from README.Rmd. Please edit that file -->

# biolutoxR package: an App-Shiny package <img src="img/logo.png" alt="biolutoxR logo" width="100" align="right"/>

<!-- badges: start -->
<!-- badges: end -->

The goal of biolutoxR is to facilitate the data analysis of toxicity
test based on bacterial bioluminescence inhibition .

## Installation & loading

You can install the biolutoxR package like so:

``` r
if (!require(devtools)) {
    install.packages("devtools")
}
devtools::install_github("bbellier/biolutoxR_package", force = TRUE)
```

You can loading the biolutoxR package like so:

``` r
library(biolutoxR)
```

## Example

An example with preloaded data is available by running the function
“example.biolutoxR()”:

``` r
# example.biolutoxR()
```

To enter your own data, you can use the “run.biolutoxR()” function:

``` r
# run.biolutoxR()
```

## Citation

Please cite this package as:

> In prep.

Also, to cite the use of this package, you can use the
“citation.biolutoxR()” function:

``` r
citation.biolutoxR()
#> [1] "To cite this document: in prep. Also, you can find all about this package in: https://bbellier.github.io/biolutoxR_website/."
```
