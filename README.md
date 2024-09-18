---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# microtoxR package: an App-Shiny package <img src="man/figures/logo.png" alt="microtoxR logo" width="100" align="right"/>

<!-- badges: start -->
<!-- badges: end -->

The goal of microtoxR is to facilitate the data analysis of the Microtox® acute toxicity test.

## Installation & loading

You can install the microtoxR package like so: 

```{r install, echo=TRUE, results='hide', message=FALSE, warning=FALSE}
# install.packages("pak")
pak::pkg_install("bbellier/microtoxR_package")
```

You can loading the microtoxR package like so: 

```{r loading}
library(microtoxR)
```


## Example

An example with preloaded data is available by running the function "example.microtoxR()":
```{r example}
# example.microtoxR()
```

To enter your own data, you can use the "run.microtoxR()" function:
```{r run}
# run.microtoxR()
```

## Citation

Please cite this package as:

> Le Picard, C. and Bellier, B., (2024). microtoxR: An R-Shiny package for easy performing data analysis of microtox® test. SoftwareX, XX, XX-XX. https://doi.org/XXXX

Also, to cite the use of this package, you can use the "citation.microtoxR()" function:
```{r citation}
citation.microtoxR()
```

