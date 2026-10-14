
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Introduction

The aim of this package is to represent phylogenetic trees in
unconstrained Euclidean space.

<!--## Dependencies
&#10;This package depends on [Stan](https://mc-stan.org/) through CmdStan. Although [other options might exist](https://mc-stan.org/docs/cmdstan-guide/cmdstan-installation.html), for instance using Conda, the simplest option is to install CmdStan via the R package [CmdStanR](https://mc-stan.org/cmdstanr/). You can do so with the commands
&#10;``` r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()
```
The last command will download and compile all the required dependencies.
-->

## Installation

You can install ContPhylo from github with:

``` r
devtools::install_github("xavierdidelot/ContPhylo")
```

The package can then be loaded using:

``` r
library(ContPhylo,quietly=T)
```

## Example

See the vignette
[here](https://github.com/xavierdidelot/ContPhylo/tree/master/vignettes).

## More information and getting help

See the help included in the package using the R command
`help(package='ContPhylo')`.

If you have any problem or question please create an issue
[here](https://github.com/xavierdidelot/ContPhylo/issues) or get in
touch by emailing `xavier.didelot@gmail.com`
