
<!-- README.md is generated from README.Rmd. Please edit that file -->

# preservedCoexp

<!-- badges: start -->

<!-- badges: end -->

preservedCoexp is an R package that quantifies the strength of preserved
co-expression at the gene-level across two co-expression networks. It
was originally designed to compute preserved co-expression between human
fetal brain and human neural organoid co-expression networks, as
implemented in Werner and Gillis 2023 (include citation). The package
provides functions to compute co-expression networks from a
gene-by-sample expression matrix. The package also provides
meta-analytic fetal cell-type markers and an aggregated fetal
co-expression network as data for users to compute preserved
co-expression within their own organoid co-expression networks.

## Installation

You can install the development version of preservedCoexp from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("JonathanMWerner/preservedCoexp")
```

## Beginning analysis

preservedCoexp takes as input a normalized gene-by-sample expression
matrix. We recommend CPM normalization. This package was originally
designed to work with single-cell expression data, but a bulk
gene-by-sample expression matrix would work as well.

``` r
library(preservedCoexp)
library(dplyr)
library(ggplot2)

data('go_genes', package = 'preservedCoexp')             #List of GO gene annotations
data('fetal_meta_markers', package = 'preservedCoexp')   #Dataframe of fetal MetaMarkers 
aggregated_fetal_network = load_fetal_coexp()            #Aggregate fetal co-expression network

#Functions to compute a rank standardized co-expression network, where exp_data is a normalized gene-by-sample expression matrix
exp_data = fit_to_GO(exp_data)              #Fit to GO gene annotations
rank_mat = get_spearman(exp_data)           #Get co-expression matrix
rank_mat = rank_coexpression(rank_mat)     #Get rank standardized co-expression matrix
```

``` r
data('meta_presCoexp_df', package = 'preservedCoexp')
plot_results = plot_meta_results(aggregated_fetal_network, rank_mat, fetal_meta_markers, meta_presCoexp_df, parallel = T)

plot_results[[1]]
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
