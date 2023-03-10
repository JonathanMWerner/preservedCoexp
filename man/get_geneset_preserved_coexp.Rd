% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preserved_coexp.R
\name{get_geneset_preserved_coexp}
\alias{get_geneset_preserved_coexp}
\title{Compute the preserved co-expression scores for genes in a geneset}
\usage{
get_geneset_preserved_coexp(
  aggregated_fetal_network,
  test_network,
  geneset,
  num_top = 10,
  parallel = F,
  cores = 10
)
}
\arguments{
\item{aggregated_fetal_network}{the aggregated fetal co-expression network to act as a reference network. Any co-expression network will work, with the same gene annotations}

\item{test_network}{the test co-expression network, we recommend using the rank standardized networks}

\item{geneset}{vector of strings}

\item{num_top}{integer, the number of top co-expressed partners of the test_gene to compute the auroc over}

\item{parallel}{boolean, if the parallelized version should be implemented, requires the parallel R package}

\item{cores}{integer, the number of computer cores for the parallelized version}
}
\value{
a named vector of gene preserved co-expression aurocs
}
\description{
Compute the preserved co-expression scores for genes in a geneset
}
