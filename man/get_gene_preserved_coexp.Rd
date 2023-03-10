% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preserved_coexp.R
\name{get_gene_preserved_coexp}
\alias{get_gene_preserved_coexp}
\title{Compute the preserved co-expression auroc for a single gene}
\usage{
get_gene_preserved_coexp(
  aggregated_fetal_network,
  test_network,
  test_gene,
  num_top = 10
)
}
\arguments{
\item{aggregated_fetal_network}{the aggregated fetal co-expression network to act as a reference network. Any co-expression network will work, with the same gene annotations}

\item{test_network}{the test co-expression network, we recommend using the rank standardized networks}

\item{test_gene}{string, the gene name}

\item{num_top}{integer, the number of top co-expressed partners of the test_gene to compute the auroc over}
}
\value{
an auroc statistic
}
\description{
Compute the preserved co-expression auroc for a single gene
}
