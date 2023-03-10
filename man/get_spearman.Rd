% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coexpression.R
\name{get_spearman}
\alias{get_spearman}
\title{Generate a spearman correlation co-expression matrix}
\usage{
get_spearman(exp_data)
}
\arguments{
\item{exp_data}{a genes x cells expression matrix, as output from fit_to_GO()}
}
\value{
a gene x gene spearman correlation matrix
}
\description{
Generate a spearman correlation co-expression matrix
}
