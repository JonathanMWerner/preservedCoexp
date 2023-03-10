% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preserved_coexp.R
\name{plot_meta_results}
\alias{plot_meta_results}
\title{Compute the preserved co-expression scores for the top 100 MetaMarker genes per cell-type for the test network and plot them in reference to our metaanalysis of fetal and organoid datasets}
\usage{
plot_meta_results(
  aggregated_fetal_network,
  test_network,
  meta_markers,
  meta_presCoexp_df,
  parallel = F
)
}
\arguments{
\item{aggregated_fetal_network}{the aggregated fetal co-expression network to act as a reference network. Any co-expression network will work, with the same gene annotations}

\item{test_network}{the test co-expression network, we recommend using the rank standardized networks}

\item{meta_presCoexp_df}{data frame of MetaMarker scores for the fetal and organoid datasets we assessed in our meta-analysis, included as data for this package}

\item{parallel}{boolean, if the parallelized version should be implemented, requires the parallel R package}

\item{metamarkers}{A datatable of metamarkers, included as data for this package}
}
\value{
a list including the MetaMarker scores for the test network and two ggplot objects plotting the results for the fetal and organoid datasets separately
}
\description{
Compute the preserved co-expression scores for the top 100 MetaMarker genes per cell-type for the test network and plot them in reference to our metaanalysis of fetal and organoid datasets
}
