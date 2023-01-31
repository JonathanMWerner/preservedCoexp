
#' Fit an expression matrix with genes as rows to the GO universe of GO annotations
#' @param dataset a matrix with genes as rows, samples (cells) as columns, and expression counts as data, we recommend a CPM normalized expression matrix
#' @return An expression matrix zero-padded for the missing GO genes in the original expression matrix
#' 
#' @export

fit_to_GO = function(dataset){

  dataset_genes = rownames(dataset)
  dataset_genes_in_go = dataset_genes[dataset_genes %in% go_genes]
  dataset_genes_missing_from_go = go_genes[! go_genes %in% dataset_genes]
  
  dataset_gene_index = dataset_genes %in% dataset_genes_in_go
  exp_data = dataset[dataset_gene_index, ]
  
  go_missing = matrix(0, nrow = length(dataset_genes_missing_from_go), ncol = ncol(exp_data), dimnames = list(dataset_genes_missing_from_go, colnames(exp_data)))
  exp_data = rbind(exp_data, go_missing)
  return(exp_data)
}
  









