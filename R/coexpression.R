
#'Imports the python numpy library using reticulate
np <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to scipy
  np <<- reticulate::import("numpy", delay_load = TRUE)

}


#' Generate a spearman correlation co-expression matrix
#'@param exp_data a genes x cells expression matrix, as output from fit_to_GO()
#'@return a gene x gene spearman correlation matrix 
#'
#'@export

get_spearman = function(exp_data){
  
  r_genes = rownames(exp_data)
  py_exp = np$array(exp_data) #Turn R matrix into numpy array
  reticulate::source_python(system.file("python_coexpression.py", package = "preservedCoexp", mustWork = T))
  python_coexp = get_coexpression(py_exp, r_genes)
  corr_matrix = as.matrix(python_coexp) #Turn to matrix
  #Convert Nans to 0, these are genes with 0 expression in any cells
  corr_matrix[!is.finite(corr_matrix)] = 0
  
  return(corr_matrix)
}


#'Generate a rank-standardized co-expression matrix
#'@param corr_matrix a gene x gene correlation matrix, as output from get_spearman()
#'@return a gene x gene correlation matrix rank standardized between 0 and 1
#'
#'@export


rank_coexpression = function(corr_matrix){
  
  rank_corr = data.table::frank(c(corr_matrix), ties.method = 'min') - 1
  #Divide by max rank to normalize between 0 and 1
  norm_rank_corr = rank_corr / max(rank_corr)
  #Put back into a matrix
  rank_corr_matrix = matrix(norm_rank_corr, nrow = nrow(corr_matrix), ncol = ncol(corr_matrix), dimnames = list(rownames(corr_matrix), colnames(corr_matrix)))
  #Set diagonal to 1
  diag(rank_corr_matrix) = 1
  return(rank_corr_matrix)
}










