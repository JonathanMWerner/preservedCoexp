

#'Compute an auroc statistic using the Mann Whitney U test
#'@param scores some metric for ranking data, must be compatible with the rank() function
#'@param label a binary numeric vector indicating data belonging to positives (1's) or negatives (0's)
#'@return an auroc statistic
#'
#'@export

compute_auc = function(scores, label){
  label = as.logical(label)
  n1 = as.numeric(sum(label))
  n2 = as.numeric(sum(!label))
  R1 = sum(rank(scores)[label])
  U1 = R1 - n1 * (n1 + 1)/2
  auc = U1/(n1 * n2)
  return(auc)
}


#'Compute the preserved co-expression auroc for a single gene
#'@param test_network the test co-expression network, we recommend using the rank standardized networks
#'@param test_gene string, the gene name
#'@param num_top integer, the number of top co-expressed partners of the test_gene to compute the auroc over
#'@return an auroc statistic
#'
#'@export

get_gene_preserved_coexp = function(test_network, test_gene, num_top = 10){
  
  #Get the top 10 coexpressed genes
  test_gene_coexp = aggregated_fetal_network[test_gene, ]
  test_gene_coexp = test_gene_coexp[names(test_gene_coexp) != test_gene]
  top_coexp = names(test_gene_coexp[order(test_gene_coexp, decreasing = T)][1:num_top])
  
  #Calculate the AUC using the organoid data
  scores = test_network[test_gene, ]
  scores = scores[names(scores) != test_gene]
  labels = names(scores) %in% top_coexp
  
  #Get the auc
  auc = compute_auc(scores, labels)
  return(auc)
}

#'Compute the preserved co-expression scores for genes in a geneset
#'@param test_network the test co-expression network, we recommend using the rank standardized networks
#'@param geneset vector of strings
#'@param num_top integer, the number of top co-expressed partners of the test_gene to compute the auroc over
#'@param parallel boolean, if the parallelized version should be implemented, requires the parallel R package
#'@param cores integer, the number of computer cores for the parallelized version
#'@return a named vector of gene preserved co-expression aurocs
#'
#'@export

get_geneset_preserved_coexp = function(test_network, geneset, num_top = 10, parallel = F, cores = 10){
  
  if(parallel){
    scores = unlist(parallel::mclapply(1:length(geneset), function(i) get_gene_preserved_coexp(test_network, geneset[i], num_top), mc.cores = cores ))
    names(scores) = geneset
    return(scores) 
    
  }else{
    scores = sapply(1:length(geneset), function(i) get_gene_preserved_coexp(test_network, geneset[i], num_top) )
    names(scores) = geneset
    return(scores)  
  }
}

#'Compute the preserved co-expression score for a geneset
#'@param test_network the test co-expression network, we recommend using the rank standardized networks
#'@param geneset vector of strings
#'@param num_top integer, the number of top co-expressed partners of the test_gene to compute the auroc over
#'@param parallel boolean, if the parallelized version should be implemented, requires the parallel R package
#'@param cores integer, the number of computer cores for the parallelized version
#'@return integer, the mean gene preserved co-expression auroc for the geneset
#'
#'@export

get_geneset_preserved_score = function(test_network, geneset, num_top = 10, parallel = F, cores = 10){
  
  if(parallel){
    return(mean(unlist(parallel::mclapply(1:length(geneset), function(i) get_gene_preserved_coexp(test_network, geneset[i], num_top), mc.cores = cores ))))
  }else{
  return(mean(sapply(1:length(geneset), function(i) get_gene_preserved_coexp(test_network, geneset[i], num_top) )))
  }
  
}

#'Plot the preserved co-expression scores for the top 100 Fetal MetaMarkers, placing the user's dataset in reference to our meta-analysis of fetal datasets
#'@param user_network gene x gene co-expression network of the user's data, we recommend using the rank standardized network
#'@return ggplot object
#'
#'@export

plot_metafetal_results = function(user_network){
  
  
}

#'Plot the preserved co-expression scores for the top 100 Fetal MetaMarkers, placing the user's dataset in reference to our meta-analysis of organoid datasets
#'@param user_network gene x gene co-expression network of the user's data, we recommend using the rank standardized network
#'@return ggplot object
#'
#'@export
plot_metaorganoid_results = function(user_network){
  
  
}







