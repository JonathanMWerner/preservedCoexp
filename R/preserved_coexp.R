
#'Load the aggregated fetal co-expression network as computed in Werner and Gillis (citation)
#'@param data_url url for the aggregated fetal co-expression network hosted on the Gillis lab public FTP server
#'@return matrix
#'@export

load_fetal_coexp = function(data_url = 'this_is_a_test_string'){
  download.file(url = data_url, destfile = 'temp_agg_fetal_network.Rdata')
  file = load('temp_agg_fetal_network.Rdata')
  aggregated_fetal_network = get(file)
  suppressWarnings(rm(file, fetal_agg_rank_mat))
  unlink('temp_agg_fetal_network.Rdata')
  gc()
  return(aggregated_fetal_network)
}


#'Load the adult co-expression network as computed in Werner and Gillis (citation)
#'@param data_url url for the adult co-expression network hosted on the Gillis lab public FTP server
#'@return matrix
#'@export

load_adult_coexp = function(data_url = 'https://labshare.cshl.edu/shares/gillislab/resource/preserved_fetal_coexp/adult_coexp_network_8_17_23.Rdata'){
  download.file(url = data_url, destfile = 'temp_adult_network.Rdata')
  file = load('temp_adult_network.Rdata')
  adult_network = get(file)
  rm(file, rank_corr_matrix) 
  unlink('temp_adult_network.Rdata')
  gc()
  return(adult_network)
}



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
#'@param aggregated_fetal_network the aggregated fetal co-expression network to act as a reference network. Any co-expression network will work, with the same gene annotations
#'@param test_network the test co-expression network, we recommend using the rank standardized networks
#'@param test_gene string, the gene name
#'@param num_top integer, the number of top co-expressed partners of the test_gene to compute the auroc over
#'@return an auroc statistic
#'
#'@export

get_gene_preserved_coexp = function(aggregated_fetal_network, test_network, test_gene, num_top = 10){
  
  if(!test_gene %in% rownames(aggregated_fetal_network)){return(NA)   #Return NA for genes not in GO
    }else{
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
}



#'Compute the preserved co-expression scores for genes in a geneset
#'@param aggregated_fetal_network the aggregated fetal co-expression network to act as a reference network. Any co-expression network will work, with the same gene annotations
#'@param test_network the test co-expression network, we recommend using the rank standardized networks
#'@param geneset vector of strings
#'@param num_top integer, the number of top co-expressed partners of the test_gene to compute the auroc over
#'@param parallel boolean, if the parallelized version should be implemented, requires the parallel R package
#'@param cores integer, the number of computer cores for the parallelized version
#'@return a named vector of gene preserved co-expression aurocs
#'
#'@export

get_geneset_preserved_coexp = function(aggregated_fetal_network,test_network, geneset, num_top = 10, parallel = F, cores = 10){
  
  if(parallel){
    scores = unlist(parallel::mclapply(1:length(geneset), function(i) get_gene_preserved_coexp(aggregated_fetal_network, test_network, geneset[i], num_top), mc.cores = cores ))
    names(scores) = geneset
    return(scores) 
    
  }else{
    scores = sapply(1:length(geneset), function(i) get_gene_preserved_coexp(aggregated_fetal_network,test_network, geneset[i], num_top) )
    names(scores) = geneset
    return(scores)  
  }
}



#'Compute the preserved co-expression score for a geneset
#'@param aggregated_fetal_network the aggregated fetal co-expression network to act as a reference network. Any co-expression network will work, with the same gene annotations
#'@param test_network the test co-expression network, we recommend using the rank standardized networks
#'@param geneset vector of strings
#'@param num_top integer, the number of top co-expressed partners of the test_gene to compute the auroc over
#'@param parallel boolean, if the parallelized version should be implemented, requires the parallel R package
#'@param cores integer, the number of computer cores for the parallelized version
#'@return integer, the mean gene preserved co-expression auroc for the geneset
#'
#'@export

get_geneset_preserved_score = function(aggregated_fetal_network,test_network, geneset, num_top = 10, parallel = F, cores = 10){
  
  if(parallel){
    return(mean(unlist(parallel::mclapply(1:length(geneset), 
                                          function(i) get_gene_preserved_coexp(aggregated_fetal_network, test_network, geneset[i], num_top), mc.cores = cores )), na.rm = T))
  }else{
  return(mean(sapply(1:length(geneset), 
                     function(i) get_gene_preserved_coexp(aggregated_fetal_network, test_network, geneset[i], num_top) ), na.rm = T))
  }
  
}


#'Compute the preserved co-expression of a cell-type marker gene set from our list of fetal MetaMarkers
#'@param aggregated_fetal_network the aggregated fetal co-expression network to act as a reference network. Any co-expression network will work, with the same gene annotations
#'@param test_network the test co-expression network, we recommend using the rank standardized networks
#'@param metamarkers A datatable of metamarkers, included as data for this package
#'@param celltype string belonging to one of the celltypes included in the MetaMarkers. Must be one of c('Dividing_Progenitor','Neural_Progenitor','Intermediate_Progenitor','GABAergic','Glutamatergic','Non-neuronal')
#'@param num_markers int, number of MetaMarkers, ranked by recurrent differential expression 
#'@param parallel boolean, if the parallelized version should be implemented, requires the parallel R package
#'@param cores integer, the number of computer cores for the parallelized version
#'@return a named vector of gene preserved co-expression aurocs
#'
#'@export

get_metaMarker_preserved_coexp = function(aggregated_fetal_network, test_network, metamarkers, celltype, num_markers, parallel = F, cores = 10){
  
  if(!celltype %in% c('Dividing_Progenitor','Neural_Progenitor','Intermediate_Progenitor','GABAergic','Glutamatergic','Non-neuronal')){
    print("Please pick a cell type in c('Dividing_Progenitor','Neural_Progenitor','Intermediate_Progenitor','GABAergic','Glutamatergic','Non-neuronal') ")
    return(NA)
  }else{
    
    celltype_markers = dplyr::filter(metamarkers, cell_type == celltype & rank <= num_markers)
    
    if(parallel){
      marker_aurocs = get_geneset_preserved_coexp(aggregated_fetal_network, test_network, celltype_markers$gene, parallel = T)
      return(marker_aurocs)
    }else{
      marker_aurocs = get_geneset_preserved_coexp(aggregated_fetal_network, test_network, celltype_markers$gene)
      return(marker_aurocs)
    }
  }
}



#'Compute the preserved co-expression of a cell-type marker gene set from our list of fetal MetaMarkers
#'@param aggregated_fetal_network the aggregated fetal co-expression network to act as a reference network. Any co-expression network will work, with the same gene annotations
#'@param test_network the test co-expression network, we recommend using the rank standardized networks
#'@param metamarkers A datatable of metamarkers, included as data for this package
#'@param celltype string belonging to one of the celltypes included in the MetaMarkers. Must be one of c('Dividing_Progenitor','Neural_Progenitor','Intermediate_Progenitor','GABAergic','Glutamatergic','Non-neuronal')
#'@param num_markers int, number of MetaMarkers, ranked by recurrent differential expression 
#'@param parallel boolean, if the parallelized version should be implemented, requires the parallel R package
#'@param cores integer, the number of computer cores for the parallelized version
#'@return integer, the mean gene preserved co-expression auroc for the geneset
#'
#'@export

get_metaMarker_preserved_score = function(aggregated_fetal_network, test_network, metamarkers, celltype, num_markers, parallel = F, cores = 10){
  
  if(!celltype %in% c('Dividing_Progenitor','Neural_Progenitor','Intermediate_Progenitor','GABAergic','Glutamatergic','Non-neuronal')){
    print("Please pick a cell type in c('Dividing_Progenitor','Neural_Progenitor','Intermediate_Progenitor','GABAergic','Glutamatergic','Non-neuronal') ")
    return(NA)
  }else{
    
    celltype_markers = dplyr::filter(metamarkers, cell_type == celltype & rank <= num_markers)
    
    if(parallel){
      marker_aurocs = get_geneset_preserved_score(aggregated_fetal_network, test_network, celltype_markers$gene, parallel = T)
      return(marker_aurocs)
    }else{
      marker_aurocs = get_geneset_preserved_score(aggregated_fetal_network, test_network, celltype_markers$gene)
      return(marker_aurocs)
    }
  }
}

#'Get a p-value for all GO terms quantifying whether each GO term has less (left-sided p-value) or more (right-sided p-value) preserved co-expression than expected by the length of the GO term
#'@param aggregated_fetal_network the aggregated fetal co-expression network to act as a reference network. Any co-expression network will work, with the same gene annotations
#'@param test_network the test co-expression network, we recommend using the rank standardized networks
#'@param go_matrix one of the 3 binary matrices provided encoding the mappings between all genes and their GO terms
#'#'@param parallel boolean, if the parallelized version should be implemented, requires the parallel R package
#'@return data.frame of the resulting statistics per GO term
#'
#'@export

get_GO_term_scores = function(aggregated_fetal_network, test_network, go_matrix, parallel = F){
  
  #Preserved coexpression auroc for each gene
  if(parallel){all_gene_aurocs = get_geneset_preserved_coexp(aggregated_fetal_network, test_network, go_genes, parallel = T)
  }else{all_gene_aurocs = get_geneset_preserved_coexp(aggregated_fetal_network, test_network, go_genes)}
  
  #Gene set lengths for all go terms
  gene_set_lengths = colSums(go_matrix)
  #Get the description for each GO term
  gene_set_descriptions = sapply(1:length(gene_set_lengths), function(i) 
    ifelse(is.null(GO_descriptions[[names(gene_set_lengths)[i]]]), NA, Term(GO_descriptions[[names(gene_set_lengths)[i]]]) ) )
  
  
  #Average preserved co-expression auroc for each gene set
  gene_set_scores = apply(go_matrix, 2, crossprod, all_gene_aurocs) / gene_set_lengths
  
  #Population mean and sd of aurocs
  pop_mu = mean(all_gene_aurocs)
  pop_sd = sd(all_gene_aurocs)
  
  #Get a p-value through the mean sample error for each GO term
  #Left sided p-value, GO terms with less preserved co-expression than expected
  se_pval_vec_left = vector(mode = 'numeric', length = length(gene_set_scores))
  for(i in 1:length(gene_set_scores)){
    
    sample_se = pop_sd / sqrt(gene_set_lengths[i])
    z_score = (gene_set_scores[i] - pop_mu) / sample_se
    se_pval_vec_left[i] = pnorm(z_score)
  }
  #Right sided p-value, GO terms with more preserved coexpression than expected
  se_pval_vec_right = 1 - se_pval_vec_left
  
  go_term_presCoexp_df = data.frame(go_term = names(gene_set_scores),description = gene_set_descriptions,
                                    gene_set_length = gene_set_lengths, presCoexp = gene_set_scores, 
                                    left_pval = se_pval_vec_left, right_pval = se_pval_vec_right,
                                    left_adj_pval = p.adjust(se_pval_vec_left, method = 'BH'), right_adj_pval = p.adjust(se_pval_vec_right, method = 'BH'))
  go_term_presCoexp_df = go_term_presCoexp_df[order(go_term_presCoexp_df$left_adj_pval), ]
  return(go_term_presCoexp_df)
}




#'Compute a percentile from ecdf
#'@param x distribution to compute ecdf
#'@param perc percentile to report
#'@return double
#'
#'@export

ecdf_fun <- function(x,perc) ecdf(x)(perc)



#'Compute the preserved co-expression scores for the top 100 MetaMarker genes per cell-type for the test network and plot them in reference to our metaanalysis of fetal and organoid datasets
#'@param aggregated_fetal_network the aggregated fetal co-expression network to act as a reference network. Any co-expression network will work, with the same gene annotations
#'@param test_network the test co-expression network, we recommend using the rank standardized networks
#'@param metamarkers A datatable of metamarkers, included as data for this package
#'@param meta_presCoexp_df data frame of MetaMarker scores for the fetal and organoid datasets we assessed in our meta-analysis, included as data for this package
#'@param parallel boolean, if the parallelized version should be implemented, requires the parallel R package
#'@return a list including the MetaMarker scores for the test network and two ggplot objects plotting the results for the fetal and organoid datasets separately
#'
#'@export

plot_meta_results = function(aggregated_fetal_network, test_network, meta_markers, meta_presCoexp_df, parallel = F){
  
  #Get the scores for the fetal metamarker sets of the test network
  test_mat_scores = vector(mode = 'numeric', length = 6)
  test_mat_celltypes = c('fetal Neural progenitor marker','fetal Dividing progenitor marker','fetal Intermediate progenitor marker',
                         'fetal GABAergic marker','fetal glutamatergic marker','fetal non-neuronal marker')
  
  if(parallel){
    test_mat_scores[1] = get_metaMarker_preserved_score(aggregated_fetal_network, test_network, meta_markers, celltype = 'Neural_Progenitor', num_markers = 100, parallel = T)
    test_mat_scores[2] = get_metaMarker_preserved_score(aggregated_fetal_network, test_network, meta_markers, celltype = 'Dividing_Progenitor', num_markers = 100, parallel = T)
    test_mat_scores[3] = get_metaMarker_preserved_score(aggregated_fetal_network, test_network, meta_markers, celltype = 'Intermediate_Progenitor', num_markers = 100, parallel = T)
    test_mat_scores[4] = get_metaMarker_preserved_score(aggregated_fetal_network, test_network, meta_markers, celltype = 'GABAergic', num_markers = 100, parallel = T)
    test_mat_scores[5] = get_metaMarker_preserved_score(aggregated_fetal_network, test_network, meta_markers, celltype = 'Glutamatergic', num_markers = 100, parallel = T)
    test_mat_scores[6] = get_metaMarker_preserved_score(aggregated_fetal_network, test_network, meta_markers, celltype = 'Non-neuronal', num_markers = 100, parallel = T)
  }else{
    test_mat_scores[1] = get_metaMarker_preserved_score(aggregated_fetal_network, test_network, meta_markers, celltype = 'Neural_Progenitor', num_markers = 100)
    test_mat_scores[2] = get_metaMarker_preserved_score(aggregated_fetal_network, test_network, meta_markers, celltype = 'Dividing_Progenitor', num_markers = 100)
    test_mat_scores[3] = get_metaMarker_preserved_score(aggregated_fetal_network, test_network, meta_markers, celltype = 'Intermediate_Progenitor', num_markers = 100)
    test_mat_scores[4] = get_metaMarker_preserved_score(aggregated_fetal_network, test_network, meta_markers, celltype = 'GABAergic', num_markers = 100)
    test_mat_scores[5] = get_metaMarker_preserved_score(aggregated_fetal_network, test_network, meta_markers, celltype = 'Glutamatergic', num_markers = 100)
    test_mat_scores[6] = get_metaMarker_preserved_score(aggregated_fetal_network, test_network, meta_markers, celltype = 'Non-neuronal', num_markers = 100)
  }
  test_df = data.frame(Var1 = rep('test_network', length = 6), gene_celltype_label = test_mat_celltypes, score = test_mat_scores)
  
  #Get the percentiles of the test network's scores, using our meta-analysis of fetal and organoid datasets
  org_fetal_percentile_vec = sapply(1:nrow(test_df), function(i) paste(as.character(round(ecdf_fun(dplyr::filter(meta_presCoexp_df, 
                                                                                                          dataset_type == 'Fetal' & gene_celltype_label == test_df$gene_celltype_label[i])$score,
                                                                                                   dplyr::filter(test_df, gene_celltype_label == test_df$gene_celltype_label[i])$score)*100)), 'percentile'))
  
  org_org_percentile_vec = sapply(1:nrow(test_df), function(i) paste(as.character(round(ecdf_fun(dplyr::filter(meta_presCoexp_df, 
                                                                                                        dataset_type == 'Organoid' & gene_celltype_label == test_df$gene_celltype_label[i])$score,
                                                                                                 dplyr::filter(test_df, gene_celltype_label == test_df$gene_celltype_label[i])$score)*100)), 'percentile'))
  
  test_df$org_fetal_percentile = org_fetal_percentile_vec 
  test_df$org_org_percentile = org_org_percentile_vec 
  #Plot
  g_fetal = ggplot2::ggplot(dplyr::filter(meta_presCoexp_df, dataset_type == 'Fetal'), ggplot2::aes(x = gene_celltype_label, y = score)) + ggplot2::geom_violin(scale = 'width') +
    ggplot2::ylim(0,1) +  ggplot2::ylab('Preserved Co-expression score') +  ggplot2::ggtitle('Fetal datasets') +  ggplot2::xlab('Top 100 Fetal MetaMarker gene sets' ) +
    ggplot2::stat_summary(data = test_df, ggplot2::aes(x = gene_celltype_label, y = score), fun = 'sum', color = 'red', geom = 'crossbar', width = .5) +
    ggplot2::geom_text(data = test_df, ggplot2::aes(label = org_fetal_percentile), y = .25) +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge=2))
  
  
  g_org =  ggplot2::ggplot( dplyr::filter(meta_presCoexp_df, dataset_type == 'Organoid'), ggplot2::aes(x = gene_celltype_label, y = score)) +  ggplot2::geom_violin(scale = 'width') +
    ggplot2::ylim(0,1) +  ggplot2::ylab('Preserved Co-expression score') +  ggplot2::ggtitle('Organoid datasets') +  ggplot2::xlab('Top 100 Fetal MetaMarker gene sets' ) +
    ggplot2::stat_summary(data = test_df, ggplot2::aes(x = gene_celltype_label, y = score), fun = 'sum', color = 'red', geom = 'crossbar', width = .5) +
    ggplot2::geom_text(data = test_df, ggplot2::aes(label = org_org_percentile), y = .25) +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge=2))
  
  return(list(test_df, g_fetal, g_org))
}









