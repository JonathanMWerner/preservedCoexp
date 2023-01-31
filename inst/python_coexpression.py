

import pandas as pd
import numpy as np
import scipy.stats as sci

def get_coexpression(gene_exp, gene_names):

  rank_test_py_exp = sci.rankdata(gene_exp, method = 'average', axis = 1)                    #Row ranks
  rank_test_py_exp = rank_test_py_exp - rank_test_py_exp.mean(axis = 1)[1]                        #Center each gene
  rank_test_py_exp = rank_test_py_exp /np.sqrt(np.square(rank_test_py_exp).sum(axis = 1))[:,None] #divide by sqrt(rowSums)
  cr_python = np.dot(rank_test_py_exp, rank_test_py_exp.T)                                        # Get correlations
  cr_python = pd.DataFrame(cr_python, columns = gene_names, index = gene_names )
  
  return cr_python
