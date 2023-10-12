
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import nsforest as ns
import pickle

adata = ad.read_h5ad("/nfs/data/tcell_deconvolution/data/bassez/nsforest_prepared_obj_sampled.h5ad")

nsf_results = ns.NSForest(adata, "groups")

filehandler = open("/nfs/data/tcell_deconvolution/data/bassez/nsforest_result.pickle", "wb")
pickle.dump(nsf_results, filehandler)
