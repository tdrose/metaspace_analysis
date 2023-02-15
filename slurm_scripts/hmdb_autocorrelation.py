import os
from tqdm import tqdm
import pickle
import random

import pandas as pd
import numpy as np
import scipy
from sklearn.metrics.pairwise import pairwise_kernels
import matplotlib.pyplot as plt

from metaspace import SMInstance
from anndata import AnnData
from metaspace2anndata import dataset_to_anndata
import scanpy as sc
import squidpy as sq


date_key = '230201'

# Store in Alexandrov g drive
data_dir = '/g/alexandr/tim/metaspace_evaluation/'

store_dir = os.path.join(data_dir, date_key)

database = ('HMDB', 'v4')

filename = 'hmdb4_autocorrelation.pickle'

# Load dataset
adata_mol = pickle.load(open( os.path.join(store_dir, 'all_datasets_mol_anndata.pickle'), "rb" ))

sc.pp.filter_genes(adata_mol, min_cells=200)
sc.pp.filter_cells(adata_mol, min_genes=50)
sc.pp.normalize_total(adata_mol, target_sum=1e4)

out_dict = {}

sm = SMInstance()
counter = 0


# Loop over all remaining datasets after filtering
for ds_id in adata_mol.obs.index:
    
    ds = sm.dataset(id=ds_id)
    
    tmp_adata = dataset_to_anndata(ds, fdr=0.5, database=database)
    
    sq.gr.spatial_neighbors(tmp_adata, coord_type='grid')
    tmp_adata.obsp['connectivities'] = tmp_adata.obsp['spatial_connectivities']
    
    out_dict[ds_id] = {
        'xdim': max(tmp_adata.obs['x']) + 1,
        'ydim': max(tmp_adata.obs['y']) + 1,
        'nfeatures': len(tmp_adata.var.index),
        'autocorrelation': {k: v for v, k in zip(sc.metrics.morans_i(tmp_adata), tmp_adata.var.index)}
    }
    
    if (counter % 100) == 0:
                pickle.dump(out_dict, 
            open( os.path.join('/scratch/trose/tmp', str(counter) + '_' + filename), "wb" ) )
    
    counter += 1
    
pickle.dump(out_dict, 
            open( os.path.join(store_dir, filename), "wb" ) )
    

# Output: {ds_id: {xdim: 0, ydim: 0, nfeatures: 0, tab: pd.DataFrame}}