import scanpy as sc
from rpy2.robjects import r, globalenv, pandas2ri
import anndata2ri
import pickle

import sys
sys.path.append("..")
import utils
from coloc_utils import *
from config import store_dir, data_dir, date_key, enrichment_dir, module_dir

pandas2ri.activate()
anndata2ri.activate()

from rpy2.robjects.packages import importr
seurat = importr('Seurat')

# Function for seurat integration
def process_data(adata, file, path, batch_key='ds', reduction='rpca'):
    
    # Prepare data for Seurat
    adata_seurat = adata.copy()
    # Data is already in log
    adata_seurat.obs = adata_seurat.obs.reset_index(drop=True)
    adata_seurat.layers["counts"] = np.exp(adata_seurat.X.copy())-1
    adata_seurat.layers["logcounts"] = adata_seurat.X.copy()
    adata_seurat.obs[batch_key] = adata_seurat.obs[batch_key].astype(str)
    
    # Remove datasets with too few cells.
    tmp = adata_seurat.obs[batch_key].value_counts()
    tmp2 = list(tmp[tmp<=adata_seurat.shape[1]].index)
    adata_seurat = adata_seurat[~adata_seurat.obs[batch_key].isin(tmp2), :]
    
    adata_seurat.obs_names_make_unique()
    del adata_seurat.uns

    # Do integration in R
    globalenv["batch_key"] = batch_key
    globalenv["adata_seurat"] = adata_seurat
    globalenv["n_features"] = adata_seurat.shape[1]
    r(f'seurat <- as.Seurat(adata_seurat, counts="counts", data="logcounts")')
    r('batch_list <- SplitObject(seurat, split.by = batch_key)')
    r('batch_list <- lapply(X = batch_list, FUN = function(x) {x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = n_features)})')
    r('anc_feats = rownames(batch_list[[1]])')
    r('batch_list <- lapply(X = batch_list, FUN = function(x) {x <- RunPCA(ScaleData(x, features = anc_feats, verbose = FALSE), features = anc_feats, verbose = FALSE, approx = F, npcs = 30)})')
    r('print(sapply(batch_list, function(x){ncol(x@reductions$pca)}))')
    # r(f'saveRDS(batch_list, file = "/g/alexandr/tim/brain_batch_list.rds")')
    r(f'anchors <- FindIntegrationAnchors(batch_list, anchor.features = anc_feats, reduction="{reduction}", k.filter = 50, scale=FALSE)')
    r('integrated <- IntegrateData(anchors)')
    # Extract the integrated expression matrix
    r('integrated_expr <- GetAssayData(integrated)')
    # Make sure the rows and columns are in the same order as the original object
    r('integrated_expr <- integrated_expr[rownames(seurat), colnames(seurat)]')
    # Transpose the matrix to AnnData format
    integrated_expr = r('t(integrated_expr)')
    
    pickle.dump(integrated_expr, open(os.path.join(path, file.split('.')[0]+'_matrix.pickle'), "wb"))
    pickle.dump(adata_seurat, open(os.path.join(path, file.split('.')[0]+'_tmp.pickle'), "wb"))
    # Do the rest in Python
    adata_seurat.X = integrated_expr
    adata_seurat.layers["seurat"] = integrated_expr
    adata_seurat.obs[batch_key] = adata_seurat.obs[batch_key].astype('category')
    
    del adata_seurat.obsm
    del adata_seurat.obsp
    
    sc.tl.pca(adata_seurat)
    sc.pp.neighbors(adata_seurat, metric='cosine')
    sc.tl.umap(adata_seurat)
    sc.tl.leiden(adata_seurat)

    return adata_seurat

file = sys.argv[1]

adc = pickle.load(open(os.path.join(store_dir, file), "rb"))

adc_seurat = process_data(adc, file=file, path=store_dir, reduction='cca')

pickle.dump(adc_seurat, 
            open(os.path.join(store_dir, file.split('.')[0]+'cca_seurat.pickle'), "wb"))