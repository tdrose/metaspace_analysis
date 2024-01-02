import os
from tqdm import tqdm
import pickle
import random
from typing import List, Dict, Tuple
from anndata import AnnData, concat
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels
import umap

from molmass import Formula
import metaspace
from scipy import stats

from collections import defaultdict

import sys
sys.path.append("..")
import utils
from coloc_utils import *
from config import store_dir, data_dir, date_key, enrichment_dir, module_dir

tissue_colocs = pickle.load(open(os.path.join(store_dir, 'pos_lip_tissue_colocs.pickle'), "rb" ))
df = pd.concat([x['c_measures'].assign(tissue=tis) for tis, x in tissue_colocs.items()])
df['mean_sig'] = df['pval_mean_corr'] <= 0.05
df['mediqr_sig'] = df['pval_mediqr_corr'] <= 0.05

pos_lip_top_datasets = pickle.load(open(os.path.join(store_dir, 'pos_lip_top_datasets_list.pickle'), "rb" ))
tissue_ads = load_alltissue_datasets(pos_lip_top_datasets)

# exclude datasets
del tissue_ads['Brain']['2022-08-24_00h20m06s']
del tissue_ads['Brain']['2022-08-23_23h48m59s']
del tissue_ads['Buccal mucosa']

# Calculate mass of every molecule
for tis, dsd in tissue_ads.items():
    for dsid, ds in dsd.items():
        tissue_ads[tis][dsid].var['mass'] = [Formula(x).mass for x in tissue_ads[tis][dsid].var['formula'].values]

# Filter by mass
for tis, dsd in tissue_ads.items():
    for dsid, ds in dsd.items():
        tmp = tissue_ads[tis][dsid]
        tmp = tmp[:, tmp.var['mass'] <= 900]
        tissue_ads[tis][dsid] = tmp[:, tmp.var['mass'] >= 400]

# Remove isobars
tissue_ads = mark_isobars(tissue_ads, ppm_threshold=3)
ll = []
for tis, dsd in tissue_ads.items():
    for dsid, ds in dsd.items():
        if ds[:, ~ds.var['has_isobar']].shape[1] < 10:
            ll.append((tis, dsid))
        else:
            # Remove too big datasets
            if ds.shape[0] > 100000:
                ll.append((tis, dsid))
for item in ll:
    del tissue_ads[item[0]][item[1]]

results = pickle.load(open(os.path.join(store_dir, 'hmdb4_results_IsobarFree.pickle'), "rb" ) )
dss = pickle.load(open(os.path.join(store_dir, 'all_datasets.pickle'), "rb" ) )

md = utils.make_metadata_dict(dss, results, only_results=True)
mdt = utils.clean_metadata_table(utils.metadata_dict_totable(md), path='../metadata_mapping/')
mdt['top_Organism'] = utils.top_feature_col(mdt['Organism'], top=6, exclusion_list=['N/A', 'cultured cells', 'Multiple'])
mdt['top_Condition'] = utils.top_feature_col(mdt['Condition'], top=10)
mdt['top_Organism_Part'] = utils.top_feature_col(mdt['Organism_Part'], top=10)
mdt['top_Polarity'] = utils.top_feature_col(mdt['Polarity'], top=10)
mdt['top_maldi_matrix'] = utils.top_feature_col(mdt['maldi_matrix'], top=8)
mdt['top_Group'] = utils.top_feature_col(mdt['Group'], top=10, exclusion_list=['not available'])





# Single pixel specific analysis
ads = {}
for tis, add in tissue_ads.items():
    # Exclude tissues
    if tis not in ['Epididymis', 'Ovary']:
        for dsid, ad in add.items():
            ads[dsid] = molecule_adata(ad, mdt)
    
# molecule frequencies
molfreq=defaultdict(float)
for ad in ads.values():
    for mol in ad.var['formula']:
        molfreq[mol] += 1
        
        
mol_cutoff = 0.2
absmc = int(mol_cutoff*len(ads))

featurelist = [key for key, val in molfreq.items() if val>=absmc ]

dssizes = []
for ad in ads.values():
    counter = 0
    for mol in ad.var['formula']:
        if mol in featurelist:
            counter +=1
    dssizes.append(counter)
    
# Subset datasets
ads_for_concat = []
for ad in ads.values():
    tmp = ad.copy()
    tmp.var['formula2'] = tmp.var['formula']
    tmp.var = tmp.var.set_index('formula')
    feats = tmp.var['formula2'].isin(featurelist)
    # Save size info
    tmp.var = tmp.var.assign(features=sum(feats))
    ads_for_concat.append(tmp[:, feats])
# concat datasets
adc = concat(ads_for_concat, join='outer')
adc.X[np.isnan(adc.X)] = 0

sc.pp.filter_cells(adc, min_genes=30)

# adc = adc[adc.obs['organism'] != 'Human', :]
adc = adc[adc.obs['organism'] != 'Multiple', :]
adc = adc[adc.obs['ds']!='2017-06-09_07h12m31s', :]
remove_ds = list(adc.obs['ds'].value_counts()[adc.obs['ds'].value_counts()<3].index)
adc = adc[~adc.obs['ds'].isin(remove_ds)]

sc.pp.normalize_total(adc, target_sum=1e4)
sc.pp.log1p(adc)

sc.pp.pca(adc)
sc.external.pp.bbknn(adc, batch_key='organ', metric='euclidean', neighbors_within_batch=3)
#sc.pp.neighbors(adc, metric='cosine')
sc.tl.umap(adc)

sc.tl.leiden(adc)
pickle.dump(adc, 
            open(os.path.join(store_dir, 'single_pixel_adata_combined_bbknn2.pickle'), "wb"))




