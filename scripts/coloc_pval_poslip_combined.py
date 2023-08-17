import os
import sys
import pickle

sys.path.append("..") # Adds higher directory to python modules path.

import utils
from coloc_utils import *
from config import store_dir, data_dir, date_key, enrichment_dir
from molmass import Formula

# Load data
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

# Combined analysis
combined_ads = {}
combined_ads['combined'] = {}

for tis, ads in tissue_ads.items():
    combined_ads['combined'].update(ads)
    
thr = 1
n_shuffles = 1000
combined_colocs = all_tissue_colocs_mp(tissue_adatas=combined_ads, min_dataset_fraction=0.2, shuffle=False, threads=thr)

tissue_colocs_shuffle_list = [all_tissue_coloc_sampling_mp(combined_colocs, min_dataset_fraction=0.2, 
                                                           threads=thr, save_coloc_dict=True) for i in range(n_shuffles)]

coloc_pvalues_mp(combined_colocs, tissue_colocs_shuffle_list, metric='mean', threads=thr)
coloc_pvalues_mp(combined_colocs, tissue_colocs_shuffle_list, metric='mediqr', threads=thr)

pickle.dump(combined_colocs, 
            open(os.path.join(store_dir, 'pos_lip_combined_colocs.pickle'), "wb"))