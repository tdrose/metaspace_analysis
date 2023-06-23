import os
import sys
import pickle

sys.path.append("..") # Adds higher directory to python modules path.

import utils
from coloc_utils import *
from config import store_dir, data_dir, date_key, enrichment_dir


# Load data
pos_lip_top_datasets = pickle.load(open(os.path.join(store_dir, 'pos_lip_top_datasets_list.pickle'), "rb" ))
tissue_ads = load_alltissue_datasets(pos_lip_top_datasets)

# exclude datasets
del tissue_ads['Brain']['2022-08-24_00h20m06s']
del tissue_ads['Brain']['2022-08-23_23h48m59s']
del tissue_ads['Buccal mucosa']

thr = 6

n_shuffles = 1000
tissue_colocs_shuffle_list = [all_tissue_colocs_mp(tissue_adatas=tissue_ads, 
                                                   min_dataset_fraction=0.2, 
                                                   shuffle=True, 
                                                   threads=thr) for x in range(n_shuffles)]

tissue_colocs = all_tissue_colocs_mp(tissue_adatas=tissue_ads, min_dataset_fraction=0.2, shuffle=False, threads=thr)

coloc_pvalues_mp(tissue_colocs, tissue_colocs_shuffle_list, metric='mean', threads=thr)
coloc_pvalues_mp(tissue_colocs, tissue_colocs_shuffle_list, metric='mediqr', threads=thr)

pickle.dump(tissue_colocs, 
            open(os.path.join(store_dir, 'pos_lip_tissue_colocs.pickle'), "wb"))