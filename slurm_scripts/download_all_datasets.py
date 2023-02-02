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


date_key = '230201'

# Store in Alexandrov g drive
data_dir = '/g/alexandr/tim/metaspace_evaluation/'

store_dir = os.path.join(data_dir, date_key)

if date_key not in os.listdir(data_dir):
    os.mkdir(store_dir)
    

dss = pickle.load(open(os.path.join(store_dir, 'all_datasets.pickle'), "rb" ) )

database = ('HMDB', 'v4')
filename = 'hmdb4_results.pickle'

all_results_dict = {}

counter = 0


for ds in dss:
    if ds.status == 'FINISHED' and ds.id not in all_results_dict.keys():
        if database in [(x.name, x.version) for x in ds.database_details]:
            
            # Download annotation table
            tmp_tab = ds.results(database=database)
            
            if tmp_tab.size > 0:
                if ds.id not in all_results_dict.keys():
                    all_results_dict[ds.id] = tmp_tab[['ionFormula', 'ion', 'fdr', 'mz', 'offSample', 'moleculeNames', 'intensity', 'moleculeIds']]
                    counter +=1
            print(counter)

            # Save intermediate results    
            if (counter % 100) == 0:
                pickle.dump(all_results_dict, 
            open( os.path.join('/scratch/trose/tmp', filename + '_' + str(counter)), "wb" ) )
                
pickle.dump(all_results_dict, 
            open( os.path.join(store_dir, filename), "wb" ) )

