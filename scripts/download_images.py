import os
from tqdm import tqdm
import pickle

from metaspace import SMInstance
from metaspace_converter import dataset_to_anndata

sys.path.append("..") # Adds higher directory to python modules path.
from config import store_dir, data_dir, date_key, enrichment_dir

# List of datasets to download
nm = pickle.load(open(os.path.join(store_dir, 'DatasetsForDownload_Neg_Met.pickle'), "rb" ) )
nl = pickle.load(open(os.path.join(store_dir, 'DatasetsForDownload_Neg_Lip.pickle'), "rb" ) )
pm = pickle.load(open(os.path.join(store_dir, 'DatasetsForDownload_Pos_Met.pickle'), "rb" ) )
pl = pickle.load(open(os.path.join(store_dir, 'DatasetsForDownload_Pos_Lip.pickle'), "rb" ) )
    
dss_names = list(set(nm+nl+pm+pl))


database = ('HMDB', 'v4')
save_path = os.path.join(store_dir, 'all_ionimages/')

counter = 0

sm = SMInstance()

for dsn in dss_names:
    
    ds = sm.dataset(id=dsn)
    
    if ds.status == 'FINISHED' and f'{dsn}.pickle' not in os.listdir(save_path):
        if database in [(x.name, x.version) for x in ds.database_details]:
            try:
                tmp_adata = dataset_to_anndata(ds, fdr=0.2, database=database)
            
                pickle.dump(tmp_adata,
                open(os.path.join(save_path, f'{dsn}.pickle'), "wb"))

            except Exception as e:
                print(f'Error in dataset {dsn}:')
                print(e)
                print('====')
    
    print(f'{counter}/{len(dss_names)}')
    
    counter += 1

