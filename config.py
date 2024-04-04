import os
import json

# Base folder where datasets will be downloaded to and results stored.
data_dir = '/g/alexandr/tim/metaspace_evaluation/'

# Key for the current analysis set. 
# Should be the date on which all datasets are downloaded on which the analysis is then performed.
date_key = '230201'

# Folder for storing results
store_dir = os.path.join(data_dir, date_key)
if date_key not in os.listdir(data_dir):
    os.mkdir(store_dir)

# Subfolder for saving enrichment related files
enrichment_key = 'enrichment'
enrichment_dir = os.path.join(store_dir, enrichment_key)
if enrichment_key not in os.listdir(store_dir):
    os.mkdir(enrichment_dir)

# Subfolder for saving module analysis related files
module_key = 'module_analysis'
module_dir = os.path.join(store_dir, module_key)
if module_key not in os.listdir(store_dir):
    os.mkdir(module_dir)
    
# Subfolder for saving module analysis related files
figuredata_key = 'figure_data'
figuredata_dir = os.path.join(store_dir, figuredata_key)
if figuredata_key not in os.listdir(store_dir):
    os.mkdir(figuredata_dir)

# For R scripts to have access to the data directories
analysis_config_dict = {'date_key': date_key,
                        'data_dir': data_dir, 
                        'store_dir': store_dir, 
                        'enrichment_key': enrichment_key, 
                        'enrichment_dir': enrichment_dir, 
                        'figuredata_key': figuredata_key,
                        'figuredata_dir': figuredata_dir,
                        'module_key': module_key,
                        'module_dir': module_dir}

with open("analysis_config.json", "w") as outfile:
    json.dump(analysis_config_dict, outfile)
