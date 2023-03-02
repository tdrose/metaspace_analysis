import os

date_key = '230201'

# Store in Alexandrov g drive
data_dir = '/g/alexandr/tim/metaspace_evaluation/'



store_dir = os.path.join(data_dir, date_key)

if date_key not in os.listdir(data_dir):
    os.mkdir(store_dir)
    
enrichment_key = 'enrichment'

enrichment_dir = os.path.join(store_dir, enrichment_key)
    
if enrichment_key not in os.listdir(store_dir):
    os.mkdir(enrichment_dir)
