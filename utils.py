import pandas as pd
import numpy as np
import matplotlib.pyplot as plt




def make_metadata_dict(dss, results_dict, only_results=False):
    
    metadata_dict = {'Organism': {}, 
                     'Condition': {}, 
                     'Organism_Part': {}, 
                     'Polarity': {},
                     'maldi_matrix': {},
                     'Group': {},
                     'mzmin': {},
                     'mzmax': {}
                    }
    
    for d in dss:
        if only_results:
            if d.id in results_dict.keys():
                metadata_dict['Organism'][d.id] = d.metadata['Sample_Information']['Organism'].strip().lower()
                metadata_dict['Condition'][d.id] = d.metadata['Sample_Information']['Condition'].strip().lower()
                metadata_dict['Organism_Part'][d.id] = d.metadata['Sample_Information']['Organism_Part'].strip().lower()
                metadata_dict['Polarity'][d.id] = d.metadata['MS_Analysis']['Polarity'].strip().lower()
                metadata_dict['maldi_matrix'][d.id] = d.metadata['Sample_Preparation']['MALDI_Matrix']

                if d.group is None:
                    metadata_dict['Group'][d.id] = "not available"
                else:
                    metadata_dict['Group'][d.id] = d.group['shortName'].strip()

                metadata_dict['mzmin'][d.id] = results_dict[d.id]['mz'].min()
                metadata_dict['mzmax'][d.id] = results_dict[d.id]['mz'].max()
                
        else:
            metadata_dict['Organism'][d.id] = d.metadata['Sample_Information']['Organism'].strip().lower()
            metadata_dict['Condition'][d.id] = d.metadata['Sample_Information']['Condition'].strip().lower()
            metadata_dict['Organism_Part'][d.id] = d.metadata['Sample_Information']['Organism_Part'].strip().lower()
            metadata_dict['Polarity'][d.id] = d.metadata['MS_Analysis']['Polarity'].strip().lower()
            metadata_dict['maldi_matrix'][d.id] = d.metadata['Sample_Preparation']['MALDI_Matrix']

            if d.group is None:
                metadata_dict['Group'][d.id] = "not available"
            else:
                metadata_dict['Group'][d.id] = d.group['shortName'].strip()

            if d.id in results_dict.keys():
                metadata_dict['mzmin'][d.id] = results_dict[d.id]['mz'].min()
                metadata_dict['mzmax'][d.id] = results_dict[d.id]['mz'].max()
            else:
                metadata_dict['mzmin'][d.id] = np.nan
                metadata_dict['mzmax'][d.id] = np.nan
        
    return metadata_dict


def metadata_dict_totable(md):
    
    convert_dict = {}
    for key in md.keys():
        convert_dict[key] = pd.Series(md[key])
        
    return pd.DataFrame(convert_dict)


def clean_metadata_table(mdt):
    
    out = mdt.copy()
    
    out['Organism'] = out['Organism'].replace({
        'homospaiens': 'homo sapiens (human)',
        'mouse': 'mus musculus (mouse)',
        'human': 'homo sapiens (human)',
        'none': 'NA',
        'n/a': 'NA',
        'rat': 'rattus norvegicus (rat)',
        'mice': 'mus musculus (mouse)'
    })
    
    out['Condition'] = out['Condition'].replace({
        'wildtype': 'normal',
        'wtype': 'normal',
        'healthy': 'normal'
    })
    
    out['maldi_matrix'] = out['maldi_matrix'].replace({
        '2,5-dihydroxybenzoic acid (DHB)': 'DHB',
        'none': 'None',
        'alpha-cyano-4-hydroxycinnamic acid (CHCA)': 'CHCA',
        'N/A': 'NA',
        'n-(1-naphthyl)ethylenediamine dihydrochloride (NEDC)' : 'NEDC',
        '9-aminoacridine (9AA)': '9AA',
        '2,5-dihydroxyacetophenone (DHA)': 'DHA',
        '1,5-diaminonaphthalene (DAN)': 'DAN',
        '1,8-bis(dimethylamino)naphthalene (DMAN)' : 'DMAN',
        'ice': 'Ice',
        'α-Cyano-4-hydroxycinnamic acid (HCCA)': 'HCCA',
        '2,5-dihydroxyacetophenone (DHAP)': 'DHAP',
        'trans-2-[3-(4-tert-Butylphenyl)-2-methyl-2-propenylidene] malononitrile (DCTB)': 'DCTB',
        '2,5-dihydroxy benzoic acid (DHB)': 'DHB',
        '2,4,6-Trihydroxyacetophenone (THAP)': 'THAP',
        '1,5-diaminonaphthalene (DAN) | DHB': 'DAN, DHB',
        '1,5-diaminonaphthalene (DAN)+HCl': 'DAN, HCl',
        '1,5-diaminonaphthalene (DAN), DHB': 'DAN, DHB',
        '2-Mercaptobenzothiazole (MBT)': 'MBT',
        '2,5-dihydroxyacetophenone (DHAP)/Gold Sputter': 'DHAP, Au',
        '2-Mercaptobenzothiazole (MBT)': 'MBT',
        '2′,6′-Dihydroxyacetophenone (DHAP)': 'DHAP',
        '5-chloro-2-Mercaptobenzothiazole (CMBT)': 'CMBT',
        '2,5-dihydroxybenzoic acid (DHB)/Fe3O4 binary mix': 'DHB, Fe3O4',
        'DHB and DAN': 'DAN, DHB',
        'DHB, CHCA': 'DHB, CHCA',
        'DHB/CHCA': 'DHB, CHCA',
        'DHB/Fe3O4': 'DHB, Fe3O4',
        'norharmane': 'Norharmane',
        'norhramne': 'Norharmane'
    })
    
    return out


def flatten(l):
    return [item for sublist in l for item in sublist]

def top_feature_col(in_col: pd.Series, top: int=20, other_key: str='Other'):
    tmp = pd.value_counts(in_col)
    replace_dict = {}
    for i in range(len(tmp)):
        if i < top:
            replace_dict[tmp.index[i]] = tmp.index[i]
        else:
            replace_dict[tmp.index[i]] = other_key
    return in_col.replace(to_replace=replace_dict)


def plot_cluster_metadata(adata, cluster='0', figsize=(9,8)):
    
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) =  plt.subplots(ncols=3, nrows=2, figsize=figsize)
    
    cols = ['top_Polarity', 'top_maldi_matrix', 'top_Group', 
                             'top_Organism', 'top_Organism_Part']
    
    color = adata.uns['leiden_colors'][int(cluster)]
    
    adata.obs[adata.obs['leiden']==cluster][cols].groupby('top_maldi_matrix')['top_maldi_matrix'].count().plot(kind='bar', stacked=True, ax=ax1, color=color)
    adata.obs[adata.obs['leiden']==cluster][cols].groupby('top_Polarity')['top_Polarity'].count().plot(kind='bar', stacked=True, ax=ax2, color=color)
    adata.obs[adata.obs['leiden']==cluster][cols].groupby('top_Group')['top_Group'].count().plot(kind='bar', stacked=True, ax=ax3, color=color)
    adata.obs[adata.obs['leiden']==cluster][cols].groupby('top_Organism')['top_Organism'].count().plot(kind='bar', stacked=True, ax=ax4, color=color)
    adata.obs[adata.obs['leiden']==cluster][cols].groupby('top_Organism_Part')['top_Organism_Part'].count().plot(kind='bar', stacked=True, ax=ax5, color=color)
    
    plt.show()
    
    