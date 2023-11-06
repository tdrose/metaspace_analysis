import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from anndata import AnnData
import os
from typing import Tuple
from scipy import ndimage
import math


def make_metadata_dict(dss, results_dict, only_results=False, fdrcutoff=0.1):
    
    metadata_dict = {'Organism': {}, 
                     'Condition': {}, 
                     'Organism_Part': {}, 
                     'Polarity': {},
                     'maldi_matrix': {},
                     'Group': {},
                     'mzmin': {},
                     'mzmax': {},
                     'Analyzer': {},
                     'Ionisation_Source': {}
                    }
    
    for d in dss:
        if only_results:
            if d.id in results_dict.keys():
                metadata_dict['Organism'][d.id] = d.metadata['Sample_Information']['Organism'].strip()
                metadata_dict['Condition'][d.id] = d.metadata['Sample_Information']['Condition'].strip()
                metadata_dict['Organism_Part'][d.id] = d.metadata['Sample_Information']['Organism_Part'].strip()
                metadata_dict['Polarity'][d.id] = d.metadata['MS_Analysis']['Polarity'].strip()
                metadata_dict['Analyzer'][d.id] = d.metadata['MS_Analysis']['Analyzer'].strip()
                metadata_dict['Ionisation_Source'][d.id] = d.metadata['MS_Analysis']['Ionisation_Source'].strip()
                metadata_dict['maldi_matrix'][d.id] = d.metadata['Sample_Preparation']['MALDI_Matrix'].strip()

                if d.group is None:
                    metadata_dict['Group'][d.id] = "not available"
                else:
                    metadata_dict['Group'][d.id] = d.group['shortName'].strip()
                    
                tmp = results_dict[d.id][results_dict[d.id]['fdr'] <= fdrcutoff]
                if tmp.shape[0] > 0:
                    metadata_dict['mzmin'][d.id] = tmp['mz'].min()
                    metadata_dict['mzmax'][d.id] = tmp['mz'].max()
                else:
                    metadata_dict['mzmin'][d.id] = np.nan
                    metadata_dict['mzmax'][d.id] = np.nan
        else:
            metadata_dict['Organism'][d.id] = d.metadata['Sample_Information']['Organism'].strip()
            metadata_dict['Condition'][d.id] = d.metadata['Sample_Information']['Condition'].strip()
            metadata_dict['Organism_Part'][d.id] = d.metadata['Sample_Information']['Organism_Part'].strip()
            metadata_dict['Polarity'][d.id] = d.metadata['MS_Analysis']['Polarity'].strip()
            metadata_dict['Analyzer'][d.id] = d.metadata['MS_Analysis']['Analyzer'].strip()
            metadata_dict['Ionisation_Source'][d.id] = d.metadata['MS_Analysis']['Ionisation_Source'].strip()
            metadata_dict['maldi_matrix'][d.id] = d.metadata['Sample_Preparation']['MALDI_Matrix']

            if d.group is None:
                metadata_dict['Group'][d.id] = "not available"
            else:
                metadata_dict['Group'][d.id] = d.group['shortName'].strip()

            if d.id in results_dict.keys():
                tmp = results_dict[d.id][results_dict[d.id]['fdr'] <= fdrcutoff]
                if tmp.shape[0] > 0:
                    metadata_dict['mzmin'][d.id] = tmp['mz'].min()
                    metadata_dict['mzmax'][d.id] = tmp['mz'].max()
                else:
                    metadata_dict['mzmin'][d.id] = np.nan
                    metadata_dict['mzmax'][d.id] = np.nan
            else:
                metadata_dict['mzmin'][d.id] = np.nan
                metadata_dict['mzmax'][d.id] = np.nan
        
    return metadata_dict


def metadata_dict_totable(md):
    
    convert_dict = {}
    for key in md.keys():
        convert_dict[key] = pd.Series(md[key])
        
    return pd.DataFrame(convert_dict)

def md_mapping_dict(path, filename):
    df = pd.read_csv(os.path.join(path, filename), na_filter = False)
    return dict(zip(df['Old'], df['New']))


def clean_metadata_table(mdt, path='./metadata_mapping/'):
    
    out = mdt.copy()
    
    out['Organism'] = out['Organism'].replace(md_mapping_dict(path, 'mapping_organism.csv'))
    out['maldi_matrix'] = out['maldi_matrix'].replace(md_mapping_dict(path, 'mapping_maldimatrix.csv'))
    out['Analyzer'] = out['Analyzer'].replace(md_mapping_dict(path, 'mapping_analyzer.csv'))
    out['Ionisation_Source'] = out['Ionisation_Source'].replace(md_mapping_dict(path, 'mapping_source.csv'))
    out['Organism_Part'] = out['Organism_Part'].replace(md_mapping_dict(path, 'mapping_organismpart.csv'))
    
    return out


def flatten(l):
    return [item for sublist in l for item in sublist]

def top_feature_col(in_col: pd.Series, top: int=20, other_key: str='Other', exclusion_list: list=None):
    tmp = pd.value_counts(in_col)
    # print(tmp[:15])
    top_list = list(tmp.index)
    if exclusion_list is None:
        exclusion_list = []
    
    replace_dict = {}
    # print(tmp[:15])
    for i in range(len(tmp)):
        if i < top and top_list[i] not in exclusion_list:
            replace_dict[top_list[i]] = top_list[i]
        else:
            replace_dict[top_list[i]] = other_key
            if top_list[i] in exclusion_list:
                top += 1
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
    
    
def make_ion_anndata(results, mdt, fdr_cutoff=0.5, only_onSample=False):
    
    features = []
    for tab in tqdm(results.values()):
        tmp_tab = tab[tab['fdr'] <= fdr_cutoff]
        
        if only_onSample:
            tmp_tab = tmp_tab[tmp_tab['offSample'] == False]
        
        for ix in tmp_tab['ion']:
                features.append(ix)
    features = list(set(features))
    
    print(len(features), ' features')
    
    fdr_data = pd.DataFrame(0, columns=list(set(features)), index=results.keys(), dtype='float64')
    
    for i in tqdm(results.keys()):
    # It is late, I lost my creativity for variable names
        tmp_tab = results[i][results[i]['fdr'] <= fdr_cutoff
                            ]
        if only_onSample:
            tmp_tab = tmp_tab[tmp_tab['offSample'] == False]

        ttt = tmp_tab.reset_index()[['ion', 'intensity']]
        ttt2 = ttt.groupby('ion').sum()

        fdr_data.loc[i, ttt2.index] = ttt2['intensity'].values
        
    return AnnData(X=fdr_data.to_numpy(), var=pd.DataFrame(features), obs=mdt.loc[fdr_data.index, :])


def make_molecule_anndata(results, mdt, fdr_cutoff=0.5, only_onSample=False):
    
    mol_features = []
    for tab in tqdm(results.values()):
        tmp_tab = tab[tab['fdr'] <= fdr_cutoff]
        
        if only_onSample:
            tmp_tab = tmp_tab[tmp_tab['offSample'] == False]
        
        for ix in tmp_tab.reset_index()['formula']:
                mol_features.append(ix)
    mol_features = list(set(mol_features))
    
    print(len(mol_features), ' features')
    
    mol_data = pd.DataFrame(0, columns=mol_features, index=results.keys(), dtype='float64')
    
    # Fill dataframe
    for i in tqdm(results.keys()):
        
        # It is late, I lost my creativity for variable names
        tmp_tab = results[i][results[i]['fdr'] <= fdr_cutoff]
        
        if only_onSample:
            tmp_tab = tmp_tab[tmp_tab['offSample'] == False]

        ttt = tmp_tab.reset_index()[['formula', 'intensity']]
        ttt2 = ttt.groupby('formula').sum()
        
        #if i =='2022-12-15_15h53m33s':
        #    print(ttt2.index)
        #    print(mol_data.loc[i, ttt2.index])

        mol_data.loc[i, ttt2.index] = ttt2['intensity'].values
        
        #if i =='2022-12-15_15h53m33s':
        #    print(mol_data.loc[i, ttt2.index])
    #print(mol_data.loc[i, 'C62H120O17P2'])        
    return AnnData(X=mol_data.to_numpy(), var=pd.DataFrame(mol_features), obs=mdt.loc[mol_data.index, :])



def get_hmdb_names(db_tab: pd.DataFrame, formula: str):
    tmp = db_tab.loc[formula, 'name']
    if type(tmp) == str:
        return [tmp]
    else:
        return list(tmp)
    
    return names


def top_annotations(tup, db, top=10, n=4, is_ion=False):
    for i in range(top):
        mol = tup[i][0]
        if is_ion:
            mol = mol.split('+')[0].split('-')[0]
        if mol in db.index:
            print(mol, ' - ', str(get_hmdb_names(db, mol)[:n]))
        else:
            print(mol)

def annotations_list(formula_list, db, n=4, is_ion=False):
    for i in range(len(formula_list)):
        mol = formula_list[i]
        if is_ion:
            mol = mol.split('+')[0].split('-')[0]
        if mol in db.index:
            print(mol, ' - ', str(get_hmdb_names(db, mol)[:n]))
        else:
            print(mol)
            
def get_sig_molecules(adata, rg='ranked_genes', max_mols=None, pval_cutoff=0.01):
    
    pvals = [x < pval_cutoff for x,y in adata.uns[rg]['pvals_adj']]
    
    names = np.array([x for x,y in adata.uns[rg]['names']])
    
    if max_mols is None:
        return names[pvals]
    else:
        return names[pvals][:max_mols]
    
def identifications(adata, sig_molecules, obsv):
    sub = adata[:, sig_molecules]
    
    l1 = []
    l2 = []
    l3 = []
    l4 = []
    
    # Loop over molecules
    for mol in sig_molecules:
        
        data_vec = sub[:, mol].X.transpose()[0]
        #print(data_vec.shape)
        
        # Loop over categories
        for cat in sub.obs[obsv].cat.categories:
            #print(sub.obs.shape)
            tmp = data_vec[sub.obs[obsv]==cat]
            zeros = sum(tmp==0)
            nonzeros = sum(tmp!=0)
            
            l1.append(mol)
            l2.append(cat)
            l3.append(zeros)
            l4.append(nonzeros)
            
    ratios = pd.DataFrame({'mol': l1, 'category': l2, 'zeros': l3, 'nonzeros': l4})
    ratios['ratio'] = ratios['nonzeros'] / ratios['zeros']
            
    return ratios


def tissue_prototyping(adat: AnnData, mol_freq_cutoff: float=0.1, top_ds_cutoff: float=0.2, mol_freq2_cutoff: float=0.1) -> Tuple[pd.Series, pd.Series, pd.Series]:
    """
    Find most characteristic datasets & molecular features per tissue.
    
    :param adat: AnnData object of one tissue
    :param mol_freq_cutoff: Cut-off for molecule frequency (top ... fraction)
    :param top_ds_cutoff: Cut-off for datasets (top ... fraction)
    :return: mol_freq2 - Dataset frequency (using top datasets [top_ds_cutoff]) for top molecules (Molecules after mol_freq2_cutoff is applied), 
             top_datasets - Molecule frequency (using top molecules [mol_freq_cutoff]) for datasets (Datasets after top_ds_cutoff is applied), 
             mol_freq - Dataset frequency (using all datasets in adat) for molecules (Molecules after mol_freq_cutoff is applied).
    """
    
    # Compute frequency of Features
    mol_freq = pd.Series((adat.X > 0).sum(axis=0), index=adat.var.index) / adat.X.shape[0]

    plt.hist(mol_freq)
    plt.show()
    
    # Find cutoff of top (mol_freq_cutoff) fraction of freatures
    com = mol_freq.sort_values(ascending=False)[int(len(mol_freq)*mol_freq_cutoff)]
    
    # With these features find which datasets have the most of these top features
    top_datasets = pd.Series((adat.X[:, mol_freq >= com] > 0).sum(axis=1), index=adat.obs.index)
    
    # Find cutoff of top (top_ds_cutoff) fraction of datasets
    co = top_datasets.sort_values(ascending=False)[int(len(top_datasets)*top_ds_cutoff)]
    
    # Compute final most characteristic features fraction (mol_freq_cutoff) of features for the subset of datasets
    mol_freq2 = pd.Series((adat.X[top_datasets>=co, :] > 0).sum(axis=0), index=adat.var.index) / adat.X[top_datasets>=co, :].shape[0]
    
    com2 = mol_freq2.sort_values(ascending=False)[int(len(mol_freq2)*mol_freq2_cutoff)]
    
    return mol_freq2[mol_freq2>=com2], top_datasets[top_datasets>=co], mol_freq[mol_freq>=com]


def update_minmax_mass(adat: AnnData):
    
    if 'mass' not in adat.var.columns:
        raise KeyError(' "mass" must be in adat.var.columns')
    
    min_mass_dict = {}
    max_mass_dict = {}

    for i in range(adat.X.shape[0]):
        ds_id = adat.obs.index[i]
        dataslice = adat.X[i, :] > 0

        min_mass_dict[ds_id] = adat.var['mass'][dataslice].min()
        max_mass_dict[ds_id] = adat.var['mass'][dataslice].max()

    adat.obs['minmass'] = pd.Series(min_mass_dict)
    adat.obs['maxmass'] = pd.Series(max_mass_dict)

def coloc_preprocessing(adata, scaling=True):
    if scaling:
        conv = (adata.X/adata.X.sum(axis=0))*1e4
    else:
        conv = adata.X
    
    tmp = conv.transpose()
    tmp = tmp.reshape((tmp.shape[0], adata.obs['y'].max()+1, -1))
    tmp2 = ndimage.median_filter(tmp, size=(1,3,3))
    
    tmp3 = tmp2.reshape((tmp2.shape[0], -1))
    mask = tmp3 < np.percentile(tmp3, q=50, axis=1)[:, np.newaxis]
    tmp3[mask] = 0
    
    return tmp3


def coloc_preprocessing_array(arr, maxy, scaling=True):
    if scaling:
        conv = (arr/arr.sum(axis=0))*1e4
    else:
        conv = arr
    
    tmp = conv.transpose()
    tmp = tmp.reshape((tmp.shape[0], maxy, -1))
    tmp2 = ndimage.median_filter(tmp, size=(1,3,3))
    
    tmp3 = tmp2.reshape((tmp2.shape[0], -1))
    mask = tmp3 < np.percentile(tmp3, q=50, axis=1)[:, np.newaxis]
    tmp3[mask] = 0
    
    return tmp3


def compute_subplot_arrangement(num_items):
    factors = []
    if num_items%2 !=0:
        num_items = num_items+1
    for i in range(2, int(math.sqrt(num_items)) + 2):
        if num_items % i == 0:
            factors.append((i, num_items // i))

    best_score = float('-inf')
    best_arrangement = None

    for rows, cols in factors:
        aspect_ratio = cols / rows
        score = aspect_ratio * (rows + cols)  # Define your scoring function here

        if score > best_score:
            best_score = score
            best_arrangement = (rows, cols)

    return best_arrangement