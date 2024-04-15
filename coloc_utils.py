import os
from tqdm import tqdm
import pickle
import random
from typing import List, Dict, Tuple
from anndata import AnnData

import pandas as pd
import numpy as np
import math
import scipy
from sklearn.metrics.pairwise import pairwise_kernels
import statsmodels
import scanpy as sc

import metaspace
import linex2metaspace as lx2m
from scipy import stats
import networkx as nx
import time
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns

import multiprocessing as mp

import utils
from config import store_dir, data_dir, date_key, enrichment_dir

# Plotting parameters
XXSMALL_SIZE = 5
XSMALL_SIZE = 6
SMALLSMALL_SIZE = 8
SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 18
cm = 1/2.54

def image_colorbar(ax, fig, im, size="5%", pad=0.1, shrink=0.6, location='bottom'):
    div = make_axes_locatable(ax)
    cax = div.append_axes(location, size=size, pad=pad)
    cbar = fig.colorbar(im, cax=cax, shrink=0.6, location=location)

def min_pixels(df, min_pixels = 20):
    ds_bool = df.groupby('ds').agg('count')['x'] > min_pixels
    
    dsl = list(ds_bool[ds_bool].index)
    
    return df[df['ds'].isin(dsl)]

def load_tissue_datasets(ds_list: pd.Series, include_offsample: bool=False, fdr: float=0.1, path=os.path.join(store_dir, 'all_ionimages/')) -> Dict[str, AnnData]:
    out_ads = {}
    for ds in ds_list:
        out_ads[ds] = pickle.load(open(os.path.join(path, f'{ds}.pickle'), "rb" ))
        if not include_offsample:
            out_ads[ds] = out_ads[ds][:, out_ads[ds].var['offSample'] == False]
        out_ads[ds] = out_ads[ds][:, out_ads[ds].var['fdr'] <= fdr]
        if out_ads[ds].shape[1] < 2:
            del out_ads[ds]
    return out_ads


def load_alltissue_datasets(tissue_ds_dict: Dict[str, pd.Series], include_offsample: bool=False, 
                            fdr: float=0.1, path=os.path.join(store_dir, 'all_ionimages/')) -> Dict[str, Dict[str, AnnData]]:
    out_dict = {}
    for k, v in tissue_ds_dict.items():
        out_dict[k] = load_tissue_datasets(v.index, include_offsample=include_offsample, fdr=fdr, path=path)
    
    return out_dict


def compute_colocs(ads: Dict[str, AnnData]) -> Tuple[Dict[str, pd.DataFrame], Dict[str, List]]:
    
    molecule_names = {}
    molecule_ids = {}
    coloc_dict = {}

    for dsid, adat in ads.items():

        # Remove isobars
        if 'has_isobar' in adat.var.columns:
            adat2 = adat[:, ~adat.var['has_isobar']]
        else:
            adat2 = adat
            
        unique_labels = np.unique(adat2.var.formula)
        sums = {}

        # Iterate over the unique labels
        for label in unique_labels:
            
            # Get the indices of rows with the current label
            indices = np.where(adat2.var.formula == label)[0]
            # Sum up the corresponding rows and store the result
            if len(indices)>1:
                sums[label] = np.sum(adat2.X[:, indices], axis=1)
            else:
                sums[label] = adat2.X[:, indices[0]]

            molecule_names[label] = adat2.var[adat2.var['formula']==label]['moleculeNames'][0]
            molecule_ids[label] = adat2.var[adat2.var['formula']==label]['moleculeIds'][0]

        tmp_array = np.stack(list(sums.values()))
        tmp_molecules = np.array(list(sums.keys()))
        tmp_ymax = adat2.obs['y'].max()+1

        # Coloc preprocessing:
        conv_data = utils.coloc_preprocessing_array(tmp_array.transpose(), tmp_ymax)

        coloc = pairwise_kernels(conv_data, metric='cosine')

        coloc_df = pd.DataFrame(coloc, columns=tmp_molecules, index=tmp_molecules)

        coloc_dict[dsid] = coloc_df
        
    return coloc_dict, molecule_names, molecule_ids


def compute_colocs_shuffle(ads: Dict[str, AnnData]) -> Tuple[Dict[str, pd.DataFrame], Dict[str, List]]:
    
    molecule_names = {}
    coloc_dict = {}

    for dsid, adat in ads.items():

        unique_labels = np.unique(adat.var.formula)
        sums = {}

        # Iterate over the unique labels
        for label in unique_labels:
            # Get the indices of rows with the current label
            indices = np.where(adat.var.formula == label)[0]
            # Sum up the corresponding rows and store the result
            if len(indices)>1:
                sums[label] = np.sum(adat.X[:, indices], axis=1)
            else:
                sums[label] = adat.X[:, indices[0]]

            molecule_names[label] = adat.var[adat.var['formula']==label]['moleculeNames'][0]

        tmp_array = np.stack(list(sums.values()))
        tmp_molecules = np.array(list(sums.keys()))
        np.random.shuffle(tmp_molecules)
        tmp_ymax = adat.obs['y'].max()+1

        # Coloc preprocessing:
        conv_data = utils.coloc_preprocessing_array(tmp_array.transpose(), tmp_ymax)

        coloc = pairwise_kernels(conv_data, metric='cosine')

        coloc_df = pd.DataFrame(coloc, columns=tmp_molecules, index=tmp_molecules)

        coloc_dict[dsid] = coloc_df
        
    return coloc_dict, molecule_names


def list_same_colocs(coloc_dict: Dict[str, pd.DataFrame]) -> Dict[Tuple[str, str], List[float]]:
    ii_dict = {}
    for ds in coloc_dict.values():
        idx = list(ds.index)
        for ion1 in range(len(idx)):
            for ion2 in range(ion1, len(idx)):
                if idx[ion1] != idx[ion2]:
                    tmp = tuple(sorted([idx[ion1], idx[ion2]]))
                    if tmp in ii_dict.keys():
                        ii_dict[tmp].append(ds.loc[idx[ion1], idx[ion2]])
                    else:
                        ii_dict[tmp] = [ds.loc[idx[ion1], idx[ion2]]]
                        
    return ii_dict


def coloc_measures(ii_dict: Dict[Tuple[str, str], List[float]],
                   num_datasets: int,
                   min_datasets: int=10,
                  ):
    mean_l = []
    median_l = []
    var_l = []
    cv_l = []
    cooc_l = []
    ion_pairs = []
    iqr_l = []
    mediqr_l = []

    for ii, x in ii_dict.items():
        if len(x) >= min_datasets:
            if not all(np.array(x)==0):
                mean_l.append(np.mean(x))
                var_l.append(np.var(x))
                cv_l.append(np.std(x)/np.mean(x))
                cooc_l.append(len(x)/num_datasets)
                ion_pairs.append(ii)
                
                med = np.median(x)
                iqr = stats.iqr(x)
                median_l.append(med)
                iqr_l.append(iqr)
                mediqr_l.append(med/iqr)
                
    
    return pd.DataFrame({'mean': mean_l, 
                         'variance': var_l, 
                         'cv': cv_l, 
                         'coocurrence': cooc_l, 
                         'ion_pairs': ion_pairs, 
                         'median': median_l,
                         'iqr': iqr_l,
                         'mediqr': mediqr_l
                        }).set_index('ion_pairs', drop=False)


def compute_lx_nets(coloc_dict: Dict[str, pd.DataFrame], molecule_names: Dict[str, List], molecule_ids, ref_lip_dict, class_reacs, bootstraps: int=30, tissue = None):
    lx_nets = {}
    lx_annotations = {}
    if tissue is not None:
            print(tissue)
            
    for dsid in coloc_dict.keys():
        
        tmp_annotations = pd.Series({x: molecule_names[x] for x in coloc_dict[dsid].columns})
        tmp_ids = pd.Series({x: molecule_ids[x] for x in coloc_dict[dsid].columns})

        parsed_lipids = lx2m.parse_annotation_series(tmp_annotations, 
                                                     ref_lip_dict, 
                                                     verbose=False,
                                                     db_ids=tmp_ids
                                                    ) # True if you want to see all lipids that were not parsed

        keep_annotations = lx2m.annotations_parsed_lipids(parsed_lipids)
        parsed_annotations = pd.DataFrame({'molecule_names': tmp_annotations, 'parsed_lipids': parsed_lipids})
        parsed_annotations = parsed_annotations.loc[keep_annotations,:]


        net = lx2m.bootstrap_networks(lx2m.unique_sum_species(parsed_annotations['parsed_lipids']), 
                                      parsed_annotations['parsed_lipids'], 
                                      n=bootstraps, 
                                      lx2_class_reacs=class_reacs, 
                                      lx2_reference_lipids=lx2m.get_lx2_ref_lips(), 
                                      return_composed=True, print_iterations=False)


        ion_net = lx2m.ion_weight_graph(net, 
                                        lx2m.unique_sum_species(parsed_annotations['parsed_lipids']), 
                                        bootstraps=bootstraps, 
                                        parsed_lipids=parsed_annotations['parsed_lipids'])

        lx_nets[dsid] = ion_net
        lx_annotations[dsid] = parsed_annotations
    
    if tissue is None:
        return lx_nets, lx_annotations
    else:
        print(f'{tissue} done')
        return lx_nets, lx_annotations, tissue


def tissue_lx_nets_mp(tissue_colocs: Dict[str, pd.DataFrame], ref_lip_dict, class_reacs, bootstraps: int=30, threads=3):
    
    pool = mp.Pool(threads)

    results = [pool.apply_async(compute_lx_nets, args=(colocs['coloc_dict'], colocs['molecule_names'], colocs['molecule_ids'],
                                                       ref_lip_dict, class_reacs, bootstraps, tissue)) for tissue, colocs in tissue_colocs.items()]
    
    # Step 3: Don't forget to close
    pool.close() 
    
    out_dict = {}
    for x in results:
        nets, annotations, tissue = x.get()
        out_dict[tissue] = {}
        out_dict[tissue]['nets'] = nets
        out_dict[tissue]['parsed_annotations'] = annotations
        
    return out_dict

    out_dict = {}
    for tissue, colocs in tissue_colocs.items():
        print(tissue)
        out_dict[tissue] = {}
        nets, annotations = compute_lx_nets(colocs['coloc_dict'], molecule_names=colocs['molecule_names'], 
                                            ref_lip_dict=ref_lip_dict, class_reacs=class_reacs, bootstraps=bootstraps)
        out_dict[tissue]['nets'] = nets
        out_dict[tissue]['parsed_annotations'] = annotations
    
    return out_dict


def catch_sp(g, source, target):
    try:
        return nx.shortest_path_length(g, source, target)
    except nx.NetworkXNoPath:
        return np.inf

    
def coloc_worker(items, min_dataset_fraction, shuffle):
    
    tissue, adatas = items
    print(tissue)
        
    tmp = int(len(adatas)*min_dataset_fraction)
    min_datasets = tmp if tmp>2 else 2
    
    out_dict = {}
    
    if shuffle:
        print('Warning! Calling legacy shuffling')
        coloc_dict, molecule_names = compute_colocs_shuffle(adatas)
    else:
        coloc_dict, molecule_names, molecule_ids = compute_colocs(adatas)
    ii_dict = list_same_colocs(coloc_dict)
    c_measures = coloc_measures(ii_dict, min_datasets=min_datasets, num_datasets=len(adatas))

    if not shuffle:
        out_dict['coloc_dict'] = coloc_dict
        out_dict['molecule_names'] = molecule_names
        out_dict['molecule_ids'] = molecule_ids
        out_dict['ii_dict'] = ii_dict
        out_dict['c_measures_min_datasets'] = min_datasets

    out_dict['c_measures'] = c_measures
    
    print(f'{tissue} finished')
    
    return (tissue, out_dict)


def all_tissue_colocs_mp(tissue_adatas: Dict[str, Dict[str, AnnData]], min_dataset_fraction: int=0.3, shuffle: bool=False, threads: int=4):
    
    pool = mp.Pool(threads)

    # Step 2: `pool.apply` the `howmany_within_range()`
    results = [pool.apply_async(coloc_worker, args=(it, min_dataset_fraction, shuffle)) for it in tissue_adatas.items()]

    # Step 3: Don't forget to close
    pool.close()    
    
    out = {}
    for x in results:
        tmp = x.get()
        out[tmp[0]] = tmp[1]
        
    return out


def flatten(l):
    return [item for subl in l for item in subl]


def coloc_pvalues(original_colocs, shuffled_colocs_list, metric: str = 'mean', var: str = 'pvalue'):
    
    for tissue in original_colocs.keys():
        pval_dict = {}
        # get set of all co-occurrences
        coloc_levels = set([x for x in original_colocs[tissue]['c_measures']['coocurrence'].values])
        for cl in coloc_levels:
            # select only pairs of the same number of colocs
            oc = original_colocs[tissue]['c_measures'][original_colocs[tissue]['c_measures']['coocurrence'] == cl]
            scl = []
            for sc_item in shuffled_colocs_list:
                scl.append(sc_item[tissue]['c_measures'][sc_item[tissue]['c_measures']['coocurrence'] == cl])
            
            sc = pd.concat(scl)
            
            # Loop over every original value:
            for ion_pair, value in oc[metric].items():
                # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC379178/
                pval_dict[ion_pair] = (sum(sc[metric] >= value)+1) / (len(sc[metric])+1)
        
        original_colocs[tissue]['c_measures'][f'{var}'] = pd.Series(pval_dict)
        original_colocs[tissue]['c_measures'][f'{var}_corr'] = pd.Series(pval_dict) * len(pval_dict)
        

def coloc_pval_worker(tissue, original_colocs, shuffled_colocs_list, metric):
    
    print(f'{tissue} pval')
    
    pval_dict = {}
    # get set of all co-occurrences
    coloc_levels = set([x for x in original_colocs[tissue]['c_measures']['coocurrence'].values])
    for cl in coloc_levels:
        # select only pairs of the same number of colocs
        oc = original_colocs[tissue]['c_measures'][original_colocs[tissue]['c_measures']['coocurrence'] == cl]
        scl = []
        for sc_item in shuffled_colocs_list:
            scl.append(sc_item[tissue]['c_measures'][sc_item[tissue]['c_measures']['coocurrence'] == cl])

        sc = pd.concat(scl)

        # Loop over every original value:
        for ion_pair, value in oc[metric].items():
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC379178/
            pval_dict[ion_pair] = (sum(sc[metric] >= value)+1) / (len(sc[metric])+1)
            
    print(f'{tissue} pval finished')
            
    return (tissue, pd.Series(pval_dict), pd.Series(pval_dict) * len(pval_dict))
    
        
def coloc_pvalues_mp(original_colocs, shuffled_colocs_list, metric: str = 'mean', threads: int=4):
    
    pool = mp.Pool(threads)

    # Step 2: `pool.apply` the `howmany_within_range()`
    results = [pool.apply_async(coloc_pval_worker, args=(tissue, original_colocs, shuffled_colocs_list, metric)) for tissue in original_colocs.keys()]

    # Step 3: Don't forget to close
    pool.close()    
    
    out = {}
    for x in results:
        (tis, pval, pval_corr) = x.get()
        
        original_colocs[tis]['c_measures'][f'pval_{metric}'] = pval
        original_colocs[tis]['c_measures'][f'pval_{metric}_corr'] = pval_corr

        
def shuffle_symmetric_array(arr):
    n = arr.shape[0]  # Number of rows/columns

    # Step 1: Extract the lower triangle
    lower_triangle = arr[np.tril_indices(n, k=-1)]

    # Step 2: Shuffle the values in the lower triangle
    np.random.shuffle(lower_triangle)

    # Step 3: Assign the shuffled values back to the lower triangle
    arr2 = np.ones(arr.shape)
    
    arr2[np.tril_indices(n, k=-1)] = lower_triangle

    # Copy the lower triangle to the upper triangle to preserve symmetry
    arr2.T[np.tril_indices(n, k=-1)] = lower_triangle

    return arr2


def coloc_sampling_worker(items, min_dataset_fraction, save_coloc_dict):
    
    tissue, measures = items
    print(tissue)
        
    tmp = int(len(measures['coloc_dict'])*min_dataset_fraction)
    min_datasets = tmp if tmp>2 else 2
    
    
    
    # shuffling
    # 1. append all lower triangle values to a long array
    all_coloc_array = np.array([])
    
    for k, v in measures['coloc_dict'].items():
        curr_array = np.array(v)
        curr_lower_triangle = curr_array[np.tril_indices(curr_array.shape[0], k=-1)]
        all_coloc_array = np.append(all_coloc_array, curr_lower_triangle)
    
    # 2. shuffle list of all colocs
    np.random.seed((os.getpid() * int(time.time())) % 123456789)
    np.random.shuffle(all_coloc_array)
    
    # 3. place them back in original coloc matrix
    coloc_dict = {}
    counter = 0
    for k, v in measures['coloc_dict'].items():
        # make new array
        curr_array = np.array(v)
        
        upper = curr_array[np.tril_indices(curr_array.shape[0], k=-1)].shape[0]
        lower_triangle = all_coloc_array[counter:(counter+upper)]
        counter += upper
        
        arr2 = np.ones(curr_array.shape)
        arr2[np.tril_indices(curr_array.shape[0], k=-1)] = lower_triangle.copy()

        # Copy the lower triangle to the upper triangle to preserve symmetry
        arr2.T[np.tril_indices(curr_array.shape[0], k=-1)] = lower_triangle.copy()
        
        coloc_dict[k] = pd.DataFrame(arr2.copy(), columns=v.columns, index=v.index)
    
    ii_dict = list_same_colocs(coloc_dict)
    c_measures = coloc_measures(ii_dict, min_datasets=min_datasets, num_datasets=len(measures['coloc_dict']))
    
    out_dict = {}
    out_dict['c_measures'] = c_measures
    if save_coloc_dict:
        out_dict['coloc_dict'] = coloc_dict
    
    print(f'{tissue} finished')
    
    return (tissue, out_dict)


def all_tissue_coloc_sampling_mp(tissue_colocs, min_dataset_fraction=0.3, threads: int=4, save_coloc_dict=False):
    
    pool = mp.Pool(threads)

    # Step 2: `pool.apply` the `howmany_within_range()`
    results = [pool.apply_async(coloc_sampling_worker, args=(it, min_dataset_fraction, save_coloc_dict)) for it in tissue_colocs.items()]

    # Step 3: Don't forget to close
    pool.close()    
    
    out = {}
    for x in results:
        tmp = x.get()
        out[tmp[0]] = tmp[1]
        
    return out


def mark_isobars(tissue_adatas, ppm_threshold=3):
    for tis in tissue_adatas.keys():
        for dsid, ds in tissue_adatas[tis].items():
            has_isobar = pd.Series(False, index=ds.var.index)

            tmp = ds.var['mz'].values
            for i in range(len(tmp)):
                for j in range(len(tmp)):
                    if i != j:
                        if abs((tmp[i]-tmp[j])/tmp[i])*1e6 < ppm_threshold:
                            has_isobar[i] = True
                            has_isobar[j] = True
            tissue_adatas[tis][dsid].var['has_isobar'] = has_isobar
    return tissue_adatas


def remove_isobars_resultstab(results_tabs):
    
    for dsid, tab in results_tabs.items():
        results_tabs[dsid] = tab[~tab['has_isobar']]
        
    return results_tabs


def mark_isobars_resultstab(results_tabs, ppm_threshold=3, remove=True):
    
    for dsid, tab in tqdm(results_tabs.items()):
        has_isobar = pd.Series(False, index=tab.index)

        for idx, val in enumerate(tab['mz'].values):
            tmp = abs((val-tab['mz'])/val)*1e6 < ppm_threshold
            tmp[idx] = False
            has_isobar = has_isobar | tmp
            
        results_tabs[dsid]['has_isobar'] = has_isobar
    
    if remove:
        return remove_isobars_resultstab(results_tabs)
    
    return results_tabs




def dict_val_min(in_dict):
    mean_values = [(key, sum(values) / len(values)) for key, values in in_dict.items()]

    # Sort the list of tuples based on the mean value (second element of the tuple)
    sorted_keys_by_mean = sorted(mean_values, key=lambda x: x[1])

    # Extract only the keys in the sorted order
    sorted_keys = [key for key, _ in sorted_keys_by_mean]
    
    return sorted_keys[0]


def avg_lipid_ranking(formulas, tissue_networks, alternative_molecule_names, alternative_ids):
    
    formula_ranking = {}
    if tissue_networks is not None:
        # Loop over all networks
        for net in tissue_networks.values():
            for formula, values in dict(net.nodes(data=True)).items():
                if formula not in formula_ranking.keys():
                    formula_ranking[formula] = {'ranking': {}, 'candidate':{}}
                for count, lipid in enumerate(values['sum_species'].index):
                    if lipid not in formula_ranking[formula]['ranking'].keys():
                        formula_ranking[formula]['ranking'][lipid] = [count]
                        for candidate in values['parsed_lipids']:
                            if candidate.sum_species_str() == lipid:
                                formula_ranking[formula]['candidate'][lipid] = candidate.get_ids()['db_id']
                    else:
                        formula_ranking[formula]['ranking'][lipid].append(count)
    
    # Loop over all formulas and get the best ranked molecule
    out_list = []
    out_ids = []
    for idx, formula in enumerate(formulas):
        if formula not in formula_ranking.keys():
            out_list.append(alternative_molecule_names[idx])
            out_ids.append(alternative_ids[idx][0]) # Just pick first one
        else:
            best_lipid = dict_val_min(formula_ranking[formula]['ranking'])
            out_list.append(best_lipid)
            out_ids.append(formula_ranking[formula]['candidate'][best_lipid])
                           
    return out_list, out_ids


def tissue_modules(coloc_df, tissue, lipid_networks, coloc_object, summary_metric='mean', scaling=False):
    tissue_df = coloc_df[coloc_df['tissue'] == tissue]
    # If only working with significant colocs:
    #tissue_df2 = tissue_df[tissue_df['mediqr_sig']]

    selected_ions = list(set([item for tuple in tissue_df.index for item in tuple]))
    coloc_matrix = pd.DataFrame(np.zeros(len(selected_ions)*len(selected_ions)).reshape((len(selected_ions), -1)), columns=selected_ions, index=selected_ions)
    for k1 in coloc_matrix.index:
        for k2 in coloc_matrix.index:
            k = tuple(sorted((k1, k2)))
            if k in tissue_df.index:
                coloc_matrix.loc[k[0], k[1]] = tissue_df[summary_metric][k]
                coloc_matrix.loc[k[1], k[0]] = tissue_df[summary_metric][k]
    # Simplify naming
    new_idx = [coloc_object[tissue]['molecule_names'][idx][0][0:30] for idx in coloc_matrix.index]
    coloc_matrix.index = new_idx
    coloc_matrix.columns = new_idx
    if lipid_networks is None:
        tis_net = None
    else:
        tis_net = lipid_networks[tissue]['nets']
    better_molecule_names, better_molecule_ids = avg_lipid_ranking(formulas=selected_ions, tissue_networks=tis_net, 
                                                                   alternative_molecule_names=new_idx,
                                                                   alternative_ids = [coloc_object[tissue]['molecule_ids'][x] for x in selected_ions]
                                                                  )
    
    
    tmp_df = pd.DataFrame({'formula': selected_ions, 
                           'molecule': better_molecule_names, 
                           'all_hmdb': [coloc_object[tissue]['molecule_ids'][x] for x in selected_ions],
                           'selected_hmdb': better_molecule_ids
                          }
                         )
    
    ad = AnnData(X=np.array(coloc_matrix), 
                 obs = tmp_df,
                 var = tmp_df,
                 dtype=np.float32)
    sc.pp.neighbors(ad)
    sc.tl.umap(ad)
    sc.tl.leiden(ad)
                 
    
    ax = sc.pl.umap(ad, color='leiden', show=False)
    
    for line in range(0,ad.shape[0]):
        ax.text(ad.obsm['X_umap'][line, 0]+np.random.normal(0, 0.05), 
                ad.obsm['X_umap'][line, 1]+np.random.normal(0, 0.05), 
                ad.var['molecule'].values[line], 
                horizontalalignment='left', size='small', color='black')
    ax.set_title(tissue)
    plt.show()
    return ad


def all_modules(coloc_df, lipid_networks, coloc_object, summary_metric='mean', scaling=False):
    # If only working with significant colocs:
    #tissue_df2 = tissue_df[tissue_df['mediqr_sig']]

    selected_ions = list(set([item for tuple in coloc_df.index for item in tuple]))
    coloc_matrix = pd.DataFrame(np.zeros(len(selected_ions)*len(selected_ions)).reshape((len(selected_ions), -1)), columns=selected_ions, index=selected_ions)
    for k1 in coloc_matrix.index:
        for k2 in coloc_matrix.index:
            k = tuple(sorted((k1, k2)))
            if k in coloc_df.index:
                coloc_matrix.loc[k[0], k[1]] = coloc_df[summary_metric][k]
                coloc_matrix.loc[k[1], k[0]] = coloc_df[summary_metric][k]
    
    # Simplify naming
    combined_molecule_names = {}
    for tis in coloc_object.keys():
        for form, names in coloc_object[tis]['molecule_names'].items():
            combined_molecule_names[form] = names
            
    combined_molecule_ids = {}
    for tis in coloc_object.keys():
        for form, names in coloc_object[tis]['molecule_ids'].items():
            combined_molecule_ids[form] = names
    
    new_idx = [combined_molecule_names[idx][0][0:30] for idx in coloc_matrix.index]
    
    coloc_matrix.index = new_idx
    coloc_matrix.columns = new_idx
    
    if lipid_networks is None:
        all_tissue_nets = None
    else:
        tissue_net_list = [lipid_networks[tis]['nets'] for tis in lipid_networks.keys()]
        all_tissue_nets = {}
        for d in tissue_net_list:
            all_tissue_nets.update(d)
    
    better_molecule_names, better_molecule_ids = avg_lipid_ranking(formulas=selected_ions, tissue_networks=all_tissue_nets, 
                                                                   alternative_molecule_names=new_idx,
                                                                   alternative_ids = [combined_molecule_ids[x] for x in selected_ions]
                                                                  )
    
    
    tmp_df = pd.DataFrame({'formula': selected_ions, 
                           'molecule': better_molecule_names, 
                           'all_hmdb': [combined_molecule_ids[x] for x in selected_ions],
                           'selected_hmdb': better_molecule_ids
                          }
                         )
    
    ad = AnnData(X=np.array(coloc_matrix), 
                 obs = tmp_df,
                 var = tmp_df,
                 dtype=np.float32)
    sc.pp.neighbors(ad)
    sc.tl.umap(ad)
    sc.tl.leiden(ad)
                 
    
    ax = sc.pl.umap(ad, color='leiden', show=False)
    
    for line in range(0,ad.shape[0]):
        ax.text(ad.obsm['X_umap'][line, 0]+np.random.normal(0, 0.05), 
                ad.obsm['X_umap'][line, 1]+np.random.normal(0, 0.05), 
                ad.var['molecule'].values[line], 
                horizontalalignment='left', size='small', color='black')
    plt.show()
    return ad


def write_cluster_sets(adata, path='unnamed'):
    for i in set(adata.obs['leiden']):
        adata.obs[adata.obs['leiden']==i][['formula']].to_csv(path+i+'.csv')

        
def get_ad_molecule_matrices(adata):
    # Molecule df summed over all adducts for the same formula
    molecule_df = pd.DataFrame(adata.X.transpose()).assign(formula=adata.var['formula'].reset_index(drop=True)).groupby('formula').sum()
    molecule_list = list(molecule_df.index)
    molecule_matrix = molecule_df.to_numpy().reshape((molecule_df.shape[0], adata.obs['y'].max()+1, -1))
    # print(molecule_matrix.shape)
    return {'molecule_list': molecule_list, 'molecule_images': molecule_matrix}

def molecule_adata(adata, mdt):
    molecule_df = pd.DataFrame(adata.X.transpose()).assign(formula=adata.var['formula'].reset_index(drop=True)).groupby('formula').sum().reset_index()
    ad= AnnData(
        X=np.array(molecule_df.drop(columns=['formula'])).transpose(), 
       obs=adata.obs.assign(
           organism=mdt.loc[adata.uns['metaspace_id'], :]['Organism'],
           organ=mdt.loc[adata.uns['metaspace_id'], :]['Organism_Part'],
           ds=adata.uns['metaspace_id']),
       var=molecule_df[['formula']]
                  )
    sc.pp.filter_genes(ad, min_cells=100)
    return ad[~((ad.X==0).all(axis=1)), :]


def tissue_mol_mat(adata_dict):
    return {key: get_ad_molecule_matrices(val) for key, val in adata_dict.items()}


def get_cluster_images(molecule_dict, cluster_assignment, q=50):
    out_dict = {}
    for cluster in set(cluster_assignment):
        molecules_in_cluster = cluster_assignment[cluster_assignment==cluster].index
        cluster_mask = [x in molecules_in_cluster for x in molecule_dict['molecule_list']]
        
        number_of_molecules = sum(cluster_mask)
        if number_of_molecules > 0:
            
            # Using percentile
            selected_images = molecule_dict['molecule_images'][cluster_mask]
            
            # Hotspot clipping
            hotspots = np.percentile(selected_images, q=99, axis=(1,2))
            for i in range(len(hotspots)):
                selected_images[i][selected_images[i] > hotspots[i]] = hotspots[i]
            
            percentile_image = np.percentile(selected_images, q=q, axis=0)
            
            centroid_image = selected_images[np.argmax(np.square((selected_images-percentile_image)).sum(axis=(1,2)))]
            
            out_dict[int(cluster)] = {'number_of_molecules': number_of_molecules, 
                                      'mean_ion_image': percentile_image}
        else:
            out_dict[int(cluster)] = {'number_of_molecules': number_of_molecules, 
                                      'mean_ion_image': np.array([])}
            
    return out_dict


def display_cluster_ion_images(cluster_assignment, adatas, adata_selection, q=50, fig=None, transpose=None):
    tmm = tissue_mol_mat(adatas)
    ds_list = []
    for dsid in adata_selection:
        ds_list.append(get_cluster_images(tmm[dsid], cluster_assignment, q=q))
    
    n_clusters = max([max(v.keys()) for v in ds_list]) + 1
    if fig is None:
        fig, axs = plt.subplots(len(adata_selection), n_clusters)
    else:
        axs = fig.subplots(len(adata_selection), n_clusters)
    for ds in range(len(ds_list)):
        for cl in range(n_clusters):
            axs[ds][cl].axis('off')
            if ds_list[ds][cl]['number_of_molecules'] > 0:
                if transpose is None:
                    axs[ds][cl].imshow(ds_list[ds][cl]['mean_ion_image'], cmap='viridis')
                else:
                    if transpose[ds]:
                        axs[ds][cl].imshow(ds_list[ds][cl]['mean_ion_image'].transpose(), cmap='viridis')
                    else:
                        axs[ds][cl].imshow(ds_list[ds][cl]['mean_ion_image'], cmap='viridis')
                nm = ds_list[ds][cl]['number_of_molecules']
            if ds == 0:
                axs[ds][cl].set_title(f'Cluster {cl}')
                
                #axs[ds][cl].set_title(f'Number of molecules: {nm}')
    plt.show()

    
def get_lipidclass_assignment(lipid_class_list, molmat, ds_lipidnetwork):
    assignment_list = []
    node_data = dict(ds_lipidnetwork.nodes(data=True))
    for mol in molmat['molecule_list']:
        if mol in node_data.keys():
            curr_class = lx2.lipid_parser(node_data[mol]['sum_species'].index[0], reference_lipids=ref_lip_dict).get_lipid_class()
            if curr_class in lipid_class_list:
                assignment_list.append(lipid_class_list.index(curr_class))
            else:
                assignment_list.append(-100)
        else:
            # We are only interested in binary comparisons here this will be ignored
            assignment_list.append(-100)
        
    return pd.Series(assignment_list, index=molmat['molecule_list'])


def display_lipidclass_ion_images(lipid_class_list, 
                                  adatas_tissue0, adatas_tissue1, 
                                  selection_tissue0, selection_tissue1,
                                  lx_tissue0, lx_tissue1
                                 ):
    tmm0 = tissue_mol_mat(adatas_tissue0)
    tmm1 = tissue_mol_mat(adatas_tissue1)
    ds_list0 = [get_cluster_images(tmm0[dsid], 
                                   get_lipidclass_assignment(lipid_class_list, 
                                                             tmm0[dsid], 
                                                             lx_tissue0[dsid])) for dsid in selection_tissue0]
    ds_list1 = [get_cluster_images(tmm1[dsid], 
                                   get_lipidclass_assignment(lipid_class_list, 
                                                             tmm1[dsid], 
                                                             lx_tissue1[dsid])) for dsid in selection_tissue1]
    
    
    n_clusters = max([max(v.keys()) for v in ds_list0] + [max(v.keys()) for v in ds_list1]) + 1
    fig, axs = plt.subplots(len(selection_tissue0), n_clusters)
    for ds in range(len(ds_list0)):
        for cl in range(n_clusters):
            axs[ds][cl].axis('off')
            if cl in ds_list0[ds].keys():
                if ds_list0[ds][cl]['number_of_molecules'] > 0:
                    axs[ds][cl].imshow(ds_list0[ds][cl]['mean_ion_image'])
                    lc = lipid_class_list[cl]
                    axs[ds][cl].set_title(f'Lipid class: {lc}')
    plt.show()
    
    fig, axs = plt.subplots(len(selection_tissue1), n_clusters)
    for ds in range(len(ds_list1)):
        for cl in range(n_clusters):
            axs[ds][cl].axis('off')
            if cl in ds_list1[ds].keys():
                if ds_list1[ds][cl]['number_of_molecules'] > 0:
                    axs[ds][cl].imshow(ds_list1[ds][cl]['mean_ion_image'])
                    lc = lipid_class_list[cl]
                    axs[ds][cl].set_title(f'Lipid class: {lc}')
    plt.show()

    
def load_enrichment(template_file, tissue, ontologies, lion_terms=None):
    # Get all files
    all_files = os.listdir(enrichment_dir)
    
    cluster_dict = {}
    # match files based on cluster
    for file in all_files:
        
        # If the file belongs to one of the ontologies
        matched_temp = [file.startswith(template_file.format(x, tissue)) for x in ontologies]
        if any(matched_temp):
            matched_ont =np.array(ontologies)[matched_temp][0]
            clusterid = int(file.split(tissue)[1].split('.')[0])
            
            if clusterid not in cluster_dict.keys():
                cluster_dict[clusterid] = []
            
            # Read table
            tab = pd.read_csv(os.path.join(enrichment_dir, file)).assign(ontology=matched_ont)
            
            if 'LION' in file:
                tab = tab.set_index('term').join(lion_terms.set_index('ID'), how='left').reset_index(drop=True).rename(columns={'name': 'term'})
            
            cluster_dict[clusterid].append(tab)
    
    # merge tables
    return {key: pd.concat(val) for key, val in cluster_dict.items()}

def plot_enrichment(tab, cluster, max_elem, ax=None, axisfont=5, textfont=4, color='ontology'):
    tab2 = tab.sort_values('ES_median')
    tab2['-log10(q-value)'] = -np.log10(tab2['q.value_median'])
    if ax is None:
        fig, ax = plt.subplots(1)
    if color in tab2.columns:
        sns.barplot(data=tab2, x='ES_median', y='term', hue=color, dodge=False, ax=ax, width=.5, saturation=1)
    else:
        sns.barplot(data=tab2, x='ES_median', y='term', color=color, dodge=False, ax=ax, width=.5, saturation=1)
    #sns.scatterplot(data=tab2, y='ES_median', x='term', hue='-log10(q-value)', size='-log10(q-value)', ax=ax, sizes=(40, 300))
    ax.set_ylim((-1, max_elem))
    ax.set_xlim((0, tab2['ES_median'].max()*2.5))
    ax.set_title(f'Cluster: {cluster}')
    
    sns.despine(offset=5, trim=False, ax=ax)
    ticks = ax.get_yticklabels()

    # for each tick write the value at the right position    
    for tick in ticks:
        term = tick.get_text()
        yval = tick.get_position()[1]
        # x position
        xval = tab2.set_index('term')['ES_median'][term]
        ax.text(xval, 
                yval, 
                term, 
                horizontalalignment='left', size=textfont, color='black', va='center')
    
    ax.set_yticklabels([])
    ax.get_yaxis().set_visible(False)
    
    
    
    

lipid_class_colors = {
    "ACoA": "#ffffff",
    "ACar": "#ffffff",
    "NAE": "#ffffff",
    "NEFA": "#ffffff",
    "MG": "#f1c232",
    "MGO": "#d0a728",
    "MGMG": "#ffffff",
    "MGMGO": "#ffffff",
    "DGMG": "#ffffff",
    "SQMG": "#ffffff",
    "DG": "#bf9000",
    "DGO": "#9d7703",
    "MGDG": "#ffffff",
    "MGDGO": "#ffffff",
    "DGDG": "#ffffff",
    "TGDG": "#ffffff",
    "SQDG": "#ffffff",
    "SQDGO": "#ffffff",
    "TG": "#7f6000",
    "TGO": "#634b01",
    "BMP": "#ffffff",
    "CL": "#e35913",
    "MLCL": "#e57e4a",
    "DLCL": "#e8a888",
    "CDPDAG": "#b8a69a",
    "LCDPDAG": "#ffffff",
    "PE": "#3c78d8",
    "PEO": "#1c4588",
    "PEP": "#1c4589",
    "LPE": "#a4c2f4",
    "LPEO": "#4f739d",
    "MMPE": "#ffffff",
    "LMMPE": "#ffffff",
    "DMPE": "#ffffff",
    "LDMPE": "#ffffff",
    "NAPE": "#ffffff",
    "NAPEO": "#ffffff",
    "NALPE": "#ffffff",
    "NALPEO": "#ffffff",
    "PC": "#6aa84f",
    "PCO": "#38761d",
    "LPC": "#b6d7a8",
    "LPCO": "#93c47d",
    "PA": "#00ffff",
    "PAO": "#009d9d",
    "LPA": "#cbffff",
    "LPAO": "#a9e5e5",
    "DGPP": "#ffffff",
    "PG": "#ff0000",
    "LPG": "#e06666",
    "LysylPG": "#ffffff",
    "PGP": "#ffffff",
    "LPGP": "#ffffff",
    "PI": "#b215da",
    "LPI": "#c470da",
    "PIM": "#ffffff",
    "LPIM": "#ffffff",
    "PIMIP": "#ffffff",
    "LPIMIP": "#ffffff",
    "PIN": "#ffffff",
    "LPIN": "#ffffff",
    "PIP": "#660d7d",
    "PIP2": "#4c0a5d",
    "PIP3": "#280530",
    "PS": "#c0e510",
    "LPS": "#dced8c",
    "Dol-": "#ffffff",
    "DolP-": "#ffffff",
    "DolPHex-": "#ffffff",
    "GD1": "#ffffff",
    "GD2": "#ffffff",
    "GD3": "#ffffff",
    "GD4": "#ffffff",
    "GM1": "#ffffff",
    "GM2": "#ffffff",
    "GM3": "#ffffff",
    "GM4": "#ffffff",
    "GT1": "#ffffff",
    "GT2": "#ffffff",
    "GT3": "#ffffff",
    "GA1": "#ffffff",
    "GQ3/2/1": "#ffffff",
    "GP3/2/1": "#ffffff",
    "LCB": "#f18ea5",
    "LCBP": "#c9909d",
    "Cer": "#f24a71",
    "CerP": "#ffffff",
    "EPC": "#ffffff",
    "HexCer": "#b03350",
    "LHexCer": "#ffffff",
    "SHexCer": "#ffffff",
    "LSHexCer": "#ffffff",
    "Hex2Cer": "#671e2f",
    "Hex3Cer": "#41121d",
    "IPC": "#ffffff",
    "LIPC": "#ffffff",
    "MIPC": "#ffffff",
    "M(IP)2C": "#ffffff",
    "SM": "#f24aa4",
    "LSM": "#d9d9d9",
    "ACer": "#d9d9d9",
    "ST": "#d9d9d9",
    "Chol": "#ffff00",
    "SE": "#d9d9d9",
    "CE": "#cfcf00"
}
