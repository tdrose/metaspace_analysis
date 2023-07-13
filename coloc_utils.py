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

import metaspace
import linex2metaspace as lx2m
from scipy import stats
import networkx as nx
import time

import multiprocessing as mp

import utils
from config import store_dir, data_dir, date_key, enrichment_dir

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

        tmp_array = np.stack(list(sums.values()))
        tmp_molecules = np.array(list(sums.keys()))
        tmp_ymax = adat2.obs['y'].max()+1

        # Coloc preprocessing:
        conv_data = utils.coloc_preprocessing_array(tmp_array.transpose(), tmp_ymax)

        coloc = pairwise_kernels(conv_data, metric='cosine')

        coloc_df = pd.DataFrame(coloc, columns=tmp_molecules, index=tmp_molecules)

        coloc_dict[dsid] = coloc_df
        
    return coloc_dict, molecule_names

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

def compute_lx_nets(coloc_dict: Dict[str, pd.DataFrame], molecule_names: Dict[str, List], ref_lip_dict, class_reacs, bootstraps: int=30, tissue = None):
    lx_nets = {}
    lx_annotations = {}
    if tissue is not None:
            print(tissue)
            
    for dsid in coloc_dict.keys():
        
        tmp_annotations = pd.Series({x: molecule_names[x] for x in coloc_dict[dsid].columns})

        parsed_lipids = lx2m.parse_annotation_series(tmp_annotations, 
                                                     ref_lip_dict, 
                                                     verbose=False) # True if you want to see all lipids that were not parsed

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

    results = [pool.apply_async(compute_lx_nets, args=(colocs['coloc_dict'], colocs['molecule_names'], 
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
        coloc_dict, molecule_names = compute_colocs_shuffle(adatas)
    else:
        coloc_dict, molecule_names = compute_colocs(adatas)
    ii_dict = list_same_colocs(coloc_dict)
    c_measures = coloc_measures(ii_dict, min_datasets=min_datasets, num_datasets=len(adatas))

    if not shuffle:
        out_dict['coloc_dict'] = coloc_dict
        out_dict['molecule_names'] = molecule_names
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
