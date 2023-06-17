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
        tmp_ymax = adat.obs['y'].max()+1

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

    for ii, x in ii_dict.items():
        if len(x) >= min_datasets:
            if not all(np.array(x)==0):
                mean_l.append(np.mean(x))
                median_l.append(np.median(x))
                var_l.append(np.var(x))
                cv_l.append(np.std(x)/np.mean(x))
                cooc_l.append(len(x)/num_datasets)
                ion_pairs.append(ii)
    
    return pd.DataFrame({'mean': mean_l, 
                         'variance': var_l, 
                         'cv': cv_l, 
                         'coocurrence': cooc_l, 
                         'ion_pairs': ion_pairs, 
                         'median': median_l}).set_index('ion_pairs', drop=False)

def compute_lx_nets(coloc_dict: Dict[str, pd.DataFrame], molecule_names: Dict[str, List], ref_lip_dict, class_reacs, bootstraps: int=30):
    lx_nets = {}
    lx_annotations = {}

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
        
    return lx_nets, lx_annotations

def tissue_lx_nets(tissue_colocs: Dict[str, pd.DataFrame], ref_lip_dict, class_reacs, bootstraps: int=30):
    
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
    
def all_tissue_colocs(tissue_adatas: Dict[str, Dict[str, AnnData]], min_dataset_fraction: int=0.3, shuffle: bool=False):
    
    out_dict = {}
    for tissue, adatas in tissue_adatas.items():
        print(tissue)
        out_dict[tissue] = {}
        
        if shuffle:
            coloc_dict, molecule_names = compute_colocs_shuffle(adatas)
        else:
            coloc_dict, molecule_names = compute_colocs(adatas)
        ii_dict = list_same_colocs(coloc_dict)
        
        tmp = int(len(adatas)*min_dataset_fraction)
        
        min_datasets = tmp if tmp>2 else 2
        
        c_measures = coloc_measures(ii_dict, min_datasets=min_datasets, num_datasets=len(adatas))
        
        if not shuffle:
            out_dict[tissue]['coloc_dict'] = coloc_dict
            out_dict[tissue]['molecule_names'] = molecule_names
            out_dict[tissue]['ii_dict'] = ii_dict
            out_dict[tissue]['c_measures_min_datasets'] = min_datasets
        
        out_dict[tissue]['c_measures'] = c_measures
    return out_dict

def flatten(l):
    return [item for subl in l for item in subl]


def coloc_pvalues(original_colocs, shuffled_colocs_list, metric: str = 'mean', ):
    
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
                
        original_colocs[tissue]['c_measures']['pvalue'] = pd.Series(pval_dict) * len(pval_dict)
