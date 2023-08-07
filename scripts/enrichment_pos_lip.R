library(bmetenrichr)
library(tidyverse)
library("rjson")
"%nin%" = Negate("%in%")


analysis_config <- fromJSON(file="../analysis_config.json")

store_path = analysis_config$store_dir
enrichment_path = analysis_config$enrichment_dir


metabo_exact_masses = readRDS("../enrichment_utils/all_exact_masses.rds")
metaspace_databases = readRDS("../enrichment_utils/all_metaspace_databases.rds")
# metabo_subclasses_univ = readRDS("enrichment_utils/Metabo_subclass_by_metabo_name.rds")
# metabo_tissue_loc = readRDS("enrichment_utils/Metabo_tissue_location_by_metabo_name.rds")
# metabo_pathway_univ = readRDS("enrichment_utils/Metabo_reduced_pathways_by_metabo_name.rds")
# metabo_classes_univ = readRDS("enrichment_utils/Metabo_class_by_metabo_name.rds")
lion_terms = readRDS("../enrichment_utils/LION.rds")

core_metab = read.delim("../enrichment_utils/core_metabolome_v3.csv")
bg_sf = read.csv(file.path(enrichment_path, "new_bg_pos.csv"))[,2] %>% unlist() %>% unique()
bg_name = metaspace_databases$name[which(metaspace_databases$formula %in% 
                                         bg_sf & metaspace_databases$db == "HMDB")] %>% unique()


all_files = list.files(enrichment_path)
files_for_enrichment = all_files[grep("^pos_lip", all_files)]

enrichment_results = list()
counter=1

for (enrichment_file in files_for_enrichment){
  curr_query = read.csv(file.path(enrichment_path, enrichment_file))[,2] %>% unlist() %>% unique()
  cat(paste0('File ', counter, ' of ', length(files_for_enrichment), ' | ',
             sum(curr_query %in% bg_sf), 
             ' of ', length(curr_query), 
             ' query molecules are in the universe.'))
  
  result_simple = bmetenrichr:::Run_bootstrap_ORA(marker_list=curr_query, 
                                                  background=lion_terms, 
                                                  custom_universe = bg_sf,
                                                  alpha_cutoff = 0.2, min_intersection = 2,
                                                  consider_isobars = T, polarization_mode = NA, 
                                                  mass_range_ppm = 3,
                                                  only_HMDB = T, annot_weights = NULL,
                                                  n_bootstraps = 50, sample_core_metab_only = F,
                                                  boot_fract_cutoff = 0.5, q.val_cutoff = 0.2,
                                                  selected_terms = NULL)
  
  enrichment_results[[enrichment_file]] = result_simple[[1]]
  counter = counter + 1
  
  outfilename = file.path(enrichment_path, paste0('enrichment_result_', enrichment_file))
  
  write_csv(result_simple[[1]]$clean_enrich_res, outfilename)
}

write_csv(bmetenrichr:::LION_LUT, file.path(enrichment_path, 'LION_ontology_translator.csv'))
