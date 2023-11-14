library(bmetenrichr)
library(tidyverse)
library("rjson")
"%nin%" = Negate("%in%")


analysis_config <- fromJSON(file="../analysis_config.json")

store_path = analysis_config$store_dir
enrichment_path = analysis_config$enrichment_dir


metabo_exact_masses = readRDS("../enrichment_utils/all_exact_masses.rds")
metaspace_databases = readRDS("../enrichment_utils/all_metaspace_databases.rds")
core_metab = read.delim("../enrichment_utils/core_metabolome_v3.csv")

metabo_classes_univ = readRDS("../enrichment_utils/Metabo_class_by_metabo_name.rds")
metabo_subclasses_univ = readRDS("../enrichment_utils/Metabo_subclass_by_metabo_name.rds")
lion_terms = readRDS("../enrichment_utils/LION.rds")


bg_sf = read.csv(file.path(enrichment_path, "new_bg_pos.csv"))[,2] %>% unlist() %>% unique()
bg_name = metaspace_databases$name[which(metaspace_databases$formula %in% 
                                         bg_sf & metaspace_databases$db == "HMDB")] %>% unique()


all_files = list.files(enrichment_path)
files_for_enrichment = all_files[grep("^single_pixel", all_files)]

enrichment_results = list()



# Metabo class enrichment
counter=1
print('Class enrichment')
for (enrichment_file in files_for_enrichment){
  
  if (paste0('enrichment_result_METABOclass_', enrichment_file)%in%all_files){
    print('skipping metaboclass file')
  } else{
    
    curr_query = read.csv(file.path(enrichment_path, enrichment_file))[,2] %>% unlist() %>% unique()
    cat(paste0('File ', counter, ' of ', length(files_for_enrichment), ' | ',
               sum(curr_query %in% bg_sf), 
               ' of ', length(curr_query), 
               ' query molecules are in the universe.'))
    
    result_simple = bmetenrichr:::Run_bootstrap_ORA(marker_list=curr_query, 
                                                    background=metabo_classes_univ, 
                                                    custom_universe = bg_sf,
                                                    alpha_cutoff = 0.3, min_intersection = 2,
                                                    consider_isobars = T, polarization_mode = NA, 
                                                    mass_range_ppm = 3,
                                                    only_HMDB = T, annot_weights = NULL,
                                                    n_bootstraps = 50, sample_core_metab_only = F,
                                                    boot_fract_cutoff = 0.333, q.val_cutoff = 0.3,
                                                    selected_terms = NULL)
    
    enrichment_results[[enrichment_file]] = result_simple[[1]]
    counter = counter + 1
    
    outfilename = file.path(enrichment_path, paste0('enrichment_result_METABOclass_', enrichment_file))
    
    write_csv(result_simple[[1]]$clean_enrich_res, outfilename)
  }
}


# Metabo subclass enrichment
counter=1
print('Subclass enrichment')
for (enrichment_file in files_for_enrichment){
  
  if (paste0('enrichment_result_METABOsubclass_', enrichment_file)%in%all_files){
    print('skipping metabosubclass file')
  } else{
    
    curr_query = read.csv(file.path(enrichment_path, enrichment_file))[,2] %>% unlist() %>% unique()
    cat(paste0('File ', counter, ' of ', length(files_for_enrichment), ' | ',
               sum(curr_query %in% bg_sf), 
               ' of ', length(curr_query), 
               ' query molecules are in the universe.'))
    
    result_simple = bmetenrichr:::Run_bootstrap_ORA(marker_list=curr_query, 
                                                    background=metabo_subclasses_univ, 
                                                    custom_universe = bg_sf,
                                                    alpha_cutoff = 0.3, min_intersection = 2,
                                                    consider_isobars = T, polarization_mode = NA, 
                                                    mass_range_ppm = 3,
                                                    only_HMDB = T, annot_weights = NULL,
                                                    n_bootstraps = 50, sample_core_metab_only = F,
                                                    boot_fract_cutoff = 0.333, q.val_cutoff = 0.3,
                                                    selected_terms = NULL)
    
    enrichment_results[[enrichment_file]] = result_simple[[1]]
    counter = counter + 1
    
    outfilename = file.path(enrichment_path, paste0('enrichment_result_METABOsubclass_', enrichment_file))
    
    write_csv(result_simple[[1]]$clean_enrich_res, outfilename)
  }
}


# LION enrichment
counter=1
print('LION enrichment')
for (enrichment_file in files_for_enrichment){
  
  if (paste0('enrichment_result_LION_', enrichment_file)%in%all_files){
    print('skipping lion file')
  } else{
    
    
    curr_query = read.csv(file.path(enrichment_path, enrichment_file))[,2] %>% unlist() %>% unique()
    cat(paste0('File ', counter, ' of ', length(files_for_enrichment), ' | ',
               sum(curr_query %in% bg_sf), 
               ' of ', length(curr_query), 
               ' query molecules are in the universe.'))
    
    result_simple = bmetenrichr:::Run_bootstrap_ORA(marker_list=curr_query, 
                                                    background=lion_terms, 
                                                    custom_universe = bg_sf,
                                                    alpha_cutoff = 0.3, min_intersection = 2,
                                                    consider_isobars = T, polarization_mode = NA, 
                                                    mass_range_ppm = 3,
                                                    only_HMDB = T, annot_weights = NULL,
                                                    n_bootstraps = 50, sample_core_metab_only = F,
                                                    boot_fract_cutoff = 0.333, q.val_cutoff = 0.3,
                                                    selected_terms = NULL)
    
    enrichment_results[[enrichment_file]] = result_simple[[1]]
    counter = counter + 1
    
    outfilename = file.path(enrichment_path, paste0('enrichment_result_LION_', enrichment_file))
    
    write_csv(result_simple[[1]]$clean_enrich_res, outfilename)
  }
}

write_csv(bmetenrichr:::LION_LUT, file.path(enrichment_path, 'LION_ontology_translator.csv'))
