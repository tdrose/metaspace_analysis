library(bmetenrichr)
library(tidyverse)
library(hypergea)
"%nin%" = Negate("%in%")

# Data --------------------------------------------------------------------
metabo_exact_masses = readRDS("enrichment_utils/all_exact_masses.rds")
metaspace_databases = readRDS("enrichment_utils/all_metaspace_databases.rds")
metabo_subclasses_univ = readRDS("enrichment_utils/Metabo_subclass_by_metabo_name.rds")
metabo_tissue_loc = readRDS("enrichment_utils/Metabo_tissue_location_by_metabo_name.rds")
metabo_pathway_univ = readRDS("enrichment_utils/Metabo_reduced_pathways_by_metabo_name.rds")
metabo_classes_univ = readRDS("enrichment_utils/Metabo_class_by_metabo_name.rds")

core_metab = read.delim("enrichment_utils/core_metabolome_v3.csv")

# Functions ----------------------------------------------------------------
get_metabo_iso = function(sf_vec, consider_isobars = T,
                          polarization_mode = NA, mass_range_ppm = 3,
                          only_HMDB = F){
  
  annotation_formulas_adduct <- gsub("\\+|\\-",".",sf_vec)
  annotation_formulas <- gsub("\\..+$","",annotation_formulas_adduct)
  annotation_adduct <- gsub("^.+\\.","",annotation_formulas_adduct)
  annotation_list <-
    lapply(annotation_formulas, function(annotation_formula_i){
      if (only_HMDB){
        res = metaspace_databases$name[which(metaspace_databases$formula == annotation_formula_i &
                                               metaspace_databases$db == "HMDB")] 
      }
      else{
        res = metaspace_databases$name[metaspace_databases$formula == annotation_formula_i]
      }
      res
    })
  is_adduct = any(c(paste0("pos", annotation_adduct), 
                    paste0("neg", annotation_adduct)) %in% colnames(metabo_exact_masses))
  if (consider_isobars & 
      is_adduct & 
      polarization_mode %in% c("negative", "positive")){
    switch(polarization_mode,
           'positive' = {
             col_name <- paste0("pos",annotation_adduct)
             
             exact_masses_slim <- metabo_exact_masses[,c(1,which(grepl("pos",colnames(metabo_exact_masses))))]
             colnames(exact_masses_slim) <- gsub("^pos","",colnames(exact_masses_slim))
             exact_masses_slim <- exact_masses_slim %>% tidyr::pivot_longer(cols = -1, values_to = "mass", names_to = "adduct")
             exact_masses_slim$formula_adduct <- paste0(exact_masses_slim$formula,".",exact_masses_slim$adduct)
           },
           'negative' = {
             col_name <- paste0("neg",annotation_adduct)
             
             exact_masses_slim <- metabo_exact_masses[,c(1,which(grepl("neg",colnames(metabo_exact_masses))))]
             colnames(exact_masses_slim) <- gsub("^neg","",colnames(exact_masses_slim))
             exact_masses_slim <- exact_masses_slim %>% tidyr::pivot_longer(cols = -1, values_to = "mass", names_to = "adduct")
             exact_masses_slim$formula_adduct <- paste0(exact_masses_slim$formula,".",exact_masses_slim$adduct)
           })
    exact_masses_annotations <-
      mapply(annotation_formulas_i = annotation_formulas,
             col_name_i = col_name,
             function(annotation_formulas_i, col_name_i) {
               mass <-
                 metabo_exact_masses[metabo_exact_masses$formula == annotation_formulas_i, col_name_i]
               if (length(mass) == 0) {
                 mass = NA
               }
               mass
               
             })
    names(exact_masses_annotations) <- annotation_formulas_adduct
    isobars_list <-
      sapply(exact_masses_annotations, function(mass_i){
        if(is.na(mass_i)){
          character(0)
        } else {
          exact_masses_slim$formula_adduct[
            between(exact_masses_slim$mass,
                    left = mass_i - (mass_range_ppm * mass_i / 1e6),
                    right = mass_i + (mass_range_ppm * mass_i / 1e6))]
        }
      }, simplify = F)
    
    isobars_list[sapply(isobars_list, length) > 1]
    
    ## remove self isobars
    isobars_list <-
      sapply(names(isobars_list), function(i){
        isobars_list[[i]][!isobars_list[[i]] %in% i]
      }, simplify = F)
    
    
  }
  else{
    isobars_list = NULL
  }
  
  if (!is.null(isobars_list)){
    annotation_list <-
      lapply(seq_along(annotation_formulas), function(i){
        c(isomer = annotation_list[[i]],
          isobar = metaspace_databases$name[metaspace_databases$formula %in% gsub("\\..+$","",isobars_list[[i]])])
        
      })
  }
  names(annotation_list) = sf_vec
  
  return(annotation_list)
}

metabo_bootstrap = function(annot_list, annot_weights = NULL, 
                            n_bootstraps = 50){
  if (!is.null(annot_weights)) {
    if (length(annot_weights) != length(annot_list)) {
      stop("annotations and annotation.weights do not have the same length")
    } else {
      ## extra check
      if (!all(sapply(annot_weights, length) == sapply(annot_list, length))) {
        stop("annotations and annotation.weights do not have the same length ")
      }
    }
  }
  bootstrapped_sublist <- pbapply::pblapply(seq(n_bootstraps),function(n_i){
    sapply(seq(length(annot_list)),function(row_number_i){
      molecules_to_sample <- annot_list[[row_number_i]]
      if(length(molecules_to_sample)==0){
          molecules_to_sample <- NA
        }
        
      if(!is.null(annot_weights)){
        weights_to_sample <- annot_weights[[row_number_i]]
      } else {
        weights_to_sample <- rep(x = 1, times = length(molecules_to_sample))
      }
      if(length(molecules_to_sample)!=1){
        sample(x = molecules_to_sample,size = 1,prob = weights_to_sample)
      } 
      else {          
        molecules_to_sample
      }
 
      })
    })
  return(bootstrapped_sublist)
}
hyper_geom_enrich = function(query,term_list, universe = NULL,
                             core_metabo_only = F){
  
  if (!is.null(universe)){
    all_bg_mols = universe
  }
  else{
    all_bg_mols = unlist(term_list) %>% unique()
  }
  
  final_res = list()
  #pb = txtProgressBar(min =0, max = length(term_list), style = 3)
  for (i in 1:length(term_list)){
    term_name = names(term_list)[i]
    term_mols = term_list[[i]]
    term_mols = term_mols[which(term_mols %in% all_bg_mols)]
    if (core_metabo_only){
      term_mols = term_mols[which(term_mols %in% core_metab$name)]
      all_bg_mols = all_bg_mols[which(all_bg_mols %in% core_metab$name)]
    }
    if (length(term_mols) == 0){
      final_res[[i]] = NULL
      #setTxtProgressBar(pb, i)
    }
    else{
      TP = length(intersect(query, term_mols))
      FP = length(which(query %nin% term_mols))
      FN = length(which(term_mols %nin% query))
      
      TN = all_bg_mols[which(all_bg_mols %nin% term_mols)]
      TN = length(TN[which(TN %nin% query)])
      
      fisher_mat = matrix(c(TP, FP, FN, TN), nrow = 2, ncol = 2, byrow = T)
      fisher_res = hypergea::hypergeom.test(fisher_mat, alternative = "greater")
      
      term_res = data.frame(term = term_name,
                            TP = TP, FP = FP, FN = FN , TN = TN,
                            OR = fisher_res$estimate %>% as.numeric(),
                            pval = fisher_res$p.value)
      final_res[[i]] = term_res
      #setTxtProgressBar(pb, i)
    }
  }
  final_res = final_res %>% dplyr::bind_rows()
  final_res$OR[which(final_res$TP == 0)] = NA
  final_res$pval[which(final_res$TP == 0)] = NA
  return(final_res)
  
  
}
simplify_hypergeom_bootstrap = function(bootstrap_list,term_list,universe = NULL,
                                        core_metabo_only = F,
                                        boot_fract_cutoff = 0.5,
                                        min_annot = 3, q.val_cutoff = 0.2){
  boot_enrich_res = pbapply::pblapply(seq(length(bootstrap_list)),function(n_i){
    enrich_res = hyper_geom_enrich(query = bootstrap_list[[n_i]],
                                   term_list = term_list, universe = universe,
                                   core_metabo_only = core_metabo_only)
    enrich_res$bootstrap = n_i
    enrich_res = enrich_res %>% dplyr::filter(!is.na(OR), !is.na(pval))
    if(nrow(enrich_res) == 0){
      return(NULL)
    }
    else{
      return(enrich_res)
    }
  })
  boot_enrich_res = boot_enrich_res %>% dplyr::bind_rows()
  boot_enrich_res <- boot_enrich_res %>%
    dplyr::group_by(term) %>%
    dplyr::mutate(fraction = length(term) / length(bootstrap_list)) %>%
    dplyr::filter(fraction > boot_fract_cutoff)
  
  final_enrich_res <- boot_enrich_res %>% group_by(bootstrap) %>%
    dplyr::mutate(q.value = p.adjust(pval, method = "fdr"))  %>%
    dplyr::group_by(term) %>%
    dplyr::summarise(n = median(TP, na.rm = T),
              ES_median = median(OR, na.rm = T),
              ES_sd = sd(OR, na.rm = T),
              p.value_median = median(pval, na.rm = T),
              p.value_sd = sd(pval, na.rm = T),
              q.value_median = median(q.value, na.rm = T),
              q.value_sd = sd(q.value, na.rm = T),
              fraction.bootstrap.presence = median(fraction, na.rm = T)) %>%
    dplyr::arrange(q.value_median)
  
  final_enrich_res = final_enrich_res %>% 
    dplyr::filter(n > min_annot, q.value_median < q.val_cutoff, term != "") %>% 
    ungroup() %>% as.data.frame()
  
  return(list("unfiltered_enrich_res" = boot_enrich_res,
              "clean_enrich_res" = final_enrich_res))
  
}


# Running enrichment for Tim's data ---------------------------------------
Tim_query_sf = read.csv("significant.csv")[,2] %>% unlist() %>% unique()
Tim_bg_sf = read.csv("background.csv")[,1] %>% unlist() %>% unique()
Tim_bg_name = metaspace_databases$name[which(metaspace_databases$formula %in% 
                                               Tim_bg_sf & metaspace_databases$db == "HMDB")] %>% unique()


x = get_metabo_iso(sf_vec = Tim_query_sf, consider_isobars = F)
y = metabo_bootstrap(annot_list = x)

tissue_final = simplify_hypergeom_bootstrap(bootstrap_list = y,
                                     term_list = metabo_tissue_loc, 
                                     universe = Tim_bg_name, core_metabo_only = F,
                                     boot_fract_cutoff = 0.5, min_annot = 3, q.val_cutoff = 0.2)

class_final = simplify_hypergeom_bootstrap(bootstrap_list = y,
                                              term_list = metabo_classes_univ, 
                                              universe = Tim_bg_name, core_metabo_only = F,
                                              boot_fract_cutoff = 0.5, min_annot = 3, q.val_cutoff = 0.2)

subclass_final = simplify_hypergeom_bootstrap(bootstrap_list = y,
                                              term_list = metabo_subclasses_univ, 
                                              universe = Tim_bg_name, core_metabo_only = F,
                                              boot_fract_cutoff = 0.5, min_annot = 3, q.val_cutoff = 0.2)
pathway_final = simplify_hypergeom_bootstrap(bootstrap_list = y,
                                              term_list = metabo_pathway_univ, 
                                              universe = Tim_bg_name, core_metabo_only = F,
                                              boot_fract_cutoff = 0.5, min_annot = 3, q.val_cutoff = 0.2)

saveRDS(list("Metabo_classes" = class_final,
             "Metabo_subclasses" = subclass_final,
             "Metabo_pathways" = pathway_final), "../Misc/Tim/enrich_results.rds")







