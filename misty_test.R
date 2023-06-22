library(mistyR)
library(future)
# data manipulation
library(dplyr)
library(purrr)
library(distances)
# plotting
library(ggplot2)
plan(multisession, workers=3)
#plan(sequential)


data("synthetic")

result.folders <- synthetic %>% imap_chr(function(sample, name) {
    
  sample.expr <- sample %>% select(-c(row, col, type, starts_with("lig")))
  sample.pos <- sample %>% select(row, col)
  
  create_initial_view(sample.expr) %>% add_paraview(sample.pos, l = 10) %>%
    run_misty(results.folder = paste0("/g/alexandr/tim/misty_test3", .Platform$file.sep, name))
})
