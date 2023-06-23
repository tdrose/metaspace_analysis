library(mistyR)
library(future)
# data manipulation
library(dplyr)
library(stringr)
library(purrr)
library(distances)
# plotting
library(ggplot2)


plan(multisession, workers=5)
#plan(sequential)

data_folder = '/scratch/trose/pos_lip_formisty/'
results_folder = '/scratch/trose/pos_lip_misty_results'


files = list.files(data_folder)

for (file in files){
    print(file)
    df = read.csv(paste0(data_folder, file))
    name = str_split(file, '\\.')[[1]][1]
    
    sample.expr <- df %>% select(-c(row, col, X))
    sample.pos <- df %>% select(row, col)

    create_initial_view(sample.expr) %>% add_paraview(sample.pos, l = 250) %>%
    run_misty(results.folder = paste0(results_folder, .Platform$file.sep, name), cached = FALSE, append = FALSE)
}
