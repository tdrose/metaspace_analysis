library(mistyR)
library(future)
# data manipulation
library(dplyr)
library(stringr)
library(purrr)
library(distances)
# plotting
library(ggplot2)


plan(multisession, workers=2)
#plan(sequential)

args <- commandArgs(trailingOnly = TRUE)

data_file = args[1]
results_folder = args[2]


print(data_file)
df = read.csv(data_file)
name = tail(str_split(str_split(data_file, '\\.')[[1]][1],'/')[[1]], n=1)

sample.expr <- df %>% select(-c(row, col, X))
sample.pos <- df %>% select(row, col)

create_initial_view(sample.expr) %>% add_paraview(sample.pos, l = 250) %>%
run_misty(results.folder = paste0(results_folder, .Platform$file.sep, name), cached = FALSE, append = FALSE)
