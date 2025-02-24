### ### ### ### ### ### ### ### ### ### ### ###
#### Script to prepare Supplementary tables ###
## Fabio Sacher, 14.08.2024 ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ###

library(tidyverse)

setwd("~/spinal_cord_paper/")

### ### ### ### ### ### ### ### ### ### ###
#### Supp table 1 / Cluster annotations ####
### ### ### ### ### ### ### ### ### ### ###

files <- list.files("annotations/", 
                    pattern = "_cluster_annotation.csv$") %>% 
  `[`(!grepl("all_int|ctrl_lumb|ctrl_poly|devel", .))

samples <- str_remove(files, "_cluster_annotation.csv$")

file_list <- list()

for (i in seq_along(samples)) {
  file_list[[i]] <- read.csv(paste0("annotations/",files[i]))
}

names(file_list) <- samples

file_df <- bind_rows(file_list, .id = "sample") %>%
  mutate(sample = str_remove(sample, "\\.\\d+")) %>% 
  mutate(sample = case_when(
    grepl("ctrl_1", sample) ~ "B10_1",
    grepl("ctrl_2", sample) ~ "B10_2",
    grepl("ctrl_int", sample) ~ "B10_int",
    grepl("D05_ctrl", sample) ~ "B05_1",
    grepl("D07_ctrl", sample) ~ "B07_1",
    grepl("lumb_1", sample) ~ "L10_1",
    grepl("lumb_2", sample) ~ "L10_2",
    grepl("lumb_int", sample) ~ "L10_int",
    grepl("poly_1", sample) ~ "Poly10_1",
    grepl("poly_2", sample) ~ "Poly10_2",
    grepl("poly_int", sample) ~ "Poly10_int"
  ))

write.csv(file_df, file = "tables/Supp_table_1.csv", row.names = FALSE)

rm(file_df, file_list, files, samples, i)

### ### ### ### ### ### ### ### ### ### ### ### ###
####  Supp table 2 / scWGCNA module annotations #### 
### ### ### ### ### ### ### ### ### ### ### ### ###

files <- list.files("annotations/", 
                    pattern = "scWGCNA_module_annotation.csv$")

samples <- str_remove(files, "_scWGCNA_module_annotation.csv$")

file_list <- list()

for (i in seq_along(files)) {
  file_list[[i]] <- read.csv(paste0("annotations/",files[i]))
}

names(file_list) <- samples

file_df <- bind_rows(file_list, .id = "sample") %>% 
  mutate(sample = str_remove(sample, "\\.\\d+")) %>% 
  mutate(sample = case_when(
    grepl("ctrl", sample) ~ "B10_int",
    grepl("lumb", sample) ~ "L10_int",
    grepl("poly", sample) ~ "Poly10_int",
    grepl("devel", sample) ~ "Devel_int"
  )) %>% 
  tidyr::separate_wider_delim(
    module,
    names = c("module_number", "module"),
    delim = "_", too_few = "align_start")

write.csv(file_df, file = "tables/Supp_table_2.csv", row.names = FALSE)

rm(file_df, file_list, files, samples, i)

### ### ### ### ### ### ### ### ###
#### Supp table 3 / GO Terms  #### 
### ### ### ### ### ### ### ### ###

files <- list.files("output/", pattern = "GOTerms_080824.rds$")

table_list <- list()

for (i in seq_along(files)) {
  
  x <- readRDS(paste0("output/",files[i]))
  
  name <- str_remove(files[i], "\\_module_all_GOTerms_080824.rds")
  
  table_list[[i]] <- do.call(rbind, x) %>% 
    rownames_to_column("tmp") %>% 
    tidyr::separate_wider_delim(
      tmp,
      names = c("module", "ID"),
      delim = ".", too_few = "align_start") %>% 
    tidyr::separate_wider_delim(
      module,
      names = c("module_number", "module"),
      delim = "_", too_few = "align_start") %>% 
    mutate(type = case_when(
      grepl("^GO", ID) ~ "GOTerms"
    )) %>% 
    mutate(sample = name)
}

Go_terms <- do.call(rbind, table_list) %>% 
  mutate(sample = case_when(
    grepl("ctrl", sample) ~ "B10_int",
    grepl("lumb", sample) ~ "L10_int",
    grepl("poly", sample) ~ "Poly10_int",
    grepl("devel", sample) ~ "Devel_int"
  )) %>% 
  relocate(sample)

write.csv(Go_terms, "tables/Supp_table_3.csv", row.names = FALSE)

rm(Go_terms, table_list, x, files, i, name)

### ### ### ### ### ### ### ### ### ###
#### Supp table 4 /  Kegg Pathways #### 
### ### ### ### ### ### ### ### ### ###

files <- list.files("output/", pattern = "KEGGPath_080824.rds$")

table_list <- list()

for (i in seq_along(files)) {
  
  x <- readRDS(paste0("output/",files[i]))
  
  name <- str_remove(files[i], "\\_module_all_KEGGPath_080824.rds$")
  
  table_list[[i]] <- do.call(rbind, x) %>% 
    rownames_to_column("tmp") %>% 
    tidyr::separate_wider_delim(
      tmp,
      names = c("module", "ID"),
      delim = ".", too_few = "align_start") %>% 
    tidyr::separate_wider_delim(
      module,
      names = c("module_number", "module"),
      delim = "_", too_few = "align_start") %>% 
    mutate(type = case_when(
      grepl("^path", ID) ~ "KEGGPath"
    )) %>% 
    mutate(sample = name)
  
}

Kegg_paths <- do.call(rbind, table_list) %>% 
  mutate(sample = case_when(
    grepl("ctrl", sample) ~ "B10_int",
    grepl("lumb", sample) ~ "L10_int",
    grepl("poly", sample) ~ "Poly10_int",
    grepl("devel", sample) ~ "Devel_int"
  )) %>% 
  relocate(sample)


write.csv(Kegg_paths, "tables/Supp_table_4.csv", row.names = FALSE)

rm(Kegg_paths, table_list, x, files, i, name)

### ### ### ### ### ### ### ### ### ### ### ###
#### Supp table 5 /  Volcano plot DE genes #### 
### ### ### ### ### ### ### ### ### ### ### ###

marker_list <- list()

marker_list[["L10int"]] <- readRDS("~/spinal_cord_paper/data/Gg_ctrl_lumb_int_markers.rds") %>% 
  mutate(clust_id = as.numeric(str_remove(cluster, "^cl-"))) %>% 
  filter(clust_id %in% c(16, 17, 18, 27)) %>%  
  mutate(cluster = fct_drop(cluster)) %>% 
  mutate(data = "B10int_vs_L10int") %>% 
  mutate(cell_type = case_when(
    clust_id == 16 ~ "MFOL_early",
    clust_id == 17 ~ "exc_neuron_1",
    clust_id == 18 ~ "MFOL_late",
    clust_id == 27 ~ "exc_neuron_2"
  ))

marker_list[["P10int"]] <- readRDS("~/spinal_cord_paper/data/Gg_ctrl_poly_int_markers.rds") %>% 
  mutate(clust_id = as.numeric(str_remove(cluster, "^cl-"))) %>% 
  filter(clust_id %in% c(24, 25)) %>%  
  mutate(cluster = fct_drop(cluster)) %>% 
  mutate(data = "B10int_vs_Poly10int") %>% 
  mutate(cell_type = case_when(
    clust_id == 24 ~ "MFOL",
    clust_id == 25 ~ "motor_neurons"
  ))



marker_df <- do.call(rbind, marker_list)

write.csv(marker_df, "tables/Supp_table_5.csv", row.names = FALSE)

rm(marker_list, marker_df)

### ### ### ### ### ### ### ### ### ### ### ###
#### Supp table 6 / DE genes all se objects #### 
### ### ### ### ### ### ### ### ### ### ### ###

files <- list.files("data/", 
                    pattern = "_fullDE_") %>% 
  `[`(!grepl("all_int|devel", .))

samples <- str_remove(files, "_fullDE_\\d{6}.txt$")

file_list <- list()

for (i in seq_along(samples)) {
  file_list[[i]] <- read.delim(paste0("data/",files[i]))
}

names(file_list) <- samples

file_df <- bind_rows(file_list, .id = "sample") %>%
  mutate(sample = str_remove(sample, "\\.\\d+")) %>% 
  mutate(sample = case_when(
    grepl("^Gg_ctrl_1", sample) ~ "B10_1",
    grepl("^Gg_ctrl_2", sample) ~ "B10_2",
    grepl("^Gg_ctrl_int", sample) ~ "B10_int",
    grepl("^Gg_D05_ctrl", sample) ~ "B05_1",
    grepl("^Gg_D07_ctrl", sample) ~ "B07_1",
    grepl("^Gg_lumb_1", sample) ~ "L10_1",
    grepl("^Gg_lumb_2", sample) ~ "L10_2",
    grepl("^Gg_lumb_int", sample) ~ "L10_int",
    grepl("^Gg_poly_1", sample) ~ "Poly10_1",
    grepl("^Gg_poly_2", sample) ~ "Poly10_2",
    grepl("^Gg_poly_int", sample) ~ "Poly10_int",
    grepl("^Gg_ctrl_lumb", sample) ~ "B_L10int",
    grepl("^Gg_ctrl_poly", sample) ~ "B_P10int"
  ))

write.csv(file_df, file = "tables/Supp_table_6.csv", row.names = FALSE)

rm(file_df, file_list, files, samples, i)
