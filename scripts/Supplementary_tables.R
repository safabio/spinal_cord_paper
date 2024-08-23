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
    grepl("D05", sample) ~ "B05_1",
    grepl("D07", sample) ~ "B07_1",
    grepl("lumb_1", sample) ~ "L10_1",
    grepl("lumb_2", sample) ~ "L10_2",
    grepl("lumb_int", sample) ~ "B10_int",
    grepl("poly_1", sample) ~ "P10_1",
    grepl("poly_2", sample) ~ "P10_2",
    grepl("poly_int", sample) ~ "P10_int"
  ))

write.csv(file_df, file = "tables/Supp_table_1.csv")

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
    grepl("devel", sample) ~ "Devel",
    grepl("ctrl", sample) ~ "B10_int",
    grepl("lumb", sample) ~ "L10_int",
    grepl("poly", sample) ~ "P10_int"
  ))

write.csv(file_df, file = "tables/Supp_table_2.csv")

### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### Supp table 3 & 4 / GO Terms and Kegg Pathways #### 
### ### ### ### ### ### ### ### ### ### ### ### ### ###

files <- list.files("output/", pattern = ".rds") %>% 
  `[`(grepl("\\_module_all_.*_080824.rds$", .))

table_list <- list()

for (i in seq_along(files)) {
  
  x <- readRDS(paste0("output/",files[i]))
  
  name <- str_remove(files[i], "\\_module_all_.*_080824.rds$")
  
  table_list[[i]] <- bind_rows(x, .id = "module") %>% 
    rownames_to_column("ID") %>% 
    mutate(ID = str_remove(ID, pattern = "\\.+\\d+$")) %>% 
    separate(
      module,
      into = c("module_number", "module"),
      sep = "_",
      fill = "left") %>% 
    mutate(type = case_when(
      grepl("^GO", ID) ~ "GOTerms",
      grepl("^path", ID) ~ "KEGGpath"
    )) %>% 
    mutate(sample = name) %>% 
    mutate(sample = case_when(
      sample == "Gg_devel_int"~ "Devel",
      sample == "Gg_ctrl_int" ~ "B10_int",
      sample == "Gg_lumb_int" ~ "L10_int",
      sample == "Gg_poly_int" ~ "P10_int"
    ))
  
  
}

Go_terms <- bind_rows(table_list[c(1,3,5,7)])

write.csv(Go_terms, "tables/Supp_table_3.csv")

Kegg_path <- bind_rows(table_list[c(2,4,6,8)])

write.csv(Kegg_path, "tables/Supp_table_4.csv")






