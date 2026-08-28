#' ---
#' title: "00. Format data for analyses"
#' authors: "Mathilde Chen & Jacques Avelino"
#' date: "September 2026"
#' ---


# CLEAR ENVIRONMENT-------------------------------------------------------------

rm(list = ls())

# > Path to data
path_project <- "D:/Mes Donnees/CIRAD/COLLABORATIONS/2026_DPCs/"

# LOAD NECESSARY PACKAGES-------------------------------------------------------

# > Open / load xlsx datasets
library(openxlsx) ; library(readxl)
library(janitor)

# > General data management 
library(tidyverse) 


# 1. NATURAL EPIDEMIC PROGRESS -------------------------------------------------
# Data on black leaf streak disease on plantain in Dominican Republic

# > Load data
load(paste0(path_project, "00_DATA/01_BANANA_BLSD_NATURAL-EPIDEMIC/data_fragment_for_analyses.rda"))

# > Full data
full_data_bsld <- data_fragment_for_analyses$ydata %>% 
  left_join(data_fragment_for_analyses$X_init, by = c(".id"="ID")) %>% 
  group_by(fragment_unique) %>% 
  mutate(year_i  = first(year(date1)),
         treat_j = compartiment,
         rep_k   = fragment_unique) %>% 
  ungroup() %>% 
  unite(".id", c(year_i, treat_j, rep_k), sep="_", remove = FALSE) %>% 
  dplyr::select(.id, year_i, treat_j, rep_k, 
                "plant"="bananier", 
                "leaf"="feuille",
                fragment,
                "day_since_leaf_rank_4"=".index", 
                "inc"=.value,
                nblesions1, date1)

# > Data for analyses 
data_for_temporal_analysis_blsd <- data_fragment_for_analyses$ydata %>% 
  left_join(data_fragment_for_analyses$X_init, 
            by = c(".id"="ID")) %>% 
  group_by(fragment_unique) %>% 
  mutate(year_i  = first(year(date1)),
         treat_j = compartiment,
         rep_k   = fragment_unique) %>% 
  ungroup() %>% 
  unite(".id", c(year_i, treat_j, rep_k), sep="_", remove = FALSE) %>%
  mutate(.id = as.character(.id)) %>% 
  # > compute % of incidence
  mutate(max_surface = as.numeric(as.character(max(.value, na.rm=T)))) %>% 
  mutate(.value = as.numeric(as.character(.value))/max_surface) %>% 
  rename(".inc"=".value") %>%
  # > replace 0 values by very small values 
  mutate(.inc = if_else(.inc == 0, 10^-6, .inc)) %>% 
  dplyr::select(.id, .index, .inc) ; str(data_for_temporal_analysis_blsd)

#tibble [2,730 × 3] (S3: tbl_df/tbl/data.frame)
#$ .id   : int [1:2730] 1 1 1 1 1 1 1 1 1 1 ...
#$ .index: num [1:2730] 3 7 14 17 21 28 31 35 38 42 ...
#$ .inc  : num [1:2730] 0.000001 0.011969 0.171679 0.177186 0.721312 ...

# > Check the number of curves
full_data_bsld %>% 
  distinct(plant, leaf, fragment) %>% 
  nrow() # 192 curves


# 2. GENOTYPE EFFECT -----------------------------------------------------------
# Effect of cacao genotype on monilia assessed in trials in Costa Rica 

# > Load data 
# one line corresponds to the 1 pod caracterized by its initiation date, diseased date (if any), and/or harvest date 
data_cacao_CR_raw <- read.xlsx(xlsxFile = paste0(path_project, "00_DATA/02_CACAO_MONILIA_GENOTYPE/BD_Moniliose.xlsx"), 
                               sheet = 2, detectDates = T) %>% 
  rename("Emergence_date"="Inicio", "Monilia_date"="Monilia", "Sporulation_date"="Esporas", "Wilt_date"="Cherele", 
         "Other_date"="Otro", "Removal_date"="Removido", "Harvest_date"="Cosechado", "Final_state"="Estado.final") 

# Identify range of possible dates accross th dataset 
min_obs <- min(data_cacao_CR_raw$Emergence_date)
max_obs <- max(data_cacao_CR_raw$Emergence_date)

# Count the number of emerged and diseased pods 
# for any of the date 
full_data_cacao_genotype <- expand_grid(Any_date = as.Date(min_obs:max_obs), 
                                        ID = unique(data_cacao_CR_raw$ID)) %>% 
  # > Add the week and year of observation 
  mutate(Any_date_week  = week(Any_date), 
         Any_date_year  = year(Any_date),
         Any_date_order = as.numeric(Any_date - min(Any_date))) %>% 
  # > Add the dates of emergence 
  left_join(., data_cacao_CR_raw %>% 
              dplyr::select(ID, Clon, Repeticion, Arbol, Emergence_date, Monilia_date), by = "ID", 
            relationship = "many-to-many") %>% 
  # > Indicate whether the pod is already emerged and diseased at date t
  mutate(Presence_at_date_t = if_else(Emergence_date <= Any_date, 1, 0),
         Monilia_at_date_t  = if_else(Monilia_date   <= Any_date, 1, 0)) %>% 
  # > Count the cumulated numbers of emerged and diseased 
  group_by(Clon, Repeticion, 
           #Arbol, 
           Any_date_order, Any_date, Any_date_week, Any_date_year) %>% 
  summarise(Cum_N_pods     = sum(Presence_at_date_t),
            Cum_N_diseased = sum(Monilia_at_date_t, na.rm = T)) %>% 
  # > Final number of pods
  group_by(Clon, 
           #Arbol, 
           Repeticion) %>%
  mutate(Max_Cum_N_pods = max(Cum_N_pods)) %>%
  # > Year and treatment
  mutate(year_i = first(year(Any_date)),
         treat_j = Clon) %>% 
  ungroup() %>% 
  mutate(cum_inc = Cum_N_diseased / Max_Cum_N_pods) %>% 
  # > Create the id variable
  mutate(Repeticion = paste0("Rep", Repeticion)) %>%
  unite(col = ".id", c("year_i", "Clon", "Repeticion"), sep = "_", remove = F)
  

# Format for temporal analysis 
data_for_temporal_analysis_cacao_genotype <- full_data_cacao_genotype %>% 
  dplyr::select(.id, ".index"="Any_date_order", ".inc" = "cum_inc")


# 3. INTERVENTION ON DISEASE ---------------------------------------------------
# Trials on treatment of initial inoculum (inoc) in Honduras (1994-1996)

# > Files containing plot data 
full_data_inoc <- NULL

# > Read and concatenate data for all years for each trial
for(year_i in 1994:1996)
{
  
  data_inoc_i <- NULL
  
  # > List all files 
  list_files_inoc_i <- list.files(path = paste0(path_project, "00_DATA/03_COFFEE_CLR_INOC/Essai inoculum initial - ", year_i), 
                                  full.names = T)
  
  # > Concatenate data in each file by year
  for(my_file in list_files_inoc_i)
  {
    
    # > Treatment (T1 to T7) and repetition (1T to 4T)
    treat_j <- substr(my_file, nchar(my_file)-8, nchar(my_file)-7)
    rep_k   <- substr(my_file, nchar(my_file)-6, nchar(my_file)-5)
    
    # > Read the file 
    full_data_ijk <- read.xlsx(xlsxFile = my_file, 
                               sheet = 1)
    
    # > Select the colomns with appropriate data 
    data_ijk <- data.frame(year_i            = year_i,
                           treat_j           = treat_j, 
                           rep_k             = rep_k, 
                           day_t             = as.numeric(as.character(full_data_ijk$X2)),
                           int_days          = as.numeric(as.character(full_data_ijk$X1)),
                           cum_total_ijkt    = as.numeric(as.character(full_data_ijk$X29)),
                           cum_diseased_ijkt = as.numeric(as.character(full_data_ijk$X30)),
                           inc_perc_ijkt     = as.numeric(as.character(full_data_ijk$X36))) %>%
      drop_na(day_t)
    
    # > Add to previous data 
    data_inoc_i <- rbind(data_inoc_i, data_ijk)
    
    rm(treat_j, rep_k, data_ijk,full_data_ijk)
    
  }
  
  full_data_inoc <- rbind(full_data_inoc, data_inoc_i)
  
  rm(data_inoc_i, my_file)
  
}

# > Format for analysis 
data_for_temporal_analysis_inoc <- full_data_inoc %>% 
  unite(col = ".id", c(year_i, treat_j, rep_k), remove = F, sep = "_") %>% 
  mutate(int_days = if_else(is.na(int_days), -1, int_days)) %>% 
  filter(int_days != 0) %>% 
  dplyr::select(".id", ".index"=day_t, ".inc"=inc_perc_ijkt) %>%
  mutate(.inc = .inc/100) %>%
  mutate(.inc = if_else(.inc == 0, 10^-6, .inc))


# 4. INTERVENTIONS ON DISEASE and HOST -----------------------------------------
# Trials fertilisation x fungicides in Honduras (1994-1996)

# > Files containing plot data 
full_data_ferti <- NULL

# > Read and concatenate data for all years for each trial
for(year_i in 1994:1996)
{
  
  data_ferti_i <- NULL
  
  # > List all files 
  list_files_ferti_i <- list.files(path = paste0(path_project, "00_DATA/04_COFFEE_CLR_FERTI-FUN/Fertilisation ", year_i), 
                                   full.names = T)
  
  # > Concatenate data in each file by year
  for(my_file in list_files_ferti_i)
  {
    
    # > Treatment (T1 to T7) and repetition (1T to 4T)
    treat_j <- substr(my_file, nchar(my_file)-8, nchar(my_file)-7)
    rep_k   <- substr(my_file, nchar(my_file)-6, nchar(my_file)-5)
    
    # > Read the file 
    full_data_ijk <- read.xlsx(xlsxFile = my_file, 
                               sheet = 1)
    
    # > Select the colomns with appropriate data 
    data_ijk <- data.frame(year_i            = year_i,
                           treat_j           = treat_j, 
                           rep_k             = rep_k, 
                           day_t             = as.numeric(as.character(full_data_ijk$X2)),
                           int_days          = as.numeric(as.character(full_data_ijk$X1)),
                           cum_total_ijkt    = as.numeric(as.character(full_data_ijk$X29)),
                           cum_diseased_ijkt = as.numeric(as.character(full_data_ijk$X30)),
                           inc_perc_ijkt     = as.numeric(as.character(full_data_ijk$X36))) %>%
      drop_na(day_t)
    
    # > Add to previous data 
    data_ferti_i <- rbind(data_ferti_i, data_ijk)
    
    rm(treat_j, rep_k, data_ijk,full_data_ijk)
    
  }
  
  full_data_ferti <- rbind(full_data_ferti, data_ferti_i)
  
  rm(data_ferti_i, my_file)
  
}

# > Format for analysis 
data_for_temporal_analysis_ferti <- full_data_ferti %>% 
  unite(col = ".id", c(year_i, treat_j, rep_k), remove = F, sep = "_") %>% 
  mutate(int_days = if_else(is.na(int_days), -1, int_days)) %>% 
  filter(int_days != 0) %>% 
  dplyr::select(".id", ".index"=day_t, ".inc"=inc_perc_ijkt) %>%
  mutate(.inc = .inc/100) %>%
  mutate(.inc = if_else(.inc == 0, 10^-6, .inc))


# 5. INTERACTION GENOTYPE x ENVIRONMENT ----------------------------------------
# Data on cacao, two clones in three shading systems, trials in Columbia (2019-2021)

# > Data
# Epidemiological, production, wilt data 
load(file = paste0(path_project, "/00_DATA/05_CACAO_MULTIPEST_GENOTYPE-ENV/2.Data_per_tree.rda"))
data_2 <- data_per_date_per_tree

full_data_cacao_genotype_shade <- data_2 %>% 
  dplyr::select(-sum_monilia, -sum_phytophtora, -sum_colletotrichum, -sum_lost_wilt,  -sum_lost) %>% 
  # correct the labels for plants 
  mutate(plant = if_else(clon == "CCN-51", as.numeric(plant)-5, as.numeric(plant)),
         plant = factor(plant)) %>% 
  # Ordering the factors
  mutate(treatment = recode(treatment, 
                            "ABARCO-CACAO"="Cacao+Mahogany", "CAUCHO-CACAO"="Cacao+Rubber tree", "MATARRATON-CACAO"="Cacao")) %>% 
  mutate(shade = factor(treatment, 
                        levels = c("Cacao+Mahogany", "Cacao+Rubber tree", "Cacao")),
         block = factor(block, 
                        levels = rev(c("3", "2", "1")),
                        labels = c("Rep 3", "Rep 2", "Rep 1")),
         clon = factor(clon, levels = c("CCN-51", "ICS-95"))) %>% 
  # nest the data per plant (no consideration of the crop cycle since the observations are continuously collected)
  group_by(block, shade, clon, plant)  %>% 
  # get the total nb of fruits from time t-1
  mutate(sum_tot_t_1    = lag(sum_tot, default = 0)) %>% 
  # compute the new fruits produced between both observations
  # total nb of fruits at time t from which we retrieve the total nb of fruits at time t-1 and those lost between both observations
  mutate(new_produced = sum_tot - (sum_tot_t_1-sum_lost_harvested),
         new_produced = if_else(new_produced < 0, 0 , new_produced),
         new_produced_t_1 = lag(new_produced, default = 0)) %>% 
  # compute the cumulative numbers of pods and diseased pods
  mutate(cumsum_np        = cumsum(new_produced),
         cumsum_diseased  = cumsum(sum_diseased)) %>% 
  # compute the dynamic cumulative % of diseased pods
  mutate(dyn_perc_cumsum_diseased  = if_else(cumsum_np==0, 0, 100*(cumsum_diseased/cumsum_np))) %>%
  # compute the total number of fruits produced by each tree
  mutate(max_cumsum_np = max(cumsum_np)) %>% 
  # compute the cumulative incidence of diseased 
  mutate(cum_inc  = cumsum_diseased/max_cumsum_np) %>% 
  group_by(clon, shade, block) %>% 
  mutate(year_i = first(year(min_date))) %>% 
  # create the .id variable
  unite(".id", c(year_i, clon, shade, block), sep = "_", remove = F) %>%
  unite("treat_j", c(clon, treatment), sep = "_", remove = F) 

# Format for temporal analysis 
data_for_temporal_analysis_genotype_shade <- full_data_cacao_genotype_shade %>% 
  ungroup() %>% 
  dplyr::select(.id, ".index"="quinzaine", ".inc" = "cum_inc")


# RESULTS AGGREGATION ----------------------------------------------------------

# > Create a list where all data will be saved

list_data_for_temporal_analysis <- list()

# > Add each dataset (full and restricted for analyses) in the complete list
list_data_for_temporal_analysis[[paste0("01_BANANA_BLSD_NAT-EPIDEMIC")]] <- list(full = full_data_bsld, 
                                                                                 ready_for_analysis = data_for_temporal_analysis_blsd)

list_data_for_temporal_analysis[[paste0("02_CACAO_MON_GENOTYPE")]] <- list(full = full_data_cacao_genotype, 
                                                                           ready_for_analysis = data_for_temporal_analysis_cacao_genotype)

list_data_for_temporal_analysis[[paste0("03_COFFEE_CLR_INOC")]] <- list(full = full_data_inoc, 
                                                                        ready_for_analysis = data_for_temporal_analysis_inoc)

list_data_for_temporal_analysis[[paste0("04_COFFEE_CLR_FERTI-FUN")]] <- list(full = full_data_ferti, 
                                                                             ready_for_analysis = data_for_temporal_analysis_ferti)

list_data_for_temporal_analysis[[paste0("05_CACAO_MULTI_GENOTYPE-ENV")]] <- list(full = full_data_cacao_genotype_shade,
                                                                                 ready_for_analysis = data_for_temporal_analysis_genotype_shade)

# > check the number of curves per dataset

list_data_for_temporal_analysis %>% 
  map_dfr(., ~{
    
    .x$ready_for_analysis %>% 
      mutate(.id = as.character(.id))
    
    
  }, .id="dataset") %>% 
  separate(dataset, c("id_dataset", "crop", "disease", "id_trial"), sep = "_", remove = FALSE) %>% 
  distinct(dataset, .id) %>% 
  group_by(dataset) %>% 
  count()

#  dataset                         n
#1 01_BANANA_BLSD_NAT-EPIDEMIC   192
#2 02_CACAO_MON_GENOTYPE          12
#3 03_COFFEE_CLR_INOC             84
#4 04_COFFEE_CLR_FERTI-FUN       108
#5 05_CACAO_MULTI_GENOTYPE-ENV    18


# SAVE RESULTS -----------------------------------------------------------------

save(list_data_for_temporal_analysis, 
     file = paste0(path_project, "00_DATA/RESULTS/00_MERGED_DATA.rda"))

