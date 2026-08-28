#' ---
#' title: "03. Comparison of parameter values across treatments"
#' author: "Mathilde Chen & Jacques Avelino, CIRAD"
#' date: "September 2026"
#' ---

# CLEAR ENVIRONMENT-------------------------------------------------------------

rm(list = ls())

# > Path of the project
path_project <- "D:/Mes Donnees/CIRAD/COLLABORATIONS/2026_DPCs/"


# LOAD NECESSARY PACKAGES-------------------------------------------------------

# > General data management 
library(tidyverse) 


# LOAD DATA --------------------------------------------------------------------

# > File containing the results of models fitting
load(paste0(path_project, "/00_DATA/RESULTS/01-1_FITS.rda"))

# > File containing the results of models fitting
load(paste0(path_project, "/00_DATA/RESULTS/01-2_FITTED_CURVES.rda"))

# EXTRACT MODELS PARAMETERS AND FPCA SCORES ------------------------------------

# > Extract parameters of growth models as estimated in the epifitter package
params_lin <- all_fits$LIN %>% 
  map_dfr(., ~{
    
    .x$Parameters %>%
      dplyr::select(.id, model, y0, r) 
    
  }, .id = "dataset") %>%
  split(.$dataset) %>%
  map_dfr(., ~{
    
    .x %>% 
      separate(col = ".id", c("year_i", "treat_j", "rep_k"), remove = F, sep = "_") %>% 
      dplyr::select(.id, model, year_i, treat_j, rep_k, y0, r) 
    
  }, .id = "dataset") 

params_nlin <- all_fits$NLIN %>% 
  map_dfr(., ~{
    
    .x$Parameters %>%
      dplyr::select(.id, model, y0, r, K) 
    
  }, .id = "dataset") %>%
  split(.$dataset) %>%
  map_dfr(., ~{
    
    .x %>% 
      separate(col = ".id", c("year_i", "treat_j", "rep_k"), remove = F, sep = "_") %>% 
      dplyr::select(.id, model, year_i, treat_j, rep_k, y0, r, K) 
    
  }, .id = "dataset") %>%
  mutate(model = paste0("nlin_", model))


# > Extract scores derived from FPCA  
params_fpca <- res_fpca_fit %>% 
  map_dfr(., ~{
    
    .x$Scores 
    
  }, .id = "dataset") %>%
  mutate(model = "FPCA") %>%
  split(.$dataset) %>%
  map_dfr(., ~{
    
    .x %>% 
      separate(col = ".id", c("year_i", "treat_j", "rep_k"), remove = F, sep = "_") %>% 
      dplyr::select(.id, model, year_i, treat_j, rep_k, starts_with("fpca_score"))
    
  }, .id = "dataset")

# > Merge both set of parameters 
params <- rbind(epi_fit_lin_params %>% 
                  gather(key = "param", value = "value", y0, r),
                epi_fit_nlin_params %>% 
                  gather(key = "param", value = "value", y0, r, K),
                fpca_fit_params %>% 
                  gather(key = "param", value = "value", starts_with("fpca_score"))) %>%
  mutate(type = case_when(
    model %in% c("Exponential", "Monomolecular", "Logistic", "Gompertz") ~ "Linear approach",
    model %in% c("nlin_Monomolecular", "nlin_Logistic", "nlin_Gompertz") ~ "Nonlinear approach",
    model == "FPCA" ~ "FPCA approach")) %>% 
  mutate(type = factor(type, levels = c("Linear approach","Nonlinear approach","FPCA approach"))) %>%
  mutate(model = recode(model,
                        "Exponential"="EXP", 
                        "Monomolecular"="MNM", 
                        "Logistic"="LOG", 
                        "Gompertz"="GOM",
                        "nlin_Monomolecular"="MNM", 
                        "nlin_Logistic"="LOG", 
                        "nlin_Gompertz"="GOM")) %>% 
  mutate(model = factor(model, levels = c("EXP", "LOG", "MNM", "GOM", "FPCA")))%>%
  mutate(param = factor(param, levels = c("y0", "r", "K", paste("fpca_score", 1:10, sep = ""))))
