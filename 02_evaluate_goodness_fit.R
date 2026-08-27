#' ---
#' title: "02. Evaluate the goodness of fit for the different models"
#' author: "Mathilde Chen & Jacques Avelino, CIRAD"
#' date: "September 2026"
#' ---

# CLEAR ENVIRONMENT-------------------------------------------------------------

rm(list = ls())

# > Path of the project
path_project <- "D:/Mes Donnees/CIRAD/COLLABORATIONS/2026_DPCs/"

# LOADING NECESSARY PACKAGES----------------------------------------------------

# > General data management 
library(tidyverse) 

# > measures of fit quality 
library(caret)  
library(DescTools)


# LOAD DATA --------------------------------------------------------------------

# > File containing the results of models fitting
load(paste0(path_project, "/00_DATA/RESULTS/FITS.rda"))

# > File containing the results of models fitting
load(paste0(path_project, "/00_DATA/RESULTS/FITTED_CURVES.rda"))

# FIT QUALITY ACROSS ALL OBSERVATIONS ------------------------------------------
#  Six metrics are used to characterize the goodness of fit based on all estimations, 
#  regardless of the curve to which they belong

# > Table with measure of fit quality
fit_quality_tot <- all_curves_fit %>% 
  group_by(dataset, type, model) %>% 
  summarise(n = n(),
            mean_O = mean(y, na.rm = T),
            mean_P = mean(predicted, na.rm = T),
            RSS    = sum(residual^2),
            # > RMSE
            RMSE   = caret::RMSE(pred = predicted, obs = y, na.rm=T),
            # > Normalized RMSE
            RRMSE  = RMSE/mean_O, 
            # > Coefficient of determination
            R2     = caret::R2(pred = predicted, obs = y, na.rm=T),
            # > Pearson
            r_Pearson = cor(y, predicted, method = "pearson"), 
            # > Lin's CCC
            CCC = as.numeric(as.character(CCC(y, predicted)$rho[1]))) %>% 
  dplyr::select(-mean_O, -mean_P) %>% 
  arrange(dataset, type)


# FIT QUALITY PER CURVE  -------------------------------------------------------
#  Six metrics are used to characterize the goodness of fit estimated using the estimations
#  of each curve

# > Extract performance of growth models as estimated in the epifitter package
perf_LIN <- all_fits$LIN %>% 
  map_dfr(., ~{
    
    .x$Parameters %>% 
      dplyr::select(.id, best_model, model, r_squared)
    
  }, .id = "dataset")

perf_NLIN <- all_fits$NLIN %>% 
  map_dfr(., ~{
    
    .x$Parameters %>% 
      mutate(model = paste0("nlin_", model)) %>% 
      dplyr::select(.id, best_model, model, r_squared)
    
  }, .id = "dataset")

# > Compute performance of FPCA 
perf_FPCA <- all_curves_fit %>% 
  filter(model == "FPCA") %>% 
  group_by(dataset, model, .id) %>% 
  summarise(r_squared = caret::R2(pred = predicted, obs = y, na.rm=T),
            CCC = as.numeric(as.character(CCC(y, predicted)$rho[1]))) %>% 
  mutate(best_model = NA) %>% 
  dplyr::select(dataset, .id, best_model, model, r_squared)

# > Table with measures of fit quality per curve and rank based on the r² (as performed in epifitter)
fit_quality_id <- rbind(perf_LIN, perf_NLIN, perf_FPCA) %>% 
  arrange(dataset, .id, desc(r_squared)) %>% 
  group_by(.id) %>% 
  mutate(best_model = 1:n()) %>% 
  mutate(type = case_when(
    model %in% c("Exponential", "Monomolecular", "Logistic", "Gompertz") ~ "Linear regression based on transformed data",
    model %in% c("nlin_Monomolecular", "nlin_Logistic", "nlin_Gompertz") ~ "Non-linear regression based on untransformed data",
    model == "FPCA" ~ "FPCA approach")) %>% 
  mutate(type = factor(type, levels = c("Linear regression based on transformed data", "Non-linear regression based on untransformed data", "FPCA approach"))) %>% 
  mutate(model = factor(model, 
                        levels = c("Exponential", "Monomolecular", "Logistic", "Gompertz",  "nlin_Monomolecular", "nlin_Logistic", "nlin_Gompertz",  "FPCA"),
                        labels = c("EXP", "MNM", "LOG", "GOM", "MNM", "LOG", "GOM", "FPCA")))


# RESULTS AGGREGATION ----------------------------------------------------------

fit_quality <- list("Total"     = fit_quality_tot, 
                    "Per_Curve" = fit_quality_id)


# SAVE RESULTS -----------------------------------------------------------------

save(fit_quality, file = paste0(path_project, "/00_DATA/RESULTS/02_FIT_QUALITY.rda"))




## Bonus: How to reproduce the output of epifitter : 
#
## lin
#res_epi_fit_lin$`02_COFFEE_CLR_INOC`$Data %>% 
#  group_by(model, .id) %>% 
#  summarise(r_squared = summary(stats::lm(linearized ~ time))$r.squared,
#            RSE       = summary(stats::lm(linearized ~ time))$sigma,
#            CCC       = DescTools::CCC(stats::lm(linearized ~ time)$fitted.values, 
#                                       linearized)$rho.c$est)  %>% head(.)
#
#res_epi_fit_lin$`02_COFFEE_CLR_INOC`$Parameters %>% arrange(model, .id) %>% dplyr::select(model, .id, r_squared, RSE, CCC) %>% head()
#
## nlin
#res_epi_fit_nlin$`02_COFFEE_CLR_INOC`$Data %>% 
#  group_by(model, .id) %>% 
#  summarise(#n = n(),
#    r_squared = caret::R2(pred = predicted, obs = y, na.rm=T),
#    #RSE       = sd(residual)/sqrt(n),
#    CCC       = DescTools::CCC(predicted, y)$rho.c$est)  %>% head(.)
#
#res_epi_fit_nlin$`02_COFFEE_CLR_INOC`$Parameters %>% arrange(model, .id) %>% dplyr::select(model, .id, r_squared, RSE, CCC) %>% head()


