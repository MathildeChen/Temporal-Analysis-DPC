#' ---
#' title: "01. Models fitting"
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

# > Temporal analysis 
library(epifitter)

# > functional PCA
library(refund)  

# LOAD DATA --------------------------------------------------------------------

# > File containing data 
load(paste0(path_project, "00_DATA/04_MERGED_DATA.rda"))

# > Aggregate all data ready for analyses but keep the id of the disease and trial
data_for_temporal_analysis <- list_data_for_temporal_analysis %>% 
  map_dfr(., ~{
    
    .x$ready_for_analysis %>% 
      mutate(.id = as.character(.id))
    
    
  }, .id="dataset") %>% 
  separate(dataset, c("id_trial", "crop", "disease", "treatment"), sep = "_", remove = FALSE) 

# MODELS FITTING ---------------------------------------------------------------

# 1. Growth models with a linear regression applied on transformed (linearized) data

# > Models fitting
fit_lin <- data_for_temporal_analysis %>% 
  # Split between each dataset 
  split(.$dataset) %>% 
  map(., ~{
    
    # Format each dataset i to allow fit_multi to work
    ydata_no_na_i <- .x %>% 
      drop_na(.inc) %>%
      # Remove missing data 
      mutate(.inc_2 = if_else(.inc == 10^-6, .inc+10^-6, .inc)) %>% 
      # rescale the time 
      group_by(.id) %>% 
      mutate(min_index = min(.index)) %>% 
      mutate(.index_2 = .index - min_index)
    
    # Apply the function which automatically fit 
    # Exponential, mono-molecular, logistic, and Gompertz (EXP, MNM, LOG, GOM)
    fit_lin_all_models_i <- fit_multi(
      data          = ydata_no_na_i, # data
      time_col      = ".index_2",    # column for the time
      intensity_col = ".inc_2",      # column for the disease intensity
      strata_cols   = ".id",         # name or vector specifying the columns for stratification
      nlin          = FALSE          # linear approach
      #estimate_K    = TRUE,       
      #maxiter       = 200 
    )
  
  })

# > Fitted curves
curves_fit_lin <- fit_lin %>% 
  map_dfr(., ~{
    
    .x$Data
    
  }, .id = "dataset")

# -------------------------------------------
# 2. Growth models with a nonlinear regression applied on non-transformed data

# > Models fitting
fit_nlin <- data_for_temporal_analysis %>% 
  # Split between each dataset 
  split(.$dataset) %>% 
  map(., ~{
    
    # Format each dataset i to allow fit_multi to work
    ydata_no_na_i <- .x %>% 
      drop_na(.inc) %>%
      # Remove missing data 
      mutate(.inc_2 = if_else(.inc == 10^-6, .inc+10^-6, .inc)) %>% 
      # rescale the time 
      group_by(.id) %>% 
      mutate(min_index = min(.index)) %>% 
      mutate(.index_2 = .index - min_index)
    
    # Apply the function which automatically fit 
    # Exponential, mono-molecular, logistic, and Gompertz (EXP, MNM, LOG, GOM)
    fit_nlin_all_models_i <- fit_multi(
      data          = ydata_no_na_i, # data
      time_col      = ".index_2",    # column for the time
      intensity_col = ".inc_2",      # column for the disease intensity
      strata_cols   = ".id",         # name or vector specifying the columns for stratification
      nlin          = TRUE,          # nonlinear approach -> epidemic levels can be lower than 100%
      estimate_K    = TRUE,          # if nlin=TRUE, estimates maximum disease intensity. 
      maxiter       = 200            # n of iterations to compute K 
    )
    
    
  })

# > Fitted curves
curves_fit_nlin <- fit_nlin %>% 
  map_dfr(., ~{
    
    .x$Data
    
  }, .id = "dataset") %>%
  mutate(model = paste0("nlin_", model))

# -------------------------------------------
# 3. Functional principal component analysis (FPCA)

# > FPCA implementation
# For the main analysis, we set the number of basis to 10. 
# 10 set as default in mgcv package, on which refund is based. 
my_nbasis <- 10 

fit_fpca <- data_for_temporal_analysis %>% 
  # Split between each dataset 
  split(.$dataset) %>% 
  map(., ~{
    
    # Change the name of the variables to be compatible with the refund's syntax 
    ydata_i <- .x %>% 
      dplyr::select(.id, .index, .inc) %>% 
      rename(".value"=".inc") %>%
      # rescale the time
      group_by(.id) %>% 
      mutate(min_index = min(.index)) %>% 
      mutate(.index = .index - min_index) %>% 
      dplyr::select(-min_index)
    
    # Apply the function which automatically fit FPCA to the whole dataset
    fpca_fit_i = fpca.sc(ydata    = ydata_i,   # data
                         nbasis   = my_nbasis, # set for sensitivity analyses
                         pve      = 0.99,      # proportion of variance explained: used to choose the number of principal components
                         var      = TRUE,      # model-based estimates for the variance should be computed
                         center   = TRUE)      # centering the data by retrieving an estimated mean function
    
    # Get the reconstructed curves 
    fitted_data_FPCA_i <- NULL
    for(j in 1:length(unique(ydata_i$.id)))
    {
      # Same format than epifitter package to allow comparison between the outputs
      fitted_data_FPCA_i <- rbind(fitted_data_FPCA_i, 
                                data.frame(.id       = rep(unique(ydata_i$.id)[j], 
                                                           length(sort(unique(ydata_i$.index)))), 
                                           time      = sort(unique(ydata_i$.index)),
                                           y         = fpca_fit_i$Y[j, ], 
                                           model     = rep("FPCA", length(sort(unique(ydata_i$.index)))),
                                           predicted = fpca_fit_i$Yhat[j, ],
                                           residual  = fpca_fit_i$Y[j, ] - fpca_fit_i$Yhat[j, ])
      ) %>%
        # remove the points where we don't have any observed data to check th fit
        drop_na()
      
      
    }
    
    # Scores derivated from the functional PCs
    X_scores_i <- cbind(ydata_i %>% distinct(.id), 
                        as.data.frame(fpca_fit_i$scores))
    colnames(X_scores_i) <- c(".id", paste("fpca_score", 1:ncol(as.data.frame(fpca_fit_i$scores)), sep = ""))
    
    
    list("Parameters" = fpca_fit_i,
         "Data"       = fitted_data_FPCA_i, 
         "Scores"     = X_scores_i)
    
  })

# > Fitted curves
curves_fpca_fit <- fit_fpca %>% 
  map_dfr(., ~{
    
    .x$Data
    
  }, .id = "dataset")


# RESULTS AGGREGATION ----------------------------------------------------------
# > Fits 
all_fits <- list("LIN"  = fit_lin,
                 "NLIN" = fit_nlin,
                 "FPCA" = fit_fpca)

# > Fitted curves
all_curves_fit <- rbind(curves_fit_lin  %>% dplyr::select(-linearized),
                        curves_fit_nlin,
                        curves_fpca_fit) %>% 
  mutate(type = case_when(
    model %in% c("Exponential", "Monomolecular", "Logistic", "Gompertz") ~ "Linear regression based on transformed data",
    model %in% c("nlin_Monomolecular", "nlin_Logistic", "nlin_Gompertz") ~ "Non-linear regression based on untransformed data",
    model == "FPCA" ~ "FPCA approach")) %>% 
  mutate(type = factor(type, levels = c("Linear regression based on transformed data", "Non-linear regression based on untransformed data", "FPCA approach"))) %>% 
  mutate(model = factor(model, 
                        levels = c("Exponential", "Monomolecular", "Logistic", "Gompertz",  "nlin_Monomolecular", "nlin_Logistic", "nlin_Gompertz",  "FPCA"),
                        labels = c("EXP", "MNM", "LOG", "GOM", "MNM", "LOG", "GOM", "FPCA")))


# SAVE RESULTS -----------------------------------------------------------------

save(all_fits, file = paste0(path_project, "/00_DATA/RESULTS/FITS.rda"))

save(all_curves_fit, file = paste0(path_project, "/00_DATA/RESULTS/FITTED_CURVES.rda"))
