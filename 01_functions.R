#' ---
#' title: "00. Functions"
#' author: "Mathilde Chen & Jacques Avelino, CIRAD"
#' date: "September 2026"
#' ---


# ----------------------------------
# Packages

# > General data management 
library(tidyverse) 

# > Temporal analysis 
library(epifitter)

# > functional PCA
library(refund)  

# > measures of fit quality 
library(caret)  
library(hydroGOF)


# -------------------------------------------
# get_density: Get density of points in 2 dimensions. Returns the density within each square.
#   - "x": A numeric vector.
#   - "y": A numeric vector.
#   - "n": Create a square n by n grid to compute density.

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
} 

# ----------------------------------
# data_to_fit: a dataset containing at least the following variables:
#   -".index",    # column for the time (can vary across the individuals)
#   -".inc",      # column for the disease intensity in proportion (ranging from 0 to 1)
#   -".id",       # name or vector specifying the columns for stratification
# fpca_nbasis: the number of basis used to fit the FPCA (need to be optimised)
# fpca_npc:    the number of functional components in the FPCA (need to be optimised)


fit_multi_fpca <- function(data_to_fit, 
                           fpca_nbasis = NULL,
                           fpca_npc = NULL)
{
  
  # ----------------------------------
  # 0. checking data and packages
  
  # > packages for temporal analysis 
  require(epifitter) ; require(refund)
  # > packages for measure of fit quality
  require(caret) ; library(hydroGOF)
  
  # > incidence variable format
  if(max(data_to_fit$.inc)>1 | min(data_to_fit$.inc)<0)
  {
    
    stop("Enter the incidence variable as proportion\n(ranging between 0 and 1)\n
         Or the wrong variable is entered as incidence. Please check.")
    
  }
  
  # > replacing 0 by very small values
  if(min(data_to_fit$.inc, na.rm=T)==0)
  {
    
    data_to_fit <- data_to_fit %>%
      # > replace 0 values by very small values 
      mutate(.inc = if_else(.inc == 0, 10^-6, .inc))
    
    warning("The 0 values in incidence variable have been replaced by 10^-6.")
    
  }
  
  # > removing any NAs in the data 
  if(anyNA(data_to_fit$.inc)==TRUE)
  {
    
    data_to_fit <- data_to_fit %>% drop_na()
    warning("NAs values in .inc have been removed.")
    
  }
  
  # > check FPCA parameters values
  if(is.na(fpca_nbasis) |  is.na(fpca_npc))
  {
    
    stop("Enter parameters to fit the FPCA approach (fpca_nbasis and fpca_npc).")
    
  }
  
  # ----------------------------------
  # 1. Modeling epidemics using the linear regression approach
  # tested approaches: Exponential, Gompterz, Logistic, and Monomolecular
  
  # > fit
  all_fit_lin <- fit_multi(
    time_col      = ".index",    # column for the time
    intensity_col = ".inc",      # column for the disease intensity
    strata_cols   = ".id",       # name or vector specifying the columns for stratification
    data          = data_to_fit, # data
    nlin          = FALSE        # estimates nonlinear approach (if FALSE = fit a linear approach)
  )
  
  # > fitted data and residuals
  data_fit_lin <- all_fit_lin$Data %>% 
    dplyr::select(.id, time, y, model, predicted, residual)
  
  
  # ----------------------------------
  # 2. Modeling epidemics using the non-linear regression approach
  # uses the Levenberg-Marquardt algorithm for least-squares estimation of nonlinear parameters
  # tested approaches: Exponential, Gompterz, Logistic, and Monomolecular
  
  # > fit
  all_fit_nlin <- fit_multi(
    time_col      = ".index",    # column for the time
    intensity_col = ".inc",      # column for the disease intensity
    strata_cols   = ".id",       # name or vector specifying the columns for stratification
    data          = data_to_fit, # data
    nlin          = TRUE         # estimates nonlinear approach (if TRUE = fit a non-linear approach)
  )
  
  # > fitted data and residuals
  data_fit_nlin <- all_fit_nlin$Data %>% 
    dplyr::select(.id, time, y, model, predicted, residual)
  
  
  # ----------------------------------
  # 3. Functional principal component analysis (FPCA)
  # > perform the functional PCA
  fit_fpca = fpca.sc(ydata  = data_to_fit,  # data
                    nbasis  = fpca_nbasis,
                    npc     = fpca_npc, 
                    simul   = T,
                    var     = TRUE,   # model-based estimates for the variance should be computed
                    center  = TRUE) 
  
  
  # > fitted data and residuals
  data_fit_FPCA <- NULL
  for(i in 1:length(unique(data_to_fit$.id)))
  {
    data_fit_FPCA <- rbind(data_fit_FPCA, 
                           data.frame(.id       = rep(unique(data_to_fit$.id)[i], length(sort(unique(data_to_fit$.index)))), 
                                      time      = sort(unique(data_to_fit$.index)),
                                      y         = fit_fpca$Y[i, ], 
                                      model     = rep("FPCA", length(sort(unique(data_to_fit$.index)))),
                                      predicted = fit_fpca$Yhat[i, ],
                                      residual  = fit_fpca$Y[i, ] - fit_fpca$Yhat[i, ]))
  }
  
  # ----------------------------------
  # 4 aggregate the different results
  
  data_fit_all <- rbind(data_fit_lin  %>% mutate(fit_method = "Linear approach"), 
                        data_fit_nlin %>% mutate(fit_method = "Non-ninear approach"),
                        data_fit_FPCA %>% mutate(fit_method = "FPCA")) 
  
  
  results <- list(fits = list("lin"  = all_fit_lin,
                              "nlin" = all_fit_nlin,
                              "fpca" = fit_fpca),
                  data_fit = data_fit_all)
  
  return(results)
  
}







