#' ---
#' title: "00. Functions"
#' author: "Mathilde Chen & Jacques Avelino, CIRAD"
#' date: "September 2026"
#' ---

# ----------------------------------
# Packages

# > General data management 
library(tidyverse) 

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







