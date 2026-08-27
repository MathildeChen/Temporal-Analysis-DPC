# A framework based on functional principal component analysis for comparing disease progress curves 

## Author and affiliations
Mathilde CHEN

CIRAD, UMR PHIM, F-97130 Capesterre-Belle-Eau, Guadeloupe, France. 

PHIM, CIRAD, INRAE, Institut Agro, IRD, Université de Montpellier, Montpellier, France. 

## Summary
Repository including the scripts used for comparison of disease progress curves (DPC) using contrasting statistical methods. 

 
5-step procedure (summarized in the attached figure):  
1. Description of the data and analyses  
2. Fitting of linear models (Exponential, Logistic, Gompertz, and Monomolecular), nonlinear models (Logistic, Gompertz, and Monomolecular), and FPCA. 
3. Evaluate the goodness of fit for the different models. Calculate various metrics (e.g., RSS, r²) to estimate the overall goodness of fit of the model (based on errors calculated point by point and curve by curve).   
4. Comparison of parameter values across treatments; ANOVA/MANOVA + linear models  
5. Biological interpretation of parameters 
<img width="4706" height="2658" alt="Figure1_workflow_V2" src="https://github.com/user-attachments/assets/db07200d-6345-41b1-8ad7-417591c9c04c" />
