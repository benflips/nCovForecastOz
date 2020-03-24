## ---------------------------
##
## Script name: parameters.R
##
## Purpose of script: load up parameters needed
##
## Author: Ben Phillips
##
## Date Created: 2020-03-24
##
## Email: phillipsb@unimelb.edu.au
##
## ---------------------------
##
## Notes: Most of these are taken from Alison Hill's curated spreadsheet (fully referenced)
##   Format will be {best, high low}
##
## --------------------------
## load up the packages we will need 

## ---------------------------

# proportion symptomatic not requiring hospitalisation
fracMild <- c(0.8, 0.27, 0.81)
# proportion requiring hospitalisation but not ICU
fracSevere <- c(0.15, 0.14, 0.24)
# proportion requiring ICU
fracCritical <- c(0.05, 0.02, 0.06)

# duration of hospitalisation severe infection
durSevInf <- c(6, 5, 14)
# duration of hospitalisation, critical infection (includes + 8 days ICU)
durCritInf <- c(14, 13, 22)

# time from symptom onset to hospitalisation (days)
t_H <- c(6, 4.5, 12)

# case fatality ratio (symptomatic only)
cfr <- c(0.033)



