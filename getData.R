## ---------------------------
##
## Script name: getData.R
##
## Purpose of script: Scrape data from gitHub repository established to track nCov20
##
## Author: Ben Phillips
##
## Date Created: 2020-03-24
##
## Email: phillipsb@unimelb.edu.au
##
## ---------------------------
##
## Notes:
##   
##
## --------------------------
## load up the packages we will need:  

library("readr")

## ---------------------------
## load up functions
source("functions.R")
focusCountry <- "Australia"
## ---------------------------

##### Case data #####
## Get case data
tsConf <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
tsDeath <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"

tsI<-read_csv(file = tsConf)
tsD<-read_csv(file = tsDeath)
#tsT<-read_csv(file = tsTesting)

## get Date range
dCols<-dateCols(tsI)
dates<-as.Date(colnames(tsI)[dCols], format = "%m/%d/%y")

## Tidy up names
names(tsI)[!dCols] <- make.names(names(tsI)[!dCols])
names(tsD)[!dCols] <- make.names(names(tsD)[!dCols])
#names(tsT)[!dCols] <- make.names(names(tsT)[!dCols])

# subset to Australia
tsI <- subset(tsI, tsI$Country.Region == focusCountry)
tsD <- subset(tsD, tsD$Country.Region == focusCountry)


## add recovery lag -- assumes all cases recover at 22 days
matI<-as.matrix(tsI[, dCols])
matD<-as.matrix(tsD[, dCols])
matA<-matI-matD #remove deaths
  matR <- cbind(matrix(0, nrow = nrow(matA), ncol = 22), matA[, -((ncol(matA)-21):ncol(matA))]) # recovered
matA <- matA - matR

tsA <- cbind(tsI[,!dCols], matA) # active cases


## Define menus
# get region names 
ddNames <- tsA$Province.State
ddReg <- ddNames
names(ddReg) <- ddNames


