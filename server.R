## ---------------------------
##
## Script name: 
##
## Purpose of script:
##
## Author: Ben Phillips
##
## Date Created: 2020-03-12
##
## Email: phillipsb@unimelb.edu.au
##
## ---------------------------
##
## Notes:
##   
##
## --------------------------
## load up the packages we will need 
library(shiny)
## ---------------------------

## source files
source("getData.R")

## ---------------------------
options(scipen=9)


# Define server logic 
shinyServer(function(input, output) {
  ##### Raw stats #####  
  output$rawStats <- renderTable({
    yA <- tsSub(tsA,tsA$Province.State %in% input$stateFinder)
    yD <- tsSub(tsD,tsD$Province.State %in% input$stateFinder)
    yI <- tsSub(tsI,tsI$Province.State %in% input$stateFinder)
    #yR <- tsSub(tsR,tsR$Province.State %in% input$stateFinder)
    nn <-length(yI)
    if (is.na(yA[nn])) nn <- nn-1
    out <- as.integer(c(yI[nn], yD[nn]))
    dim(out) <-c(1,2)
    colnames(out) <- c("Total", "Deaths")
    format(out, big.mark = ",")
  }, rownames = FALSE)
  
##### Raw plot #####  
  output$rawPlot <- renderPlot({
    yA <- tsSub(tsA,tsA$Province.State %in% input$stateFinder)
    lDat <- projSimple(yA, dates)
    yMax <- max(c(lDat$y[,"fit"], yA), na.rm = TRUE)
    yTxt <- "Confirmed active cases"
    plot(yA~dates, 
         xlim = c(min(dates), max(lDat$x)),
         ylim = c(0, yMax),
         pch = 19, 
         bty = "u", 
         xlab = "Date", 
         ylab = yTxt,
         main = input$stateFinder)
    axis(side = 4)
    lines(lDat$y[, "fit"]~lDat$x)
    lines(lDat$y[, "lwr"]~lDat$x, lty = 2)
    lines(lDat$y[, "upr"]~lDat$x, lty = 2)
  })
  
##### Log plot #####    
  output$logPlot <- renderPlot({
    yA <- tsSub(tsA,tsA$Province.State %in% input$stateFinder)
    lDat <- projSimple(yA, dates)
    yMax <- max(c(lDat$y[,"fit"], yA), na.rm = TRUE)
    yTxt <- "Confirmed active cases (log scale)"
    plot((yA+0.1)~dates, 
         xlim = c(min(dates), max(lDat$x)),
         ylim = c(1, yMax),
         log = "y",
         pch = 19, 
         bty = "u", 
         xlab = "Date", 
         ylab = yTxt,
         main = input$stateFinder)
    axis(side=4)
    lines(lDat$y[, "fit"]~lDat$x)
    lines(lDat$y[, "lwr"]~lDat$x, lty = 2)
    lines(lDat$y[, "upr"]~lDat$x, lty = 2)
  })
  
##### Detection rate #####    
  output$detRate <- renderText({
    yD <- tsSub(tsD,tsD$Province.State %in% input$stateFinder)
    yI <- tsSub(tsI,tsI$Province.State %in% input$stateFinder)
    dR<-round(detRate(yI, yD), 4)
    if (is.na(dR)) "Insufficient data for estimation" else dR
  })
  
##### Prediction table confirmed #####    
  output$tablePredConf <- renderTable({
    yA <- tsSub(tsA,tsA$Province.State %in% input$stateFinder)
    lDat <- projSimple(yA, dates)
    nowThen <- format(as.integer(c(tail(yA[!is.na(yA)], 1), tail(lDat$y[,"lwr"],1), tail(lDat$y[,"upr"],1))), big.mark = ",")
    nowThen <- c(nowThen[1], paste(nowThen[2], "-", nowThen[3]))
    dim(nowThen) <- c(1, 2)
    colnames(nowThen)<-c("Now", "In 10 days (min-max)")
    nowThen
  }, rownames = FALSE)
  
##### Prediction table true #####    
  output$tablePredTrue <- renderText({
    yA <- tsSub(tsA,tsA$Province.State %in% input$stateFinder)
    yD <- tsSub(tsD,tsD$Province.State %in% input$stateFinder)
    yI <- tsSub(tsI,tsI$Province.State %in% input$stateFinder)
    dRate <- detRate(yI, yD)
    lDat <- projSimple(yA, dates)
    now <- tail(yA[!is.na(yA)], 1)
    nowTrue <- format(round(now/dRate, 0), big.mark = ",")
    nowTrue
  })
  
##### Curve-flattenning #####    
  output$cfi <- renderPlot({
    pDat <- subset(tsA, tsA$Province.State %in% input$stateGrowthRate)
    pMat<-as.matrix(log(pDat[,-(1:4)]))
    row.names(pMat)<-pDat$Province.State
    cfiDat<-apply(pMat, MARGIN = 1, FUN = "cfi")
    cfiDat[!is.finite(cfiDat)]<-0
    clrs<-hcl.colors(length(pDat$Province.State))
    dateSub<-3:length(dates) # date tsSub
    plot(cfiDat[,1]~dates[dateSub], 
         type = "n", 
         ylim = range(c(-1.2,1.2)*sd(cfiDat)),
         bty = "l",
         xlab = "Date",
         ylab = "Curve-flatenning index")
    abline(a = 0, b = 0, lty = 2, lwd = 2)
    for (cc in 1:ncol(cfiDat)){
      cfiSmooth<-loess(cfiDat[,cc]~as.numeric(dates[dateSub]))
      lines(cfiSmooth$fitted~dates[dateSub], col = clrs[cc], lwd=2)
    }
    legend("topleft", 
           legend = pDat$Province.State,
           col = clrs,
           lty = 1,
           bty = "n")
  })
##### Growth rate #####    
  output$growthRate <- renderPlot({
    pDat <- subset(tsA[,dCols], tsA$Province.State %in% input$stateGrowthRate)
    gRate <- as.matrix(growthRate(pDat))
    clrs<-hcl.colors(length(pDat$Province.State))
    dates10 <- dates[(length(pDat)-10+1):length(pDat)]
    counts <- table(gRate)
    barplot(gRate,
            main="Growth rate",
            xlab="Date", 
            ylab="Growth rate (% per day)",
            beside=TRUE,
            col = clrs,
            legend = pDat$Province.State,
            args.legend = list(bty = "n", x = "topleft"))
  })
  
##### Doubling time ##### 
  output$doubTime <- renderText({
    pDat <- tsSub(tsA, tsA$Province.State %in% input$stateFinder)
    dTime <- round(doubTime(pDat, dates), 1)
  })
  
##### Doubling time plot #####    
  output$doubTimePlot <- renderPlot({
    pDat <- tsSub(tsA, tsA$Province.State %in% input$stateGrowthRate)
    dTime <- as.matrix(doubTime(pDat))
    dTime[!is.finite(dTime)]<-NA
    clrs<-hcl.colors(length(input$stateGrowthRate))
    dates10 <- dates[(length(pDat)-10+1):length(pDat)]
    counts <- table(dTime)
    barplot(dTime,
            main="Doubling time",
            xlab="Date", 
            ylab="Doubling time (days)",
            beside=TRUE,
            col = clrs,
            legend = input$stateGrowthRate,
            args.legend = list(bty = "n", x = "topleft"))
  })

})
