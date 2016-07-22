####################################
#
# Nest success code
# Created by: Max Henschell
# Modified: 2 June 2015
#
####################################

#srm(list=ls()) #clean out the system

####################################
#
# Install required packages
#
####################################

packages.list<-c('MuMIn','BMA','xtable','RCurl', 'reshape', 'ggplot2', 'corrplot')

required.packages <- packages.list
new.packages <- required.packages[!(required.packages %in% installed.packages()[,'Package'])]
if(length(new.packages)) install.packages(new.packages, repos = 'http://cran.rstudio.com/')
lapply(packages.list, require, character.only=T)
# 
# library(BMA)
# library(xtable)
# Logistical Exposure Link Function
# See Shaffer, T.  2004. A unifying approach to analyzing nest success. 
# Auk 121(2): 526-540.
# library(MASS)
# logexp <- function(exposure = 1)
# {
#   linkfun <- function(mu) qlogis(mu^(1/exposure))
#   linkinv <- function(eta)  plogis(eta)^exposure
#   mu.eta <- function(eta) exposure * plogis(eta)^(exposure-1) *
#     .Call(stats:::C_logit_mu_eta, eta, PACKAGE = 'stats')
#   valideta <- function(eta) TRUE
#   link <- paste('logexp(', deparse(substitute(exposure)), ')',
#                 sep='')
#   structure(list(linkfun = linkfun, linkinv = linkinv,
#                  mu.eta = mu.eta, valideta = valideta, 
#                  name = link),
#             class = 'link-glm')
# }
# 
# panel.cor <- function(x, y, digits = 2, cex.cor, ...)
# {
#   usr <- par("usr"); on.exit(par(usr))
#   par(usr = c(0, 1, 0, 1))
#   # correlation coefficient
#   r <- cor(x, y)
#   txt <- format(c(r, 0.123456789), digits = digits)[1]
#   txt <- paste("r= ", txt, sep = "")
#   text(0.5, 0.6, txt)
#   
#   # p-value calculation
#   p <- cor.test(x, y)$p.value
#   txt2 <- format(c(p, 0.123456789), digits = digits)[1]
#   txt2 <- paste("p= ", txt2, sep = "")
#   if(p<0.01) txt2 <- paste("p= ", "<0.01", sep = "")
#   text(0.5, 0.4, txt2)
# }

####################################
#
# Source logistic exposure and correlation plot code from github
#
####################################

source_https <- function(url, ...) {
  # load package
  require(RCurl)
  # parse and evaluate each .R script
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file('CurlSSL', 'cacert.pem', package = 'RCurl'))), envir = .GlobalEnv)
  })
}

source_https('https://raw.github.com/maxhenschell/R-functions/master/logexp.R', 
             'https://raw.github.com/maxhenschell/R-functions/master/CorrPlots.R')

####################################
#
# WOTH
#
####################################

WOTH <- read.csv('DSR-WOTH.csv', header = T, na.strings = 'NA')
WOTH <- WOTH[WOTH$EXP != 0,]#Remove exp = 0 (first visit or visits after termination)
WOTH <- WOTH[complete.cases(WOTH),]
WOTH.MJP <- WOTH[WOTH$Study %in% "MJP",]
WOTH.MJP$Study <- factor(WOTH.MJP$Study)
WOTH.MJP$Year <- factor(WOTH.MJP$Year)
WOTH.MJP$PlotType <- factor(WOTH.MJP$PlotType)
WOTH.MJP$PatchType <- factor(WOTH.MJP$PatchType)
WOTH.MJP$GPH[WOTH.MJP$GPH == 0.00] <- (min(WOTH.MJP$GPH[WOTH.MJP$GPH > 0]))/2
WOTH.MJP$OrdDate <- WOTH.MJP$OrdDate-min(WOTH.MJP$OrdDate)+1
WOTH.MJP$OrdDate2 <- WOTH.MJP$OrdDate^2
WOTH.MJP$ParaStat <- as.factor(WOTH.MJP$ParaStat)
WOTH.MJP$Year <- as.factor(WOTH.MJP$Year)
WOTH.MJP1 <- WOTH.MJP
WOTH.MJP.nl <- WOTH.MJP1[!duplicated(WOTH.MJP[,1]),]

#Center and scale
head(WOTH.MJP[,c(12:26)])
str(WOTH.MJP)
for (i in 12:26){
  WOTH.MJP[,i] <- (WOTH.MJP[,i] - mean(WOTH.MJP[,i]))/sd(WOTH.MJP[,i])
}

####################################
#
# Bayesian Model Averaging
# Develop a full model, then run code
#
####################################
str(WOTH.MJP)
WOTH.MJP.BMA <- bic.glm(FATE ~ Core10km + Edge60m + OrdDate + CoreArea + ParaStat + TrailDist+ WBH + GPH + EdgeDist + PlotType, glm.family = binomial(logexp(exposure = WOTH.MJP$EXP)), data = WOTH.MJP)
imageplot.bma(WOTH.MJP.BMA)
summary(WOTH.MJP.BMA)
MJP.BMA.df <- data.frame(summary(WOTH.MJP.BMA))
