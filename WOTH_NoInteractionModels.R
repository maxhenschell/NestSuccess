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

packages.list<-c('MuMIn','BMA','xtable','RCurl', 'reshape', 'ggplot2')

required.packages <- packages.list
new.packages <- required.packages[!(required.packages %in% installed.packages()[,'Package'])]
if(length(new.packages)) install.packages(new.packages, repos = 'http://cran.rstudio.com/')
# packages.list2<-c('coefplot2')
# required.packages2 <- packages.list2
# new.packages2 <- required.packages2[!(required.packages2 %in% installed.packages()[,'Package'])]
# if(length(new.packages2)) install.packages(new.packages2, ,
#                                            repos='http://www.math.mcmaster.ca/bolker/R',
#                                            type='source')
# packages.list<-c('BMA','xtable','RCurl', 'bbmle', 'reshape', 'ggplot2', 'BEST', 'arm', 'MCMCpack', 'popbio', 'extrafont', 'coefplot2')
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
WOTH <- WOTH[complete.cases(WOTH),]#[,c('FATE','EXP','Year','TrailDist', 'WTR', 'WBH', 'GPH', 'ParaStat','RoadDist', 'EdgeDist')]),]#remove nests with missing data
#GPH.m <- mean(WOTH$GPH[WOTH$GPH > 0])
# WOTH$GPH.m <- WOTH$GPH
# GPH.m <- mean(WOTH$GPH[WOTH$GPH > 0])
#WOTH$GPH.m[WOTH$GPH.m == 0.00] <- GPH.m 
WOTH$GPH[WOTH$GPH == 0.00] <- (min(WOTH$GPH[WOTH$GPH > 0]))/2
WOTH$GPH[WOTH$Study == 'MP'] <- 0

WOTH$OrdDate <- WOTH$OrdDate-min(WOTH$OrdDate)+1
WOTH$OrdDate2 <- WOTH$OrdDate^2
WOTH <- WOTH[WOTH$Study %in% c('MAH', 'MJP'),]
WOTH$Study <- factor(WOTH$Study)
#WOTH <- WOTH[WOTH$TrailDist < 400,]
WOTH$ParaStat <- as.factor(WOTH$ParaStat)
WOTH$Year <- as.factor(WOTH$Year)
#WOTH$CoreArea <- (pi*10^2)*WOTH$Core10km
#WOTH$CoreArea <- WOTH$CoreArea/100
#WOTH$EdgeArea <- (pi*.5^2)*WOTH$Edge60m
#WOTH$EdgeArea <- WOTH$EdgeArea/100
WOTH1 <- WOTH
WOTH.nl <- WOTH1[!duplicated(WOTH[,1]),]

#Center and scale
str(WOTH)
head(WOTH[,c(12:26)])
for (i in 12:26){
  WOTH[,i] <- (WOTH[,i] - mean(WOTH[,i]))/sd(WOTH[,i])
}

####################################
#
# Bayesian Model Averaging
# Develop a full model, then run code
#
####################################


#Marty
WOTH.MJP <- WOTH[WOTH$Study %in% "MJP",]
WOTH.MJP$Study <- factor(WOTH.MJP$Study)
WOTH.MJP$Year <- factor(WOTH.MJP$Year)
WOTH.MJP$PlotType <- factor(WOTH.MJP$PlotType)
WOTH.MJP.BMA <- bic.glm(FATE ~ Core10km + Edge60m + OrdDate + CoreArea + ParaStat + TrailDist+ WBH + GPH + EdgeDist + PlotType, glm.family = binomial(logexp(exposure = WOTH.MJP$EXP)), data = WOTH.MJP)
imageplot.bma(WOTH.MJP.BMA)
summary(WOTH.MJP.BMA)
MJP.BMA.df <- data.frame(summary(WOTH.MJP.BMA))

#Max
WOTH.MAH <- WOTH[WOTH$Study %in% "MAH",]
WOTH.MAH$Study <- factor(WOTH.MAH$Study)
WOTH.MAH$Year <- factor(WOTH.MAH$Year)
WOTH.MAH.BMA <- bic.glm(FATE ~ Year + Core10km + Edge60m + OrdDate + CoreArea + PatchType + ParaStat + TrailDist+ WBH + GPH + EdgeDist + PlotType, glm.family = binomial(logexp(exposure = WOTH.MAH$EXP)), data = WOTH.MAH)
imageplot.bma(WOTH.MAH.BMA)
summary(WOTH.MAH.BMA)
MAH.BMA.df <- data.frame(summary(WOTH.MAH.BMA))