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
# AMRO
#
####################################

AMRO <- read.csv('DSR-AMRO.csv', header = T, na.strings = 'NA')
# AMRO$PatchArea <- as.numeric(levels(AMRO$PatchArea))[AMRO$PatchArea]
# AMRO$CoreArea <- as.numeric(levels(AMRO$CoreArea))[AMRO$CoreArea]
AMRO <- AMRO[AMRO$EXP != 0,]#Remove exp = 0 (first visit or visits after termination)
AMRO <- AMRO[complete.cases(AMRO[,c('FATE','EXP','Year','TrailDist', 'WTR', 'WBH', 'GPH','CoreArea', 
                                    'RoadDist', 'EdgeDist')]),]
AMRO.MJP <- AMRO[AMRO$Study %in% "MJP",]
AMRO.MJP$Study <- factor(AMRO.MJP$Study)
AMRO.MJP$Year <- factor(AMRO.MJP$Year)
AMRO.MJP$PlotType <- factor(AMRO.MJP$PlotType)
AMRO.MJP$PatchType <- factor(AMRO.MJP$PatchType)
AMRO.MJP$GPH[AMRO.MJP$GPH == 0.00] <- (min(AMRO.MJP$GPH[AMRO.MJP$GPH > 0]))/2
AMRO.MJP$OrdDate <- AMRO.MJP$OrdDate-min(AMRO.MJP$OrdDate)+1
AMRO.MJP$OrdDate2 <- AMRO.MJP$OrdDate^2
AMRO.MJP$ParaStat <- as.factor(AMRO.MJP$ParaStat)
AMRO.MJP$Year <- as.factor(AMRO.MJP$Year)
AMRO.MJP1 <- AMRO.MJP
AMRO.MJP.nl <- AMRO.MJP1[!duplicated(AMRO.MJP[,1]),]

#Center and scale
head(AMRO.MJP[,c(13:25)])
str(AMRO.MJP)
for (i in 13:25){AMRO.MJP[,i] <- (AMRO.MJP[,i] - mean(AMRO.MJP[,i]))/sd(AMRO.MJP[,i])}

####################################
#
# Bayesian Model Averaging
# Develop a full model, then run code
#
####################################
str(AMRO.MJP)
AMRO.MJP.BMA <- bic.glm(FATE ~ Core10km + Edge60m + OrdDate + CoreArea + TrailDist+ WBH + GPH + EdgeDist +PatchType + PlotType, glm.family = binomial(logexp(exposure = AMRO.MJP$EXP)), data = AMRO.MJP)
imageplot.bma(AMRO.MJP.BMA)
summary(AMRO.MJP.BMA)
MJP.BMA.df <- data.frame(summary(AMRO.MJP.BMA))
