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
# ACFL
#
####################################

ACFL <- read.csv('DSR-ACFL.csv', header = T, na.strings = 'NA')

ACFL <- ACFL[ACFL$EXP != 0,]#Remove exp = 0 (first visit or visits after termination)
ACFL <- ACFL[complete.cases(ACFL),]
ACFL.MJP <- ACFL[ACFL$Study %in% "MJP",]
ACFL.MJP$Study <- factor(ACFL.MJP$Study)
ACFL.MJP$Year <- factor(ACFL.MJP$Year)

ACFL.MJP$GPH[ACFL.MJP$GPH == 0.00] <- (min(ACFL.MJP$GPH[ACFL.MJP$GPH > 0]))/2
ACFL.MJP$OrdDate <- ACFL.MJP$OrdDate-min(ACFL.MJP$OrdDate)+1
ACFL.MJP$OrdDate2 <- ACFL.MJP$OrdDate^2
ACFL.MJP$ParaStat <- as.factor(ACFL.MJP$ParaStat)
ACFL.MJP$Year <- as.factor(ACFL.MJP$Year)
ACFL.MJP1 <- ACFL.MJP
ACFL.MJP.nl <- ACFL.MJP1[!duplicated(ACFL.MJP[,1]),]

#Center and scale
for (i in 12:26){
  ACFL.MJP[,i] <- (ACFL.MJP[,i] - mean(ACFL.MJP[,i]))/sd(ACFL.MJP[,i])
}

####################################
#
# Bayesian Model Averaging
# Develop a full model, then run code
#
####################################

ACFL.MJP.BMA <- bic.glm(FATE ~ 
                          CoreArea +
                          OrdDate*GPH + 
                          ParaStat*WBH + 
                          ParaStat*OrdDate + 
                           WBH*GPH, glm.family = binomial(logexp(exposure = ACFL.MJP$EXP)), data = ACFL.MJP)
imageplot.bma(ACFL.MJP.BMA)
summary(ACFL.MJP.BMA)

MJP.BMA.df <- data.frame(summary(ACFL.MJP.BMA))