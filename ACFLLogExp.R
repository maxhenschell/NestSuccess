rm(list=ls())

packages.list<-c("MASS","MuMIn","lattice","ggplot2", "xlsx","RCurl","repmis","BMA","xtable","lme4")
required.packages <- packages.list
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.rstudio.com/")
lapply(packages.list, require, character.only=T)

# #Get logexp function from dropbox
# setwd(tempdir())
# destfile = "logexp.R"
# x = getBinaryURL("https://dl.dropboxusercontent.com/u/54200206/R/functions/logexp.R", followlocation = TRUE, ssl.verifypeer = FALSE)
# writeBin(x, destfile, useBytes = TRUE)
# source(paste(tempdir(), "/logexp.R", sep = ""))

# Logistical Exposure Link Function
# See Shaffer, T.  2004. A unifying approach to analyzing nest success. 
# Auk 121(2): 526-540.
library(MASS)
logexp <- function(exposure = 1)
{
  linkfun <- function(mu) qlogis(mu^(1/exposure))
  linkinv <- function(eta)  plogis(eta)^exposure
  mu.eta <- function(eta) exposure * plogis(eta)^(exposure-1) *
    .Call(stats:::C_logit_mu_eta, eta, PACKAGE = "stats")
  valideta <- function(eta) TRUE
  link <- paste("logexp(", deparse(substitute(exposure)), ")",
                sep="")
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, 
                 name = link),
            class = "link-glm")
}

#setwd("D:\\Dropbox\\R\\RStudioProjects\\NestSuccess")#Widnows
#setwd("/mnt/Removable/Dropbox/R/RStudioProjects/NestSuccess")#Maxtop ubuntu
setwd("/home/max/Dropbox/R/RStudioProjects/NestSuccess")
ACFL <- read.csv("LogExpACFL.csv", header = T, na.strings = "NA")

str(ACFL)
str(ACFL <- ACFL[ACFL$EXP != 0,])#Remove exp = 0 (first visit)

ACFL <- ACFL[complete.cases(ACFL[,"DIST"]),] #remove missing distance
ACFL <- ACFL[complete.cases(ACFL[,"HEIGHT"]),] #remove missing height
ACFL <- ACFL[complete.cases(ACFL[,"OrdDate"]),] #remove missing Ordinal date
ACFL <- ACFL[complete.cases(ACFL[,c("W1","W2")]),] #remove missing width
ACFL <- ACFL[complete.cases(ACFL[,c("EXP")]),] #remove missing exposure dat
ACFL <- ACFL[complete.cases(ACFL[,c("FATE")]),] #remove missing fates (checks after completion)
str(ACFL)

source("/home/max/Dropbox/R/R-cheatScripts/CorrPlotRVals.R")
pdf("VariableCorrelationPlot.pdf")
pairs(ACFL[,c("DIST","HEIGHT","TYPE","W1","W2")], upper.panel = panel.cor)
dev.off()
write(cor(ACFL[,c("DIST","HEIGHT","W1","W2")]), "VariableCorrelation.txt", sep = ",")
ACFL.fate <- glm(FATE ~ OrdDate + OrdDate^2 + HEIGHT*DIST +  W1*DIST + W2*DIST + TYPE*DIST + W1*TYPE + W2*TYPE, family = binomial(logexp(exposure = as.numeric(ACFL$EXP))), data = ACFL)
summary(ACFL.fate)

ACFL.mods <- dredge(ACFL.fate)
sink("ACFLmods.txt")
ACFL.mods
sink()

#Build some a priori, ecologically meaningful models

ACFL.mods <- list() #list of models
ACFL.sum <- list() # list of model summaries

ACFL.mods[1] <- glm(FATE ~ OrdDate + OrdDate^2 + 
                      HEIGHT*DIST + + DIST*TYPE +
                      W1*DIST + W1*TYPE, 
                    family = binomial(logexp(exposure = as.numeric(ACFL$EXP))), data = ACFL)#FULL MODEL with tread width

ACFL.mods[1] <- glm(FATE ~ OrdDate + OrdDate^2 + 
                      HEIGHT*DIST + + DIST*TYPE +
                      W2*DIST + W2*TYPE, 
                    family = binomial(logexp(exposure = as.numeric(ACFL$EXP))), data = ACFL)#FULL MODEL with WBH

ACFL.mods[2] <- glm(FATE ~ OrdDate + OrdDate^2,
                    family = binomial(logexp(exposure = as.numeric(ACFL$EXP))), data = ACFL)#Ordinal date

ACFL.mods[2] <- glm(FATE ~ W1*DIST,
                    family = binomial(logexp(exposure = as.numeric(ACFL$EXP))), data = ACFL)#tread width + distance

ACFL.mods[2] <- glm(FATE ~ W2*DIST,
                    family = binomial(logexp(exposure = as.numeric(ACFL$EXP))), data = ACFL)#WBH + distnace

