####################################
#
# Nest success for Acadian Flycatchers 
# in the Baraboo Hills
# Created by: Max Henschell
# Modified: 20 March 2014 
#
####################################

rm(list=ls()) #clean out the system

####################################
#
# Install required packages
#
####################################

packages.list<-c("popbio","ResourceSelection","pROC","MASS","MuMIn","lattice","ggplot2","BMA","xtable","lme4","RCurl", "repmis", "xtable", "knitr", "aod")
required.packages <- packages.list
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.rstudio.com/")
lapply(packages.list, require, character.only=T)

# Logistical Exposure Link Function
# See Shaffer, T.  2004. A unifying approach to analyzing nest success. 
# Auk 121(2): 526-540.
# library(MASS)
# logexp <- function(exposure = 1)
# {
#   linkfun <- function(mu) qlogis(mu^(1/exposure))
#   linkinv <- function(eta)  plogis(eta)^exposure
#   mu.eta <- function(eta) exposure * plogis(eta)^(exposure-1) *
#     .Call(stats:::C_logit_mu_eta, eta, PACKAGE = "stats")
#   valideta <- function(eta) TRUE
#   link <- paste("logexp(", deparse(substitute(exposure)), ")",
#                 sep="")
#   structure(list(linkfun = linkfun, linkinv = linkinv,
#                  mu.eta = mu.eta, valideta = valideta, 
#                  name = link),
#             class = "link-glm")
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
    eval(parse(text = getURL(u, followlocation = TRUE, cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}

source_https("https://raw.github.com/maxhenschell/R-functions/master/logexp.R", 
             "https://raw.github.com/maxhenschell/R-functions/master/CorrPlots.R")

####################################
#
# Download data from dropbox
#
####################################

#ACFL <- source_data("https://www.dropbox.com/s/wv2c57ketiek6ac/DSR-ACFL.csv")
ACFL <- read.csv("DSR-ACFL.csv", header = T, na.strings = "NA")

####################################
#
#Clean up data
#
####################################

ACFL <- ACFL[ACFL$EXP != 0,]#Remove exp = 0 (first visit or visits after termination)
ACFL <- ACFL[complete.cases(ACFL[,c("FATE","EXP","NestHeight", "TrailDist", "WTR", "WBH", "GPH", 
                                    "RoadDist", "EdgeDist", "RoadDensity", "TrailType")]),]
        #remove nests with missing data
ACFL$GPH[ACFL$GPH == 0] <- min(ACFL$GPH)/2
ACFL$OrdDate <- ACFL$OrdDate-min(ACFL$OrdDate)+1
####################################
#
#Correlation
#
####################################

pairs(ACFL[,c("NestHeight", "TrailDist", "TrailType", "WTR", "WBH", "GPH", 
              "RoadDist", "EdgeDist", "RoadDensity")], upper.panel = panel.cor)

dat.mah <- ACFL[ACFL$Study == "MAH",]
dat.mah <- dat.mah[complete.cases(dat.mah[,c("NestHeight", "TrailDist", "TrailType", "WTR", "WBH", "GPH", 
                                             "RoadDist", "EdgeDist", "RoadDensity""TPH", "TrailType")]),]
pairs(dat.mah[,c("NestHeight", "TrailDist", "TrailType", "WTR", "WBH", "GPH", 
                 "RoadDist", "EdgeDist", "RoadDensity""TPH", "BA_ha", "Stem_m2", "TrailType")], upper.panel = panel.cor)

dat.mjp <- ACFL[ACFL$STUDY == "MJP",]
pairs(dat.mjp[,c("DIST", "HEIGHT", "WTR", "WBH", "GPH", "BA_ha", "TYPE")], upper.panel = panel.cor)

dat.al <- ACFL[ACFL$STUDY == "AL",]
pairs(dat.al[,c("DIST", "HEIGHT", "WTR", "WBH", "GPH")], upper.panel = panel.cor)


####################################
#
# Bayesian Model Averaging
# Develop a full model, then run code
#
####################################

# variables:
# Forest1km: proportion forest cover within 1 km of nest -> all studies
# Road1km: km of roads within 1 km of nest -> all studies
# Homesteads1km : number of homesteads within 1 km of nest -> all studies
# RoadDist: distance to closest road -> all studies
# W2: trail width at 1.3 m from ground -> all studies
# GPH: groups per hour -> MAH,MJP
# TPH: trees pre hectare around nest -> MAH
# BA_ha: basal area per hectare around nest -> MAH,MJP
# TYPE: trail type (constructed, old road, (social?)) -> all studies
# PLOT: trail name -> all studies 
# HEIGHT: nest height -> all studies
# DIST: distance from nest to trail -> all studies
# ORIENT: Nest orientation from trail (90 = opposite side of tree from trail)
# OrdDate: Ordinal date starting from Jan 1 (=1) -> all studies
# STUDY: study
# YEAR: year of nest
# NEST: nest ID

# f.ls <- formula(FATE ~ Forest1km + Road1km + Homesteads1km + RoadDist)
# f.hab <- formula(FATE ~ W2*GPH + TPH + BA_ha + T_TYPE + PLOT)
# f.nest <- formula(FATE ~ HEIGHT + DIST + ORIENT + OrdDate + OrdDate^2)
# f.ect <- formula(FATE ~ NEST + STUDY + YEAR)

#Adrian, Marty, and Max
str(ACFL)
#ACFL$FATE=as(ACFL$FATE)
ACFL.AMM <- bic.glm(FATE ~ OrdDate + OrdDate^2 + NestHeight + TrailDist + WBH + GPH + RoadDist + EdgeDist + RoadDensity + as.factor(Year),
                    glm.family = binomial(logexp(exposure = as.numeric(ACFL$EXP))), data = ACFL, na.rm = TRUE)
imageplot.bma(ACFL.AMM)

ACFL.AMM.glm <- glm(FATE ~ OrdDate + OrdDate^2 + NestHeight + TrailDist + WBH
                    + GPH + RoadDist + EdgeDist + RoadDensity + as.factor(Year), 
                    family = binomial(logexp(exposure = as.numeric(ACFL$EXP))), data = ACFL)
summary(ACFL.AMM.glm)

confint(ACFL.AMM.glm)#logLikelihood
confint.default(ACFL.AMM.glm)#SE
#pdf("ACFL_diagnostics.pdf")
plot(ACFL.AMM.glm)
#dev.off()

plot(predict(ACFL.AMM.glm), residuals(ACFL.AMM.glm))
abline(h=0,lty=2,col="grey")
plot(predict(ACFL.AMM.glm), residuals(ACFL.AMM.glm), col=c("blue","red")[1+ACFL$FATE])
abline(h=0,lty=2,col="grey")
lines(lowess(predict(ACFL.AMM.glm),residuals(ACFL.AMM.glm)),col="black",lwd=2)

#Diagnostics
roc(FATE ~ OrdDate + OrdDate^2 + NestHeight + TrailDist, data = ACFL3)
hoslem.test(m$y, fitted(m))
hoslem.test(ACFL$FATE, fitted(ACFL.AMM.glm))
dat <- ACFL

dat <- ACFL
pdf("ACFL_diagnostics.pdf")
for (i in 10:17){
  #png(file=paste('ACFL_plots_cont_',names(dat)[i],'.png', sep=''), height=8, width=8, res=500, units='in')
  par(mfrow=c(2,2))    
  #plot(histogram(~dat[,i]|resp, data=dat, xlab=paste(names(dat)[i]), main="Boulder", col="orange", pch=16))       
  #ttest=t.test(dat[,i]~dat$resp)
  wilcox=wilcox.test(dat[,i],dat$FATE)
  pvalue=wilcox$p.value
  boxplot(dat[,i]~dat$FATE, xlab="FATE", main=paste("ACFL",names(dat)[i]," Wilcox p-value=",round(pvalue,3)), col=rainbow(10))
  #bwplot( dat[,i]~dat$resp, main=paste("Boulder",names(dat)[i]," p-value=",round(pvalue,3)), col=rainbow(10),data=dat)   
  shapiro=shapiro.test(dat[,i])
  p.value=round(shapiro$p.value, 3)
  qqnorm(dat[,i],ylab=paste("Sample Quantiles, ",names(dat)[i]), main=paste("ACFL,Shapiro test, p-value=",round(p.value,3)), col=2)
  qqline(dat[,i], lty=5, col=3)
  g=glm(dat$FATE~dat[,i],family=binomial(logexp(exposure = as.numeric(dat$EXP))),dat) 
  plot(jitter(dat[,i]),dat$FATE,xlab=paste(names(dat)[i]),ylab="Probability of Nest Survival",main="ACFL",)
  points(dat[,i],fitted(g),pch=16, col="orange")
  logi.hist.plot(dat[,i],dat$FATE,boxp=FALSE,resp="hist",col="gray", main=paste("ACFL",names(dat)[i]))
  #dev.off()
}
dev.off()

dat$FATE <- as.factor(dat$FATE)#, levels=c('0','1'), labels=c('Survived','Failed')) #I'm labeling my numbers without actually changing the table.
pdf("ACFL_diagnostics_QQ.pdf")
for (i in 10:17){
 # png(file=paste('ACFL_QQplots_cont_',names(dat)[i],'.png', sep=''), height=8, width=8, res=500, units='in')
  plot <- qqmath(~dat[,i]|FATE, aspect="fill", data=dat,ylab=paste("Sample Quantiles, ",names(dat)[i]),main="ACFL",
         prepanel= prepanel.qqmathline,
         panel = function(x, ...) {
           panel.qqmathline(x, ...)
           panel.qqmath(x, ...)
         })
  print(plot)
  #dev.off()
}
dev.off()










#Standardize coefficients
hist(resid(ACFL.AMM.glm))
str(ACFL[10:17])
ACFL2 <- ACFL
ACFL2[10:17] <- scale(ACFL2[10:17], center = TRUE, scale = TRUE)

ACFL2.AMM.glm <- glm(FATE ~ OrdDate + OrdDate^2 + NestHeight + TrailDist + WBH + GPH + RoadDist + EdgeDist + RoadDensity + as.factor(Year), family = binomial(logexp(exposure = as.numeric(ACFL2$EXP))), data = ACFL2)
summary(ACFL2.AMM.glm)
plot(ACFL2.AMM.glm)
hist(resid(ACFL2.AMM.glm))

ACFL3 <- ACFL
for (i in 10:17){
  ACFL3[i] <- log(ACFL[i]+1)
}
ACFL3.AMM.glm <- glm(FATE ~ OrdDate + OrdDate^2 + NestHeight + TrailDist + WBH + GPH + RoadDist + EdgeDist + RoadDensity + as.factor(Year), family = binomial(logexp(exposure = as.numeric(ACFL3$EXP))), data = ACFL3)
summary(ACFL3.AMM.glm)
plot(ACFL3.AMM.glm)
hist(resid(ACFL3.AMM.glm))












#Max
ACFL.MAH <- bic.glm(FATE ~ OrdDate + OrdDate^2 + HEIGHT + DIST + WTR + WBH + GPH + TYPE + BA_ha + TPH + as.factor(Year),
                    glm.family = binomial(logexp(exposure = as.numeric(ACFL[ACFL$STUDY == "MAH",]$EXP))), data = ACFL[ACFL$STUDY == "MAH",], na.rm = TRUE)
imageplot.bma(ACFL.MAH)

#Marty
ACFL.MJP <- bic.glm(FATE ~ OrdDate + OrdDate^2 + HEIGHT + DIST + WTR + WBH + GPH + TYPE + BA_ha,
                    glm.family = binomial(logexp(exposure = as.numeric(ACFL[ACFL$STUDY == "MJP",]$EXP))), data = ACFL[ACFL$STUDY == "MJP",], na.rm = TRUE)
imageplot.bma(ACFL.MJP)

#Marty and Max
ACFL.MM <- bic.glm(FATE ~ OrdDate + OrdDate^2 + HEIGHT + DIST + WTR + WBH + GPH + TYPE + BA_ha + as.factor(Year),
                   glm.family = binomial(logexp(exposure = as.numeric(ACFL[ACFL$STUDY != "AL",]$EXP))), data = ACFL[ACFL$STUDY != "AL",], na.rm = TRUE)
imageplot.bma(ACFL.MM)

#Adrian
# ACFL.AL <- bic.glm(FATE ~ OrdDate + OrdDate^2 + HEIGHT*DIST + HEIGHT*WBH + TYPE*GPH,
#                    glm.family = binomial(logexp(exposure = as.numeric(ACFL[ACFL$STUDY == "AL",]$EXP))), data = ACFL[ACFL$STUDY == "AL",], na.rm = TRUE)
# imageplot.bma(ACFL.AL)





####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################