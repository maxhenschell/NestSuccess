####################################
#
# BMA for ACFL
# in the Baraboo Hills
# Created by: Max Henschell
# Modified: 17 September 2014 
#
####################################

lapply(packages.list, require, character.only=T)

####################################
#
# Correlation
#
####################################

#pdf(file = paste(D.home,'ACFL_CorrelationPlots.pdf', sep = ''), width = 10, height = 8)
pairs(ACFL[,c('FATE','ParaStat', 'OrdDate','TrailDist', 'WTR', 'WBH', 'GPH', 
              'Core10km', 'Edge60m')], upper.panel = panel.cor, main = 'MAH + MJP')
#dev.off()

####################################
#
# Marginal plots
#
####################################

ACFL.M <- ACFL1
ACFL.M <- as.data.frame(ACFL)
#str(ACFL.M[,c(9:16)])
attr(ACFL.M,'scaled:center')<-NULL 
attr(ACFL.M,'scaled:scale')<-NULL 
#str(ACFL.M) 

mdat <- melt(ACFL.M, id.var = c('NID','FATE','EXP'), measure.vars = c('OrdDate', 'TrailDist', 'WTR', 'WBH', 'GPH', 
                                                                      'RoadDist', 'EdgeDist', 'Core10km', 'Edge60m'))

####################################
#
# BMA
#
####################################

#Centered and scaled
# attr(ACFL,'scaled:center')<-NULL 
# attr(ACFL,'scaled:scale')<-NULL 

ACFL.mah <- ACFL[ACFL$Study %in% "MAH",]
ACFL.mah$Study <- factor(ACFL.mah$Study)
ACFL.mah$Year <- factor(ACFL.mah$Year)
str(ACFL.mah)

ACFL.mjp <- ACFL[ACFL$Study %in% "MJP",]
ACFL.mjp$Study <- factor(ACFL.mjp$Study)
ACFL.mjp$Year <- factor(ACFL.mjp$Year)
str(ACFL.mjp)

##From meeting with Murray 7/13

ACFL.MAH.glm <- glm(FATE ~ Year + Core10km + Edge60m
    + ParaStat*TrailDist  
    + ParaStat*WBH + GPH*WBH
    + ParaStat*WTR + GPH*WTR 
    + GPH*TrailDist  + GPH*OrdDate,
    family = binomial(logexp(exposure = ACFL.mah$EXP)), data = ACFL.mah)

summary(ACFL.MAH.glm)
require(MuMIn)
options(na.action = "na.fail")
ACFL.MAH.dredge <- dredge(ACFL.MAH.glm, rank = AIC)
subset(ACFL.MAH.dredge, subset = delta < 4)

summary(glm(FATE ~ OrdDate + WTR,
    family = binomial(logexp(exposure = ACFL.mah$EXP)), data = ACFL.mah))
summary(glm(FATE ~ OrdDate + WBH,
            family = binomial(logexp(exposure = ACFL.mah$EXP)), data = ACFL.mah))
summary(glm(FATE ~ OrdDate,
            family = binomial(logexp(exposure = ACFL.mah$EXP)), data = ACFL.mah))

#Just Marty
ACFL.MJP.glm <- glm(FATE ~ Core10km + Edge60m
                    + ParaStat*TrailDist  
                    + ParaStat*WBH + GPH*WBH
                    + ParaStat*WTR + GPH*WTR 
                    + GPH*TrailDist  + GPH*OrdDate,
                    family = binomial(logexp(exposure = ACFL.mjp$EXP)), data = ACFL.mjp)

summary(ACFL.MJP.glm)
require(MuMIn)
options(na.action = "na.fail")
ACFL.MJP.dredge <- dredge(ACFL.MJP.glm, rank = BIC)
subset(ACFL.MJP.dredge, delta < 4)

ACFL.MJP.glm2 <- glm(FATE ~ Core10km + Edge60m + OrdDate
                    + ParaStat*TrailDist  
                    + ParaStat*WBH + #GPH*WBH
                    + ParaStat*WTR #+ GPH*WTR 
                   # + GPH*TrailDist  + GPH*OrdDate,
                    ,family = binomial(logexp(exposure = ACFL.mjp$EXP)), data = ACFL.mjp)

summary(ACFL.MJP.glm2)
require(MuMIn)
options(na.action = "na.fail")
ACFL.MJP.dredge2 <- dredge(ACFL.MJP.glm2, rank = BIC)
subset(ACFL.MJP.dredge2, delta < 4)


summary(bic.glm(FATE ~ Core10km + Edge60m + OrdDate
                     + ParaStat*TrailDist  
                     + ParaStat*WBH + GPH*WBH
                       + ParaStat*WTR + GPH*WTR 
                      + GPH*TrailDist  + GPH*OrdDate,
                     ,glm.family = binomial(logexp(exposure = ACFL.mjp$EXP)), data = ACFL.mjp))




## Just Max
## GPH and WBH
summary(bic.glm(FATE ~ Year + Core10km + Edge60m
        + ParaStat*TrailDist  
        + ParaStat*WBH + GPH*WBH
       + ParaStat*WTR + GPH*WTR 
        + GPH*TrailDist  + GPH*OrdDate,
        glm.family = binomial(logexp(exposure = ACFL.mah$EXP)), data = ACFL.mah))

## Just Max
## GPH and WTR
summary(bic.glm(FATE ~ Year + Core10km + Edge60m
                + ParaStat*TrailDist  
               # + ParaStat*WBH + GPH*WBH
                + ParaStat*WTR + GPH*WTR 
                + GPH*TrailDist  + GPH*OrdDate,
                glm.family = binomial(logexp(exposure = ACFL.mah$EXP)), data = ACFL.mah))


## Matry and Max
## GPH and WTR
summary(bic.glm(FATE ~ Year + Core10km + Edge60m
                + ParaStat*TrailDist  
                # + ParaStat*WBH + GPH*WBH
                + ParaStat*WTR + GPH*WTR 
                + GPH*TrailDist  + GPH*OrdDate,
                    glm.family = binomial(logexp(exposure = ACFL$EXP)), data = ACFL))



# imageplot.bma(ACFL.bic, color = c('green', 'red', 'white'))
# summary(ACFL.bic, n.models = 6)



summary(ACFL.bic.m, n.models = 6)


##Combine variable P(ne0) > .5 with study interactions

summary(bic.glm(FATE ~
                ParaStat*WTR*Study + GPH*WTR*Study
                + OrdDate*GPH*Study,
                glm.family = binomial(logexp(exposure = ACFL$EXP)), data = ACFL))


summary(bic.glm(FATE ~#Year +
                Core10km*Study + Edge60m*Study
                + ParaStat*Study 
                + TrailDist*Study  
                + OrdDate*Study
                + WBH*Study + GPH*WBH*Study
                + GPH*Study
               # + WTR*Study + GPH*WTR*Study 
                +ParaStat*TrailDist*Study
                +ParaStat*WBH*Study
                +GPH*OrdDate*Study
                +GPH*TrailDist*Study,
                glm.family = binomial(logexp(exposure = ACFL$EXP)), data = ACFL))

####################################
#
# Save results
#
####################################

pdf(file = paste(D.home,'ACFL_bma_plots_',Sys.Date(),'.pdf', sep = ''), width = 10, height = 8)

#plot(ParaStat~TrailDist, data = ACFL.nl, ylab = 'Parasitized', xlab = 'Distance from trail')

pairs(ACFL[,c('FATE','ParaStat', 'OrdDate','TrailDist', 'WTR', 'WBH', 'GPH', 
              'CoreArea', 'EdgeArea')], upper.panel = panel.cor, main = 'MAH + MJP')
#marginal plots
print(
  ggplot(mdat,aes(x=value,y=FATE))+
    geom_point(alpha=0.5,aes(size=EXP))+geom_smooth(method='loess')+
    facet_wrap(~variable,scale='free')+
    coord_cartesian(ylim=c(-0.05,1.05))+xlab('')+ylab('Survival')
)
imageplot.bma(ACFL.bic, color = c('green', 'red', 'white'))

plot(ACFL.bic)
dev.off()
dev.off()

write.table(summary(ACFL.bic, n.models = 6), file = paste(D.home,'ACFL_bma_summary_',Sys.Date(),'.csv', sep = ''), sep = ',')

write.table(summary(ACFL.bic.m, n.models = 6), file = paste(D.home,'ACFL_bma.m_summary_',Sys.Date(),'.csv', sep = ''), sep = ',')
#####################   END   #####################
#####################   END   #####################
#####################   END   #####################
#####################   END   #####################
#####################   END   #####################
#####################   END   #####################



##Old Code
# ####################################
# #
# # TrailDist Histogram (effort)
# #
# ####################################
# 
# #pdf(file = paste(D.home,'ACFL_hist.pdf', sep = ''), width = 10, height = 8)
# hist(ACFL.nl$TrailDist, main = 'ACFL nests', ylab = '# nests', xlab = 'Distance from trail')
# dev.off()
# 
# ####################################
# #
# # Brood Parasitism
# #
# ####################################
# 
# #pdf(file = paste(D.home, 'ACFL_ParaStat.pdf', sep = ''), width = 10, height = 8)
# plot(ParaStat~TrailDist, data = ACFL.nl, ylab = 'Parasitized', xlab = 'Distance from trail')