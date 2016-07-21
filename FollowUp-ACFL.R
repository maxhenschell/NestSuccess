library(popbio)
require("coefplot2")

####################################
#
# Full model
#
####################################

acfl.tmp <- read.csv("DSR-ACFL.csv", header = T)

acfl.fm <- glm(FATE ~ OrdDate + Year 
               + RoadDist + Core10km + Edge60m
               + ParaStat*TrailDist  + ParaStat*WBH + ParaStat*WTR
               + GPH*TrailDist + GPH*WTR + GPH*WBH + GPH*OrdDate,
            family = binomial(logexp(exposure = ACFL$EXP)), data = ACFL)


coefplot2(acfl.fm)
summary(acfl.fm)

acad.rm <- glm(FATE ~ OrdDate + ParaStat*WBH + ParaStat*WTR + GPH*WBH,
                family = binomial(logexp(exposure = ACFL$EXP)), data = ACFL)
summary(acad.rm)
coefplot2(acad.rm)

# summary(acad.fm)
mean(acad.rm$fitted)

# library(hier.part)
# 
# Fate <- ACFL$FATE
# Vars <- ACFL[c("OrdDar")]
# 
# hier.part(ACFL$FATE, )

####################################
#
# Null model
#
####################################

acad.n <- glm(FATE~ 1,
            family = binomial(logexp(exposure = ACFL$EXP)), data = ACFL)
summary(acad.n)
mean(acad.n$fitted)
#ACFL5 <- as.data.frame(ACFL2)
#acad.n <- glm(FATE~ 1,
 #             family = binomial(logexp(exposure = ACFL$EXP)), data = ACFL)
# str(ACFL2)
# attr(ACFL,"scaled:center")<-NULL 
# attr(ACFL,"scaled:scale")<-NULL 
# str(ACFL)
# mean(acad$fitted)
# OrdDate <- ACFL2$OrdDate
# plot(FATE~OrdDate, xlab = "OrdDate", ylab = "Fate", data = ACFL)
# curve(x, add = TRUE)
# plot(predict(acad))
# curve(predict(acad), add = TRUE)
# x <- function(x) {predict(acad, type = "response")}

predict(acad.n)


logi.hist.plot(as.vector(ACFL$OrdDate), as.vector(ACFL$FATE), boxp = FALSE, type = "hist", col ="gray")

# library(ggplot2)
# ggplot(ACFL, aes(x=OrdDate, y=FATE)) + geom_point() + stat_smooth(method = "glm", family =binomial, se = FALSE)
#                                                                   
#                                                                   
# plot(acad.fit)
#                                                                   
# exp(3.57)


####################################
#
# Controls v Trails
#
####################################

max.m <- ACFL[ACFL$Study %in% c('MAH'),]

MP <- ACFL[ACFL$Study %in% c('MJP'),]

summary(glm(FATE~as.factor(PlotType),  
    family = binomial(logexp(exposure = max.m$EXP)), data = max.m))

controls <- max.m[max.m$PlotType %in% 'I',]
constructed <- ACFL[ACFL$PlotType %in% 'C',]
roads <- ACFL[ACFL$PlotType %in% 'R',]

controls.n <- glm(FATE~ 1,
              family = binomial(logexp(exposure = controls$EXP)), data = controls)
summary(controls.n)
mean(controls.n$fitted)
unique(controls$NID)

controls.MP <- glm(FATE~ 1,
                   family = binomial(logexp(exposure = MP$EXP)), data = MP)
unique(MP$NID)
summary(controls.MP)
mean(controls.MP$fitted)


constructed.n <- glm(FATE~ 1,
                  family = binomial(logexp(exposure = constructed$EXP)), data = constructed)
summary(constructed.n)
mean(constructed.n$fitted)

roads.n <- glm(FATE~ 1,
                     family = binomial(logexp(exposure = roads$EXP)), data = roads)
summary(roads.n)
mean(roads.n$fitted)

ACFL.c <- ACFL[ACFL$Study %in% c('MAH', 'MP'),]
ACFL.c <- ACFL.c[ACFL.c$PlotType %in% 'I',]

c.glm <- glm(FATE~ 1,
             family = binomial(logexp(exposure = MP$EXP)), data = ACFL.c)
unique(ACFL.c$NID)
summary(c.glm)
mean(c.glm$fitted)



constructed.ac <- ACFL[ACFL$TrailType == 'Constructed',]
roads.ac <- ACFL[ACFL$TrailType %in% 'Old Road',]


ACFL.TT <- glm(FATE~as.factor(TrailType),  
            family = binomial(logexp(exposure = ACFL$EXP)), data = ACFL)
summary(ACFL.TT)
mean(ACFL.TT$fitted.values)

constructed.n <- glm(FATE~ 1,
                     family = binomial(logexp(exposure = constructed$EXP)), data = constructed)
summary(constructed.n)
mean(constructed.n$fitted)

roads.n <- glm(FATE~ 1,
               family = binomial(logexp(exposure = roads$EXP)), data = roads)
summary(roads.n)
mean(roads.n$fitted)

#Ordinal Date

max.ac <- ACFL[ACFL$Study %in% 'MAH',]
controls.ac <- max.ac[max.ac$PlotType %in% 'I',]
constructed.ac <- max.ac[max.ac$PlotType %in% 'C',]
roads.ac <- max.ac[max.ac$PlotType %in% 'R',]
str(max.ac)
max.plots.ac <- glm(FATE~ PlotType,family = binomial(logexp(exposure = max.ac$EXP)), data = max.ac)
summary(max.plots.ac)


controls.OD <- glm(FATE~ OrdDate, family = binomial(logexp(exposure = controls.ac$EXP)), data = controls.ac)
summary(controls.OD)
mean(controls.OD$fitted)
logi.hist.plot(as.vector(controls$OrdDate), as.vector(controls$FATE), boxp = FALSE, type = "hist", col ="gray")

OD.predict <- data.frame( OrdDate = seq(min(ACFL$OrdDate),max(ACFL$OrdDate), 0.01))
controls.pr <- predict(controls.OD, type = 'link', se = TRUE)
controls.resp <- predict(controls.OD, type = "response", se = TRUE)


require(ggplot2)

binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

ggplot(controls.ac, aes(x = OrdDate, y = FATE)) + geom_point() +
  binomial_smooth(se = FALSE) + coord_cartesian(xlim = c(-2, 2))

ggplot(constructed.ac, aes(x = OrdDate, y = FATE)) + geom_point() +
  binomial_smooth(se = FALSE)

ggplot(roads.ac, aes(x = OrdDate, y = FATE)) + geom_point() +
  binomial_smooth(se = FALSE)

plot(controls.ac$OrdDate, controls.resp$fit, type = "l", ylim = c(0,1))

lines(OD.predict, controls.resp)
plot(controls.ac$OrdDate, controls.ac$FATE)
curve(predict(controls.OD, data.frame(OrdDate = x), type = "response"), add = TRUE)

newdata.control <- cbind(OD.predict, controls.resp)



constructed.ac <- data.frame(constructed.ac)
constructed.OD <- glm(FATE ~ as.vector(OrdDate), family = binomial(logexp(exposure = constructed.ac$EXP)), data = constructed.ac)
summary(constructed.OD)
mean(constructed.OD$fitted)

OD.predict <- seq(min(constructed.ac$OrdDate), max(constructed.ac$OrdDate), 0.01)
fate.constructed <- predict(constructed.OD, newdata = data.frame(OrdDate = OD.predict), type="response")

plot(predict(constructed.OD, newdata = OD.predict), type="l")

dev.off()
plot(constructed.ac$OrdDate, constructed.ac$FATE, type = "n")
curve(OD.predict, fate.constructed)


logi.hist.plot(as.vector(constructed$OrdDate), as.vector(constructed$FATE), boxp = FALSE, type = "hist", col ="gray")

roads.OD <- glm(FATE~ OrdDate, family = binomial(logexp(exposure = roads.ac$EXP)), data = roads.ac)
summary(roads.OD)
mean(roads.OD$fitted)




pdf("ACFL_OrdDateByPlotType.pdf", height = 10, width = 8)
par(mfrow=c(3,1)) 
logi.hist.plot(as.vector(controls$OrdDate), as.vector(controls$FATE), boxp = FALSE, type = "hist", col ="gray", main = "Interior", ylabel = "Fate")
logi.hist.plot(as.vector(constructed$OrdDate), as.vector(constructed$FATE), boxp = FALSE, type = "hist", col ="gray", main = "Constructed Trails", ylabel = "Fate")
logi.hist.plot(as.vector(roads$OrdDate), as.vector(roads$FATE), boxp = FALSE, type = "hist", col ="gray", main = "Old Roads", xlab = "Ordinal Date", ylabel = "Fate")
dev.off()

#init date is 145 (5/23)


####################################
#
# Controls v Trails
#
####################################

Use <- read.csv("UseTime.csv", header = T, na.string = "NA")

pdf(file = paste(D.home,"UsePlots.pdf", sep = ""), width = 8, height = 10)
par(mfrow = c(2,1))
plot(Busy~J_Date, data = Use, ylab = "Groups per hour", xlab = " ", pch = 1, ylim = c(0,15))
#points(Use$Not.Busy, pch = 2)

plot(Not.Busy~J_Date, data = Use, ylab = "Groups per hour", xlab = "Ordinal Date", pch = 2, ylim = c(0,1))
dev.off()

points()


##
# controls.OD <- glm(FATE~ OrdDate,
#                    family = binomial(logexp(exposure = controls$EXP)), data = controls)
# summary(controls.OD)
# mean(controls.OD$fitted)
# logi.hist.plot(as.vector(controls$OrdDate), as.vector(controls$FATE), boxp = FALSE, type = "hist", col ="gray")
# 
# 
# constructed.OD <- glm(FATE~ OrdDate,
#                       family = binomial(logexp(exposure = constructed$EXP)), data = constructed)
# summary(constructed.OD)
# mean(constructed.OD$fitted)
# logi.hist.plot(as.vector(constructed$OrdDate), as.vector(constructed$FATE), boxp = FALSE, type = "hist", col ="gray")
# 
# 
# roads.OD <- glm(FATE~ OrdDate,
#                 family = binomial(logexp(exposure = roads$EXP)), data = roads)
# summary(roads.OD)
# mean(roads.OD$fitted)