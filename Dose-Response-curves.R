# https://rstats4ag.org/dose-response-curves.html#which-dose-response-model
library(drc)
head(ryegrass,3)
op <- par(mfrow = c(1, 2), mar=c(3.2,3.2,2,.5), mgp=c(2,.7,0)) #make two plots in two columns 
plot(rootl ~ conc, data = ryegrass, main="Original Dose Scale")
plot(rootl ~ log(conc+.1), data = ryegrass, main="Logarithmic Dose Scale")

# drm(): Fitting dose-response models 
# LL.4: Log-logistic (ED50 as parameter) (4 parms)
#  and the Weibull model W1.4. Use getMeanFunctions for a full list.
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, 
                   fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(ryegrass.m1)

op <- par(mfrow = c(1, 2), mar=c(3.2,3.2,.5,.5), mgp=c(2,.7,0))
#Plot of regression and averages of observations within each concentration.
#aggregate(ryegrass[, 1], list(ryegrass[, 2]), mean)
plot(ryegrass.m1, broken=TRUE, bty="l",
     xlab="Concentration of Ferulic Acid", ylab="Length of Roots")
#Plot all observation
plot(ryegrass.m1, broken=TRUE, bty="l",
     xlab="Concentration of Ferulic Acid", ylab="Length of Roots",type="all")
#Test for lack of fit
modelFit(ryegrass.m1)
#Assumption: Correct regression model, Variance Homogeneity,Normally distributed measurement errors, Mutually independent measurement error ε
#Graphical analysis of residuals
op <- par(mfrow = c(1, 2), mar=c(3.2,3.2,2,.5), mgp=c(2,.7,0)) #put two graphs together
plot(residuals(ryegrass.m1) ~ fitted(ryegrass.m1), main="Residuals vs Fitted")
abline(h=0)
qqnorm(residuals(ryegrass.m1))
qqline(residuals(ryegrass.m1))

#liner with R^2, root~conc
op <- par(mfrow = c(1, 2), mar=c(3.2,3.2,0,.5), mgp=c(2,.7,0)) #put two graphs together
plot(rootl ~ conc, data = ryegrass, bty="l")
linear.m1 <- lm(rootl ~ conc, data = ryegrass)
summary(linear.m1)
abline(linear.m1)
plot(linear.m1, which=1, bty="l")

# which model better, small better.
AIC(ryegrass.m1, linear.m1)

#find ED 10, 50, 90
ED(ryegrass.m1, c(10, 50, 90), interval = "delta")

# 11.4 More Dose-Response Curves
head(S.alba,n=2)
## Fitting a log-logistic model with
S.alba.m1 <- drm(DryMatter ~ Dose, Herbicide, data=S.alba, fct = LL.4())
modelFit(S.alba.m1)
summary(S.alba.m1)

##  common lower and upper limits
S.alba.m2 <- drm(DryMatter ~ Dose, Herbicide, data=S.alba, fct = LL.4(),
                 pmodels=data.frame(Herbicide, 1, 1, Herbicide)) 
summary(S.alba.m2)
par(mfrow=c(1,1), mar=c(3.2,3.2,.5,.5), mgp=c(2,.7,0))
plot(S.alba.m1, broken=TRUE, bty="l")
plot(S.alba.m2, col="red", broken=TRUE, add=TRUE, legend=FALSE)
# Test for lack of fit not signifiancant
anova(S.alba.m1, S.alba.m2) 
#A test for lack of fit was non-significant and, 
#therefore, we presume that the curves have similar upper an lower limits


##  common lower and upper limits and slope
S.alba.m3 <- drm(DryMatter~Dose, Herbicide, data=S.alba, fct = LL.4(),
                 pmodels=data.frame(1, 1, 1, Herbicide)) 
anova(S.alba.m2, S.alba.m3)
#the test for lack of fit now is significant, 
#so the assumption of similar curves except for the ED50 is not supported
par(mfrow=c(1,1), mar=c(3.2,3.2,.5,.5), mgp=c(2,.7,0))
plot(S.alba.m3, broken=TRUE, bty="l")

#different relative potencies:  ED50
EDcomp(S.alba.m2, c(50, 50), interval="delta", reverse=TRUE)

#relationship between changes of the relative potency and the predicted responses.
par(mfrow=c(1,1), mar=c(3.2,3.2,.5,.5), mgp=c(2,.7,0))
relpot(S.alba.m2, interval = "delta", bty="l")

#11.5 When Upper and Lower Limits are not similar
op <- par(mfrow = c(1, 1), mar=c(3.2,3.2,.5,.5), mgp=c(2,.7,0))
TwoCurves<-read.csv("/Users/Hui/Hui2020/Julie/Drug_synergy/twoCurves.csv")
head(TwoCurves,3)
Bispyribac.sodium.m1 <- drm(Dry.weight.per.pot ~ Dose, Accesion,
                            data = TwoCurves, fct=LL.4())
modelFit(Bispyribac.sodium.m1)
summary(Bispyribac.sodium.m1)
plot(Bispyribac.sodium.m1, broken=TRUE, bp=.1, ylab="Dry Matter(g/pot)",
     xlab="Bispyribac sodium (g/ha)", bty="l")

# drawing lines on the plot
ED.SAM.E.8 <- 0.24591 + (2.24201 - 0.24591) * 0.5
ED.TEK.E.10 <- -0.05308 + (1.82755 +- 0.05308) * 0.5
arrows(0.1, ED.SAM.E.8, 21.7, ED.SAM.E.8, code=0)
arrows(21.7, ED.SAM.E.8, 21.7, 0, code=0)
arrows(0.1, ED.TEK.E.10, 2.67, ED.TEK.E.10, code=0, lty=2)
arrows(2.67, ED.TEK.E.10, 2.67, 0, code=0, lty=2)

Bispyribac.sodium.m2 <- drm(Dry.weight.per.pot ~ Dose, Accesion,
                            data = TwoCurves, fct=LL.4(),
                            lowerl=c(NA,0,NA,NA))
summary(Bispyribac.sodium.m2)
EDcomp(Bispyribac.sodium.m2, c(50, 50),interval="delta")
ED(Bispyribac.sodium.m2, 1.24, type="absolute", interval="delta")
par(mfrow = c(1, 1), mar=c(3.2,3.2,.5,.5), mgp=c(2,.7,0))
plot(Bispyribac.sodium.m2, broken=TRUE, bp=.1, bty="l",
     ylab="g Dry Matter/pot",
     xlab="Dose (g Bispyribac sodium/ha)")
arrows(  0.1, 1.24, 21.79, 1.24, code=0, lty=1, col="red")
arrows(21.79, 1.24, 21.79,    0, code=0, lty=1, col="red")
arrows(  1.31, 1.24,   1.31,    0, code=0, lty=2, col="red")
SAM.E.8.C<-coef(Bispyribac.sodium.m1)[3]
TEK.E.10.C<-coef(Bispyribac.sodium.m1)[4]
SAM.E.8.d<-coef(Bispyribac.sodium.m1)[5]
TEK.E.10.d<-coef(Bispyribac.sodium.m1)[6]
Rel.Y<-with(TwoCurves,
            ifelse(Accesion=="SAM-E-8", 
                   ((Dry.weight.per.pot-SAM.E.8.C)/(SAM.E.8.d-SAM.E.8.C)),
                   ((Dry.weight.per.pot-TEK.E.10.C)/(TEK.E.10.d-TEK.E.10.C))))

Scaled.m1<-drm(Rel.Y ~ Dose, Accesion, pmodels=data.frame(Accesion, 1, Accesion),
               data=TwoCurves, fct=LL.3())
par(mfrow = c(1, 1), mar=c(3.2,3.2,.5,.5), mgp=c(2,.7,0))
plot(Scaled.m1, broken=TRUE, bp=.1, bty="l", 
     ylab="Relative Dry Matter", ylim=c(0,1.2),
     xlab="Bispyribac-sodium dose (g/ha)")
arrows(  0.1, .5, 21.7,.5, code=0, lty=1, col="red")
arrows(21.79, 0.5, 21.79,    0, code=0, lty=1, col="red")
arrows(  2.66, .5,   2.66,    0, code=0, lty=2, col="red")
summary(Scaled.m1)
EDcomp(Scaled.m1, c(50, 50),interval="delta")

# Remedy for heterogeneous variance.
par(mfrow = c(1, 1), mar=c(3.2,3.2,.5,.5), mgp=c(2,.7,0))
ryegrass.BX <- boxcox(ryegrass.m1,method="anova")
summary(ryegrass.BX)
ED(ryegrass.BX,50,interval="delta")
#Graphical analysis of residuals
op <- par(mfrow = c(1, 2), mar=c(3.2,3.2,2,.5), mgp=c(2,.7,0)) #put two graphs together
plot(residuals(ryegrass.BX) ~ fitted(ryegrass.BX), main="Residuals vs Fitted")
abline(h=0)
qqnorm(residuals(ryegrass.BX))
qqline(residuals(ryegrass.BX))

##############  DND41  ###########################
par(mfrow = c(2, 1), mar=c(3.2,3.2,.5,.5), mgp=c(2,.7,0))
#A1: GSK4716
G<- drm(A1.vib ~ A1.conc, fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
ED(G, c(30, 50, 70), interval = "delta")
modelFit(G)
summary(G)
plot(G, broken=TRUE, bp=1, ylab="Response",
     xlab="GSK", bty="l", type = "all", xlim = c(0, 20000 ), ylim = c(-0.1, 1))

#A2:Dexamethasone
D<- drm(A2.vib ~ A2.conc, fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
ED(D, c(30, 50, 70), interval = "delta")
modelFit(D)
summary(D)
plot(D, broken=TRUE, bp=1, ylab="Response",
     xlab="Dex", bty="l", type = "all", xlim = c(0, 1500 ), ylim = c(-0.1, 1))
