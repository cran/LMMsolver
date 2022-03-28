## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.dim = c(7, 4)
)

## ----setup--------------------------------------------------------------------
library(LMMsolver)
library(ggplot2)

## ----oatsdata-----------------------------------------------------------------
## Load data.
data(john.alpha, package = "agridat")
head(john.alpha)

## ----Pspline1D----------------------------------------------------------------
## Fit mixed model with fixed and spline part.
obj1 <- LMMsolve(fixed = yield ~ rep + gen,
                 spline = ~spl1D(x = plot, nseg = 100),
                 data = john.alpha)

## ----deviance-----------------------------------------------------------------
round(deviance(obj1), 2)

## ----summaryED----------------------------------------------------------------
summary(obj1)

## ----summaryED2, echo=FALSE---------------------------------------------------
tbl <- summary(obj1)
EDs <- round(tbl$Effective[5], 1)

## ----coeff_example------------------------------------------------------------
genEff <- coef(obj1)$gen
head(genEff, 4)

## ----ObtainSmoothTrend1-------------------------------------------------------
## Extract smooth trend from mixed model.
plotDat1 <- obtainSmoothTrend(obj1, grid = 1000, includeIntercept = TRUE)
head(plotDat1)

## ----plot1D-------------------------------------------------------------------
## Plot smooth trend.
ggplot(data = plotDat1, aes(x = plot, y = ypred)) +
  geom_line(color = "red", size = 1) +
  labs(title = "Smooth spatial trend oats data", x = "plotnr", y = "yield") +
  theme(panel.grid = element_blank())

## ----define LVmodel-----------------------------------------------------------
## Add plot as factor.
john.alpha$plotF <- as.factor(john.alpha$plot)
## Define the precision matrix, see eqn (A2) in Boer et al (2020).
N <- nrow(john.alpha)
cN <- c(1 / sqrt(N - 1), rep(0, N - 2), 1 / sqrt(N - 1))
D <- diff(diag(N), diff = 1)
Delta <- 0.5 * crossprod(D)
LVinv <- 0.5 * (2 * Delta + cN %*% t(cN))
## Add LVinv to list, with name corresponding to random term.
lGinv <- list(plotF = LVinv)

## ----modelLV------------------------------------------------------------------
## Fit mixed model with first degree B-splines and first order differences.
obj2 <- LMMsolve(fixed = yield ~ rep + gen,
                 random = ~plotF, 
                 ginverse = lGinv, 
                 data = john.alpha)

## ----summaryVAR---------------------------------------------------------------
summary(obj2, which = "variances")

## ----APSIMdat-----------------------------------------------------------------
data(APSIMdat)
head(APSIMdat)

## ----APSIMmodel---------------------------------------------------------------
obj2 <- LMMsolve(biomass ~ 1,
                 spline = ~spl1D(x = das, nseg = 200), 
                 data = APSIMdat)

## ----APSIMmodelSummary--------------------------------------------------------
summary(obj2)

## ----APSIMplot----------------------------------------------------------------
plotDat2 <- obtainSmoothTrend(obj2, grid = 1000, includeIntercept = TRUE)

ggplot(data = APSIMdat, aes(x = das, y = biomass)) +
  geom_point(size = 1.2) +
  geom_line(data = plotDat2, aes(y = ypred), color = "red", size = 1) +
  geom_line(data = plotDat2, aes(y = ypred-2*se), col='blue', size=1) +
  geom_line(data = plotDat2, aes(y = ypred+2*se), col='blue', size=1) +
  labs(title = "APSIM biomass as function of time", 
       x = "days after sowing", y = "biomass (kg)") +
  theme(panel.grid = element_blank())

## ----APSIMDeriv---------------------------------------------------------------
plotDatDt <- obtainSmoothTrend(obj2, grid = 1000, deriv = 1)

ggplot(data = plotDatDt, aes(x = das, y = ypred)) +
  geom_line(color = "red", size = 1) +
  labs(title = "APSIM growth rate as function of time", 
       x = "days after sowing", y = "growth rate (kg/day)") +
  theme(panel.grid = element_blank())

## ----USprecip data------------------------------------------------------------
## Get precipitation data from spam
data(USprecip, package = "spam")

## Only use observed data
USprecip <- as.data.frame(USprecip)
USprecip <- USprecip[USprecip$infill == 1, ]

## ----runobj3------------------------------------------------------------------
obj3 <- LMMsolve(fixed = anomaly ~ 1,
                 spline = ~spl2D(x1 = lon, x2 = lat, nseg = c(41, 41)),
                 data = USprecip)

## ----ED_USprecip--------------------------------------------------------------
summary(obj3)

## ----Plot_USprecip------------------------------------------------------------
plotDat3 <- obtainSmoothTrend(obj3, grid = c(200, 300), includeIntercept = TRUE)
usa = maps::map("usa", regions = "main", plot = FALSE)
v <- sp::point.in.polygon(plotDat3$lon, plotDat3$lat, usa$x, usa$y)
plotDat3 <- plotDat3[v == 1, ]

ggplot(plotDat3, aes(x = lon, y = lat, fill = ypred)) +
  geom_tile(show.legend = TRUE) +
  scale_fill_gradientn(colors = topo.colors(100))+
  labs(title = "Precipitation (anomaly)", x = "Longitude", y = "Latitude") +
  coord_fixed() +
  theme(panel.grid = element_blank())

## ----newdataARG---------------------------------------------------------------
## Predictions for new data, using city coordinates from maps package.
data(us.cities, package = "maps")
## Column names have to match column names used for fitting the model.
colnames(us.cities)[colnames(us.cities) == "long"] <- "lon"
## Select columns name, lat and lon
us.cities <- us.cities[, c(1,4,5)]
head(us.cities)

pred3 <- obtainSmoothTrend(obj3, newdata = us.cities, includeIntercept = TRUE)
head(pred3)

## ----multipop-----------------------------------------------------------------
## Load data for multiparental population.
data(multipop)
head(multipop)

## ----residualARG--------------------------------------------------------------
## Fit null model.
obj4 <- LMMsolve(fixed = pheno ~ cross, 
                 residual = ~cross, 
                 data = multipop)
dev4 <- deviance(obj4)

## ----groupOPTION--------------------------------------------------------------
## Fit alternative model - include QTL with probabilities defined in columns 3:5 
lGrp <- list(QTL = 3:5)
obj5 <- LMMsolve(fixed = pheno ~ cross, 
                 group = lGrp,
                 random = ~grp(QTL),
                 residual = ~cross,
                 data = multipop) 
dev5 <- deviance(obj5)

## ----approxPvalue-------------------------------------------------------------
## Deviance difference between null and alternative model.
dev <- dev4 - dev5
## Calculate approximate p-value. 
minlog10p <- -log10(0.5 * pchisq(dev, 1, lower.tail = FALSE))
round(minlog10p, 2)

## ----QTLeffects---------------------------------------------------------------
coef(obj5)$QTL

