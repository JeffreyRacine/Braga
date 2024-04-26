## ----global_options-----------------------------------------------------------
#| include: false
library(ggplot2)
library(kableExtra)
library(fda)
library(gamair)
library(np)
library(latex2exp)
source("R_scripts_data/func_lib.R")


## ----sampleelements-----------------------------------------------------------
#| label: fig-sampleelements
#| fig-cap: "Classical Versus Functional Sample Elements ($N=25$)"
## Simple function to demonstrate classical versus FDA sample elements

set.seed(42)

## Set parameters here

## Regularity parameters for generating data

h_first = 0.5
h_second = 0.5
h_slope = 0
change_point_t = 0.5

## Parameters for generating data from a Brownian motion

sigma = 0.5
tau = 1
L = 1

## Set size of batch (N), (average) number of observations for each online curve
## (M), number of online curves, and grid of points on which to estimate curves

N = 25
M = 1

t0.grid = seq(0,1,length=100)

## Set the mean function

#mu.func <- function(x) { 1 + x }
#mu.func <- function(x) { sin(2*pi*x) }
mu.func <- function(t) {t+exp(-200*(t-0.5)^2)}

## Generate batch data

points_dist <- functional::Curry(runif, min = 0, max = 1)

points_list <- generates_equal_length_points(N = N, m = M,
                                             distribution = points_dist)

batch <- generate_curves(points_list,
                         hurst = hurst_logistic,
                         grid = t0.grid,
                         sigma = sigma,
                         tau = tau,
                         L = L)

curves_batch <- lapply(batch$curves,
                       function(x) list(t = x$observed$t,
                                        x = x$observed$x + mu.func(x$observed$t),
                                        grid_true = x$ideal$t,
                                        x_true = x$ideal$x + mu.func(x$ideal$t)))

x <- runif(N)
y <- mu.func(x) + rnorm(N,sd=sigma)
ylim.xy <- range(y)

ylim.curve <- range(curves_batch[[1]]$x_true)
if(N>1) {
  for(n in 2:N) {
    ylim.curve <- range(ylim.curve,curves_batch[[n]]$x_true)
  }
}

par(mfrow=c(1,2))
plot(x[1:N],y[1:N],pch=19,
     xlab="X",
     ylab="Y",
     ylim=ylim.xy,
     xlim=c(0,1),
     main="Classical Sample Elements",
     panel.first=grid(lty=1))

plot(t0.grid,curves_batch[[1]]$x_true,
     type="l",
     ylim=ylim.curve,
     xlim=c(0,1),
     main="FDA Sample Elements",
     xlab="t",
     ylab="Curve",
     panel.first=grid(lty=1))

if(N>1) {
  for(j in 2:N) {
    lines(t0.grid,curves_batch[[j]]$x_true,col=j,lty=j)
  }
}



## ----noisyfuncdata------------------------------------------------------------
## Generate data for plotting, plot in next code chunk
set.seed(42)

## Set parameters here

## Regularity parameters for generating data

h_first = 0.05
h_second = 0.05
h_slope = 0
change_point_t = 0.5

## Parameters for generating data from a multifractional Brownian motion

sigma = 0.0
tau = 0.0
L = 0.25

## Set size of batch (N), (average) number of observations for each online curve
## (M), number of online curves, and grid of points on which to estimate curves

N = 1
M = 100

t0.grid = seq(0,1,length=2000)

## Set the mean function

mu.func <- function(x) { 1 + x }
#mu.func <- function(x) { -cos(2*pi*x) }
#mu.func <- function(t) {t+exp(-200*(t-0.5)^2)}

## Generate batch data

points_dist <- functional::Curry(runif, min = 0, max = 1)

points_list <- generates_equal_length_points(N = N, m = M,
                                             distribution = points_dist)

batch <- generate_curves(points_list,
                         hurst = hurst_logistic,
                         grid = t0.grid,
                         sigma = sigma,
                         tau = tau,
                         L = L)

curves_irregular <- lapply(batch$curves,
                       function(x) list(t = x$observed$t,
                                        x = x$observed$x + mu.func(x$observed$t),
                                        grid_true = x$ideal$t,
                                        x_true = x$ideal$x + mu.func(x$ideal$t)))

ylim.irregular <- range(curves_irregular[[1]]$x_true,curves_irregular[[1]]$x)
if(N>1) {
  for(n in 2:N) {
    ylim.irregular <- range(ylim.irregular,curves_irregular[[n]]$x_true)
  }
}

sigma = 0.25
tau = 0.0
L = 0

points_dist <- functional::Curry(runif, min = 0, max = 1)

points_list <- generates_equal_length_points(N = N, m = M,
                                             distribution = points_dist)

batch <- generate_curves(points_list,
                         hurst = hurst_logistic,
                         grid = t0.grid,
                         sigma = sigma,
                         tau = tau,
                         L = L)

curves_regular <- lapply(batch$curves,
                       function(x) list(t = x$observed$t,
                                        x = x$observed$x + mu.func(x$observed$t),
                                        grid_true = x$ideal$t,
                                        x_true = x$ideal$x + mu.func(x$ideal$t)))

ylim.regular <- range(curves_regular[[1]]$x_true,curves_regular[[1]]$x)
if(N>1) {
  for(n in 2:N) {
    ylim.regular <- range(ylim.regular,curves_regular[[n]]$x_true)
  }
}


## ----nonnoisyfuncall----------------------------------------------------------
#| label: fig-nonnoisyfuncall
#| fig-cap: "Irregular Function, Data Measured Without Error"
par(mfrow=c(1,2))
plot(curves_irregular[[1]]$t,curves_irregular[[1]]$x,
     type="p",
     ylim=ylim.irregular,
     xlab=TeX("$t$"),
     main="Discrete Measurements for 1 Curve",
     ylab=TeX("$Y$"),
     pch=19,
     panel.first=grid(lty=1))
plot(t0.grid,curves_irregular[[1]]$x_true,
     type="l",
     ylim=ylim.irregular,
     xlab=TeX("$t$"),
     main="Irregular Function, Data Measured Without Error",
     col="grey",
     ylab=TeX("$Y$"),
     panel.first=grid(lty=1))
points(curves_irregular[[1]]$t,curves_irregular[[1]]$x,pch=19)


## ----noisyfuncall-------------------------------------------------------------
#| label: fig-noisyfuncall
#| fig-cap: "Regular Function, Noisy Data"
par(mfrow=c(1,2))
plot(curves_regular[[1]]$t,curves_regular[[1]]$x,
     type="p",
     ylim=ylim.regular,
     xlab=TeX("$t$"),
     main="Discrete Measurements for 1 Curve",
     ylab=TeX("$Y$"),
     pch=19,
     panel.first=grid(lty=1))
plot(t0.grid,curves_regular[[1]]$x_true,
     type="l",
     ylim=ylim.regular,
     xlab=TeX("$t$"),
     main="Smooth Function, Data Measured With Error",
     ylab=TeX("$Y$"),
     panel.first=grid(lty=1))
points(curves_regular[[1]]$t,curves_regular[[1]]$x,pch=19)


## ----mfbm---------------------------------------------------------------------
#| label: fig-mfbm
#| fig-cap: "Irregular Function, Varying Regularity, Noisy Data"
set.seed(42)
require(GA)
require(MASS)

## Set parameters here

## Regularity parameters for generating data

h_first = 0
h_second = 1
h_slope = 10
change_point_t = 0.5

## Parameters for generating data from a multifractional Brownian motion

sigma = 0.25
tau = 0
L = 1

## Set size of batch (N), (average) number of observations for each online curve
## (M), number of online curves, and grid of points on which to estimate curves

N = 1
M = 50

t0.grid = seq(0,1,length=1000)

## Set the mean function

mu.func <- function(x) { 1 + x }
#mu.func <- function(x) { sin(2*pi*x) }
#mu.func <- function(t) {t+exp(-200*(t-0.5)^2)}

## Generate batch data

points_dist <- functional::Curry(runif, min = 0, max = 1)

points_list <- generates_equal_length_points(N = N, m = M,
                                             distribution = points_dist)

batch <- generate_curves(points_list,
                         hurst = hurst_logistic,
                         grid = t0.grid,
                         sigma = sigma,
                         tau = tau,
                         L = L)

curve.Brownian <- lapply(batch$curves,
                         function(x) list(t = x$observed$t,
                                          x = x$observed$x + mu.func(x$observed$t),
                                          grid_true = x$ideal$t,
                                          x_true = x$ideal$x + mu.func(x$ideal$t)))

## Some plots and true values

true_H <- hurst_logistic(t0.grid)
C1 <- 3*L**2/((2*true_H+1)*(2*true_H+3))
C2 <- (3/5)*sigma**2
h.MSE <- (C2/(2*true_H*C1*M))**(1/(2*true_H+1))

par(mfrow=c(1,1))
plot(t0.grid,curve.Brownian[[1]]$x_true,
     type="l",
     ylim=range(curve.Brownian[[1]]$x_true,curve.Brownian[[1]]$x),
     main=if(N>1){paste("True Curves ($M$ = ",M,", $N$ = ",N,")",sep="")}else{NULL},
     sub=TeX("1 draw from $X_t$ discretely sampled with error ($Y^{(i)}(T_i)=X^{(i)}(T_i)+\\epsilon(T_i)$) at randomly spaced $T_i$"),
     xlab=TeX("$t$"),
     ylab=TeX("$X^{(i)}(t)$, $Y^{(i)}(t)$"),
     col="grey",
     panel.first=grid(lty=1))
points(curve.Brownian[[1]]$t,curve.Brownian[[1]]$x,
       col=1,
       pch=19)
rug(curve.Brownian[[1]]$t,lwd=2)
legend("topright",c(TeX("Sample Curve (random draw from $X_t$)"),TeX("Sample Data (curve measured with error at discrete points)")),
       lty=c(1,NA),
       pch=c(NA,19),
       pt.cex=c(NA,.5),
       col=c("grey",1),
       bty="n")


## ----growth-------------------------------------------------------------------
#| label: fig-growth
#| fig-cap: "Berkeley Growth Study Data"
par(mfrow=c(1,1))
with(growth, matplot(age, hgtf[, 1:NCOL(hgtf)], ylab="Height (cm)", xlab = "Age (years)",type="l",
sub="Heights of 54 girls from age 1 to 18 (ages not equally spaced)",panel.first=grid(lty=1)))
with(growth,rug(age))


## ----canWeather---------------------------------------------------------------
#| label: fig-canWeather
#| fig-cap: "Canadian Weather Study Data"
require(gamair)
data(canWeather)
reg <- unique(CanWeather$region)
place <- unique(CanWeather$place)
col <- 1:4;names(col) <- reg
par(mfrow=c(1,1))
for (k in 1:35) {
  if (k==1) plot(1:365,CanWeather$T[CanWeather$place==place[k]],
            ylim=range(CanWeather$T),type="l",
            sub="Data on daily mean temperature (C) throughout the year at 35 Canadian locations",
            col=col[CanWeather$region],
            xlab="Day of the year",
            ylab="Temperature (C)",
            panel.first=grid(lty=1)) else lines(1:365,CanWeather$T[CanWeather$place==place[k]],
            col=k)

}


## ----mfbm100------------------------------------------------------------------
set.seed(123)
require(GA)
require(MASS)

## Set parameters here

## Regularity parameters for generating data

h_first = 0.25
h_second = 0.75
h_slope = 10
change_point_t = 0.5

## Parameters for generating data from a Multifractional Brownian motion

sigma = 0.50
tau = 1
L = 1

## Set size of batch (N), (average) number of observations for each online curve
## (M), number of online curves, and grid of points on which to estimate curves

N = 100
M = 100

t0.grid = seq(0,1,length=50)

## Set the mean function

#mu.func <- function(x) { 1 + x }
mu.func <- function(x) { sin(2*pi*x) }
#mu.func <- function(t) {t+exp(-200*(t-0.5)^2)}

## Generate batch data

points_dist <- functional::Curry(runif, min = 0, max = 1)

points_list <- generates_equal_length_points(N = N, m = M,
                                             distribution = points_dist)

batch <- generate_curves(points_list,
                         hurst = hurst_logistic,
                         grid = t0.grid,
                         sigma = sigma,
                         tau = tau,
                         L = L)

curves.Brownian <- lapply(batch$curves,
                       function(x) list(t = x$observed$t,
                                        x = x$observed$x + mu.func(x$observed$t),
                                        grid_true = x$ideal$t,
                                        x_true = x$ideal$x + mu.func(x$ideal$t)))

## Some plots and true values

true_H <- hurst_logistic(t0.grid)
C1 <- 3*L**2/((2*true_H+1)*(2*true_H+3))
C2 <- (3/5)*sigma**2
h.MSE <- (C2/(2*true_H*C1*M))**(1/(2*true_H+1))

## Call the function on the batch, use preliminary starting values for mean
## function

mu.starting <- mu.pooled.func(curves.Brownian,t0.grid,method="np")

fda.la.Brownian <- fda.la(curves_batch=curves.Brownian,
                          t.grid=t0.grid,
                          mu.hat=mu.starting)

H.Brownian <- fda.la.Brownian$H.updated


## ----mfbmcurves---------------------------------------------------------------
#| label: fig-mfbmcurves
#| fig-cap: "Multifractional Brownian Motion - True Curves"
ylim <- range(curves.Brownian[[1]]$x_true)
if(N>1) {
  for(n in 2:N) {
    ylim <- range(ylim,curves.Brownian[[n]]$x_true)
  }
}
par(mfrow=c(1,1))
plot(t0.grid,curves.Brownian[[1]]$x_true,
     type="l",
     ylim=ylim,
     main=TeX(paste("True Curves ($M$ = ",M,", $N$ = ",N,")",sep="")),
     xlab=TeX("$t$"),
     ylab="True Curves",
     panel.first=grid(lty=1))
if(N>1) {
  for(n in 2:N) {
    lines(t0.grid,curves.Brownian[[n]]$x_true,col=n,lty=n)
  }
}


## ----mfbmdata-----------------------------------------------------------------
#| label: fig-mfbmdata
#| fig-cap: "Multifractional Brownian Motion - Data"
ylim <- range(curves.Brownian[[1]]$x)
if(N>1) {
  for(n in 2:N) {
    ylim <- range(ylim,curves.Brownian[[n]]$x)
  }
}
par(mfrow=c(1,1))
plot(curves.Brownian[[1]]$t,curves.Brownian[[1]]$x,
     type="p",
     ylim=ylim,
     main=TeX(paste("Data ($M$ = ",M,", $N$ = ",N,")",sep="")),
     xlab=TeX("$t$"),
     ylab="Data",
     panel.first=grid(lty=1))
if(N>1) {
  for(n in 2:N) {
    points(curves.Brownian[[n]]$t,curves.Brownian[[n]]$x,col=n)
  }
}


## ----mfbmiterhHL--------------------------------------------------------------
#| label: fig-mfbmiterhHL
#| fig-cap: "Multifractional Brownian Motion - Bandwidth Iteration"
par(mfrow=c(1,1))
matplot(rbind(fda.la.Brownian$h.HL.starting,fda.la.Brownian$h.HL.updated.mat),
        xlab="Iteration",
        ylab=TeX("Iterated bandwidths $h_t$"),
        main=TeX(paste("Evolution of Bandwidths $h_t$ for $H_t$ and $L_t$ ($M$ = ",M,", $N$ = ",N,", $t$ = ",length(t0.grid)," grid points)",sep="")),
        type="l",
        panel.first=grid(lty=1))


## ----mfbmiterH----------------------------------------------------------------
#| label: fig-mfbmiterH
#| fig-cap: "Multifractional Brownian Motion - Local Hölder Exponents Iteration"
par(mfrow=c(1,1))
matplot(fda.la.Brownian$H.updated.mat,
        xlab="Iteration",
        ylab=TeX("Iterated $H_t$"),
        main=TeX(paste("Evolution of Local Hölder Exponents $H_t$ ($M$ = ",M,", $N$ = ",N,", $t$ = ",length(t0.grid)," grid points)",sep="")),
        type="l",
        panel.first=grid(lty=1))


## ----mfbmHtrue----------------------------------------------------------------
#| label: fig-mfbmHtrue
#| fig-cap: "Multifractional Brownian Motion - Estimated and True Regularity"
par(mfrow=c(1,1))
plot(t0.grid,fda.la.Brownian$H.updated,
     ylim=c(0,1),
     ylab=TeX("$H_t$"),
     xlab=TeX("t"),
     main=TeX(paste("Estimated and True $H_t$ ($M$ = ",M,", $N$ = ",N,", $t$ = ",length(t0.grid)," grid points)",sep="")),
     type="p",
     pch=19,
     panel.first=grid(lty=1))
lines(t0.grid,true_H)
rug(t0.grid)
legend("topleft",c("True","Estimated"),lty=c(1,NA),pch=c(NA,19),bty="n")


## ----growthH------------------------------------------------------------------
#| label: fig-growthH
#| fig-cap: "Berkeley Growth Example - Estimated Regularity"
require(fda)
curves.growth <- with(growth,vector(mode="list",length=NCOL(hgtf)))
t0.grid <- with(growth,age)
for(j in 1:NCOL(with(growth,hgtf))) {
  curves.growth[[j]]$x <- with(growth,hgtf[,j])
  curves.growth[[j]]$t <- t0.grid
}

mu.starting <- mu.pooled.func(curves.growth,t0.grid,method="np")

fda.la.growth <- fda.la(curves_batch=curves.growth,
                        t.grid=t0.grid,
                        t.lb=min(t0.grid),
                        t.ub=max(t0.grid),
                        common.design=TRUE,
                        mu.hat=mu.starting)

H.growth <- fda.la.growth$H.updated

par(mfrow=c(1,1))
plot(t0.grid,fda.la.growth$H.updated,
     ylim=c(0,1),
     ylab=TeX("$H_t$"),
     xlab=TeX("t"),
     main=TeX(paste("Estimated $H_t$ (",length(t0.grid)," grid points)",sep="")),
     type="p",
     pch=19,
     panel.first=grid(lty=1))
rug(t0.grid)


## ----canWeatherH--------------------------------------------------------------
#| label: fig-canWeatherH
#| fig-cap: "Canadian Weather Example - Estimated Regularity"
require(gamair)
data(canWeather)
place <- unique(CanWeather$place)
curves.canWeather <- vector(mode="list",length=length(place))
t0.grid <- 1:365
for(j in 1:length(place)) {
  curves.canWeather[[j]]$x <- CanWeather$T[CanWeather$place==place[j]]
  curves.canWeather[[j]]$t <- t0.grid
}

mu.starting <- mu.pooled.func(curves.canWeather,t0.grid,method="np")

fda.la.canWeather <- fda.la(curves_batch=curves.canWeather,
                            t.grid=t0.grid,
                            t.lb=min(t0.grid),
                            t.ub=max(t0.grid),
                            common.design=TRUE,
                            mu.hat=mu.starting)

H.canWeather <- fda.la.canWeather$H.updated

par(mfrow=c(1,1))
plot(t0.grid,fda.la.canWeather$H.updated,
     ylim=c(0,1),
     ylab=TeX("$H_t$"),
     xlab=TeX("t"),
     main=TeX(paste("Estimated $H_t$ (",length(t0.grid)," grid points)",sep="")),
     type="p",
     pch=19,
     panel.first=grid(lty=1))
rug(t0.grid)


## ----boxplotH-----------------------------------------------------------------
#| label: fig-boxplotH
#| fig-cap: "Estimated Regularity Comparison"
par(mfrow=c(1,1))
boxplot(list(growth=H.growth,canWeather=H.canWeather,MfBm=H.Brownian),
        notch=TRUE,
        outline=FALSE,
        ylab=TeX("Regularity $H_t$"),
        border = NA,
        xaxt='n',
        yaxt = "n",
        frame = FALSE)
        grid(lty=1)
boxplot(list(growth=H.growth,canWeather=H.canWeather,MfBm=H.Brownian),
        notch=TRUE,
        outline=FALSE,
        ylab=TeX("Regularity $H_t$"),
        add = TRUE,
        ann = FALSE)


## ----growthgamma--------------------------------------------------------------
#| label: fig-growthgamma
#| fig-cap: "Berkeley Growth"
GA::persp3D(fda.la.growth$t.grid[order(fda.la.growth$t.grid)],
            fda.la.growth$t.grid[order(fda.la.growth$t.grid)],
            fda.la.growth$Gamma.hat[order(fda.la.growth$t.grid),order(fda.la.growth$t.grid)],
            xlab="s",
            ylab="t",
            zlab="Estimated Gamma")


## ----canweathergamma----------------------------------------------------------
#| label: fig-canweathergamma
#| fig-cap: "Canadian Weather"
GA::persp3D(fda.la.canWeather$t.grid[order(fda.la.canWeather$t.grid)],
            fda.la.canWeather$t.grid[order(fda.la.canWeather$t.grid)],
            fda.la.canWeather$Gamma.hat[order(fda.la.canWeather$t.grid),order(fda.la.canWeather$t.grid)],
            xlab="s",
            ylab="t",
            zlab="Estimated Gamma")


## ----mfbmgamma----------------------------------------------------------------
#| label: fig-mfbmgamma
#| fig-cap: "Multifractional Brownian Motion"
GA::persp3D(fda.la.Brownian$t.grid[order(fda.la.Brownian$t.grid)],
            fda.la.Brownian$t.grid[order(fda.la.Brownian$t.grid)],
            fda.la.Brownian$Gamma.hat[order(fda.la.Brownian$t.grid),order(fda.la.Brownian$t.grid)],
            xlab="s",
            ylab="t",
            zlab="Estimated Gamma")


## ----Brownianrecon------------------------------------------------------------
#| label: fig-Brownianrecon
#| fig-cap: "MfBm Curve Reconstruction"
par(mfrow=c(1,1))
plot(fda.la.Brownian$t.grid,fda.la.Brownian$X.hat.lp[100,],
     xlab=TeX("$t$"),
     ylab=TeX("$X^{(i)}(t)$, $Y^{(i)}(t)$"),
     type="l",
     panel.first=grid(lty=1),
     col=3)
lines(fda.la.Brownian$t.grid,fda.la.Brownian$X.hat[100,],lty=2,col=2,type="l")
lines(curves.Brownian[[100]]$grid_true,curves.Brownian[[100]]$x_true,lty=1,col="grey",lwd=1.5,type="l")
points(curves.Brownian[[100]]$t, curves.Brownian[[100]]$x,pch=19)
legend("topright",c("Reconstructed Estimator","Local Estimator","True Curve","Data"),lty=c(1,2,1,NA),lwd=c(1,1,1.5,NA),col=c(3,2,"grey",1),pch=c(NA,NA,NA,19),bty="n")


## ----rmseMtabler--------------------------------------------------------------
#| label: tbl-rmseMtable
#| tbl-cap: Mean RMSE after final recursion
load("R_scripts_data/table_skip.RData")
knitr::kable(table.skip.out[-1,c(1,3,6,8,10,11)],
             digits=4,
             row.names=FALSE,
             escape=FALSE,
             format="markdown",
             align='c')


## ----rmseM200Ntabler----------------------------------------------------------
#| label: tbl-rmseM200Ntable
#| tbl-cap: Mean RMSE after final recursion
load("R_scripts_data/table_M_200.RData")
knitr::kable(table.foo.200[,c(-2,-4)],
             digits=4,
             row.names=FALSE,
             escape=FALSE,
             align='c',
             format="markdown")

