## Simple function to demonstrate classical versus FDA sample elements

source("func_lib.R")

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

par(mfrow=c(2,1))

for(n in 1:N) {

  plot(x[1:n],y[1:n],pch=19,
       xlab="X",
       ylab="Y",
       ylim=ylim.xy,
       xlim=c(0,1),
       main=paste("Classical Sample elements (independent random pairs, N = ",n,")",sep=""),
       panel.first=grid(lty=1))
  
  plot(t0.grid,curves_batch[[1]]$x_true,
       type="l",
       ylim=ylim.curve,
       xlim=c(0,1),
       main=paste("FDA Sample elements (independent random curves, N = ",n,")",sep=""),
       xlab="t",
       ylab="Curve",
       panel.first=grid(lty=1))
  
  if(n>1) {
    for(j in 2:n) {
      lines(t0.grid,curves_batch[[j]]$x_true,col=j,lty=j)
    }
  }
  
  readline(prompt="Press [enter] to continue\r")
  
}
