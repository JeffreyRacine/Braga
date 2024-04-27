## Assess estimator of sigma.hat

source("func_lib.R")

#set.seed(42)

## Set parameters here

## Regularity parameters for generating data

h_first = 0.25
h_second = 0.75
h_slope = 10
change_point_t = 0.5

sigma = 1
tau = 1
L = 10

## Set size of batch (N), (average) number of observations for each online curve
## (M), number of online curves, and grid of points on which to estimate curves

N = 100
M = 100

t0.grid = seq(0,1,length=100)

## Set the mean function

#mu.func <- function(x) { 1 + x }
#mu.func <- function(x) { -cos(2*pi*x) }
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


## See how sigma estimator performs... try M=N=100, grid of, say, 25, sigma=.1,.5,1

sigma.sq.summand <- matrix(NA,N,length(t0.grid))
Tmi.pooled <- numeric()
Ymi.pooled <- numeric()

for(n in 1:N) {
  Tmi <- curves_batch[[n]]$t
  Ymi <- curves_batch[[n]]$x
  Tmi.pooled <- c(Tmi.pooled,Tmi)
  Ymi.pooled <- c(Ymi.pooled,Ymi)
  for(j in 1:length(t0.grid)) {
    ## For computing sigma.hat - need at least 2 sample realizations to
    ## compute. If element is empty will be set to NA and ignored when
    ## computing mean over all N curves. Summand for sigma.sq involves two
    ## closest Ymi values
    if(length(Ymi)>1) {
      close.1 <- which.min(abs(Tmi-t0.grid[j]))
      sigma.sq.summand[n,j] <-  0.5*(Ymi[close.1]-(Ymi[-close.1])[which.min(abs(Tmi[-close.1]-t0.grid[j]))])^2
    }
  }
}

## Estimate sigma (no smoothing required, does not change with
## bandwidth) Compute average sample size, call this Mi.mean

par(mfrow=c(3,1))

plot.true(curves_batch)

plot.data(curves_batch)

sigma.hat <- sqrt(colMeans(sigma.sq.summand, na.rm=TRUE))

print(summary(sigma.hat))
plot(t0.grid,sigma.hat,
     xlab="t",
     ylab="sigma.hat",
     main=paste("True sigma = ",sigma," (homoskedastic)",sep=""),
     sub=paste("Mean sigma.hat = ",formatC(mean(sigma.hat),format="f",digits=4),sep=""),
     ylim=range(0,2*sigma,sigma.hat),
     panel.first=grid(lty=1))
abline(h=sigma)

