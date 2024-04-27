## This function uses the R manipulate package to create a dynamic menu that
## enables one to modify various parameters for a multifractional Brownian
## motion to be used to illustrate the impact of modifying the process. In
## particular, we are interested in demonstration how the local regularity
## parameter Ht (which can vary across the domain) and local Holder constant L
## affect the resulting sample curve element and also the discretely sampled
## points, measured with error (controlled by sigma).

## Comments to racinej@mcmaster.ca (Jeffrey S. Racine)

## Points to note:

## We include the random variable tau in order to avoid all the curves starting
## from the origin. Almost the same effect can be obtained if we generate the
## process on, say, [0,2] and we only retain the sample path for t in [1,2],
## which we next rename [0,1]. It is like in time series when one wants to
## generate a stationary AR(1) and starts from the origin. Then one has to let
## the series "burn in" a bit before starting to look at it. Here there is no
## stationarity of stationary increments requirement, but the idea is similar:
## one wants to avoid the effect of the starting value, hence one has simply to
## let the process "burn in" a bit.

## Setting (unrealistically) L=1 and tau=0 results in a sample curve that is
## equal to the mean curve. In some sense we can demonstrate a more "classical"
## sample drawn from some unknown conditional mean (which can be modified to
## include a constant function, linear function, sin function, etc.)

suppressPackageStartupMessages(require(manipulate))

manipulate.plot <- function(M,
                            L,
                            h_first,
                            h_second,
                            length_grid,
                            sigma,
                            tau,
                            mufunc,
                            plot.mufunc,
                            plot.curve,
                            plot.data,
                            plot.regularity,
                            setseed){
  
  ################################################################################
  ### Data generation process
  ###
  ### Generation of MFBM
  ### See the following paper for the covariance structure
  ### https://doi.org/10.3390/fractalfract6020074
  ### Equations (6) and (7)
  ################################################################################
  
  hurst_logistic <- function(
    t_vec,
    h_left = h_first,
    h_right = h_second,
    slope = h_slope,
    change_point_position = change_point_t
  ) {
    change_point <- change_point_position * (max(t_vec) + min(t_vec))
    u <- (t_vec - change_point) / (max(t_vec) - min(t_vec))
    (h_right - h_left) / (1 + exp(-slope * u)) + h_left
  }
  
  constant_d <- function(x, y) {
    a <- gamma(2 * x + 1) * gamma(2 * y + 1) * sin(pi * x) * sin(pi * y)
    b <- 2 * gamma(x + y + 1) * sin(pi * (x + y) / 2)
    sqrt(a) / b
  }
  
  covariance_mfbm <- function(points, hurst, ...) {
    tmp <- expand.grid(s = points, t = points)
    s <- tmp$s
    t <- tmp$t
    hs <- hurst(s, ...)
    ht <- hurst(t, ...)
    hh <- hs + ht
    values <- constant_d(hs, ht) * (t**hh + s**hh - abs(t - s)**hh)
    matrix(values, ncol = length(points))
  }
  
  ##
  ## Data generating process
  ##
  
  generates_random_length_points <- function(N, m, distribution = runif) {
    M <- rpois(N, m)
    lapply(
      M,
      function(x) {
        sort(distribution(x))
      }
    )
  }
  
  generates_equal_length_points <- function(N, m, distribution = runif) {
    M <- rep(m,N)
    lapply(
      M,
      function(x) {
        sort(distribution(x))
      }
    )
  }
  
  generates_increasing_length_points <- function(N, m, distribution = runif) {
    M <- round(seq(m,m*log(N),length=N))
    lapply(
      M,
      function(x) {
        sort(distribution(x))
      }
    )
  }
  
  generate_curves <- function(
    points_list,
    hurst,
    sigma = 0.1,
    grid = seq(0,1,length=100),
    tau = 0,
    L = 1,
    shift = 0,
    scale_mu = 1,
    ...
  ) {
    
    t_min <- min(grid)
    t_max <- max(grid)
    covariance_list <- lapply(
      points_list,
      function(point,...) {
        pp <- sort(c(point, grid))
        covariance_mfbm(pp, hurst, ...)
      }
    )
    curves <- lapply(
      covariance_list,
      function(covariance) {
        MASS::mvrnorm(
          1,
          mu = rep(tau * rnorm(1), ncol(covariance)),
          Sigma = L * covariance
        )
      }
    )
    x_min <- min(sapply(curves, min))
    x_max <- max(sapply(curves, max))
    #
    curves_random <- list()
    curves_grid <- list()
    for (i in 1:length(curves)) {
      pp <- sort(c(points_list[[i]], grid))
      idx <- which(pp %in% grid)
      curves_grid[[i]] <- curves[[i]][idx]
      curves_random[[i]] <- curves[[i]][-idx]
    }
    
    out_curves <- list()
    for (i in seq_along(points_list)) {
      ideal <- list(
        t = grid,
        x = curves_grid[[i]]
      )
      observed <- list(
        t = points_list[[i]],
        x = curves_random[[i]] +
          rnorm(length(curves_random[[i]]), 0, sigma)
      )
      class(ideal) <- "curves"
      class(observed) <- "curves"
      out_curves[[i]] <- list(
        ideal = ideal,
        observed = observed
      )
    }
    
    out <- list()
    out$t_min <- t_min
    out$t_max <- t_max
    out$x_min <- x_min
    out$x_max <- x_max
    out$curves <- out_curves
    class(out) <- "curves_ideal_observed"
    out
  }
  
  generate_mean_curve <- function(
    grid = seq(0, 1, length.out = 101 ), shift = 0, scale_mu = 1,
    K = 100,
    alpha = 1
  ) {
    K_seq <- 1:K
    Z <- rnorm(K)
    xi <- 1/((K_seq - 1/2)**alpha * pi**2)
    
    mu_kt <- sapply(grid, function(t) {
      sqrt(2) * Z * sqrt(xi) * sin((K_seq - 0.5) * pi * t)
      
    })
    
    mu_t <- scale_mu * colSums(mu_kt) +  rep(shift, length(grid))
    list(t = grid, x = mu_t)
  }
  
  add_mean_curve <- function(data, mu) {
    lapply(data, function(curve) {
      idx <- map_dbl(curve$t, ~which.min(abs(.x - mu$t)))
      
      list(
        t = curve$t,
        x = curve$x + mu$mu[idx]
      )
    })
  }
  
  sample_curve <- function(curve, grid) {
    idx <- map_dbl(grid, ~which.min(abs(.x - curve$t)))
    curve$x[idx]
  }
  
  plot.curves_ideal_observed <- function(
    curves,
    which_curves = 1:length(curves$curves),
    type = "observed", #c("true", "observed")
    ...
  ) {
    courbes <- curves$curves[which_curves]
    my_colors <- rainbow(length(courbes), alpha = 0.5)
    plot(NULL,
         type = "n",
         xlim = c(curves$t_min, curves$t_max),
         ylim = c(curves$x_min, curves$x_max),
         ...
    )
    for (i in seq_along(courbes)) {
      courbe <- courbes[[i]]
      mc <- my_colors[i]
      if (type == "true") {
        with(
          courbe,
          {
            points(observed$t, observed$x, col = mc)
            lines(ideal$t, ideal$x, col = mc)
          }
        )
      } else {
        with(
          courbe,
          {
            lines(observed$t, observed$x, col = mc)
          }
        )
      }
      
    }
  }
  
  ############################ The manipulate main code ###########################################
  
  if(setseed) set.seed(42)
  
  ## Set size of batch (N), changepoint, and slope (these are not meant to be
  ## modified and there should be no reason to do so for the purposes of this
  ## demo)
  
  N = 1
  change_point_t = 0.5
  h_slope = 10
  
  t0.grid = seq(0,1,length=length_grid)
  
  ## Set the mean function
  
  if(mufunc=="constant") {
    mu.func <- function(x) { rep(0,length(x)) }
  } else if (mufunc=="linear") {
    mu.func <- function(x) { 1 + x }
  } else if (mufunc=="sin") {
    mu.func <- function(x) {  sin(2*pi*x) }
  } else if (mufunc=="hardle") {
    mu.func <- function(t) {t+exp(-200*(t-0.5)^2)}
  } else if (mufunc=="wiggly") {  
    mu.func <- function(x){
      #set.seed(12)
      #xi.coef <-rnorm(20)
      xi.coef <- c(-1.48056759,1.57716947,-0.95674448,-0.92000525,-1.99764210,-0.27229604,-0.31534871,-0.62825524,-0.10646388,0.42801480,-0.77771958,-1.29388230,-0.77956651,0.01195176,-0.15241624,-0.70346425,1.18887916,0.34051227,0.50696817,-0.29330515)
      mu <- rep(0, length = length(x))
      for(j in 1:length(xi.coef))
      {mu <- mu + xi.coef[j]*sin(pi*j*x)/(j-1/2) }
      return(mu)
    }
  } else {
    stop("You must enter a mean function type (constant, linear, sin, hardle, wiggly)")
  }
  
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
  
  
  
  
  ## Some plots and true values
  
  true_H <- hurst_logistic(t0.grid)
  
  if(plot.regularity) {
    par(mfrow=c(2,1),cex=.65)
  } else {
    par(mfrow=c(1,1),cex=.65)
  }
  if(plot.curve & plot.data) {
    plot(t0.grid,curves_batch[[1]]$x_true,
         type="l",
         ylim=range(curves_batch[[1]]$x_true,curves_batch[[1]]$x,mu.func(curves_batch[[1]]$t)),
         #main="Multifractional Brownian Motion",
         sub=paste("sigma = ",sigma,", tau = ",tau,", L = ",L,", M = ",M,", curve grid = ",length_grid,sep=""),
         xlab="t",
         ylab="X, Y",
         col="grey",
         panel.first=grid(lty=1))
    points(curves_batch[[1]]$t,curves_batch[[1]]$x,
           col="blue",
           pch=19)
    rug(curves_batch[[1]]$t,col="blue")
  } else if(plot.curve & !plot.data) {
    plot(t0.grid,curves_batch[[1]]$x_true,
         type="l",
         ylim=range(curves_batch[[1]]$x_true,curves_batch[[1]]$x,mu.func(curves_batch[[1]]$t)),
         #main="Multifractional Brownian Motion",
         sub=paste("sigma = ",sigma,", tau = ",tau,", L = ",L,", M = ",M,", curve grid = ",length_grid,sep=""),
         xlab="t",
         ylab="X, Y",
         col="grey",
         panel.first=grid(lty=1))
  } else if(!plot.curve & plot.data) {
    plot(curves_batch[[1]]$t,curves_batch[[1]]$x,
         type="p",
         pch=19,
         ylim=range(curves_batch[[1]]$x_true,curves_batch[[1]]$x,mu.func(curves_batch[[1]]$t)),
         #main="Multifractional Brownian Motion",
         sub=paste("sigma = ",sigma,", tau = ",tau,", L = ",L,", M = ",M,", curve grid = ",length_grid,sep=""),
         xlab="t",
         ylab="X, Y",
         col="blue",
         panel.first=grid(lty=1))
    rug(curves_batch[[1]]$t,col="blue")
  }
  
  if(plot.curve & plot.mufunc & plot.data) {
    lines(t0.grid,mu.func(t0.grid),lty=2,col=2)
    legend("topright",c("Sample Curve (random draw from X)","Mean Curve","Sample Data (curve measured with error at discrete points)"),
           lty=c(1,2,NA),
           pch=c(NA,NA,19),
           pt.cex=c(NA,NA,.5),
           col=c("grey","red","blue"),
           bty="n")   
  } else if(!plot.curve & plot.mufunc & plot.data) {
    lines(t0.grid,mu.func(t0.grid),lty=2,col=2)
    legend("topright",c("Mean Curve","Sample Data (curve measured with error at discrete points)"),
           lty=c(2,NA),
           pch=c(NA,19),
           pt.cex=c(NA,.5),
           col=c("red","blue"),
           bty="n")   
  } else if(plot.curve & plot.mufunc & !plot.data) {
    lines(t0.grid,mu.func(t0.grid),lty=2,col=2)
    legend("topright",c("Sample Curve (random draw from X)","Mean Curve"),
           lty=c(1,2),
           col=c("grey","red"),
           bty="n")       
  } else if(plot.curve & !plot.mufunc & plot.data)  {
    legend("topright",c("Sample Curve (random draw from X)","Sample Data (curve measured with error at discrete points)"),
           lty=c(1,NA),
           pch=c(NA,19),
           pt.cex=c(NA,.5),
           col=c("grey","blue"),
           bty="n")
  } else if(!plot.curve & !plot.mufunc & plot.data)  {
    legend("topright",c("Sample Data (curve measured with error at discrete points)"),
           pch=c(19),
           pt.cex=c(.5),
           col=c("blue"),
           bty="n")
  }  else if(plot.curve & !plot.mufunc & !plot.data)  {
    legend("topright",c("Sample Curve (random draw from X)"),
           lty=c(1),
           col=c("grey"),
           bty="n")
  }
  if(plot.regularity) {
    plot(t0.grid,true_H,
         main="Regularity",
         ylim=c(0,1),
         ylab="True H",
         xlab="t",
         type="l",
         panel.first=grid(lty=1))
    legend("topright",paste("H", " (L = ",L,")",sep=""),lty=1,bty="n")
  }
  
}

manipulate(manipulate.plot(M,
                           L,
                           h_first,
                           h_second,
                           #h_slope,
                           #change_point_t,
                           length_grid,
                           sigma,
                           tau,
                           mufunc,
                           plot.mufunc,
                           plot.curve,
                           plot.data,
                           plot.regularity,
                           setseed),
           h_first=slider(0.05,1,0.75,step=0.05,label="H low (leftmost Holder exponent Ht)"),
           h_second=slider(0.05,1,0.75,step=0.05,label="H high (rightmost Holder exponent Ht)"),
           L=slider(0,2,1,step=0.1,label="L (Holder constant L)"),
           sigma=slider(0,1,0,step=0.05,label="Noise dispersion (sigma)"),
           tau=slider(0,2,1,step=0.1,label="tau"),
           M=slider(20,500,100,step=20,label="M (discrete sample points)"),
           length_grid=slider(100,2000,500,step=100,label="No. equal spaced curve points"),
           mufunc=picker("sin", "hardle", "wiggly", "constant", "linear",label="Mean Function"),
           plot.mufunc=checkbox(FALSE, "Plot Mean Function"),
           plot.curve=checkbox(TRUE, "Plot Curve"),
           plot.data=checkbox(FALSE, "Plot Data"),
           plot.regularity=checkbox(TRUE, "Plot Regularity Function"),
           setseed=checkbox(TRUE, "Set Seed (42)"))


