## This contains a set of functions that I wrote (Jeffrey Racine,
## racinej@mcmaster.ca), and supporting functions for data generation those
## Sunny wrote (Sunny Wang <sunnywanggw@gmail.com>). Sunny relies on the
## functional library via :: so make sure it is installed, and this load will
## fail if it is not.

## Feb 18 2023, adding Mi[n]/Mi.mean and Mi/Mi.mean to computation of
## X.hat.gamma which is used to compute covariance matrix. For mu it is correct
## in the batch since we pool, while it is correct in the online since we use
## Robbins-Munro with weights adjusted for different Mi

suppressPackageStartupMessages(library(functional))

## The function fda.la() is designed to be used with a) batches of
## curves, b) online curves one-at-a-time, c) many online curves, and d) batches
## and many online curves.  When processing batches it uses an iterative
## procedure over the batch files to determine H and L. Then it pools the data
## to estimate the mean function mu(t) and covariance function Gamma(s,t). When
## processing online curves it uses recursion using only the new information in
## the curve plus simple input from its previous invocation (or invocation of a
## batch). In this way it is designed for cases where the batch is simply too
## large to contemplate processing in real-time hence an online approach (one
## new online curve at a time) has obvious appeal. Unlike processing a batch of
## files, recursion is used to incorporate information from batch (existing)
## curves.

fda.la <- function(common.design=FALSE,
                   compute.lp.estimator=TRUE,
                   curves_batch=NULL,
                   curves_online=NULL,
                   delta=NULL,
                   degenerate.method.HL=c("skip","1nn","lp"),
                   degenerate.method.mu=c("skip","1nn","lp"),
                   degenerate.method.Gamma=c("skip","1nn","lp"),
                   f.hat.grid=NULL,
                   fda.la.previous=NULL,
                   Gamma.bandwidth=c("logNMsq","NMsq"),
                   Gamma.hat.colSums.vec=NULL,
                   Gamma.hat.crossprod.mat=NULL,
                   Gamma.hat=NULL,
                   h.Gamma=NULL,
                   h.mu=NULL,
                   h.HL.starting=NULL,
                   h.HL=NULL,
                   H.lb=0.1,
                   H.ub=1.0,
                   issue.warnings=TRUE,
                   L.sq.lb=.Machine$double.eps,
                   mean.correct.HL=TRUE,
                   Mi.mean=NULL,
                   Mi.sq.mean=NULL,
                   mu.bandwidth=c("logNM","M","NM"),
                   mu.hat=NULL,
                   mu.method=c("default","Bspline.cv","kernel.cv"),
                   N=NULL,
                   num.iter=25,
                   plot.batch.results=FALSE,
                   plot.online.results=FALSE,
                   sigma.hat=NULL,
                   t.grid=NULL,
                   t.lb=0,
                   t.ub=1,
                   theta.hat.t1.t2=NULL,
                   theta.hat.t1.t3=NULL,
                   theta.hat.t2.t3=NULL,
                   X.hat.count=NULL,
                   X.hat.mu.count=NULL,
                   X.hat.Gamma.count=NULL,
                   X.hat.Gamma.count.mat=NULL){
  
  ## Note - we are overloading this function so that if curves_batch is not
  ## passed, then we presume N, Mi.mean, sigma.hat, f.hat.grid, h, and the
  ## theta.hat.*.* have all been passed correctly (they can be extracted
  ## automatically from a previous invocation of the function that was assigned
  ## to an object via fda.la.previous=foo where foo was generated via
  ## foo <- fda.la(...)). The intent is to be able to use this function
  ## purely for online purposes (get one curve, update existing objects, push
  ## them out). So, this function could be used with a) set of batch and set of
  ## online curves, b) no batch and set of online, or c) no batch and 1 online.
  
  ## Basic checking for valid input is done in two places, a) some initial
  ## checks (did you even pass batch or online curves?), then b) at the online
  ## stage that appears further down in the code.
  
  ## mu.method is used for set of batch curves and computes mean on pooled data
  ## using h.mu computed in the batch. degenerate.method determines how
  ## degenerate curve points are handled.
  
  mu.method <- match.arg(mu.method)
  degenerate.method.HL <- match.arg(degenerate.method.HL)
  degenerate.method.mu <- match.arg(degenerate.method.mu)
  degenerate.method.Gamma <- match.arg(degenerate.method.Gamma)
  
  ## If, in the batch stage only, you wish to compute the pooled mean using
  ## instead splines or global bandwidth kernel regression, we load the
  ## appropriate libraries and use their functions.
  
  if(mu.method=="Bspline.cv") { suppressPackageStartupMessages(library(crs));options(crs.messages=FALSE) }
  if(mu.method=="kernel.cv") { suppressPackageStartupMessages(library(np));options(np.tree=TRUE,np.messages=FALSE) }
  if(plot.batch.results | plot.online.results) { suppressPackageStartupMessages(library(GA)) }
  
  ## mu.bandwidth and Gamma.bandwidth determine essentially whether to use Mi
  ## (MSE bandwidth for one online curve), N x Mi, or N x Mi^2 for mu and Gamma
  ## estimation.
  
  mu.bandwidth <- match.arg(mu.bandwidth)
  Gamma.bandwidth <- match.arg(Gamma.bandwidth)
  
  if(is.null(t.grid) & is.null(fda.la.previous)) stop("You must provide t.grid")
  length.grid <- length(t.grid)
  
  if(is.null(curves_batch) & is.null(fda.la.previous)) {
    if(is.null(N)) stop("Must provide N")
    if(is.null(Mi.mean)) stop("Must provide Mi.mean")
    if(is.null(sigma.hat)) stop("Must provide sigma.hat")
    if(is.null(f.hat.grid)) stop("Must provide f.hat.grid")
    if(is.null(h.HL)) stop("Must provide h.HL")
    if(is.null(theta.hat.t1.t3)) stop("Must provide theta.hat.t1.t3")
    if(is.null(theta.hat.t1.t2)) stop("Must provide theta.hat.t1.t2")
    if(is.null(theta.hat.t2.t3)) stop("Must provide theta.hat.t2.t3")
  } else {
    N <- length(curves_batch)
  }
  
  if(!is.null(curves_online)) {
    N.online <- length(curves_online)
  } else {
    N.online <- 0
    i.online <- N
  }
  
  ## Further checking in case invalid values are provided
  
  if(!is.null(Mi.mean)&any(is.na(Mi.mean))|any(is.infinite(Mi.mean))|any(is.nan(Mi.mean))) stop("Mi.mean contains invalid entries (i.e., NA, Inf, or NaN)")
  if(!is.null(mu.hat)&any(is.na(mu.hat))|any(is.infinite(mu.hat))|any(is.nan(mu.hat))) stop("mu.hat contains invalid entries (i.e., NA, Inf, or NaN)")
  if(!is.null(f.hat.grid)&any(is.na(f.hat.grid))|any(is.infinite(f.hat.grid))|any(is.nan(f.hat.grid))) stop("f.hat.grid contains invalid entries (i.e., NA, Inf, or NaN)")  
  if(!is.null(h.HL)&any(is.na(h.HL))|any(is.infinite(h.HL))|any(is.nan(h.HL))) stop("h.HL contains invalid entries (i.e., NA, Inf, or NaN)")  
  if(!is.null(sigma.hat)&any(is.na(sigma.hat))|any(is.infinite(sigma.hat))|any(is.nan(sigma.hat))) stop("sigma.hat contains invalid entries (i.e., NA, Inf, or NaN)")
  if(!is.null(theta.hat.t1.t3)&any(is.na(theta.hat.t1.t3))|any(is.infinite(theta.hat.t1.t3))|any(is.nan(theta.hat.t1.t3))) stop("theta.hat.t1.t3 contains invalid entries (i.e., NA, Inf, or NaN)")  
  if(!is.null(theta.hat.t1.t2)&any(is.na(theta.hat.t1.t2))|any(is.infinite(theta.hat.t1.t2))|any(is.nan(theta.hat.t1.t2))) stop("theta.hat.t1.t2 contains invalid entries (i.e., NA, Inf, or NaN)")  
  if(!is.null(theta.hat.t2.t3)&any(is.na(theta.hat.t2.t3))|any(is.infinite(theta.hat.t2.t3))|any(is.nan(theta.hat.t2.t3))) stop("theta.hat.t2.t3 contains invalid entries (i.e., NA, Inf, or NaN)")  
  
  if(is.null(curves_online)&is.null(curves_batch)) stop("Must provide either online or batch curve(s)")
  
  ## Remove any "Null" curves (i.e., zero observation members of the list of
  ## curves), issue warning, also retrieve min/max of Tmi for setting some
  ## sensible defaults (starting bandwidth, delta, grid) if not provided
  
  if(!is.null(curves_batch)){
    
    list.obs <- numeric()
    min.Tmi <- numeric()
    max.Tmi <- numeric()
    
    for(j in 1:N) {
      Tmi <- curves_batch[[j]]$t
      list.obs[j] <- length(Tmi)
      min.Tmi <- min(min.Tmi,Tmi)
      max.Tmi <- max(max.Tmi,Tmi)
    }
    
    curves_batch <- curves_batch[which(list.obs > 0)]
    if(issue.warnings & any(list.obs<1)) warning("NULL curve(s) encountered (i.e., containing zero observations), removed from batch",immediate.=TRUE)
    
    ## If the user neglected to set t.lb and t.ub to be appropriate
    ## for their data, correct and optionally issue a warning if
    ## requested (correct means use t.lb=min() and t.ub=max(), i.e.,
    ## empirical support)
    
    if(min.Tmi < t.lb) { t.lb <- min.Tmi; if(issue.warnings) warning("t.lb lies above min(Tmi), resetting",immediate.=TRUE) }
    if(max.Tmi > t.ub) { t.ub <- max.Tmi; if(issue.warnings) warning("t.ub lies below max(Tmi), resetting",immediate.=TRUE) }
    
    ## If no grid is passed use a default grid of 100 points.
    
    if(is.null(t.grid)) t.grid <- seq(min.Tmi,max.Tmi,length=100)
    length.grid <- length(t.grid)
    
    X.hat.t1 <- matrix(NA,N,length.grid)
    X.hat.t2 <- matrix(NA,N,length.grid)
    X.hat.t3 <- matrix(NA,N,length.grid)
    X.hat.Gamma <- matrix(NA,N,length.grid)
    X.hat <- matrix(NA,N,length.grid)
    X.hat.lp <- matrix(NA,N,length.grid)
    
    ## If no mean curve is provided, start with constant vector
    
    if(is.null(mu.hat)) mu.hat <- rep(0,length.grid)
    
    H.updated.mat <- matrix(NA,num.iter,length.grid)
    L.sq.updated.mat <- matrix(NA,num.iter,length.grid)
    h.HL.updated.mat <- matrix(NA,num.iter,length.grid)
    X.hat.mu.updated.mat <- NULL
    
    sigma.sq.summand <- matrix(NA,N,length.grid)
    
    Mi <- numeric()
    
    ## Tmi.pooled needed for histogram density estimator used in denominator of
    ## MSE-optimal bandwidth formula, Ymi.pooled and Tmi.pooled for the pooled
    ## batch estimator of mu.
    
    Tmi.pooled <- numeric()
    Ymi.pooled <- numeric()
    
    for(n in 1:N) {
      Tmi <- curves_batch[[n]]$t
      Ymi <- curves_batch[[n]]$x
      Tmi.pooled <- c(Tmi.pooled,Tmi)
      Ymi.pooled <- c(Ymi.pooled,Ymi)
      Mi[n] <- length(Tmi)
      for(j in 1:length.grid) {
        ## For computing sigma.hat - need at least 2 sample realizations to
        ## compute. If element is empty will be set to NA and ignored when
        ## computing mean over all N curves. Summand for sigma.sq involves two
        ## closest Ymi values
        if(length(Ymi)>1) {
          close.1 <- which.min(abs(Tmi-t.grid[j]))
          sigma.sq.summand[n,j] <-  0.5*(Ymi[close.1]-(Ymi[-close.1])[which.min(abs(Tmi[-close.1]-t.grid[j]))])^2
        }
      }
    }
    
    ## Estimate sigma (no smoothing required, does not change with
    ## bandwidth) Compute average sample size, call this Mi.mean
    
    sigma.hat <- sqrt(colMeans(sigma.sq.summand, na.rm=TRUE))
    
    ## Compute the average sample size in the batch of N curves (here Mi is a
    ## vector of length N) and the average of the squared sample sizes
    
    ## Feb 20 - use harmonic.mean()
    
    Mi.mean <- mean(Mi)
    Mi.sq.mean <- mean(Mi^2)
    
    ## Rule-of-thumb for sample information...
    
    if(issue.warnings & Mi.mean < 16) warning(paste("the average number of points per curve (",format(Mi.mean,digits=2),") is unreasonably small for estimating regularity, proceed with extreme caution",sep=""),immediate.=TRUE)
    
    ## Seat of the pants check... if data is drawn from continuous
    ## distribution and _not_ common for all curves (i.e., no "common
    ## design") then there should be no repeated values (unless there is
    ## substantial rounding applied to the data), so issue a simple
    ## heads up if less than 1/2 of the support are unique (better than
    ## ignoring)
    
    if(issue.warnings & identical(unique(as.numeric(Tmi.pooled)),unique(as.numeric(t.grid))) & common.design==TRUE & length(unique(Tmi.pooled))==Mi.mean) {
      warning("Common design encountered, inputs appear correct (all input checks passed)")
    } else if(issue.warnings & length(unique(Tmi.pooled))==Mi.mean & N > 1) {
      warning("Common design encountered (number of unique design points equals mean number of observations over all N curves)\n    1. ensure you have set the option common.design=TRUE in your function call\n    2. also ensure your grid points exist exclusively of the common design points",immediate.=TRUE)
    }
    
    ## Histogram density, trivial hack for closed right to get density at right
    ## bound and left, we are presuming N x M data points is large so this is
    ## truly negligible. We estimate the density on the grid points, this appears
    ## in the denominator of the MSE-optimal bandwidth. We avoid standard smooth
    ## estimators due to a) boundary effects, and b) large sample
    ## sizes/computational overhead. The non-smooth nonparametric histogram
    ## density estimator should be appropriate in this instance. Note f.hat is
    ## based on pooled Tmi.
    
    f.hat.hist <- hist(Tmi.pooled,plot=FALSE)
    f.hat.grid <- f.hat.hist$density[findInterval(t.grid, f.hat.hist$breaks,all.inside=TRUE)]
    
    ## Delta.star is the only "nuisance parameter" to be set, uses theoretical
    ## justification. Distance between t1 and t2 of 3 consecutive equally spaced
    ## points (t1,t2,t3) is Delta/2, so between any 2 consecutive points is
    ## Delta/4.  We compute H at t2 +- delta (i.e., at t1=t2-delta and
    ## t3=t1+delta). Set Delta.star if it is not passed.
    
    ## This calculation is for [0,1] support - multiply delta times
    ## empirical range otherwise
    
    if(is.null(delta)) {
      Delta.star = exp(-log(Mi.mean)^{1/3})
      delta <- (t.ub-t.lb)*Delta.star/4
    }
    
    if(issue.warnings & delta > (t.ub-t.lb)/2) warning("unreasonably large delta, recommended delta < (t.ub-t.lb)/2",immediate.=TRUE)
    
    ## Starting vector of bandwidths for initial estimation of H, all equal to
    ## input starting scalar, will be updated/overwritten once iteration
    ## commences
    
    if(is.null(h.HL.starting)) {
      h.HL <- h.Gamma <- h.HL.starting <- rep(0.001*IQR(Tmi.pooled),length.grid)
    } else {
      h.HL <- h.Gamma <- h.HL.starting <- rep(h.HL.starting,length.grid)
    }
    h.HL[h.HL>0.5*(t.ub-t.lb)-delta] <- 0.5*(t.ub-t.lb)-delta
    
    ## Iterate the batch-based procedure beginning with the starting bandwidths
    ## provided - you can check visually whether convergence of the vector of
    ## bandwidths converged by invoking the function with the argument
    ## plot.batch.results=TRUE and examining the figures produced.
    
    for(cc in 1:num.iter) {
      
      ## Compute the kernel smooths of the N batch curves on
      ## (t-delta,t,t+delta), i.e., on (t1,t2,t3) evaluated at each grid point
      
      t.grid.boundary <- t.grid
      grid.corr <- (delta+h.HL)
      t.grid.boundary[t.grid<(t.lb+grid.corr)] <- (t.lb+grid.corr)[t.grid<(t.lb+grid.corr)]
      t.grid.boundary[t.grid>(t.ub-grid.corr)] <- (t.ub-grid.corr)[t.grid>(t.ub-grid.corr)]
      
      for(n in 1:N) {
        Tmi <- curves_batch[[n]]$t
        Ymi <- Ymi.demeaned <- curves_batch[[n]]$x
        X.hat[n,] <- NW.estimator(Ymi,Tmi,h.HL,t.grid,"1nn")$Y.hat
        if(mean.correct.HL) Ymi.demeaned <- Ymi-mu.hat[sapply(Tmi, function(x) which.min(abs(t.grid - x)))]
        X.hat.t1[n,] <- NW.estimator(Ymi.demeaned,Tmi,h.HL,t.grid.boundary-delta,"1nn")$Y.hat
        X.hat.t2[n,] <- NW.estimator(Ymi.demeaned,Tmi,h.HL,t.grid.boundary,"1nn")$Y.hat
        X.hat.t3[n,] <- NW.estimator(Ymi.demeaned,Tmi,h.HL,t.grid.boundary+delta,"1nn")$Y.hat
        ## Adjust for unequal sample sizes when computing Gamma crossproducts etc.
        X.hat.Gamma[n,] <- (Mi[n]/Mi.mean)*NW.estimator(Ymi,Tmi,h.Gamma,t.grid,"1nn")$Y.hat
      }
      
      ## With the N estimates of the curves at the grid points, compute
      ## the theta values needed to construct \hat H_t and \hat L_t
      
      theta.hat.t1.t3 <- colMeans((X.hat.t1-X.hat.t3)^2)
      theta.hat.t1.t2 <- colMeans((X.hat.t1-X.hat.t2)^2)
      theta.hat.t2.t3 <- colMeans((X.hat.t2-X.hat.t3)^2)
      
      ## Trap zeros (feed instead .Machine$double.eps hence numerator is
      ## zero for log diff)
      
      theta.hat.t1.t3[theta.hat.t1.t3==0] <- .Machine$double.eps
      theta.hat.t1.t2[theta.hat.t1.t2==0] <- .Machine$double.eps
      theta.hat.t2.t3[theta.hat.t2.t3==0] <- .Machine$double.eps
      
      H.updated.mat[cc,] <- (log(theta.hat.t1.t3)-log((theta.hat.t1.t2+theta.hat.t2.t3)/2))/(2*log(2))
      
      ## Trap cases where estimated H is negative (cannot occur
      ## theoretically) OR exceeds 1
      
      H.updated.mat[cc,H.updated.mat[cc,] <= H.lb] <- H.lb
      H.updated.mat[cc,H.updated.mat[cc,] > H.ub] <- H.ub
      
      ## Note we take the mean of estimates on |t1-t2|=|t2-t3|=delta. If instead
      ## you used theta.hat.t1.t3 then you need (2*delta)^{} in the denominator
      ## where 2*delta=|t1-t3|
      
      L.sq.updated.mat[cc,] <- 0.5*(theta.hat.t1.t2+theta.hat.t2.t3)/(delta^{2*H.updated.mat[cc,]})
      
      ## Trap case where estimated L^2 evolves to NA
      
      L.sq.updated.mat[cc,L.sq.updated.mat[cc,]<=.Machine$double.eps] <- L.sq.lb
      
      ## Update MSE optimal h using iterated values for the unknown constants C1
      ## and C2, H, and L, will be used as the starting value in the next batch
      ## iteration. Note N*Mi.mean is the sum of the sample sizes, N*Mi.sq.mean
      ## the sum of the squared sample sizes.
      
      C1.updated <- 3*L.sq.updated.mat[cc,]/((2*H.updated.mat[cc,]+1)*(2*H.updated.mat[cc,]+3))
      C2.updated <- (3/5)*sigma.hat^2
      
      h.HL <- (C2.updated/(2*H.updated.mat[cc,]*C1.updated*Mi.mean*f.hat.grid))**(1/(2*H.updated.mat[cc,]+1))
      h.HL[h.HL>0.5*(t.ub-t.lb)-delta] <- 0.5*(t.ub-t.lb)-delta
      h.HL.updated.mat[cc,] <- h.HL
      
      if(mu.bandwidth=="M") {
        h.mu <- h.HL
      } else if(mu.bandwidth=="NM"){
        h.mu <- (C2.updated/(2*H.updated.mat[cc,]*C1.updated*N*Mi.mean*f.hat.grid))**(1/(2*H.updated.mat[cc,]+1))
      } else if(mu.bandwidth=="logNM"){
        h.mu <- (C2.updated*log(N*Mi.mean)/(2*H.updated.mat[cc,]*C1.updated*N*Mi.mean*f.hat.grid))**(1/(2*H.updated.mat[cc,]+1))
      }
      if(Gamma.bandwidth=="NMsq") {
        h.b1 <- (C2.updated/(4*H.updated.mat[cc,]*C1.updated*N*Mi.mean*f.hat.grid))**(1/(2*H.updated.mat[cc,]+1))
        h.b2 <- (C2.updated/(4*H.updated.mat[cc,]*C1.updated*N*Mi.sq.mean*f.hat.grid))**(1/(2*H.updated.mat[cc,]+2))
        h.Gamma <- pmax(h.b1,h.b2)
      } else if(Gamma.bandwidth=="logNMsq"){
        h.b1 <- (C2.updated*log(N*Mi.mean)/(4*H.updated.mat[cc,]*C1.updated*N*Mi.mean*f.hat.grid))**(1/(2*H.updated.mat[cc,]+1))
        h.b2 <- (C2.updated*log(N*Mi.sq.mean)/(4*H.updated.mat[cc,]*C1.updated*N*Mi.sq.mean*f.hat.grid))**(1/(2*H.updated.mat[cc,]+2))
        h.Gamma <- pmax(h.b1,h.b2)
      }
      
    }
    
    ## Generate mu and Gamma for the batch once we have iterated and updated the
    ## vectors C1, C2, H, and L. mu is computed on pooled batch data, Gamma on
    ## the batch of estimated curves.
    
    ## Don't need to adjust for different Mi from curve to curve since we don't
    ## take the mean of the curves. If we did then adjusting would produce
    ## identical results when instead using the mean of the curves at the grid
    ## points.
    
    if(mu.method=="default") {
      
      ## Compute the kernel smooth on the pooled data with the computed iterated bandwidths
      
      mu.hat <- NW.estimator(Ymi.pooled,Tmi.pooled,h.mu,t.grid,"1nn")$Y.hat
      
    } else if(mu.method == "Bspline.cv") {
      
      mu.hat <- suppressWarnings(predict(crs(Ymi.pooled~Tmi.pooled,cv="exhaustive"),
                                         newdata=data.frame(Tmi.pooled=t.grid)))
      
    } else if(mu.method == "kernel.cv") {
      
      mu.hat <- predict(npreg(Ymi.pooled~Tmi.pooled,ckertype="epanechnikov"),
                        newdata=data.frame(Tmi.pooled=t.grid))
    }
    
    ## Compute the batch covariance matrix Gamma at the grid points, pass
    ## forward the crossproducts and column sums for dynamic updating of newly
    ## acquired online curves.
    
    Gamma.hat.crossprod.mat <- crossprod(X.hat.Gamma,X.hat.Gamma)
    Gamma.hat.colSums.vec <- colSums(X.hat.Gamma)
    Gamma.hat <- Gamma.hat.crossprod.mat/N-tcrossprod(Gamma.hat.colSums.vec/N,Gamma.hat.colSums.vec/N)
    
    ## Given Gamma.hat is now constructed, we revisit the N batch curves and compute the
    ## "reconstructed" (linear-projection) estimator
    
    if(compute.lp.estimator) {
      for(n in 1:N) {
        Tmi <- curves_batch[[n]]$t
        Ymi <- curves_batch[[n]]$x
        X.hat.lp[n,] <- lp.estimator(Ymi,
                                     Tmi,
                                     t.grid,
                                     Gamma.hat,
                                     sigma.hat**2,
                                     mu.hat)
      }
    }
    
    if(plot.batch.results) {
      
      # Some simple plotting which can be helpful (say, for confirming
      # convergence of bandwidths)
      
      par(mfrow=c(2,3))
      
      matplot(rbind(h.HL.starting,h.HL.updated.mat),
              xlab="Iteration",
              ylab="Iterated bandwidth & starting value",
              main="Evolution of bandwidths",
              sub=paste("Mean bandwidth: ", format(mean(h.HL.updated.mat[num.iter,]),digits=3),sep=""),
              type="l",
              panel.first=grid(lty=1))
      
      matplot(H.updated.mat,
              xlab="Iteration",
              ylab="Iterated H",
              main="Evolution of H",
              sub=paste("Mean H: ", format(mean(H.updated.mat[num.iter,]),digits=3),sep=""),
              type="l",
              panel.first=grid(lty=1))
      
      matplot(sqrt(L.sq.updated.mat),
              xlab="Iteration",
              ylab="Iterated L",
              main="Evolution of L",
              sub=paste("Mean L.sq: ", format(mean(L.sq.updated.mat[num.iter,]),digits=3),sep=""),
              type="l",
              panel.first=grid(lty=1))
      
      matplot(t.grid[order(t.grid)], t(X.hat[,order(t.grid)]),
              xlab="grid",
              ylab="Estimated Curves",
              main="Estimated Curves",
              type="l",
              panel.first=grid(lty=1))
      
      plot(t.grid[order(t.grid)],
           mu.hat[order(t.grid)],
           xlab="grid",
           ylab="Estimated mean curve",
           main="Estimated mean curve",
           sub=paste("Mean curve estimation method: ",mu.method,sep=""),
           type="b",
           panel.first=grid(lty=1))
      
      GA::persp3D(t.grid[order(t.grid)],
                  t.grid[order(t.grid)],
                  Gamma.hat[order(t.grid),order(t.grid)],
                  xlab="s",
                  ylab="t",
                  zlab="Estimated Gamma",
                  main="Estimated Gamma",
                  sub=paste("Mean Gamma: ",format(mean(Gamma.hat),digits=3),sep=""))
      
    }
    
    ## Pass forward zero counts to be used when object fda.la.previous
    ## is fed upon subsequent invocation of the function for online data. For
    ## Gamma counts, we need to construct counts for colSums() and crossprod() -
    ## former generated by call to NW function, latter computed manually by
    ## necessity. When processing curves one-at-a-time need to keep track of
    ## this and read it in and push it out.
    
    X.hat.count <- rep(0,length.grid)
    X.hat.mu.count <- rep(0,length.grid)
    X.hat.Gamma.count <- rep(0,length.grid)
    X.hat.Gamma.count.mat <- matrix(0,length.grid,length.grid)
    
  }
  
  ## Here in the code the batch block above either completed running or was not
  ## called (i.e., the function was invoked with online data, which relies on an
  ## initial estimate on batch data typically). Now we head into the recursion
  ## code block for processing the online curves, either one-at-a-time or a
  ## series jointly. Typically, at this stage we have our initial estimates of H
  ## & L based on the batch, which we update as we process the online curves
  
  ## NOTE - Valentin and I are in total agreement to not throw out curves in
  ## batch mode but instead use the 1nn estimator for the degenerate case. The
  ## batch is intended to "get things going" or "bootstrap" in the
  ## non-statistical sense of the term ("pull yourself up by your own
  ## bootstraps"). So above we always use the 1nn estimator (this is hard coded
  ## so unaffected by the option degenerate.method=). Below we branch depending
  ## on whether we "skip" or "1nn" in the degenerate case.
  
  if(is.null(curves_batch) & !is.null(curves_online) & !is.null(fda.la.previous)) {
    
    ## If a previous fda.la object was passed in (useful for online
    ## processing and not having to remember which entries are required to be
    ## passed), we extract all necessary information from the previous
    ## invocation (options, counts etc.)
    
    compute.lp.estimator <- fda.la.previous$compute.lp.estimator
    t.grid <- fda.la.previous$t.grid
    length.grid <- fda.la.previous$length.grid
    degenerate.method.HL <- fda.la.previous$degenerate.method.HL
    degenerate.method.mu <- fda.la.previous$degenerate.method.mu
    degenerate.method.Gamma <- fda.la.previous$degenerate.method.Gamma
    delta <- fda.la.previous$delta
    plot.online.results <- fda.la.previous$plot.online.results
    N <- fda.la.previous$N
    mu.hat <- fda.la.previous$mu.hat
    mean.correct.HL <- fda.la.previous$mean.correct.HL
    Gamma.hat <- fda.la.previous$Gamma.hat
    Gamma.hat.flattop <- fda.la.previous$Gamma.hat.flattop.correction
    Gamma.hat.crossprod.mat <- fda.la.previous$Gamma.hat.crossprod.mat
    Gamma.hat.colSums.vec <- fda.la.previous$Gamma.hat.colSums.vec
    mu.bandwidth <- fda.la.previous$mu.bandwidth
    Gamma.bandwidth <- fda.la.previous$Gamma.bandwidth
    Mi.mean <- fda.la.previous$Mi.mean
    Mi.sq.mean <- fda.la.previous$Mi.sq.mean
    sigma.hat <- fda.la.previous$sigma.hat
    f.hat.grid <- fda.la.previous$f.hat.grid
    h.HL <- fda.la.previous$h.HL
    h.mu <- fda.la.previous$h.mu
    h.Gamma <- fda.la.previous$h.Gamma
    theta.hat.t1.t3 <- fda.la.previous$theta.hat.t1.t3
    theta.hat.t1.t2 <- fda.la.previous$theta.hat.t1.t2
    theta.hat.t2.t3 <- fda.la.previous$theta.hat.t2.t3
    X.hat.count <- fda.la.previous$X.hat.count
    X.hat.mu.count <- fda.la.previous$X.hat.mu.count
    X.hat.Gamma.count <- fda.la.previous$X.hat.Gamma.count
    X.hat.Gamma.count.mat <- fda.la.previous$X.hat.Gamma.count.mat
    issue.warnings <- fda.la.previous$issue.warnings
    
  }
  
  if(is.null(curves_batch) & !is.null(curves_online) & is.null(fda.la.previous)) {
    
    ## If previous fda.la object was not passed in, we check (or double
    ## check) that necessary non-null objects were passed
    
    if(is.null(t.grid)) stop("must pass t.grid")
    if(is.null(length.grid)) stop("must pass length.grid")
    if(is.null(degenerate.method.HL)) stop("must pass dengenerate.method.HL")
    if(is.null(degenerate.method.mu)) stop("must pass dengenerate.method.mu")
    if(is.null(degenerate.method.Gamma)) stop("must pass dengenerate.method.Gamma")
    if(is.null(delta)) stop("must pass delta")
    if(is.null(N)) stop("must pass N")
    if(is.null(mu.hat)) stop("must pass mu.hat")
    if(is.null(mean.correct.HL)) stop("must pass mean.correct.HL")
    if(is.null(Gamma.hat)) stop("must pass Gamma.hat")
    if(is.null(Gamma.hat.crossprod.mat)) stop("must pass Gamma.hat.crossprod.mat")
    if(is.null(Gamma.hat.colSums.vec)) stop("must pass Gamma.hat.colSums.vec")
    if(is.null(mu.bandwidth)) stop("must pass mu.bandwidth")
    if(is.null(Gamma.bandwidth)) stop("must pass Gamma.bandwidth")
    if(is.null(Mi.mean)) stop("must pass Mi.mean")
    if(is.null(Mi.sq.mean)) stop("must pass Mi.sq.mean")
    if(is.null(sigma.hat)) stop("must pass sigma.hat")
    if(is.null(f.hat.grid)) stop("must pass f.hat.grid")
    if(is.null(h.HL)) stop("must pass h.HL")
    if(is.null(h.mu)) stop("must pass h.mu")
    if(is.null(h.Gamma)) stop("must pass h.Gamma")
    if(is.null(theta.hat.t1.t3)) stop("must pass theta.hat.t1.t3")
    if(is.null(theta.hat.t1.t2)) stop("must pass theta.hat.t1.t2")
    if(is.null(theta.hat.t2.t3)) stop("must pass theta.hat.t2.t3")
    if(is.null(X.hat.count)) stop("must pass X.hat.count")
    if(is.null(X.hat.mu.count)) stop("must pass X.hat.mu.count")
    if(is.null(X.hat.Gamma.count)) stop("must pass X.hat.Gamma.count)")
    if(is.null(X.hat.Gamma.count.mat)) stop("must pass X.hat.Gamma.count.mat")
    
  }
  
  ## Now we process the online curves
  
  if(!is.null(curves_online)) {
    
    if(issue.warnings & delta > (t.ub-t.lb)/2) warning("unreasonably large delta, recommended delta < (t.ub-t.lb)/2",immediate.=TRUE)
    
    sigma.sq.summand <- numeric()
    
    H.updated.mat <- matrix(NA,N.online,length.grid)
    L.sq.updated.mat <- matrix(NA,N.online,length.grid)
    h.HL.updated.mat <- matrix(NA,N.online,length.grid)
    X.hat.mu.updated.mat <- matrix(NA,N.online,length.grid)
    X.hat <- matrix(NA,N.online,length.grid)
    X.hat.lp <- matrix(NA,N.online,length.grid)
    
    ## Counts for X.hat and X.hat.mu used to adjust weights. If not passed in
    ## start from zero.
    
    if(is.null(X.hat.count)) X.hat.count <- rep(0,length.grid)
    if(is.null(X.hat.mu.count)) X.hat.mu.count <- rep(0,length.grid)
    
    ## Counts for colSums() and crossprod() - former generated by call to NW
    ## function, latter computed manually by necessity. When processing curves
    ## one-at-a-time need to keep track of this and read it in and push it out.
    
    if(is.null(X.hat.Gamma.count)) X.hat.Gamma.count <- rep(0,length.grid)
    if(is.null(X.hat.Gamma.count.mat)) X.hat.Gamma.count.mat <- matrix(0,length.grid,length.grid)
    
    for(n in 1:N.online) {
      
      ## i.online is the number of online curves already processed. N is the number
      ## in the initial batch (passed in if batch called with previous
      ## invocation of this function), i.online gets augmented if you pass in a set
      ## of > 1 online curves at once (N.online curves as once).
      
      i.online <- N + n
      
      Tmi <- curves_online[[n]]$t
      Ymi <- Ymi.demeaned <- curves_online[[n]]$x
      Mi <- length(Ymi)
      
      ## Here we are computing one online curve, save in matrix with each row
      ## the X.hat for each online curve if a set of online curves is fed in
      
      X.hat[n,] <- NW.estimator(Ymi,Tmi,h.HL,t.grid,degenerate.method.HL,Gamma.hat,sigma.hat**2,mu.hat)$Y.hat
      
      is.na.X.hat.n <- is.na(X.hat[n,])
      X.hat.count[is.na.X.hat.n] <- X.hat.count[is.na.X.hat.n] + 1
      
      ## Create weights for gamma & (1-gamma) for updating the mean online
      ## sample size
      
      gamma.fullsample <- (i.online/(i.online+1))
      
      ## Update mean sample and square mean sample size
      
      Mi.mean.previous <- Mi.mean
      Mi.mean <- gamma.fullsample*Mi.mean + (1-gamma.fullsample)*Mi
      Mi.sq.mean <- gamma.fullsample*Mi.sq.mean + (1-gamma.fullsample)*Mi^2
      
      if(min(Tmi) < t.lb) { t.lb <- min.Tmi; if(issue.warnings) warning("t.lb lies above min(Tmi), resetting",immediate.=TRUE) }
      if(max(Tmi) > t.ub) { t.ub <- max.Tmi; if(issue.warnings) warning("t.ub lies below max(Tmi), resetting",immediate.=TRUE) }
      
      ## Take care of boundary concerns when computing H & L
      
      t.grid.boundary <- t.grid
      grid.corr <- (delta+h.HL)
      t.grid.boundary[t.grid<(t.lb+grid.corr)] <- (t.lb+grid.corr)[t.grid<(t.lb+grid.corr)]
      t.grid.boundary[t.grid>(t.ub-grid.corr)] <- (t.ub-grid.corr)[t.grid>(t.ub-grid.corr)]
      
      ## mean correcting - create demeaned vector (otherwise set to Ymi, no
      ## demeaning)
      
      if(mean.correct.HL) Ymi.demeaned <- Ymi-mu.hat[sapply(Tmi, function(x) which.min(abs(t.grid - x)))]
      
      X.hat.t1 <- NW.estimator(Ymi.demeaned,Tmi,h.HL,t.grid.boundary-delta,degenerate.method.HL,Gamma.hat,sigma.hat**2,mu.hat)$Y.hat
      X.hat.t2 <- NW.estimator(Ymi.demeaned,Tmi,h.HL,t.grid.boundary,degenerate.method.HL,Gamma.hat,sigma.hat**2,mu.hat)$Y.hat
      X.hat.t3 <- NW.estimator(Ymi.demeaned,Tmi,h.HL,t.grid.boundary+delta,degenerate.method.HL,Gamma.hat,sigma.hat**2,mu.hat)$Y.hat
      ## X.hat.mu will be adjusted for unequal sample sizes, don't do it here
      ## (in batch pooled hence correct)
      X.hat.mu <- NW.estimator(Ymi,Tmi,h.mu,t.grid,degenerate.method.mu,Gamma.hat,sigma.hat**2,mu.hat)$Y.hat
      is.na.X.hat.mu <- is.na(X.hat.mu)
      X.hat.mu.count[is.na.X.hat.mu] <- X.hat.mu.count[is.na.X.hat.mu] + 1
      ## Adjust for unequal sample sizes when computing Gamma crossproducts etc.
      X.hat.Gamma <- (Mi/Mi.mean.previous)*NW.estimator(Ymi,Tmi,h.Gamma,t.grid,degenerate.method.Gamma,Gamma.hat,sigma.hat**2,mu.hat)$Y.hat
      is.na.X.hat.Gamma <- is.na(X.hat.Gamma)
      is.na.tc.X.hat.Gamma <- is.na(tcrossprod(X.hat.Gamma,X.hat.Gamma))
      X.hat.Gamma.count[is.na.X.hat.Gamma] <- X.hat.Gamma.count[is.na.X.hat.Gamma] + 1
      X.hat.Gamma.count.mat[is.na.tc.X.hat.Gamma] <- X.hat.Gamma.count.mat[is.na.tc.X.hat.Gamma]+1
      
      ## Summand for sigma.sq involves two closest Ymi values
      
      for(j in 1:length.grid) {
        if(Mi>1) {
          close.1 <- which.min(abs(Tmi-t.grid[j]))
          sigma.sq.summand[j] <-  0.5*(Ymi[close.1]-(Ymi[-close.1])[which.min(abs(Tmi[-close.1]-t.grid[j]))])^2
        }
      }
      
      ## With new online Tmi data, updated design density estimates at grid points
      
      f.hat.hist <- hist(Tmi,plot=FALSE)
      f.hat.grid <- gamma.fullsample*f.hat.grid + (1-gamma.fullsample)*f.hat.hist$density[findInterval(t.grid,f.hat.hist$breaks,all.inside=TRUE)]
      
      ## In what follows (don't want to do this for average of Mi nor fhat
      ## above), the weights for updating theta, gamma & (1-gamma) now adjust
      ## for larger/smaller sample size in online sample i, Mi. Note with
      ## X.hat.count subtracted i.online is now a vector (compatible with all
      ## objects it is used in) reflecting skipped grid point counts (subtracts
      ## # times each grid point is skipped).
      
      gamma <- (i.online-X.hat.count)/(i.online-X.hat.count+1)*Mi.mean.previous/Mi.mean
      gamma.mu <- (i.online-X.hat.mu.count)/(i.online-X.hat.mu.count+1)*Mi.mean.previous/Mi.mean
      
      ## With the new online curve data, update estimates of the theta values
      ## needed to construct \hat H_t and \hat L. Here we skip all estimates at
      ## grid points where the kernel estimator is degenerate if
      ## degenerate.method=="skip" which produces NA values in the skipped
      ## X.hat.* element.
      
      is.na.X.hat.t1 <- is.na(X.hat.t1)
      is.na.X.hat.t2 <- is.na(X.hat.t2)
      is.na.X.hat.t3 <- is.na(X.hat.t3)
      theta.hat.t1.t3[!is.na.X.hat.t1 & !is.na.X.hat.t3] <- (gamma*theta.hat.t1.t3 + (1-gamma)*(X.hat.t1-X.hat.t3)^2)[!is.na.X.hat.t1 & !is.na.X.hat.t3]
      theta.hat.t1.t2[!is.na.X.hat.t1 & !is.na.X.hat.t2] <- (gamma*theta.hat.t1.t2 + (1-gamma)*(X.hat.t1-X.hat.t2)^2)[!is.na.X.hat.t1 & !is.na.X.hat.t2]
      theta.hat.t2.t3[!is.na.X.hat.t2 & !is.na.X.hat.t3] <- (gamma*theta.hat.t2.t3 + (1-gamma)*(X.hat.t2-X.hat.t3)^2)[!is.na.X.hat.t2 & !is.na.X.hat.t3]
      
      ## Trap zeros (feed instead .Machine$double.eps hence numerator is
      ## zero for log diff)
      
      theta.hat.t1.t3[theta.hat.t1.t3==0] <- .Machine$double.eps
      theta.hat.t1.t2[theta.hat.t1.t2==0] <- .Machine$double.eps
      theta.hat.t2.t3[theta.hat.t2.t3==0] <- .Machine$double.eps
      
      ## Update the estimate of mu, the mean curve. Again, we skip
      ## degenerate grid point estimates if degenerate.method=="skip".
      
      mu.hat[!is.na.X.hat.mu] <- (gamma.mu*mu.hat + (1-gamma.mu)*X.hat.mu)[!is.na.X.hat.mu]
      X.hat.mu.updated.mat[n,] <- X.hat.mu
      
      ## Update the estimate of sigma and H
      
      sigma.hat <- sqrt(gamma.fullsample*sigma.hat^2 + (1-gamma.fullsample)*sigma.sq.summand)
      
      H.updated.mat[n,] <- (log(theta.hat.t1.t3)-log((theta.hat.t1.t2+theta.hat.t2.t3)/2))/(2*log(2))
      
      ## Trap cases where estimated H is negative (cannot occur theoretically)
      ## OR exceeds 1
      
      H.updated.mat[n,H.updated.mat[n,] <= H.lb] <- H.lb
      H.updated.mat[n,H.updated.mat[n,] > H.ub] <- H.ub
      
      ## Note we take the mean of estimates on |t1-t2|=|t2-t3|=delta. If instead
      ## you used theta.hat.t1.t3 then you need (2*delta)^{} in the denominator
      ## where 2*delta=|t1-t3|
      
      L.sq.updated.mat[n,] <- 0.5*(theta.hat.t1.t2+theta.hat.t2.t3)/(delta^{2*H.updated.mat[n,]})
      
      ## Trap case where estimated L^2 evolves to NA
      
      L.sq.updated.mat[n,L.sq.updated.mat[n,]<=.Machine$double.eps] <- L.sq.lb
      
      ## Update cross product matrix & colSums with updated online estimates,
      ## recompute Gamma.hat. Note that X.hat.Gamma that contains NA elements at
      ## degenerate grid points when degenerate.method=="skip" (no Tmi/Ymi in
      ## +/- h.Gamma) which is used to compute mean vector via X.hat.Gamma and
      ## crossproduct matrix via tcrossprod(X.hat.Gamma,X.hat.Gamma). We know
      ## the location and count of the NAs for X.hat.Gamma from
      ## X.hat.Camma.count but must construct them for the matrix and do so via
      ## tc.X.hat.Gamma and update X.hat.Gamma.count.mat. We skip NA elements
      ## for crossprod and colSums and use the previous values of each when
      ## degenerate.method=="skip" (i.e.., Gamma.hat.crossprod.mat and
      ## Gamma.hat.colSums.vec don't get updated at elements for which the
      ## current matrix/vector to be added contains NA elements)
      
      ## Update non-degenerate crossprod and colSums entries based on tcrossprod
      ## of X.hat.Gamma and X.hat.Gamma NA entries
      
      Gamma.hat.crossprod.mat[!is.na.tc.X.hat.Gamma] <- (Gamma.hat.crossprod.mat + tcrossprod(X.hat.Gamma,X.hat.Gamma))[!is.na.tc.X.hat.Gamma]
      Gamma.hat.colSums.vec[!is.na.X.hat.Gamma] <- (Gamma.hat.colSums.vec + X.hat.Gamma)[!is.na.X.hat.Gamma]
      Gamma.hat <- Gamma.hat.crossprod.mat/(i.online-X.hat.Gamma.count.mat)-tcrossprod(Gamma.hat.colSums.vec/(i.online-X.hat.Gamma.count),Gamma.hat.colSums.vec/(i.online-X.hat.Gamma.count))
      
      ## Given Gamma.hat is now constructed, we compute the "reconstructed"
      ## (linear-projection) estimator, save the online kernel estimator as well
      
      if(compute.lp.estimator) {
        X.hat.lp[n,] <- lp.estimator(Ymi,
                                     Tmi,
                                     t.grid,
                                     Gamma.hat,
                                     sigma.hat**2,
                                     mu.hat)
      }
      
      ## Update MSE optimal h using recursive values for the unknown constants
      ## C1 and C2, H, and L, will be used as the starting value in the next
      ## recursion. When degenerate.method=="skip" we update the number of
      ## curves appearing in the bandwidth to reflect the index i.online not
      ## equalling the number of curves N+n rather it is N_n-#skipped curves for
      ## each curve (X.hat.*.count is a vector)
      
      C1.updated <- 3*L.sq.updated.mat[n,]/((2*H.updated.mat[n,]+1)*(2*H.updated.mat[n,]+3))
      C2.updated <- (3/5)*sigma.hat^2
      
      h.HL <- (C2.updated/(2*H.updated.mat[n,]*C1.updated*Mi.mean*f.hat.grid))**(1/(2*H.updated.mat[n,]+1))
      h.HL[h.HL>0.5*(t.ub-t.lb)-delta] <- 0.5*(t.ub-t.lb)-delta
      h.HL.updated.mat[n,] <- h.HL
      
      if(mu.bandwidth=="M") {
        h.mu <- h.HL
      } else if(mu.bandwidth=="NM"){
        h.mu <- (C2.updated/(2*H.updated.mat[n,]*C1.updated*(i.online-X.hat.mu.count)*Mi.mean*f.hat.grid))**(1/(2*H.updated.mat[n,]+1))
      } else if(mu.bandwidth=="logNM"){
        h.mu <- (C2.updated*log((i.online-X.hat.mu.count)*Mi.mean)/(2*H.updated.mat[n,]*C1.updated*(i.online-X.hat.mu.count)*Mi.mean*f.hat.grid))**(1/(2*H.updated.mat[n,]+1))
      }
      if(Gamma.bandwidth=="NMsq") {
        h.b1 <- (C2.updated/(4*H.updated.mat[n,]*C1.updated*(i.online-X.hat.Gamma.count)*Mi.mean*f.hat.grid))**(1/(2*H.updated.mat[n,]+1))
        h.b2 <- (C2.updated/(4*H.updated.mat[n,]*C1.updated*(i.online-X.hat.Gamma.count)*Mi.sq.mean*f.hat.grid))**(1/(2*H.updated.mat[n,]+2))
        h.Gamma <- pmax(h.b1,h.b2)
      } else if(Gamma.bandwidth=="logNMsq"){
        h.b1 <- (C2.updated*log((i.online-X.hat.Gamma.count)*Mi.mean)/(4*H.updated.mat[n,]*C1.updated*(i.online-X.hat.Gamma.count)*Mi.mean*f.hat.grid))**(1/(2*H.updated.mat[n,]+1))
        h.b2 <- (C2.updated*log((i.online-X.hat.Gamma.count)*Mi.sq.mean)/(4*H.updated.mat[n,]*C1.updated*(i.online-X.hat.Gamma.count)*Mi.sq.mean*f.hat.grid))**(1/(2*H.updated.mat[n,]+2))
        h.Gamma <- pmax(h.b1,h.b2)
      }
      
      if(plot.online.results) {
        
        # Some simple plotting of online curves which can be helpful
        
        par(mfrow=c(2,2))
        
        plot(t.grid[order(t.grid)], H.updated.mat[n,],
             xlab="grid",
             ylab="Estimated H",
             main="Estimated H",
             ylim=c(H.lb-0.1,H.ub+0.1),
             type="b",
             panel.first=grid(lty=1))
        
        ylim <- range(Ymi,(X.hat[n,])[!is.na(X.hat[n,])])
        
        plot(t.grid[order(t.grid)], X.hat[n,order(t.grid)],
             xlab="grid",
             ylab="Estimated online curve",
             ylim=ylim,
             main=paste("Estimated online curve ",i.online,sep=""),
             type="b",
             panel.first=grid(lty=1))
        points(Tmi,Ymi,cex=.5,col=2)
        legend("topleft",c("Curve","Data"),col=c(1,2),lty=c(1,NA),pch=c(NA,1),bty="n")
        
        plot(t.grid[order(t.grid)],
             mu.hat[order(t.grid)],
             xlab="grid",
             ylab="Estimated mean curve",
             main="Estimated mean curve",
             type="b",
             panel.first=grid(lty=1))
        
        GA::persp3D(t.grid[order(t.grid)],
                    t.grid[order(t.grid)],
                    Gamma.hat[order(t.grid),order(t.grid)],
                    xlab="s",
                    ylab="t",
                    zlab="Estimated Gamma",
                    main="Estimated Gamma")
        
      }
      
    }
    
  }
  
  ## Some helpful warnings related to potentially smooth curves, constant
  ## curves, and oversmoothed estimates
  
  prop.large.h <- length(which(h.HL.updated.mat[NROW(h.HL.updated.mat),]>0.9*(0.5*(t.ub-t.lb)-delta)))/length(h.HL.updated.mat[NROW(h.HL.updated.mat),])
  prop.large.H <- length(which(H.updated.mat[NROW(H.updated.mat),]>0.9))/length(H.updated.mat[NROW(H.updated.mat),])
  prop.small.L <- length(which(sqrt(L.sq.updated.mat[NROW(L.sq.updated.mat),])<0.1))/length(L.sq.updated.mat[NROW(L.sq.updated.mat),])
  if(issue.warnings & prop.large.h > 0.05) warning(paste("a substantial proportion of h_t values (",format(prop.large.h*100,digits=1),"%) are quite large (i.e., individual curves, mu, and Gamma may be oversmoothed, proceed with caution)",sep=""),immediate.=TRUE)
  if(issue.warnings & prop.large.H > 0.05) warning(paste("a substantial proportion of H_t values (",format(prop.large.H*100,digits=1),"% > 0.9) are close to 1 (i.e., curves may be continuously differentiable)",sep=""),immediate.=TRUE)
  if(issue.warnings & prop.small.L > 0.05) warning(paste("a substantial proportion of L_t values (",format(prop.small.L*100,digits=1),"% < 0.1) are close to 0 (i.e., curves may be constant functions)",sep=""),immediate.=TRUE)
  
  return(list(compute.lp.estimator=compute.lp.estimator,
              degenerate.method.HL=degenerate.method.HL,
              degenerate.method.mu=degenerate.method.mu,
              degenerate.method.Gamma=degenerate.method.Gamma,
              delta=delta,
              f.hat.grid=f.hat.grid,
              Gamma.bandwidth=Gamma.bandwidth,
              Gamma.hat.colSums.vec=Gamma.hat.colSums.vec,
              Gamma.hat.crossprod.mat=Gamma.hat.crossprod.mat,
              Gamma.hat=Gamma.hat,
              Gamma.hat.sigma.correction=Gamma.hat-diag(sigma.hat^2),
              Gamma.hat.flattop.correction=flattop.Gamma(Gamma=Gamma.hat,h.vec=h.Gamma,t.vec=t.grid),
              h.Gamma=h.Gamma,
              h.mu=h.mu,
              h.HL.updated.mat=h.HL.updated.mat,
              h.HL.updated=h.HL.updated.mat[NROW(h.HL.updated.mat),],
              h.HL=h.HL,
              h.HL.starting=h.HL.starting,
              H.updated.mat=H.updated.mat,
              H.updated=H.updated.mat[NROW(H.updated.mat),],
              issue.warnings=issue.warnings,
              length.grid=length.grid,
              L.sq.updated.mat=L.sq.updated.mat,
              L.sq.updated=L.sq.updated.mat[NROW(L.sq.updated.mat),],
              mean.correct.HL=mean.correct.HL,
              Mi.mean=Mi.mean,
              Mi.sq.mean=Mi.sq.mean,
              mu.hat.mat=X.hat.mu.updated.mat,
              mu.hat=mu.hat,
              mu.bandwidth=mu.bandwidth,
              N=i.online,
              plot.online.results=plot.online.results,
              sigma.hat=sigma.hat,
              t.grid=t.grid,
              theta.hat.t1.t2=theta.hat.t1.t2,
              theta.hat.t1.t3=theta.hat.t1.t3,
              theta.hat.t2.t3=theta.hat.t2.t3,
              X.hat=X.hat,
              X.hat.lp=X.hat.lp,
              X.hat.count=X.hat.count,
              X.hat.mu.count=X.hat.mu.count,
              X.hat.Gamma.count=X.hat.Gamma.count,
              X.hat.Gamma.count.mat=X.hat.Gamma.count.mat))
  
}

## Function for "smoothing" the diagonal band of a covariance matrix to deal
## with overlapping elements adding sigma^2 to the diagonal and nearby
## off-diagonal elements. This is an alternative to naively subtracting an
## estimate of sigma2 from just the diagonal elements

flattop.Gamma <- function(Gamma=NULL,
                          h.vec=NULL,
                          t.vec=NULL,
                          step.vector=NULL,
                          increase.step=0,
                          issue.warnings=FALSE) {
  
  if(is.null(Gamma)) stop("Must input Gamma")
  if(is.null(h.vec)) stop("Must input h.vec")
  if(is.null(t.vec)) stop("Must input t.vec")
  
  K <- NROW(Gamma)
  
  ## Stepping off the diagonal at 90 degrees (anti diagonal) limits the number
  ## of steps to at most dim(Gamma)/2 which, since Gamma is square is
  ## nrow(Gamma)/2
  
  K.one.half <- round(K/2)
  
  ## Allow the user to pass in an integer 1,2,... K/2-1 that widens the diagonal
  ## "flattop" beyond what the bandwidths and grid suggests
  
  if(increase.step < 0 || increase.step >= K.one.half) stop("increase.step must be a non-negative integer less than 1/2 the dimension of Gamma")
  
  ## Storage vector for number of steps taken off main diagonal k for each
  ## bandwidth/grid combination
  
  l.vec <- numeric()
  
  ## Since we take the anti diagonal band 0 and the one to its "upper left",
  ## band 1, we start from the lower right corner of Gamma and proceed (hence "k
  ## in K:1)
  
  seq.1.K.one.half <- 1:K.one.half
  
  for(k in K:1) {
    
    ## Adaptive l, step-wise blocking with anti diagonal band 0 and band 1.
    ## First, identify the appropriate covariance to copy into the bands (first
    ## appearance of a |t_i-t_j|>h_i+h_j with no overlap, minimum 1 so minimum
    ## band 0 has 3 anti diagonal elements, band 1 2). The integer l the becomes
    ## the offset from the diagonal of Gamma for a given k.
    
    R.idx <- pmax(1,k-seq.1.K.one.half)
    C.idx <- pmin(k+seq.1.K.one.half,K)
    
    ## L is a vector of logicals, TRUE = no overlap (distance between grid
    ## points exceeds sum of bandwidths), FALSE = overlap. Presuming both TRUE
    ## and FALSE values exist, obtain integer l representing how many steps away
    ## from the diagonal the first non-overlapping grid point lies. If all
    ## points overlap (all FALSE) l is the max value possible K.one.half.
    
    L <- (abs(t.vec[C.idx]-t.vec[R.idx])>=(h.vec[C.idx]+h.vec[R.idx]))
    if(any(L)) {l <- min(which(L)) } else {l <- K.one.half}
    
    ## Can manually increase l by a set number of steps if desired.
    
    if(increase.step) l <- l + increase.step
    
    ## Allow the user to instead provide a step vector (one for each k) which
    ## can contain different integers for each k, say, mean values from previous
    ## recursions on smaller amounts of data (will override increase.step)
    
    if(!is.null(step.vector)) l <- step.vector[k]
    
    l.vec[k] <- l
    
    ## Next, identify band 0 and band 1 indices and ensure they are valid matrix
    ## entries
    
    I <- (k-l):(k+l)
    I <- I[I > 0 & I <= K]
    I.band.0 <- cbind((k+l):(k-l),(k-l):(k+l))
    I.band.0 <- I.band.0[apply(I.band.0, 1, function(row) all(row > 0 & row <= K )),,drop=FALSE]
    I.band.1 <- cbind(c((k+(l-1)):k,(k-1):(k-l)),c((k-l):(k-1),k:(k+(l-1))))
    I.band.1 <- I.band.1[apply(I.band.1, 1, function(row) all(row > 0 & row <= K )),,drop=FALSE]
    
    ## Finally, copy the appropriate covariance entry into anti diagonal band 0
    ## and band 1 of the kth (2l+1)x(2l+1) matrix in Gamma centred on (k,k)
    
    Gamma[rbind(I.band.0,I.band.1)] <- Gamma[cbind(R.idx,C.idx)][l]
    
  }
  
  ## Assign attribute to the l vector of offsets from the diagonal k
  
  attr(Gamma,"k.offset") <- l.vec
  
  if(issue.warnings & max(l.vec)==1) warning("max. diagonal offset is 1, consider a denser grid",immediate.=TRUE)
  
  return(Gamma)
  
}

colMedians <- function(data) {
  data[which(sapply(data, is.infinite),arr.ind=TRUE)] <- NA
  ## Remove runs with NAs
  data <- na.omit(data)
  colmed <- numeric(ncol(data))
  for(i in 1:ncol(data)) {
    colmed[i] <- median(data[,i])
  }
  return(colmed)
}

## Global pooled smooth using flexible smoother, 3 choices

mu.pooled.func <- function(curves,
                           grid,
                           method=c("np","crs","loess")){
  
  ## One-dimensional smooth on all curves (pooled data) evaluated on grid points
  ## - note loess() can produce NA values if grid points lie beyond sample, note
  ## ... gets passed to npreg only
  
  method <- match.arg(method)
  N <- length(curves)
  TN <- NULL
  YN <- NULL
  for(n in 1:N) {
    TN <- c(TN,curves[[n]]$t)
    YN <- c(YN,curves[[n]]$x)
  }
  if(method=="crs") {
    suppressPackageStartupMessages(require(crs))
    options(crs.messages=FALSE)
    as.numeric(predict(crs(YN~TN,cv="exhaustive"),
                       newdata=data.frame(TN=grid)))
  } else if(method=="loess") {
    predict(loess(YN~TN),
            newdata=grid)
  } else if(method=="np") {
    suppressPackageStartupMessages(require(np))
    options(np.messages=FALSE,
            np.tree=TRUE)
    ## No need to anal-retentive precision on optimization, override
    ## anal-retentive search precision defaults
    predict(npreg(YN~TN,ckertype="epanechnikov",tol=1e-04,ftol=1e-04,nmulti=1),
            newdata=data.frame(TN=grid))
  }
}

## Some plotting routines to be deployed on a set of curves

plot.true <- function(curves) {
  N <- length(curves)
  ylim <- range(curves[[1]]$x_true)
  if(N>1) {
    for(n in 2:N) {
      ylim <- range(ylim,curves[[n]]$x_true)
    }
  }
  plot(curves[[1]]$grid_true[order(curves[[1]]$grid_true)],curves[[1]]$x_true[order(curves[[1]]$grid_true)],type="l",ylim=ylim,xlab="t",ylab="True Curves",panel.first=grid(lty=1))
  if(N>1) {
    for(n in 2:N) {
      lines(curves[[n]]$grid_true[order(curves[[n]]$grid_true)],curves[[n]]$x_true[order(curves[[n]]$grid_true)],col=n,lty=n)
    }
  }
}

plot.data <- function(curves) {
  N <- length(curves)
  ylim <- range(curves[[1]]$x)
  if(N>1) {
    for(n in 2:N) {
      ylim <- range(ylim,curves[[n]]$x)
    }
  }
  plot(curves[[1]]$t[order(curves[[1]]$t)],curves[[1]]$x[order(curves[[1]]$t)],type="b",ylim=ylim,xlab="Index",ylab="Data",cex=.5,panel.first=grid(lty=1))
  if(N>1) {
    for(n in 2:N) {
      lines(curves[[n]]$t[order(curves[[n]]$t)],curves[[n]]$x[order(curves[[n]]$t)],col=n,lty=n,type="b",cex=.5)
    }
  }
}

summary.data <- function(curves,kable=FALSE) {
  N <- length(curves)
  M <- length(curves[[1]]$x)
  ylim <- range(curves[[1]]$x)
  for(n in 2:N) {
    M <- c(M,length(curves[[n]]$x))
    ylim <- range(ylim,curves[[n]]$x)
  }
  if(kable) {
    kableExtra::footnote(knitr::kable(table(M),caption="Summary of curve object (M = obs. per curve, Freq = no. curves having M obs.)"),
                         alphabet=c(paste("Number of curves: ",N,sep=""),
                                    paste("Mean sample size: ",format(mean(M),digits=2),sep=""),
                                    paste("Total number of observations by curve: ",sum(M),sep="")))
  } else {
    print(paste("Number of curves: ",N,sep=""))
    print(paste("Mean sample size: ",format(mean(M),digits=2),sep=""))
    print(paste("Total number of observations: ",sum(M),sep=""))
    print("Summary of unique number of observations by curve:")
    print(table(M,dnn=c("M")))
  }
  
  
}

## The epanechnikov kernel

ker.epa <- function(t,Ti,h){
  u <- (t-Ti)/h
  ifelse(abs(u)<=1,0.75*(1-u**2),.Machine$double.eps)
}

## Nadaraya-Watson (local constant) estimator that if degenerate (i.e., no local
## obs in t+- h) returns the closest Yi,Ti to t0 if method=="1nn", NA if
## method="skip", and the reconstructed estimator if method="lp", along with the
## vector of fitted values, vector of logicals (whether degeneracy occurred at
## an evaluation point or not), and a vector of local observations in +-h of
## each evaluation point.

NW.estimator <- function(Yi=NULL,
                         Ti=NULL,
                         h=NULL,
                         t=NULL,
                         method=c("skip","1nn","lp"),
                         Gamma=NULL,
                         sigma.sq=NULL,
                         mu=NULL) {
  
  ## Some basic input checking
  
  method <- match.arg(method)
  if(is.null(Yi)) stop("must pass curve vector Yi")
  if(is.null(Ti)) stop("must pass curve vector Ti")
  if(is.null(h)) stop("must pass bandwidth vector h")
  if(is.null(t)) stop("must pass grid vector t")
  
  degenerate.logical <- logical()
  Y.hat <- numeric()
  L <- numeric()
  
  ## Reconstructed estimator needs to invert Gamma + Sigma, and match closest
  ## values on grid with sample Tmi. There is also an issue of multiple matches
  ## which appears to require the corresponding value of sigma to increase.
  
  if(method=="lp") {
    if(is.null(Gamma)) stop("must pass Gamma")
    if(is.null(sigma.sq)) stop("must pass sigma.sq")
    if(is.null(mu)) stop("must pass mu")
    idx <- sapply(Ti, function(x) which.min(abs(t - x)))
    ## Use of count.idx correction here damages RMSE
    ## significantly. Use of count.idx in lp.estimator() significantly
    ## improves RMSE.
    G.S.inv <- solve(Gamma[idx,idx]+diag(sigma.sq[idx]))
  }
  
  for(j in 1:length(t)) {
    
    ## Compute kernel vector for point t[j] in the grid vector t
    
    ker.mu <- ker.epa(t[j],Ti,h[j])
    
    ## Check whether grid point t[j] is a degenerate case, i.e., there are zero
    ## observations in +- h[j] of t[j] (happens when all kernel values are equal
    ## to the value assigned by ker.eps() when all sample Ti lie outside t[j]
    ## +/- h[j], which  is .Machine$double.eps). If the case is degenerate,
    ## handle the case using method=="..."
    
    if((L[j] <- length(ker.mu[ker.mu !=.Machine$double.eps]))==0) {
      if(method=="1nn") {
        ## Degenerate case returns closest Yi (1st nearest neighbour) if no obs
        ## in +- h of t
        Y.hat[j] <- Yi[which.min(abs(Ti-t[j]))]
      } else if(method=="skip") {
        ## Degenerate case returns NA if no obs in +- h of t
        Y.hat[j] <- NA
      } else if(method=="lp") {
        ## Degenerate case returns reconstructed estimator based on estimate of
        ## Gamma if no obs in +- h of t
        Y.hat[j] <- mu[j] + (Yi - mu[idx]) %*% G.S.inv %*% Gamma[j,idx]
      }
      degenerate.logical[j] <- TRUE
    } else {
      ## Non-degenerate case
      Y.hat[j] <- sum(Yi*ker.mu)/sum(ker.mu)
      degenerate.logical[j] <- FALSE
    }
  }
  
  return(list(Y.hat=Y.hat,degenerate=degenerate.logical,local.obs=L))
  
}

## Reconstructed curve using "snap-to-grid" (knn=1) approach.

lp.estimator <- function(Yi=NULL,
                         Ti=NULL,
                         t=NULL,
                         Gamma=NULL,
                         sigma.sq=NULL,
                         mu=NULL) {
  
  ## Some basic input checking
  
  if(is.null(Yi)) stop("must pass curve vector Yi")
  if(is.null(Ti)) stop("must pass curve vector Ti")
  if(is.null(t)) stop("must pass grid vector t")
  if(is.null(Gamma)) stop("must pass Gamma")
  if(is.null(sigma.sq)) stop("must pass sigma.sq")
  if(is.null(mu)) stop("must pass mu")
  
  if(length(Yi) < 3) stop("linear projection for a curve must be based on at least 3 sample points:\n    1. restart using option compute.lp.estimator=FALSE, or\n    2. discard any such curves")
  
  ## First find index of grid values (t) closest to the Ti for the sample,
  ## correct for presence of multiple repeated indices, then compute the
  ## reconstructed curve using the 1nn Gamma values.
  
  idx <- sapply(Ti, function(x) which.min(abs(t - x)))
  count.idx <- sapply(idx, function(x) sum(idx==x))
  G.S.inv <- solve(Gamma[idx,idx]+diag(count.idx*sigma.sq[idx]))
  
  X.hat.lp <- numeric()
  
  for(j in 1:length(t)) {
    X.hat.lp[j] <- mu[j] + (Yi - mu[idx]) %*% G.S.inv %*% Gamma[j,idx]
  }
  
  return(X.hat.lp)
  
}

############# Sunny's code - Sunny Wang <sunnywanggw@gmail.com> ################

################################################################################
### Data generation process
###
### Generation of MFBM
### See the following paper for the covariance structure
### https://doi.org/10.3390/fractalfract6020074
### Equations (6) and (7)
################################################################################

# Load libraries

##
## Usual Hurst functions
##

hurst_atan <- function(t_vec) {
  atan(t_vec) / pi + 0.5
}

hurst_linear <- function(t_vec, h_left = 0.2, h_right = 0.8) {
  t1 <- max(t_vec)
  t0 <- min(t_vec)
  a <- (h_right - h_left) / (t1 - t0)
  b <- h_right - a * t1
  pmin(a * t_vec + b, 1)
}

## Supplants function defined in library/dgp_old.R, modified by Sunny but can be
## overwritten by code below depending on where it is place, so commented out
## code below to avoid issues such as those I encountered...

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

#hurst_logistic <- function(
    #    t_vec,
#    h_left = 0.5,
#    h_right = 0.8,
#    slope = 30,
#    change_point_position = 0.5
#) {
#    change_point <- change_point_position * (max(t_vec) + min(t_vec))
#    u <- (t_vec - change_point) / (max(t_vec) - min(t_vec))
#    (h_right - h_left) / (1 + exp(-slope * u)) + h_left
#}

# not totally in agreement with the theory
hurst_piecewise <- function(
    t_vec,
    h_left = 0.2,
    h_right = 0.8,
    change_point_position = 0.5
) {
  change_point <- change_point_position * (max(t_vec) + min(t_vec))
  h_left * (t_vec <= change_point) + h_right * (t_vec > change_point)
}

##
## Covariance structure
##

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
  #    grid <- seq(t_min, t_max, length.out = add_regular_grid_of_size)
  covariance_list <- lapply(
    points_list,
    function(point) {
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

##
## Class "curves_ideal_observed": methods for generic functions
##

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

