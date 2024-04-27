\name{fda.la}
\alias{fda.la}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Locally Adaptive Functional Data Analysis}
\description{
\code{fda.la()} uses a variety of nonparametric methods to adaptively learn key unknown local parameters for a class of irregular functions in a Functional Data Analysis (FDA) framework. Local kernel smoothing methods are used, and key unknown parameters in the pointwise-optimal bandwidth formulas are learned adaptively from supplied data. The class of irregular functions lie in the so-called local Hölder class (the typical Hölder class \eqn{f:\mathbb{R}\to\mathbb{R}}, and \eqn{f} is Hölder continuous if there exist constant \eqn{C>0} and exponent \eqn{\alpha>0} such that \eqn{|f(x)-f(y)|\le C|x-y|^\alpha} for all \eqn{x} and \eqn{y} - note that when \eqn{\alpha=1} we have the class of Lipschitz functions). The approach is fully data-driven and automatic, requiring nothing more than an object containing a set of functional data.

Smoothing adapts to the unknown regularity of the object being estimated which is characterized by \eqn{H_t} and \eqn{L_t}, the regularity parameters associated with the fractional derivative \eqn{(X_u-X_v)/|u-v|^{H_t}} where \eqn{\mathbb{E}((X_u-X_v)^2)\simeq L_t^2|u-v|^{2H_t}} and where \eqn{X} is functional data. 

In particular, for a pointwise-optimal bandwidth for estimating the unknown \eqn{H_t}, we use data-driven methods for obtaining key constants based on the fractional derivative and kernel function adopted (here the Epanechnikov kernel) which deliver estimates of \eqn{H_t} and \eqn{L_t} (these are obtained by a convergent iterative procedure). These same constants appear in the pointwise MSE-optimal bandwidths for estimating a) each online functional curve \eqn{X^{(i)}}, \eqn{i=1,\dots,N}, b) the mean functional curve \eqn{\mu(t)}, and c) the covariance matrix for the online functional curves \eqn{\Gamma(s,t)}. The difference between the smoothing parameters used for estimating \eqn{H_t} and \eqn{L_t}, on one hand, and \eqn{\mu(t)} and \eqn{\Gamma(s,t)}, on the other, are the factors of proportionality and rate of convergence. The factors of proportionality can include \eqn{1/M_i}, the number of observations for the \eqn{i}th online curve, \eqn{1/\sum_{i=1}^N M_i}, or \eqn{\log(\sum_{i=1}^N M_i)/\sum_{i=1}^N M_i}, or \eqn{1/N} where \eqn{N} is the number of curves at hand, while the rates of convergence can include \eqn{1/(2H_t+1)} or \eqn{1/(2H_t+2)} (the former for estimating \eqn{X^{(i)}} and \eqn{\mu(t)}, both the former and latter for estimating \eqn{\Gamma(s,t)}).

This function can be used for optimal smoothing for a batch of curves simultaneously, or for optimal smoothing of a sequence of new online curves using a dynamic recursive updating algorithm. When used for estimating newly acquired online curve data, the method is computationally efficient from both CPU time and memory requirements.
}
\usage{
fda.la(common.design = FALSE, 
       compute.lp.estimator = TRUE,
       curves_batch = NULL, 
       curves_online = NULL, 
       delta = NULL, 
       degenerate.method.HL=c("skip","1nn","lp"),
       degenerate.method.mu=c("skip","1nn","lp"),
       degenerate.method.Gamma=c("skip","1nn","lp"),
       f.hat.grid = NULL, 
       fda.la.previous = NULL,
       Gamma.bandwidth = c("logNMsq", "NMsq"), 
       Gamma.hat.colSums.vec = NULL, 
       Gamma.hat.crossprod.mat = NULL, 
       Gamma.hat = NULL, 
       h.Gamma = NULL, 
       h.mu = NULL, 
       h.HL.starting = NULL, 
       h.HL = NULL, 
       H.lb = 0.1, 
       H.ub = 1, 
       issue.warnings = TRUE, 
       L.sq.lb = .Machine$double.eps, 
       mean.correct.HL=TRUE,
       Mi.mean = NULL, 
       Mi.sq.mean = NULL
       mu.bandwidth = c("logNM", "M", "NM"), 
       mu.hat = NULL, 
       mu.method = c("default", "Bspline.cv", "kernel.cv"), 
       N = NULL, 
       num.iter = 25, 
       plot.batch.results = FALSE, 
       plot.online.results = FALSE, 
       sigma.hat = NULL, 
       t.grid = NULL, 
       t.lb = 0, 
       t.ub = 1, 
       theta.hat.t1.t2 = NULL, 
       theta.hat.t1.t3 = NULL, 
       theta.hat.t2.t3 = NULL,
       X.hat.count=NULL,
       X.hat.mu.count=NULL,
       X.hat.Gamma.count=NULL,
       X.hat.Gamma.count.mat=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{common.design}{
A logical value indicating whether the curves have common design (\code{common.design=TRUE} or independent design \code{common.design=FALSE})
}
  \item{compute.lp.estimator}{
A logical value that determines whether to compute the linear projection estimator for each curve using the estimated \eqn{\Gamma(s,t)} and sample points for each curve
}
  \item{curves_batch}{
A list, each element containing elements \code{curves_batch[[i]]$t} and \code{curves_batch[[i]]$x}, \eqn{i=1,...,N}, where \eqn{N} is the number of curves in the batch, \code{t} the design points for the \eqn{i}th curve, and \code{x} the sample points
}
  \item{curves_online}{
A list, each element containing elements \code{curves_online[[i]]$t} and \code{curves_online[[i]]$x}, \eqn{i=1,...,N.online}, where \eqn{N.online} is the number of online curves, \code{t} the design points for the \eqn{i}th curve, and \code{x} the sample points
}
  \item{delta}{
A parameter used to estimate \eqn{H_t} that defaults to \eqn{\Delta^*/4} where \eqn{\Delta=e^{-\log(M)^{1/3}}}
}
  \item{degenerate.method.HL}{
A character string that determines how degenerate points (i.e., no observations in \eqn{t\pm h}) are handled when computing \eqn{H_t} and \eqn{L_t} for the online curves: \code{degenerate.method="1nn"} indicates 1st nearest neighbor, \code{degenerate.method="skip"} indicates to discard the online curve's contribution at each degenerate point (in the "vertical" sense) and return \code{NA}, and code{degenerate.method="lp"} indicates to use a linear curve reconstruction method based on the estimated \eqn{\Gamma(s,t)}
}
  \item{degenerate.method.mu}{
A character string that determines how degenerate points (i.e., no observations in \eqn{t\pm h}) are handled when computing \eqn{\mu(t)} for the online curves: \code{degenerate.method="1nn"} indicates 1st nearest neighbor, \code{degenerate.method="skip"} indicates to discard the online curve's contribution at each degenerate point (in the "vertical" sense) and return \code{NA}, and code{degenerate.method="lp"} indicates to use a linear curve reconstruction method based on the estimated \eqn{\Gamma(s,t)}
}
  \item{degenerate.method.Gamma}{
A character string that determines how degenerate points (i.e., no observations in \eqn{t\pm h}) are handled when computing \eqn{\Gamma(s,t)} for the online curves: \code{degenerate.method="1nn"} indicates 1st nearest neighbor, \code{degenerate.method="skip"} indicates to discard the online curve's contribution at each degenerate point (in the "vertical" sense) and return \code{NA}, and code{degenerate.method="lp"} indicates to use a linear curve reconstruction method based on the estimated \eqn{\Gamma(s,t)}
}
  \item{f.hat.grid}{
When called recursively, the density estimates for the design points \eqn{t} evaluated on the grid of points \eqn{t_0}
}
\item{fda.la.previous}{
When called recursively, output from previous invocation of \code{fda.la} stored in object via 

\code{fda.la.previous <- fda.la(...)} 

is used to retrieve all necessary objects for performing recursion, alternately one can provide them manually
}
  \item{Gamma.bandwidth}{
When called recursively, the vector of bandwidths used to compute \eqn{\Gamma(s,t)} evaluated on the grid of points \eqn{t_0}
}
  \item{Gamma.hat.colSums.vec}{
When called recursively, the vector of column sums used to compute \eqn{\Gamma(s,t)} evaluated on the grid of points \eqn{t_0}
}
  \item{Gamma.hat.crossprod.mat}{
When called recursively, the matrix of cross product terms used to compute \eqn{\Gamma(s,t)} evaluated on the grid of points \eqn{t_0}
}
  \item{Gamma.hat}{
The estimate of the covariance matrix \eqn{\Gamma(s,t)} evaluated on the grid of points \eqn{t_0}
}
  \item{h.Gamma}{
The vector of bandwidths used to compute the covariance matrix \eqn{\Gamma(s,t)} evaluated on the \eqn{2\times 2} matrix of grid points \eqn{t_0}
}
  \item{h.mu}{
The vector of bandwidths used to compute the mean function \eqn{\mu(t)} evaluated on the grid of points \eqn{t_0}
}
  \item{h.HL.starting}{
  A scalar for uninitialized bandwidths used to generate a vector of starting values to compute the regularization vectors \eqn{H_t} and \eqn{L_t} evaluated on the grid of points \eqn{t_0}
}
  \item{h.HL}{
The vector of bandwidths used to compute the regularization vectors \eqn{H_t} and \eqn{L_t} evaluated on the grid of points \eqn{t_0}
}
  \item{H.lb}{
Lower bound (scalar) provided for the regularization vector \eqn{H_t} evaluated on the grid of points \eqn{t_0}
}
  \item{H.ub}{
Upper bound (scalar) provided for the regularization vector \eqn{H_t} evaluated on the grid of points \eqn{t_0}
}
  \item{issue.warnings}{
Logical that determines whether to issue warnings as they occur or not (\code{TRUE} issues them immediately)
}
  \item{L.sq.lb}{
Lower bound (scalar) provided for the regularization vector \eqn{L^2_t} evaluated on the grid of points \eqn{t_0}
}
  \item{mean.correct.HL}{
A logical value that determines whether to mean correct individual curves used to construct \eqn{H_t} and \eqn{L_t}
}
  \item{Mi.mean}{
When called recursively, the mean sample size for all curves called previously, but also returned by invocation of the function
}
  \item{Mi.sq.mean}{
When called recursively, the mean of the squared sample size for all curves called previously, but also returned by invocation of the function
}
  \item{mu.bandwidth}{
The vector of bandwidths used to compute the mean function \eqn{\mu(t)} evaluated on the grid of points \eqn{t_0}
}
  \item{mu.hat}{
The estimated mean function \eqn{\mu(t)} evaluated on the grid of points \eqn{t_0}
}
  \item{mu.method}{
When invoked on a batch of curves, whether to use the default kernel smoothing method to estimate the mean function \eqn{\mu(t)} or, alternatively, to use cross-validated B-splines or cross-validated kernel smoothing with a global (scalar) smoothing parameter 
}
  \item{N}{
The number of curves involved from previous invocations of the function (used to inform subsequent online curve estimation)
}
  \item{num.iter}{
An integer provided that determines the number of iterations of the batch of curves used to determine the vector of smoothing parameters \code{h.HL} for determining the regularization parameter vectors \eqn{H_t} and \eqn{L_t} that are evaluated on the grid of points \eqn{t_0}
}
  \item{plot.batch.results}{
A logical value that determines whether to produce some simple graphical summaries when the function is invoked on a batch of curves
}
  \item{plot.online.results}{
A logical value that determines whether to produce some simple graphical summaries when the function is invoked on a set of online curves
}
  \item{sigma.hat}{
When called recursively, the estimated vector of \eqn{\sigma_t} evaluated on the grid of points \eqn{t_0}
}
  \item{t.grid}{
The grid of points \eqn{t_0} provided by the user
}
  \item{t.lb}{
The lower bound of the support of the online and batch curves \eqn{t_{lb}} (if not provided will use minimum of the empirical support)
}
  \item{t.ub}{
The upper bound of the support of the online and batch curves \eqn{t_{ub}} (if not provided will use maximum of the empirical support)
}
  \item{theta.hat.t1.t2}{
When called recursively, a set of estimates used to compute the regularization parameters \eqn{H_t} and \eqn{L_t}
}
  \item{theta.hat.t1.t3}{
When called recursively, a set of estimates used to compute the regularization parameters \eqn{H_t} and \eqn{L_t}
}
  \item{theta.hat.t2.t3}{
When called recursively, a set of estimates used to compute the regularization parameters \eqn{H_t} and \eqn{L_t}
}
  \item{X.hat.count}{
  A vector of counts representing the number of skipped curves (degenerate cases poitwise) at each grid point for the computation of \eqn{H_t} and \eqn{L_t}, which is non-zero when \code{degenerate.method.HT="skip"}
}
  \item{X.hat.mu.count}{
  A vector of counts representing the number of skipped curves (degenerate cases poitwise) at each grid point for the computation of \eqn{\mu(t)}, which is non-zero when \code{degenerate.method.mu="skip"}
}
  \item{X.hat.Gamma.count}{
  A vector of counts representing the number of skipped curves (degenerate cases poitwise) at each grid point for the computation of the mean vector used to compute \eqn{\Gamma(s,t)}, which is non-zero when \code{degenerate.method.Gamma="skip"}
}
  \item{X.hat.Gamma.count.mat}{
  A matrix of counts representing the number of skipped curves (degenerate cases pointwise) at each grid point for the computation of the crossproduct vector used to compute \eqn{\Gamma(s,t)}, which is non-zero when \code{degenerate.method.Gamma="skip"}
}
}
\details{
  This function contains a variety of optional arguments which allow for a range of tasks to be accomplished, and there are two ways to think about using this function depending on whether you are interested in conducting analysis on a \dQuote{batch} of curves or, instead, conducting \dQuote{online} analysis of one (or more) newly arrived curves having already processed a batch of curves. 
  
  In batch mode, a set of \eqn{N} existing curves are used to compute various objects. The process is self-starting and uses sensible defaults for all objects (bandwidths etc.) save for the \dQuote{mean curve} (it is recommended that the user provide some initial mean curve based on the pooled data, or use the function \code{mu.pooled.func()} - see the example below for an illustration). In batch mode all data is used, i.e., for curves with no sample points within a bandwidth's reach of a grid point the first nearest neighbour is used (\code{degenerate.method.*="1nn"}) rather than, say, discarding that curve's information. It is intended that batch mode is used for a relatively small number of existing curves where batch computation can be done within the existing memory footprint of your computer.
  
  Essentially, batch mode proceeds from arbitrary starting values for all objects and then iterates the process to produce counterparts to the arbitrary starting values that are obtained from optimally smoothed curves where the optimal bandwidths are obtained via an interative procedure that converges quite quickly (after roughly 10 or so iterations, experience would indicate).
  
  Typically we envision that batch mode is used for feasibly computing a set of existing curves that will be used for subsequent updating as, for instance, when new curves become available either one at a time or as a set.
  
    The salient difference between the two modes of calling this function is that in online mode one takes the results from the initial batch run and updates these results using recursion based on a newer set of curves drawn from the same stochastic process as that for the original batch. This way, one can call the function each time a new curve arises and recursively update the prior estimates. In this manner, the estimates are refined and as the number of curves processed increases the estimates will themselves stabilize. Alternatively, one could take a set of new curves and process them using online mode, and they will be treated recursively as if each curve was appearing in succession. This method of processing is computationally expedient requiring minimal memory and computation even when when the number of previously processed curves increases without limit.
  
  Some very basic plotting can be done for batch or online invocation with plots being specific to each.

}
\value{

  This function returns...
  
  \item{ degenerate.method.HL }{ per above }
  \item{ degenerate.method.mu }{ per above }
  \item{ degenerate.method.Gamma }{ per above }
  \item{ delta }{ per above }
  \item{ f.hat.grid }{ per above post recursion }
  \item{ Gamma.bandwidth }{ per above post recursion }
  \item{ Gamma.hat.colSums.vec }{ per above post recursion }
  \item{ Gamma.hat.crossprod.mat }{ per above post recursion }
  \item{ Gamma.hat }{ per above post recursion }
  \item{ Gamma.hat.sigma.correction }{ Gamma.hat with diagonal corrected by substraction of sigma.hat estimate post recursion }
  \item{ Gamma.hat.flattop.correction }{ Gamma.hat with diagonal band corrected by flattop method post recursion }
  \item{ h.Gamma }{ per above post recursion }
  \item{ h.mu }{ per above post recursion }
  \item{ h.HL.updated.mat }{ matrix of updated values of h.HL post-recursion over all recursions }
  \item{ h.HL.updated }{ updated value of h.HL post-recursion for last recursion }
  \item{ h.HL }{ per above post recursion }
  \item{ h.HL.starting }{ per above }
  \item{ H.updated.mat }{ matrix of updated values of H post-recursion over all recursions }
  \item{ H.updated }{ updated value of H post-recursion for last recursion }
  \item{ issue.warnings }{ per above }
  \item{ length.grid }{ per above }
  \item{ L.sq.updated.mat }{ matrix of updated values of \eqn{L^2} post-recursion over all recursions }
  \item{ L.sq.updated }{ updated value of \eqn{L^2} post-recursion for last recursion }
  \item{ mean.correct.HL }{ per above }
  \item{ Mi.mean }{ per above post recursion }
  \item{ Mi.sq.mean }{ per above post recursion }
  \item{ mu.hat.mat }{ per above post-recursion over all recursions }
  \item{ mu.hat }{ per above post recursion }
  \item{ mu.bandwidth }{ per above post recursion }
  \item{ N }{ per above }
  \item{ plot.online.results }{ per above }
  \item{ sigma.hat }{ per above post recursion }
  \item{ t.grid }{ per above }
  \item{ theta.hat.t1.t2 }{ per above post recursion }
  \item{ theta.hat.t1.t3 }{ per above post recursion }
  \item{ theta.hat.t2.t3 }{ per above post recursion }
  \item{ X.hat }{ estimated curve matrix }
  \item{ X.hat.lp }{ estimated curve matrix using lp reconstruction }
  \item{ X.hat.count }{ per above post recursion }
  \item{ X.hat.mu.count }{ per above post recursion }
  \item{ X.hat.Gamma.count }{ per above post recursion }
  \item{ X.hat.Gamma.count.mat }{ per above post recursion }

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
 Valentin Patilea \email{valentin.patilea@gmail.com}, Jeffrey S. Racine \email{racinej@mcmaster.ca}
}
\note{
  There are some functions that operate on the list of curves provided: \code{plot.data(curves)}, \code{plot.curves(curves)} (for simulated data containing \code{curves$x_true} and \code{curves$grid_true}), and \code{summary.data(curves,kable=TRUE/FALSE)} that provides summary information on the curves with regards to sample points etc.

  Note that the option \code{common.design} triggers...
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Require libraries for 3D plotting and multivariate random number generation

require(GA)
require(MASS)

## Set parameters here

## Regularity parameters for generating data

h_first = 0.25
h_second = 0.75
h_slope = 10
change_point_t = 0.5

## Parameters for generating data from a Brownian motion

sigma = 0.25
tau = 1
L = 1

## Set size of batch (N), (average) number of observations for each online curve
## (M), number of online curves, and grid of points on which to estimate curves

N = 100
M = 100
N.online = 1
t0.grid = seq(0,1,length=25)

## Set the mean function

mu.func <- function(x) { sin(2*pi*x) }

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

## Call the function on the batch using a preliminary estimator of the pooled mean
## (the preliminary estimator is replaced by the proposed estimator of mu(t), it is
## only used for computing H_t and L_t for the batch)

mu.starting <- mu.pooled.func(curves_batch,t0.grid,method="np")

fda.la.out <- fda.la(curves_batch=curves_batch,
                     t.grid=t0.grid,
                     mu.hat=mu.starting,
                     plot.batch.results=TRUE,
                     plot.online.results=TRUE)

## Generate online data

points_dist <- functional::Curry(runif, min = 0, max = 1)

points_list <- generates_increasing_length_points(N = N.online, m = M,
                                                  distribution = points_dist)

online <- generate_curves(points_list,
                          hurst = hurst_logistic,
                          grid = t0.grid,
                          sigma = sigma,
                          tau = tau,
                          L = L)

curves_online <- lapply(online$curves,
                        function(x) list(t = x$observed$t,
                                         x = x$observed$x + mu.func(x$observed$t),
                                         grid_true = x$ideal$t))

## Call the function on the online data (call with plot.online.results=TRUE to
## create a simple plot that overwrites the batch plot generated above)

fda.la.out <- fda.la(curves_online=curves_online,
                     fda.la.previous=fda.la.out)
                                  
## Or you can call the function on the online data and pass all required objects 
## instead of simply inputting fda.la.out (previous invocation on batch)

## fda.la.out <- fda.la(curves_online=curves_online, 
##                      plot.online.results=TRUE,
##                      t.grid=t0.grid,
##                      delta=fda.la.out$delta,
##                      N=fda.la.out$N,
##                      mu.hat=fda.la.out$mu.hat,
##                      Gamma.hat=fda.la.out$Gamma.hat,
##                      Gamma.hat.crossprod.mat=fda.la.out$Gamma.hat.crossprod.mat,
##                      Gamma.hat.colSums.vec=fda.la.out$Gamma.hat.colSums.vec,
##                      Mi.mean=fda.la.out$Mi.mean,
##                      Mi.sq.mean=fda.la.out$Mi.sq.mean,
##                      sigma.hat=fda.la.out$sigma.hat,
##                      f.hat.grid=fda.la.out$f.hat.grid,
##                      h.HL=fda.la.out$h.HL,
##                      h.mu=fda.la.out$h.mu,
##                      h.Gamma=fda.la.out$h.Gamma,
##                      theta.hat.t1.t3=fda.la.out$theta.hat.t1.t3,
##                      theta.hat.t1.t2=fda.la.out$theta.hat.t1.t2,
##                      theta.hat.t2.t3=fda.la.out$theta.hat.t2.t3,
##                      X.hat.count=fda.la.out$X.hat.count,
##                      X.hat.mu.count=fda.la.out$X.hat.mu.count,
##                      X.hat.Gamma.count=fda.la.out$X.hat.Gamma.count, 
##                      X.hat.Gamma.count.mat=fda.la.out$X.hat.Gamma.count.mat) 

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
