\name{dream}
\alias{dream}
\title{DiffeRential Evolution Adaptive Metropolis (DREAM)}
\usage{dream(FUN, func.type, pars, FUN.pars=list(), INIT=LHSInit,
INIT.pars=list(), control=list(), measurement)}
\description{
  Efficient global MCMC even in high-dimensional spaces.
  From J.A. Vrugt, C.J.F. ter Braak et al.
}
\arguments{
  \item{FUN}{model function with first argument a vector of parameter values of length ndim}
  \item{func.type}{type of value FUN returns.
    one of: posterior.density, logposterior.density,calc.loglik,calc.rmse,calc.weighted.rmse}
  \item{pars}{a list of variable ranges. Any names will be propagated to output}
  \item{FUN.pars}{A named list of extra arguments to FUN}
  \item{INIT}{A function f(pars,nseq,...) returning an nseq x ndim
    matrix of initial parameter values}
  \item{INIT.pars}{A named list of extra arguments to INIT}
  \item{control}{list of settings for dream. see below}
  \item{measurement}{Required parameter for func.types starting with calc. list with
    element data, a vector containing observed dependent variable values.
  }
}
\details{
  Elements of control are:
  \tabular{lll}{
    \emph{Element} \tab \emph{Default} \tab \emph{Description} \cr
    \code{ndim} \tab \tab Number of parameters. Calculated from
    parameters \code{pars} \cr
    \code{nseq} \tab \code{ndim} \tab Number of parallel chains to evolve \cr
    \code{DEpairs} \tab \code{(nseq-1)/2} \tab Number of pairs to evolve at each
    generation \cr
    \code{nCR} \tab 3 \tab Crossover values used to generate proposals (geometric
    series). scalar \cr
  \code{gamma} \tab 0 \tab Kurtosis parameter Bayesian Inference Scheme \cr
  \code{steps} \tab 10 \tab Number of steps in sem \cr
  \code{eps} \tab 5e-2 \tab Random error for ergodicity \cr
  \code{outlierTest} \tab 'IQR_test' \tab Test used to detect outlier
    chains.
    One off: IQR_test, Grubbs_test,Mahal_test \cr
  \code{pCR.Update} \tab \code{TRUE} \tab Whether to use adaptive tuning of crossover
    values \cr
  \code{boundHandling} \tab 'reflect' \tab Method used to handle parameter
    values outside of parameter bounds.
    One of: "reflect", "bound", "fold", "none","rand"
  \cr
  \code{thin.t} \tab NA \tab Thinning interval. NA if thinning is not
  desired \cr
  \code{REPORT} \tab 1000 \tab Approximate number of function evaluations
    between calculation and reporting of convergence diagnostics. 0 if
    no reporting or calculation is desired.
    Frequent calculation is likely to slow performance.
    \code{\link{gelman.diag}} can calculate convergence statistic after completion
  \cr
  \code{ndraw} \tab  1e5 \tab Maximum number of function evaluations. May
    terminate before convergence. \cr
  \code{maxtime} \tab Inf \tab Maximum duration of optimization in seconds. May
    terminate before convergence.\cr
  \code{Rthres} \tab 1.01 \tab Value of Gelman & Rubin's convergence
    diagnostic R value below which the sequences are considered to have
    converged, and execution is terminated. Vrugt suggests 1.2
  \cr
  }
  Execution terminates when Gelman-Rubin statistic < \code{control$Rthres}, or
  \code{control$ndraw} or \code{control$maxtime} are reached

  Two options for INIT are provided: latin hypercube sampling
  \code{\link{LHSInit}} and covariance-based sampling \code{\link{CovInit}}.
  
  Default settings may not be sufficient for model to run. For few
  parameters, nseq may need to be increased.
  
  There are methods for print, summary, plot, coef

  coef uses \code{\link{maxLikPars}}, a naive approach to extract the most likely parameter value.

}
\value{
  A list with elements:
  \item{X}{converged nseq points in parameter space. matrix nseq x ndim}
  \item{Sequences}{ \code{\link{mcmc.list}}. nseq mcmc elements of ndim variables}
  \item{Reduced.Seq}{ \code{\link{mcmc.list}}. nseq mcmc elements of ndim
    variables, if \code{control$thin.t!=NA}}
  \item{AR}{Acceptance rate for each draw. matrix n.elem x 2}
  \item{outlier}{Iterations at which chains were replaced. vector of variable length}
  \item{R.stat}{gelman.diag statistic for each variable at each step. matrix n.elem/steps x 1+ndim}
  \item{CR}{Probability of crossover at each step. matrix n.elem/steps x 1+length(pCR)}
}

\seealso{
  \code{\link{possibility.envelope}}
  
  Examples in demo folder:
  \itemize{
    \item{FME non linear model: }{Calibrating the non-linear model shown in the FME vignette}
    \item{example1: }{Fitting a banana shaped distribution - the first example
  in Vrugt matlab code}
    }
  }

\references{

  Accelerating Markov Chain Monte Carlo Simulation by Differential Evolution with Self-Adaptive Randomized Subspace Sampling
Jasper A. Vrugt, C.J.F. ter Braak, C.G.H. Diks, Bruce A. Robinson, James M. Hyman and Dave Higdon
 International Journal of Nonlinear Sciences and Numerical Simulation
 (in press 2009)

 Equifinality of formal (DREAM) and informal (GLUE)
Bayesian approaches in hydrologic modeling?
Jasper A. Vrugt, Cajo J. F. ter Braak, Hoshin V. Gupta, Bruce A. Robinson
Stoch Environ Res Risk Assess (2009) 23:1011\2261026
DOI 10.1007/s00477-008-0274-y
  }