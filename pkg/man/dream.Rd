\name{dream}
\alias{dream}
\title{DiffeRential Evolution Adaptive Metropolis (DREAM)}
\usage{

dream(FUN, func.type, pars, FUN.pars = list(), INIT = LHSInit,
      INIT.pars = list(), control = list(), measurement)

}
\description{
  Efficient global MCMC even in high-dimensional spaces.
  From J.A. Vrugt, C.J.F. ter Braak et al.
}
\arguments{
  \item{FUN}{
    model function with first argument a vector of parameter values of
    length ndim.
  }
  \item{func.type}{
    type of value \code{FUN} returns. One of:
    \describe{
      \item{\code{"posterior.density"}}{
        \code{FUN} returns the likelihood.
      }
      \item{\code{"logposterior.density"}}{
        \code{FUN} returns the log-likelihood.
      }
      \item{\code{"calc.loglik"}}{
        \code{FUN} returns the residual vector to be compared with
        \code{measurement}.
        The exponential power density function of Box and Tiao (1973) is
        used, based on the parameter \code{control$gamma}, and
        \code{measurement$N} and \code{measurement$sigma}.
      }
      \item{\code{"calc.rmse"}}{
        \code{FUN} returns the residual vector to be compared with
        \code{measurement}. The likelihood is calculated as the sum of
        squared residuals to a power \code{-N}, adjusted by
        \code{control$gamma}. 
      }
      \item{\code{"calc.weighted.rmse"}}{
        similar to \code{"calc.rmse"} but weighted with
        \code{measurement$sigma}. 
      }
    }
  }
  \item{pars}{
    a list of variable ranges. Any names will be propagated to output.
  }
  \item{FUN.pars}{
    a list of any extra arguments to be passed to \code{FUN}.
  }
  \item{INIT}{
    A function f(pars,nseq,...) returning an \code{nseq x ndim}
    matrix of initial parameter values. 
  }
  \item{INIT.pars}{
    a list of any extra arguments to be passed to the \code{INIT} function.
  }
  \item{control}{
    list of settings for the DREAM algorithm. see below.
  }
  \item{measurement}{
    Required parameter for \code{func.type}s starting with \code{calc}.
    A list with element \code{data}, a vector containing observed
    dependent variable values. 
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
    One of: IQR_test, Grubbs_test, Mahal_test \cr
    \code{pCR.Update} \tab \code{TRUE} \tab Whether to use adaptive tuning of crossover
    values \cr
    \code{boundHandling} \tab 'reflect' \tab Method used to handle parameter
    values outside of parameter bounds.
    One of: "reflect", "bound", "fold", "none","rand"
    \cr
    
    \code{burnin.length} \tab 0.1 \tab If \code{<1} proportion of
    function evaluations to use as burn-in period. If \code{>1} number
    of function evaluations to use.
    In the burnin-period, adaptive tuning of crossover values is
    performed if \code{pCR.Update=TRUE} and the method specified by
    \code{outlierTest} is used to remove outliers. 
    If outliers are detected outside the burn-in period, dream emits a
    warning and returns to the burn-in period for another \code{burnin.length} iterations.
    If \code{=0} there is no burn-in period, and all outlier warnings
    are suppressed. \cr
    
    \code{thin.t} \tab NA \tab MCMC chain thinning interval. NA if thinning is not
    desired \cr
    \code{REPORT} \tab 1000 \tab Approximate number of function evaluations
    between calculation and reporting of convergence diagnostics. 0 if
    no reporting or calculation is desired.
    Frequent calculation is likely to slow performance. The summary
    function can use \code{\link{gelman.diag}} to calculate convergence statistic after completion
    \cr
    \code{parallel} \tab \code{"none"} \tab Character vector
    of packages to
    use for parallelisation of function evaluations, in order of
    preference. Options are \code{"multicore", "snow", "foreach"}. For foreach, a backend must be registered from one of the doMC, doMPI,
    doSNOW packages. \cr
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

  If dream frequently warns that it is reentering the burn-period,
  \code{control$burnin.length} may need to be increased.
  
  There are S3 methods for print, summary,
  \code{\link[=plot.dream]{plot}},
  \code{\link[=coef.dream]{coef}},
  \code{\link[=predict.dream]{predict}},
  simulate

  coef uses \code{\link{maxLikPars}}, a naive approach to extracting the
  most likely parameter value, using the bandwith parameters used for
  CODA density plots.

}
\value{
  A list with elements:
  \item{X}{converged nseq points in parameter space. matrix \code{nseq x ndim}.}
  \item{Sequences}{ \code{\link{mcmc.list}}. \code{nseq} mcmc elements of \code{ndim} variables.}
  \item{Reduced.Seq}{ \code{\link{mcmc.list}}. \code{nseq} mcmc elements of \code{ndim}
    variables, if \code{control$thin.t!=NA}.}
  \item{hist.logp}{History of log(\var{p}) values. matrix \code{nseq x n.elem}.}
  \item{AR}{Acceptance rate for each draw. matrix \code{n.elem x 2}.}
  \item{outlier}{Iterations at which chains were replaced. vector of variable length}
  \item{R.stat}{gelman.diag statistic for each variable at each step. matrix \code{n.elem/steps x 1+ndim}.}
  \item{CR}{Probability of crossover at each step. matrix \code{n.elem/steps x
    1+length(pCR)}.}
  \item{in.burnin}{Boolean showing whether dream was in the burn-in
    period at termination (see description above).}
  \item{fun.evals}{Total number of function evaluations.}
  \item{time}{running time (wall time) in seconds.}
  \item{EXITMSG}{a message about the stopping condition.}
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
  
  Vrugt, J. A., ter Braak, C. J. F., Diks, C. G. H., Robinson, B. A.,
  Hyman, J. M., Higdon, D., 2009. Accelerating Markov chain Monte Carlo
  simulation by differential evolution with self-adaptive randomized
  subspace sampling. \emph{International Journal of Nonlinear Sciences
  and Numerical Simulation} 10 (3), 273-290. 
 
  Vrugt, J. A., Braak, C. J. F., Gupta, H. V., Robinson, B. A.,
  2009. Equifinality of formal (DREAM) and informal (GLUE) Bayesian
  approaches in hydrologic modeling?
  \emph{Stochastic Environmental Research and Risk Assessment} 23 (7), 1011--1026. 
  \url{http://dx.doi.org/10.1007/s00477-008-0274-y}
  
}
