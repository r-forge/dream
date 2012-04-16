\name{snow.chains}
\alias{snowchains}
\alias{snow.chains}
\title{DREAM with external batch files in parallel folders}
\description{
  N.B. This functionality is only provided on request for expert
  users of external model software. It is assumed that users are capable
  of setting up batch files, executing them from R, and reading in the
  results. The authors can only provide limited support.

  With the setting \code{control$parallel="snow.chains"}, dream
  interprets \code{FUN} to be of the form \code{logp=f(chain.id, list of
  pars)}. By passing the \code{chain.id} to the function, the user can
  run external programs in parallel (using SNOW) within separate
  folders.

  In addition:
  \itemize{
   \item \code{func.type} must be "logposterior.density",
   i.e. \code{FUN} must return the log-likelihood
   \item \code{measurement} should be left at its default value of
   \code{NULL} or dream will stop with an error.
   \item \code{FUN.pars} should be left at its
   default value of \code{list()} or dream will stop with an error.
   \item dreamCalibrate cannot be used
   }
 }
 \seealso{\code{\link{dream}} for documentation of other arguments, and
  \code{\link{dreamCalibrate}} for the recommended way to calibrate a
  function using dream.
  demo("parallelisation_chain_id")
}
