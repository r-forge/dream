\name{likelihood_functions}
\alias{calc.weighted.rmse}
\alias{calc.rmse}
\alias{calc.loglik}
\title{
Likelihood functions
}
\description{
Calculates likelihood or log likelihood given predicted and observed values
}
\usage{
calc.rmse(predicted, observed, control = list(gamma = 0))
calc.weighted.rmse(predicted, observed, control = list(gamma = 0))
calc.loglik(predicted, observed, control = list(gamma = 0))
}
\arguments{
  \item{predicted}{
    output of model
}
  \item{observed}{
    observed values expected, to compare with output of model
}
  \item{control}{
    list of options
}
}
\details{
  \code{calc.rmse} is the default for \code{\link{dreamCalibrate}} but it is recommended to use an objective
  function better suited to the particular problem

  \itemize{
    \item Sigma-weighted RMSE log likelihood function
    \item Calculate log RMSE using formulas used by Vrugt
    }
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  log-likelihood
}
