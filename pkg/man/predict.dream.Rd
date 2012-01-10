\name{predict.dream_model}
\alias{predict.dream_model}
\title{Predict values using dream_model object}
\usage{
\S3method{predict}{dream_model}(object, newdata=NULL,method="sample.ml", level=0.99, \dots)
}
\description{
  Predict values using function calibrated by \code{\link{dreamCalibrate}}, optionally with new data,
  using various methods of summarising the posterior parameter and
  output distributions
}
\arguments{
  \item{object}{dream_model object, as returned by \code{\link{dreamCalibrate}}}
  \item{newdata}{new FUN.pars list. If NULL, use object's.}
  \item{method}{CI or a \code{method} of \code{\link{coef.dream}}}
  \item{level}{Requested two-sided level of confidence. For CI method.}
  \item{...}{Extra arguments for \code{\link{window.dream}}}
}
\value{the same as the return value of \code{FUN}, provided to
  \code{\link{dreamCalibrate}}  (either numeric, ts or list). For
  \code{method="CI"}, either a matrix with upper and lower bound or list
  of matrices.
  }