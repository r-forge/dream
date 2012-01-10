\name{compareToMatlab}
\alias{compareToMatlab}
\title{
Compare dream and matlab results
}
\description{
  Compare results from dream run to results from matlab
}
\usage{
compareToMatlab(mat.file, dream.obj)
}
\arguments{
  \item{mat.file}{
    Path to matlab .mat file
}
  \item{dream.obj}{
    object returned by dream
}
}
\details{
  Used in demos
}
\value{
  Prints diagnostics. \cr
  Did matlab and R use same inputs? \cr
  Do outputs have same dimensions? \cr
  ks.test that samples are from same distribution for each variable \cr
  returns result of \code{\link{plotMCMCQQ}}
}
