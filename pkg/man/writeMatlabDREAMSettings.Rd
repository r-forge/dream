\name{writeMatlabDREAMSettings}
\alias{writeMatlabDREAMSettings}
\title{
  Write settings for dream in matlab
}
\description{
  Write matlab .mat file with settings and .m file to run dream in Matlab
}
\usage{
writeMatlabDREAMSettings(dream.obj, ModelName, InitPopulation, matlab.dream.dir, in.mat.file = "in.mat", out.mat.file = "out.mat", run.m.file = "run_once.m")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dream.obj}{
    \code{\link{dream}} object
}
  \item{ModelName}{
    name of matlab .m file and function containing model to be run
}
  \item{InitPopulation}{
    Setting to generate InitPopulation in dream matlab version
}
  \item{matlab.dream.dir}{
    Directory in which to write files
}
  \item{in.mat.file}{
    name of .mat file to which to write dream settings
}
  \item{out.mat.file}{
    name of .mat file to which to write output of dream when run (used
  in \code{run.m.file}.
}
  \item{run.m.file}{
    Name of Matlab .m file to write which will run dream
}
}
\details{
  Used to create demo
}
\value{
  Creates files to allow running in Matlab dream version
}
