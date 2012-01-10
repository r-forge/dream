## Example of parallelisation
## Uses FME data (see other example)

obs.all <- data.frame(x=c(   28,  55,   83,  110,  138,  225,  375),   # mg COD/l
                  y=c(0.053,0.06,0.112,0.105,0.099,0.122,0.125))   # 1/hour


unloadNamespace("dream")
library(dream)

control <- list(
                nseq=4
                )
pars <- list(p1=c(0,1),p2=c(0,100))


##########################
## Demonstration of parallelisation packages
##  Comparison of run times for known fixed termination problem

Model.y <- function(p,x) {
  p[1]*x/(x+p[2])

}

system.time(sapply(1:5e4,function(x) Model.y(c(0,0),obs.all$x)))[["elapsed"]]/5e4

## Optional parallelisation of function evaluations
## Uses foreach, which allows Rmpi, SNOW or multicore to be registered,
##   or just normal sequential evaluation
## Example here is for SNOW with a socket cluster, which doesn't require any other software or setup
## Parallelisation incurs an overhead and is not worthwhile for fast functions
if (require(doSNOW)){
  cl <- makeCluster(2, type = "SOCK")
  registerDoSNOW(cl)
}

for (p in c("none","snow","foreach")){

  print(p)
  control$parallel <- p

  set.seed(456)
  dd <- dreamCalibrate(
                       FUN=Model.y,
                       pars = pars,
                       obs=obs.all$y,
                       FUN.pars=list(
                         x=obs.all$x
                         ),
                       control = control
                       )


  print("Coefficients:")
  print(coef(dd))
  print(sprintf("Elapsed time %f seconds",dd$time))
}
if (require(doSNOW))  stopCluster(cl)

## Seconds taken:
## none: 6.234, 6.219
## snow: 12.234, 12.125
## foreach: 107.999, 107.893

## For function with 1e-3 sleep
## [1] "none"
##         p1         p2
##  0.1511961 50.7765344
## [1] 109.483
## [1] "snow"
##         p1         p2
##  0.1511961 50.7765344
## [1] 54.858
## [1] "foreach"
##         p1         p2
##  0.1511961 50.7765344
## [1] 109.593


##################
## Test of number of function evaluations in fixed time for a time-expensive function

Model.y <- function(p,x) {
  Sys.sleep(0.1)
  p[1]*x/(x+p[2])
}

system.time(sapply(1:5e1,function(x) Model.y(c(0,0),obs.all$x)))[["elapsed"]]/5e1

if (require(doSNOW)){
  cl <- makeCluster(2, type = "SOCK")
  registerDoSNOW(cl)
}

for (p in c("none","snow","foreach")){
  print(p)
  set.seed(456)
  control$parallel <- p
  control$maxtime <- 20

  dd <- dreamCalibrate(
                       FUN=Model.y,
                       pars = pars,
                       obs=obs.all$y,
                       FUN.pars=list(
                         x=obs.all$x
                         ),
                       control = control
                       )
  print(sprintf("Number of function evaluations: %f",dd$fun.evals))
}

if (require(doSNOW))  stopCluster(cl)

## Number of function evaluations 10/01/2012
## Note foreach has significant overhead
## [1] "none"
## [1] "Number of function evaluations: 200.000000"
## [1] "snow"
## [1] "Number of function evaluations: 400.000000"
## [1] "foreach"
## [1] "Number of function evaluations: 280.000000"
