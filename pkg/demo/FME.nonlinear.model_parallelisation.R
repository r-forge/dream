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

  set.seed(456)
  control$parallel <- p
  
  dd <- dream(
              FUN=Model.y, func.type="calc.rmse",
              pars = pars,
              FUN.pars=list(
                x=obs.all$x
                ),
              INIT = LHSInit,
              measurement=list(data=obs.all$y),
              control = control
              )

  print(dd$control$parallel)
  print(coef(dd))
  print(dd$time)
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
  Sys.sleep(1e-3)
  p[1]*x/(x+p[2])

}

system.time(sapply(1:5e1,function(x) Model.y(c(0,0),obs.all$x)))[["elapsed"]]/5e1


if (require(doSNOW)){
  cl <- makeCluster(2, type = "SOCK")
  registerDoSNOW(cl)
}

for (p in c("none","snow","foreach")){
  set.seed(456)
  control$parallel <- p
  control$maxtime <- 20
  
  dd <- dream(
              FUN=Model.y, func.type="calc.rmse",
              pars = pars,
              FUN.pars=list(
                x=obs.all$x
                ),
              INIT = LHSInit,
              measurement=list(data=obs.all$y),
              control = control
              )
  print(dd$control$parallel)
  print(dd$fun.evals)
}

if (require(doSNOW))  stopCluster(cl)

## Number of function evaluations
## TODO: appears to be an error in foreach evaluation
## [1] "none"
## [1] 1280
## [1] "snow"
## [1] 2560
## [1] "foreach"
## [1] 1280
