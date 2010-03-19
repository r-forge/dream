
## Example from FME vignette, p21, section 7. Soetaert & Laine
## Nonlinear model parameter estimation
## TODO: document more clearly, better outputs
## 

## http://r-forge.r-project.org/plugins/scmsvn/viewcvs.php/pkg/FME/inst/doc/FMEmcmc.Rnw?rev=96&root=fme&view=markup

Obs <- data.frame(x=c(   28,  55,   83,  110,  138,  225,  375),   # mg COD/l
                  y=c(0.053,0.06,0.112,0.105,0.099,0.122,0.125))   # 1/hour
Obs2<- data.frame(x=c(   20,  55,   83,  110,  138,  240,  325),   # mg COD/l
                   y=c(0.05,0.07,0.09,0.10,0.11,0.122,0.125))   # 1/hour
obs.all <- rbind(Obs,Obs2)

##########################
## DREAM results

unloadNamespace("dream")
library(dream)

Model.y <- function(p,x) p[1]*x/(x+p[2])

control <- list(
                nseq=4
                )

pars <- list(p1=c(0,1),p2=c(0,100))

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
