## Example of parallelisation
## Uses FME data (see other example)

obs.all <- data.frame(x=c(   28,  55,   83,  110,  138,  225,  375),   # mg COD/l
                  y=c(0.053,0.06,0.112,0.105,0.099,0.122,0.125))   # 1/hour


control <- list(
                nseq=4
                )
pars <- list(p1=c(0,1),p2=c(0,100))


##########################
## Demonstration of parallelisation packages
##  Comparison of run times for known fixed termination problem

lik <- function(Qobs,Qs,...){
  res<-Qobs-Qs
  sd.res<-sd(res)
  var.res<-var(res)
  logp<-sum(sapply(res,function(res_k)
                   log(1/sqrt(2*pi)*sd.res)*exp(-res_k^2/(2*var.res))))
}

Model.split <- function(id,pars){

  ##Make sure that you have as many batch scripts as you have parameter chains
  ## placeholder for this example

  ##Select out this parameter set
  this.par.set=pars[id,]

  ##Do something to set the parameters in SWAT for this instance
  ## Because they're running on other nodes, they don't have access to data or functions in R
  x <- c(   28,  55,   83,  110,  138,  225,  375)
  Model.y <- function(p,x) {
    p[1]*x/(x+p[2])
  }

  ##Run the model
  ans <- Model.y(this.par.set,x)

  ##Do something to collect the answer, in appropriate form for lik.fun in dreamCalibrate
  return(ans)
}

Model.fit <- function(id,pars){
  ##Make sure that you have as many batch scripts as you have parameter chains
  ## placeholder for this example

  ##Select out this parameter set
  this.par.set=pars[id,]

  ##Do something to set the parameters in SWAT for this instance
  ## Because they're running on other nodes, they don't have access to data or functions in R
  x <- c(   28,  55,   83,  110,  138,  225,  375)
  Model.y <- function(p,x) {
    p[1]*x/(x+p[2])
  }

  ##Run the model
  ans <- Model.y(this.par.set,x)

  ##Do something to collect the answer as logp
  Qs <- ans
  Qobs <- c(0.053,0.06,0.112,0.105,0.099,0.122,0.125)

  res<-Qobs-Qs
  sd.res<-sd(res)
  var.res<-var(res)
  logp<-sum(sapply(res,function(res_k)
                   log(1/sqrt(2*pi)*sd.res)*exp(-res_k^2/(2*var.res))))
  return(logp)
}

## Example here is for SNOW with a socket cluster, which doesn't require any other software or setup
## Parallelisation incurs an overhead and is not worthwhile for fast functions
library(doSNOW)
cl <- makeCluster(2, type = "SOCK")
registerDoSNOW(cl)

## Using dream directly, with Model.fit, which returns a likelihood
set.seed(456)
control$parallel <- "snow.chains"
dd <- dream(FUN = Model.fit,
            pars = pars,
            func.type = "logposterior.density",
            control=control
            )
summary(dd)

## TODO: dreamCalibrate cannot be used because of the way it calls FUN
## ## Using dreamCalibrate, which allows the likelihood function to be changed separately from the model function
## set.seed(456)
## dd <- dreamCalibrate(
##                      FUN=Model.split,
##                      pars = pars,
##                      obs=obs.all$y,
##                      lik.fun=lik,
##                      control = control
##                      )
## summary(dd)

stopCluster(cl)
