
## Example from FME vignette, p21, section 7. Soetaert & Laine
## Nonlinear model parameter estimation
## TODO: document more clearly, better outputs
## 


unloadNamespace("dream")
library(dream)

## http://r-forge.r-project.org/plugins/scmsvn/viewcvs.php/pkg/FME/inst/doc/FMEmcmc.Rnw?rev=96&root=fme&view=markup

Obs <- data.frame(x=c(   28,  55,   83,  110,  138,  225,  375),   # mg COD/l
                  y=c(0.053,0.06,0.112,0.105,0.099,0.122,0.125))   # 1/hour
Obs2<- data.frame(x=c(   20,  55,   83,  110,  138,  240,  325),   # mg COD/l
                   y=c(0.05,0.07,0.09,0.10,0.11,0.122,0.125))   # 1/hour
obs.all <- rbind(Obs,Obs2)

Model <- function(p,x) return(data.frame(x=x,y=p[1]*x/(x+p[2])))
##Model(c(0.1,1),obs.all$x)

Model.y <- function(p,x) p[1]*x/(x+p[2])


control <- list(
                nseq=5,
                gamma=0,
                nCR=3,
                ndraw=1e5,
                steps=10,
                eps=5e-2,
                outlierTest='IQR_test',
                pCR.Update=TRUE,
                thin=TRUE,
                thin.t=10,
                boundHandling="fold"
                )

pars <- list(p1=c(0,1),p2=c(0,100))

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

cat(sprintf("
Exit message:  %s
Num fun evals: %d
Time (secs):   %.1f
",
            dd$EXITMSG,
            dd$fun.evals,
            dd$time
            ))
tail(dd$R.stat,1)
maxLikPars(window(dd$Sequences, start = 0.75*end(dd$Sequences) + 1))

ss <- window(dd$Sequences, start = end(dd$Sequences)/2 + 1)
pars.maxp <- maxLikPars(ss)
print(pars.maxp)

plot(ss)
gelman.plot(ss)


### FME

library(FME)
Residuals  <- function(p) {
   cost<-modCost(model=Model(p,Obs$x),obs=Obs,x="x")
   modCost(model=Model(p,Obs2$x),obs=Obs2,cost=cost,x="x")
}

P      <- modFit(f=Residuals,p=c(0.1,1))
print(P$par)

## rbind(
##       pars.maxp-maxp.res,
##       P$par,
##       pars.maxp+maxp.res
##       )

plot(Obs,xlab="mg COD/l",ylab="1/hour", pch=16, cex=1.5)
points(Obs2,pch=18,cex=1.5, col="red")
lines(Model(p=P$par,x=0:375))
lines(Model(p=pars.maxp,x=0:375),col="green")
