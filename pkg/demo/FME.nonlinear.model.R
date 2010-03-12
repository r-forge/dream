
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


###########################
### FME results

Model <- function(p,x) return(data.frame(x=x,y=p[1]*x/(x+p[2])))
##Model(c(0.1,1),obs.all$x)


library(FME)
Residuals  <- function(p) {
   cost<-modCost(model=Model(p,Obs$x),obs=Obs,x="x")
   modCost(model=Model(p,Obs2$x),obs=Obs2,cost=cost,x="x")
}

P      <- modFit(f=Residuals,p=c(0.1,1))
print(P$par)

plotFME <- function(){
plot(Obs,xlab="mg COD/l",ylab="1/hour", pch=16, cex=1.5,
     xlim=c(25,400),ylim=c(0,0.15))
points(Obs2,pch=18,cex=1.5, col="red")
lines(Model(p=P$par,x=0:375))
}


##########################
## DREAM results


unloadNamespace("dream")
library(dream)


Model.y <- function(p,x) as.ts(p[1]*x/(x+p[2]))

set.seed(456)

control <- list(
                nseq=4,
                use.multicore=FALSE
                ##                REPORT=0
                ##                ndraw=1000
                ##                Rthres=1+1e-3
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

dd

summary(dd)

coef(dd)

plot(dd)

plotFME()
lines(Model(p=coef(dd),x=0:375),col="green")


########################
## Calculate bounds around output estimates

plotCIs <- function(x,cis,...){
  em <- strwidth("M")/2
  segments(x0=x,y0=cis[,1],
           x1=x,y1=cis[,2],
           ...
           )
  segments(x0=x-em,y0=cis[,1],
           x1=x+em,y1=cis[,1],
           ...
           )
  segments(x0=x-em,y0=cis[,2],
           x1=x+em,y1=cis[,2],
           ...
           )
}##plotCIs

## Calibrate with Obs
dd <- dream(
            FUN=Model.y, func.type="calc.rmse",
            pars = pars,
            FUN.pars=list(
              x=Obs$x
              ),
            INIT = LHSInit,
            measurement=list(data=Obs$y),
            control = control
            )

##Obs1
plotFME()
lines(Obs$x,predict(dd),col="blue")
lines(Obs$x,predict(dd,method="mean"),col="red")
lines(Obs$x,predict(dd,method="median"),col="orange")
plotCIs(Obs$x,predict(dd,method="CI"),col="black")

##Obs2
plotFME()
lines(Obs2$x,predict(dd,list(x=Obs2$x)),col="blue")
lines(Obs2$x,predict(dd,list(x=Obs2$x),method="mean"),col="red")
lines(Obs2$x,predict(dd,list(x=Obs2$x),method="median"),col="orange")
plotCIs(Obs2$x,predict(dd,list(x=Obs2$x),method="CI"),col="red")

### Example with new sample
dd.sim <- simulate(dd)
predict(dd.sim)
plotFME()
lines(Obs2$x,predict(dd,list(x=Obs2$x)),col="blue")
lines(Obs2$x,predict(dd.sim,list(x=Obs2$x)),col="purple")


########################
## Legacy examples
plotFME()
lines(Model(p=coef(dd),x=0:375),col="green")

## Naive 95% bounds from residuals
resid <- Model.y(p=coef(dd),x=Obs$x)-Obs$y
##densityplot(resid)
qq <- quantile(resid,c(0.005,.995))
gg <- t(sapply(Model.y(p=coef(dd),x=Obs$x),function(v) v+qq))
plotCIs(Obs$x,gg,col="grey")

## Test on Obs2
cis.2 <- predict(dd,list(x=Obs2$x),out="CI")
plotCIs(Obs2$x,cis.2,col="red")

## TODO: add residual error, using method p6, vrugt. equifinality
