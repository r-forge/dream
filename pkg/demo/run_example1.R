
source("example1.R")

set.seed(11)

dd <- dream(
            FUN=Banshp, func.type="logposterior.density",
            pars = pars,
            FUN.pars=list(
              imat=solve(cmat),
              bpar=0.1),
            INIT = CovInit,
            INIT.pars=list(
              muX=muX,
              qcov=qcov,
              bound.handling=control$boundHandling
              ),
            control = control
            )

plot(dd)
