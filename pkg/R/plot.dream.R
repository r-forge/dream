##' Plot characteristics of dream object
##' @param dream object
##' @param interactive - stop for each plot

##' Uses second half of sequences

plot.dream <- function(dream.obj,interactive=TRUE){
  devAskNewPage(interactive)
  
  ss <- window(dream.obj$Sequences, start = end(dd$Sequences)/2 + 1)

  ## Trace and parameter density
  
  plot(ss)

  xyplot(ss)
  densityplot(ss)

  ## Acceptance rate
  plot(table(dd$AR[,2]),main="Distribution of % acceptance rate")
  
  ##Convergence
  
  try(gelman.plot(ss))

  plot(dd$R.stat[,1],dd$R.stat[,2],type="l",ylim=c(0,2))
  for (i in 2:dd$control$ndim) lines(dd$R.stat[,1],dd$R.stat[,i+1],ylim=c(0,2))
  title(main="Evolution of R.stat",sub="Equivalent to gelman.plot")

  ## Multi-variate density for first chain
  
  splom(as.data.frame(dd$Sequences[[1]]),
      upper.panel = panel.smoothScatter, nrpoints = 0,
      lower.panel = function(x, y, ...) {
          panel.grid(-1, -1)
          panel.loess(x, y, span = 1/3, lwd = 1)
          panel.loess(x, y, span = 2/3, lwd = 2)
          grid::grid.text(paste("cor =", round(cor(x, y),2)),
                          y = 0.1)
      })
  
}##plot.dream
