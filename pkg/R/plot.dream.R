##' Plot characteristics of dream object
##' @param dream object
##' @param interactive - stop for each plot

##' Uses second half of sequences

plot.dream <- function(x,interactive=TRUE,...){
  devAskNewPage(interactive)
  
  ss <- window(x$Sequences, start = end(x$Sequences)/2 + 1)

  ## Trace and parameter density
  
  plot(ss)

  print(xyplot(ss))
  print(densityplot(ss))

  ## Acceptance rate
  plot(table(x$AR[,2]),main="Distribution of % acceptance rate")
  
  ##Convergence
  
  try(gelman.plot(ss))

  plot(x$R.stat[,1],x$R.stat[,2],type="l",ylim=c(0,2))
  for (i in 2:x$control$ndim) lines(x$R.stat[,1],x$R.stat[,i+1],ylim=c(0,2))
  title(main="Evolution of R.stat",sub="Equivalent to gelman.plot")

  ## Multi-variate density for first chain
  
  print(splom(as.data.frame(x$Sequences[[1]]),
      upper.panel = panel.smoothScatter, nrpoints = 0,
      lower.panel = function(x, y, ...) {
          panel.grid(-1, -1)
          panel.loess(x, y, span = 1/3, lwd = 1)
          panel.loess(x, y, span = 2/3, lwd = 2)
          grid::grid.text(paste("cor =", round(cor(x, y),2)),
                          y = 0.1)
      }))
  
}##plot.dream
