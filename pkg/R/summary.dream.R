
summary.dream <- function(dream.obj){

  cat(sprintf("
Exit message:  %s
Num fun evals: %d
Time (secs):   %.1f
Final R.stats:
",
              dream.obj$EXITMSG,
              dream.obj$fun.evals,
              dream.obj$time
              ))

  R.stat.last <- tail(dream.obj$R.stat,1)
  for (i in 2:ncol(dream.obj$R.stat)){
    cat(sprintf("\t%s:\t%f\n",
                colnames(dream.obj$R.stat)[i],
                R.stat.last[i]
                ))
  } ##for

  cat("
CODA summary for last 50% of MCMC chains:
")
  print(summary(window(dream.obj$Sequences, start = end(dd$Sequences)/2 + 1)))

  cat("
Acceptance Rate
")
  summary(dream.obj$AR[,2])
  
} ##summary.dream
