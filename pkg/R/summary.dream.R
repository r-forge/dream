
summary.dream <- function(object,...){

  cat(sprintf("
Exit message:  %s
Num fun evals: %d
Time (secs):   %.1f
Final R.stats:
",
              object$EXITMSG,
              object$fun.evals,
              object$time
              ))

  R.stat.last <- tail(object$R.stat,1)
  for (i in 2:ncol(object$R.stat)){
    cat(sprintf("\t%s:\t%f\n",
                colnames(object$R.stat)[i],
                R.stat.last[i]
                ))
  } ##for

  cat("
CODA summary for last 50% of MCMC chains:
")
  print(summary(window(object$Sequences, start = end(object$Sequences)/2 + 1)))

  cat("
Acceptance Rate
")
  summary(object$AR[,2])
  
} ##summary.dream
