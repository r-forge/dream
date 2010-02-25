print.dream <- function(dream.obj){
  cat("
Call:
")
  
  print(dream.obj$call)
  
  cat("
Control:
")
  for (i in names(dream.obj$control)){
    v <- dream.obj$control[[i]]
    cat(sprintf("%15s: ",i))
    cat(v,"\n")
  } 

  cat("\nExit condition:",dream.obj$EXITMSG,"\n")
}
