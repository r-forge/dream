

fitted.dream <-
    function(object,
             start = 1+(end(object$Sequences)-1)*(1-fraction),
             fraction = 0.5,
             ...)
{
    window(object$Sequences, start = start, ...)
}
