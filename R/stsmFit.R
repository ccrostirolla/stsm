
stsmFit <- function(x, stsm.method = c("maxlik.fd.scoring", "maxlik.td.scoring", 
  "maxlik.fd.optim", "maxlik.td.optim"), ...)
{
  method <- match.arg(stsm.method)
  
  do.call(stsm.method, args = c(list(m = x), list(...)))
}
