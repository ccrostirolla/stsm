
##FIXME TODO for other methods (time domain, concentrated likelihood, optim, EM)

##FIXME if xreg is not null (in ellipsis, list(...)) and derivatives is "analytical"
#change to numerical (at least for maxlik.fd/td.optim functions)

#stsmFit <- function(x, xreg = NULL, method = c("maxlik.fd.scoring"), args = NULL)
#stsmFit <- function(x, xreg = NULL, method = c("maxlik.fd.scoring"), ...)
stsmFit <- function(x, stsm.method = c("maxlik.fd.scoring", "maxlik.td.scoring", 
  "maxlik.fd.optim", "maxlik.td.optim"), ...)
{
  method <- match.arg(stsm.method)
  
  #do.call(method, args = c(list(m = x, xreg = xreg), args))
  do.call(stsm.method, args = c(list(m = x), list(...)))
}

# do.call("stsmFit", args = list(x = x, args = list()))
# fcn <- function(x, xreg = NULL, method = c("maxlik.fd.scoring"), ...)
# {
#   #match.call(expand.dots = TRUE)
#   res <- list(...)
#   res
# }
