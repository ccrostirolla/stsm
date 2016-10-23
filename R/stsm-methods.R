
setGeneric("get.pars", function(x, rescale = FALSE, gradient = FALSE){
  standardGeneric("get.pars") })
setMethod("get.pars", "stsm", function(x, rescale = FALSE, gradient = FALSE)
{
  tmp <- transPars(x, gradient = gradient, hessian = FALSE)
  p <- tmp$pars

  if (rescale)
  {
    idvars <- grep("^var\\d{1,2}$", names(x@pars))

    if (!is.null(x@cpar)) # 'x@pars' are relative variances
    {
      # absolute variances
      p[idvars] <- p[idvars] * x@cpar
    } else {
      warning("'rescale = TRUE' had no effect since slot 'cpar' is NULL.")
    }
  }

  # gradient=TRUE is used by char2numeric (see note there)
  if (gradient)
    return(list(pars = p, gradient = tmp$gradient))

  p  
})

setGeneric("get.nopars", function(x, rescale = FALSE){ 
  standardGeneric("get.nopars") })
setMethod("get.nopars", "stsm", function(x, rescale = FALSE)
{
  # relative variances if 'cpar' is not null, 
  # otherwise absolute variances
  p <- x@nopars
 
  if (rescale)
  {
    idvars <- grep("^var\\d{1,2}$", names(x@pars))

    if (!is.null(x@cpar)) # 'x@pars' are relative variances
    {
      # absolute variances
      p[idvars] <- p[idvars] * x@cpar 
    } else {
      warning("'rescale = TRUE' had no effect since slot 'cpar' is NULL.")
    }
  }

  p
})

setGeneric("get.cpar", function(x, rescale = FALSE){
  standardGeneric("get.cpar") })
setMethod("get.cpar", "stsm", function(x, rescale = FALSE)
{
  p <- x@cpar

##FIXME change use in stsm.sgf() (rescale=FALSE)
  if (!rescale && !is.null(x@cpar))
  {
    p[] <- 1 # p / p
  }
  
  # if 'rescale' is TRUE and x@cpar is NULL
  # NULL (x@cpar) is returned
  # no warning is given, output NULL is revealing enough

  p
})

setGeneric("set.pars", function(x, v, check = TRUE, inplace = FALSE){
  standardGeneric("set.pars") })
setMethod("set.pars", "stsm", function(x, v, check = TRUE, inplace = FALSE)
{
  nms <- names(v)  
  if (is.null(nms)) {
    nullnms <- TRUE
    nms <- seq_along(x@pars)
  } else nullnms <- FALSE

  if (is.numeric(v))
  {
    if (inplace) {
      #http://tolstoy.newcastle.edu.au/R/help/04/02/0966.html
      eval(eval(substitute(expression(x@pars[nms] <<- v))))      
    } else x@pars[nms] <- v
  } else
    stop("Wrong class of argument 'v'.")

  if (check)
  {
    if (nullnms) {
      if (length(x@pars) > 1 && length(v) == 1)
        warning("number of items to replace is not a multiple of replacement length")
    } #else # run 'validObject' anyway to check further issues (lower, upper bounds)
    validObject(x)
  }

  ifelse(inplace, return(invisible()), return(x))
})

setGeneric("set.nopars", function(x, v, check = TRUE, inplace = FALSE){ 
  standardGeneric("set.nopars") })
setMethod("set.nopars", "stsm", function(x, v, check = TRUE, inplace = FALSE)
{
  nms <- names(v)  
  if (is.null(nms)) {
    nullnms <- TRUE
    nms <- seq_along(x@nopars)
  } else nullnms <- FALSE

  if (is.numeric(v))
  {
    if (inplace) {
      #http://tolstoy.newcastle.edu.au/R/help/04/02/0966.html
      eval(eval(substitute(expression(x@nopars[nms] <<- v))))      
    } else x@nopars[nms] <- v
  } else
    stop("Wrong class of argument 'v'.")

  if (check)
  {
    if (nullnms) {
      if (length(x@pars) > 1 && length(v) == 1)
        warning("number of items to replace is not a multiple of replacement length")
    } #else # run 'validObject' anyway to check further issues (lower, upper bounds)
    validObject(x)
  }

  ifelse(inplace, return(invisible()), return(x))
})

setGeneric("set.cpar", function(x, value, check = TRUE, inplace = FALSE){
  standardGeneric("set.cpar") })
setMethod("set.cpar", "stsm", function(x, value, check = TRUE, inplace = FALSE)
{
  if (is.null(x@cpar))
  {
    return(warning("Nothing done. ", 
      "The model is not defined in terms of relative variances ",
      "(slot 'cpar' is null)."))
  }

  if (is.numeric(value))
  {    
    if (inplace) {
      #http://tolstoy.newcastle.edu.au/R/help/04/02/0966.html
      eval(eval(substitute(expression(x@cpar[] <<- value))))
    } else x@cpar[] <- value
  } else
    stop("Wrong class of argument 'value'.")

  if (check)
    validObject(x)

  ifelse(inplace, return(invisible()), return(x))
})

setGeneric("set.sgfc", function(x, inplace = FALSE){ 
  standardGeneric("set.sgfc") })
setMethod("set.sgfc", "stsm", function(x, inplace = FALSE)
{
  if (!is.null(x@sgfc))
  {
    return(warning("Nothing done, x@sgfc was not empty.", 
      "Set it to null and retry to overwrite."))
  }

  if (x@model %in%  c("local-level", "local-trend", "BSM", "llm+seas"))
  {
    tmp <- stsm.sgf(x, FALSE, FALSE, FALSE)$constants

    if (inplace) {
      #http://tolstoy.newcastle.edu.au/R/help/04/02/0966.html
      eval(eval(substitute(expression(x@sgfc <<- tmp))))
      return(invisible())
    } else {
      x@sgfc <- tmp
      return(x)
    }

  } else
    stop("'set.sgfc' is not implemented for model ", sQuote(x@model), ".")
})

setGeneric("set.xreg", function(x, xreg, coefs = NULL){
  standardGeneric("set.xreg") })
setMethod("set.xreg", "stsm", function(x, xreg, coefs = NULL)
{
  ss.xreg <- x@ss$xreg
  isnotnull.xreg <- !is.null(xreg)

  if (is.null(dim(xreg)))
  {
    if (isnotnull.xreg)
    {
      xreg <- cbind(xreg = xreg)
      x@ss$xreg <- xregnms <- "xreg"
    } else {
      x@ss$xreg <- NULL
    }
  } else {
  if (is.null(xregnms <- colnames(xreg)))
    xregnms <- paste("xreg", seq.int(ncol(xreg)), sep = "")
    colnames(xreg) <- xregnms
    x@ss$xreg <- xregnms
  }

  x@xreg <- xreg

  if (!is.null(ss.xreg))
  {
    id <- match(ss.xreg, names(x@pars))
    x@pars <- x@pars[-id]
    id <- match(ss.xreg, names(x@lower))
    x@lower <- x@lower[-id]
    x@upper <- x@upper[-id]
  }

  # coefficients of regressors are attached to slot "pars",
  # otherwise the regressors would be an offset and a model 
  # without regressors could be created for y - xreg %*% xregcoefs

  if (isnotnull.xreg)
  {
    xregcoefs <- rep(1, length(xregnms))
    names(xregcoefs) <- xregnms
    if (!is.null(coefs))
    {
      #if (anyNA(match(names(coefs), xregnms))))
      if (!all(names(coefs) %in% xregnms))
        stop("some of the names in ", dQuote("coefs"), 
          " do not match the column names in ", dQuote("xreg"), ".")
      xregcoefs[names(coefs)] <- coefs
    }

    x@pars <- c(x@pars, xregcoefs)
    xregcoefs[] <- Inf
    x@upper <- c(x@upper, xregcoefs)
    x@lower <- c(x@lower, -xregcoefs)
  }

  x
})

valid.stsmObject <- function(object)
{
  allnames <- names(c(object@pars, object@nopars, object@cpar))
  stopifnot(!any(duplicated(allnames)))
  ss <- object@ss

  uss <- unlist(ss)
  ref <- which(lapply(uss, FUN = class) == "call")
  if (length(ref) > 0)
  {
    uss <- c(uss, all.vars(as.expression(uss[ref])))
    uss[ref] <- NULL
  }
  if (!all(allnames %in% uss))
    stop("Wrong definition of 'pars', 'nopars'. (e1)")
  if (!all(uss %in% c("0", "1", "-1", allnames)))
    stop("Wrong definition of 'pars', 'nopars'. (e2)")

  if (!is.null(object@lower))
  {
    if (any(!names(object@lower) %in% names(object@pars)))
      stop("Parameters in 'lower' are not in 'pars'.")

    if (any(!names(object@pars) %in% names(object@lower)))
      stop("Parameters in 'pars' are not in 'lower'.")
  }
  if (!is.null(object@upper))
  {
    if (any(!names(object@upper) %in% names(object@pars)))
      stop("Parameters in 'upper' are not in 'pars'.")
    if (any(!names(object@pars) %in% names(object@upper)))
      stop("Parameters in 'pars' are not in 'upper'.")
  }

  if (!is.null(object@cpar))
  {
    if (any(names(object@cpar) %in% names(object@pars)))
      stop(paste("the parameter to be concentrated out of the", 
        "likelihood function is also defined in 'pars'."))
    if (any(names(object@cpar) %in% names(object@nopars)))
      stop(paste("the parameter to be concentrated out of the",
        "likelihood function is also defined in 'pars'."))
  }

  if (!is.null(object@xreg))
  {
    if (nrow(object@xreg) != length(object@y))
      stop("lengths of ", sQuote("y"), " and ", sQuote("xreg"), " do not match.")
  }

  if (any(is.finite(c(object@lower, object@upper))))
    check.bounds(object)

  TRUE
}

setValidity("stsm", valid.stsmObject)

setGeneric("check.bounds", function(x){ standardGeneric("check.bounds") })
setMethod("check.bounds", "stsm", function(x)
{
  ref <- cbind(x@pars, x@lower, x@upper)
  ref <- apply(ref, MARGIN = 1, FUN = function(x)
    findInterval(x[1], c(x[2], x[3]), rightmost.closed = TRUE))
  if (!all(ref) == 1) {
    stop("some parameter values are outside lower-upper bounds.")
  } else 
  invisible(TRUE)
})

setGeneric("char2numeric", function(x, P0cov = FALSE, rescale = FALSE){
  standardGeneric("char2numeric") })
setMethod("char2numeric", "stsm", function(x, P0cov = FALSE, rescale = FALSE)
{
  ss.fill <- function(M, allpars, allnames)
  {
    if (is.numeric(M)) {
      m <- M
    } else if (is.character(M))
    {
      dimM <- dim(M)
      if (is.null(dimM))
        dimM <- c(1, 1)
      m <- array(0, dim = dimM)
      ref <- which(M != "0" & M != "1")
      m[ref] <- allpars[match(M[ref], allnames)]
      m[which(M == "1")] <- 1
    } else
      stop("Wrong class in x@ss.")
    m
  }

  # gradient=TRUE is helpful when "char2numeric" is used to 
  # create "ss" to be passed is called in "KF.deriv";
  # returning here the gradient avoids a second call to this 
  # function in order to get the gradient
  
  getgrad <- !is.null(x@transPars)
  pars <- get.pars(x, rescale, getgrad)
  if (getgrad)
  {
    tpGradient <- pars$gradient
    pars <- pars$pars    
  } #else # this case is arranged at the end of this function
    #tpGradient <- NULL

##FIXME check (changed after renamed to 'rescale')
  allpars <- c(pars, get.cpar(x, rescale), get.nopars(x, rescale))
  allnames <- names(allpars)
  ss <- x@ss

  if (x@model == "local-level")
  {
    ss$Vid <- NULL
    ss$Qid <- NULL

    ss$H <- allpars["var1"]
    ss$V <- ss$Q <- rbind(allpars["var2"])
    ss$a0 <- allpars["a01"]
    ss$P0 <- rbind(allpars["P01"])
  } else 
  if (x@model == "local-trend")
  {
    ss$Vid <- c(0, 3)
    ss$Qid <- c(0, 3)

    ss$H <- allpars["var1"]
    ss$V <- ss$Q <- diag(c(allpars["var2"], allpars["var3"]))
    ss$a0 <- c(allpars["a01"], allpars["a02"])
    ss$P0 <- diag(c(allpars["P01"], allpars["P02"]))
  } else 
  if (x@model == "llm+seas" && frequency(x@y) == 4)
  {
    ss$Vid <- c(0, 3)
    ss$Qid <- c(0, 5)

    ss$H <- allpars["var1"]
    ss$V <- diag(c(allpars["var2"], allpars["var3"]))
    ss$Q <- rbind(cbind(ss$V,0,0),0,0)
    ss$a0 <- c(allpars["a01"], allpars["a02"], 
      allpars["a03"], allpars["a04"])
    ss$P0 <- diag(c(allpars["P01"], allpars["P02"], 
      allpars["P03"], allpars["P04"]))
  } else 
  if (x@model == "BSM" && frequency(x@y) == 4)
  {
    ss$Vid <- NULL
    ss$Qid <- NULL

    ss$H <- allpars["var1"]
    ss$V <- diag(c(allpars["var2"], allpars["var3"], allpars["var4"]))
    ss$Q <- rbind(cbind(ss$V,0,0),0,0)
    ss$a0 <- c(allpars["a01"], allpars["a02"], 
      allpars["a03"], allpars["a04"], allpars["a05"])
    ss$P0 <- diag(c(allpars["P01"], allpars["P02"], 
      allpars["P03"], allpars["P04"], allpars["P05"]))

  } else 
  if (x@model == "BSM" && frequency(x@y) == 12)
  {
    ss$Vid <- NULL
    ss$Qid <- NULL

    ss$H <- allpars["var1"]
    ss$V <- diag(c(allpars["var2"], allpars["var3"], allpars["var4"]))
    ss$Q <- rbind(cbind(ss$V,0,0,0,0,0,0,0,0,0,0),0,0,0,0,0,0,0,0,0,0)
    ss$a0 <- c(allpars["a01"], allpars["a02"], 
      allpars["a03"], allpars["a04"], allpars["a05"],
      allpars["a06"], allpars["a07"], allpars["a08"],
      allpars["a09"], allpars["a010"], allpars["a011"],
      allpars["a012"], allpars["a013"])
    ss$P0 <- diag(c(allpars["P01"], allpars["P02"], 
      allpars["P03"], allpars["P04"], allpars["P05"],
      allpars["P06"], allpars["P07"], allpars["P08"],
      allpars["P09"], allpars["P010"], allpars["P011"],
      allpars["P012"], allpars["P013"]))

  } else
  if (x@model == "trend+ar2")
  {
    ss$Vid <- NULL
    ss$Qid <- NULL

    p <- 2
    ss$T <- rbind(c(1, 1, rep(0, p)), c(0, 1, rep(0, p)),
      c(0, 0, allpars[c("phi1", "phi2")]), 
      cbind(0, 0, diag(p - 1), 0))
    ss$H <- allpars["var1"]
    ss$V <- diag(c(allpars["var2"], allpars["var3"], allpars["var4"]))
    ss$Q <- rbind(cbind(ss$V, 0), 0)
    ss$a0 <- c(allpars["a01"], allpars["a02"], 
      allpars["a03"], allpars["a04"])
    ss$P0 <- diag(c(allpars["P01"], allpars["P02"], 
      allpars["P03"], allpars["P04"]))

#NOTE allpars["P0i"] is ignored
#P0 <- diag(4)
#P0[1,1] <- 100000000
#tmp <- solve(diag(4) - kronecker(ss$T[c(3,4),c(3,4)], ss$T[c(3,4),c(3,4)])) %*% as.vector(ss$Q[c(3,4),c(3,4)])
#P0[c(3,4),c(3,4)] <- matrix(tmp, nrow=2, ncol=2)
#NOTE non-zero values outside the diagonal ??? maybe, because cov23 != 0
#ss$P0 <- P0

  } else 
  if (x@model == "level+AR2")
  {
    ss$Vid <- NULL
    ss$Qid <- NULL

    p <- 2
    ss$T <- rbind(c(1, 0, 0), c(0, allpars[c("phi1", "phi2")]), c(0,1,0))      
    ss$H <- allpars["var1"]
    ss$V <- diag(c(allpars["var2"], allpars["var3"]))
    ss$Q <- rbind(cbind(ss$V, 0), 0)
    ss$a0 <- c(allpars["a01"], allpars["a02"], allpars["a03"])
    ss$P0 <- diag(c(allpars["P01"], allpars["P02"], allpars["P03"]))

  } else
  if (x@model == "level+drift+AR2")
  {
    ss$Vid <- NULL
    ss$Qid <- NULL

    #NOTE code of original paper uses allpars["vari"]^2 
    #i.e., the model is parameterized in terms of standar deviations
    #here in term of variances, also relevant in transPars

    p <- 2
    ss$T <- rbind(c(1, 0, 0), c(0, allpars[c("phi1", "phi2")]), c(0,1,0))      
    ss$H <- allpars["var1"]
    ss$V <- diag(c(allpars["var2"], allpars["var3"]))
    ss$V[1,2] <- ss$V[2,1] <- allpars["cov23"]
    ss$Q <- rbind(cbind(ss$V, 0), 0)
    ss$Q[1,2] <- ss$Q[2,1] <- allpars["cov23"]
    ss$a0 <- c(allpars["a01"], allpars["a02"], allpars["a03"])
    #ss$P0 <- diag(c(allpars["P01"], allpars["P02"], allpars["P03"]))

    ss$mu <- c(allpars["drift"], 0, 0)

##NOTE allpars["P0i"] are ignored
#http://www.trigconsulting.co.uk/gauss/man_algebra.html
P0 <- diag(3)
P0[1,1] <- 100000000
tmp <- solve(diag(4) - kronecker(ss$T[-1,-1], ss$T[-1,-1])) %*% as.vector(ss$Q[-1,-1])
P0[-1,-1] <- matrix(tmp, nrow=2, ncol=2)
##FIXME non-zero values outside the diagonal ??? maybe, because cov23 != 0
#but the non-zero off-diagonal element seems to be located in a place other that than the cell related to cov23
ss$P0 <- P0

  } else ##NOTE currently this option is not used
  {
    ss$Z <- ss.fill(ss$Z, allpars, allnames)
    if (is.list(ss$T)) 
    {
      #mT <- lapply(ss$T, FUN = eval, envir = as.list(allpars))
      if (is.matrix(ss$T[[1]])) {
        mT1 <- ss$T[[1]]
        ss$T[[1]] <- NULL
        mT2 <- lapply(ss$T, FUN = eval, envir = as.list(allpars))
        mT2 <- matrix(unlist(mT2), nrow = length(mT2)/2)
        aux <- diag(0, nrow(mT1))
        ss$T <- rbind(cbind(mT1, aux), cbind(aux, mT2))
      } else {
          mT <- lapply(ss$T, FUN = eval, envir = as.list(allpars))
          ss$T <- matrix(unlist(mT), nrow = ncol(ss$Z))
      }
    } else
    ss$T <- ss.fill(ss$T, allpars, allnames)
    ss$H <- ss.fill(ss$H, allpars, allnames)
    ss$R <- ss.fill(ss$R, allpars, allnames)
    ss$V <- ss.fill(ss$V, allpars, allnames)
    ss$Q <- ss$R %*% ss$V %*% t(ss$R)

    ss$a0 <- ss.fill(ss$a0, allpars, allnames)
    ss$P0 <- ss.fill(ss$P0, allpars, allnames)
  }

  if (P0cov)
    ss$P0[] <- ss$P0[1]

  pnms <- names(pars)
  ss$parsNames <- pnms
  ss$id <- list(
    var = grep("^var\\d{1,2}$", pnms),
    phi = grep("^phi\\d{1,2}$", pnms),
    a0 = grep("^a0\\d{1,2}$", pnms))

  if (!getgrad)
  {
    # define a vector arranged to be used by "KF.deriv"
    # this output is intended to save doing the arrangements "KF.deriv",
    # which is likely to be called several times and is already more demanding
    #
    # when model@transPars is NULL 
    # the derivative with respect to the variance parameter 
    # is retured, 1, (not the derivative of transPars with respect the auxiliary
    # parameter, which would be 0 if no transformation or reparameterization is used);
    # this way is more convenient when this output used by "KF.deriv"
    
    tpg.var <- rep(1, length(ss$id$var))
    names(tpg.var) <- pnms[ss$idvar]
    #ss$tpGradient <- list(var = tpg.var, phi1 = NULL, phi2 = NULL)
    ss$tpGradient <- list(var = tpg.var)

  } else {
    #ss$tpGradient <- tpGradient
    if (is.matrix(tpGradient)) # then "phi1" and "phi2" are defined in "pars"
    {
      ss$tpGradient <- list(var = diag(tpGradient)[ss$id$var],
        phi1 = tpGradient["phi1",c("phi1","phi2")], 
        phi2 = tpGradient["phi2",c("phi1","phi2")])
    } else 
      ss$tpGradient <- list(var = tpGradient)
  }

  class(ss) <- "stsmSS"
  ss
})

setMethod("show", "stsm",
  function(object) {
    cat(paste("Model:", object@model, "\n\n"))
    p <- get.pars(object, rescale = !is.null(object@cpar))
    print(p)
    invisible()
  }
)
