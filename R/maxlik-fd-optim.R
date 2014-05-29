
maxlik.fd.optim <- function(m, xreg = NULL,
  barrier = list(type = c("1", "2"), mu = 0), inf = 99999, 
  method = c("BFGS", "L-BFGS-B", "Nelder-Mead", "CG", "SANN"), 
  gr = c("analytical", "numerical"), optim.control = list(), hessian = TRUE)
{
  mcall <- match.call()
  method <- match.arg(method)

##FIXME pass xreg already differenced to optim to reduce computations

  bargs <- list(type = "1", mu = 0)
  if (barrier$mu > 0) #&& !is.null(barrier$mu)
  {
    nbargs0 <- names(bargs)
    bargs[nbargs <- names(barrier)] <- barrier
    if (length(nobargs <- nbargs[!nbargs %in% nbargs0]))
      warning("unknown names in 'barrier': ", paste(nobargs, collapse = ", "))
  }

  gr0 <- match.arg(gr)[1]
  if(bargs$mu > 0 && gr0 == "analytical")
    warning("barrier term was ignored.")

  if (is.null(m@cpar))
  {
    fcn <- mloglik.fd

    if(gr0 == "analytical") {
      gr <- function(x, ...) { mloglik.fd.grad(x = x, ...) }
    } else gr <- NULL
  } else {
    fcn <- mcloglik.fd
##FIXME TODO mcloglik with mcloglik.fd

    if(gr0 == "analytical") {
      gr <- function(x, ...) { mcloglik.fd.grad(x = x, ...) }
    } else gr <- NULL
  }

  if (!is.null(xreg))
  {
##FIXME see svd()

    stopifnot(NROW(xreg) == length(m@y))

    #if (length(dim(xreg)) == 1L)
    #  xreg <- cbind(xreg)
    #
    # if length(dim(xreg)) == 1L colnames will be null
    if (is.null(xregnms <- colnames(xreg)))
    {
      # column names are necessary to match coefficients in slots "pars" and/or "nopars"
      # which may follow a different order than the arrangement of variables 
      # by columns in "xreg"
      stop("column names must be given to regressor variables ",
        "defined in argument ", sQuote("xreg")) 
    }

    #xregcoefs <- model@pars[grepl("^xreg", names(model@pars))]
    allpars <- c(m@pars, m@nopars)
    id.isna <- is.na(match(xregnms, names(allpars)))
    if (all(id.isna))
    { # xreg coefficients are not defined in "stsm" object "m"
      # define them here as parameters to be estimated (slot "pars")
      ncxreg <- ncol(xreg)
      xregpars <- rep(1, ncxreg)
      names(xregpars) <- xregnms
      m@pars <- c(m@pars, xregpars)
      tmp <- rep(Inf, ncxreg)
      names(tmp) <- xregnms
      m@upper <- c(m@upper, tmp)
      m@lower <- c(m@lower, -tmp)      
    } else 
    if (all(!id.isna))
    { # nothing to do since
      # xreg coefficients are defined in "stsm" object "m" and their names match
      # the column names of the regressor variables passed in "xreg";
      # it is checked whether lower and upper bounds to be used by 
      # method "L-BFGS-B" are defined
      if (method == "L-BFGS-B")
      {
        parsnms <- names(m@pars)
        tmpnms <- intersect(xregnms, parsnms)
        if (length(tmpnms) > 0)
        {
          ##NOTE
          # if no initial values for the coefficients were not specified 
          # then they will be initialized above and the lower and upper bounds
          # would also be defined;
          # unless specific initial values are desired, it is more convenient
          # not to define the parameters in the slot "pars" and let this function 
          # do all the arrangements if "xreg" is not NULL
          if (!all(tmpnms %in% names(m@lower)) || !all(tmpnms %in% names(m@upper)))
            stop("method ", sQuote("L-BFGS-B"), " requires lower and upper bounds.\n",
            "Lower and upper bounds for the coefficients of the regressor variables ",
            "must be defined in the slots ", sQuote("lower"), " and ", sQuote("upper"), 
            " of the ", sQuote("stsm"), " object passed as input in argument ", 
            sQuote("m"), ".")
        }
      }
    } else 
    if (any(id.isna))
    { # some of the xreg coefficients seem to be defined because their names
      # match the column names of "xreg", but for some others there is no match
      # do not make a guess about the input, stop for safety
      stop("some, but not all, names of the parameters in object ", sQuote("m"), 
        " match the column names of ", sQuote("xreg"), 
        ". Either all or neither of the coefficients should be defined in the ", 
        sQuote("stsm"), "input model ", sQuote("m"), ".")
    }

    ##NOTE
    #xreg is overwritten 
    xreg <- m@fdiff(xreg, frequency(m@y))
  }
  
  if (method == "L-BFGS-B")
  {
    res <- optim(par = m@pars, fn = fcn, gr = gr, 
      model = m, xreg = xreg, barrier = bargs, inf = inf, method = method, 
      lower = m@lower, upper = m@upper, 
      control = optim.control, hessian = hessian)
  } else
  {
    res <- optim(par = m@pars, fn = fcn, gr = gr, 
      model = m, xreg = xreg, barrier = bargs, inf = inf, method = method, 
      control = optim.control, hessian = hessian)
  }
##FIXME see add if (method == "AB-NM") as in maxlik.td.optim()

  pars0 <- m@pars
  m@pars <- res$par

  if (!is.null(m@cpar))
  {
    if (is.null(m@ssd)) {
      pg <- Mod(fft(m@diffy))^2 / (2 * pi * length(m@diffy))
    } else
      pg <- m@ssd

    n <- length(m@diffy)
    nh <- n / 2
    pi2 <- 2 * pi
    sgf <- stsm.sgf(m, FALSE, FALSE, FALSE)$sgf
    m@cpar[] <- (pi2 / n) * sum(pg / sgf)
    res$value <- -nh * (log(pi2) + 1) - nh * log(m@cpar) - 0.5 * sum(log(sgf))
  } else {
    # the barrier is not added to the final likelihood value
    if (barrier$mu > 0)
      res$value <- logLik(object = m, domain = "frequency", barrier = list(mu = 0))
  }

  # "xregcoefs" are not considered part of a "stsm" object;
  # isolate the coefficients related to regressor variables
  # (only those in slot "pars)
  # for convenience they are treated separately by other functions  

  if (!is.null(xreg))
  {
    id <- na.omit(match(xregnms, names(m@pars)))
    if (length(id) > 0) 
    {
      xregcoefs <- m@pars[id]
      m@pars <- m@pars[-id]
      if (length(id <- match(xregnms, names(m@lower))) > 0)
        m@lower <- m@lower[-id]
      if (length(id <- match(xregnms, names(m@upper))) > 0)
        m@upper <- m@upper[-id]
    }
    xreg <- list(coefs = xregcoefs, stde = NULL)
  } #else 
    #xregcoefs <- NULL

  res2 <- c(list(call = mcall, model = m, 
    init = pars0, pars = m@pars, xreg = xreg, loglik = -res$value,
    convergence = res$convergence, iter = res$counts, message = res$message,
    Mpars = NULL, steps = NULL), list(ls.iter = NULL, ls.counts = NULL),
    list(hessian = res$hessian))
  class(res2) <- "stsmFit"
  res2
}
