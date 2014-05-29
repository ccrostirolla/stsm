
##FIXME TODO analytical derivatives when xreg!=NULL (see vcov expression in H89)

mloglik.fd <- function(x, model, xreg = NULL,
  barrier = list(type = c("1", "2"), mu = 0), inf = 99999)
{
  if (!missing(x)) if (!is.null(x))
    model@pars <- x

  n <- length(model@diffy)
  pi2 <- 2 * pi

  if (!is.null(xreg))
  {
    #stopifnot((NROW(xreg) == length(model@y)) || (NROW(xreg) == length(model@diffy)))

    #xregcoefs <- model@pars[grepl("^xreg", names(model@pars))]
    allpars <- c(model@pars, model@nopars)
    id <- match(colnames(xreg), names(allpars))
    if (any(is.na(id)))
      stop("column names of ", sQuote("xreg"), " do not match the names of the parameters")
    xregcoefs <- allpars[id]

    ##NOTE
    #xreg is overwritten 
    #an object called for instance "dxreg" could be created but it would involve 
    #two objects with the same content;
    #no additional arguments added to function, xreg or differenced xreg is 
    #deduced from the length/nrows; e.g., xreg is passed already appropriately 
    #differenced by function maxlik.fd.optim
    if (NROW(xreg) == length(model@y)) {
      xreg <- model@fdiff(xreg, frequency(model@y))
    }

    model@diffy <- model@diffy - xreg %*% cbind(xregcoefs)
    #in general update only if (!is.null(model@ssd)
    model@ssd <- as.vector(Mod(fft(model@diffy))^2 / (pi2 * n))
  }

  if (is.null(model@ssd)) {
    pg <- Mod(fft(model@diffy))^2 / (pi2 * n)
  } else
    pg <- model@ssd

  # spectral generating function of the stationary differenced data

  sgf <- stsm.sgf(model, FALSE, FALSE, FALSE)$sgf

  # minus log-likelihood

  if (sgf[1] == 0)
  {
    sgf <- sgf[-1]
    pg <- pg[-1]
  }

  mll <- 0.5 * (n * log(pi2) + sum(log(sgf))) + pi * sum(pg/sgf)

  # barrier term (optional)

  if (barrier$mu != 0)
  {
    bar <- barrier.eval(model, barrier$type, barrier$mu, 
      FALSE, FALSE)$barrier
  } else bar <- 0

  mll <- mll + bar

  if (!is.finite(mll))
    mll <- sign(mll) * inf

  if (is.na(mll))
    mll <- abs(inf) #minus log-lik, penalize with a large value

  mll
}

##FIXME TODO for non null xreg or at least adapt vcov() to it

mloglik.fd.deriv <- function(model, 
  gradient = TRUE, hessian = TRUE, infomat = TRUE, modcovgrad = TRUE,
  barrier = list(type = c("1", "2"), mu = 0),
  version = c("2", "1"))
{
  version <- match.arg(version)[1]

  dtpars <- if (version == "1") TRUE else FALSE
  if (!is.null(model@transPars) && version == "2" && any(gradient, hessian, modcovgrad))
    dtrans <- transPars(model, gradient = gradient, hessian = hessian)

  pi2 <- 2 * pi

  tmp <- stsm.sgf(model, TRUE, (hessian || modcovgrad), dtpars)
  sgf <- tmp$sgf
  sgf.d1 <- tmp$gradient
  sgf.d2 <- tmp$hessian

  if (any(c(gradient, hessian, infomat, modcovgrad)) && barrier$mu != 0)
    bar <- barrier.eval(model, barrier$type, barrier$mu, 
      gradient, (hessian || infomat))

  if (is.null(model@ssd)) {
    pg <- Mod(fft(model@diffy))^2 / (pi2 * length(model@diffy))
  } else
    pg <- model@ssd

  if (sgf[1] == 0)
  {
    sgf <- sgf[-1]
    sgf.d1 <- sgf.d1[-1,]
    sgf.d2 <- sgf.d2[,,-1]
    pg <- pg[-1]
  }

  if (hessian || infomat || modcovgrad)
    sgf.sq <- sgf^2
  if (gradient || hessian)
  {
    pgog <- as.vector(pg/sgf)
    pipgog <- pi * pgog
    if (gradient)
      gd1og <- as.matrix(sgf.d1 / sgf)
  }

  # first order derivatives

  if (gradient || (!is.null(model@transPars) && any(hessian, infomat, modcovgrad)))
  {
    d10 <- -0.5 * colSums((2 * pipgog - 1) * gd1og)

    if (!is.null(model@transPars) && version == "2") {
      d1 <- d10 * dtrans$gradient
    } else d1 <- d10

    if (barrier$mu != 0)
    {      
      if (!is.null(model@lower)) {
        d1 <- d1 + barrier$mu * bar$dl1
      }
      if (!is.null(model@upper)) {
        d1 <- d1 + barrier$mu * bar$du1
      }
    }
  } else d1 <- NULL

  # second order derivatives

  if (hessian)
  {
    tmp1 <- (0.5 - 2 * pipgog) / sgf.sq
    gd1cp <- apply(sgf.d1, MARGIN = 1, FUN = tcrossprod)
    d2 <- matrix(-colSums(tmp1 * t(gd1cp)), nrow = length(model@pars))
    tmp2 <- (-0.5 + pipgog) / sgf

if (model@model %in% c("cycle", "trend-cycle")) {
  tmp <- matrix(0, nrow(d2), nrow(d2))
  for(i in seq(dim(sgf.d2)[3]))
    tmp <- tmp + tmp2[i] * sgf.d2[,,i]
  d2 <- d2 - tmp
} else {
    diag(d2) <- diag(d2) - colSums(tmp2 * sgf.d2)
}

    if (!is.null(model@transPars) && version == "2")
      d2 <- tcrossprod(dtrans$gradient) * d2 + d10 * dtrans$hessian

    if (barrier$mu != 0)
    {
      if (!is.null(model@lower))
        diag(d2) <- diag(d2) + barrier$mu * bar$dl2
      if (!is.null(model@upper))
        diag(d2) <- diag(d2) + barrier$mu * bar$du2
    }
  } else d2 <- NULL

  if (infomat || modcovgrad)
  {
    if (hessian)
    {
      gcov <- matrix(0.5 * colSums(t(gd1cp) / sgf.sq), nrow = ncol(d2))
    } else
    {
      gcov <- 0
      for (i in seq(along = sgf))
        gcov <- gcov + tcrossprod(sgf.d1[i,]) / sgf.sq[i]
      gcov <- 0.5 * gcov
    }

    if (!is.null(model@transPars) && version == "2") {
      if (modcovgrad) 
        gcovmod <- gcov
      gcov <- tcrossprod(dtrans$gradient) * gcov
    } else if (modcovgrad) gcovmod <- gcov

    if (barrier$mu != 0)
    {
      if (!is.null(model@lower))
        diag(gcov) <- diag(gcov) + barrier$mu * bar$dl2
      if (!is.null(model@upper))
        diag(gcov) <- diag(gcov) + barrier$mu * bar$du2
    }
  } else gcov <- NULL

  if (modcovgrad)
  {
    if (!hessian)
      tmp2 <- (-0.5 + pipgog) / sgf

    if (model@model %in% c("cycle", "trend-cycle"))
    { # for devel version
      tmp <- matrix(0, nrow(d2), nrow(d2))
      for(i in seq(dim(sgf.d2)[3]))
        tmp <- tmp + tmp2[i] * sgf.d2[,,i]
      gcovmod <- gcovmod + tmp
    } else {
        diag(gcovmod) <- diag(gcovmod) + colSums(tmp2 * sgf.d2)
    }

    if (!is.null(model@transPars) && version == "2")
      gcovmod <- tcrossprod(dtrans$gradient) * gcovmod - d10 * dtrans$hessian

    if (barrier$mu != 0)
    {
      if (!is.null(model@lower))
        diag(gcovmod) <- diag(gcovmod) + barrier$mu * bar$dl2
      if (!is.null(model@upper))
        diag(gcovmod) <- diag(gcovmod) + barrier$mu * bar$du2
    }
  } else gcovmod <- NULL

  list(gradient = d1, hessian = d2, infomat = gcov, modcovgrad = gcovmod)
}

mloglik.fd.grad <- function(x, model,
  barrier = list(type = c("1", "2"), mu = 0),
  inf, xreg)
#arguments 'inf' and 'xreg' are not used here but it is needed when used within optim 
#'maxlik.fd.optim' where this function is passed as the gradient
{
  if (!missing(x)) if (!is.null(x))
    model@pars <- x

  mloglik.fd.deriv(model = model, 
    gradient = TRUE, hessian = FALSE, infomat = FALSE, modcovgrad = FALSE,
    barrier = barrier, version = "2")$gradient
}
