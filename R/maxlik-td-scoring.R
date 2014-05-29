
maxlik.td.scoring <- function(m, xreg = NULL, step = NULL, 
  KF.args = list(), check.KF.args = TRUE,
  ls = list(type = "optimize", tol = .Machine$double.eps^0.25, cap = 1),
  control = list(maxit = 100, tol = 0.001, trace = FALSE, silent = FALSE), 
  debug = FALSE)
{
  barrier <- list(type = "1", mu = 0)
  mcall <- match.call()

  if (!is.null(step) && (step <= 0 || !is.numeric(step)))
    stop("'step' must be a numeric positive value.")

  mfcn.step <- function(x, m, pd, barrier)
  {
    m <- set.pars(m, get.pars(m) + x * pd)

    -logLik(object = m, domain = "time", xreg = NULL,
      td.args = list(KF.version = "KFKSDS", P0cov = FALSE, t0 = 1), FALSE,
      barrier = list(mu = 0), inf = 99999)
  }

  if (check.KF.args)
    KF.args <- make.KF.args(char2numeric(m), "KFKSDS", KF.args) 

  cargs <- list(maxit = 100, tol = 0.001, trace = FALSE, silent = FALSE)
  ncargs0 <- names(cargs)
  cargs[ncargs <- names(control)] <- control
  if (length(nocargs <- ncargs[!ncargs %in% ncargs0]))
    warning("unknown names in 'control': ", paste(nocargs, collapse = ", "))
  control <- cargs

  info <- "expected"
  use.IM <- TRUE

  lsres <- list(ls.iter = NULL, ls.counts = c("fcn" = 0, "grd" = 0))
  lsintv <- c(0, ls$cap)

  pars0 <- m@pars
  convergence <- FALSE
  conv.e1 <- FALSE
  iter <- 0

  # regressor variables (define constant variables)

  ##NOTE
  # this block is based on same block in function maxlik.fd.scoring()

  if (!is.null(xreg))
  {
    # check for safety if a single regressor is passed as a vector
    if (is.null(dim(xreg)))
      xreg <- data.matrix(xreg)
    stopifnot(nrow(xreg) == length(m@y))

    dxreg <- m@fdiff(xreg, frequency(m@y))
    y0 <- m@y

    xregcoefs <- coef(lm.fit(x = dxreg, y = m@diffy))

    m@y <- y0 - xreg %*% xregcoefs

  } else
    xregcoefs <- NULL

  # storage matrices for tracing information

  Mpars <- if (control$trace) rbind(c(m@pars, xregcoefs), 
    matrix(nrow = control$maxit + 1, 
      ncol = length(m@pars) + max(0, ncol(xreg)))) else NULL
  steps <- if (control$trace) rep(NA, control$maxit + 2) else NULL

if (debug)
{
  val <- logLik(object = m, domain = "time", 
    td.args = list(P0cov = FALSE, t0 = 1, KF.version = "KFKSDS"), 
    barrier = list(mu = 0), inf = 99999)
  cat(paste("\niter =", iter, "logLik =", round(val, 4), "\n"))
  print(get.pars(m))
}

  # begin iterative process

j1 <- 1

  while (!(convergence || iter > control$maxit))
  {
    tmp <- mloglik.td.deriv(m, TRUE, use.IM, KF.args, "1", NULL)

    # gradient

    G <- -tmp$gradient

    # information matrix

    IM <- tmp$infomat

if (debug)
{
  eg <- eigen(IM, only.values = TRUE)$values
  if (min(eg) <= 0) {
    warning("IM is not a positive definite matrix.")
    print(eg)
  }
}

    # direction vector

    pd <- drop(solve(IM) %*% G)

    # step size (choose the step that maximizes the increase in 
    # the log-likelihood function given the direction vector 'pd')

    if (is.null(step))
    {
      lsintv[2] <- step.maxsize(m@pars, m@lower, m@upper, pd, ls$cap)

      lsout <- switch(ls$type,

      "optimize" = optimize(f = mfcn.step, interval = lsintv, 
        maximum = FALSE, tol = ls$tol, m = m, pd = pd, 
        barrier = barrier),

      "brent.fmin" = Brent.fmin(a = 0, b = lsintv[2], 
        fcn = mfcn.step, tol = ls$tol, m = m, pd = pd, 
        barrier = barrier))

      lambda <- lsout$minimum
      lsres$ls.iter <- c(lsres$ls.iter, lsout$iter)
      lsres$ls.counts <- lsres$ls.counts + lsout$counts

    } else
    if (is.numeric(step))
      lambda <- step

    # update (new guess)

    pars.old <- m@pars
    pars.new <- pars.old + lambda * pd
    m <- set.pars(m, pars.new)

    # stopping criteria

    if (sqrt(sum((pars.old - m@pars)^2)) < control$tol)
    {
      convergence <- TRUE
    } else
    if (sum(abs(pars.old - m@pars) > control$tol) == 1)
      j1 <- j1 + 1
    if (j1 == 10) {
      convergence <- TRUE
      conv.e1 <- TRUE
      if (!control$silent)
        warning(paste("Possible convergence problem.\n",
        "Over 10 iterations, failure to concergence was caused by just 1 parameter",
        "(not necessarily the same parameter)."))
    }

    # coefficients of regressor variables

    if (!is.null(xreg))
    {
      ss <- char2numeric(m)
      kf <- KF(y0, ss)
      ks <- KS(y0, ss, kf)
      dy <- m@fdiff(rowSums(ks$ahat), frequency(m@y))
      xregcoefs <- coef(lm.fit(x = dxreg, y = dy))

      if (!convergence)
      {
        m@y <- y0 - xreg %*% xregcoefs
      } 
    }

    # trace

    iter <- iter + 1

    if (control$trace)
    {
      Mpars[iter+1,] <- c(m@pars, xregcoefs)
      steps[iter] <- lambda
    }

if (debug)
{
  val <- logLik(object = m, domain = "time", 
    td.args = list(P0cov = FALSE, t0 = 1, KF.version = "KFKSDS"), 
    barrier = list(mu = 0), inf = 99999)
  cat(paste("\niter =", iter, "logLik =", round(val, 4), "\n"))
  print(get.pars(m))
}

if (debug && !is.null(m@lower) && !is.null(m@upper))
{
  check.bounds(m)
}

  }

  if (!control$silent && !convergence)
    warning(paste("Possible convergence problem.",
    "Maximum number of iterations reached."))

  if (control$trace)
  {
    Mpars <- na.omit(Mpars)
    attr(Mpars, "na.action") <- NULL
    steps <- na.omit(steps)
    attr(steps, "na.action") <- NULL
  }

  # output

  # the barrier term is not added to the final likelihood value
  # xreg is set to NULL and m@y <- y0 - xreg %*% xregcoefs (done above)
  # is used, xreg coefficients are not part of the "stsm" model
  val <- logLik(object = m, domain = "time", xreg = NULL,
    td.args = c(KF.args, KF.version = "KFKSDS"), check.td.args = TRUE,
    barrier = list(mu = 0), inf = 99999)

  if (!is.null(xreg))
  {
    tmp <- coef(summary(lm(dy ~ 0 + dxreg)))

    xregcoefs <- tmp[,"Estimate"]
    xregstde <- tmp[,"Std. Error"]
##FIXME return xregtstat also in function maxlik.fd.scoring
    xregtstat <- tmp[,"t value"]
    names(xregcoefs) <- names(xregstde) <- names(xregtstat) <- colnames(xreg)
    xreg <- list(coef = xregcoefs, stde = xregstde, tstat = xregtstat)

    m@y <- y0
  }

  if (conv.e1)
    convergence <- FALSE

  res <- c(list(call = mcall, model = m, 
    init = pars0, pars = m@pars, xreg = xreg, loglik = val,
    convergence = convergence, iter = iter, message = "",
    Mpars = Mpars, steps = steps), lsres)
  class(res) <- "stsmFit"
  res
}
