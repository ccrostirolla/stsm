
maxlik.fd.optim <- function(m,  
  barrier = list(type = c("1", "2"), mu = 0), inf = 99999, 
  method = c("BFGS", "L-BFGS-B", "Nelder-Mead", "CG", "SANN"), 
  gr = c("analytical", "numerical"), optim.control = list(), hessian = TRUE)
{
  mcall <- match.call()

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

    if(gr0 == "analytical") {
      gr <- function(x, ...) { mcloglik.fd.grad(x = x, ...) }
    } else gr <- NULL
  }

  method <- match.arg(method)[1]
  if (method == "L-BFGS-B")
  {
    res <- optim(par = m@pars, fn = fcn, gr = gr, 
      model = m, barrier = bargs, inf = inf, method = method, 
      lower = m@lower, upper = m@upper, 
      control = optim.control, hessian = hessian)
  } else
  {
    res <- optim(par = m@pars, fn = fcn, gr = gr, 
      model = m, barrier = bargs, inf = inf, method = method, 
      control = optim.control, hessian = hessian)
  }

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

  res2 <- c(list(call = mcall,
    init = pars0, pars = m@pars, model = m, loglik = -res$value,
    convergence = res$convergence, iter = res$counts, message = res$message,
    Mpars = NULL, steps = NULL), list(ls.iter = NULL, ls.counts = NULL),
    hessian = res$hessian)
  class(res2) <- "stsmFit"
  res2
}
