
maxlik.td.optim <- function(m,
  KF.version = eval(formals(KFKSDS::KalmanFilter)$KF.version),
  KF.args = list(), check.KF.args = TRUE,
  barrier = list(type = c("1", "2"), mu = 0), inf = 99999, 
  method = c("BFGS", "L-BFGS-B", "Nelder-Mead", "CG", "SANN", "AB-NM"),
  gr = c("numerical", "analytical"), optim.control = list(), hessian = FALSE)
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

  KF.version <- match.arg(KF.version)[1]
  gr0 <- match.arg(gr)[1]
  if (gr0 == "analytical" && KF.version != "KFKSDS")
  {
    warning("analytical gradient is not available for 'KF.version = ", 
      KF.version, "'. Changed to 'KF.version = KFKSDS'.")
    KF.version <- "KFKSDS"
  }

  if (check.KF.args)
    KF.args <- make.KF.args(char2numeric(m), KF.version, KF.args) 

  if (gr0 == "analytical")
  {
    if (!is.null(m@cpar)) {
      warning("numerical gradient was used (analytical gradient is not available for ", 
        sQuote("KFconvar"), ").")
      gr <- NULL
    } else gr <- function(x, ...) { mloglik.td.grad(x = x, ...) }
  } else gr <- NULL

  if (is.null(m@cpar))
  {
    if(gr0 == "analytical") {
      gr <- function(x, ...) { mloglik.fd.grad(x = x, ...) }
    } else gr <- NULL
  } else {

    if(gr0 == "analytical") {
      gr <- function(x, ...) { mcloglik.fd.grad(x = x, ...) }
    } else gr <- NULL
  }

  method <- match.arg(method)

  if (method == "L-BFGS-B")
  {
    res <- optim(par = m@pars, fn = mloglik.td, gr = gr, 
      model = m, KF.version = KF.version, KF.args = KF.args, 
      barrier = bargs, inf = inf, #load.package = FALSE, 
      method = method, lower = m@lower, upper = m@upper, 
      control = optim.control, hessian = hessian)

  } else 
  if (method == "BFGS") {
    res <- optim(par = m@pars, fn = mloglik.td, gr = gr, 
      model = m, KF.version = KF.version, KF.args = KF.args, 
      barrier = bargs, inf = inf, #load.package = FALSE, 
      method = method, 
      control = optim.control, hessian = hessian)
  } else
  if (method == "AB-NM")
  {
    res <- constrOptim(theta = m@pars, f = mloglik.td, grad = NULL,
      model = m, KF.version = KF.version, KF.args = KF.args,
      ui = diag(length(m@pars)), ci = m@lower, mu = 1e-04, 
      method = "Nelder-Mead")
  }

  pars0 <- c(m@pars, m@cpar)
  m@pars <- res$par

  if (!is.null(m@cpar))
  {
    # the barrier is not added to the final likelihood value
    tmp <- KFconvar(m, P0cov = KF.args$P0cov, debug = TRUE)
    m@cpar[] <- tmp$cpar
##FIXME ver si tmp$mll concide con el valor final de 'optim' o se actualiza algo
    res$value <- -tmp$mll
  } else {
    # the barrier is not added to the final likelihood value
    val <- logLik(object = m, domain = "time", 
      td.args = c(KF.version = KF.version, KF.args), 
      barrier = list(mu = 0), inf = 99999)
  }

  res2 <- c(list(call = mcall,
    init = pars0, pars = m@pars, model = m, loglik = res$value,
    convergence = res$convergence, iter = res$counts, message = res$message,
    Mpars = NULL, steps = NULL), list(ls.iter = NULL, ls.counts = NULL),
    hessian = res$hessian)
  class(res2) <- "stsmFit"
  res2
}
