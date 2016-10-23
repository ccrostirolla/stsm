
setGeneric(name = "transPars", def = function(x, 
##FIXME see rename as var.type
  type = c("square", "StructTS", "exp", "exp10sq", 
    "square+v1", "StructTS+v1", "exp+v1", "exp10sq+v1",
    "square+v2", "StructTS+v2", "exp+v2", "exp10sq+v2",
    "square+stats", "StructTS+stats", "exp+stats", "exp10sq+stats",
    "square+none", "StructTS+none", "exp+none", "exp10sq+none",
    "level+drift+AR2"), 
  #ar.type = c("v1", "v2", "stats", "none"),
  gradient = FALSE, hessian = FALSE, 
  rp, sclrho = 1.7, sclomega = 1.7, ftrans = NULL, ...) {
    standardGeneric("transPars") })

setMethod(f = "transPars", signature = "numeric", 
  definition = function(x, 
    type = eval(formals(transPars)$type),
    #ar.type = eval(formals(transPars)$ar.type), 
    gradient = FALSE, hessian = FALSE,
    rp, sclrho = 1.7, sclomega = 1.7, ftrans = NULL, ...) {

  tpStructTS <- function(pars, rp)
  {
    if (length(idvar) > 0) {
      tpars[idvar] <- pars[idvar] * rp
    } #else warning("No changes done by 'transPars'.")
    # changes may be done in other parameters, e.g. 'phi1', 'phi2'

    if (gradient)
    {
      if (length(idvar) > 0)
        d1[idvar] <- rp
    }
    #if (hessian)
    #  d2[] <- 0 #array(0, dim = c(p, p))
    list(pars = tpars, d1 = d1, d2 = d2)
  }

  tpsq <- function(pars)
  {
    if (length(idvar) > 0) {
      tpars[idvar] <- pars[idvar]^2
    } #else warning("No changes done by 'transPars'.")

    if (gradient)
    {
      if (length(idvar) > 0)
        d1[idvar] <- 2 * pars[idvar]
    }
    if (hessian) {
      #d2[] <- 0 #array(0, dim = c(np, np)) #diag(2, np, np)
      diag(d2)[idvar] <- 2
    }
    list(pars = tpars, d1 = d1, d2 = d2)
  }

  tpexp <- function(pars)
  {
    if (length(idvar) > 0) {
      tpars[idvar] <- exp(pars[idvar])
    } #else warning("No changes done by 'transPars'.")

    if (gradient)
    {
      if (length(idvar) > 0)
        d1[idvar] <- tpars[idvar]
    }
    if (hessian) {
      #d2[] <- 0 #array(0, dim = c(np, np)) #diag(2, np, np)
      diag(d2)[idvar] <- tpars[idvar]
    }
    list(pars = tpars, d1 = d1, d2 = d2)
  }

  tpexp10sq <- function(pars)
  {
    if (length(idvar) > 0) {
      tpars[idvar] <- (exp(-pars[idvar])/10)^2
    } #else warning("No changes done by 'transPars'.")

##FIXME TODO
    if (gradient)
    {
      d1 <- NULL
    }
    if (hessian) {
      d2 <- NULL
    }
    list(pars = tpars, d1 = d1, d2 = d2)
  }

  if (is.function(ftrans))
    return(ftrans(x = x, gradient = gradient, hessian = hessian, ...))

  type <- match.arg(type)[1]
  tpars <- pars <- x

  label <- "var"
  np <- length(pars)
  pnms <- names(pars)
  idvar <- grep("^var|P0\\d{1,2}$", pnms, value = FALSE)
  idphi <- grep("^phi\\d", pnms)

  tmp <- strsplit(type, "\\+")[[1]]
  if (length(idphi) > 0)
  {    
    type.phi <- tmp[2]
    if (is.na(type.phi)) {
      type.phi <- "none"
      # give warning because the user did not explicity specify the option "+none"
      warning("AR parameters are *not* transformed to remain in the region of stationarity.")
    } else 
      type <- tmp[1]
  } else 
  if (!is.na(tmp[2]))
    warning(sQuote(tmp[2]), " in argument ", sQuote("type"), " is ignored ", 
    "(AR coefficients are not defined in the slot ", sQuote("pars"), ")")

  if (type == "StructTS")
    stopifnot(!is.null(rp))

  if (gradient)
  {
    d1 <- rep(NA, np)
    names(d1) <- pnms
    d1[grep("^a0\\d", pnms)] <- 0
  } else
    d1 <- NULL
  
  if (hessian)
  {
    d2 <- matrix(0, np, np)
    rownames(d2) <- colnames(d2) <- pnms
  } else
    d2 <- NULL

  if (type == "level+drift+AR2")
  {
##FIXME devel
##FIXME vcov and all vars except var1 must be defined in "pars", not in "nopars"

##FIXME see allpars (see if some vars are in pars and others in nopars

    if ("cov23" %in% pnms)
    {
      if ("var1" %in% pnms)
        tpars["var1"] <- exp(-pars["var1"])
        #tpars["var1"] <- exp(-pars["var1"])^2

      varcov <- crossprod(rbind(c(exp(-pars["var2"]), 0), c(pars["cov23"], exp(-pars["var3"]))))
      #NOTE original code parameterizes the model in terms of standar deviations
      #here in terms of variance parameters (do not take square root)
      #tpars["var2"] <- sqrt(varcov[1,1])
      #tpars["var3"] <- sqrt(varcov[2,2])
      tpars["var2"] <- varcov[1,1]
      tpars["var3"] <- varcov[2,2]
      tpars["cov23"] <- varcov[2,1]
    } else {
      tpars[idvar] <- exp(-pars[idvar])
      #tpars[idvar] <- exp(-pars[idvar])^2
    }

    aaa <- pars[idphi[1]] / (1 + abs(pars[idphi[1]]))
    ccc <- (1 - abs(aaa))*pars[idphi[2]] / (1 + abs(pars[idphi[2]])) + abs(aaa) - aaa^2
    tpars[idphi[1]] <- 2*aaa
    tpars[idphi[2]] <- -(aaa^2 + ccc)

#opabsx <- 1 + abs(pars[idphi])
#aux <- pars[idphi] / opabsx
#print(c(sum(aux), -prod(aux)))
#tpars[idphi] <- c(sum(aux), -prod(aux))

    ##NOTE finish here
    return(list(pars = tpars, gradient = d1, hessian = d2))
  }

##FIXME TODO none (none transformation for vars, return here the list)
#list(pars = pars, gradient = d1, hessian = d2)
#i.e., possibility not to transform vars but transform phi

  res <- switch(type, "StructTS" = tpStructTS(pars, rp), "square" = tpsq(pars), 
    "exp" = tpexp(pars), "exp10sq" = tpexp10sq(pars))

  if (length(idphi) > 0)
  {
    stopifnot(length(idphi) == 2)

##NOTE the ransformed parameters depend on each other, cross-derivatives are not zero

    if (type.phi == "v1") # Kim and Nelson (1999)
    {
      opabsx <- 1 + abs(pars[idphi])
      aux <- pars[idphi] / opabsx
      res$pars[idphi] <- c(sum(aux), -prod(aux))

      if (gradient)
      {
        res$d1["phi1"] <- (opabsx["phi1"] - 
          pars["phi1"] * sign(pars["phi1"])) / opabsx["phi1"]^2
        res$d1["phi2"] <- -aux["phi1"] * (opabsx["phi2"] - 
          pars["phi2"] * sign(pars["phi2"])) / opabsx["phi2"]^2
      }

      if (hessian)
      {
##FIXME TODO
      }
    } else
    if (type.phi == "v2") # Morley, Nelson and Zivot (2003)
    {
      aaa <- pars["phi1"] / (1 + abs(pars["phi1"]))
      ccc <- (1 - abs(aaa))*pars["phi2"] / (1 + abs(pars["phi2"])) + abs(aaa) - aaa^2
      res$pars["phi1"] <- 2*aaa
      res$pars["phi2"] <- -(aaa^2 + ccc)

      if (gradient)
      {
##FIXME TODO
        #res$d1["phi1"] <- 
        #res$d1["phi2"] <- 
      }

      if (hessian)
      {
##FIXME TODO
      }
    } else
    if (type.phi == "stats") # stats package
    {
      #NOTE in this model, AR(2), the periodicity of the data 
      #has no effect here, 1 is passed for convenience
      #.Call(stats:::C_ARIMA_transPars, c(pars["phi1"], pars["phi2"]), 
      #  c(2L,0L,0L,0L,as.integer(frequency(model$y)),0L,0L), TRUE)[[1]]      
      #tmp <- .Call(stats:::C_ARIMA_transPars, c(pars["phi1"], pars["phi2"]), 
      #  c(2L,0L,0L,0L,1L,0L,0L), TRUE)[[1]]
      #res$pars["phi1"] <- tmp[1]
      #res$pars["phi2"] <- tmp[2]
      
      #based on source code of "stats:::C_ARIMA_transPars"
      tanhp <- tanh(c(pars["phi1"], pars["phi2"]))
      res$pars["phi1"] <- tanhp[1] - tanhp[2] * tanhp[1]
      res$pars["phi2"] <- tanhp[2]

      if (gradient)
      {
        #tmp <- .Call(stats:::C_ARIMA_Gradtrans, c(pars["phi1"], pars["phi2"]), 
        #  c(2L,0L,0L,0L,1L,0L,0L))
        #res$d1["phi1"] <- tmp[1,1]
        #res$d1["phi2"] <- tmp[2,2]

        #analytical derivatives of the equations inferred from the 
        #source code of "stats:::C_ARIMA_transPars"
        gp <- matrix(0, nrow=2, ncol=2)
        tmp <- 1/cosh(c(pars["phi1"], pars["phi2"]))^2
        gp[1,1] <- tmp[1] - tanhp[2] * tmp[1]
        gp[2,2] <- tmp[2]
        gp[2,1] <- -tanhp[1] * tmp[2]

        # return a matrix (derivative with respect "phi2" depends on "phi1")
        d1 <- array(0, c(np, np), dimnames=list(pnms, pnms))
        d1[cbind(idvar,idvar)] <- res$d1[idvar]
        d1[idphi,idphi] <- gp
        # location or other possible parameters (currently "a0")
        id <- seq_along(pnms)[-c(idvar,idphi)]
        if (length(id) > 0)
          d1[id,id] <- res$d1[id]
        res$d1 <- d1
      }

      if (hessian)
      {
##FIXME TODO
      }
    } else    
# if (type.phi == "stats")
#this transformation of phi (phi^2) returns the right derivatives of mloglik;
#this shows that the problem arise from the fact that the transformation of phi depends on phi2
#not on having phi1 and phi2 in matrix T, 
# {
#   res$pars["phi1"] <- res$pars["phi1"]^2
#   res$pars["phi2"] <- res$pars["phi2"]^2
#   if (gradient)
#   {
#     res$d1[idphi] <- 2 * pars[idphi]
#   }
#   
# } else
    if (type.phi == "none")
    {
      if (gradient)
        res$d1["phi1"] <- res$d1["phi2"] <- 0
    }
  }

  list(pars = res$pars, gradient = res$d1, hessian = res$d2)
})

setMethod(f = "transPars", signature = "stsm", 
  definition = function(x, type = NULL,
    gradient = FALSE, hessian = FALSE,     
    rp, sclrho = 1.7, sclomega = 1.7, ftrans = NULL, ...) {

  if (is.null(x@transPars))
    return(list(pars = x@pars))

  if (is.function(x@transPars)) {
    type <- NULL
    ftrans <- x@transPars
    rp <- NULL
  } else {
    #type <- match.arg(x@transPars, c("StructTS", "square", "exp", "4", "nk"))
    type <- match.arg(x@transPars, eval(formals(transPars)$type))

    rp <- switch(type, 
      "StructTS" = var(x@y, na.rm = TRUE) / 100,
      NULL) # default
  }
    
  transPars(x@pars, type = type, gradient = gradient, hessian = hessian,
    rp = rp, sclrho = sclrho, sclomega = sclomega, ftrans = ftrans, ...)
})
