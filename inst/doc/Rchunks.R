###################################################
### .Rprofile
###################################################
library("stsm.class", quietly = TRUE)
library("stats", quietly = TRUE)
library("stsm", quietly = TRUE)
library("KFKSDS", quietly = TRUE)

###################################################
### chunk 'stsm-class-1'
###################################################
library("stsm.class")
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson)

###################################################
### chunk 'stsm-class-2'
###################################################
m <- set.pars(m, c("var1" = 2, "var2" = 15, "var3" = 30))
get.pars(m)

###################################################
### chunk 'stsm-class-3'
###################################################
ss <- char2numeric(m)
ss$T
ss$V

###################################################
### chunk 'stsm-class-4'
###################################################
m <- stsm.model(model = "llm+seas", y = JohnsonJohnson, 
  pars = c(var1 = 2, var2 = 3, var3 = 4), transPars = "square")
m@pars

###################################################
### chunk 'stsm-class-5'
###################################################
get.pars(m)

###################################################
### chunk 'airp-stsm-ml-td-1'
###################################################
library("stsm")
res0 <- StructTS(log(AirPassengers), type = "BSM")$coef[c(4,1:3)]
mairp <- stsm.model(model = "BSM", y = log(AirPassengers), 
  transPars = "StructTS")
res1 <- maxlik.td.optim(mairp, KF.version = "KFKSDS", 
  KF.args = list(P0cov = TRUE), method = "L-BFGS-B", gr = "numerical")
mairp <- set.pars(mairp, pmax(res1$par, .Machine$double.eps))
round(get.pars(mairp), 6)
all.equal(get.pars(mairp), res0, tol = 1e-05, check.att = FALSE)

###################################################
### chunk 'airp-stsm-ml-td-2'
###################################################
mairp <- set.pars(mairp, c(1,1,1,1))
res2 <- maxlik.td.optim(mairp, KF.version = "KFKSDS", 
  KF.args = list(P0cov = FALSE), method = "L-BFGS-B", gr = "numerical")
mairp <- set.pars(mairp, pmax(res2$par, .Machine$double.eps))
round(get.pars(mairp), 6)

###################################################
### chunk 'airp-stsm-ml-td-3'
###################################################
mairp <- stsm.model(model = "BSM", y = log(AirPassengers), 
  transPars = "square")
res3 <- maxlik.td.optim(mairp, KF.version = "KFKSDS", 
  KF.args = list(P0cov = FALSE), method = "BFGS", gr = "numerical")
round(res3$par^2, 6)

###################################################
### chunk 'airp-stsm-ml-td-4'
###################################################
mairp <- stsm.model(model = "BSM", y = log(AirPassengers), 
  transPars = "exp")
res4 <- maxlik.td.optim(mairp, KF.version = "KFKSDS", 
  KF.args = list(P0cov = FALSE), method = "BFGS", gr = "numerical")
round(exp(res4$par), 6)

###################################################
### chunk 'airp-stsm-ml-td-5'
###################################################
mairp <- stsm.model(model = "BSM", y = log(AirPassengers), 
  transPars = "StructTS", nopars = c(var3 = 0))
res5 <- maxlik.td.optim(mairp, KF.version = "KFKSDS", 
  KF.args = list(P0cov = TRUE), method = "L-BFGS-B", gr = "numerical")
mairp <- set.pars(mairp, res5$par)
round(get.pars(mairp), 6)

###################################################
### chunk 'airp-stsm-ml-td-6'
###################################################
mairp <- stsm.model(model = "llm+seas", y = log(AirPassengers), 
  transPars = "StructTS")
res6 <- maxlik.td.optim(mairp, KF.version = "KFKSDS",
  KF.args = list(P0cov = FALSE), method = "L-BFGS-B", gr = "numerical")
mairp <- set.pars(mairp, res6$par)
round(get.pars(mairp), 6)

###################################################
### chunk 'airp-stsm-ml-td-7'
###################################################
mairp <- stsm.model(model = "BSM", y = log(AirPassengers), 
  nopars = c(var3 = 0), cpar = c(var2 = 1), transPars = "StructTS")
res7 <- maxlik.td.optim(mairp, KF.version = "KFKSDS",
  KF.args = list(P0cov = FALSE), method = "L-BFGS-B", gr = "numerical")
mairp <- res7$model
round(coef(res7), 6)

###################################################
### chunk 'fig-draft-airp-ml-td'
###################################################
library("KFKSDS")
plot(tsSmooth(StructTS(mairp@y, type = "BSM"))[,1:3])
ss <- char2numeric(mairp)
kf <- KF(mairp@y, ss)
ks <- KS(mairp@y, ss, kf)
plot(ks$ahat[,c(1,2,3)])

###################################################
### chunk 'fig-airp-ml-td'
###################################################
mairp1 <- stsm.model(model = "BSM", y = log(AirPassengers), 
  transPars = "StructTS")
mairp1 <- set.pars(mairp1, pmax(res1$par, .Machine$double.eps))
ss <- char2numeric(mairp1, P0cov = TRUE)
kf1 <- KF(mairp1@y, ss)
ks1 <- KS(mairp1@y, ss, kf1)

mairp2 <- stsm.model(model = "BSM", y = log(AirPassengers), 
  nopars = c(var3 = 0), cpar = c(var2 = 1), transPars = "StructTS")
mairp2 <- set.pars(mairp2, res7$pars[1:2])
mairp2 <- set.cpar(mairp2, tail(res7$pars, 1))
ss <- char2numeric(mairp2, P0cov = FALSE)
kf2 <- KF(mairp2@y, ss)
ks2 <- KS(mairp2@y, ss, kf2)
lks <- list(ks1$ahat[,c(1,3)], ks2$ahat[,c(1,3)])

lairp <- log(AirPassengers)
nf <- 4
titles <- c("Based on results stored in `res1'", 
  "Based on results stored in `res7'")

postscript(file = file.path("figures", "fig-airp-ml-td.eps"), horizontal = FALSE, 
  family = "Times", paper = "special", width = 5.9, height = 2.5)
par(mfcol = c(1, 2), mar = c(1, 1.1, 1, 1.1), las = 1, cex.main = 0.75)
for (i in c(1, 2))
{
  aux <- cbind(lairp, lks[[i]])
  tmp <- apply(aux, 2, scale)
  tmp[,1] <- tmp[,1] + nf
  tmp[,2] <- tmp[,2] + nf
  tmp <- ts(tmp)
  attributes(tmp) <- attributes(aux)

  plot(tmp, plot.type = "single", type = "n", 
    xaxt = "n", yaxt = "n", bty = "n", ylab = "", xlab = "")
  mtext(titles[i], side = 3, adj = 0, cex = 0.75)

  ax <- seq(1948, 1962, 2)
  abline(v = ax, lty = 3, col = "gray75")
  axis(side = 1, at = ax, labels = FALSE, tcl = 0.5, lwd = 0, lwd.tick = 1)
  axis(side = 1, at = ax + 1, labels = FALSE, tcl = 0.25, lwd = 0, lwd.tick = 1)
  axis(side = 1, at = ax, labels = ax, lwd = 0, lwd.tick = 0, line = -1.1, cex.axis = 0.7)

  ay1 <- seq(5, 6.5, 0.5)
  ay1labs <- c("5.0", "5.5", "6.0", "6.5")
  ay2 <- (ay1 - mean(lairp)) / sd(lairp) + nf
  abline(h = ay2, lty = 3, col = "gray75")
  axis(side = 2, at = ay2, labels = FALSE, tcl = 0.25, lwd = 0, lwd.tick = 1)
  axis(side = 2, at = ay2, labels = ay1labs, lwd = 0, lwd.tick = 0, line = -0.8, cex.axis = 0.7)

  ay1 <- seq(-2, 2, 1)
  ay2 <- (ay1 - mean(tmp[,3])) / sd(tmp[,3])
  abline(h = ay2, lty = 3, col = "gray75")
  axis(side = 4, at = ay2, labels = FALSE, tcl = 0.25, lwd = 0, lwd.tick = 1)
  axis(side = 4, at = ay2, labels = ay1, lwd = 0, lwd.tick = 0, line = -0.8, cex.axis = 0.7)

  lines(tmp[,1], col = "gray70")
  lines(tmp[,2])
  lines(tmp[,3])
  box()
}
invisible(dev.off())

###################################################
### chunk 'airp-stsm-ml-td-8'
###################################################
a0 <- c(mairp@y[1], rep(0, 12))
names(a0) <- paste("a0", 1:13, sep = "")
mairp <- stsm.model(model = "BSM", y = log(AirPassengers), pars = a0,
  nopars = c("var1" = 0, "var2" = 0.000772, "var3" = 0, "var4" = 0.001397))
res8 <- maxlik.td.optim(mairp, KF.version = "KFKSDS",
  KF.args = list(P0cov = TRUE), 
  method = "L-BFGS-B", gr = "numerical", hessian = FALSE)
round(res8$par, 6)

###################################################
### chunk 'mlfd-scoring-1'
###################################################
data("llmseas")
m <- stsm.model(model = "llm+seas", y = llmseas[,1], 
  pars = c(var1 = 1, var2 = 1, var3 = 1), 
  ssd = TRUE, sgfc = TRUE)
sgf.d1 <- m@sgfc

###################################################
### chunk 'mlfd-scoring-2'
###################################################
convergence <- FALSE
tol <- 0.001
maxit <- 100
step <- 1
iter <- 0
while (!(convergence || iter > maxit))
{
  sgf <- drop(m@sgfc %*% m@pars)
  pipgog <- pi * as.vector(m@ssd / sgf)
  gd1og <- as.matrix(sgf.d1 / sgf)
  G <- 0.5 * colSums((2 * pipgog - 1) * gd1og)
  IM <- 0
  for (i in seq(along = sgf))
    IM <- IM + 0.5 * tcrossprod(sgf.d1[i,]) / sgf[i]^2
  pd <- drop(solve(IM) %*% G)
  pars.old <- m@pars
  pars.new <- pars.old + step * pd
  m <- set.pars(m, pars.new)
  if (sqrt(sum((pars.old - pars.new)^2)) < tol)
    convergence <- TRUE
  iter <- iter + 1
}

###################################################
### chunk 'mlfd-scoring-airp'
###################################################
mairp <- stsm.model(model = "BSM", y = log(AirPassengers))
res8 <- maxlik.fd.scoring(m = mairp, step = NULL, information = "expected",
  control = list(maxit = 50, tol = 0.001))
round(res8$pars, 6)

###################################################
### chunk 'confint-td-1'
###################################################
ses <- sqrt(diag(vcov(res8, type = "infomat")))
pmax(0, res8$pars - 1.96 * ses)
res8$pars + 1.96 * ses

###################################################
### chunk 'confint-td-2'
###################################################
civcov <- confint(res8, type = "vcov", vcov.type = "infomat", 
  level = 0.95)
round(civcov, 6)

###################################################
### chunk 'confint-fd'
###################################################
set.seed(123)
ciboot <- confint(res8, type = "bootstrap", breps = 100)
round(ciboot, 6)

###################################################
### chunk 'em-1'
###################################################
m <- stsm.model(model = "llm+seas", y = llmseas[,1])
resem1 <- maxlik.em(model = m, type = "standard", 
  tol = 0.001, maxiter = 250)
round(resem1$pars, 6)
paste("number of iterations:", resem1$iter)

###################################################
### chunk 'em-2'
###################################################
resem2 <- maxlik.em(model = m, type = "modified", 
  tol = 0.001, maxiter = 250)
round(resem2$pars, 6)
paste("number of iterations:", resem2$iter)

###################################################
### chunk 'em-3'
###################################################
resem3 <- maxlik.em(model = m, type = "mix", 
  tol = 0.001, maxiter = 250, mod.steps = seq(3, 250, 5))
round(resem3$pars, 6)
paste("number of iterations:", resem3$iter)

###################################################
### chunk 'em-4'
###################################################
mairp <- stsm.model(model = "llm+seas", y = log(AirPassengers))
resem.airp <- maxlik.em(mairp, type = "mix", tol = 0.001, maxiter = 250, 
  mod.steps = seq(3, 250, 5))
round(resem.airp$pars, 6)
paste("number of iterations:", resem.airp$iter)

###################################################
### chunk 'clark-1'
###################################################
data("gdp4795")
gdp <- log(gdp4795)
pars <- c("var2" = 2, "var3" = 2, "var4" = 2, "phi1" = 2, "phi2" = 1)
nopars <- c("var1" = 0, "a01" = 0, "a02" = 0, "a03" = 0, "a04" = 0,
  "P01" = 100, "P02" = 100, "P03" = 100, "P04" = 100)
mgdp <- stsm.model(model = "trend+ar2", y = gdp, 
  pars = pars, nopars = nopars, transPars = "exp10sq")

###################################################
### chunk 'clark-2'
###################################################
res.gdp <- maxlik.td.optim(mgdp, KF.version = "KFKSDS", 
  KF.args = list(P0cov = FALSE, t0 = 21), method = "BFGS")
mgdp <- set.pars(mgdp, res.gdp$par)
optpars <- c(sqrt(get.pars(mgdp)[1:3]), get.pars(mgdp)[4:5])
names(optpars)[1:3] <- paste("stdv", 2:4, sep = "")
round(optpars, 6)

###################################################
### chunk 'clark-3'
###################################################
ss <- char2numeric(mgdp)
kf.mgdp <- KF(mgdp@y, ss)
plot(window(kf.mgdp$a.upd[,3], start = time(mgdp@y)[21]))

###################################################
### chunk 'clark-4'
###################################################
mgdp2 <- stsm.model(model = "local-trend", y = gdp, transPars = "square")
res2.gdp <- maxlik.td.optim(mgdp2, KF.version = "KFKSDS", 
  KF.args = list(P0cov = FALSE, t0 = 1), method = "BFGS")
mgdp2 <- set.pars(mgdp2, get.pars(res2.gdp$model))

###################################################
### chunk 'fig-clark-cycle'
###################################################
postscript(file = file.path("figures", "fig-clark-cycle.eps"), horizontal = FALSE, 
  family = "Times", paper = "special", width = 5.9, height = 2.2)
par(mfcol = c(1, 2), mar = c(1, 1.5, 1, 0.8), las = 1, cex.main = 0.75)
#
# left hand side plot: filtered cycle in the 'trend+ar2' model
#
ss <- char2numeric(mgdp)
kf <- KF(mgdp@y, ss)
wgdp <- window(kf$a.upd[,3], start = time(mgdp@y)[21])
plot(wgdp, xlab = "", ylab = "", 
  #ylim = c(-0.06, 0.04), type = "n",
  ylim = c(-0.042, 0.03), type = "n",  
  xaxt = "n", yaxt = "n", bty = "n")
  #main = "Trend-cycle model"
mtext("Trend-cycle model", side = 3, adj = 0, cex = 0.75)

#turning points dated by the NBER
#Source: http://www.nber.org/cycles/
peaks <- c(1953, 1957, 1960, 1969, 1973, 1980, 1981, 1990) + c(2, 3, 2, 4, 4, 1, 3, 3)/4
troughs <- c(1954, 1958, 1961, 1970, 1975, 1980, 1982, 1991) + c(2, 2, 1, 4, 1, 3, 4, 1)/4
for(j in seq(along = troughs))
{
  rect(troughs[j], -0.1, peaks[j],0.1, border=FALSE, col="gray90")
}

lines(wgdp)
ax <- seq(1945, 1995, 5)
axis(side = 1, at = ax, labels = FALSE, tcl = 0.25, lwd = 0, lwd.tick = 1)
axis(side = 1, at = ax, labels = ax, lwd = 0, lwd.tick = 0, line = -1.1, cex.axis = 0.7)
ay <- seq(-0.06, 0.04, 0.02)
axis(side = 2, at = ay, labels = FALSE, tcl = 0.25, lwd = 0, lwd.tick = 1)
axis(side = 2, at = ay, labels = ay, lwd = 0, lwd.tick = 0, line = -0.8, cex.axis = 0.7)

box()
#
# right hand side plot: smoothed slope in the local trend model
#
mgdp2 <- stsm.model(model = "local-trend", y = gdp, transPars = "square")
res2.gdp <- maxlik.td.optim(mgdp2, KF.version = "KFKSDS", 
  KF.args = list(P0cov = FALSE, t0 = 1), method = "BFGS")
mgdp2 <- set.pars(mgdp2, get.pars(res2.gdp$model))
ss2 <- char2numeric(mgdp2)
kf2 <- KF(mgdp2@y, ss2)
ks2 <- KS(mgdp2@y, ss2, kf2)
wgdp <- window(ks2$ahat[,2], start = time(mgdp@y)[21])

plot(wgdp, xlab = "", ylab = "", ylim = c(0.0025, 0.0125), type = "n", 
  xaxt = "n", yaxt = "n", bty = "n")
mtext("Local trend model", side = 3, adj = 0, cex = 0.75)

for(j in seq(along = troughs))
{
  rect(troughs[j], -0.1, peaks[j],0.1, border=FALSE, col="gray90")
}

lines(wgdp)
ax <- seq(1945, 1995, 5)
axis(side = 1, at = ax, labels = FALSE, tcl = 0.25, lwd = 0, lwd.tick = 1)
axis(side = 1, at = ax, labels = ax, lwd = 0, lwd.tick = 0, line = -1.1, cex.axis = 0.7)
ay <- seq(0, 0.02, len = 6)
axis(side = 2, at = ay, labels = FALSE, tcl = 0.25, lwd = 0, lwd.tick = 1)
axis(side = 2, at = ay, labels = ay, lwd = 0, lwd.tick = 0, line = -0.8, cex.axis = 0.7)
box()
invisible(dev.off())
