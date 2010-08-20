disperse.M <- function(qtlnet.object, k = 1)
{
  Mav <- qtlnet.object$Mav
  Msd <- qtlnet.object$Msd
  runs <- length(attr(qtlnet.object, "nSamples"))

  ravel <- row(as.matrix(Mav)) > col(as.matrix(Mav))
  Mavr <- Mav[ravel] - t(Mav)[ravel]
  
  M1 <- Mavr * runs

  ## This captures some but not all of it.
  s1 <- sqrt((runs - 1) ^ 2 + 1) / (runs ^ 2 * sqrt(runs - 1))
  sk <- sqrt(k * (runs * (runs + k) - 2 * k) / (runs - 1)) / runs
  SD <- abs(M1) * s1

  M2 <- sign(M1) * Msd / s1
  data.frame(Mavr = Mavr, Msd = Msd, M1 = M1, M2 = M2, SD = SD)
}
zero.M <- function(qtlnet.object, run = which.min(mbic),
                   burnin = attr(qtlnet.object, "burnin"))
{
  nSamples <- attr(qtlnet.object, "nSamples")
  runs <- length(nSamples)
  run.id <- rep(seq(runs), nSamples)
  M0 <- attr(qtlnet.object, "M0")
  ravel <- row(as.matrix(M0)) > col(as.matrix(M0))

  tmpfn <- function(post.model, burnin, ravel) {
    M <- apply(model2M(post.model), 1:2, sum)
    (M[ravel] + t(M)[ravel]) == 0
  }
  cat("Extracting network matrices...\n")
  out <- matrix(unlist(tapply(qtlnet.object$post.model, run.id, tmpfn,
                              burnin, ravel)),
                ncol = runs)

  mbic <- meanbic(qtlnet.object, burnin)
  wh <- which.min(mbic)

  data.frame(nonzero = apply(out, 2, sum),
             agree = apply(out, 2, function(x,y) sum(x == y & y > 0),
               out[,run]),
             mean.bic = mbic)
}
best.qtlnet <- function(x, burnin = attr(x, "burnin"))
{
  nSamples <- attr(x, "nSamples")
  runs <- length(nSamples)
  if(runs == 1)
    stop("Need more than one run to find best.")
  
  run.id <- rep(seq(runs), nSamples)

  tmpfn <- function(x, burnin) {
    tmp <- which(seq(x) >= burnin * length(x))
    range(x[tmp])
  }
  range.bic <- range(unlist(tapply(x$post.bic, run.id,
                                   tmpfn, burnin)))
  
  plot(c(1,max(nSamples)),range.bic, type = "n",
       xlab = "Sample Index", ylab = "BIC")
  title(paste("BIC samples for", runs, "MCMC runs"))

  tmpfn <- function(x, burnin) {
    tmp <- which(seq(x) >= burnin * length(x))
    lines(tmp, x[tmp])
  }
  tapply(x$post.bic, run.id, tmpfn, burnin)

  mbic <- meanbic(x, burnin)
  
  wh <- which.min(mbic)
  cat("Best sample is", wh, "\n")

  subset(x, wh)
}
meanbic <- function(qtlnet.object, burnin = attr(qtlnet.object, "burnin"))
{
  nSamples <- attr(qtlnet.object, "nSamples")
  runs <- length(nSamples)
  run.id <- rep(seq(runs), nSamples)
  
  tmpfn <- function(x, burnin) {
    tmp <- which(seq(x) >= burnin * length(x))
    mean(x[tmp])
  }
  mbic <- tapply(qtlnet.object$post.bic, run.id, tmpfn, burnin)
}


  
