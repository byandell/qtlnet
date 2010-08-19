get.model.average <- function(qtlnet.object)
{
  mav <- qtlnet.object$Mav
  if(is.null(mav)) {
    post.burnin <- get.post.burnin(qtlnet.object)

    mav <- apply(qtlnet.object$post.net.str[,,post.burnin], c(1,2), mean)
    pheno.names <- attr(qtlnet.object, "pheno.names")
    dimnames(mav) <- list(pheno.names, pheno.names)
  }
  mav
}
##########################################################################
get.post.burnin <- function(qtlnet.object)
{
  burnin <- attr(qtlnet.object, "burnin")
  nSamples <- attr(qtlnet.object, "nSamples")
  which(apply(matrix(nSamples), 1,
              function(x,b) seq(x) >= b * x, burnin))
}
##########################################################################
get.posterior.prob <- function(qtlnet.object)
{
  post.burnin <- get.post.burnin(qtlnet.object)
  post.model <- qtlnet.object$post.model[post.burnin]
  post.bic <- qtlnet.object$post.bic[post.burnin]
  
  out <- data.frame(post.prob = as.vector(table(post.model)) / length(post.model),
               BIC  = tapply(post.bic, post.model, mean))
  out[order(-out$post.prob, out$BIC),]
}
##########################################################################
## Need to condense qtlnet object.
## First, get rid of post.net.str.
## Could have utility to translate post.model to post.net.str.
##########################################################################
model2M <- function(post.model)
{
  ## Convert post.model into 3-D M array.
  
  ## Strip out parentheses and split by node.
  a <- strsplit(substring(post.model, 2, nchar(post.model) - 1), ")(", fixed = TRUE)
  nSamples <- length(a)
  n.pheno <- length(a[[1]])
  phenos <- rep(0, n.pheno)
  asplit <- function(a, phenos = length(a)) {
    a1 <- strsplit(a, "|", fixed = TRUE)
    a2 <- sapply(a1, function(x, phenos)
                 {
                   if(length(x) == 2) {
                     xx <- as.numeric(strsplit(x[2], ",", fixed = TRUE)[[1]])
                     phenos[xx] <- 1
                   }
                   phenos
                 }, phenos)
    a2
  }
  array(unlist(lapply(a, asplit, phenos)), c(n.pheno, n.pheno, nSamples))
}
##########################################################################
subset.qtlnet <- function(x, run, ...)
{
  if(missing(run))
    stop("run must be specified")

  nSamples <- attr(x, "nSamples")
  runs <- length(nSamples)
  if(any(run < 0) | any(run > runs))
    stop(paste("run must be integer between 1 and", runs))
  
  run.id <- rep(seq(runs), nSamples)

  out <- x
  wh <- which(run.id %in% run)
  for(i in c("post.model","post.bic","all.bic"))
    out[[i]] <- out[[i]][wh]
  M <- model2M(out$post.model)
  out$Msd <- NULL
  attr(out, "M0") <- M[,,1]
  attr(out, "nSamples") <- attr(out, "nSamples")[run]

  ## These are approximate as only using saved MCMC samples.
  out$Mav <- apply(M, 1:2, mean)
  out$freq.accept <- mean(out$all.bic == out$post.bic)
  
  out
}
