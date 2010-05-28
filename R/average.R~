qtlnet.average <- function(qtlnet.object, burnin = 0)
{
  n <- length(qtlnet.object[[1]])
  if(burnin > 0)
    b <- floor(burnin * n)
  else
    b <- 1
  nSamples <- n - b + 1
  
  pp <- get.posterior.prob(qtlnet.object, b, n)
  aux.pp <- rep(NA,nSamples)
  le <- nrow(pp)
  np <- ncol(qtlnet.object[[3]][,,1])
  unique.str <- array(dim=c(np,np,le))

  ## Somewhere this is not working.
  mav <- matrix(NA,np,np)
  pheno.names <- attr(qtlnet.object, "pheno.names")
  dimnames(mav) <- list(pheno.names, pheno.names)
  for(i in 1:le){
    aux1 <- pp[i,1]
    aux2 <- (which(qtlnet.object[[1]][b:n] == aux1)[1])+b-1
    unique.str[,,i] <- qtlnet.object[[3]][,,aux2]*pp[i,2]
  } 
  out <- data.frame(matrix(NA,(np^2)-np, 2))
  names(out) <- c("arc","Posterior")
  k <- 1
  for(i in 1:np){
    for(j in 1:np){
      mav[i,j] <- sum(unique.str[i,j,])
      if(i != j){
        out[k,1] <- paste(i, paste("-->",j,sep=" "), sep=" ")
        out[k,2] <- mav[i,j]
        k <- k + 1
      }
    }
  }
  aux.o <- order(out[,2],decreasing=TRUE)
  out <- out[aux.o,]
  row.names(out) <- c(1:((np^2)-np))  
  list(mav=mav,out=out,pp=pp)
}
######################################################################
get.posterior.prob <- function(qtlnet.object, b, n)
{
  nSamples <- n - b + 1

  sampled.models <- unique(qtlnet.object[[1]][b:n])
  le <- length(sampled.models)
  post.prob <- data.frame(matrix(NA,le,3))
  names(post.prob) <- c("Model","Posterior prob","BIC")
  post.prob[,1] <- sampled.models
  for(i in 1:le){
    wh <- which(qtlnet.object$post.model[b:n] == sampled.models[i])
    post.prob[i,2] <- length(wh)/nSamples
    post.prob[i,3] <- mean(qtlnet.object$post.bic[b:n][wh])
  }
  aux <- order(post.prob[,2],decreasing=TRUE)
  post.prob <- post.prob[aux,]
  row.names(post.prob) <- c(1:le)
  post.prob
}
