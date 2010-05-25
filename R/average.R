qtlnet.average <- function(metroObject, burnin = 0)
{
  n <- length(metroObject[[1]])
  if(burnin > 0)
    b <- floor(burnin * n)
  else
    b <- 1
  nSamples <- n - b + 1
  
  pp <- get.posterior.prob(metroObject, b, n)
  aux.pp <- rep(NA,nSamples)
  le <- nrow(pp)
  np <- ncol(metroObject[[3]][,,1])
  unique.str <- array(dim=c(np,np,le))
  mav <- matrix(NA,np,np)
  pheno.names <- attr(metroObject, "pheno.names")
  dimnames(mav) <- list(pheno.names, pheno.names)
  for(i in 1:le){
    aux1 <- pp[i,1]
    aux2 <- (which(metroObject[[1]][b:n] == aux1)[1])+b-1
    unique.str[,,i] <- metroObject[[3]][,,aux2]*pp[i,2]
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
get.posterior.prob <- function(metroObject, b, n)
{
  nSamples <- n - b + 1

  sampled.models <- unique(metroObject[[1]][b:n])
  le <- length(sampled.models)
  post.prob <- data.frame(matrix(NA,le,2))
  names(post.prob) <- c("Model","Posterior prob")
  post.prob[,1] <- sampled.models
  for(i in 1:le){
    post.prob[i,2] <- length(which(metroObject[[1]][b:n] == sampled.models[i]))/nSamples
  }
  aux <- order(post.prob[,2],decreasing=TRUE)
  post.prob <- post.prob[aux,]
  row.names(post.prob) <- c(1:le)
  post.prob
}
