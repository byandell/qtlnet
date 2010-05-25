print.qtlnet <- function(x, cutoff = 0.01, digits = 3, ...)
{
  cat("\nModel averaged probabilities for edge direction (row -> col):\n")
  print(round(x$model.average$mav, digits))

  cat("\nPosterior probabilities by causal model:\n")
  wh <- which(x$model.average$pp[,2] >= cutoff)
  if(!length(wh))
    wh <- which.max(x$model.average$pp[,2])
  print(x$model.average$pp[wh,, drop = FALSE])
  
  invisible(x$model.average)
}  
######################################################################
summary.qtlnet <- function(object, ...)
{
  mav <- object$model.average[[1]]
  node.nms <- attr(object, "pheno.names")

  n <- nrow(mav)
  new <- mav
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      m1 <- mav[i,j]
      m2 <- mav[j,i]
      m3 <- 1 - m1 - m2
      aux <- which.max(c(m1,m2,m3))
      if(aux == 1) 
        new[j,i] <- 0
      if(aux == 2) 
        new[i,j] <- 0
      if(aux == 3){ 
        new[j,i] <- 0
        new[i,j] <- 0
      }
    }
  }  
  out <- list(freq.accept = object$freq.accept,
              averaged.net = M.2.lista(new, node.nms),
              posterior.table = averaged.posterior.table(mav, node.nms),
              parent.patterns = parent.qtlnet(object))
  class(out) <- c("summary.qtlnet", "list")
  out
}
######################################################################
print.summary.qtlnet <- function(x, ...)
{
  cat("\nModel-averaged network:\n")
  print(x$averaged.net)

  cat("\nPosterior probabilities by direction:\n")
  print(x$posterior.table)

  cat("\nAcceptance frequency for MCMC:", x$freq.accept, "\n")

  print(x$parent.patterns, ...)
  
  invisible(x)
}
######################################################################
parent.qtlnet <- function(qtlnet.object)
{
  ## Split model into trait|depends), with "|" and "depends" optional.

  ## Model patterns.
  model <- qtlnet.object$post.model
  model <- unlist(strsplit(model, "(", fixed = TRUE))
  model <- model[model != ""]
  
  ## Parent patterns.
  parents <- sapply(strsplit(model, "|", fixed = TRUE), function(x) x[2])

  ## Size of parent patterns.
  n.parents <- sapply(strsplit(parents, ",", fixed = TRUE),
                      function(x) ifelse(is.na(x[1]), 0, length(x)))

  ## How many different traits share the same parent set?
  tmp <- sapply(tapply(model, parents, table), length)
  tmp2 <- sapply(strsplit(names(tmp), ",", fixed = TRUE), length)
  scans <- sapply(tapply(tmp, tmp2,
                       function(x) c(solo = sum(x == 1), more = sum(x > 1))),
                function(x) x)
  scans <- rbind(scans, total = apply(scans, 2, sum))
  scans <- cbind(scans, total = apply(scans, 1, sum))
  
  out <- list(model = model, parents = parents, n.parents = n.parents,
              scans = scans)
  class(out) <- c("parent.qtlnet", "list")
  out
}
######################################################################
print.parent.qtlnet <- function(x, freq.max = 10, ...)
{
  mytable <- function(parents, freq.max = 10) {
    tmp <- table(parents)
    tmp2 <- "rest" ## paste(">=", freq.max, sep = "")
    tmp[tmp >= freq.max] <- tmp2
    tmp <- ordered(tmp, c(seq(freq.max - 1), tmp2))
    tmp <- c(table(tmp))
    tmp["total"] <- sum(tmp)
    tmp
  }
  
  cat("\nHow frequently are model patterns sampled?\n")
  print(mytable(x$model, freq.max = freq.max))
  
  cat("\nHow frequently are parent patterns sampled?\n")
  print(mytable(x$parents, freq.max = freq.max))
  
  cat("\nWhat are sizes of parent pattern sets?\n")
  tmp <- table(x$n.parents)
  tmp["total"] <- sum(tmp)
  print(tmp)
  
  cat("\nHow frequently are parent sets sampled (by set size)?\n")
  parents.size <- tapply(x$parents[x$n.parents > 0],
                         x$n.parents[x$n.parents > 0],
                         mytable, freq.max = freq.max)
  print(parents.size)

  cat("\nWhat are frequencies of each single-parent pattern?\n")
  tmp <- table(x$parents[x$n.parents == 1])
  names(tmp) <- substring(names(tmp), 1, nchar(names(tmp)) - 1)
  tmp <- tmp[order(as.numeric(names(tmp)))]
  print(tmp)
  
  cat("\nHow many scans are solo or more for each size parent pattern?\n")
  print(x$scans)
  
  ## Evidence suggest it is not worth doing multiple traits
  ## for n > le.pheno / 2.
}
######################################################################
######################################################################
averaged.posterior.table <- function(maM,nms)
{
  nn <- nrow(maM)
  np <- choose(nn,2)
  out <- data.frame(matrix(NA,np,5))
  k <- 1
  for(i in 1:(nn-1)){
    for(j in (i+1):nn){
      out[k,1] <- i
      out[k,2] <- j
      out[k,3] <- round(maM[i,j],3)
      out[k,4] <- round(maM[j,i],3)
      out[k,5] <- round(1 - out[k,3] - out[k,4], 3)
      k <- k + 1
    }
  }
  new.out <- out
  for(k in 1:np){
    new.out[k,1] <- nms[out[k,1]]
    new.out[k,2] <- nms[out[k,2]]
  }
  names(new.out) <- c("node1","node2","-->","<--","no")
  new.out
}
################################################################
M.2.lista <- function(M,node.nms)
{
  le <- nrow(M)
  k <- 1
  le2 <- le*(le-1)/2
  out <- data.frame(cause = factor(rep(NA, le2), node.nms),
                    effect = factor(rep(NA, le2), node.nms),
                    prob = rep(NA, le2))

  for(i in 1:(le-1)){
    for(j in (i+1):le){
      if(M[i,j] != 0){
        out[k,1] <- node.nms[i]
##        out[k,2] <- "---->"
        out[k,2] <- node.nms[j]
        out[k,3] <- M[i,j]
        k <- k + 1
      }
      if(M[j,i] != 0){
        out[k,2] <- node.nms[i]
        out[k,1] <- node.nms[j]
##        out[k,1] <- node.nms[i]
##        out[k,2] <- "<----"
##        out[k,3] <- node.nms[j]
        out[k,3] <- M[j,i]
        k <- k + 1
      }
    }
  }
  out[1:(k-1),]
}
