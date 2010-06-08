######################################################################
mcmc.qtlnet <- function(cross, pheno.col, threshold,
                        addcov=NULL, intcov=NULL,
                        nSamples = 1000, thinning=1,
                        max.parents = 3,
                        M0 = init.qtlnet(pheno.col, max.parents, n.edge),
                        burnin = 0.1, method = "hk", random.seed = NULL,
                        n.edge = 0,
                        saved.scores = NULL,
                        verbose = FALSE, ...)
{
  ## Check input parameters.

  ## Random number generator seed.
  if(!is.null(random.seed)) {
    if(!is.numeric(random.seed))
      stop("random seed must be numeric")
    set.seed(random.seed)
  }

  ## Burnin must be between 0 and 1.
  if(is.logical(burnin))
    burnin <- ifelse(burnin, 0.1, 0)
  if(burnin < 0 | burnin > 1)
    stop("burnin must be between 0 and 1")

  ## Adjust phenotypes and covariates to be numeric.
  cross <- adjust.pheno(cross, pheno.col, addcov, intcov)
  pheno.col <- cross$pheno.col
  pheno.names <- cross$pheno.names
  addcov <- cross$addcov
  intcov <- cross$intcov
  cross <- cross$cross

  n.pheno <- length(pheno.col)

  ## Initial network matrix.
  if(nrow(M0) != n.pheno | ncol(M0) != n.pheno)
    paste("M0 must be square matrix the size of pheno.col")

  ## LOD threshold by phenotype.
  if(length(threshold) == 1)
    threshold <- rep(threshold, n.pheno)
  if(length(threshold) != n.pheno)
    stop("threshold must have same length as pheno.col")

  ## Calculate genotype probabilities if not already done.
  if (!("prob" %in% names(cross$geno[[1]]))) {
    warning("First running calc.genoprob.")
    cross <- calc.genoprob(cross)
  }
  
  post.net.str <- array(dim=c(n.pheno,n.pheno,nSamples))
  post.bic <- rep(NA,nSamples)
  post.model <- rep(NA,nSamples)
  all.bic <- rep(NA,nSamples)

  saved.scores <- make.saved.scores(pheno.names, max.parents,
                                    saved.scores = saved.scores,
                                    verbose = verbose, ...)
        
  M.old <- M0
  ne.old <- nbhd.size(M.old, max.parents)[[1]]
  aux.new <- score.model(M.old, saved.scores, cross, addcov, intcov,
                         threshold, verbose, method = method, ...)
  if(verbose) cat("\n")
  
  bic.old <- aux.new$model.score
  tmp <- aux.new$update.scores
  codes <- dimnames(saved.scores)[[1]]
  if(!is.null(tmp)) {
    index <- nrow(saved.scores)
    index <- match(tmp[,1], codes) + index * (tmp[,2] - 1)
    saved.scores[index] <- tmp[,3]
  }
  model.old <- aux.new$model.name

  k <- 0
  if(thinning <= 1) {
    post.bic[1] <- all.bic[1] <- bic.old
    post.net.str[,,1] <- M.old
    post.model[1] <- model.old
    k <- 1
  }
  cont.accept <- 0
  for(i in 2:(nSamples*thinning)){
    M.new <- propose.new.structure(M.old, max.parents,
                                   verbose = (verbose > 1))
    ne.new <- nbhd.size(M.new, max.parents)[[1]]
    aux.new <- score.model(M.new, saved.scores, cross, addcov, intcov,
                           threshold, verbose, method = method, ...)

    ## Typically only a few update.scores.
    tmp <- aux.new$update.scores
    if(!is.null(tmp)) {
      index <- nrow(saved.scores)
      tmp2 <- match(tmp[,1], codes)
      if(any(is.na(tmp2)))
         browser()
      index <- match(tmp[,1], codes) + index * (tmp[,2] - 1)
      saved.scores[index] <- tmp[,3]
    }

    bic.new <- aux.new$model.score
    model.new <- aux.new$model.name


    mr <- exp(-0.5*(bic.new - bic.old))*(ne.old/ne.new)
    if(runif(1) <= min(1,mr)){
      M.old <- M.new
      bic.old <- bic.new
      ne.old <- ne.new
      model.old <- model.new
      cont.accept <- cont.accept + 1
    }

    ## Bookkeeping to save sample.
    if((i %% thinning) == 0){      
      k <- k + 1
      all.bic[k] <- bic.new
      post.bic[k] <- bic.old
      post.net.str[,,k] <- M.old
      post.model[k] <- model.old
      if(verbose) print(c(i,k)) 
    }
  }           
  out <- list(post.model=post.model,
              post.bic=post.bic, 
              post.net.str=post.net.str, 
              freq.accept = cont.accept / (nSamples * thinning), 
              saved.scores=saved.scores,
              all.bic=all.bic,
              cross = cross)

  ## Attributes of qtlnet object.
  attr(out, "M0") <- M0
  attr(out, "threshold") <- threshold
  attr(out, "nSamples") <- nSamples
  attr(out, "thinning") <- thinning
  attr(out, "pheno.col") <- pheno.col
  attr(out, "pheno.names") <- pheno.names
  attr(out, "addcov") <- addcov
  attr(out, "intcov") <- intcov
  attr(out, "burnin") <- burnin
  attr(out, "method") <- method
  attr(out, "random.seed") <- random.seed
  attr(out, "random.kind") <- RNGkind()

  ## Add model average calculation.
  out$model.average <- qtlnet.average(out, burnin = burnin)

  class(out) <- c("qtlnet","list")
  
  out
}
######################################################################
c.qtlnet <- function(...)
{
  ## Combine qtlnet objects.
  ## Might be useful to have summary that just pulls out the model averages.
  ## Also need to think carefully about burnin--what does it mean here?
  
  netlist <- list(...)
  out <- netlist[[1]]
  if(!inherits(out, "qtlnet")) {
    netlist <- out
    out <- netlist[[1]]
  }
  if(!inherits(out, "qtlnet"))
    stop("argument must be list of qtlnet objects")
  
  if(length(netlist) > 1) {
    ## Need to do some checking here that qtlnet objects match up.
    ## Minimal for now.

    ## Assumes that parameters stay the same, although nSamples could change.
    
    n.pheno <- length(attr(out, "pheno.col"))
    if(n.pheno != mean(sapply(netlist, function(x) length(attr(x, "pheno.col")))))
      stop("different numbers of phenotypes not allowed")
    if(attr(out, "thinning") != mean(sapply(netlist, function(x) attr(x, "thinning"))))
      stop("thinning values differ")
    ## check threshold vector.
    if(!all(attr(out, "threshold") == apply(sapply(netlist, function(x) attr(x, "threshold")), 1, mean)))
      stop("threshold values differ")

    for(i in seq(2, length(netlist))) {
      ## Per step summaries.
      for(j in c("post.model","post.bic","all.bic"))
        out[[j]] <- c(out[[j]], netlist[[i]][[j]])

      ## Attributes.
      n1 <- sum(attr(out, "nSamples"))
      n2 <- attr(netlist[[i]], "nSamples")
      attr(out, "nSamples") <- c(attr(out, "nSamples"), n2)

      ## Matrices of network structure (3-D array).
      out$post.net.str <- array(c(out$post.net.str, netlist[[i]]$post.net.str),
                                c(n.pheno, n.pheno, n1 + n2))
      ## Acceptance frequency.
      out$freq.accept <- (n1 * out$freq.accept + n2 * netlist[[i]]$freq.accept) / (n1 + n2)
    }
    out$model.average <- qtlnet.average(out, burnin = attr(out, "burnin"))
  }
  out
}
######################################################################
nbhd.size <- function(M, max.parents = 3)
{
  n.deletions <- sum(M)
  n.additions <- 0 
  n.reversions <- 0 
  le <- ncol(M)
  for(j in 1:le){
        add.forbid <- forbidden.additions(M, j, max.parents)
        n.additions <- n.additions + (le - 1 - length(c(add.forbid$upf, add.forbid$downf)))
        rev.allow <- check.reversions(M, j, max.parents)$allowed
        if(!is.null(rev.allow))
          n.reversions <- n.reversions + nrow( rev.allow )
  }
  nbhd.size <- n.deletions + n.additions + n.reversions
  list(nbhd.size=nbhd.size, 
       n.deletions=n.deletions,
       n.additions=n.additions, 
       n.reversions=n.reversions)
}######################################################################
init.qtlnet <- function(pheno.col, max.parents = 3, init.edges = NULL)
{
  n.pheno <- length(pheno.col)
  if(is.null(init.edges)) {
    n.edges <- n.pheno * (n.pheno - 1) / 2
    init.edges <- sample(seq(0, n.edges), 1)
  }
  M <- matrix(0, n.pheno, n.pheno)
  if(init.edges > 0)
    for(i in seq(init.edges)) {
      cause <- sample(n.pheno, 1)
      if(i == 1)
        effect <- sample(seq(n.pheno)[-cause], 1)
      else
        effect <- propose.add(M, cause, max.parents)
      if(effect == 0)
        break
      M[cause, effect] <- 1
    }
  M
}
######################################################################
propose.add <- function(M, node, max.parents)
{
  ## Add an edge with a direction.
  aux1 <- forbidden.additions(M, node, max.parents)
  ## Is there a down option?
  aux3 <- unique(c(node, aux1$upf, aux1$downf))
  aux3 <- seq(ncol(M))[-aux3]
  
  if(length(aux3)) {
    ## Check if any of aux3 is at or above max.parents.
    wh <- which(apply(M[, aux3, drop = FALSE], 2, sum) >= max.parents)
    if(length(wh))
      aux3 <- aux3[-wh]
  }
  if(length(aux3) == 0)
    aux3 <- 0
  else {
    if(length(aux3) > 1)
      aux3 <- sample(aux3, 1)
  }
  aux3
}  
######################################################################
propose.new.structure <- function(M, max.parents = 3, verbose = FALSE)
{
  ## Acceptance rate is 20%. Could we improve here?

  ## To speed up, use M as logical matrix. Still return 0/1 matrix for now.
  le.nodes <- ncol(M)
  flag <- TRUE
  while(flag){
    ## Pick node and decide on add/delete/reverse.
    ## Keep doing this until successful.
    
    node <- sample(seq(le.nodes),1)
    move <- sample(c("add","delete","reverse"),1)
    if(verbose)
      cat(node, move, "")

    switch(move,
           add = {
             aux3 <- propose.add(M, node, max.parents)
             if(!(flag <- (aux3 == 0))) {
               M[node,aux3] <- 1
               if(verbose)
                 cat(aux3, "\n")
             } 
           },
           reverse = {
             ## Reverse direction of an edge.
             aux1 <- check.reversions(M, node, max.parents)$allowed
             if(!(flag <- is.null(aux1))) {
               le.rev <- nrow(aux1)
               aux3 <- sample(seq(le.rev), 1)
               aux3 <- aux1[aux3, ]
               M[aux3[1], aux3[2]] <- 0
               M[aux3[2], aux3[1]] <- 1
               if(verbose)
                 cat(aux3[1], aux3[2], "\n")
             }
           },
           delete = {
             ## Delete an existing edge.
             aux1 <- which(M[, node] == 1)
             if(!(flag <- (length(aux1) == 0))) {
               if(length(aux1) > 1)
                 aux1 <- sample(aux1, 1)
               M[aux1, node] <- 0
               if(verbose)
                 cat(aux1, "\n")
             }
           })
  }
  M
}
######################################################################
forbidden.additions <- function(M, node, max.parents = 3)
{
  ## upf are forbidden upstream nodes (already present downstream).
  ## downf are forbidden downstream nodes (already present upstream).
  ## Check on cycles is one step. See check.reversions() for longer cycles.

  ## Forbidden upstream additions
  upf <- which(M[node,] == 1)
  if(length(upf)==0)
    upf <- NULL

  ## Forbidden downstream additions. More complicated.
  downf <- which(M[,node] == 1)
  le <- length(downf)
  
  ## If at least max.parents are causal for node, forbid any more.
  if(le >= max.parents)
    downf <- seq(nrow(M))[-node]
  else {
    if(le > 0)
      downf <- check.downstream(M, downf)
    else
      downf <- NULL
  }
  
  list(upf=upf, downf=downf)
} 
######################################################################
check.downstream <- function(M, downf)
{
  ## After accounting for scanone, 85% of time is nbhd.size.
  ## Of that, 85% (75% overall) is in check.downstream.

  count.downok <- nrow(M) - 1 - length(downf)
  
  ## Check downstream to see if a cycle would be created.
  tmpfn <- function(x) sum(x)
  is.down <- apply(M[, downf, drop = FALSE], 1, tmpfn)
  flag <- TRUE
  while(flag){
    aux1 <- which(is.down > 0)
    if(flag <- (length(aux1) > 0 & count.downok > 0)) {
      aux1 <- aux1[!(aux1 %in% downf)]
      if(flag <- length(aux1)) {
        is.down <- apply(M[, aux1, drop = FALSE], 1, tmpfn)
        downf <- c(downf, aux1)
        count.downok <- count.downok - flag
      }
    }
  }
  downf
}
######################################################################
check.reversions <- function(M, node, max.parents = 3)
{
  ## Check on possible reversals of directions.
  ## allowed if no cycles produced.
  ## forbidden if cycles result.
  down <- which(M[,node] == 1)
  le <- length(down)

  allowed <- forbidden <- NULL
  
  if(le == 0)
    forbidden <- cbind(seq(nrow(M))[-node], node)
  else {
    if(le == 1)
      allowed <- cbind(down, node)
    else { ## le > 1
      ## Multiple colliders.
      forbid.down <- rep(FALSE, le)
      for(k in 1:le)
        forbid.down[k] <- down[k] %in% check.downstream(M, down[-k])
      if(any(forbid.down)){
        if(any(!forbid.down))
          allowed <- cbind(down[!forbid.down], node)
        forbidden <- cbind(down[forbid.down], node)
      }
      else
        allowed <- cbind(down, node)
    }

    if(!is.null(allowed)) {
      ## Final check of max.parents.
      wh <- which(apply(M[, allowed[,1], drop = FALSE], 2, sum) >= max.parents)
      if(length(wh)) {
        ## Move pairs from allowed to forbidden.
        forbidden <- rbind(forbidden, allowed[wh,])
        allowed <- allowed[-wh,, drop = FALSE]
        if(nrow(allowed) == 0)
          allowed <- NULL
      }
    }
  }
    
  list(allowed = allowed, forbidden = forbidden)   
}
