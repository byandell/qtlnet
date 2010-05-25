######################################################################
mcmc.qtlnet <- function(cross, pheno.col, threshold,
                        addcov=NULL, intcov=NULL,
                        nSamples = 1000, thinning=1,
                        M0 = matrix(0, le.pheno, le.pheno),
                        burnin = 0.1, method = "hk", random.seed = NULL,
                        verbose = FALSE)
{
  le.pheno <- length(pheno.col)

  ## Check input parameters.

  ## LOD threshold by phenotype.
  if(length(threshold) == 1)
    threshold <- rep(threshold, le.pheno)
  threshold <- as.list(threshold)
  if(length(threshold) != le.pheno)
    stop("threshold must have same length as pheno.col")

  ## Random number generator seed.
  if(!is.null(random.seed)) {
    if(!is.numeric(random.seed))
      stop("random seed must be numeric")
    set.seed(random.seed)
  }

  ## Burnin must be between 0 and 1.
  if(is.logical(burnin)) {
    if(burnin)
      burnin <- 0.1
    else
      burnin <- 0
  }
  if(burnin < 0 | burnin > 1)
    stop("burnin must be between 0 and 1")

  ## Initial network matrix.
  if(nrow(M0) != le.pheno | ncol(M0) != le.pheno)
    paste("M0 must be square matrix the size of pheno.col")

  ## Adjust phenotype and covariate names and columns.
  ## Make sure addcov and intcov are NULL or lists.
  pheno.names <- find.pheno.names(cross, pheno.col)
  make.list <- function(x, namex, len) {
    out <- x
    if(!is.null(x)) {
      if(is.list(x)) {
        if(length(x) != len)
          stop(paste(namex, "is list but not of same length as pheno.col"))
      }
      else {
        out <- vector(mode="list", length = len)
        for(i in seq(len))
          out[[i]] <- x
      }
    }
    out
  }
  addcov.names <- make.list(find.pheno.names(cross, addcov))
  intcov.names <- make.list(find.pheno.names(cross, intcov))
  cross$pheno <-
    cross$pheno[, unique(c(pheno.names, unlist(addcov.names), unlist(intcov.names)))]
  pheno.col <- find.pheno(cross, pheno.names)
  if(!is.null(addcov))
    for(i in seq(length(addcov)))
        addcov[[i]] <- find.pheno(cross, addcov.names[[i]])
  if(!is.null(intcov))
    for(i in seq(length(intcov)))
        intcov[[i]] <- find.pheno(cross, intcov.names[[i]])


  ## Calculate genotype probabilities if not already done.
  if (!("prob" %in% names(cross$geno[[1]]))) {
    warning("First running calc.genoprob.")
    cross <- calc.genoprob(cross)
  }
  
  post.net.str <- array(dim=c(le.pheno,le.pheno,nSamples))
  post.bic <- rep(NA,nSamples)
  post.model <- rep(NA,nSamples)
  all.bic <- rep(NA,nSamples)
  
  saved.scores <- list()
  for(i in 1:le.pheno){
    saved.scores[[i]] <- create.saved.score(node=pheno.col[i], node.nms=pheno.col)
  }
  names(saved.scores) <- pheno.names
        
  M.old <- M0
  ne.old <- nbhd.size(M=M.old)[[1]]
  aux.old <- score.model(M=M.old, saved.scores=saved.scores, cross=cross,
                         addcov=addcov, intcov=intcov, thr=threshold, method = method,
                         verbose = TRUE)
  if(verbose) cat("\n")
  
  bic.old <- aux.old[[1]]
  saved.scores <- aux.old[[2]]
  model.old <- aux.old[[3]]
  post.bic[1] <- all.bic[1] <- bic.old
  post.net.str[,,1] <- M.old
  post.model[1] <- model.old
  cont.accept <- 0
  k <- 2
  aux.thin <- c(1:nSamples)*thinning
  for(i in 2:(nSamples*thinning)){
    M.new <- propose.new.structure(M0=M.old)
    ne.new <- nbhd.size(M=M.new)[[1]]
    aux.new <- score.model(M=M.new, saved.scores=saved.scores, cross=cross,
                           addcov=addcov, intcov=intcov, thr=threshold, method = method,
                           verbose = verbose)
    bic.new <- aux.new[[1]]
    if(i == aux.thin[k]) 
      all.bic[k] <- aux.new[[1]]
    saved.scores <- aux.new[[2]]
    model.new <- aux.new[[3]]
    mr <- exp(-0.5*(bic.new - bic.old))*(ne.old/ne.new)
    if(runif(1) <= min(1,mr)){
      if(i == aux.thin[k]){
        post.bic[k] <- bic.new
        post.net.str[,,k] <- M.new
        post.model[k] <- model.new
        k <- k + 1
      }
      M.old <- M.new
      bic.old <- bic.new
      ne.old <- ne.new
      model.old <- model.new
      cont.accept <- cont.accept + 1
    }
    else{
      if(i == aux.thin[k]){      
        post.bic[k] <- bic.old
        post.net.str[,,k] <- M.old
        post.model[k] <- model.old
        k <- k + 1
      }
    }
    if(verbose) print(c(i,k-1)) 
  }           
  out <- list(post.model=post.model,
              post.bic=post.bic, 
              post.net.str=post.net.str, 
              freq.accept=cont.accept/(nSamples*thinning), 
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
create.saved.score <- function(node,node.nms)
{
  n <- length(node.nms)
  out <- nms <- rep(NA, 2^(n-1))
  nms[1] <- node
  condset <- node.nms[-which(node.nms == node)]
  le <- length(condset)
  k <- 2  
  for(i in 1:le){
    aux <- combn(x=condset,m=i)
    n.comb <- ncol(aux)
    for(j in 1:n.comb){
      nms[k] <- paste(node, paste(aux[,j],collapse=","), sep="|")
      k <- k + 1 
    }
  }
  names(out) <- nms
  out
} 
######################################################################
propose.new.structure <- function(M0)
{
  M <- M0
  le.nodes <- ncol(M0)
  flag <- 0
  while(flag == 0){
    node <- sample(c(1:le.nodes),1)
    move <- sample(c("add","delete","reverse"),1)
    if(move == "add"){
      aux1 <- forbidden.additions(M=M0,node=node)
      if(!is.null(c(aux1$upf,aux1$downf))){
        aux2 <- match(sort(c(aux1$upf,aux1$downf)),c(1:le.nodes))
        aux3 <- c(1:le.nodes)[-c(aux2,node)]
      }
      else{
        aux3 <- c(1:le.nodes)[-node]
      }
      if(length(aux3) == 1){
        M[node,aux3] <- 1
        flag <- 1
      } 
      if(length(aux3) > 1){
        head <- sample(aux3,1)
        M[node,head] <- 1
        flag <- 1
      }
    }
    if(move == "reverse"){
      aux1 <- check.reversions(M=M0,node=node)$allowed
      if(!is.null(aux1)){
        le.rev <- nrow(aux1)
        aux2 <- sample(c(1:le.rev),1)
        aux3 <- aux1[aux2,]
        M[aux3[1],aux3[2]] <- 0
        M[aux3[2],aux3[1]] <- 1
        flag <- 1
      }
    }
    if(move == "delete"){
      aux1 <- which(M0[,node] == 1)
      if(length(aux1) == 1){
        M[aux1,node] <- 0
        flag <- 1
      }
      if(length(aux1) > 1){
        aux2 <- sample(aux1,1)
        M[aux2,node] <- 0
        flag <- 1
      }
    }
  } 
  M
}
######################################################################
nbhd.size <- function(M)
{
  n.deletions <- sum(M)
  n.additions <- 0 
  n.reversions <- 0 
  le <- ncol(M)
  for(j in 1:le){
        add.forbid <- forbidden.additions(M=M,node=j)
        n.additions <- n.additions + (le - 1 - length(c(add.forbid$upf, add.forbid$downf)))
        rev.allow <- check.reversions(M=M,node=j)$allowed
        if(!is.null(rev.allow))
          n.reversions <- n.reversions + nrow( rev.allow )
  }
  nbhd.size <- n.deletions + n.additions + n.reversions
  list(nbhd.size=nbhd.size, 
       n.deletions=n.deletions,
       n.additions=n.additions, 
       n.reversions=n.reversions)
}
######################################################################
######################################################################
forbidden.additions <- function(M,node)
{
  upf <- which(M[node,] == 1)
  downf <- which(M[,node] == 1)
  le <- length(downf)
  if(le > 0){  
    flag <- 0
    while(flag == 0){
      aux1 <- c()
      for(i in 1:le){
        aux1 <- c(aux1,which(M[,downf[i]] == 1))
        aux1 <- unique(aux1)
      }
      new.downf <- sort(unique(c(aux1,downf)))
      if(identical(new.downf,downf)) flag <- 1
      downf <- new.downf
      le <- length(downf)
    }
   }
   if(length(upf)==0) upf <- NULL
   if(length(downf)==0) downf <- NULL   
   list(upf=upf, downf=downf)
} 
######################################################################
check.reversions <- function(M,node)
{
  down <- which(M[,node] == 1)
  le <- length(down)
  if(le == 0){
    allowed <- NULL 
    forbidden <- matrix( c( c(1:nrow(M))[-node],rep(node,nrow(M)-1)), 
                        nrow(M)-1, 2, byrow=FALSE )
  }
  if(le == 1){
    allowed <- matrix(c(down,node),1,2)
    forbidden <- NULL
  }
  if(le > 1){
    forbidden <- matrix(NA,le,2)
    for(k in 1:le){
      tail <- down[k]
      down1 <- down[-k]
      flag <- 0
      while(flag == 0){
        aux1 <- c()
        le <- length(down1)
        for(i in 1:le){
          aux1 <- c(aux1,which(M[,down1[i]] == 1))
          aux1 <- unique(aux1)
        }
        new.down1 <- sort(unique(c(aux1,down1)))
        if(identical(new.down1,down1)) flag <- 1
        down1 <- new.down1
      }
      if(!is.na(match(tail,down1))) 
        forbidden[k,] <- c(tail,node)
    }  
    forbidden <- na.omit(forbidden)
    attr(forbidden,"na.action") <- NULL
    aux2 <- match(forbidden[,1],down)
    if(length(aux2) > 0){
      allowed <- matrix(c(down[-aux2], rep(node,length(down[-aux2]))),
                        length(down[-aux2]),2,byrow=FALSE)
    }
    else{
      allowed <- matrix(c(down, rep(node,length(down))),
                        length(down),2,byrow=FALSE)
      forbidden <- NULL
    }
  } 
  list(allowed=allowed, forbidden=forbidden)   
}
######################################################################
node.parents <- function(M, node, node.nms)
{
  aux <- which(M[,node] == 1)
  if(length(aux) == 0){ 
    parents <- NULL
    identifier <- as.character(node)
  }
  else{
    parents <- aux
    identifier <- paste(node,paste(parents,collapse=","),sep="|")
  }
  list(parents=parents, identifier=identifier)
}
######################################################################
hk.design.matrix <- function(qtlo, cross.type="f2")
{
  nr <- nrow(qtlo$prob[[1]])
  ng <- length(qtlo$prob)
  if(cross.type == "f2"){
    tmp <- unlist(lapply(qtlo$prob, function(x) cbind(x[,1]-x[,3],x[,2])))
    hkm <- matrix(tmp,nr,2*ng)
  }
  if(cross.type == "bc"){
    tmp <- unlist(lapply(qtlo$prob, function(x) x[,1]-x[,2]))
    hkm <- matrix(tmp,nr,ng)
  }
  cbind(rep(1,nr),hkm)
}
######################################################################
find.pheno.names <- function(cross, pheno.col)
{
  if(!is.null(pheno.col)) {
    if(is.list(pheno.col)) {
      for(i in seq(length(pheno.col)))
        if(!is.character(pheno.col[[i]]))
          pheno.col[[i]] <- names(cross$pheno)[pheno.col[[i]]]
    }
    else {
      if(!is.character(pheno.col))
        pheno.col <- names(cross$pheno)[pheno.col]
    }
  }
  pheno.col
}
