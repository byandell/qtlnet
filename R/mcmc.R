pheno.stuff <- function(cross, pheno.col)
{
  ## I would like to see something like this in R/qtl.
  if (is.character(pheno.col)) {
    num <- find.pheno(cross, pheno.col)
    if (any(is.na(num))) {
      if (sum(is.na(num)) > 1) 
        stop("Couldn't identify phenotypes ",
             paste(paste("\"", pheno.col[is.na(num)], "\"", sep = ""), collapse = " "))
      else stop("Couldn't identify phenotype \"", pheno.col[is.na(num)], "\"")
    }
    pheno.names <- pheno.col
    pheno.col <- num
  }
  else {
    if(is.logical(pheno.col))
      pheno.col <- seq(ncol(cross$pheno))[pheno.col]
    if(!is.numeric(pheno.col))
      stop("Cannot interpret pheno.col:", pheno.col)
    pheno.names <- names(cross$pheno)[pheno.col]
  }
  list(pheno.col = pheno.col, pheno.names = pheno.names)
}
######################################################################
qtlnet.mcmc <- function(M0, cross, thr, nSamples, thinning=1, pheno.col, addcov=NULL,
                        intcov=NULL, burnin=TRUE, method = "hk", random.seed = NULL,
                        verbose=FALSE)
{
  if(!is.null(random.seed)) {
    if(!is.numeric(random.seed))
      stop("random seed must be numeric")
    set.seed(random.seed)
  }
    
  tmp <- pheno.stuff(cross, pheno.col)
  pheno.col <- tmp$pheno.col
  pheno.names <- tmp$pheno.names
  le.pheno <- length(pheno.col)
  
  post.net.str <- array(dim=c(le.pheno,le.pheno,nSamples))
  post.bic <- rep(NA,nSamples)
  post.model <- rep(NA,nSamples)
  all.bic <- rep(NA,nSamples)
  
  saved.scores <- list()
  for(i in 1:le.pheno){
    saved.scores[[i]] <- create.save.score.matrix(node=pheno.col[i], node.nms=pheno.col)
  }
  names(saved.scores) <- pheno.names
        
  M.old <- M0
  ne.old <- nbhd.size(M=M.old)[[1]]
  aux.old <- score.model(M=M.old, saved.scores=saved.scores, cross=cross,
                         addcov=addcov, intcov=intcov, thr=thr, method = method)
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
                           addcov=addcov, intcov=intcov, thr=thr, method = method)
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
       freq.accept=cont.accept/nSamples, 
       saved.scores=saved.scores, 
       all.bic=all.bic)
  attr(out, "M0") <- M0
  attr(out, "threshold") <- thr
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

  out$model.average <- qtlnet.average(out, burnin = burnin)

  ## Attributes of qtlnet object.
  class(out) <- c("qtlnet","list")
  
  out
}
######################################################################
create.save.score.matrix <- function(node,node.nms)
{
  n <- length(node.nms)
  out <- matrix(NA, 2^(n-1), 1)
  nms <- rep(NA, 2^(n-1))
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
  dimnames(out) <- list(nms, "score")
  out
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
score.model <- function(M, saved.scores, cross, addcov, intcov, thr, method = "hk")
{
  ## addcov, intcov, thr, are lists
  le <- ncol(M)
  mod.score <- 0
  mymodel <- rep(NA,le)
  for(i in 1:le){
    pheno <- node.parents(M=M,node=i)
    mymodel[i] <- paste("(",paste(pheno$identifier,")",sep=""),sep="")
    score.pointer <- match(pheno$identifier, dimnames(saved.scores[[i]])[[1]])
    aux.score <- saved.scores[[i]][score.pointer,1]
    if(is.na(aux.score)){
      addcovM.dat <- create.cov.matrix(cross,cov.names=addcov[[i]])
      intcovM.dat <- create.cov.matrix(cross,cov.names=intcov[[i]])
      addintcovM.dat <- create.cov.matrix(cross,cov.names=unique(c(addcov[[i]],intcov[[i]])))
      addcov.dat <- data.frame(cross$pheno[,addcov[[i]]])
      names(addcov.dat) <- addcov[[i]]
      intcov.dat <- data.frame(cross$pheno[,intcov[[i]]])
      names(intcov.dat) <- intcov[[i]]      
      addintcov.dat <- data.frame(cross$pheno[,unique(c(addcov[[i]],intcov[[i]]))])
      names(addintcov.dat) <- unique(c(addcov[[i]],intcov[[i]]))      
      y <- cross$pheno[,i]
      covM.dat <- NULL
      if(!is.null(pheno$parents)){
        covM.dat <- data.frame(cross$pheno[,pheno$parents])
        names(covM.dat) <- names(cross$pheno)[pheno$parents]
      }
      if(!is.null(covM.dat) & !is.null(addcovM.dat)) 
        aux.addcov <- cbind(covM.dat,addcovM.dat)
      if(!is.null(covM.dat) & is.null(addcovM.dat)) 
        aux.addcov <- covM.dat
      if(is.null(covM.dat) & !is.null(addcovM.dat)) 
        aux.addcov <- addcovM.dat
      if(is.null(covM.dat) & is.null(addcovM.dat)) 
        aux.addcov <- NULL
      scan <- scanone(cross, pheno.col=i, addcov=aux.addcov, intcov=intcovM.dat,
                      method = method)
      ss <- summary(scan,thr=thr[[i]])
      markers <- row.names(ss)
      le.markers <- length(markers)
      if(le.markers > 0){ 
        qtlo <- makeqtl(cross,chr=ss[,1],pos=ss[,2],what="prob")
        geno.dat <- hk.design.matrix(qtlo=qtlo)[,-1]
        dimnames(geno.dat)[[2]] <- as.vector(rbind(paste("add",1:le.markers,sep=""),
          paste("dom",1:le.markers,sep="")))
        if(!is.null(covM.dat) & !is.null(addcov.dat) & !is.null(intcov.dat)){
          dat <- data.frame(y, covM.dat, addintcov.dat, geno.dat)
          form <- myformula(addcov=c(names(covM.dat),addcov[[i]]), intcov=intcov[[i]], nQ=le.markers)
        } 
        if(!is.null(covM.dat) & !is.null(addcov.dat) & is.null(intcov.dat)){
          dat <- data.frame(y, covM.dat, addcov.dat, geno.dat)
          form <- myformula(addcov=c(names(covM.dat),addcov[[i]]), intcov=NULL, nQ=le.markers)
        } 
        if(!is.null(covM.dat) & is.null(addcov.dat) & !is.null(intcov.dat)){
          dat <- data.frame(y, covM.dat, intcov.dat, geno.dat)
          form <- myformula(addcov=names(covM.dat), intcov=intcov[[i]], nQ=le.markers)
        } 
        if(is.null(covM.dat) & !is.null(addcov.dat) & !is.null(intcov.dat)){
          dat <- data.frame(y, addintcov.dat, geno.dat)
          form <- myformula(addcov=addcov[[i]], intcov=intcov[[i]], nQ=le.markers)
        } 
        if(!is.null(covM.dat) & is.null(addcov.dat) & is.null(intcov.dat)){
          dat <- data.frame(y, covM.dat, geno.dat)
          form <- myformula(addcov=names(covM.dat), intcov=NULL, nQ=le.markers)
        } 
        if(is.null(covM.dat) & !is.null(addcov.dat) & is.null(intcov.dat)){
          dat <- data.frame(y, addcov.dat, geno.dat)
          form <- myformula(addcov=addcov[[i]], intcov=NULL, nQ=le.markers)
        } 
        if(is.null(covM.dat) & is.null(addcov.dat) & !is.null(intcov.dat)){
          dat <- data.frame(y, intcov.dat, geno.dat)
          form <- myformula(addcov=NULL, intcov=intcov[[i]], nQ=le.markers)
        } 
        if(is.null(covM.dat) & is.null(addcov.dat) & is.null(intcov.dat)){
          dat <- data.frame(y, geno.dat)
          form <- myformula(addcov=NULL, intcov=NULL, nQ=le.markers)
        } 
      }
      else{
        if(!is.null(covM.dat) & !is.null(addcov.dat) & !is.null(intcov.dat)){
          dat <- data.frame(y, covM.dat, addintcov.dat)
          form <- myformula(addcov=c(names(covM.dat),addcov[[i]]), intcov=intcov[[i]], nQ=0)
        } 
        if(!is.null(covM.dat) & !is.null(addcov.dat) & is.null(intcov.dat)){
          dat <- data.frame(y, covM.dat, addcov.dat)
          form <- myformula(addcov=c(names(covM.dat),addcov[[i]]), intcov=NULL, nQ=0)
        } 
        if(!is.null(covM.dat) & is.null(addcov.dat) & !is.null(intcov.dat)){
          dat <- data.frame(y, covM.dat, intcov.dat)
          form <- myformula(addcov=names(covM.dat), intcov=intcov[[i]], nQ=0)
        } 
        if(is.null(covM.dat) & !is.null(addcov.dat) & !is.null(intcov.dat)){
          dat <- data.frame(y, addintcov.dat)
          form <- myformula(addcov=addcov[[i]], intcov=intcov[[i]], nQ=0)
        } 
        if(!is.null(covM.dat) & is.null(addcov.dat) & is.null(intcov.dat)){
          dat <- data.frame(y, covM.dat)
          form <- myformula(addcov=names(covM.dat), intcov=NULL, nQ=0)
        } 
        if(is.null(covM.dat) & !is.null(addcov.dat) & is.null(intcov.dat)){
          dat <- data.frame(y, addcov.dat)
          form <- myformula(addcov=addcov[[i]], intcov=NULL, nQ=0)
        } 
        if(is.null(covM.dat) & is.null(addcov.dat) & !is.null(intcov.dat)){
          dat <- data.frame(y, intcov.dat)
          form <- myformula(addcov=NULL, intcov=intcov[[i]], nQ=0)
        } 
        if(is.null(covM.dat) & is.null(addcov.dat) & is.null(intcov.dat)){
          dat <- data.frame(y)
          form <- myformula(addcov=NULL, intcov=NULL, nQ=0)
        } 
      }
      fm <- lm(form, dat)
      aux.score <- saved.scores[[i]][score.pointer,1] <- AIC(fm,k=log(length(y)))[1]
    }
    mod.score <- mod.score + aux.score
  }
  mymodel <- paste(mymodel,collapse="")
  list(mod.score=mod.score,saved.scores=saved.scores,mymodel=mymodel)
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
