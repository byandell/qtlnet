##############################################################################
##
## $Id: codeQDG.R,v 2007/11/28 byandell Exp $
##
##     Copyright (C) 2007 Elias Chaibub Neto and Brian S. Yandell
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the
## Free Software Foundation; either version 2, or (at your option) any
## later version.
##
## These functions are distributed in the hope that they will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## The text of the GNU General Public License, version 2, is available
## as http://www.gnu.org/copyleft or by writing to the Free Software
## Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
##
## Routines: qdgAlgo, summary.qdgAlgo, print.qdgAlgo
##           qdg.perm.test
##           summary.qdg.sem, print.qdg.sem, summary.qdg.perm.test, print.qdg.perm.test
##############################################################################

## Results of qdgAlgo can be plotted using graph.qdg (using igraph).

qdgAlgo <- function(cross, 
                     phenotype.names, 
                     marker.names, 
                     QTL, 
                     alpha, 
                     n.qdg.random.starts, 
                     addcov = NULL, 
                     intcov = NULL,
                     skel.method="pcskel",
                     udg.order = 2)
{
  if(!inherits(cross, "cross"))
    stop("cross must be an object of class cross")

  ################################
  transformPCtoUDG <- function(PC) {
    edges <- graph::edges(PC@graph)  
    tmp1 <- rep(names(edges), sapply(edges, length))
    tmp2 <- unlist(edges)
    UDG <- data.frame(matrix(NA,length(tmp1), 2))
    UDG[,1] <- as.numeric(tmp1)
    UDG[,2] <- as.numeric(tmp2)
    UDG <- UDG[as.numeric(tmp1) < as.numeric(tmp2), ]
    names(UDG) <- paste("node", 1:2, sep = "")
    UDG$edge <- 1
    class(UDG) <- c("qdg", "data.frame")
    attr(UDG, "edgemode") <- "undirected"
    attr(UDG, "message") <- ""
    attr(UDG, "cont") <- 0
    UDG
  }
  ###################################
  renameUDG <- function(selpheno,UDG) {
    rUDG <- UDG
    n <- length(UDG[,1])
    for(i in 1:n){
      rUDG[i,1] <- selpheno[UDG[i,1]] 
      rUDG[i,2] <- selpheno[UDG[i,2]]
    }
    rUDG
  }
  ########################################################
  myformula <- function(addcov=NULL, intcov=NULL, nQ, dat) {
    if(!is.null(addcov) & is.null(intcov)){ 
      mycovs <- dat[,addcov]
      form <- as.formula(paste(" ~ ", paste(addcov,collapse="+")))
      mycovs <- data.frame(model.matrix(form,dat)[,-1])
      addcov <- names(mycovs)
      if(nQ > 0){
        Qnames <- paste("Q", 1:nQ, sep = "")
        myform <- as.formula(paste("y ~ ", paste(c(addcov,Qnames), collapse = "+")))
      }
      else{
        myform <- as.formula(paste("y ~ ", paste(addcov, collapse = "+")))
      }
    }
    if(!is.null(intcov)){
      le <- length(intcov)	
      intaddcov <- unique(c(intcov,addcov))
      mycovs <- dat[,c(intaddcov)]
      form <- as.formula(paste(" ~ ", paste(intaddcov,collapse="+")))
      mycovs <- data.frame(model.matrix(form,dat)[,-1])
      intaddcov <- names(mycovs)
      if(nQ > 0){
        Qnames <- paste("Q", 1:nQ, sep = "")
        intQnames <- c()
        for(i in 1:le){
          intQnames <- c(intQnames,paste(intaddcov[i], Qnames, sep=":"))
        }
        myform <- as.formula(paste("y ~ ", paste(c(intaddcov,Qnames,intQnames), collapse = "+")))
      }
      else{
        myform <- as.formula(paste("y ~ ", paste(intaddcov, collapse = "+")))
      }
    }
    if(is.null(addcov) & is.null(intcov)){
      if(nQ > 0){
        Qnames <- paste("Q", 1:nQ, sep = "")
        myform <- as.formula(paste("y ~ ", paste(Qnames, collapse = "+")))
        mycovs <- NULL
      }
      else{
        myform <- as.formula("y ~ 1")
        mycovs <- NULL
      }
    }	
    list(myform,mycovs)
  }
  ##########################################################################################################################
  lod.score <- function(cross, node1, node2, qtl.node1, qtl.node2, cov.node1=NULL, cov.node2=NULL,
                        intcov = NULL, artfact.qtl) {
    nQ1 <- length(qtl.node1$chr)
    nQ2 <- length(qtl.node2$chr)
    cov.node1 <- unique(c(intcov,cov.node1))
    cov.node2 <- unique(c(intcov,cov.node2))
    if(nQ1 == 0 & nQ2 > 0) qtl.node1 <- qtl.node2
    if(nQ1 > 0 & nQ2 == 0) qtl.node2 <- qtl.node1
    if(nQ1 == 0 & nQ2 == 0) qtl.node1 <- qtl.node2 <- artfact.qtl
    
    old.fitqtl <- compareVersion(qtlversion(), "1.08-43") < 0
    if(old.fitqtl)
      myfitqtl <- function(cross, pheno.col, ...)
        fitqtl(cross$pheno[[pheno.col]], ...)
    else
      myfitqtl <- function(cross, pheno.col, ...)
        fitqtl(cross, pheno.col, ...)
    
    if( !is.null(cov.node1) & !is.null(cov.node2) ){
      mycovs1 <- data.frame(cross$pheno[,cov.node1])
      names(mycovs1) <- cov.node1
      mycovs2 <- data.frame(cross$pheno[,cov.node2])
      names(mycovs2) <- cov.node2
      node1.col <- find.pheno(cross, pheno=node1)
      node1 <- data.frame(cross$pheno[,node1.col])
      names(node1) <- "node1"
      node2.col <- find.pheno(cross, pheno=node2)
      node2 <- data.frame(cross$pheno[,node2.col])
      names(node2) <- "node2"
      auxf1 <- myformula(addcov=cov.node1, intcov=intcov, nQ=nQ1, dat=mycovs1) 
      
      lf1 <- myfitqtl(cross, node1.col, qtl.node1, 
                      formula = auxf1[[1]], cov=auxf1[[2]],
                      dropone=FALSE)$result.full[10]
      auxf2 <- myformula(addcov=c("node1",cov.node2), intcov=intcov, nQ=nQ2, 
                         dat=data.frame(node1,mycovs2))
      lf2 <- myfitqtl(cross, node2.col, qtl.node2, 
                      formula = auxf2[[1]], cov=auxf2[[2]], 
                      dropone=FALSE)$result.full[10]
      auxg2 <- myformula(addcov=cov.node2, intcov=intcov, nQ=nQ2, dat=mycovs2)
      lg2 <- myfitqtl(cross, node2.col, qtl.node2, 
                      formula = auxg2[[1]], cov=auxg2[[2]], 
                      dropone=FALSE)$result.full[10]
      auxg1 <- myformula(addcov=c("node2",cov.node1), intcov=intcov, nQ=nQ1, 
                         dat=data.frame(node2,mycovs1))
      lg1 <- myfitqtl(cross, node1.col, qtl.node1, 
                      formula = auxg1[[1]], cov=auxg1[[2]], 
                      dropone=FALSE)$result.full[10]
    }
    if( !is.null(cov.node1) & is.null(cov.node2) ){
      mycovs1 <- data.frame(cross$pheno[,cov.node1])
      names(mycovs1) <- cov.node1
      node1.col <- find.pheno(cross, pheno=node1)
      node1 <- data.frame(cross$pheno[,node1.col])
      names(node1) <- "node1"
      node2.col <- find.pheno(cross, pheno=node2)
      node2 <- data.frame(cross$pheno[,node2.col])
      names(node2) <- "node2"
      auxf1 <- myformula(addcov=cov.node1, intcov=intcov, nQ=nQ1, dat=mycovs1) 
      lf1 <- myfitqtl(cross, node1.col, qtl.node1, 
                      formula = auxf1[[1]], cov=auxf1[[2]],
                      dropone=FALSE)$result.full[10]
      auxf2 <- myformula(addcov="node1", intcov=intcov, nQ=nQ2, dat=node1)
      lf2 <- myfitqtl(cross, node2.col, qtl.node2, 
                      formula = auxf2[[1]], cov=auxf2[[2]], 
                      dropone=FALSE)$result.full[10]
      auxg2 <- myformula(addcov=NULL, intcov=intcov, nQ=nQ2, dat=NULL)
      lg2 <- myfitqtl(cross, node2.col, qtl.node2, 
                      formula = auxg2[[1]], cov=auxg2[[2]], 
                      dropone=FALSE)$result.full[10]
      auxg1 <- myformula(addcov=c("node2",cov.node1), intcov=intcov, nQ=nQ1, 
                         dat=data.frame(node2,mycovs1))
      lg1 <- myfitqtl(cross, node1.col, qtl.node1, 
                      formula = auxg1[[1]], cov=auxg1[[2]], 
                      dropone=FALSE)$result.full[10]
    }
    if( is.null(cov.node1) & !is.null(cov.node2) ){
      mycovs2 <- data.frame(cross$pheno[,cov.node2])
      names(mycovs2) <- cov.node2
      node1.col <- find.pheno(cross, pheno=node1)
      node1 <- data.frame(cross$pheno[,node1.col])
      names(node1) <- "node1"
      node2.col <- find.pheno(cross, pheno=node2)
      node2 <- data.frame(cross$pheno[,node2.col])
      names(node2) <- "node2"
      auxf1 <- myformula(addcov=NULL, intcov=intcov, nQ=nQ1, dat=NULL) 
      lf1 <- myfitqtl(cross, node1.col, qtl.node1, 
                      formula = auxf1[[1]], cov=auxf1[[2]],
                      dropone=FALSE)$result.full[10]
      auxf2 <- myformula(addcov=c("node1",cov.node2), intcov=intcov, nQ=nQ2, 
                         dat=data.frame(node1,mycovs2))
      lf2 <- myfitqtl(cross, node2.col, qtl.node2, 
                      formula = auxf2[[1]], cov=auxf2[[2]], 
                      dropone=FALSE)$result.full[10]
      auxg2 <- myformula(addcov=cov.node2, intcov=intcov, nQ=nQ2, dat=mycovs2)
      lg2 <- myfitqtl(cross, node2.col, qtl.node2, 
                      formula = auxg2[[1]], cov=auxg2[[2]], 
                      dropone=FALSE)$result.full[10]
      auxg1 <- myformula(addcov="node2", intcov=intcov, nQ=nQ1, dat=node2)
      lg1 <- myfitqtl(cross, node1.col, qtl.node1, 
                      formula = auxg1[[1]], cov=auxg1[[2]], 
                      dropone=FALSE)$result.full[10]
    }
    if( is.null(cov.node1) & is.null(cov.node2) ){
      node1.col <- find.pheno(cross, pheno=node1)
      node1 <- data.frame(cross$pheno[,node1.col])
      names(node1) <- "node1"
      node2.col <- find.pheno(cross, pheno=node2)
      node2 <- data.frame(cross$pheno[,node2.col])
      names(node2) <- "node2"
      auxf1 <- myformula(addcov=NULL, intcov=intcov, nQ=nQ1, dat=NULL) 
      lf1 <- myfitqtl(cross, node1.col, qtl.node1, 
                      formula = auxf1[[1]], cov=auxf1[[2]],
                      dropone=FALSE)$result.full[10]
      auxf2 <- myformula(addcov="node1", intcov=intcov, nQ=nQ2, dat=node1)
      lf2 <- myfitqtl(cross, node2.col, qtl.node2, 
                      formula = auxf2[[1]], cov=auxf2[[2]], 
                      dropone=FALSE)$result.full[10]
      auxg2 <- myformula(addcov=NULL, intcov=intcov, nQ=nQ2, dat=NULL)
      lg2 <- myfitqtl(cross, node2.col, qtl.node2, 
                      formula = auxg2[[1]], cov=auxg2[[2]], 
                      dropone=FALSE)$result.full[10]
      auxg1 <- myformula(addcov="node2", intcov=intcov, nQ=nQ1, dat=node2)
      lg1 <- myfitqtl(cross, node1.col, qtl.node1, 
                      formula = auxg1[[1]], cov=auxg1[[2]], 
                      dropone=FALSE)$result.full[10]
    }
    lf1+lf2-lg1-lg2
  }
  ######################################################################
  orient.graph.edges <- function(cross,UDG,QTLs,addcov=NULL,intcov=NULL) {
    UDG <- subset(UDG,UDG[,3]==1)
    le <- length(UDG[,1])
    DG <- data.frame(matrix(0,le,4))
    for(i in 1:le){
      node1 <- DG[i,1] <- UDG[i,1]
      node2 <- DG[i,3] <- UDG[i,2]
      s <- lod.score(cross=cross, node1=node1, node2=node2, 
                     qtl.node1=QTLs[[node1]], qtl.node2=QTLs[[node2]], 
                     cov.node1=addcov, cov.node2=addcov, intcov=intcov, artfact.qtl=QTLs[[1]])
      DG[i,2] <- ifelse(s >= 0, "---->", "<----")
      DG[i,4] <- s
    }
    names(DG) <- c("node1", "direction", "node2", "lod score")
    class(DG) <- c("qdg", "data.frame")
    attr(DG, "edgemode") <- "directed"
    DG
  }
  ########################################################################
  recheck.directions <- function(cross,QTLs,oldDG,addcov=NULL,intcov=NULL) {
    DG1 <- DG2 <- oldDG
    aux1 <- "not equal"
    cont <- 0
    le <- length(DG1[,1])
    while((aux1 == "not equal") & (cont < 30)){
      for(i in 1:le){
        DG2 <- orient.graph.edges.cov(cross=cross, QTLs=QTLs, oldDG=DG2, i=i, addcov=addcov, intcov=intcov)
      }
      aux1 <- check.DG(newDG=DG2,oldDG=DG1)
      DG1 <- DG2 
      cont <- cont + 1
    }
    if(aux1 == "equal"){
      newDG <- DG1
      message <- "algorithm converged"
    }
    if(aux1 != "equal"){
      newDG <- oldDG
      message <- "algorithm didn't converge"
    }
    attr(newDG, "edgemode") <- "directed"
    attr(newDG, "message") <- message
    attr(newDG, "cont") <- cont
      newDG
  }
  ##############################################################################
  orient.graph.edges.cov <- function(cross,QTLs,oldDG,i,addcov=NULL,intcov=NULL) {
    newDG <- oldDG
    cov <- get.covariates(pair=i,DG=oldDG)
    node1 <- oldDG[i,1]
    node2 <- oldDG[i,3]
    ls <- lod.score(cross=cross, node1=node1, node2=node2, 
                    qtl.node1=QTLs[[node1]], qtl.node2=QTLs[[node2]], 
                    cov.node1=c(cov[[1]],addcov), cov.node2=c(cov[[2]],addcov), 
                    intcov=intcov, artfact.qtl=QTLs[[1]])
    if(ls >= 0){newDG[i,2] <- "---->"}
    else{newDG[i,2] <- "<----"}
    newDG[i,4] <- ls
    names(newDG) <- c("node1", "direction", "node2", "lod score")
      return(newDG)
  }
  #################################
  check.DG <- function(newDG,oldDG) {
    aux <- all.equal(newDG[,2],oldDG[,2])
    return( ifelse(aux == TRUE, "equal", "not equal") )
  }
  ###################################
  get.covariates <- function(pair,DG) {
    nDG <- length(DG[,1])
    node1 <- DG[pair,1]
    node2 <- DG[pair,3]
    cov1 <- c()
    cov2 <- c()
    for(i in 1:nDG) {
      if((i != pair) & (DG[i,1] == node1) & (DG[i,2] == "<----"))
        cov1 <- c(cov1,DG[i,3])
      if((i != pair) & (DG[i,3] == node1) & (DG[i,2] == "---->"))
        cov1 <- c(cov1,DG[i,1])
      if((i != pair) & (DG[i,1] == node2) & (DG[i,2] == "<----"))
        cov2 <- c(cov2,DG[i,3])
      if((i != pair) & (DG[i,3] == node2) & (DG[i,2] == "---->"))
        cov2 <- c(cov2,DG[i,1])
    }
    return(list(cov1,cov2))
  }
  ##########################
  shuffle.DG <- function(DG) {
    le <- length(DG[,1])
    aux <- sample(c(1:le),le,replace=FALSE)
    return(DG[aux,]) 
  }
  ###################################
  pull.geno.argmax <- function(cross) {
    le <- length(cross$geno)
    all <- cross$geno[[1]]$argmax
    for(i in 2:le){
      aux <- cross$geno[[i]]$argmax
      all <- cbind(all,aux)
    }
    all
  }
  ###############################################################################################################
  get.all.solutions <- function(DG, rc, n.shuffles, cross, QTLs, markers, phenotypes, genotypes,
                                addcov = NULL, intcov = NULL) {
    mylist <- list()
    rc <- order.as(rc)
    mylist[[1]] <- rc
    myloglik <- c()
    myBIC <- c()
    n.arrows <- length(DG[,1])
    aux <- log.likelihood(cross = cross, DG = rc, markers = markers, phenotypes = phenotypes,
                          genotypes = genotypes, addcov = addcov, intcov = intcov) 
    myloglik[1] <- aux[[1]]
    myBIC[1] <- aux[[2]]
    for(i in 2:n.shuffles){
      DG <- shuffle.DG(DG)
      rc <- recheck.directions(cross=cross,QTLs=QTLs,oldDG=DG,addcov=addcov,intcov=intcov)
      if(attr(rc, "message") == "algorithm converged"){
        rc <- order.as(rc)
        mylist[[i]] <- rc
        aux <- log.likelihood(cross = cross, DG = rc, markers = markers, phenotypes = phenotypes,
                              genotypes = genotypes, addcov = addcov, intcov = intcov)	
        myloglik[i] <- aux[[1]]
        myBIC[i] <- aux[[2]]
      }
    }
    ed <- unique(myloglik)
    le <- length(ed)
    newlist <- list()
    loglikelihood <- bic <- rep(0,le)
    for(i in 1:le){
      aux.pos <- which(myloglik==ed[i])
      pos <- aux.pos[1]
      newlist[[i]] <- mylist[[pos]]
      loglikelihood[i] <- myloglik[pos]
      bic[i] <- myBIC[pos]
    }
    outlist <- list(newlist,loglikelihood,bic)
    names(outlist) <- c("solutions","loglikelihood","BIC")
    outlist
  }
  ########################
  order.as <- function(as) {
    aux1 <- row.names(as)
    n <- length(aux1)
    ordered <- data.frame(matrix(0,n,4))
    for(i in 1:n){
      aux2 <- which(aux1==i)
      ordered[i,] <- as[aux2,]
    }
    names(ordered) <- c("node1","direction","node2","lod")
    ordered
  }
  #########################################################################################
  log.likelihood <- function(cross,DG,markers,phenotypes,genotypes,addcov=NULL,intcov=NULL) {
    aux1 <- 0
    aux2 <- 0
    n.phe <- length(phenotypes)
    for(i in 1:n.phe){
      covs <- get.cov.loglik(resp=phenotypes[i],DG=DG)
      aux <- log.likelihood.s(cross=cross, 
                              node=phenotypes[i], markers=markers[[phenotypes[i]]], 
                              cov.node=covs, genotypes=genotypes, 
                              addcov=addcov, intcov=intcov)
      aux1 <- aux1 + aux[[1]]
      aux2 <- aux2 + aux[[2]]
    }
    list(aux1,aux2)
  }
  ###################################
  get.cov.loglik <- function(resp,DG) {
    nDG <- length(DG[,1])
    cov1 <- c()
    for(i in 1:nDG){	
      if((DG[i,1] == resp) & (DG[i,2] == "<----")){
        cov1 <- c(cov1,DG[i,3])
      }
      if((DG[i,3] == resp) & (DG[i,2] == "---->")){
        cov1 <- c(cov1,DG[i,1])
      }
    }
    return(cov1)
  }
  #################################################################################################
  log.likelihood.s <- function(cross, node, markers, cov.node, genotypes, addcov=NULL, intcov=NULL) {
    node.col <- find.pheno(cross, pheno=node)
    y <- data.frame(cross$pheno[,node.col])
    names(y) <- "y"
    nQ <- length(markers)
    if(nQ > 0){
      mygeno <- data.frame(genotypes[,markers])
      nQ <- length(mygeno[1,])
      for(i in 1:nQ){
        mygeno[,i] <- as.factor(mygeno[,i])
      }
      names(mygeno) <- paste("Q",1:nQ,sep="")
      if(!is.null(cov.node)){
        covar.node <- data.frame(cross$pheno[,c(intcov,addcov)],cross$pheno[,find.pheno(cross, pheno=cov.node)])
        names(covar.node) <- c(intcov,addcov,cov.node)
        aux <- myformula(addcov=names(covar.node), intcov=intcov, nQ=nQ, dat=covar.node)
        mylm <- lm(aux[[1]],data.frame(y,aux[[2]],mygeno))
        mylist <- list(logLik(mylm)[[1]],AIC(mylm,k=log(length(y[,1]))))
      }
      else{
        covar.node <- data.frame(cross$pheno[,c(intcov,addcov)])
        names(covar.node) <- c(intcov,addcov)
        aux <- myformula(addcov=names(covar.node), intcov=intcov, nQ=nQ, dat=covar.node)
        if(!is.null(aux[[2]])) mylm <- lm(aux[[1]],data.frame(y,aux[[2]],mygeno))
        if(is.null(aux[[2]])) mylm <- lm(aux[[1]],data.frame(y,mygeno))
        mylist <- list(logLik(mylm)[[1]],AIC(mylm,k=log(length(y[,1]))))
      }
    }
    else{
      if(!is.null(cov.node)){
        covar.node <- data.frame(cross$pheno[,c(intcov,addcov)],cross$pheno[,find.pheno(cross, pheno=cov.node)])
        names(covar.node) <- c(intcov,addcov,cov.node)
        aux <- myformula(addcov=names(covar.node), intcov=intcov, nQ=0, dat=covar.node)
        mylm <- lm(aux[[1]],data.frame(y,aux[[2]]))
        mylist <- list(logLik(mylm)[[1]],AIC(mylm,k=log(length(y[,1]))))
      }
      else{
        covar.node <- data.frame(cross$pheno[,c(intcov,addcov)])
        names(covar.node) <- c(intcov,addcov)
        aux <- myformula(addcov=names(covar.node), intcov=intcov, nQ=0, dat=covar.node)
        if(!is.null(aux[[2]])) mylm <- lm(aux[[1]],data.frame(y,aux[[2]]))
        if(is.null(aux[[2]])) mylm <- lm(aux[[1]],y)
        mylist <- list(logLik(mylm)[[1]],AIC(mylm,k=log(length(y[,1]))))
      }
    }
    mylist
  }
  ##############################################################
  approximate.UDG <- function(Data, alpha, fixed.order = 2) {
    partial.correlation <- function(i, j, k, comb, R){
      RR <- R[c(i, j, comb[, k]), c(i, j, comb[, k])]
      RRinv <- solve(RR)
      D <- diag(1/sqrt(diag(RRinv)))
      return(-D%*%RRinv%*%D)
    }
    pvalue <- function(correlation, n, np){
      tobs <- correlation/sqrt((1-correlation^2)/(n-2-np))
      return(2*pt(abs(tobs), df = n-2-np, lower.tail = FALSE))	
    }
    n <- length(Data[, 1])
    nv <- length(Data[1, ])
    aux.comb <- c(1:nv)
    R <- cor(Data, method = "spearman")
    UDG <- data.frame(matrix(1, nv*(nv-1)/2, 3))
    names(UDG) <- c("node1", "node2", "edge")
    cp <- 1
    for(i in 1:(nv-1)){
      for(j in (i+1):nv){
        UDG[cp, "node1"] <- names(Data)[i]
        UDG[cp, "node2"] <- names(Data)[j]
        order0.cor <- R[i, j]
        pv <- pvalue(correlation = order0.cor, n = n, np = 0)
        order0.ht <- ifelse(pv > alpha, 0, 1)
        if(order0.ht == 0) UDG[cp, "edge"] <- 0
        else{
          if(nv > 2) end <- 0 
          else end <- 1
          order <- 1
          while((order <= fixed.order) & (end == 0)){
            comb <- combn(x = aux.comb[-c(i, j)], m = order)
            nc <- length(comb[1, ])
            k <- 1
            while((k <= nc) & (end == 0)){
              PR <- partial.correlation(i = i, j = j, k = k, comb = comb, R = R)
              pv <- pvalue(correlation = PR[1, 2], n = n, np = order)
              ht <- ifelse(pv > alpha, 0, 1)
              if(ht == 0){
                UDG[cp, "edge"] <- 0
                end <- 1
              }
              k <- k+1
            }
            order <- order+1
            if(nv <= order+1){end <-1}				
          }
        }
        cp <- cp + 1
      }
    }
    class(UDG) <- c("qdg", "data.frame")
    attr(UDG, "edgemode") <- "undirected"
    attr(UDG, "message") <- ""
    attr(UDG, "cont") <- 0
    
    UDG
  }
  #################################################
  pheno.data <- cross$pheno[,phenotype.names]
  if(skel.method == "pcskel"){
    pcskel <- function(pheno.data, alpha) {
      n <- nrow(pheno.data)
      p <- ncol(pheno.data)
      ## define independence test (partial correlations)
      indepTest <- gaussCItest
      ## define sufficient statistics
      suffStat <- list(C = cor(pheno.data), n = n)
      ## estimate Skeleton
      skeleton(suffStat, indepTest, p, alpha)
    }
    
    ## pcskeleton <- pcAlgo(pheno.data, alpha = alpha)
    pcskeleton <- pcskel(pheno.data, alpha)
    UDG <- transformPCtoUDG(pcskeleton)
    UDG <- renameUDG(selpheno=phenotype.names,UDG=UDG)
  }
  else if (skel.method == "udgskel") 
    UDG <- approximate.UDG(Data = pheno.data, alpha = alpha, fixed.order = udg.order)
  else stop("skel.method must either pcskel or udgskel")
  DG <- orient.graph.edges(cross=cross,UDG=UDG,QTLs=QTL,addcov=addcov,intcov=intcov)
  rc <- recheck.directions(cross=cross,QTLs=QTL,oldDG=DG,addcov=addcov,intcov=intcov)
  aux.cross <- argmax.geno(cross)
  genotypes <- pull.geno.argmax(aux.cross)
  as <- get.all.solutions(DG=DG, rc=rc, n.shuffles=n.qdg.random.starts, cross=cross,
                          QTLs=QTL, markers=marker.names, phenotypes=phenotype.names,
                          genotypes=genotypes, addcov=addcov, intcov=intcov)
  best <- which(as$BIC == min(as$BIC))
  mylist <- list(UDG, DG, best, as)
  names(mylist) <- c("UDG","DG","best.lm","Solutions")
  mylist$marker.names <- marker.names
  mylist$phenotype.names <- phenotype.names
  mylist$addcov <- addcov
  class(mylist) <- c("qdgAlgo", "qdg", "list")
  
  mylist
}

summary.qdgAlgo <- function(object, ...)
{
  cat("\n Number of solutions:\n")
  print(length(object$Solutions$BIC))
  cat("\nBest solution:\n")
  print(object$Solutions$solutions[[object$best.lm]])
  bic.lm <- object$Solutions$BIC[object$best.lm]
  cat("\nBIC:\n")
  print(c(lm = bic.lm))
  cat("\nBest solution is solution number:\n")
  print(object$best.lm)
  cat("\nCaution:\n")
  print("If one of the solutions is a cyclic graph you should run qdg.sem in order to score the networks using SEM.")
  invisible()
}

print.qdgAlgo <- function(x, ...) summary(x, ...)

#################################################################################################
qdg.sem <- function(qdgAlgoObject, cross) 
{
  #################################################################################
  score.sem.models <- function(cross,pheno.names,all.solutions,steptol,addcov=NULL) {
    n.sol <- length(all.solutions[[1]])
    mypheno <- cross$pheno[,pheno.names]
    np <- length(mypheno[1,])
    n.paths <- nrow(all.solutions[[1]][[1]])
    semBIC <- rep(NA,n.sol)
    path.coeffs <- matrix(NA,n.paths,n.sol)
    if(!is.null(addcov)){
      addcov <- paste("cross$pheno$",addcov,sep="")
      myresid <- matrix(0,nind(cross),np)
      for(i in 1:np){
        fm <- lm(as.formula(paste("mypheno[,i] ~ ", paste(addcov, collapse = "+"))))
        myresid[,i] <- fm$resid
      }
      mycov <- cov(myresid)
      for(i in 1:n.sol){
        ramMatrix <- create.sem.model(DG=all.solutions[[1]][[i]],pheno.names=pheno.names)	
        mysem <- try(sem(ramMatrix, S = mycov, N = nind(cross), var.names = pheno.names,
                         steptol = steptol, analytic.gradient = FALSE), silent = TRUE)
        if(class(mysem) != "try-error"){
          aux.summary <- try(summary(mysem),silent=TRUE)
          if(class(aux.summary) != "try-error"){ 
            semBIC[i] <- aux.summary$BIC
            path.coeffs[,i] <- include.path.coefficients(sem.summary=aux.summary,output=all.solutions[[1]][[i]])
          }
        }
      }
    }
    else {
      mycov <- cov(mypheno)
      for(i in 1:n.sol){
        ramMatrix <- create.sem.model(DG=all.solutions[[1]][[i]],pheno.names=pheno.names)	
        mysem <- try(sem(ramMatrix, S = mycov, N = nind(cross), var.names = pheno.names,
                         steptol = steptol, analytic.gradient = FALSE), silent = TRUE)
        if(class(mysem) != "try-error"){
          aux.summary <- try(summary(mysem),silent=TRUE)
          if(class(aux.summary) != "try-error"){ 
            semBIC[i] <- aux.summary$BIC
            path.coeffs[,i] <- include.path.coefficients(sem.summary=aux.summary,output=all.solutions[[1]][[i]])
          } 
        }
      }
    }
    
    ## Drop solutions that did not work with sem().
    tmp <- !is.na(semBIC)
    if(!any(tmp)) {
      stop("No qdgAlgo solutions could be fit with sem().")
    }
    if(any(!tmp)) {
      warning(paste(sum(!tmp), "qdgAlgo solutions could not be fit with sem() and were dropped."))
      semBIC <- semBIC[tmp]
      path.coeffs <- path.coeffs[, tmp, drop = FALSE]
      n.sol <- sum(tmp)
      dropped <- which(!tmp)
    }
    else
      dropped <- NULL
    
    output <- data.frame(cbind(semBIC,approx.posterior(semBIC)))
    names(output) <- c("sem.BIC","posterior prob")
    row.names(output) <- paste("model.",1:n.sol,sep="")
    ## if there are ties, returns the first.
    best <- which(output[,2] == max(output[,2]))[1]	
    list(output,path.coeffs[,best], dropped)
  }
  #########################################################
  include.path.coefficients <- function(sem.summary,output) {
    ne <- length(output[,1])
    mypathcoef <- rep(NA,ne)
    aux <- sem.summary$coeff
    aux <- aux[1:ne,]
    for(i in 1:ne){
      if(output[i,2] == "---->") aux1 <- paste(output[i,3], output[i,1], sep=" <--- ")
      if(output[i,2] == "<----") aux1 <- paste(output[i,1], output[i,3], sep=" <--- ")
      aux2 <- match(aux1,aux[,5])
      mypathcoef[i] <- aux[aux2,1]
    }
    mypathcoef
  }
  ############################################
  create.sem.model <- function(DG,pheno.names) {
    n <- length(DG[,1])
    myvector <- c()
    for(i in 1:n){
      aux1 <- which(DG[i,1]==pheno.names)
      aux2 <- which(DG[i,3]==pheno.names)
      if(DG[i,2] == "---->"){
        aux.vector <- c(1,aux2,aux1,i,NA)
      }
      else{aux.vector <- c(1,aux1,aux2,i,NA)}
      myvector <- c(myvector,aux.vector)
    }
    for(i in 1:length(pheno.names)){
      aux.vector <- c(2,i,i,n+i,NA)
      myvector <- c(myvector,aux.vector)
      }
    matrix(myvector,ncol=5,byrow=TRUE)
  }
  ##################################
  approx.posterior <- function(bics) {
    aux <- min(bics)
    round(exp(-0.5*(bics-aux))/sum(exp(-0.5*(bics-aux))),6)
  }

  #################################################
  ss <- score.sem.models(cross = cross,
                         pheno.names = qdgAlgoObject$phenotype.names,
                         all.solutions = qdgAlgoObject$Solutions,
                         steptol = 1 / 100000,
                         addcov = qdgAlgoObject$addcov)
  best <- which(ss[[1]][,1] == min(ss[[1]][,1]))
  mylist <- list(best, ss[[1]], ss[[2]])
  names(mylist) <- c("best.SEM","BIC.SEM","path.coeffs")
  mylist$Solutions <- qdgAlgoObject$Solutions
  mylist$marker.names <- qdgAlgoObject$marker.names
  mylist$phenotype.names <- qdgAlgoObject$phenotype.names
  mylist$dropped <- ss[[3]]
  class(mylist) <- c("qdg.sem", "qdg", "list")
  
  mylist
}

summary.qdg.sem <- function(object, ...)
{
  cat("\nBest SEM solution:\n")
  print(object$Solution$solution[[object$best.SEM]])
  bic.sem <- object$BIC.SEM[object$best.SEM, "sem.BIC"]
  cat("\nBIC:\n")
  print(c(sem = bic.sem))
  cat("\nBest SEM solution is solution number:\n")
  print(object$best.SEM)
  if(!is.null(object$dropped)) {
    cat(length(object$dropped), "qdg.sem solution were dropped; sem() failed for graphs",
        paste(object$dropped, collapse = ","))
  }
  invisible()
}

print.qdg.sem <- function(x, ...) summary(x, ...)

#################################################################################################
qdg.perm.test <- function(cross,nperm,node1,node2,common.cov=NULL,DG,QTLs,addcov=NULL,intcov=NULL)
{
  ###################################
  get.covariates <- function(pair,DG)
    {
      nDG <- length(DG[,1])
      node1 <- DG[pair,1]
      node2 <- DG[pair,3]
      cov1 <- c()
      cov2 <- c()
      for(i in 1:nDG){	
        if((i != pair) & (DG[i,1] == node1) & (DG[i,2] == "<----")){
          cov1 <- c(cov1,DG[i,3])
        }
        if((i != pair) & (DG[i,3] == node1) & (DG[i,2] == "---->")){
          cov1 <- c(cov1,DG[i,1])
        }
        if((i != pair) & (DG[i,1] == node2) & (DG[i,2] == "<----")){
          cov2 <- c(cov2,DG[i,3])
        }
        if((i != pair) & (DG[i,3] == node2) & (DG[i,2] == "---->")){
          cov2 <- c(cov2,DG[i,1])
        }
      }
      list(cov1,cov2)
    }
  ##########################################################################################################################
  lod.score <- function(cross, node1, node2, qtl.node1, qtl.node2, cov.node1=NULL, cov.node2=NULL, intcov=NULL, artfact.qtl)
    {
      nQ1 <- length(qtl.node1$chr)
      nQ2 <- length(qtl.node2$chr)
      cov.node1 <- unique(c(intcov,cov.node1))
      cov.node2 <- unique(c(intcov,cov.node2))
      if(nQ1 == 0 & nQ2 > 0) qtl.node1 <- qtl.node2
      if(nQ1 > 0 & nQ2 == 0) qtl.node2 <- qtl.node1
      if(nQ1 == 0 & nQ2 == 0) qtl.node1 <- qtl.node2 <- artfact.qtl

      old.fitqtl <- compareVersion(qtlversion(), "1.08-43") < 0
      if(old.fitqtl)
        myfitqtl <- function(cross, pheno.col, ...)
          fitqtl(cross$pheno[[pheno.col]], ...)
      else
        myfitqtl <- function(cross, pheno.col, ...)
          fitqtl(cross, pheno.col, ...)

      if( !is.null(cov.node1) & !is.null(cov.node2) ){
        mycovs1 <- data.frame(cross$pheno[,cov.node1])
        names(mycovs1) <- cov.node1
        mycovs2 <- data.frame(cross$pheno[,cov.node2])
        names(mycovs2) <- cov.node2
        node1.col <- find.pheno(cross, pheno=node1)
        node1 <- data.frame(cross$pheno[,node1.col])
        names(node1) <- "node1"
        node2.col <- find.pheno(cross, pheno=node2)
        node2 <- data.frame(cross$pheno[,node2.col])
        names(node2) <- "node2"
        auxf1 <- myformula(addcov=cov.node1, intcov=intcov, nQ=nQ1, dat=mycovs1) 

        lf1 <- myfitqtl(cross, node1.col, qtl.node1, 
                      formula = auxf1[[1]], cov=auxf1[[2]],
                      dropone=FALSE)$result.full[10]
        auxf2 <- myformula(addcov=c("node1",cov.node2), intcov=intcov, nQ=nQ2, 
                           dat=data.frame(node1,mycovs2))
        lf2 <- myfitqtl(cross, node2.col, qtl.node2, 
                      formula = auxf2[[1]], cov=auxf2[[2]], 
                      dropone=FALSE)$result.full[10]
        auxg2 <- myformula(addcov=cov.node2, intcov=intcov, nQ=nQ2, dat=mycovs2)
        lg2 <- myfitqtl(cross, node2.col, qtl.node2, 
                      formula = auxg2[[1]], cov=auxg2[[2]], 
                      dropone=FALSE)$result.full[10]
        auxg1 <- myformula(addcov=c("node2",cov.node1), intcov=intcov, nQ=nQ1, 
                           dat=data.frame(node2,mycovs1))
        lg1 <- myfitqtl(cross, node1.col, qtl.node1, 
                      formula = auxg1[[1]], cov=auxg1[[2]], 
                      dropone=FALSE)$result.full[10]
      }
      if( !is.null(cov.node1) & is.null(cov.node2) ){
        mycovs1 <- data.frame(cross$pheno[,cov.node1])
        names(mycovs1) <- cov.node1
        node1.col <- find.pheno(cross, pheno=node1)
        node1 <- data.frame(cross$pheno[,node1.col])
        names(node1) <- "node1"
        node2.col <- find.pheno(cross, pheno=node2)
        node2 <- data.frame(cross$pheno[,node2.col])
        names(node2) <- "node2"
        auxf1 <- myformula(addcov=cov.node1, intcov=intcov, nQ=nQ1, dat=mycovs1) 
        lf1 <- myfitqtl(cross, node1.col, qtl.node1, 
                      formula = auxf1[[1]], cov=auxf1[[2]],
                      dropone=FALSE)$result.full[10]
        auxf2 <- myformula(addcov="node1", intcov=intcov, nQ=nQ2, dat=node1)
        lf2 <- myfitqtl(cross, node2.col, qtl.node2, 
                      formula = auxf2[[1]], cov=auxf2[[2]], 
                      dropone=FALSE)$result.full[10]
        auxg2 <- myformula(addcov=NULL, intcov=intcov, nQ=nQ2, dat=NULL)
        lg2 <- myfitqtl(cross, node2.col, qtl.node2, 
                      formula = auxg2[[1]], cov=auxg2[[2]], 
                      dropone=FALSE)$result.full[10]
        auxg1 <- myformula(addcov=c("node2",cov.node1), intcov=intcov, nQ=nQ1, 
                           dat=data.frame(node2,mycovs1))
        lg1 <- myfitqtl(cross, node1.col, qtl.node1, 
                      formula = auxg1[[1]], cov=auxg1[[2]], 
                      dropone=FALSE)$result.full[10]
      }
      if( is.null(cov.node1) & !is.null(cov.node2) ){
        mycovs2 <- data.frame(cross$pheno[,cov.node2])
        names(mycovs2) <- cov.node2
        node1.col <- find.pheno(cross, pheno=node1)
        node1 <- data.frame(cross$pheno[,node1.col])
        names(node1) <- "node1"
        node2.col <- find.pheno(cross, pheno=node2)
        node2 <- data.frame(cross$pheno[,node2.col])
        names(node2) <- "node2"
        auxf1 <- myformula(addcov=NULL, intcov=intcov, nQ=nQ1, dat=NULL) 
        lf1 <- myfitqtl(cross, node1.col, qtl.node1, 
                      formula = auxf1[[1]], cov=auxf1[[2]],
                      dropone=FALSE)$result.full[10]
        auxf2 <- myformula(addcov=c("node1",cov.node2), intcov=intcov, nQ=nQ2, 
                           dat=data.frame(node1,mycovs2))
        lf2 <- myfitqtl(cross, node2.col, qtl.node2, 
                      formula = auxf2[[1]], cov=auxf2[[2]], 
                      dropone=FALSE)$result.full[10]
        auxg2 <- myformula(addcov=cov.node2, intcov=intcov, nQ=nQ2, dat=mycovs2)
        lg2 <- myfitqtl(cross, node2.col, qtl.node2, 
                      formula = auxg2[[1]], cov=auxg2[[2]], 
                      dropone=FALSE)$result.full[10]
        auxg1 <- myformula(addcov="node2", intcov=intcov, nQ=nQ1, dat=node2)
        lg1 <- myfitqtl(cross, node1.col, qtl.node1, 
                      formula = auxg1[[1]], cov=auxg1[[2]], 
                      dropone=FALSE)$result.full[10]
      }
      if( is.null(cov.node1) & is.null(cov.node2) ){
        node1.col <- find.pheno(cross, pheno=node1)
        node1 <- data.frame(cross$pheno[,node1.col])
        names(node1) <- "node1"
        node2.col <- find.pheno(cross, pheno=node2)
        node2 <- data.frame(cross$pheno[,node2.col])
        names(node2) <- "node2"
        auxf1 <- myformula(addcov=NULL, intcov=intcov, nQ=nQ1, dat=NULL) 
        lf1 <- myfitqtl(cross, node1.col, qtl.node1, 
                      formula = auxf1[[1]], cov=auxf1[[2]],
                      dropone=FALSE)$result.full[10]
        auxf2 <- myformula(addcov="node1", intcov=intcov, nQ=nQ2, dat=node1)
        lf2 <- myfitqtl(cross, node2.col, qtl.node2, 
                      formula = auxf2[[1]], cov=auxf2[[2]], 
                      dropone=FALSE)$result.full[10]
        auxg2 <- myformula(addcov=NULL, intcov=intcov, nQ=nQ2, dat=NULL)
        lg2 <- myfitqtl(cross, node2.col, qtl.node2, 
                      formula = auxg2[[1]], cov=auxg2[[2]], 
                      dropone=FALSE)$result.full[10]
        auxg1 <- myformula(addcov="node2", intcov=intcov, nQ=nQ1, dat=node2)
        lg1 <- myfitqtl(cross, node1.col, qtl.node1, 
                      formula = auxg1[[1]], cov=auxg1[[2]], 
                      dropone=FALSE)$result.full[10]
      }
      lf1+lf2-lg1-lg2
    }
  ########################################################
  myformula <- function(addcov=NULL, intcov=NULL, nQ, dat)
    {
      if(!is.null(addcov) & is.null(intcov)){ 
        mycovs <- dat[,addcov]
        form <- as.formula(paste(" ~ ", paste(addcov,collapse="+")))
        mycovs <- data.frame(model.matrix(form,dat)[,-1])
        addcov <- names(mycovs)
        if(nQ > 0){
          Qnames <- paste("Q", 1:nQ, sep = "")
          myform <- as.formula(paste("y ~ ", paste(c(addcov,Qnames), collapse = "+")))
        }
        else{
          myform <- as.formula(paste("y ~ ", paste(addcov, collapse = "+")))
        }
      }
      if(!is.null(intcov)){
        le <- length(intcov)	
        intaddcov <- unique(c(intcov,addcov))
        mycovs <- dat[,c(intaddcov)]
        form <- as.formula(paste(" ~ ", paste(intaddcov,collapse="+")))
        mycovs <- data.frame(model.matrix(form,dat)[,-1])
        intaddcov <- names(mycovs)
        if(nQ > 0){
          Qnames <- paste("Q", 1:nQ, sep = "")
          intQnames <- c()
          for(i in 1:le){
            intQnames <- c(intQnames,paste(intaddcov[i], Qnames, sep=":"))
          }
          myform <- as.formula(paste("y ~ ", paste(c(intaddcov,Qnames,intQnames), collapse = "+")))
        }
        else{
          myform <- as.formula(paste("y ~ ", paste(intaddcov, collapse = "+")))
        }
      }
      if(is.null(addcov) & is.null(intcov)){
        if(nQ > 0){
          Qnames <- paste("Q", 1:nQ, sep = "")
          myform <- as.formula(paste("y ~ ", paste(Qnames, collapse = "+")))
          mycovs <- NULL
        }
        else{
          myform <- as.formula("y ~ 1")
          mycovs <- NULL
        }
      }	
      list(myform,mycovs)
    }
  ##################################################
  permuta.block.pheno <- function(cross,node1,node2,common.cov,addcov,intcov)
    {
      pheno <- cross$pheno
      le <- length(pheno[,1])
      aux <- sample(c(1:le),le,replace=FALSE)
      perm.pheno <- pheno
      block <- c(node1,node2,common.cov,addcov,intcov)
      perm.pheno[,block] <- pheno[aux,block]
      perm.cross <- cross
      perm.cross$pheno <- perm.pheno
      perm.cross
    }
  #######################################################
  intersect <- function(x, y) y[match(x, y, nomatch = 0)]
  ############################################################
  pair <- intersect(x = which(DG[,1] == node1), y = which(DG[,3] == node2))
  cov <- get.covariates(pair=pair, DG=DG)
  obs.lod <- DG[pair,4]
  ps <- rep(0,nperm)
  for(i in 1:nperm){
    perm.cross <- permuta.block.pheno(cross,node1,node2,common.cov,addcov,intcov)
    ps[i] <- lod.score(cross=perm.cross, node1=node1, node2=node2, 
                       qtl.node1=QTLs[[node1]], qtl.node2=QTLs[[node2]], 
                       cov.node1=c(cov[[1]],addcov), 
                       cov.node2=c(cov[[2]],addcov),
                       intcov=intcov,artfact.qtl=QTLs[[1]])
  }
  if(obs.lod > 0) pvalue <- length(which(ps >= obs.lod))/nperm
  else pvalue <- length(which(ps <= obs.lod))/nperm	
  mylist <- list(pvalue,obs.lod,ps,node1,node2)
  names(mylist) <- c("pvalue","obs.lod","permSample","node1","node2")
  class(mylist) <- c("qdg.perm.test", "list")

  mylist
}

summary.qdg.perm.test <- function(object, ...)
{
  cat("\nNodes:\n")
  print(c(object$node1,object$node2))
  cat("\nPermutation p-value for direction:\n")
  print(object$pvalue)
  cat("\nObserved direction LOD score:\n")
  print(object$obs.lod)
  invisible()
}

print.qdg.perm.test <- function(x, ...) summary(x, ...)
