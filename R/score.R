score.model <- function(M, saved.scores, cross, addcov, intcov, thr, method = "hk",
                        verbose = TRUE)
{
  ## addcov, intcov, thr, are lists
  le <- ncol(M)
  model.score <- 0
  mymodel <- rep(NA,le)
  count.score <- 0
  update.scores <- NULL
  saved.patterns <- dimnames(saved.scores)[[1]]
  for(i in 1:le){
    pheno <- node.parents(M, i)
    mymodel[i] <- paste("(",paste(pheno$identifier,")",sep=""),sep="")
    score.pointer <- match(pheno$code, saved.patterns)
    if(is.na(score.pointer)) {
      ## This should not happen.
      cat(score.pointer,i, pheno$code, "\n")
      browser()
    }
    bic <- saved.scores[score.pointer, i]
    if(is.na(bic)){
      ## Phenotype i is response.
      y <- cross$pheno[,i]

      ## Phenotype matrices of actual data. Can be mix of numeric and factor.
      intcov.dat <- pull.pheno.null(cross, intcov[[i]])
      addcov.dat <- pull.pheno.null(cross, unique(c(addcov[[i]],intcov[[i]])))

      ## Matrix for parent phenotypes of phenotype i.
      covM.dat <- pull.pheno.null(cross, pheno$parents)

      ## Design matrices for qtl/scanone. Must be numeric entries.
      ## Additive covariates: 
      addcovM.dat <- covM.dat
      tmp <- unique(c(addcov[[i]], intcov[[i]]))
      if(!is.null(tmp))
        addcovM.dat <- cbind(addcovM.dat, create.cov.matrix(cross, cov.names = tmp))
      ## Interactive covariates.
      intcovM.dat <- create.cov.matrix(cross,cov.names=intcov[[i]])

      ## This is the big time commitment.
      count.score <- count.score + 1
      if(verbose) {
        if(count.score == 1) cat("\nscan ")
        cat(mymodel[i], "")
      }
      scan <- scanone(cross, pheno.col=i, addcov=addcovM.dat, intcov=intcovM.dat,
                      method = method)

      ## Mainly intersted in this summary to determine QTLs.
      ss <- summary(scan,thr=thr[[i]])

      ## Build data.frame and formula for linear model fit.
      le.markers <- length(row.names(ss))
      if(le.markers > 0){ 
        qtlo <- makeqtl(cross, chr=ss[,1], pos=ss[,2], what="prob")
        geno.dat <- hk.design.matrix(qtlo=qtlo)[,-1]

        tmp <- paste("add", 1:le.markers, sep = "")
        if(inherits(newcross, "f2"))
          tmp <- as.vector(rbind(tmp, paste("dom", 1:le.markers, sep = "")))
        dimnames(geno.dat)[[2]] <- tmp

        form <- set.dat.form(y, covM.dat, addcov.dat, intcov.dat, geno.dat, le.markers)
        dat <- form$dat
        form <- form$form
      }
      else{
        form <- set.dat.form(y, covM.dat, addcov.dat, intcov.dat)
        dat <- form$dat
        form <- form$form
      }

      ## Fit linear model.
      fm <- lm(form, dat)

      ## Update saved scores.
      bic <- AIC(fm, k = log(length(y)))[1]
      update.scores <- rbind(update.scores, c(score.pointer, i, bic))
    }
    ## Accumulate model score.
    model.score <- model.score + bic
  }
  
  list(model.score = model.score,
       update.scores = update.scores,
       model.name = paste(mymodel,collapse=""))
}
###########################################################################################
set.dat.form <- function(y, covM.dat=NULL, addcov.dat=NULL, intcov.dat=NULL,
                         geno.dat=NULL, le.markers = 0)
{
  ## Set up data.frame
  dat <- data.frame(y = y)
  if(!is.null(covM.dat))
    dat <- cbind.data.frame(dat, covM.dat)
  if(!is.null(addcov.dat))
    dat <- cbind.data.frame(dat, addcov.dat)
  if(!is.null(intcov.dat))
    dat <- cbind.data.frame(dat, intcov.dat)
  if(!is.null(geno.dat))
    dat <- cbind.data.frame(dat, geno.dat)

  ## Set up formula.
  form <- myformula(c(names(covM.dat),names(addcov.dat)),
                    names(intcov.dat),
                    le.markers)

  list(dat = dat, form = form)
}
######################################################################
node.parents <- function(M, node)
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
  code <- sum(M[-node, node] * 2 ^ seq(0, nrow(M) - 2))
  list(parents=parents, identifier=identifier, code = code)
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
