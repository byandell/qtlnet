score.model <- function(M, saved.scores, cross, addcov, intcov, threshold, method = "hk",
                        scan.parents = le - 1, verbose = TRUE)
{
  le <- ncol(M)
  model.score <- 0
  mymodel <- rep(NA,le)
  count.score <- 0
  update.scores <- NULL
  
  for(i in 1:le){
    pheno <- node.parents(M, i)
    mymodel[i] <- paste("(",paste(pheno$identifier,")",sep=""),sep="")
    score.pointer <- match(pheno$code, dimnames(saved.scores)[[1]])

    ## Find saved BIC score if already computed.
    bic <- find.score(saved.scores, update.scores, score.pointer, i)

    if(is.na(bic)){
      if(length(pheno$parents) <= scan.parents) {
        ## Find all pheno.col that have the same parents
        ## and have the same addcov, intcov.
        pheno.cols <- which(!(seq(le) %in% pheno$parents))
        if(length(pheno.cols) > 1 & !is.null(addcov))
          pheno.cols <- pheno.cols[sapply(addcov[pheno.cols], agree.covs, addcov[[i]])]
        if(length(pheno.cols) > 1 & !is.null(intcov))
          pheno.cols <- pheno.cols[sapply(intcov[pheno.cols], agree.covs, intcov[[i]])]
      }
      else
        pheno.cols <- i
        
      ## Run pheno.cols through scan.genome.
      bic <- scan.genome(cross, pheno.cols, pheno$parents, addcov[[i]], intcov[[i]],
                         threshold[pheno.cols], method)

      ## Need to get score.pointer for all phenos.
      if(length(pheno.cols) > 1) {
        code <- apply(outer(pheno.cols, pheno$parents, function(x,y) y - (x < y)),
                      1, function(x) sum(2 ^ (x - 1)))
        score.pointer <- match(code, dimnames(saved.scores)[[1]])
      }
      
      ## Update saved scores.
      update.scores <- rbind(update.scores, cbind(score.pointer, pheno.cols, bic))

      ## Print scan info if verbose.
      count.score <- count.score + 1
      if(verbose) {
        if(count.score == 1) cat("\nscan ")
        cat(mymodel[i], "")
      }
      bic <- bic[i == pheno.cols]
    }
    ## Accumulate model score.
    model.score <- model.score + bic
  }
    
  list(model.score = model.score,
       update.scores = update.scores,
       model.name = paste(mymodel,collapse=""))
}
###########################################################################################
agree.covs <- function(x,y) {
  out <- length(x) == length(y)
  if(out)
    out <- all(sort(x) == sort(y))
  out
}
###########################################################################################
find.score <- function(saved.scores, update.scores, score.pointer, i)
{
  bic <- saved.scores[score.pointer, i]
  if(is.na(bic) & !is.null(update.scores)) {
    wh <- which(update.scores[,1] == score.pointer & update.scores[,2] == i)
    if(length(wh))
      bic <- update.scores[wh[1], 3]
  }
  bic 
}
###########################################################################################
scan.genome <- function(cross, pheno.col, pheno.parents, addcov, intcov, threshold, method)
{
  ## This currently scans one phenotype at a time.
  ## Want to extend to all phenotypes not in names(covM.dat).
  ## To do that, assume covariates are same for all i.
  
  ## Design matrices for qtl/scanone. Must be numeric entries.

  ## Design matrix for parent phenotypes
  covM.dat <- pull.pheno.null(cross, pheno.parents)
      
  ## Design matrix for additive covariates.
  addcovM.dat <- covM.dat
  tmp <- unique(c(addcov, intcov))
  if(!is.null(tmp))
    addcovM.dat <- cbind(addcovM.dat, create.cov.matrix(cross, cov.names = tmp))
  ## Design matrix for interactive covariates.
  intcovM.dat <- create.cov.matrix(cross,cov.names=intcov)
  
  ## This is the big time commitment.
  scan <- scanone(cross, pheno.col = pheno.col, addcov=addcovM.dat, intcov=intcovM.dat,
                  method = method)
  
  ## Mainly intersted in this summary to determine QTLs.
  ss <- summary(scan, format = ifelse(length(pheno.col) == 1, "onepheno", "allpeaks"),
                threshold = threshold)

  ## Build data.frame and formula for linear model fit.
  n.pheno <- length(pheno.col)
  bic <- rep(NA, n.pheno)

  ## Phenotype matrices of actual data. Can be mix of numeric and factor.
  addcov.dat <- pull.pheno.null(cross, unique(c(addcov, intcov)))
  intcov.dat <- pull.pheno.null(cross, intcov)
  
  if(nrow(ss)) {
    ## Model may be different for each trait, depending on QTL.
    ## For loop is clumsy, but calculations are pretty fast.
    for(i in seq(n.pheno)) {
      ## Response for linear model.
      y <- cross$pheno[, pheno.col[i]]
      
      signif.lods <- (ss[ , 1 + 2 * i] >= threshold[i])
      le.markers <- sum(signif.lods)
      if(le.markers > 0){
        ## Need to extend this to multiple phenotypes.
        qtlo <- makeqtl(cross, chr = ss[signif.lods, 1], pos = ss[signif.lods, 2 * i], what="prob")
        geno.dat <- hk.design.matrix(qtlo=qtlo)[,-1]
        
        tmp <- paste("add", 1:le.markers, sep = "")
        if(inherits(cross, "f2"))
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

      ## Record BIC.
      bic[i] <- AIC(fm, k = log(length(y)))[1]
    }
  }
  else {
    ## Same model for all traits (no QTL).
    for(i in seq(n.pheno)) {
      if(i == 1) {
        y <- cross$pheno[, pheno.col[1]]

        ## Need only do this once if no QTL.
        form <- set.dat.form(y, covM.dat, addcov.dat, intcov.dat)
        dat <- form$dat
        form <- form$form
        
        ## Fit linear model.
        fm <- lm(form, dat)
      }
      else
        dat$y <- cross$pheno[, pheno.col[i]]
        
      ## Record BIC.
      bic[i] <- AIC(update(fm, data = dat), k = log(length(y)))[1]
    }
  }
  bic
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
