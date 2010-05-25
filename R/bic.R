######################################################################
bic.qtlnet <- function(cross, pheno.col, threshold,
                       addcov=NULL, intcov=NULL,
                       max.parents = 3,
                       parents = parents.qtlnet(pheno.col, max.parents),
                       verbose = TRUE,
                       ...)
{
  ## Pre-compute BIC for many patterns.

  ## Calculate genotype probabilities if not already done.
  if (!("prob" %in% names(cross$geno[[1]]))) {
    warning("First running calc.genoprob.")
    cross <- calc.genoprob(cross)
  }

  ## Adjust phenotypes and covariates to be numeric.
  cross <- adjust.pheno(cross, pheno.col, addcov, intcov)
  pheno.col <- cross$pheno.col
  pheno.names <- cross$pheno.names
  addcov <- cross$addcov
  intcov <- cross$intcov
  cross <- cross$cross

  n.pheno <- length(pheno.col)

  ## LOD threshold by phenotype.
  if(length(threshold) == 1)
    threshold <- rep(threshold, n.pheno)
  threshold <- as.list(threshold)
  if(length(threshold) != n.pheno)
    stop("threshold must have same length as pheno.col")

  n.scores <- sum(n.pheno - sapply(parents, length))
  update.scores <- matrix(NA, n.scores, 3)
  dimnames(update.scores)[[2]] <- c("code","pheno.col","bic")
  k <- 0

  for(parent.code in names(parents)) {
    if(verbose)
      print(parent.code)
    ## Compute BIC ahead of time.
    ## This method just marches through based on parent pattern.
    ## To do this in parallel, 
    pheno.cols <- which(!(seq(n.pheno) %in% parents[[parent.code]]))
    code <- find.code(pheno.cols, parents[[parent.code]])
    for(i in seq(length(pheno.cols))) {
      bic <- find.bic(code[i], pheno.cols[i], update.scores, ...)
      if(is.na(bic)){
        run <- calc.bic(cross, code[i], pheno.cols[i], parents[[parent.code]],
                        addcov, intcov, threshold, n.pheno,
                        update.scores = update.scores, ...)
        n.run <- nrow(run)
        update.scores[k + seq(n.run), ] <- run
        k <- k + n.run
      }
    }
  }
  update.scores
}
######################################################################
bic.join <- function(cross, pheno.col, ..., max.parents = 3)
{
  pheno.names <- find.pheno.names(cross, pheno.col)
  
  scores <- list(...)
  if(length(scores) == 0) return(NULL)
  if(length(scores) == 1 & is.list(scores[[1]]))
    scores <- scores[[1]]
  
  n.pheno <- length(pheno.names)
  
  codes <- parents.qtlnet(pheno.col[-1], max.parents, TRUE)
  n.codes <- length(codes)
  saved.scores <- matrix(NA, n.codes, n.pheno)
  dimnames(saved.scores) <- list(codes, pheno.names)

  for(i in seq(length(scores))) {
    index <- match(scores[[i]][, "code"], codes) +
      n.codes * (scores[[i]][, "pheno.col"] - 1)
    saved.scores[index] <- scores[[i]][,"bic"]
  }
  saved.scores
}
######################################################################
make.saved.scores <- function(pheno.names, max.parents = 3, saved.scores = NULL,
                              verbose = TRUE, ...)
{
  ## This routine 
  n.pheno <- length(pheno.names)
  max.parents <- max(1, min(max.parents, n.pheno))

  codes <- parents.qtlnet(seq(n.pheno - 1), max.parents, TRUE)

  if(verbose)
    cat("Making saved.scores matrix of size", length(codes), "by", n.pheno, "\n")
  
  out <- matrix(NA, length(codes), n.pheno)
  dimnames(out) <- list(codes, pheno.names)

  if(!is.null(saved.scores)) {
    tmp <- (ncol(saved.scores) == n.pheno)
    if(tmp & n.pheno == 3)
      tmp <- tmp & !all(dimnames(saved.scores)[[2]] == c("code","pheno.col","bic"))
    if(tmp) {
      dim.scores <- dimnames(saved.scores)
      if(!all(pheno.names == dim.scores[[2]]))
        stop("all pheno.names in user-supplied saved.scores must agree")
      wh <- match(dim.scores[[1]], codes, nomatch = 0)
      tmp <- (wh == 0)
      if(any(tmp))
        warning("some codes in user-supplied saved.scores not found")
      if(!all(tmp))
        out[wh,] <- saved.scores[!tmp,]
    }
    else {
      ## The saved.scores are from bic.qtlnet.
      index <- nrow(out)
      index <- match(saved.scores[,1], codes) + index * (saved.scores[,2] - 1)
      out[index] <- saved.scores[,3]
    }
  }
  out
}
######################################################################
parents.qtlnet <- function(pheno.col, max.parents, codes.only = FALSE)
{
  n.pheno <- length(pheno.col)
  
  ## Create possible patterns. Must be a faster way to do this, but so what.
  if(codes.only)
    parents <- "0"
  else
    parents <- list("0" = numeric(0))
  for(i in seq(1, max.parents)) {
    combs <- combn(n.pheno, i)
    name.combs <- as.character(apply(2 ^ (combs - 1), 2, sum))
    if(codes.only)
      parents <- c(parents, name.combs)
    else {
      combs <- as.list(as.data.frame(combs))
      names(combs) <- name.combs
      parents <- c(parents, combs)
    }
  }
  class(parents) <- c("parents.qtlnet", "list")
  attr(parents, "n.pheno") <- n.pheno
  parents
}
################################################################################
print.parents.qtlnet <- function(x, ...) print(summary(x, ...))
################################################################################
summary.parents.qtlnet <- function(object, ...)
{
  n.pheno <- attr(object, "n.pheno")
  
  out <- data.frame(parents = sapply(object, paste, collapse = ","))
  out$n.child <- n.pheno - sapply(object, length)
  out
}
################################################################################
calc.bic <- function(cross, code, pheno.col, parents, addcov, intcov, threshold,
                     n.pheno, method = "hk", min.nonmissing = 50, ...)
{
  
  ## Match phenotypes with same parents to run scan across genome.
  run <- match.parents(cross, code, pheno.col, parents, addcov, intcov,
                       n.pheno, ...)
  
  ## Run pheno.cols through scanone in scan.genome.
  if(sum(!run$na) < min.nonmissing) {
    warning(paste("BIC set to 0: only ", sum(!run$na),
                  " non-missing values for combination: (",
                  paste(run$pheno.cols, collapse = ","), "|",
                  paste(parents, collapse = ","), ")", sep = ""))
    bic <- rep(0, length(run$pheno.cols))
  }
  else
    bic <- scan.genome(subset(cross, ind = !run$na),
                       run$pheno.cols, parents,
                       addcov[[pheno.col]], intcov[[pheno.col]],
                       threshold[run$pheno.cols], method)
  run$bic <- bic
  cbind(code = as.numeric(run$code),
        pheno.col = run$pheno.cols,
        bic = bic)
}
###########################################################################################
find.code <- function(pheno.cols, pheno.parents)
{
  apply(outer(pheno.cols, pheno.parents, function(x,y) y - (x < y)),
        1, function(x) sum(2 ^ (x - 1)))
}
###########################################################################################
match.parents <- function(cross, code, i, parents, addcov, intcov,
                          n.pheno,
                          scan.parents = n.pheno - 1, ...)
{
  ## Find which individuals have missing phenotype for i, parents, covariates.
  ## na.pheno should not be all TRUE.
  na.pheno <- is.na(cross$pheno[[i]])
  na.covar <- rep(FALSE, nind(cross))
  if(length(parents))
    na.covar <- {na.covar |
                 apply(cross$pheno[, parents, drop = FALSE], 1,
                       function(x) any(is.na(x)))}
  if(!is.null(addcov[[i]]))
    na.covar <- {na.covar |
                 apply(cross$pheno[, addcov[[i]], drop = FALSE], 1,
                       function(x) any(is.na(x)))}
  if(!is.null(intcov[[i]]))
    na.covar <- {na.covar |
                 apply(cross$pheno[, intcov[[i]], drop = FALSE], 1,
                       function(x) any(is.na(x)))}

  scan.parents <- max(1, min(scan.parents, n.pheno - 1))
  if(length(parents) <= scan.parents) {
    ## Find all pheno.col that have the same parents
    pheno.cols <- which(!(seq(n.pheno) %in% parents))

    ## Need to eliminate any that already have saved scores.
    code <- find.code(pheno.cols, parents)
    bic <- find.bic(code, pheno.cols, ...)
    pheno.cols <- pheno.cols[is.na(bic)]

    ## Make sure they have the same addcov, intcov.
    if(length(pheno.cols) > 1 & !is.null(addcov))
      pheno.cols <- pheno.cols[sapply(addcov[pheno.cols], agree.covs, addcov[[i]])]
    if(length(pheno.cols) > 1 & !is.null(intcov))
      pheno.cols <- pheno.cols[sapply(intcov[pheno.cols], agree.covs, intcov[[i]])]

    if(length(pheno.cols) > 1) {

      ## Find out if i and other pheno.cols overlap with NA pattern.
      na.other <- apply(cross$pheno[, pheno.cols, drop = FALSE], 2,
                        function(x, na.pheno, na.covar) {
                          na.x <- is.na(x)
                          any((na.x != na.pheno)[!na.covar])
                        },
                        na.pheno, na.covar)
      ## Should pick up at least i, maybe others.
      pheno.cols <- pheno.cols[!na.other]
    }
  }
  else
    pheno.cols <- i
  
  ## Need to get code for all phenos if more than one..
  if(length(pheno.cols) > 1)
    code <- find.code(pheno.cols, parents)

  list(code = code, pheno.cols = pheno.cols, na = na.pheno | na.covar)
}
######################################################################
adjust.pheno <- function(cross, pheno.col, addcov, intcov)
{
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

  list(cross = cross, pheno.col = pheno.col, pheno.names = pheno.names,
       addcov = addcov, intcov = intcov)
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
