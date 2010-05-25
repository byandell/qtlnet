loci.qtlnet <- function(qtlnet.object, chr.pos=FALSE, ...)
{
  cross <- qtlnet.object$cross
  ## Make sure cross object has genotype probabilities.
  if (!("prob" %in% names(cross$geno[[1]]))) {
      warning("First running calc.genoprob.")
    cross <- calc.genoprob(cross)
  }

  pheno.net.str <- summary(qtlnet.object)$averaged.net

  ## Extract needed attributes from qtlnet.object.
  pheno.nms <- attr(qtlnet.object, "pheno.names")
  addcov <- attr(qtlnet.object, "addcov") 
  intcov <- attr(qtlnet.object, "intcov") 
  thr <- attr(qtlnet.object, "threshold")
  method <- attr(qtlnet.object, "method")

  le <- length(pheno.nms)
  QTLnodes <- list()
  for(i in 1:le){
    pa <- get.parents(pheno=pheno.nms[i], pheno.net.str=pheno.net.str)
    covM.dat <- NULL
    if(!is.null(pa)){
      covM.dat <- data.frame(cross$pheno[,pa])
      names(covM.dat) <- names(cross$pheno)[pa]
    }

    addcovM.dat <- create.cov.matrix(cross,cov.names=addcov[[i]])
    intcovM.dat <- create.cov.matrix(cross,cov.names=intcov[[i]])
    addintcovM.dat <- create.cov.matrix(cross,cov.names=unique(c(addcov[[i]],intcov[[i]])))
    addcov.dat <- data.frame(cross$pheno[,addcov[[i]]])
    names(addcov.dat) <- addcov[[i]]
    intcov.dat <- data.frame(cross$pheno[,intcov[[i]]])
    names(intcov.dat) <- intcov[[i]]      
    addintcov.dat <- data.frame(cross$pheno[,unique(c(addcov[[i]],intcov[[i]]))])
    names(addintcov.dat) <- unique(c(addcov[[i]],intcov[[i]]))
    if(!is.null(covM.dat) & !is.null(addcovM.dat)) 
      aux.addcov <- cbind(covM.dat,addcovM.dat)
    if(!is.null(covM.dat) & is.null(addcovM.dat)) 
      aux.addcov <- covM.dat
    if(is.null(covM.dat) & !is.null(addcovM.dat)) 
      aux.addcov <- addcovM.dat
    if(is.null(covM.dat) & is.null(addcovM.dat)) 
      aux.addcov <- NULL
    scan <- scanone(cross, pheno.col=i, addcov=aux.addcov, intcov=intcovM.dat, method = method)
    ss <- summary(scan,thr=thr[[i]])
    markers <- row.names(ss)
    le.markers <- length(markers)
    if(le.markers > 0){ 
      if(chr.pos){
        QTLnodes[[i]] <- paste(paste(paste("chr",ss[,1],sep=""),"@",sep=""),round(ss[,2],2),sep="")
      }
      else{
        QTLnodes[[i]] <- markers
      }
    }
    if(le.markers == 0){
      QTLnodes[[i]] <- markers
    }   
  }
  names(QTLnodes) <- pheno.nms
  QTLnodes
}
######################################################################
get.parents <- function(pheno, pheno.net.str)
{
  aux1 <- which(pheno.net.str[,2] == pheno)
  pa <- pheno.net.str[pheno.net.str[,2] == pheno, 1]
  if(length(pa) == 0) return(NULL)
  else return(pa)
}
