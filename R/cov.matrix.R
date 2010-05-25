myformula <- function(addcov=NULL, intcov=NULL, nQ)
{
  if(!is.null(addcov) & is.null(intcov)){ 
    if(nQ > 0){
      Qnames <- as.vector(rbind(paste("add",1:nQ,sep=""),paste("dom",1:nQ,sep="")))
      myform <- as.formula(paste("y ~ ", paste(c(addcov,Qnames), collapse = "+")))
    }
    else{
      myform <- as.formula(paste("y ~ ", paste(addcov, collapse = "+")))
    }
  }
  if(!is.null(intcov)){
    le <- length(intcov)	
    intaddcov <- unique(c(intcov,addcov))
    if(nQ > 0){
      Qnames <- as.vector(rbind(paste("add",1:nQ,sep=""),paste("dom",1:nQ,sep="")))
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
      Qnames <- as.vector(rbind(paste("add",1:nQ,sep=""),paste("dom",1:nQ,sep="")))
      myform <- as.formula(paste("y ~ ", paste(Qnames, collapse = "+")))
    }
    else{
      myform <- as.formula("y ~ 1")
    }
  }	
  myform
}

create.cov.matrix <- function(cross, cov.names)
{
  if(!is.null(cov.names)){
    myformula <- formula(paste("~",paste(cov.names,collapse="+")))
    out <- model.matrix(myformula,cross$pheno)[,-1]
  }
  else{
    out <- NULL
  }
  out
}
