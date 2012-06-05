\name{generate.qtl}
\alias{generate.qtl.pheno}
\alias{generate.qtl.markers}
\title{Generate QTLs ane phenotypes from cross object}
\description{
Generate QTLs ane phenotype data for individual examples from cross. These are utility routines to
illustrate the examples. They are not meant for users per se.
}
\usage{
generate.qtl.markers(cross, n.phe, nqtl = 3)
generate.qtl.pheno(name, cross, bp, bq, stdev, allqtl,
    burnin = 2000, geno)
}
\arguments{
  \item{cross}{object of class \code{cross}; see
    \code{\link[qtl]{read.cross}}}
  \item{name}{character string for example name}
  \item{bp}{vector or matrix of coefficients for dependencies between
    phenotypes; see cyclic and acyclic examples}
  \item{bq}{vector or matrix of coefficients for QTL effects on
    phenotypes; see cyclic and acyclic examples}
  \item{stdev}{vector of standard deviations per phenotype}
  \item{allqtl}{list of objects of class \code{qtl} produced by
    \code{generate.qtl.sample}}
  \item{burnin}{number of burnin cycles for MCMC; default is 2000}
  \item{geno}{genotypes at markers, typically extracted with
    \code{\link[qtl]{pull.geno}}}
  \item{n.phe}{number of phenotypes}
  \item{nqtl}{number of QTL}
}
\seealso{
\code{\link{acyclic}}, 
\code{\link{cyclica}}, 
\code{\link{cyclicb}}, 
\code{\link{cyclicc}}
}
\examples{
\dontrun{
example(acyclic)
example(cyclica)
example(cyclicb)
example(cyclicc)
}
}
\keyword{datagen}
