\name{generate.data}
\alias{generate.data}
\alias{generate.data.graph.a}
\alias{generate.data.graph.b}
\alias{generate.data.graph.c}
\alias{generate.data.2or3}
\alias{produce.qtl.sample}
\title{Generate graph data for examples}
\description{
Generate data for individual examples. These are utility routines to

illustrate the examples. They are not meant for users per se.
}
\usage{
generate.data(cross, bp, bq, stdev, allqtl)
generate.data.graph.a(cross, burnin, bq, bp, stdev, geno)
generate.data.graph.b(cross, burnin, bq, bp, stdev, geno)
generate.data.graph.c(cross, burnin, bq, bp, stdev, geno)
generate.data.2or3(cross, bp, bq, stdev, allqtl)
produce.qtl.sample(cross, n.phe, nqtl = 3)
}
\arguments{
  \item{cross}{object of class \code{cross}; see \code{\link[qtl]{read.cross}}}
  \item{bp}{vector or matrix of coefficients for dependencies between
    phenotypes; see cyclic and acyclic examples}
  \item{bq}{vector or matrix of coefficients for QTL effects on
    phenotypes; see cyclic and acyclic examples}
  \item{stdev}{vector of standard deviations per phenotype}
  \item{allqtl}{list of objects of class \code{qtl} produced by
    \code{produce.qtl.sample}}
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
