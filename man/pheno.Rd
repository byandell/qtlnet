\name{loci.qtlnet}
\alias{loci.qtlnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{QTL architecture per node as list}
\description{
Determines QTL that affect each phenotype conditional on the model-averaged
network and on covariates.
}
\usage{
loci.qtlnet(qtlnet.object, chr.pos = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qtlnet.object}{
Object of class \code{qtlnet}.
}
  \item{chr.pos}{
Include chromsome and position if \code{TRUE}.
}
}
\value{
List containing, for each phenotype in the network, a character vector
of the QTL names as \code{chr@pos}, or pseudomarker name if
\code{chr.pos} is \code{FALSE}.
}
\author{
Brian S. Yandell and Elias Chaibub Neto
}
\references{
Chaibub Neto E, Keller MP, Attie AD, Yandell BS (2010)
Causal Graphical Models in Systems Genetics: a unified
framework for joint inference of causal network and
genetic architecture for correlated phenotypes.
Ann Appl Statist 3: 000-000.
TR 1146R, Department of Statistics, UW-Madison.
\url{http://www.stat.wisc.edu/~yandell/doc/2010/92.AnnApplStat.pdf}
}
\seealso{
\code{\link{newcross.mcmc}}
}
\examples{
loci.qtlnet(newcross.mcmc)
}
\keyword{ ~kwd1 }
