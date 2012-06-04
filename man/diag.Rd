\name{dist.qtlnet}
\alias{dist.qtlnet}
\alias{edgematch.qtlnet}
\alias{mds.qtlnet}
\alias{plotbic.qtlnet}
\title{QTL network diagnostic routines}
\description{
Various QTLnet diagnostic routines.
}
\usage{
dist.qtlnet(qtlnet.object, min.prob = 0.9, method = "manhattan", cex = 5)
edgematch.qtlnet(qtlnet.object, min.prob = 0.9, method = "manhattan", cex = 5)
mds.qtlnet(qtlnet.object, min.prob = 0.9, method = "manhattan", cex = 5)
plotbic.qtlnet(x, ..., smooth = TRUE)
}
\arguments{
  \item{qtlnet.object,x}{Object of class \code{qtlnet}.}
  \item{min.prob}{Minimum probability to include edge in network.}
  \item{method}{Distance method to be used between columns of connection
    matrix. Used by \code{\link[stats]{dist}}. (Only used for \code{mds.qtlnet}.)}
  \item{cex}{Character expansion. (Only used for \code{mds.qtlnet},
    scaled by range of BIC values.)}
  \item{smooth}{Use \code{\link[stats]{lowess}} smoother if \code{TRUE}.}
  \item{\dots}{Additional unused arguments.}
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
Ann Appl Statist 4: 320-339.
\url{http://www.stat.wisc.edu/~yandell/doc/2010/92.AnnApplStat.pdf}
}
\seealso{
\code{\link{mcmc.qtlnet}}
}
\examples{
loci.qtlnet(Pscdbp.qtlnet)
}
\keyword{ ~kwd1 }
