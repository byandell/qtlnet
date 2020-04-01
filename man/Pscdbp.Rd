\name{Pscdbp}
\alias{Pscdbp}
\alias{Pscdbp.qtlnet}
\docType{data}
\title{
  Cross and qtlnet objects with Ghazalpour et al. (2006) data.
  Only 13 phenotypes are included.
}
\description{
The R/qtl \code{cross} object was created from data at source. The
\code{qtlnet} object was created using \code{\link{mcmc.qtlnet}}.
}
\usage{
data(Pscdbp)
data(Pscdbp.qtlnet)
}
\seealso{
  \code{\link[qtl]{read.cross}},
  \code{\link{mcmc.qtlnet}}
}
\source{
\url{https://horvath.genetics.ucla.edu/coexpressionnetwork/}
}
\references{
Ghazalpour A, Doss S, Zhang B, Wang S, Plaisier C, Castellanos R,
Brozell A, Schadt EE, Drake TA, Lusis AJ, Horvath S (2006) Integrating
genetic and network analysis to characterize genes related to mouse
weight. PLoS Genetics 2: e130-NA.
\url{http://dx.doi.org/10.1371/journal.pgen.0020130}

Chaibub Neto E, Keller MP, Attie AD, Yandell BS (2010)
Causal Graphical Models in Systems Genetics: a unified
framework for joint inference of causal network and
genetic architecture for correlated phenotypes.
Ann Appl Statist 4: 320-339.
\url{http://dx.doi.org/10.1214/09-AOAS288}
}
\examples{
summary(Pscdbp)
\dontrun{
summary(Pscdbp.qtlnet)
}
}
\keyword{datasets}
