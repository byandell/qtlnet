\name{summary.qtlnet}
\alias{print.qtlnet}
\alias{summary.qtlnet}
\alias{print.summary.qtlnet}
\title{summary of model average network and posterior table}
\description{
Print and summary methods for qtlnet objects.
}
\usage{
\method{print}{qtlnet}(x, cutoff = 0.01, digits = 3, ...)
\method{summary}{qtlnet}(object, ...)
\method{print}{summary.qtlnet}(x, ...)
}
\arguments{
  \item{x,object}{
Object of class \code{qtlnet}.
}
  \item{cutoff}{
Frequency cutoff for model patterns to be displayed. Always shows at
least the most common pattern.
}
  \item{digits}{
Number of digits to display for posterior probabilities on directed edges.
}
  \item{\dots}{
Other arguments currently ignored.
}
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
data(Pscdbp.qtlnet)
print(Pscdbp.qtlnet)
summary(Pscdbp.qtlnet)
}
\keyword{utilities}
