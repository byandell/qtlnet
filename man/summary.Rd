\name{summary.qtlnet}
\alias{print.qtlnet}
\alias{summary.qtlnet}
\alias{print.summary.qtlnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{summary of model average network and posterior table}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{print}{qtlnet}(x, ...)
\method{summary}{qtlnet}(object, ...)
\method{print}{summary.qtlnet}(x, ...)
}
\arguments{
  \item{x,object}{
Object of class \code{qtlnet}.
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
Ann Appl Statist 3: 000-000.
TR 1146R, Department of Statistics, UW-Madison.
\url{http://www.stat.wisc.edu/~yandell/doc/2010/92.AnnApplStat.pdf}
}
\seealso{
\code{\link{qtlnet.mcmc}}
}
\examples{
data(newcross.mcmc)
print(newcross.mcmc)
summary(newcross.mcmc)
}
\keyword{utilities}
