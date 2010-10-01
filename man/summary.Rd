\name{summary.qtlnet}
\alias{check.qtlnet}
\alias{print.qtlnet}
\alias{summary.qtlnet}
\alias{print.summary.qtlnet}
\title{summary of model average network and posterior table}
\description{
Print and summary methods for qtlnet objects.
}
\usage{
\method{print}{qtlnet}(x, cutoff = 0.01, digits = 3, \dots)
\method{summary}{qtlnet}(object, parent.patterns = FALSE, \dots)
\method{print}{summary.qtlnet}(x, \dots)
check.qtlnet(object, min.prob = 0.9, correct = TRUE, verbose = FALSE, \dots)
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
  \item{parent.patterns}{
    Include summary of parent patterns if \code{TRUE}.
  }
  \item{min.prob}{
    Set the minimum posterior probability for inclusion of an edge.
  }
  \item{correct}{
    Correct \code{min.prob} if \code{TRUE}.
  }
  \item{verbose}{
    Print forbidden edges in model-averaged solution if \code{TRUE}.
  }
  \item{\dots}{
    Other hidden arguments. These include \code{min.prob}, which
    can also be passed to other \code{\link{igraph.qtlnet}} and
    \code{\link{plot.qtlnet}} routines.
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
