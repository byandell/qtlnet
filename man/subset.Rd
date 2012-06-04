\name{subset.qtlnet}
\alias{subset.qtlnet}
\alias{c.qtlnet}
\alias{best.qtlnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Catenate or subset qtlnet object(s).
}
\description{
Multiple qtlnet objects can be catenated together or subsetted by run.
}
\usage{
\method{subset}{qtlnet}(x, run, \dots)
\method{c}{qtlnet}(\dots)
best.qtlnet(x, burnin = attr(x, "burnin"), wh = which.min(meanbic(x, burnin)))
}
\arguments{
  \item{x}{
Object of class \code{qtlnet}. See \code{\link{mcmc.qtlnet}}.
}
  \item{run}{
Numeric index to desired run. Must be between 0 and number of runs.
}
  \item{burnin}{
Proportion of MCMC samples to be considered as burnin. Taken from
\code{qtlnet} object usually.
}
\item{wh}{Number identifying which model is best.}
  \item{\dots}{
For \code{c.qtlnet}, objects of class \code{qtlnet} to be
joined. Ignored for \code{subset.qtlnet}.
}
}
\value{
Both return an object of class \code{qtlnet}.
}
\details{
The catenation is used by \code{\link{parallel.qtlnet}} in phase 5 to
join together multiple independent MCMC runs. Note that the averaged
network and the frequency of acceptance for a derived subset are only
based on the saved samples, while the original \code{qtlnet} objects
used all samples. Thus catenation and subset are not strictly reversible
functions.

The \code{best.qtlnet} routine picks the run with the best (lowest) BIC
score on average and returns that run as a qtlnet object. It also
produces a trace plot of BIC for all the runs.
}
\author{
Brian Yandell
}
\seealso{
\code{\link{mcmc.qtlnet}}
}
\examples{
\dontrun{
joined <- c(qtlnet1, qtlnet2)
sub1 <- subset(joined, 1)
best <- best.qtlnet(joined)
## qtlnet1 and sub1 should be nearly identical.
}
}
\keyword{utilities}
