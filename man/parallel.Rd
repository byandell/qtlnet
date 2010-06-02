\name{parallel.qtlnet}
\alias{parallel.qtlnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Code to parallelize use of qtlnet}
\description{
This routine calls one of four phases in a parallelized version of qtlnet.
}
\usage{
parallel.qtlnet(dirpath, phase, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dirpath}{
Character string for directory were user can read and write files.
}
  \item{phase}{
Phase of parallelization as number 1 through 4. See details.
}
  \item{\dots}{
    Additional arguments for phases. See details.
  }
}
\details{
See \url{http://www.stat.wisc.edu/~yandell/sysgen/qtlnet} for details of
implementation in progress. The plan is to run qtlnet via Condor
(\url{http://www.condor.org}) to scale up to larger networks, say up to
100 nodes. Phase 1 sets up the problem and determines how many
\code{\link[qtl]{scanone}} calls (via \code{\link{bic.qtlnet}} will be
needed. Phase 2 runs \code{\link[qtl]{scanone}} in a parallelized
fashion. Phase 3 runs the MCMC sequence multiple times when done in
parallel. Phase 4 is the wrapup. There is a Phase 0 to install R/qtl and
R/qtlnet on cluster machines.
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

\url{http://www.stat.wisc.edu/~yandell/sysgen/qtlnet}
}
\seealso{
\code{\link{mcmc.qtlnet}},
\code{\link{bic.qtlnet}}
}
\examples{
  \dontrun{
    parallel.qtlnet("/u/y/a/yandell/public/html/sysgen/qtlnet/condor", 1)
  }
}
\keyword{ ~kwd1 }
