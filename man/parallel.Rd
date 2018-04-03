\name{parallel.qtlnet}
\alias{parallel.qtlnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Code to parallelize use of qtlnet}
\description{
This routine calls one of four phases in a parallelized version of qtlnet.
}
\usage{
parallel.qtlnet(phase, index = 1, \dots, dirpath = ".")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phase}{
Phase of parallelization as number 1 through 4. See details.
}
  \item{index}{
    Index for phase. Used in phases 2 and 4, and for error codes saved
    in \code{RESULT.phase.index} file.
  }
  \item{\dots}{
    Additional arguments for phases. See details.
  }
  \item{dirpath}{
Character string for directory were user can read and write files. When
submitting to a cluster, this should remain the default.
}
}
\details{
See \url{http://www.stat.wisc.edu/~yandell/sysgen/qtlnet} for details of
implementation in progress. The plan is to run qtlnet via Condor
(\url{https://research.cs.wisc.edu/htcondor/}) to scale up to larger networks, say up to
100 nodes. Most important information is passed in files. Phase 1
imports arguments from the \code{params.txt} file, which must have
parse-able assignments to the arguments of
\code{qtlnet:::qtlnet.phase1}. This first phase produces file
\code{Phase1.RData}, which included objects used by all other phases.

Phase 1 also creates file \code{groups.txt}, which for each line has
begin and end indices for the \code{parents} that would result from a
call to \code{\link{parents.qtlnet}}. Phase 2 should be run the same
number of times as the number of lines in file \code{groups.txt}. Each
run produces a \code{bicN.RData} file containing BIC computations. These
computations are aggregated in Phase 3 to create \code{Phase3.RData},
which contains the \code{saved.scores} used for
\code{\link{mcmc.qtlnet}} runs in Phase 4, which each produce an
\code{mcmcN.RData} file. The number of runs of Phase 4
is an argument \code{nruns} stored in the \code{params.txt} file
processed in Phase 1. Finally, Phase 5 aggregates the MCMC results from
multiple independent runs into one \code{qtlnet} object.
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
