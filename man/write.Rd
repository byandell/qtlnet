\name{write.qtlnet}
\alias{write.qtlnet}
\title{write qtlnet as text file}
\description{
Write resulting graph as text file
}
\usage{
write.qtlnet(x, filename, edges, loci.list,
  include.qtl = TRUE, digits = 3, \dots)
}
\arguments{
  \item{x}{
Object of class \code{qtlnet}.
}
\item{filename}{
  Character string with name of text file (usually ends in \code{.txt}).
}
  \item{edges}{
Data frame with first two columns being \code{cause} and \code{effect}
directed phenotype pairs. Typically determined as \code{averaged.net}
element from call to \code{\link{summary.qtlnet}}.
}
  \item{loci.list}{
List of character names of loci by phenotype. Typically determined by
call to \code{\link{loci.qtlnet}}.
}
  \item{include.qtl}{
Include QTL in graph if \code{TRUE} and \code{loci.list} is not \code{NULL}.
}
  \item{digits}{
Number of significant digits for \code{width}.
}
  \item{\dots}{
Additional arguments passed to called routines.
}
}
\details{
Simple write of causal network, for instance to use with Cytoscape.
}
\value{
Invisibly returns data frame that corresponds to saved file.
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
  \code{\link{igrap.qtlnet}}
}
\examples{
write.qtlnet(Pscdbp.qtlnet, "Pscdbp.txt")
}
\keyword{ ~kwd1 }
