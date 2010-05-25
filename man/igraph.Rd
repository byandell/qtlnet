\name{igraph.qtlnet}
\alias{igraph.qtlnet}
\alias{graph.qtlnet}
\alias{plot.qtlnet}
\title{qtlnet plot using igraph}
\description{
Plot inferred causal network using igraph package.
}
\usage{
igraph.qtlnet(x, edges, loci.list,
  pheno.color = "green", qtl.color = "red", vertex.color,
  include.qtl = TRUE, \dots)
\method{plot}{qtlnet}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class \code{qtlnet}.
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
  \item{pheno.color,qtl.color}{
Name of color to use for phenotypes and QTLs, respectively.
}
  \item{vertex.color}{
Vertex colors in order of pheno-pheno edged augmented by
\code{loci.list}, by default determined by \code{pheno.color} and \code{qtl.color}.
}
  \item{include.qtl}{
Include QTL in graph if \code{TRUE} and \code{loci.list} is not \code{NULL}.
}
  \item{\dots}{
Additional arguments passed to called routines.
}
}
\details{
Uses the igraph package to create graph objects. These can be exported
to a variety of other modern graphics packages. \code{graph.qtlnet} is
synonymous with \code{igraph.qtlnet}.
}
\value{
  Object of class \code{\link[igraph]{igraph}}.
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
  \code{\link{summary.qtlnet}},
  \code{\link{loci.qtlnet}},
  \code{\link[igraph]{graph.data.frame}},
  \code{\link[igraph]{tkplot}}
}
\examples{
Pscdbp.graph <- igraph.qtlnet(Pscdbp.qtlnet)
Pscdbp.graph
\dontrun{
tkplot(Pscdbp.graph)
}
}
\keyword{ ~kwd1 }
