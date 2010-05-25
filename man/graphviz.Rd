\name{graphviz.qtlnet}
\alias{Rgraphviz.qtlnet}
\title{qtlnet plot using Rgraphviz}
\description{
Plot inferred causal network using Rgraphviz package.
}
\usage{
Rgraphviz.qtlnet(graph, layout.method = "dot", \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph}{
Object of class \code{\link[igraph]{igraph}}, as produced by \code{\link{graph.qtlnet}}.
}
  \item{layout.method}{
    Character string indicating \code{Rgraphviz} layout method. Choices are
    \code{"dot"}, \code{"neato"}, \code{"twopi"}, \code{"circo"} and
    \code{"fdp"}.
    See \code{y} argument to
\code{\link[Rgraphviz]{plot.graphNEL}}.
}
  \item{\dots}{
Additional arguments passed to internal routines for properties of edges
(\code{lwd=3}) and vertices (\code{fontsize=18}, \code{width=3}, \code{height=1}).
}
}
\details{
Rgraphviz package is very flexible, but it depend on the external
package graphviz. The package igraph is offered as default for plotting
(see \code{\link{igraph.qtlnet}}),
but users can invoke \code{Rgraphviz.qtlnet} if they have the Rgraphviz
package and its external routine properly installed on their system.
}
\value{
Whatever is returned by \code{\link[Rgraphviz]{plot.graph}}.
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
  \code{\link{graph.qtlnet}},
  \code{\link[igraph]{igraph}},
  \code{\link[Rgraphviz]{plot.graphNEL}}.
}
\examples{
\dontrun{
gr <- igraph.qtlnet(Pscdbp.qtlnet)
Rgraphviz.qtlnet(gr)
}
}
\keyword{ ~kwd1 }
