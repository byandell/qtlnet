\name{plot.qdg}
\alias{plot.qdg}
\title{ Plot a directed graph }
\description{
  Plot the outputs of the qdgAlgo and qdgSEM functions.
}
\usage{
\method{plot}{qdg}(x, simple = FALSE, breaks = c(1, 3, 10, 20),
  col = c(pos.color = "green", neg.color = "red",
  pheno.color = "yellow", qtl.color = "magenta"),
  include.qtl = TRUE, \dots)
}
\arguments{
  \item{x}{ list containing the output of the qdgAlgo or qdgSEM functions. }
  \item{simple}{ plot a simple graph (no QTL nodes, color or varying node size or 
        arrow thickness). }
  \item{breaks}{ control the thickness of the arrows according to the values of the
        LOD scores. Must be a vector of length four.}
  \item{col}{ color of the arrows. When plotting from qdgSEM green and red arrows 
        represent positive and negative path coefficients, respectively. Arrows are 
        black by default when plotting from qdgAlgo.}
  \item{include.qtl}{ If TRUE, the QTL nodes are included in the plot. }
  \item{\dots}{ Plot arguments (not used). }
}
\details{
  Create directed graph as an object of Class
 \code{\link[graph]{graphNEL-class}} and 
  plot it. Size of a node is proportional to its degree (number of connections).
}
\seealso{ \code{\link[graph]{graphNEL-class}},
 \code{\link{qdgAlgo}}, \code{\link{qdgSEM}}
 }
\examples{
\dontrun{
example(acyclic)
example(cyclica)
example(cyclicb)
example(cyclicc)}
}
\keyword{ hplot }
