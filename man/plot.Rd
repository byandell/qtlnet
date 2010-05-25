\name{plot.qtlnet}
\alias{plot.qtlnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot using Rgraphviz}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.qtlnet(pheno.output, marker.list, simple = FALSE, pheno.color = "yellow", qtl.color = "red", include.qtl = TRUE, thick = 3, myfontsize, mynodewidth, mynodeheight, plotstyle = "dot", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pheno.output}{
%%     ~~Describe \code{pheno.output} here~~
}
  \item{marker.list}{
%%     ~~Describe \code{marker.list} here~~
}
  \item{simple}{
%%     ~~Describe \code{simple} here~~
}
  \item{pheno.color}{
%%     ~~Describe \code{pheno.color} here~~
}
  \item{qtl.color}{
%%     ~~Describe \code{qtl.color} here~~
}
  \item{include.qtl}{
%%     ~~Describe \code{include.qtl} here~~
}
  \item{thick}{
%%     ~~Describe \code{thick} here~~
}
  \item{myfontsize}{
%%     ~~Describe \code{myfontsize} here~~
}
  \item{mynodewidth}{
%%     ~~Describe \code{mynodewidth} here~~
}
  \item{mynodeheight}{
%%     ~~Describe \code{mynodeheight} here~~
}
  \item{plotstyle}{
%%     ~~Describe \code{plotstyle} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
  \code{\link[igraph]{tkplot}}
}
\examples{
## Load object newcross with 5 traits.
data(newcross)

## Load MCMC samples created by qtlnet.mcmc.
data(newcross.mcmc)

plot.qtlnet(newcross.mcmc, newcross)
}
\keyword{ ~kwd1 }
