\name{qtlnet.pheno}
\alias{qtlnet.pheno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{QTL architecture per node as list}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qtlnet.pheno(pheno.net.str, cross, pheno.nms, addcov, intcov, thr, chr.pos = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pheno.net.str}{
%%     ~~Describe \code{pheno.net.str} here~~
}
  \item{cross}{
%%     ~~Describe \code{cross} here~~
}
  \item{pheno.nms}{
%%     ~~Describe \code{pheno.nms} here~~
}
  \item{addcov}{
%%     ~~Describe \code{addcov} here~~
}
  \item{intcov}{
%%     ~~Describe \code{intcov} here~~
}
  \item{thr}{
%%     ~~Describe \code{thr} here~~
}
  \item{chr.pos}{
%%     ~~Describe \code{chr.pos} here~~
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
\code{\link{newcross.mcmc}}
}
\examples{
## Load object newcross with 5 traits.
data(newcross)

## Load MCMC samples created by qtlnet.mcmc.
data(newcross.mcmc)

qtlnet.pheno(newcross.mcmc, newcross, chr.pos=TRUE)
}
\keyword{ ~kwd1 }
