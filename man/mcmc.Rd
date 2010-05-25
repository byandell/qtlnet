\name{qtlnet.mcmc}
\alias{qtlnet}
\alias{qtlnet.mcmc}
\alias{newcross.mcmc}
\alias{newcross}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MCMC across DAGs}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qtlnet.mcmc(M0, cross, thr, nSamples, thinning = 1, pheno.nms,
  addcov = NULL, intcov = NULL, verbose = FALSE)
data(newcross.mcmc)
data(newcross)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M0}{
%%     ~~Describe \code{M0} here~~
}
  \item{cross}{
%%     ~~Describe \code{cross} here~~
}
  \item{thr}{
%%     ~~Describe \code{thr} here~~
}
  \item{nSamples}{
%%     ~~Describe \code{nSamples} here~~
}
  \item{thinning}{
%%     ~~Describe \code{thinning} here~~
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
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
\examples{
\dontrun{
  ## Load object mycross.
  data(mycross)

  ## Subset to five traits to reduce calculation time.
  newcross <- mycross
  newcross$pheno <- mycross$pheno[,c(1,2,4,5,6)]
  newcross <- calc.genoprob(newcross,step=1)
  save(newcross, file = "newcross.RData", compress = TRUE)
}
data(newcross)

## Parameters for qtlnet.mcmc.
M <- matrix(0,5,5)
thr <- vector(mode="list", length=5)
thr[1:5] <- 3.83
random.seed <- 92387475

\dontrun{
  newcross.mcmc <- qtlnet.mcmc(M0=M, cross=newcross, thr=thr, 
                               addcov=NULL, intcov=NULL, nSamples=1000, thinning=20, 
                               pheno.col = c(1:5),
                               random.seed = random.seed, verbose=TRUE)
  save(newcross.mcmc, file = "newcross.mcmc.RData", compress = TRUE)
}
data(newcross.mcmc)
}
\keyword{ ~kwd1 }
