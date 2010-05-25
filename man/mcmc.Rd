\name{mcmc.qtlnet}
\alias{qtlnet}
\alias{mcmc.qtlnet}
\alias{newcross.mcmc}
\alias{newcross}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MCMC across DAGs}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mcmc.qtlnet(cross, pheno.col, threshold, addcov = NULL, intcov = NULL,
  nSamples = 1000, thinning = 1, M0, burning = TRUE, method = "hk",
  random.seed = NULL, verbose = FALSE)
data(newcross.mcmc)
data(newcross)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{
Object of class \code{cross}. See \code{\link[qtl]{read.cross}}.
}
  \item{pheno.col}{
Phenotype identifiers from \code{cross} object. May be numeric, logical
or character.
}
  \item{threshold}{
Scalar or list of thresholds, one per each node.
}
  \item{addcov}{
    Additive covariates for each phenotype (\code{NULL} if not used).
    If entered as scalar or vector (same format as \code{pheno.col}),
    then the same \code{addcov} is used for all
    phenotypes. Altenatively, may be a list of additive covariate identifiers.
}
  \item{intcov}{
Interactive covariates, entered in the same manner as \code{addcov}.
}
  \item{nSamples}{
Number of samples to record.
}
  \item{thinning}{
Thinning rate. Number of MCMC samples is \code{nSamples*thinning}.
}
  \item{M0}{
Matrix of 0s and 1s with initial directed graph of row->col if (row,col)
entry is 1. Cycles are forbidden (e.g. 1s on diagonal or symmetric 1s
across diagonal).
}
\item{burnin}{
  Proportion of MCMC samples to use as burnin. Default is 0.1 if burnin
  is \code{TRUE}. Must be between 0 and 1.
}
\item{method}{
Model fitting method for \code{\link[qtl]{scanone}}.
}
\item{random.seed}{
  Initialization seed for random number generator. Must be \code{NULL}
  (no reset) or positive numeric. Used in \code{\link[base]{set.seed}}.
}
  \item{verbose}{
Print iteration and number of models fit.
  }
}
\details{
  Models are coded compactly as \code{(1)(2|1)(3|1,2,4,5)(4|2)(5|2)}. Each
  parenthetical entry is a of form \code{(node|parents)}.
  }
\value{
  List of class \code{qtlnet}
  \item{post.model}{Model code (see details).}
  \item{post.bic}{Posterior BIC}
  \item{post.net.str}{3-D array with \code{M} matrix for each of the
    \code{nSamples} values of third dimension.}
  \item{freq.accept}{Frequency of acceptance M-H proposals.}
  \item{saved.scores}{Saved LOD score for each phenotype and all
    possible sets of the other phenotypes as parent nodes.}
  \item{all.bic}{}
  \item{cross}{The \code{cross} object with calculated genotype
    probabilities.}
  In addition, a number of attributes are recorded:
  \item{M0}{Initial network matrix.}
  \item{threshold}{threshold list}
  \item{nSamples}{Number of samples saved}
  \item{thinning}{Thinning rate}
  \item{pheno.col}{Phenotype columns.}
  \item{pheno.names}{Phenotype names}
  \item{addcov}{Additive covariate columns.}
  \item{intcov}{Interactive covariate columns.}
  \item{burnin}{Burnin proportion}
  \item{method}{Method used for \code{\link[qtl]{scanone}}.}
  \item{random.seed}{Initial random number generator seed.}
  \item{random.kind}{Random number generator kind from \code{\link[base]{RNGkind}}.}
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
  \code{\link[qtl]{read.cross}},  \code{\link[qtl]{scanone}},  \code{\link[base]{set.seed}}
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

\dontrun{
  newcross.mcmc <- mcmc.qtlnet(newcross, pheno.col = c(1:5),
                               threshold = 3.83,
                               nSamples = 1000, thinning = 20, 
                               random.seed = 92387475, verbose = TRUE)
  save(newcross.mcmc, file = "newcross.mcmc.RData", compress = TRUE)
}
data(newcross.mcmc)
}
\keyword{ ~kwd1 }
