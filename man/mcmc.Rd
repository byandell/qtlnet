\name{mcmc.qtlnet}
\alias{qtlnet}
\alias{mcmc.qtlnet}
\alias{init.qtlnet}
\title{Sample genetic architecture and QTL network}
\description{
Use MCMC to alternatively sample genetic architecture and QTL network as
directed acyclic graphs (DAGs).
}
\usage{
mcmc.qtlnet(cross, pheno.col, threshold, addcov = NULL, intcov = NULL,
  nSamples = 1000, thinning = 1, max.parents = 3, M0 = NULL,
  burnin = 0.1, method = "hk", random.seed = NULL, init.edges = 0,
  saved.scores = NULL, rev.method = c("nbhd", "single"), verbose = FALSE, \dots)
init.qtlnet(pheno.col, max.parents, init.edges)
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
  \item{max.parents}{
Maximum number of parents to a node. This reduces the complexity of
graphs and shortens run time. Probably best to consider values of 3-5.
}
  \item{M0}{
Matrix of 0s and 1s with initial directed graph of row->col if (row,col)
entry is 1. Cycles are forbidden (e.g. 1s on diagonal or symmetric 1s
across diagonal). Default (if \code{NULL}) is sampled by a call to \code{init.qtlnet};
all 0s if \code{init.edges} = 0 (default).
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
\item{init.edges}{
  Initial number of edges for \code{M0}, to be sampled using
  \code{{init.qtlnet}}. Chosen uniformly from 0 to the number of
  possible edges if set to \code{NULL}.
}
\item{saved.scores}{
  Updated scores, typically pre-computed by
  \code{\link{bic.qtlnet}}.
}
  \item{rev.method}{
Method to use for reversing edges. See details.
  }
  \item{verbose}{
Print iteration and number of models fit.
  }
  \item{\dots}{
Additional arguments. Advanced users may want to supply pre-computed
\code{saved.scores} to speed up calculations.
  }
}
\details{
  Models are coded compactly as \code{(1)(2|1)(3|1,2,4,5)(4|2)(5|2)}. Each
  parenthetical entry is a of form \code{(node|parents)}; these each
  require a model fit, for now with \code{\link[qtl]{scanone}}.

  The \code{\link[qtl]{scanone}} routine is run on multiple
  phenotypes in the network that could all have the same parents. For
  instance, for 5 phenotypes, if \code{(1|2,4)} is sampled, then do
  scanone of this model as well as \code{(3|2,4)} and \code{(5|2,4)}.
  Setting the hidden parameter \code{scan.parents} to a value smaller
  than \code{length(pheno.col) - 1} (default) disallows multiple
  trait scanning with more than that number of parents.

  The \code{saved.scores} parameter can greatly reduce MCMC run time,
  by supplying pre-computed BIC scores. See
  \code{\link{bic.qtlnet}}. Another option is to capture
  \code{saved.scores} from a previous \code{mcmc.qtlnet} run with the
  same phenotypes (and covariates). Caution is advised as only a modest
  amount of checking can be done.

  The \code{init.qtlnet} routine can be used to randomly find an initial
  causal network \code{M0} with up to \code{init.edges} edges.

  MCMC updates include delete, add or reverse edge direction. The early
  version of this method only considered the edge on its own
  (\code{rev.method = "single"}), while the neighborhood method
  (\code{rev.method = "nbhd"}) uses the update
  }
\value{
  List of class \code{qtlnet}
  \item{post.model}{Model code (see details).}
  \item{post.bic}{Posterior BIC}
  \item{Mav}{Model average of \code{M} across MCMC samples.}
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
Ann Appl Statist 4: 320-339.
\url{http://dx.doi.org/10.1214/09-AOAS288}

Grzegorczyk and Husmeier (2008) Improving the structure MCMC sampler for
Bayesian networks by introducing a new edge reversal move.
Mach Learn 71: 265–305.
\url{http://dx.doi.org/10.1007/s10994-008-5057-7}
}
\seealso{
  \code{\link[qtl]{read.cross}},  \code{\link[qtl]{scanone}},
  \code{\link[base]{set.seed}},
  \code{\link{bic.qtlnet}}.
}
\examples{
data(Pscdbp)
\dontrun{
  ## Run of subset of traits. Still takes some time.
  Pscdbp.qtlnet <- mcmc.qtlnet(Pscdbp, pheno.col = c(1,2,4,5,6),
                               threshold = 3.83,
                               nSamples = 1000, thinning = 20, 
                               random.seed = 92387475, verbose = TRUE)
  save(Pscdbp.qtlnet, file = "Pscdbp.qtlnet.RData", compress = TRUE)
}
data(Pscdbp.qtlnet)

\dontrun{
  out.qtlnet <- mcmc.qtlnet(Pscdbp, pheno.col = 1:13,
                            threshold = 3.83,
                            nSamples = 1000, thinning = 20, 
                            random.seed = 92387475, verbose = TRUE,
                            saved.scores = Pscdbp.bic)
}
}
\keyword{ ~kwd1 }
