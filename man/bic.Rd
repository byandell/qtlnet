\name{bic.qtlnet}
\alias{bic.qtlnet}
\alias{parents.qtlnet}
\alias{bic.join}
\alias{Pscdbp.bic}
\title{
Pre-compute BIC values for qtlnet sampling.
}
\description{
Pre-compute BIC values for qtlnet sampling to speed up MCMC sampling.
}
\usage{
bic.qtlnet(cross, pheno.col, threshold, addcov = NULL, intcov = NULL,
  max.parents = 3, parents, verbose = TRUE, \dots)
parents.qtlnet(pheno.col, max.parents, codes.only = FALSE)
bic.join(cross, pheno.col, \dots, max.parents = 3)
data(Pscdbp.bic)
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
  \item{max.parents}{
    Maximum number of parents per node. This reduces the complexity of
    graphs and shortens run time. Probably best to consider values of 3-5.
  }
  \item{parents}{
    List containing all possible parents up to \code{max.parents} in
    size. May be a subset
  }
  \item{verbose}{
    Print iteration and number of models fit.
  }
  \item{\dots}{
    Additional arguments passed to internal routines. In the case of
    \code{bic.join}, these are a list of objects produced by
    \code{bic.qtlnet} (see example below).
  }
  \item{codes.only}{
    Return only codes of parents if \code{TRUE}.
  }
}
\details{
The most expensive part of calculations is running
\code{\link[qtl]{scanone}} on each phenotype with parent phenotypes as
covariates. One strategy is to pre-compute the BIC contributions using a
cluster and save them for later use.
}
\value{
  Matrix with columns:
  \item{code}{Binary code as decimal for the parents of a phenotype
    node, excluding the phenotype. Value is between 0 (no parents) and
    \code{2 ^ (length(pheno.col) - 1)}.}
  \item{pheno.col}{Phenotype column in reduced set, in range
    \code{1:length(pheno.col)}.}
  \item{bic}{BIC score for phenotype conditional on parents (and covariates).}
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
}
\seealso{
\code{\link{mcmc.qtlnet}}
}
\examples{
\dontrun{
## Compute all phenotype/parent combinations.
## This takes several minutes.
pheno.col <- 1:13
threshold <- 3.83
Pscdbp.bic <- bic.qtlnet(Pscdbp, pheno.col, threshold)
save(Pscdbp.bic, file = "Pscdbp.bic.RData", compress = TRUE)
}
dim(Pscdbp.bic)

\dontrun{
## Break up into pieces to run on several machines.
max.parents <- 3
parents <- parents.qtlnet(pheno.col, max.parents)
## Here is the number of BIC terms to compute (3083).
sum.parents <- summary(parents)
sum.parents
## Index so that approximately 40 groups of 100 per group.
groups <- 1 + floor(cumsum(sum.parents$n.child) / 100)
tapply(sum.parents$n.child, groups, sum)

## Run 30 smaller calculations (say on separate machines).
bic.group <- list()
for(i in unique(groups)) {
  print(c(i, sum(sum.parents$n.child[groups == i])))

  bic.group[[i]] <- bic.qtlnet(Pscdbp, pheno.col, threshold,
    parents = parents[groups == i])
}

## The saved.scores should be saved to use for subsequent call.
saved.scores <- bic.join(cross, pheno.col, bic.group)
out.qtlnet <- mcmc.qtlnet(cross, pheno.col, threshold = 3.83,
  saved.scores = saved.scores)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
