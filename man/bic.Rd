\name{bic.qtlnet}
\alias{bic.qtlnet}
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
}
\details{
  The most expensive part of calculations is running
  \code{\link[qtl]{scanone}} on each phenotype with parent phenotypes as
  covariates. One strategy is to pre-compute the BIC contributions using a
  cluster and save them for later use.
  
  We divide the job into three steps: 1) determine parents and divide into
  reasonable sized groups; 2) compute BIC scores using scanone on a grid
  of computers; 3) compute multiple MCMC runs on a grid of computers. See
  the example for details.
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
  \code{\link{mcmc.qtlnet}},
  \code{\link{parents.qtlnet}}
}
\examples{
pheno.col <- 1:13
max.parents <- 12
size.qtlnet(pheno.col, max.parents)

\dontrun{
## Compute all phenotype/parent combinations.
## This shows how to break up into many smaller jobs.

#########################################
## STEP 1: Preparation. Fast. Needed in steps 2 and 3.

library(qtlnet) ## depends on R/qtl, R/igraph.

pheno.col <- 1:13
max.parents <- 12
threshold <- 3.83

## Load cross object. Here we use internal object.
data(Pscdbp)
## or: load("Pscdbp.RData")
cross <- Pscdbp
## or: cross <- read.cross("Pscdbp.csv", "csv")

## Break up into groups to run on several machines.
## ~53 groups of ~1000, for a total of 53248 scanone runs.
parents <- parents.qtlnet(pheno.col, max.parents)
groups <- group.qtlnet(parents = parents, group.size = 1000)

## Save all relevant objects for later steps.
save(cross, pheno.col, max.parents, threshold, parents, groups,
  file = "Step1.RData", compress = TRUE)

#########################################
## STEP 2: Compute BIC scores. Parallelize.

## NB: Configuration of parallelization determined using Step 1 results.

library(qtlnet) ## depends on R/qtl, R/igraph.

## Load Step 1 computations.
load("Step1.RData")

## Parallelize this:
for(i in seq(nrow(groups)))
{
  ## Pre-compute BIC scores for selected parents.
  bic <- bic.qtlnet(cross, pheno.col, threshold,
    max.parents = max.parents,
    parents = parents[seq(groups[i,1], groups[i,2])])

  save(bic, file = paste("bic", i, ".RData", sep = ""), compress = TRUE)
}

#########################################
## STEP 3: Sample Markov chain (MCMC). Parallelize.

## NB: n.runs sets the number of parallel runs.

n.runs <- 100

library(qtlnet) ## depends on R/qtl, R/igraph.

## Load Step 1 computations.
load("Step1.RData")

## Read in saved BIC scores and combine into one object.
bic.group <- list()
for(i in seq(nrow(groups)))
{
  load(paste("bic", i, ".RData", sep = ""))
  bic.group[[i]] <- bic
}
saved.scores <- bic.join(cross, pheno.col, bic.group)


## Parallelize this:
for(i in seq(n.runs))
{
  ## Run MCMC with randomized initial network.
  mcmc <- mcmc.qtlnet(cross, pheno.col, threshold = threshold,
    max.parents = max.parents, saved.scores = saved.scores, init.edges = NULL)

  save(mcmc, file = paste("mcmc", i, ".RData", sep = ""), compress = TRUE)
}

#########################################
## STEP 4: Combine results for post-processing.

## NB: n.runs needed from Step 3.

n.runs <- 100

library(qtlnet) ## depends on R/qtl, R/igraph.

## Combine outputs together.
outs.qtlnet <- list()
for(i in seq(n.runs))
{
  load(paste("mcmc", i, ".RData", sep = ""))
  outs.qtlnet[[i]] <- mcmc
}

out.qtlnet <- c.qtlnet(outs.qtlnet)
summary(out.qtlnet)
print(out.qtlnet)

## End of parallel example.
#########################################
}

dim(Pscdbp.bic)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line