\name{parents.qtlnet}
\alias{parents.qtlnet}
\alias{summary.parents.qtlnet}
\alias{size.qtlnet}
\alias{group.qtlnet}
\title{
Determine and group node-parent combinations.
}
\description{
Routines useful for examining the size of node-parent combinations.
}
\usage{
parents.qtlnet(pheno.col, max.parents = 3, codes.only = FALSE)
\method{summary}{parents}(object, \dots)
size.qtlnet(pheno.col, max.parents = 3)
group.qtlnet(pheno.col, max.parents = 3, n.groups,
  group.size = 1000, parents = parents.qtlnet(pheno.col, max.parents))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pheno.col}{
Phenotype identifiers from \code{cross} object. May be numeric, logical
or character.
}
  \item{max.parents}{
    Maximum number of parents per node. This reduces the complexity of
    graphs and shortens run time. Probably best to consider values of 3-5.
  }
  \item{parents}{
    List containing all possible parents up to \code{max.parents} in
    size. May be a subset
  }
  \item{codes.only}{
    Return only codes of parents if \code{TRUE}.
  }
  \item{n.groups}{
    Number of groups for parallel computation. Determined from
    \code{group.size} if missing.
  }
  \item{group.size}{
    Size of groups for parallel computation. See details.
  }
  \item{object}{
    Object of class \code{parent.qtlnet}.
  }
  \item{\dots}{
    Additional arguments ignored.
  }
}
\details{
  The most expensive part of calculations is running
  \code{\link[qtl]{scanone}} on each phenotype with parent phenotypes as
  covariates. One strategy is to pre-compute the BIC contributions using a
  cluster and save them for later use. The \code{parents.qtlnet} routine
  creates a list of all possible parent sets (up to \code{max.parents} in
  size). The \code{size.qtlnet} determines the number of
  \code{\link[qtl]{scanone}} calculations possible for a network with
  nodes \code{pheno.col} and maximum parent size \code{max.parents}. The
  \code{group.qtlnet} groups the parent sets into roughly equal size
  groups for parallel computations. See \code{\link{bic.qtlnet}} for
  further details.
}
\value{
  The \code{size.qtlnet} returns the number of possible
  \code{\link[qtl]{scanone}} computations needed for BIC scores.
  
  The \code{group.qtlnet} produces and index into the parents list
  created by \code{parents.qtlnet}. See details.
  
  The \code{parents.qtlnet} creates a list object with names being the
  \code{code}.
  
  The \code{summary} method for such an object is a data
  frame with row.names being the \code{code}, a binary code as decimal
  for the parents of a phenotype node, excluding the phenotype. Value is
  between 0 (no parents) and \code{2 ^ (length(pheno.col) - 1)}. The
  columns are 
  \item{parents}{Comma-separated string of parents to potential child node.}
  \item{n.child}{Number of possible child nodes to this parent set.}
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
\code{\link{bic.qtlnet}}
}
\examples{
## Restrict to at most 3 parents per node.
pheno.col <- 1:6
max.parents <- 3
size.qtlnet(pheno.col, max.parents)
parents <- parents.qtlnet(pheno.col, max.parents)
summary(parents)

## Allow an arbitrary number (up to 12) of parents per node.
pheno.col <- 1:13
max.parents <- 12
size.qtlnet(pheno.col, max.parents)

## Make ~53 groups of ~1000, for a total of 53248 scanone runs.
parents <- parents.qtlnet(pheno.col, max.parents)
n.child <- summary(parents)$n.child
groups <- group.qtlnet(parents = parents, group.size = 1000)
tapply(n.child, groups, sum)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line