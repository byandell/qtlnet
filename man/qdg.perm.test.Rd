\name{qdg.perm.test}
\alias{qdg.perm.test}
\alias{summary.qdg.perm.test}
\alias{print.qdg.perm.test}
\title{ Conduct permutation test for LOD score of edge direction on
 directed graph }
\description{
  Conduct permutation test for LOD score of edge direction on
  directed graph.
}
\usage{
qdg.perm.test(cross, nperm, node1, node2, common.cov = NULL,
  DG, QTLs, addcov = NULL, intcov = NULL)
\method{summary}{qdg.perm.test}(object, \dots)
\method{print}{qdg.perm.test}(x, \dots)
}
\arguments{
  \item{cross}{ Object of class \code{cross} (see \code{\link[qtl]{read.cross}}). }
  \item{nperm}{ Number of permutations. }
  \item{node1}{ Character string with name of a phenotype nodes. }
  \item{node2}{ Character string with name of a phenotype nodes. }
  \item{common.cov}{ Character string with name of common phenotype covariates. }
  \item{DG}{ Directed graph of class \code{QDG} }
  \item{QTLs}{ List of objects of class \code{qtl}. }
  \item{addcov}{ Names of additive covariates. Must be valid phenotype
    names in \code{cross}. Expanded to include all \code{intcov} names. }
  \item{intcov}{ Names of additive covariates. Must be valid phenotype
    names in \code{cross}. }
  \item{x,object}{Object of class \code{qdg.perm.test}.}
  \item{\dots}{Additional arguments ignored.}
}


\details{
  \code{qdg.perm.test} performs \code{nperm} permutation-based test
 of LOD score for an 
   edge of a directed graph.

}


\value{List composed by:

  \item{pvalue}{Permutation p-value.}

  \item{obs.lod}{Observed LOD score.}

  \item{PermSample}{Permutation LOD scores sample.}

  \item{node1}{Character string with name of a phenotype nodes.}

  \item{node2}{Character string with name of a phenotype nodes.}

}


\references{ Chaibub Neto et al. (2008) Inferring causal phenotype networks from 
            segregating populations. Genetics 179: 1089-1100.}



\examples{
\dontrun{
example(glxnet)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
