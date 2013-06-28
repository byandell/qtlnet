\name{qdg}
\alias{qdg}
\alias{summary.qdg}
\alias{print.qdg}
\alias{graph.qdg}
\title{ Produces a directed graph using QDG algorithm }
\description{
  This function implements the QDG algorithm described in Chaibub Neto et al 2008. It
  creates and scores QDGs. The computed scores (log-likelihood and BIC) are only 
  valid for acyclic graphs. For cyclic networks qdgSEM should be used to compute the 
  scores.
}
\usage{
qdg(cross, phenotype.names, marker.names, QTL, alpha, 
    n.qdg.random.starts, addcov = NULL, intcov = NULL, 
    skel.method = c("pcskel","udgskel"), udg.order = 2)
graph.qdg(x, \dots)
\method{print}{qdg}(x, \dots)
\method{summary}{qdg}(object, \dots)
}
\arguments{
  \item{cross}{ object of class \code{cross} (see \code{\link[qtl]{read.cross}}). }
  \item{phenotype.names}{ character string with names of phenotype nodes 
        corresponding to phenotypes in \code{cross}. }
  \item{marker.names}{ list of character strings, one for each of
        \code{phenotype.names}. Each character string has the marker names
        for that phenotype. }
  \item{QTL}{ object of class \code{qtl} (see \code{\link[qtl]{makeqtl}}). }
  \item{alpha}{ significance level threshold for PC or UDG algorithms (for the 
        inference of the graph skeleton. See step 1 of the QDG algorithm). Must be
        between 0 and 1. }
  \item{n.qdg.random.starts}{ number of random starts for the QDG algorithm (see 
        step 3 of the QDG algorithm). }
  \item{addcov}{ names of additive covariates. Must be valid phenotype names in 
        \code{cross}. Expanded to include all \code{intcov} names. }
  \item{intcov}{ names of additive covariates. Must be valid phenotype names in 
        \code{cross}. }
  \item{skel.method}{ Either "pcskel" for the PC skeleton algorithm
    (\code{\link[pcalg]{skeleton}}) or "udgskel" for the UDG algorithm
    (\code{approximate.UDG} routine defined internal to \code{qdg}). }
  \item{udg.order}{ maximum allowed order of the UDG algorithm. Must be between zero 
        and the number of variables minus 2. }
  \item{x,object}{object of class \code{qdg}.}
  \item{\dots}{additional arguments (ignored).}
}
\details{
  The log-likelihood and BIC scores are computed based in the factorization of the 
  joint distribution, and hence are only valid for acyclic networks. For cyclic 
  networks these scores are relative to the unnormalized likelihoods. Models include
  phenotypes and QTLs. The 'udgskel' method for the computation of the skeleton of 
  the causal model should be used for small networks only (the UDG algorithm quickly 
  becomes computationally infeasible as the number of nodes increases).
}
\value{List object that inherits class "qdg" and "qdg" with components:
  \item{UDG}{Undirected dependency graph from PC skeleton or UDG algorithms.}
  \item{DG}{Directed dependency graph before recheck step (output of the step 2 of the        QDG algorithm).} 
  \item{best.lm}{Solution with lowest BIC (best fit to the data).}
  \item{Solutions}{Solutions of dependency graph after recheck step (output of steps 3, 4
        and 5 of the QDG algorithm.)}
  \item{marker.names}{List of character strings, one for each of
 \code{phenotype.names}. 
        Each character string has the marker names
 for that phenotype.}
  \item{phenotype.names}{Character string with names of phenotype nodes
        corresponding to phenotypes in \code{cross}.}
}
\references{Chaibub Neto et al. (2008) Inferring causal phenotype networks from 
            segregating populations. Genetics 179: 1089-1100.}
\seealso{ \code{\link[pcalg]{skeleton}} }
\examples{
## simulate a genetic map (20 autosomes, 10 not equaly spaced markers per 
## chromosome)
mymap <- sim.map(len=rep(100,20), n.mar=10, eq.spacing=FALSE, include.x=FALSE)

## simulate an F2 cross object with n.ind (number of individuals)
n.ind <- 200
mycross <- sim.cross(map=mymap, n.ind=n.ind, type="f2")

## produce multiple imputations of genotypes using the 
## sim.geno function. The makeqtl function requires it,
## even though we are doing only one imputation (since 
## we don't have missing data and we are using the 
## genotypes in the markers, one imputation is enough)
mycross <- sim.geno(mycross,n.draws=1)

## sample markers (2 per phenotype)
genotypes <- pull.geno(mycross)
geno.names <- dimnames(genotypes)[[2]]
m1 <- sample(geno.names,2,replace=FALSE)
m2 <- sample(geno.names,2,replace=FALSE)
m3 <- sample(geno.names,2,replace=FALSE)
m4 <- sample(geno.names,2,replace=FALSE)

## get marker genotypes
g11 <- genotypes[,m1[1]]; g12 <- genotypes[,m1[2]]
g21 <- genotypes[,m2[1]]; g22 <- genotypes[,m2[2]]
g31 <- genotypes[,m3[1]]; g32 <- genotypes[,m3[2]]
g41 <- genotypes[,m4[1]]; g42 <- genotypes[,m4[2]]

## generate phenotypes
y1 <- runif(3,0.5,1)[g11] + runif(3,0.5,1)[g12] + rnorm(n.ind)
y2 <- runif(3,0.5,1)[g21] + runif(3,0.5,1)[g22] + rnorm(n.ind)
y3 <- runif(1,0.5,1) * y1 +  runif(1,0.5,1) * y2 + runif(3,0.5,1)[g31] +
      runif(3,0.5,1)[g32] + rnorm(n.ind)
y4 <- runif(1,0.5,1) * y3 + runif(3,0.5,1)[g41] + runif(3,0.5,1)[g42] +
      rnorm(n.ind)

## incorporate phenotypes to cross object
mycross$pheno <- data.frame(y1,y2,y3,y4)

## create markers list
markers <- list(m1,m2,m3,m4)
names(markers) <- c("y1","y2","y3","y4")

## create qtl object
allqtls <- list()
m1.pos <- find.markerpos(mycross, m1)
allqtls[[1]] <- makeqtl(mycross, chr = m1.pos[,"chr"], pos = m1.pos[,"pos"])
m2.pos <- find.markerpos(mycross, m2)
allqtls[[2]] <- makeqtl(mycross, chr = m2.pos[,"chr"], pos = m2.pos[,"pos"])
m3.pos <- find.markerpos(mycross, m3)
allqtls[[3]] <- makeqtl(mycross, chr = m3.pos[,"chr"], pos = m3.pos[,"pos"])
m4.pos <- find.markerpos(mycross, m4)
allqtls[[4]] <- makeqtl(mycross, chr = m4.pos[,"chr"], pos = m4.pos[,"pos"])
names(allqtls) <- c("y1","y2","y3","y4")

## infer QDG 
out <- qdg(cross=mycross, 
		phenotype.names = c("y1","y2","y3","y4"), 
		marker.names = markers, 
		QTL = allqtls, 
		alpha = 0.005, 
		n.qdg.random.starts=10, 
		skel.method="pcskel")

out
gr <- graph.qdg(out)
gr
plot(gr)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
