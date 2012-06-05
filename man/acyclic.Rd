\name{acyclic}
\alias{acyclic}
\alias{acyclic.data}
\alias{acyclic.DG}
\alias{acyclic.qdg}
\alias{acyclic.qtl}
\title{Acyclic graph example}
\description{
  We generate synthetic data (sample size 300) according to a DAG composed
  by 100 nodes and 107 edges (exactly as in Figure 1). 
  Each phenotype node is affected by three 
  QTLs, and we allow only additive genetic effects. The QTLs for each phenotype are 
  randomly
  selected among 200 markers, with 10 markers unevenly distributed on each of 20 autosomes. 
  We allowed different 
  phenotypes to potentially share common QTLs. For each phenotype, 
  the regression coefficients with other phenotypes are 
  chosen uniformly between 0.5 and 1; 
  QTL effects are chosen between 0.2 to 0.6; and residual standard deviations are chosen 
  from 0.1 to 0.5. For each realization we apply the QDG algorithm to infer causal 
  directions for the edges of the skeleton 
  obtained by the PC-skeleton algorithm.}
\details{For cyclic graphs, the output of the qdgAlgo function computes the 
log-likelihood 
up to the normalization constant (un-normalized log-likelihood). We can 
use the un-normalized 
log-likelihood to compare cyclic graphs with reversed directions (since they 
have the same normalization constant). However we cannot compare cyclic and 
acyclic graphs.
}
\seealso{
\code{\link[qtl]{sim.cross}}, 
\code{\link[qtl]{sim.geno}},
\code{\link[qtl]{sim.map}}, 
\code{\link[pcalg]{skeleton}},
\code{\link{qdgAlgo}},
\code{\link{graph.qdg}},
\code{\link{generate.qtl.pheno}}
}
\references{Chaibub Neto et al. (2008) Inferring causal phenotype networks from 
            segregating populations. Genetics 179: 1089-1100.}
\usage{data(acyclic)}
\examples{
\dontrun{
## This reproduces Figure 1 exactly.
set.seed(3456789)

tmp <- options(warn=-1)
acyclic.DG <- randomDAG(n = 100, prob = 2 / 99)

options(tmp)

## Simulate cross object using R/qtl routines.
n.ind <- 300
mymap <- sim.map(len=rep(100,20), n.mar=10, eq.spacing=FALSE, include.x=FALSE)
mycross <- sim.cross(map=mymap, n.ind=n.ind, type="f2")
summary(mycross)
mycross <- sim.geno(mycross,n.draws=1)


## Produce 100 QTL at three markers apiece.
acyclic.qtl <- generate.qtl.markers(cross=mycross,n.phe=100)

## Generate data from directed graph.
bp <- runif(100,0.5,1)
stdev <- runif(100,0.1,0.5)
bq <- matrix(0,100,3)
bq[,1] <- runif(100,0.2,0.4)
bq[,2] <- bq[,1]+0.1
bq[,3] <- bq[,2]+0.1
## Generate phenotypes.
acyclic.data <- generate.qtl.pheno("acyclic", cross = mycross,
  bp = bp, bq = bq, stdev = stdev, allqtl = acyclic.qtl$allqtl)

acyclic.qdg <- qdgAlgo(cross=acyclic.data, 
		phenotype.names=paste("y",1:100,sep=""),
		marker.names=acyclic.qtl$markers, 
		QTL=acyclic.qtl$allqtl, 
		alpha=0.005, 
		n.qdg.random.starts=1,
		skel.method="pcskel")
save(acyclic.DG, acyclic.qtl, acyclic.data, acyclic.qdg,
  file = "acyclic.RData", compress = TRUE)
}

data(acyclic)

dims <- dim(acyclic.data$pheno)
SuffStat <- list(C = cor(acyclic.data$pheno), n = dims[1])
pc <- skeleton(SuffStat, gaussCItest, p = dims[2], alpha = 0.005)
summary(pc)

summary(graph.qdg(acyclic.qdg))
gr <- graph.qdg(acyclic.qdg, include.qtl = FALSE)
plot(gr)
}
\keyword{datagen}
