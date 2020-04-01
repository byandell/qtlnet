\name{cyclicb}
\alias{cyclicb.data}
\alias{cyclicb.qtl}
\alias{cyclicb}
\title{Cyclic graph (b) example}
\description{We use a Gibbs sampling scheme to generate a data-set with
  200 individuals (according with cyclic graph (b)). Each phenotype is
  affected by 3 QTLs. We fixed the regression coefficients at 0.5, 
  error variances at 0.025 and the QTL effects at 0.2, 0.3 
  and 0.4 for the three F2 genotypes. We used 
  a burn-in of 2000 for the Gibbs sampler.}
\details{For cyclic graphs, the output of the qdg function
  computes the log-likelihood up to the normalization constant
  (un-normalized log-likelihood). We can use the un-normalized 
  log-likelihood to compare cyclic graphs with reversed directions 
  (since they have the same normalization constant). However we cannot 
  compare cyclic and acyclic graphs.}
\references{Chaibub Neto et al. (2008) Inferring causal phenotype networks from 
  segregating populations. Genetics 179: 1089-1100.}
\usage{data(cyclicb)}
\seealso{
\code{\link[qtl]{sim.cross}}, 
\code{\link[qtl]{sim.geno}},
\code{\link[qtl]{sim.map}}, 
\code{\link[pcalg]{skeleton}},
\code{\link{qdg}},
\code{\link{graph.qdg}},
\code{\link{generate.qtl.pheno}}
}
\examples{
\dontrun{
bp <- matrix(0, 6, 6)
bp[2,1] <- bp[1,5] <- bp[3,1] <- bp[4,2] <- bp[5,4] <- bp[5,6] <- bp[6,3] <- 0.5
stdev <- rep(0.025, 6)

## Use R/qtl routines to simulate.
set.seed(3456789)
mymap <- sim.map(len = rep(100,20), n.mar = 10, eq.spacing = FALSE,
  include.x = FALSE)
mycross <- sim.cross(map = mymap, n.ind = 200, type = "f2")
mycross <- sim.geno(mycross, n.draws = 1)

cyclicb.qtl <- generate.qtl.markers(cross = mycross, n.phe = 6)
mygeno <- pull.geno(mycross)[, unlist(cyclicb.qtl$markers)]

cyclicb.data <- generate.qtl.pheno("cyclicb", cross = mycross, burnin = 2000,
  bq = c(0.2,0.3,0.4), bp = bp, stdev = stdev, geno = mygeno)
save(cyclicb.qtl, cyclicb.data, file = "cyclicb.RData", compress = TRUE)

data(cyclicb)
out <- qdg(cross=cyclicb.data, 
		phenotype.names=paste("y",1:6,sep=""),
		marker.names=cyclicb.qtl$markers, 
		QTL=cyclicb.qtl$allqtl, 
		alpha=0.005, 
		n.qdg.random.starts=10,
		skel.method="pcskel")

gr <- graph.qdg(out)
gr
plot(gr)
}
}
\keyword{datagen}
