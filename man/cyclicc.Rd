\name{cyclicc}
\alias{cyclicc.data}
\alias{cyclicc.qtl}
\alias{cyclicc}
\title{Cyclic graph (c) example}
\description{
  We use a Gibbs sampling scheme to generate a data-set with 200
  individuals (according with cyclic graph (c)). Each phenotype is
  affected by 3 QTLs. We fixed the regression coefficients at 0.5, (except
  for beta[5,2]=0.8) error variances at 0.025 and the QTL effects at 0.2,
  0.3 and 0.4 for the three F2 genotypes. We used a burn-in of 2000 for 
  the Gibbs sampler. This example illustrates that even though our method
  cannot detect reciprocal interactions (e.g. between phenotypes 2 and 5
  in cyclic graph (c)), it can still infer the stronger direction, that
  is, the direction corresponding to the higher regression
  coefficient. Since beta[5,2] is greater than beta[2,5], the QDG method
  should infer the direction from 2 to 5.
}
\details{For cyclic graphs, the output of the qdgAlgo function computes the 
log-likelihood up to the normalization constant (un-normalized log-likelihood). We can use the un-normalized log-likelihood to compare cyclic graphs with reversed directions (since they have the same normalization constant). However we 
cannot compare cyclic and acyclic graphs.}
\references{Chaibub Neto et al. (2008) Inferring causal phenotype networks from 
            segregating populations. Genetics 179: 1089-1100.}
\usage{data(cyclicc)}
\seealso{
\code{\link[qtl]{sim.cross}}, 
\code{\link[qtl]{sim.geno}},
\code{\link[qtl]{sim.map}}, 
\code{\link[pcalg]{skeleton}},
\code{\link{qdgAlgo}},
\code{\link{graph.qdg}},
\code{\link{generate.qtl.pheno}}
}
\examples{
\dontrun{
bp <- matrix(0, 6, 6)
bp[2,5] <- 0.5
bp[5,2] <- 0.8
bp[2,1] <- bp[3,2] <- bp[5,4] <- bp[6,5] <- 0.5
stdev <- rep(0.025, 6)

## Use R/qtl routines to simulate map and genotypes.
set.seed(34567899)
mymap <- sim.map(len = rep(100,20), n.mar = 10, eq.spacing = FALSE,
  include.x = FALSE)
mycross <- sim.cross(map = mymap, n.ind = 200, type = "f2")
mycross <- sim.geno(mycross, n.draws = 1)

## Use R/qdg routines to produce QTL sample and generate phenotypes.
cyclicc.qtl <- generate.qtl.markers(cross = mycross, n.phe = 6)
mygeno <- pull.geno(mycross)[, unlist(cyclicc.qtl$markers)]

cyclicc.data <- generate.qtl.pheno("cyclicc", cross = mycross, burnin = 2000,
  bq = c(0.2,0.3,0.4), bp = bp, stdev = stdev, geno = mygeno)
save(cyclicc.qtl, cyclicc.data, file = "cyclicc.RData", compress = TRUE)
}

data(cyclicc)
out <- qdgAlgo(cross=cyclicc.data, 
		phenotype.names=paste("y",1:6,sep=""),
		marker.names=cyclicc.qtl$markers, 
		QTL=cyclicc.qtl$allqtl, 
		alpha=0.005, 
		n.qdg.random.starts=1,
		skel.method="pcskel")

gr <- graph.qdg(out)
plot(gr)
}
\keyword{datagen}
