\name{glxnet}
\alias{glxnet}
\alias{glxnet.cross}
\alias{glxnet.qtl}
\alias{glxnet.qdg}
\title{Generate and graph Glx network}

\references{Chaibub Neto et al. 2008 Inferring causal phenotype networks 
from segregating populations. Genetics 179: 1089-1100.

Ferrara et al. 2008 Genetic networks of liver metabolism revealed by 
integration of metabolomic and transcriptomic profiling. PLoS Genetics 4:
e1000034. 
}
\description{
This is the Glx network reported in Chaibub Neto et al 2008 and in Ferrara et
al 2008. Age was used as an additive covariate and we allowed for sex by 
genotype interaction. The network differs slightly from the published network 
due to improved code.
}
\seealso{\code{\link{qdg}}}
\examples{
data(glxnet)
glxnet.cross <- calc.genoprob(glxnet.cross)
set.seed(1234)
glxnet.cross <- sim.geno(glxnet.cross)

n.node <- nphe(glxnet.cross) - 2 ## Last two are age and sex.
markers <- glxnet.qtl <- vector("list", n.node)
for(i in 1:n.node) {
    ac <- model.matrix(~ age + sex, glxnet.cross$pheno)[, -1]
    ss <- summary(scanone(glxnet.cross, pheno.col = i,
                          addcovar = ac, intcovar = ac[,2]),
                  threshold = 2.999)
    glxnet.qtl[[i]] <- makeqtl(glxnet.cross, chr = ss$chr, pos = ss$pos)
    markers[[i]] <- find.marker(glxnet.cross, chr = ss$chr, pos = ss$pos)
}
names(glxnet.qtl) <- names(markers) <- names(glxnet.cross$pheno)[seq(n.node)]

glxnet.qdg <- qdg(cross=glxnet.cross, 
		phenotype.names = names(glxnet.cross$pheno[,seq(n.node)]), 
		marker.names = markers, 
		QTL = glxnet.qtl, 
		alpha = 0.05, 
		n.qdg.random.starts=10, 
		addcov="age", 
		intcov="sex", 
		skel.method="udgskel",
		udg.order=6)
glxnet.qdg
gr <- graph.qdg(glxnet.qdg)
plot(gr)
## Or use tkplot().
\dontrun{
glxnet.cross <- clean(glxnet.cross)
save(glxnet.cross, glxnet.qdg, glxnet.qtl, file = "glxnet.RData", compress = TRUE)
}
}
\keyword{datagen}
