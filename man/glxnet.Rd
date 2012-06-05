\name{glxnet}
\alias{glxnet}
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
\examples{
data(glxnet)
glxnet <- calc.genoprob(glxnet)
set.seed(1234)
glxnet <- sim.geno(glxnet)

allqtls <- list()
markers <- list()
for(i in 1:8){
    age <- glxnet$pheno$age
    sex <- as.numeric(glxnet$pheno$sex == "M")
    ac <- glxnet$pheno[,c("age","sex")]
    ic <- glxnet$pheno[,c("sex")]
    tmp1 <- scanone(glxnet, pheno.col = i,
                    addcovar = data.frame(age,sex), 
                    intcovar = sex)
    tmp2 <- summary(tmp1, threshold = 2.999)
    chr <- as.numeric(as.vector(tmp2[[1]]))
    pos <- tmp2[[2]]
    allqtls[[i]] <- makeqtl(glxnet, chr = chr, pos = pos)
    markers[[i]] <- find.marker(glxnet, chr = chr, pos = pos)
}
names(allqtls) <- names(markers) <- names(glxnet$pheno[,1:8])

out <- qdgAlgo(cross=glxnet, 
		phenotype.names = names(glxnet$pheno[,1:8]), 
		marker.names = markers, 
		QTL = allqtls, 
		alpha = 0.05, 
		n.qdg.random.starts=10, 
		addcov="age", 
		intcov="sex", 
		skel.method="udgskel",
		udg.order=6)
out
gr <- graph.qdg(out)
plot(gr)
## Or use tkplot().
}
\keyword{datagen}
