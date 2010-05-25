##########
## Call ##
##########
library(qtlnet)

## Load object mycross.
data(mycross)
summary(mycross)

newcross <- mycross
newcross$pheno <- mycross$pheno[,c(1,2,4,5,6)]

newcross <- calc.genoprob(newcross,step=1)

summary(newcross)

M <- matrix(0,5,5)

thr <- vector(mode="list", length=5)
thr[1:5] <- 3.83

random.seed <- 92387475

\dontrun{
  newcross.mcmc <- qtlnet.mcmc(M0=M, cross=newcross, thr=thr, 
                               addcov=NULL, intcov=NULL, nSamples=1000, thinning=20, 
                               pheno.col=c(1:5), random.seed = random.seed, verbose=TRUE)
  save(newcross.mcmc, file = "newcross.mcmc.RData", compress = TRUE)
}
data(newcross.mcmc)

## example of how to create addcov and intcov
##addcov <- vector(mode="list", length=10)
##addcov[1:10] <- "adipose.batch"
##intcov <- vector(mode="list", length=10)
##intcov[1:10] <- "Sex"

## Model average (average.net).
print(newcross.mcmc)

## averaged.net and averaged.posterior.table
summary(newcross.mcmc)

pns.qtl <- qtlnet.pheno(newcross.mcmc, cross=newcross, chr.pos=TRUE)
pns.qtl

##library(graph)
##library(Rgraphviz)
## OR
libary(igraph)

plot.qtlnet(newcross.mcmc, newcross, marker.list = pns.qtl, simple = FALSE,
            pheno.color="transparent", qtl.color="lightgrey", include.qtl=TRUE)


## Plot BIC over MCMC samples.
n <- length(newcross.mcmc[[2]])
## b is burnin
b <- floor(0.1*n)
xaxis <- seq(b,n,by=1)
plot(xaxis,newcross.mcmc[[2]][b:n],type="n",ylab="BIC",xlab="")
points(xaxis,newcross.mcmc[[2]][b:n], cex = 0.5)

library(lattice)
stripplot(newcross.mcmc$post.model[b:n] ~ newcross.mcmc$post.bic[b:n], horizontal = TRUE)

plot(as.numeric(table(newcross.mcmc$post.model[b:n])),
     unlist(tapply(newcross.mcmc$post.bic[b:n],newcross.mcmc$post.model[b:n],mean)),
     xlab = "model frequency", ylab = "model BIC")




















