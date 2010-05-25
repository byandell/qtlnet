library(qtlnet)

## How many scanone runs might be needed?
n.pheno <- 13
max.parents <- 12
for(i in 3:max.parents)
  print(c(i, dim(make.saved.scores(1:n.pheno,i))))
## Guess: maximum feasible size is ~50K.

##############################################################
                                               ## pa size min
load("~/Rlib/qtlnet/data/mycross.mcmc.RData")  ##  3  299  16
load("~/Rlib/qtlnet/data/mycross4.mcmc.RData") ##  4  794  15*
load("~/Rlib/qtlnet/data/mycross5.mcmc.RData") ##  5 1586  28
load("~/Rlib/qtlnet/data/mycross6.mcmc.RData") ##  6 2510  24*
                                               ##  7 3302  ??
load("~/Rlib/qtlnet/data/mycrossall.mcmc.RData") ## 12 4096  26*

g3 <- graph.qtlnet(mycross.mcmc)
g4 <- graph.qtlnet(mycross4.mcmc)
g5 <- graph.qtlnet(mycross5.mcmc)
g6 <- graph.qtlnet(mycross6.mcmc)
g12 <- graph.qtlnet(mycrossall.mcmc)

tkplot(g3)
tkplot(g4)
tkplot(g5)
tkplot(g6)
tkplot(g12)

## What is complexity of graphs? Here are counts of colliders.
table(table(g3[[4]][g3[[3]] < 13]))
table(table(g4[[4]][g4[[3]] < 13]))
table(table(g5[[4]][g5[[3]] < 13]))
table(table(g6[[4]][g6[[3]] < 13]))
table(table(g12[[4]][g12[[3]] < 13]))
######################################################################
tmp <- make.saved.scores(attr(mycross6.mcmc, "pheno.names"), 12, mycross6.mcmc$saved.scores)
tmp3 <- make.saved.scores(attr(mycross6.mcmc, "pheno.names"), 12, mycross.mcmc$saved.scores)
tmp1 <- is.na(tmp)
tmp[tmp1] <- tmp3[tmp1]
tmp4 <- make.saved.scores(attr(mycross6.mcmc, "pheno.names"), 12, mycross4.mcmc$saved.scores)
tmp1 <- is.na(tmp)
tmp[tmp1] <- tmp4[tmp1]
tmp5 <- make.saved.scores(attr(mycross6.mcmc, "pheno.names"), 12, mycross5.mcmc$saved.scores)
tmp1 <- is.na(tmp)
tmp[tmp1] <- tmp5[tmp1]
