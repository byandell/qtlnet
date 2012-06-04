\name{qdgSEM}
\alias{qdgSEM}
\alias{summary.qdgSEM}
\alias{print.qdgSEM}
\title{Score directed graphs outputed by qdgAlgo using structural equation models (SEM)}
\description{
  Score directed graphs (cyclic or acyclic) outputed by qdgAlgo function using the sem R
  package.     
}
\usage{
qdgSEM(qdgAlgoObject, cross)
\method{print}{qdgSEM}(x, \dots)
\method{summary}{qdgSEM}(object, \dots)
}
\arguments{
  \item{qdgAlgoObject}{list containing the output of \code{\link{qdgAlgo}}.}
  \item{cross}{object of class \code{cross} (see \code{\link[qtl]{read.cross}}).}
  \item{x,object}{object of class \code{qdgAlgo}.}
  \item{\dots}{extra arguments to print or summary (ignored).}
}
\details{
  Fits a SEM to the phenotypes network. QTLs are not included as
  variables in the model. When additive covariates are used in qdgAlgo, qdgSEM
  fits a SEM model to the residuals of the variables after adjustment of the 
  additive covariates.
}
\value{List object that inherits class "qdgSEM" and "qdg" composed by:
  \item{best.SEM}{Solution with lowest SEM BIC (best fit to the data).}
  \item{BIC.SEM}{Vector with the BIC values of all solutions from qdgAlgo.}
  \item{path.coeffs}{Path coefficients associated with the best SEM solution.}
  \item{Solutions}{Solutions of dependency graph after recheck step
    (output of steps 3, 4 and 5 of the QDG algorithm.)}
  \item{marker.names}{List of character strings, one for each of
    \code{phenotype.names}. 
    Each character string has the marker names
    for that phenotype.}
  \item{phenotype.names}{Character string with names of phenotype nodes
    corresponding to phenotypes in \code{cross}.}
  \item{dropped}{Indexes of solutions that were dropped (\code{NULL} if
    none dropped).}
}


\seealso{ \code{\link[sem]{sem}} }


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
y3 <- runif(1,0.5,1) * y1 +  runif(1,0.5,1) * y2 + runif(3,0.5,1)[g31] + runif(3,0.5,1)[g32] + rnorm(n.ind)
y4 <- runif(1,0.5,1) * y3 + runif(3,0.5,1)[g41] + runif(3,0.5,1)[g42] + rnorm(n.ind)

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
out <- qdgAlgo(cross=mycross, 
		phenotype.names = c("y1","y2","y3","y4"), 
		marker.names = markers, 
		QTL = allqtls, 
		alpha = 0.005, 
		n.qdg.random.starts=10, 
		skel.method="pcskel")


out2 <- qdgSEM(out, cross=mycross)
out2
plot(out2)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
