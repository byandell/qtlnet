#####################################################################
##
## $Id: generate.R,v 2007/11/28 byandell Exp $
##
##     Copyright (C) 2007 Elias Chaibub Neto and Brian S. Yandell
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the
## Free Software Foundation; either version 2, or (at your option) any
## later version.
##
## These functions are distributed in the hope that they will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## The text of the GNU General Public License, version 2, is available
## as http://www.gnu.org/copyleft or by writing to the Free Software
## Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
##
## Routines: generate.qtl.markers, generate.qtl.pheno
##############################################################################
#########################################
#########################################
## Generating data for some DAGs examples
## This code is specific for particular graphs, and is not meant
## to provide general tools for generating directed graphs.
#########################################
#########################################

#########################################
generate.qtl.markers <- function(cross, n.phe, nqtl = 3)
{
  ## randomly selects 2 or 3 markers (per phenotype) 
  nqtl <- array(nqtl, n.phe)
  allqtls <- list()
  markers <- list()
  for(i in 1:n.phe){
    if(nqtl[i]==2){
      chrm <- sample(c(1:20), 2, replace = TRUE)
      position <- sample(c(1:10), 2, replace = FALSE)
      position <- c(cross$geno[[chrm[1]]]$map[position[1]],
                    cross$geno[[chrm[2]]]$map[position[2]])
    }
    else{
      chrm <- sample(c(1:20), 3, replace = TRUE)
      position <- sample(c(1:10), 3, replace = FALSE)
      position <- c(cross$geno[[chrm[1]]]$map[position[1]],
                    cross$geno[[chrm[2]]]$map[position[2]],
                    cross$geno[[chrm[3]]]$map[position[3]])	
    }
    allqtls[[i]] <- makeqtl(cross, chr = chrm, pos = position)
    markers[[i]] <- find.marker(cross, chr = chrm, pos = position)
  }
  names(allqtls) <- paste("y", 1:n.phe, sep = "")
  names(markers) <- paste("y", 1:n.phe, sep = "")
  list(allqtl = allqtls, markers = markers)
}

##################################################################
generate.qtl.pheno <- function(name = c("acyclic","acyc2or3","cyclica","cyclicb","cyclicc"),
                               cross,
                               bp, bq, stdev, allqtl,
                               burnin = 2000, geno)
{
  name <- match.arg(name)
  switch(name,
         acyclic = { generate.data(cross, bp, bq, stdev, allqtl) },
         acyc2or3 = { generate.data.2or3(cross, bp, bq, stdev, allqtl) },
         cyclica = { generate.data.graph.a(cross, burnin, bq, bp, stdev, geno) },
         cyclicb = { generate.data.graph.b(cross, burnin, bq, bp, stdev, geno) },
         cyclicc = { generate.data.graph.c(cross, burnin, bq, bp, stdev, geno) })
}
         
##################################################################
## Acyclic example (100 phenotypes network)
##
generate.data <- function(cross,bp,bq,stdev,allqtl)
{
  n <- length(cross$pheno[,1])
  y <- matrix(0,n,100)
  y[,6] <- bq[6,allqtl[[6]]$geno[,1,]]+bq[6,allqtl[[6]]$geno[,2,]]+bq[6,allqtl[[6]]$geno[,3,]]+rnorm(n,0,stdev[6])
  y[,10] <- bq[10,allqtl[[10]]$geno[,1,]]+bq[10,allqtl[[10]]$geno[,2,]]+bq[10,allqtl[[10]]$geno[,3,]]+rnorm(n,0,stdev[10])
  y[,21] <- bq[21,allqtl[[21]]$geno[,1,]]+bq[21,allqtl[[21]]$geno[,2,]]+bq[21,allqtl[[21]]$geno[,3,]]+rnorm(n,0,stdev[21])
  y[,22] <- bq[22,allqtl[[22]]$geno[,1,]]+bq[22,allqtl[[22]]$geno[,2,]]+bq[22,allqtl[[22]]$geno[,3,]]+rnorm(n,0,stdev[22])
  y[,24] <- bq[24,allqtl[[24]]$geno[,1,]]+bq[24,allqtl[[24]]$geno[,2,]]+bq[24,allqtl[[24]]$geno[,3,]]+rnorm(n,0,stdev[24])
  y[,27] <- bq[27,allqtl[[27]]$geno[,1,]]+bq[27,allqtl[[27]]$geno[,2,]]+bq[27,allqtl[[27]]$geno[,3,]]+rnorm(n,0,stdev[27])
  y[,38] <- bq[38,allqtl[[38]]$geno[,1,]]+bq[38,allqtl[[38]]$geno[,2,]]+bq[38,allqtl[[38]]$geno[,3,]]+rnorm(n,0,stdev[38])
  y[,48] <- bq[48,allqtl[[48]]$geno[,1,]]+bq[48,allqtl[[48]]$geno[,2,]]+bq[48,allqtl[[48]]$geno[,3,]]+rnorm(n,0,stdev[48])
  y[,89] <- bq[89,allqtl[[89]]$geno[,1,]]+bq[89,allqtl[[89]]$geno[,2,]]+bq[89,allqtl[[89]]$geno[,3,]]+rnorm(n,0,stdev[89])
  y[,98] <- bq[98,allqtl[[98]]$geno[,1,]]+bq[98,allqtl[[98]]$geno[,2,]]+bq[98,allqtl[[98]]$geno[,3,]]+rnorm(n,0,stdev[98])
  y[,100] <- bq[100,allqtl[[100]]$geno[,1,]]+bq[100,allqtl[[100]]$geno[,2,]]+bq[100,allqtl[[100]]$geno[,3,]]+rnorm(n,0,stdev[100])
  y[,1] <- bq[1,allqtl[[1]]$geno[,1,]]+bq[1,allqtl[[1]]$geno[,2,]]+bq[1,allqtl[[1]]$geno[,3,]]+rnorm(n,0,stdev[1])
  y[,16] <- bq[16,allqtl[[16]]$geno[,1,]]+bq[16,allqtl[[16]]$geno[,2,]]+bq[16,allqtl[[16]]$geno[,3,]]+rnorm(n,0,stdev[16])
  y[,53] <- bq[53,allqtl[[53]]$geno[,1,]]+bq[53,allqtl[[53]]$geno[,2,]]+bq[53,allqtl[[53]]$geno[,3,]]+rnorm(n,0,stdev[53])
  y[,43] <- bq[43,allqtl[[43]]$geno[,1,]]+bq[43,allqtl[[43]]$geno[,2,]]+bq[43,allqtl[[43]]$geno[,3,]]+rnorm(n,0,stdev[43])
  y[,11] <- bq[11,allqtl[[11]]$geno[,1,]]+bq[11,allqtl[[11]]$geno[,2,]]+bq[11,allqtl[[11]]$geno[,3,]]+rnorm(n,0,stdev[11])
  y[,8] <- bq[8,allqtl[[8]]$geno[,1,]]+bq[8,allqtl[[8]]$geno[,2,]]+bq[8,allqtl[[8]]$geno[,3,]]+rnorm(n,0,stdev[8])
  y[,18] <- bq[18,allqtl[[18]]$geno[,1,]]+bq[18,allqtl[[18]]$geno[,2,]]+bq[18,allqtl[[18]]$geno[,3,]]+rnorm(n,0,stdev[18])
  y[,26] <- bq[26,allqtl[[26]]$geno[,1,]]+bq[26,allqtl[[26]]$geno[,2,]]+bq[26,allqtl[[26]]$geno[,3,]]+rnorm(n,0,stdev[26])
  y[,14] <- bq[14,allqtl[[14]]$geno[,1,]]+bq[14,allqtl[[14]]$geno[,2,]]+bq[14,allqtl[[14]]$geno[,3,]]+rnorm(n,0,stdev[14])
  y[,80] <- bq[80,allqtl[[80]]$geno[,1,]]+bq[80,allqtl[[80]]$geno[,2,]]+bq[80,allqtl[[80]]$geno[,3,]]+rnorm(n,0,stdev[80])
  y[,3] <- bq[3,allqtl[[3]]$geno[,1,]]+bq[3,allqtl[[3]]$geno[,2,]]+bq[3,allqtl[[3]]$geno[,3,]]+rnorm(n,0,stdev[3])
  y[,2] <- bq[2,allqtl[[2]]$geno[,1,]]+bq[2,allqtl[[2]]$geno[,2,]]+bq[2,allqtl[[2]]$geno[,3,]]+rnorm(n,0,stdev[2])
  y[,4] <- bq[4,allqtl[[4]]$geno[,1,]]+bq[4,allqtl[[4]]$geno[,2,]]+bq[4,allqtl[[4]]$geno[,3,]]+rnorm(n,0,stdev[4])
  y[,85] <- bq[85,allqtl[[85]]$geno[,1,]]+bq[85,allqtl[[85]]$geno[,2,]]+bq[85,allqtl[[85]]$geno[,3,]]+rnorm(n,0,stdev[85])
  y[,50] <- bq[50,allqtl[[50]]$geno[,1,]]+bq[50,allqtl[[50]]$geno[,2,]]+bq[50,allqtl[[50]]$geno[,3,]]+rnorm(n,0,stdev[50])
  y[,9] <- bq[9,allqtl[[9]]$geno[,1,]]+bq[9,allqtl[[9]]$geno[,2,]]+bq[9,allqtl[[9]]$geno[,3,]]+rnorm(n,0,stdev[9])
  y[,44] <- bq[44,allqtl[[44]]$geno[,1,]]+bq[44,allqtl[[44]]$geno[,2,]]+bq[44,allqtl[[44]]$geno[,3,]]+rnorm(n,0,stdev[44])
  y[,57] <- bq[57,allqtl[[57]]$geno[,1,]]+bq[57,allqtl[[57]]$geno[,2,]]+bq[57,allqtl[[57]]$geno[,3,]]+rnorm(n,0,stdev[57])
  y[,7] <- bq[7,allqtl[[7]]$geno[,1,]]+bq[7,allqtl[[7]]$geno[,2,]]+bq[7,allqtl[[7]]$geno[,3,]]+rnorm(n,0,stdev[7])
  y[,49] <- bq[49,allqtl[[49]]$geno[,1,]]+bq[49,allqtl[[49]]$geno[,2,]]+bq[49,allqtl[[49]]$geno[,3,]]+rnorm(n,0,stdev[49])
  y[,64] <- bq[64,allqtl[[64]]$geno[,1,]]+bq[64,allqtl[[64]]$geno[,2,]]+bq[64,allqtl[[64]]$geno[,3,]]+rnorm(n,0,stdev[64])
  y[,5] <- bq[5,allqtl[[5]]$geno[,1,]]+bq[5,allqtl[[5]]$geno[,2,]]+bq[5,allqtl[[5]]$geno[,3,]]+rnorm(n,0,stdev[5])
  y[,32] <- bq[32,allqtl[[32]]$geno[,1,]]+bq[32,allqtl[[32]]$geno[,2,]]+bq[32,allqtl[[32]]$geno[,3,]]+rnorm(n,0,stdev[32])
  y[,47] <- bq[47,allqtl[[47]]$geno[,1,]]+bq[47,allqtl[[47]]$geno[,2,]]+bq[47,allqtl[[47]]$geno[,3,]]+rnorm(n,0,stdev[47])
  y[,59] <- bq[59,allqtl[[59]]$geno[,1,]]+bq[59,allqtl[[59]]$geno[,2,]]+bq[59,allqtl[[59]]$geno[,3,]]+rnorm(n,0,stdev[59])
  y[,61] <- bp[61]*y[,6]+bq[61,allqtl[[61]]$geno[,1,]]+bq[61,allqtl[[61]]$geno[,2,]]+bq[61,allqtl[[61]]$geno[,3,]]+rnorm(n,0,stdev[61])
  y[,13] <- bp[13]*y[,6]+bq[13,allqtl[[13]]$geno[,1,]]+bq[13,allqtl[[13]]$geno[,2,]]+bq[13,allqtl[[13]]$geno[,3,]]+rnorm(n,0,stdev[13])
  y[,34] <- bp[34]*y[,27]+bq[34,allqtl[[34]]$geno[,1,]]+bq[34,allqtl[[34]]$geno[,2,]]+bq[34,allqtl[[34]]$geno[,3,]]+rnorm(n,0,stdev[34])
  y[,87] <- bp[87]*y[,1]+bq[87,allqtl[[87]]$geno[,1,]]+bq[87,allqtl[[87]]$geno[,2,]]+bq[87,allqtl[[87]]$geno[,3,]]+rnorm(n,0,stdev[87])
  y[,31] <- bp[31]*y[,1]+bq[31,allqtl[[31]]$geno[,1,]]+bq[31,allqtl[[31]]$geno[,2,]]+bq[31,allqtl[[31]]$geno[,3,]]+rnorm(n,0,stdev[31])
  y[,29] <- bp[29]*y[,16]+bq[29,allqtl[[29]]$geno[,1,]]+bq[29,allqtl[[29]]$geno[,2,]]+bq[29,allqtl[[29]]$geno[,3,]]+rnorm(n,0,stdev[29])
  y[,23] <- bp[23]*y[,16]+bq[23,allqtl[[23]]$geno[,1,]]+bq[23,allqtl[[23]]$geno[,2,]]+bq[23,allqtl[[23]]$geno[,3,]]+rnorm(n,0,stdev[23])
  y[,70] <- bp[70]*y[,1]+bp[70]*y[,53]+bp[70]*y[,43]+bq[70,allqtl[[70]]$geno[,1,]]+bq[70,allqtl[[70]]$geno[,2,]]+bq[70,allqtl[[70]]$geno[,3,]]+rnorm(n,0,stdev[70])
  y[,19] <- bp[19]*y[,11]+bq[19,allqtl[[19]]$geno[,1,]]+bq[19,allqtl[[19]]$geno[,2,]]+bq[19,allqtl[[19]]$geno[,3,]]+rnorm(n,0,stdev[19])
  y[,55] <- bp[55]*y[,11]+bp[55]*y[,13]+bq[55,allqtl[[55]]$geno[,1,]]+bq[55,allqtl[[55]]$geno[,2,]]+bq[55,allqtl[[55]]$geno[,3,]]+rnorm(n,0,stdev[55])
  y[,60] <- bp[60]*y[,34]+bq[60,allqtl[[60]]$geno[,1,]]+bq[60,allqtl[[60]]$geno[,2,]]+bq[60,allqtl[[60]]$geno[,3,]]+rnorm(n,0,stdev[60])
  y[,17] <- bp[17]*y[,8]+bq[17,allqtl[[17]]$geno[,1,]]+bq[17,allqtl[[17]]$geno[,2,]]+bq[17,allqtl[[17]]$geno[,3,]]+rnorm(n,0,stdev[17])
  y[,88] <- bp[88]*y[,8]+bq[88,allqtl[[88]]$geno[,1,]]+bq[88,allqtl[[88]]$geno[,2,]]+bq[88,allqtl[[88]]$geno[,3,]]+rnorm(n,0,stdev[88])
  y[,81] <- bp[81]*y[,8]+bp[81]*y[,31]+bq[81,allqtl[[81]]$geno[,1,]]+bq[81,allqtl[[81]]$geno[,2,]]+bq[81,allqtl[[81]]$geno[,3,]]+rnorm(n,0,stdev[81])
  y[,40] <- bp[40]*y[,29]+bq[40,allqtl[[40]]$geno[,1,]]+bq[40,allqtl[[40]]$geno[,2,]]+bq[40,allqtl[[40]]$geno[,3,]]+rnorm(n,0,stdev[40])
  y[,28] <- bp[28]*y[,23]+bq[28,allqtl[[28]]$geno[,1,]]+bq[28,allqtl[[28]]$geno[,2,]]+bq[28,allqtl[[28]]$geno[,3,]]+rnorm(n,0,stdev[28])
  y[,75] <- bp[75]*y[,16]+bp[75]*y[,70]+bp[75]*y[,11]+bq[75,allqtl[[75]]$geno[,1,]]+bq[75,allqtl[[75]]$geno[,2,]]+bq[75,allqtl[[75]]$geno[,3,]]+rnorm(n,0,stdev[75])
  y[,65] <- bp[65]*y[,29]+bp[65]*y[,18]+bq[65,allqtl[[65]]$geno[,1,]]+bq[65,allqtl[[65]]$geno[,2,]]+bq[65,allqtl[[65]]$geno[,3,]]+rnorm(n,0,stdev[65])
  y[,25] <- bp[25]*y[,18]+bp[25]*y[,19]+bq[25,allqtl[[25]]$geno[,1,]]+bq[25,allqtl[[25]]$geno[,2,]]+bq[25,allqtl[[25]]$geno[,3,]]+rnorm(n,0,stdev[25])
  y[,46] <- bp[46]*y[,19]+bq[46,allqtl[[46]]$geno[,1,]]+bq[46,allqtl[[46]]$geno[,2,]]+bq[46,allqtl[[46]]$geno[,3,]]+rnorm(n,0,stdev[46])
  y[,33] <- bp[33]*y[,26]+bq[33,allqtl[[33]]$geno[,1,]]+bq[33,allqtl[[33]]$geno[,2,]]+bq[33,allqtl[[33]]$geno[,3,]]+rnorm(n,0,stdev[33])
  y[,37] <- bp[37]*y[,26]+bq[37,allqtl[[37]]$geno[,1,]]+bq[37,allqtl[[37]]$geno[,2,]]+bq[37,allqtl[[37]]$geno[,3,]]+rnorm(n,0,stdev[37])
  y[,56] <- bp[56]*y[,1]+bp[56]*y[,40]+bq[56,allqtl[[56]]$geno[,1,]]+bq[56,allqtl[[56]]$geno[,2,]]+bq[56,allqtl[[56]]$geno[,3,]]+rnorm(n,0,stdev[56])
  y[,30] <- bp[30]*y[,17]+bq[30,allqtl[[30]]$geno[,1,]]+bq[30,allqtl[[30]]$geno[,2,]]+bq[30,allqtl[[30]]$geno[,3,]]+rnorm(n,0,stdev[30])
  y[,91] <- bp[91]*y[,14]+bp[91]*y[,80]+bp[91]*y[,65]+bq[91,allqtl[[91]]$geno[,1,]]+bq[91,allqtl[[91]]$geno[,2,]]+bq[91,allqtl[[91]]$geno[,3,]]+rnorm(n,0,stdev[91])
  y[,39] <- bp[39]*y[,29]+bq[39,allqtl[[39]]$geno[,1,]]+bq[39,allqtl[[39]]$geno[,2,]]+bq[39,allqtl[[39]]$geno[,3,]]+rnorm(n,0,stdev[39])
  y[,71] <- bp[39]*y[,4]+bq[71,allqtl[[71]]$geno[,1,]]+bq[71,allqtl[[71]]$geno[,2,]]+bq[71,allqtl[[71]]$geno[,3,]]+rnorm(n,0,stdev[71])
  y[,15] <- bp[15]*y[,3]+bq[15,allqtl[[15]]$geno[,1,]]+bq[15,allqtl[[15]]$geno[,2,]]+bq[15,allqtl[[15]]$geno[,3,]]+rnorm(n,0,stdev[15])
  y[,63] <- bp[63]*y[,40]+bp[63]*y[,2]+bp[63]*y[,4]+bq[63,allqtl[[63]]$geno[,1,]]+bq[63,allqtl[[63]]$geno[,2,]]+bq[63,allqtl[[63]]$geno[,3,]]+rnorm(n,0,stdev[63])
  y[,94] <- bp[94]*y[,65]+bp[94]*y[,4]+bp[94]*y[,85]+bp[94]*y[,19]+bq[94,allqtl[[94]]$geno[,1,]]+bq[94,allqtl[[94]]$geno[,2,]]+bq[94,allqtl[[94]]$geno[,3,]]+rnorm(n,0,stdev[94])
  y[,54] <- bp[54]*y[,50]+bq[54,allqtl[[54]]$geno[,1,]]+bq[54,allqtl[[54]]$geno[,2,]]+bq[54,allqtl[[54]]$geno[,3,]]+rnorm(n,0,stdev[54])
  y[,74] <- bp[74]*y[,50]+bq[74,allqtl[[74]]$geno[,1,]]+bq[74,allqtl[[74]]$geno[,2,]]+bq[74,allqtl[[74]]$geno[,3,]]+rnorm(n,0,stdev[74])
  y[,51] <- bp[51]*y[,18]+bp[51]*y[,46]+bq[51,allqtl[[51]]$geno[,1,]]+bq[51,allqtl[[51]]$geno[,2,]]+bq[51,allqtl[[51]]$geno[,3,]]+rnorm(n,0,stdev[51])
  y[,68] <- bp[68]*y[,46]+bq[68,allqtl[[68]]$geno[,1,]]+bq[68,allqtl[[68]]$geno[,2,]]+bq[68,allqtl[[68]]$geno[,3,]]+rnorm(n,0,stdev[68])
  y[,82] <- bp[82]*y[,4]+bp[82]*y[,9]+bp[82]*y[,44]+bq[82,allqtl[[82]]$geno[,1,]]+bq[82,allqtl[[82]]$geno[,2,]]+bq[82,allqtl[[82]]$geno[,3,]]+rnorm(n,0,stdev[82])
  y[,45] <- bp[45]*y[,44]+bq[45,allqtl[[45]]$geno[,1,]]+bq[45,allqtl[[45]]$geno[,2,]]+bq[45,allqtl[[45]]$geno[,3,]]+rnorm(n,0,stdev[45])
  y[,42] <- bp[42]*y[,33]+bq[42,allqtl[[42]]$geno[,1,]]+bq[42,allqtl[[42]]$geno[,2,]]+bq[42,allqtl[[42]]$geno[,3,]]+rnorm(n,0,stdev[42])
  y[,41] <- bp[41]*y[,37]+bq[41,allqtl[[41]]$geno[,1,]]+bq[41,allqtl[[41]]$geno[,2,]]+bq[41,allqtl[[41]]$geno[,3,]]+rnorm(n,0,stdev[41])
  y[,90] <- bp[90]*y[,57]+bp[90]*y[,56]+bq[90,allqtl[[90]]$geno[,1,]]+bq[90,allqtl[[90]]$geno[,2,]]+bq[90,allqtl[[90]]$geno[,3,]]+rnorm(n,0,stdev[90])
  y[,36] <- bp[36]*y[,30]+bq[36,allqtl[[36]]$geno[,1,]]+bq[36,allqtl[[36]]$geno[,2,]]+bq[36,allqtl[[36]]$geno[,3,]]+rnorm(n,0,stdev[36])
  y[,96] <- bp[96]*y[,7]+bq[96,allqtl[[96]]$geno[,1,]]+bq[96,allqtl[[96]]$geno[,2,]]+bq[96,allqtl[[96]]$geno[,3,]]+rnorm(n,0,stdev[96])
  y[,35] <- bp[35]*y[,7]+bq[35,allqtl[[35]]$geno[,1,]]+bq[35,allqtl[[35]]$geno[,2,]]+bq[35,allqtl[[35]]$geno[,3,]]+rnorm(n,0,stdev[35])
  y[,83] <- bp[83]*y[,30]+bp[83]*y[,7]+bp[83]*y[,39]+bp[83]*y[,71]+bp[83]*y[,49]+bp[83]*y[,54]+bq[83,allqtl[[83]]$geno[,1,]]+bq[83,allqtl[[83]]$geno[,2,]]+bq[83,allqtl[[83]]$geno[,3,]]+rnorm(n,0,stdev[83])
  y[,20] <- bp[20]*y[,3]+bq[20,allqtl[[20]]$geno[,1,]]+bq[20,allqtl[[20]]$geno[,2,]]+bq[20,allqtl[[20]]$geno[,3,]]+rnorm(n,0,stdev[20])
  y[,78] <- bp[78]*y[,15]+bp[78]*y[,43]+bp[78]*y[,51]+bq[78,allqtl[[78]]$geno[,1,]]+bq[78,allqtl[[78]]$geno[,2,]]+bq[78,allqtl[[78]]$geno[,3,]]+rnorm(n,0,stdev[78])
  y[,73] <- bp[73]*y[,63]+bq[73,allqtl[[73]]$geno[,1,]]+bq[73,allqtl[[73]]$geno[,2,]]+bq[73,allqtl[[73]]$geno[,3,]]+rnorm(n,0,stdev[73])
  y[,86] <- bp[86]*y[,63]+bp[86]*y[,64]+bq[86,allqtl[[86]]$geno[,1,]]+bq[86,allqtl[[86]]$geno[,2,]]+bq[86,allqtl[[86]]$geno[,3,]]+rnorm(n,0,stdev[86])
  y[,62] <- bp[62]*y[,51]+bq[62,allqtl[[62]]$geno[,1,]]+bq[62,allqtl[[62]]$geno[,2,]]+bq[62,allqtl[[62]]$geno[,3,]]+rnorm(n,0,stdev[62])
  y[,52] <- bp[52]*y[,49]+bq[52,allqtl[[52]]$geno[,1,]]+bq[52,allqtl[[52]]$geno[,2,]]+bq[52,allqtl[[52]]$geno[,3,]]+rnorm(n,0,stdev[52])
  y[,12] <- bp[12]*y[,7]+bp[12]*y[,5]+bq[12,allqtl[[12]]$geno[,1,]]+bq[12,allqtl[[12]]$geno[,2,]]+bq[12,allqtl[[12]]$geno[,3,]]+rnorm(n,0,stdev[12])
  y[,72] <- bp[72]*y[,5]+bq[72,allqtl[[72]]$geno[,1,]]+bq[72,allqtl[[72]]$geno[,2,]]+bq[72,allqtl[[72]]$geno[,3,]]+rnorm(n,0,stdev[72])
  y[,67] <- bp[67]*y[,54]+bp[67]*y[,19]+bq[67,allqtl[[67]]$geno[,1,]]+bq[67,allqtl[[67]]$geno[,2,]]+bq[67,allqtl[[67]]$geno[,3,]]+rnorm(n,0,stdev[67])
  y[,77] <- bp[77]*y[,45]+bq[77,allqtl[[77]]$geno[,1,]]+bq[77,allqtl[[77]]$geno[,2,]]+bq[77,allqtl[[77]]$geno[,3,]]+rnorm(n,0,stdev[77])
  y[,66] <- bp[66]*y[,36]+bq[66,allqtl[[66]]$geno[,1,]]+bq[66,allqtl[[66]]$geno[,2,]]+bq[66,allqtl[[66]]$geno[,3,]]+rnorm(n,0,stdev[66])
  y[,97] <- bp[97]*y[,41]+bp[97]*y[,40]+bp[97]*y[,96]+bp[97]*y[,83]+bq[97,allqtl[[97]]$geno[,1,]]+bq[97,allqtl[[97]]$geno[,2,]]+bq[97,allqtl[[97]]$geno[,3,]]+rnorm(n,0,stdev[97])
  y[,84] <- bp[84]*y[,40]+bp[84]*y[,20]+bp[84]*y[,73]+bq[84,allqtl[[84]]$geno[,1,]]+bq[84,allqtl[[84]]$geno[,2,]]+bq[84,allqtl[[84]]$geno[,3,]]+rnorm(n,0,stdev[84])
  y[,93] <- bp[93]*y[,39]+bp[93]*y[,78]+bq[93,allqtl[[93]]$geno[,1,]]+bq[93,allqtl[[93]]$geno[,2,]]+bq[93,allqtl[[93]]$geno[,3,]]+rnorm(n,0,stdev[93])
  y[,79] <- bp[79]*y[,32]+bp[79]*y[,62]+bq[79,allqtl[[79]]$geno[,1,]]+bq[79,allqtl[[79]]$geno[,2,]]+bq[79,allqtl[[79]]$geno[,3,]]+rnorm(n,0,stdev[79])
  y[,76] <- bp[76]*y[,47]+bp[76]*y[,52]+bq[76,allqtl[[76]]$geno[,1,]]+bq[76,allqtl[[76]]$geno[,2,]]+bq[76,allqtl[[76]]$geno[,3,]]+rnorm(n,0,stdev[76])
  y[,58] <- bp[58]*y[,12]+bq[58,allqtl[[58]]$geno[,1,]]+bq[58,allqtl[[58]]$geno[,2,]]+bq[58,allqtl[[58]]$geno[,3,]]+rnorm(n,0,stdev[58])
  y[,99] <- bp[99]*y[,4]+bp[99]*y[,79]+bq[99,allqtl[[99]]$geno[,1,]]+bq[99,allqtl[[99]]$geno[,2,]]+bq[99,allqtl[[99]]$geno[,3,]]+rnorm(n,0,stdev[99])
  y[,95] <- bp[95]*y[,76]+bp[95]*y[,50]+bq[95,allqtl[[95]]$geno[,1,]]+bq[95,allqtl[[95]]$geno[,2,]]+bq[95,allqtl[[95]]$geno[,3,]]+rnorm(n,0,stdev[95])
  y[,69] <- bp[69]*y[,58]+bp[69]*y[,50]+bp[69]*y[,43]+bp[69]*y[,59]+bq[69,allqtl[[69]]$geno[,1,]]+bq[69,allqtl[[69]]$geno[,2,]]+bq[69,allqtl[[69]]$geno[,3,]]+rnorm(n,0,stdev[69])
  y[,92] <- bp[92]*y[,59]+bq[92,allqtl[[92]]$geno[,1,]]+bq[92,allqtl[[92]]$geno[,2,]]+bq[92,allqtl[[92]]$geno[,3,]]+rnorm(n,0,stdev[92])

  y <- data.frame(y)
  names(y) <- paste("y",1:100,sep="")
  cross$pheno <- y
  return(cross)
}

###########################################################
## Actual example used in paper with two or three QTL per node.
##
generate.data.2or3 <- function(cross, bp, bq, stdev, allqtl)
{
  n <- length(cross$pheno[,1])
  y <- matrix(0,n,100)

  ## Yes, the following is clumsy, but the ordering is important
  ## in terms of the dependencies in the graph.
  ## This is set up for a specific graph used in the paper.
  
  y[,6] <- bq[6,allqtl[[6]]$geno[,1,]]+bq[6,allqtl[[6]]$geno[,2,]]+rnorm(n,0,stdev[6])

  y[,10] <- bq[10,allqtl[[10]]$geno[,1,]]+bq[10,allqtl[[10]]$geno[,2,]]+rnorm(n,0,stdev[10])

  y[,21] <- bq[21,allqtl[[21]]$geno[,1,]]+bq[21,allqtl[[21]]$geno[,2,]]+bq[21,allqtl[[21]]$geno[,3,]]+rnorm(n,0,stdev[21])

  y[,22] <- bq[22,allqtl[[22]]$geno[,1,]]+bq[22,allqtl[[22]]$geno[,2,]]+bq[22,allqtl[[22]]$geno[,3,]]+rnorm(n,0,stdev[22])

  y[,24] <- bq[24,allqtl[[24]]$geno[,1,]]+bq[24,allqtl[[24]]$geno[,2,]]+bq[24,allqtl[[24]]$geno[,3,]]+rnorm(n,0,stdev[24])

  y[,27] <- bq[27,allqtl[[27]]$geno[,1,]]+bq[27,allqtl[[27]]$geno[,2,]]+bq[27,allqtl[[27]]$geno[,3,]]+rnorm(n,0,stdev[27])

  y[,38] <- bq[38,allqtl[[38]]$geno[,1,]]+bq[38,allqtl[[38]]$geno[,2,]]+rnorm(n,0,stdev[38])

  y[,48] <- bq[48,allqtl[[48]]$geno[,1,]]+bq[48,allqtl[[48]]$geno[,2,]]+rnorm(n,0,stdev[48])

  y[,89] <- bq[89,allqtl[[89]]$geno[,1,]]+bq[89,allqtl[[89]]$geno[,2,]]+bq[89,allqtl[[89]]$geno[,3,]]+rnorm(n,0,stdev[89])

  y[,98] <- bq[98,allqtl[[98]]$geno[,1,]]+bq[98,allqtl[[98]]$geno[,2,]]+rnorm(n,0,stdev[98])

  y[,100] <- bq[100,allqtl[[100]]$geno[,1,]]+bq[100,allqtl[[100]]$geno[,2,]]+bq[100,allqtl[[100]]$geno[,3,]]+rnorm(n,0,stdev[100])

  y[,1] <- bq[1,allqtl[[1]]$geno[,1,]]+bq[1,allqtl[[1]]$geno[,2,]]+rnorm(n,0,stdev[1])

  y[,16] <- bq[16,allqtl[[16]]$geno[,1,]]+bq[16,allqtl[[16]]$geno[,2,]]+bq[16,allqtl[[16]]$geno[,3,]]+rnorm(n,0,stdev[16])

  y[,53] <- bq[53,allqtl[[53]]$geno[,1,]]+bq[53,allqtl[[53]]$geno[,2,]]+bq[53,allqtl[[53]]$geno[,3,]]+rnorm(n,0,stdev[53])

  y[,43] <- bq[43,allqtl[[43]]$geno[,1,]]+bq[43,allqtl[[43]]$geno[,2,]]+rnorm(n,0,stdev[43])

  y[,11] <- bq[11,allqtl[[11]]$geno[,1,]]+bq[11,allqtl[[11]]$geno[,2,]]+bq[11,allqtl[[11]]$geno[,3,]]+rnorm(n,0,stdev[11])

  y[,8] <- bq[8,allqtl[[8]]$geno[,1,]]+bq[8,allqtl[[8]]$geno[,2,]]+bq[8,allqtl[[8]]$geno[,3,]]+rnorm(n,0,stdev[8])

  y[,18] <- bq[18,allqtl[[18]]$geno[,1,]]+bq[18,allqtl[[18]]$geno[,2,]]+rnorm(n,0,stdev[18])

  y[,26] <- bq[26,allqtl[[26]]$geno[,1,]]+bq[26,allqtl[[26]]$geno[,2,]]+bq[26,allqtl[[26]]$geno[,3,]]+rnorm(n,0,stdev[26])

  y[,14] <- bq[14,allqtl[[14]]$geno[,1,]]+bq[14,allqtl[[14]]$geno[,2,]]+bq[14,allqtl[[14]]$geno[,3,]]+rnorm(n,0,stdev[14])

  y[,80] <- bq[80,allqtl[[80]]$geno[,1,]]+bq[80,allqtl[[80]]$geno[,2,]]+rnorm(n,0,stdev[80])

  y[,3] <- bq[3,allqtl[[3]]$geno[,1,]]+bq[3,allqtl[[3]]$geno[,2,]]+rnorm(n,0,stdev[3])

  y[,2] <- bq[2,allqtl[[2]]$geno[,1,]]+bq[2,allqtl[[2]]$geno[,2,]]+bq[2,allqtl[[2]]$geno[,3,]]+rnorm(n,0,stdev[2])

  y[,4] <- bq[4,allqtl[[4]]$geno[,1,]]+bq[4,allqtl[[4]]$geno[,2,]]+bq[4,allqtl[[4]]$geno[,3,]]+rnorm(n,0,stdev[4])

  y[,85] <- bq[85,allqtl[[85]]$geno[,1,]]+bq[85,allqtl[[85]]$geno[,2,]]+rnorm(n,0,stdev[85])

  y[,50] <- bq[50,allqtl[[50]]$geno[,1,]]+bq[50,allqtl[[50]]$geno[,2,]]+bq[50,allqtl[[50]]$geno[,3,]]+rnorm(n,0,stdev[50])

  y[,9] <- bq[9,allqtl[[9]]$geno[,1,]]+bq[9,allqtl[[9]]$geno[,2,]]+bq[9,allqtl[[9]]$geno[,3,]]+rnorm(n,0,stdev[9])

  y[,44] <- bq[44,allqtl[[44]]$geno[,1,]]+bq[44,allqtl[[44]]$geno[,2,]]+rnorm(n,0,stdev[44])

  y[,57] <- bq[57,allqtl[[57]]$geno[,1,]]+bq[57,allqtl[[57]]$geno[,2,]]+rnorm(n,0,stdev[57])

  y[,7] <- bq[7,allqtl[[7]]$geno[,1,]]+bq[7,allqtl[[7]]$geno[,2,]]+bq[7,allqtl[[7]]$geno[,3,]]+rnorm(n,0,stdev[7])

  y[,49] <- bq[49,allqtl[[49]]$geno[,1,]]+bq[49,allqtl[[49]]$geno[,2,]]+rnorm(n,0,stdev[49])

  y[,64] <- bq[64,allqtl[[64]]$geno[,1,]]+bq[64,allqtl[[64]]$geno[,2,]]+rnorm(n,0,stdev[64])

  y[,5] <- bq[5,allqtl[[5]]$geno[,1,]]+bq[5,allqtl[[5]]$geno[,2,]]+bq[5,allqtl[[5]]$geno[,3,]]+rnorm(n,0,stdev[5])

  y[,32] <- bq[32,allqtl[[32]]$geno[,1,]]+bq[32,allqtl[[32]]$geno[,2,]]+bq[32,allqtl[[32]]$geno[,3,]]+rnorm(n,0,stdev[32])

  y[,47] <- bq[47,allqtl[[47]]$geno[,1,]]+bq[47,allqtl[[47]]$geno[,2,]]+rnorm(n,0,stdev[47])

  y[,59] <- bq[59,allqtl[[59]]$geno[,1,]]+bq[59,allqtl[[59]]$geno[,2,]]+bq[59,allqtl[[59]]$geno[,3,]]+rnorm(n,0,stdev[59])

  y[,61] <- bp[61]*y[,6]+bq[61,allqtl[[61]]$geno[,1,]]+bq[61,allqtl[[61]]$geno[,2,]]+bq[61,allqtl[[61]]$geno[,3,]]+rnorm(n,0,stdev[61])

  y[,13] <- bp[13]*y[,6]+bq[13,allqtl[[13]]$geno[,1,]]+bq[13,allqtl[[13]]$geno[,2,]]+bq[13,allqtl[[13]]$geno[,3,]]+rnorm(n,0,stdev[13])

  y[,34] <- bp[34]*y[,27]+bq[34,allqtl[[34]]$geno[,1,]]+bq[34,allqtl[[34]]$geno[,2,]]+bq[34,allqtl[[34]]$geno[,3,]]+rnorm(n,0,stdev[34])

  y[,87] <- bp[87]*y[,1]+bq[87,allqtl[[87]]$geno[,1,]]+bq[87,allqtl[[87]]$geno[,2,]]+bq[87,allqtl[[87]]$geno[,3,]]+rnorm(n,0,stdev[87])

  y[,31] <- bp[31]*y[,1]+bq[31,allqtl[[31]]$geno[,1,]]+bq[31,allqtl[[31]]$geno[,2,]]+bq[31,allqtl[[31]]$geno[,3,]]+rnorm(n,0,stdev[31])

  y[,29] <- bp[29]*y[,16]+bq[29,allqtl[[29]]$geno[,1,]]+bq[29,allqtl[[29]]$geno[,2,]]+rnorm(n,0,stdev[29])

  y[,23] <- bp[23]*y[,16]+bq[23,allqtl[[23]]$geno[,1,]]+bq[23,allqtl[[23]]$geno[,2,]]+bq[23,allqtl[[23]]$geno[,3,]]+rnorm(n,0,stdev[23])

  y[,70] <- bp[70]*y[,1]+bp[70]*y[,53]+bp[70]*y[,43]+bq[70,allqtl[[70]]$geno[,1,]]+bq[70,allqtl[[70]]$geno[,2,]]+rnorm(n,0,stdev[70])

  y[,19] <- bp[19]*y[,11]+bq[19,allqtl[[19]]$geno[,1,]]+bq[19,allqtl[[19]]$geno[,2,]]+rnorm(n,0,stdev[19])

  y[,55] <- bp[55]*y[,11]+bp[55]*y[,13]+bq[55,allqtl[[55]]$geno[,1,]]+bq[55,allqtl[[55]]$geno[,2,]]+bq[55,allqtl[[55]]$geno[,3,]]+rnorm(n,0,stdev[55])

  y[,60] <- bp[60]*y[,34]+bq[60,allqtl[[60]]$geno[,1,]]+bq[60,allqtl[[60]]$geno[,2,]]+rnorm(n,0,stdev[60])

  y[,17] <- bp[17]*y[,8]+bq[17,allqtl[[17]]$geno[,1,]]+bq[17,allqtl[[17]]$geno[,2,]]+rnorm(n,0,stdev[17])

  y[,88] <- bp[88]*y[,8]+bq[88,allqtl[[88]]$geno[,1,]]+bq[88,allqtl[[88]]$geno[,2,]]+bq[88,allqtl[[88]]$geno[,3,]]+rnorm(n,0,stdev[88])

  y[,81] <- bp[81]*y[,8]+bp[81]*y[,31]+bq[81,allqtl[[81]]$geno[,1,]]+bq[81,allqtl[[81]]$geno[,2,]]+rnorm(n,0,stdev[81])

  y[,40] <- bp[40]*y[,29]+bq[40,allqtl[[40]]$geno[,1,]]+bq[40,allqtl[[40]]$geno[,2,]]+rnorm(n,0,stdev[40])

  y[,28] <- bp[28]*y[,23]+bq[28,allqtl[[28]]$geno[,1,]]+bq[28,allqtl[[28]]$geno[,2,]]+bq[28,allqtl[[28]]$geno[,3,]]+rnorm(n,0,stdev[28])

  y[,75] <- bp[75]*y[,16]+bp[75]*y[,70]+bp[75]*y[,11]+bq[75,allqtl[[75]]$geno[,1,]]+bq[75,allqtl[[75]]$geno[,2,]]+rnorm(n,0,stdev[75])

  y[,65] <- bp[65]*y[,29]+bp[65]*y[,18]+bq[65,allqtl[[65]]$geno[,1,]]+bq[65,allqtl[[65]]$geno[,2,]]+bq[65,allqtl[[65]]$geno[,3,]]+rnorm(n,0,stdev[65])

  y[,25] <- bp[25]*y[,18]+bp[25]*y[,19]+bq[25,allqtl[[25]]$geno[,1,]]+bq[25,allqtl[[25]]$geno[,2,]]+bq[25,allqtl[[25]]$geno[,3,]]+rnorm(n,0,stdev[25])

  y[,46] <- bp[46]*y[,19]+bq[46,allqtl[[46]]$geno[,1,]]+bq[46,allqtl[[46]]$geno[,2,]]+rnorm(n,0,stdev[46])

  y[,33] <- bp[33]*y[,26]+bq[33,allqtl[[33]]$geno[,1,]]+bq[33,allqtl[[33]]$geno[,2,]]+bq[33,allqtl[[33]]$geno[,3,]]+rnorm(n,0,stdev[33])

  y[,37] <- bp[37]*y[,26]+bq[37,allqtl[[37]]$geno[,1,]]+bq[37,allqtl[[37]]$geno[,2,]]+bq[37,allqtl[[37]]$geno[,3,]]+rnorm(n,0,stdev[37])

  y[,56] <- bp[56]*y[,1]+bp[56]*y[,40]+bq[56,allqtl[[56]]$geno[,1,]]+bq[56,allqtl[[56]]$geno[,2,]]+rnorm(n,0,stdev[56])

  y[,30] <- bp[30]*y[,17]+bq[30,allqtl[[30]]$geno[,1,]]+bq[30,allqtl[[30]]$geno[,2,]]+rnorm(n,0,stdev[30])

  y[,91] <- bp[91]*y[,14]+bp[91]*y[,80]+bp[91]*y[,65]+bq[91,allqtl[[91]]$geno[,1,]]+bq[91,allqtl[[91]]$geno[,2,]]+rnorm(n,0,stdev[91])

  y[,39] <- bp[39]*y[,29]+bq[39,allqtl[[39]]$geno[,1,]]+bq[39,allqtl[[39]]$geno[,2,]]+rnorm(n,0,stdev[39])

  y[,71] <- bp[39]*y[,4]+bq[71,allqtl[[71]]$geno[,1,]]+bq[71,allqtl[[71]]$geno[,2,]]+bq[71,allqtl[[71]]$geno[,3,]]+rnorm(n,0,stdev[71])

  y[,15] <- bp[15]*y[,3]+bq[15,allqtl[[15]]$geno[,1,]]+bq[15,allqtl[[15]]$geno[,2,]]+rnorm(n,0,stdev[15])

  y[,63] <- bp[63]*y[,40]+bp[63]*y[,2]+bp[63]*y[,4]+bq[63,allqtl[[63]]$geno[,1,]]+bq[63,allqtl[[63]]$geno[,2,]]+rnorm(n,0,stdev[63])

  y[,94] <- bp[94]*y[,65]+bp[94]*y[,4]+bp[94]*y[,85]+bp[94]*y[,19]+bq[94,allqtl[[94]]$geno[,1,]]+bq[94,allqtl[[94]]$geno[,2,]]+bq[94,allqtl[[94]]$geno[,3,]]+rnorm(n,0,stdev[94])

  y[,54] <- bp[54]*y[,50]+bq[54,allqtl[[54]]$geno[,1,]]+bq[54,allqtl[[54]]$geno[,2,]]+rnorm(n,0,stdev[54])

  y[,74] <- bp[74]*y[,50]+bq[74,allqtl[[74]]$geno[,1,]]+bq[74,allqtl[[74]]$geno[,2,]]+rnorm(n,0,stdev[74])

  y[,51] <- bp[51]*y[,18]+bp[51]*y[,46]+bq[51,allqtl[[51]]$geno[,1,]]+bq[51,allqtl[[51]]$geno[,2,]]+rnorm(n,0,stdev[51])

  y[,68] <- bp[68]*y[,46]+bq[68,allqtl[[68]]$geno[,1,]]+bq[68,allqtl[[68]]$geno[,2,]]+bq[68,allqtl[[68]]$geno[,3,]]+rnorm(n,0,stdev[68])

  y[,82] <- bp[82]*y[,4]+bp[82]*y[,9]+bp[82]*y[,44]+bq[82,allqtl[[82]]$geno[,1,]]+bq[82,allqtl[[82]]$geno[,2,]]+bq[82,allqtl[[82]]$geno[,3,]]+rnorm(n,0,stdev[82])

  y[,45] <- bp[45]*y[,44]+bq[45,allqtl[[45]]$geno[,1,]]+bq[45,allqtl[[45]]$geno[,2,]]+rnorm(n,0,stdev[45])

  y[,42] <- bp[42]*y[,33]+bq[42,allqtl[[42]]$geno[,1,]]+bq[42,allqtl[[42]]$geno[,2,]]+rnorm(n,0,stdev[42])

  y[,41] <- bp[41]*y[,37]+bq[41,allqtl[[41]]$geno[,1,]]+bq[41,allqtl[[41]]$geno[,2,]]+rnorm(n,0,stdev[41])

  y[,90] <- bp[90]*y[,57]+bp[90]*y[,56]+bq[90,allqtl[[90]]$geno[,1,]]+bq[90,allqtl[[90]]$geno[,2,]]+rnorm(n,0,stdev[90])

  y[,36] <- bp[36]*y[,30]+bq[36,allqtl[[36]]$geno[,1,]]+bq[36,allqtl[[36]]$geno[,2,]]+rnorm(n,0,stdev[36])

  y[,96] <- bp[96]*y[,7]+bq[96,allqtl[[96]]$geno[,1,]]+bq[96,allqtl[[96]]$geno[,2,]]+rnorm(n,0,stdev[96])

  y[,35] <- bp[35]*y[,7]+bq[35,allqtl[[35]]$geno[,1,]]+bq[35,allqtl[[35]]$geno[,2,]]+rnorm(n,0,stdev[35])

  y[,83] <- bp[83]*y[,30]+bp[83]*y[,7]+bp[83]*y[,39]+bp[83]*y[,71]+bp[83]*y[,49]+bp[83]*y[,54]+bq[83,allqtl[[83]]$geno[,1,]]+bq[83,allqtl[[83]]$geno[,2,]]+rnorm(n,0,stdev[83])

  y[,20] <- bp[20]*y[,3]+bq[20,allqtl[[20]]$geno[,1,]]+bq[20,allqtl[[20]]$geno[,2,]]+bq[20,allqtl[[20]]$geno[,3,]]+rnorm(n,0,stdev[20])

  y[,78] <- bp[78]*y[,15]+bp[78]*y[,43]+bp[78]*y[,51]+bq[78,allqtl[[78]]$geno[,1,]]+bq[78,allqtl[[78]]$geno[,2,]]+bq[78,allqtl[[78]]$geno[,3,]]+rnorm(n,0,stdev[78])

  y[,73] <- bp[73]*y[,63]+bq[73,allqtl[[73]]$geno[,1,]]+bq[73,allqtl[[73]]$geno[,2,]]+bq[73,allqtl[[73]]$geno[,3,]]+rnorm(n,0,stdev[73])

  y[,86] <- bp[86]*y[,63]+bp[86]*y[,64]+bq[86,allqtl[[86]]$geno[,1,]]+bq[86,allqtl[[86]]$geno[,2,]]+rnorm(n,0,stdev[86])

  y[,62] <- bp[62]*y[,51]+bq[62,allqtl[[62]]$geno[,1,]]+bq[62,allqtl[[62]]$geno[,2,]]+rnorm(n,0,stdev[62])

  y[,52] <- bp[52]*y[,49]+bq[52,allqtl[[52]]$geno[,1,]]+bq[52,allqtl[[52]]$geno[,2,]]+rnorm(n,0,stdev[52])

  y[,12] <- bp[12]*y[,7]+bp[12]*y[,5]+bq[12,allqtl[[12]]$geno[,1,]]+bq[12,allqtl[[12]]$geno[,2,]]+rnorm(n,0,stdev[12])

  y[,72] <- bp[72]*y[,5]+bq[72,allqtl[[72]]$geno[,1,]]+bq[72,allqtl[[72]]$geno[,2,]]+bq[72,allqtl[[72]]$geno[,3,]]+rnorm(n,0,stdev[72])

  y[,67] <- bp[67]*y[,54]+bp[67]*y[,19]+bq[67,allqtl[[67]]$geno[,1,]]+bq[67,allqtl[[67]]$geno[,2,]]+bq[67,allqtl[[67]]$geno[,3,]]+rnorm(n,0,stdev[67])

  y[,77] <- bp[77]*y[,45]+bq[77,allqtl[[77]]$geno[,1,]]+bq[77,allqtl[[77]]$geno[,2,]]+rnorm(n,0,stdev[77])

  y[,66] <- bp[66]*y[,36]+bq[66,allqtl[[66]]$geno[,1,]]+bq[66,allqtl[[66]]$geno[,2,]]+rnorm(n,0,stdev[66])

  y[,97] <- bp[97]*y[,41]+bp[97]*y[,40]+bp[97]*y[,96]+bp[97]*y[,83]+bq[97,allqtl[[97]]$geno[,1,]]+bq[97,allqtl[[97]]$geno[,2,]]+bq[97,allqtl[[97]]$geno[,3,]]+rnorm(n,0,stdev[97])

  y[,84] <- bp[84]*y[,40]+bp[84]*y[,20]+bp[84]*y[,73]+bq[84,allqtl[[84]]$geno[,1,]]+bq[84,allqtl[[84]]$geno[,2,]]+bq[84,allqtl[[84]]$geno[,3,]]+rnorm(n,0,stdev[84])

  y[,93] <- bp[93]*y[,39]+bp[93]*y[,78]+bq[93,allqtl[[93]]$geno[,1,]]+bq[93,allqtl[[93]]$geno[,2,]]+rnorm(n,0,stdev[93])

  y[,79] <- bp[79]*y[,32]+bp[79]*y[,62]+bq[79,allqtl[[79]]$geno[,1,]]+bq[79,allqtl[[79]]$geno[,2,]]+rnorm(n,0,stdev[79])

  y[,76] <- bp[76]*y[,47]+bp[76]*y[,52]+bq[76,allqtl[[76]]$geno[,1,]]+bq[76,allqtl[[76]]$geno[,2,]]+rnorm(n,0,stdev[76])

  y[,58] <- bp[58]*y[,12]+bq[58,allqtl[[58]]$geno[,1,]]+bq[58,allqtl[[58]]$geno[,2,]]+bq[58,allqtl[[58]]$geno[,3,]]+rnorm(n,0,stdev[58])

  y[,99] <- bp[99]*y[,4]+bp[99]*y[,79]+bq[99,allqtl[[99]]$geno[,1,]]+bq[99,allqtl[[99]]$geno[,2,]]+rnorm(n,0,stdev[99])

  y[,95] <- bp[95]*y[,76]+bp[95]*y[,50]+bq[95,allqtl[[95]]$geno[,1,]]+bq[95,allqtl[[95]]$geno[,2,]]+rnorm(n,0,stdev[95])

  y[,69] <- bp[69]*y[,58]+bp[69]*y[,50]+bp[69]*y[,43]+bp[69]*y[,59]+bq[69,allqtl[[69]]$geno[,1,]]+bq[69,allqtl[[69]]$geno[,2,]]+bq[69,allqtl[[69]]$geno[,3,]]+rnorm(n,0,stdev[69])

  y[,92] <- bp[92]*y[,59]+bq[92,allqtl[[92]]$geno[,1,]]+bq[92,allqtl[[92]]$geno[,2,]]+bq[92,allqtl[[92]]$geno[,3,]]+rnorm(n,0,stdev[92])

  y <- data.frame(y)
  names(y) <- paste("y",1:100,sep="")
  cross$pheno <- y
  cross
}

###############
###############
## cyclic graphs 
###############
###############

###################################################
## compute the means (that depend on the genotypes) 
## of each individual. These means are used in the
## generating of the phenotype data 
##
compute.mu <- function(ind.geno,bq){
  mu <- rep(0, 6)
  for(i in 1:6)
    mu[i] <- sum(bq[ind.geno[(3 * (i - 1)) + (1:3)]])
  mu
}

############################################################
## generate the phenotype data. Each data point is generated 
## by a separate Markov chain
##
generate.data.graph.a <- function(cross,burnin,bq,bp,stdev,geno)
{
  ## Gibbs sampler for graph (a) 
  gibbs.graph.a <- function(n,burnin,bp,stdev,mu){
    phi6 <- stdev[6]
    aux2 <- stdev[4]+stdev[2]*bp[4,2]^2
    phi2 <- stdev[2]*stdev[4]/aux2
    aux4 <- stdev[5]+stdev[4]*bp[5,4]^2
    phi4 <- stdev[4]*stdev[5]/aux4
    aux5 <- stdev[2]*stdev[6]+(bp[2,5]^2)*stdev[5]*stdev[6]+(bp[6,5]^2)*stdev[2]*stdev[5]
    phi5 <- stdev[2]*stdev[5]*stdev[6]/aux5
    aux1 <- stdev[2]+stdev[1]*bp[2,1]^2
    phi1 <- stdev[1]*stdev[2]/aux1
    aux3 <- stdev[4]+stdev[3]*bp[4,3]^2
    phi3 <- stdev[3]*stdev[4]/aux3
    y1 <- y2 <- y3 <- y4 <- y5 <- y6 <- rep(0,n+burnin)
    for(i in 2:(n+burnin)){
      m1 <- (stdev[2]*mu[1]+stdev[1]*bp[2,1]*(y2[i-1]-mu[2]-bp[2,5]*y5[i-1]))/aux1
      y1[i] <- rnorm(1,m1,sqrt(phi1))

      m2 <- (stdev[4]*(mu[2]+bp[2,1]*y1[i]+bp[2,5]*y5[i-1])+stdev[2]*bp[4,2]*(y4[i-1]-mu[4]-bp[4,3]*y3[i-1]))/aux2
      y2[i] <- rnorm(1,m2,sqrt(phi2))

      m3 <- (stdev[4]*mu[3]+stdev[3]*bp[4,2]*(y4[i-1]-mu[4]-bp[4,2]*y2[i]))/aux3
      y3[i] <- rnorm(1,m3,sqrt(phi3))

      m4 <- (stdev[5]*(mu[4]+bp[4,2]*y2[i]+bp[4,3]*y3[i])+stdev[4]*bp[5,4]*(y5[i-1]-mu[5]))/aux4
      y4[i] <- rnorm(1,m4,sqrt(phi4))

      m5 <- (stdev[2]*stdev[6]*(mu[5]+bp[5,4]*y4[i])+bp[2,5]*stdev[5]*stdev[6]*(y2[i]-mu[2]-bp[2,1]*mu[1])+
             bp[6,5]*stdev[2]*stdev[5]*(y6[i-1]-mu[6]))/aux5
      y5[i] <- rnorm(1,m5,sqrt(phi5))

      m6 <- mu[6]+bp[6,5]*y5[i]
      y6[i] <- rnorm(1,m6,sqrt(phi6))
    }	
    y <- cbind(y1,y2,y3,y4,y5,y6)
    return(y[-c(1:burnin),])
  }

  n <- length(cross$pheno[,1])
  y <- matrix(0,n,6)
  for(i in 1:n){
    mu <- compute.mu(ind.geno=geno[i,],bq=bq)
    y[i,] <- as.vector(gibbs.graph.a(n=1,burnin=burnin,bp,stdev,mu=mu))
  }
  y <- data.frame(y)
  names(y) <- c("y1","y2","y3","y4","y5","y6")
  cross$pheno <- y
  return(cross)
} 


###############################
## Cyclic toy example (graph b)
generate.data.graph.b <- function(cross,burnin,bq,bp,stdev,geno)
{
  gibbs.graph.b <- function(n,burnin,bp,stdev,mu){
    aux1 <- stdev[2]*stdev[3]+stdev[1]*stdev[3]*(bp[2,1]^1)+stdev[1]*stdev[2]*(bp[3,1]^2)
    phi1 <- stdev[1]*stdev[2]*stdev[3]/aux1
    aux2 <- stdev[4]+stdev[2]*bp[4,2]^2
    phi2 <- stdev[2]*stdev[4]/aux2
    aux3 <- stdev[6]+stdev[3]*bp[6,3]^2
    phi3 <- stdev[3]*stdev[6]/aux3
    aux4 <- stdev[5]+stdev[4]*bp[5,4]^2
    phi4 <- stdev[4]*stdev[5]/aux4
    aux5 <- stdev[1]+stdev[5]*bp[1,5]^2
    phi5 <- stdev[1]*stdev[5]/aux5
    aux6 <- stdev[5]+stdev[6]*bp[5,6]^2
    phi6 <- stdev[5]*stdev[6]/aux6
    y1 <- y2 <- y3 <- y4 <- y5 <- y6 <- rep(0,n+burnin)
    for(i in 2:(n+burnin)){
      m1 <- (stdev[2]*stdev[3]*(mu[1]+bp[1,5]*y5[i-1])+stdev[1]*stdev[3]*bp[2,1]*(y2[i-1]-mu[2])+stdev[1]*stdev[2]*bp[3,1]*(y3[i-1]-mu[3]))/aux1
      y1[i] <- rnorm(1,m1, sqrt(phi1))

      m2 <- (stdev[4]*(mu[2]+bp[2,1]*y1[i])+stdev[2]*bp[4,2]*(y4[i-1]-mu[4]))/aux2
      y2[i] <- rnorm(1, m2, sqrt(phi2))

      m3 <- (stdev[6]*(mu[3]+bp[3,1]*y1[i])+stdev[3]*bp[6,3]*(y6[i-1]-mu[6]))/aux3
      y3[i] <- rnorm(1, m3, sqrt(phi3))

      m4 <- (stdev[5]*(mu[4]+bp[4,2]*y2[i])+stdev[4]*bp[5,4]*(y5[i-1]-mu[5]-bp[5,6]*y6[i-1]))/aux4
      y4[i] <- rnorm(1, m4, sqrt(phi4))

      m5 <- (stdev[1]*(mu[5]+bp[5,4]*y4[i]+bp[5,6]*y6[i-1])+stdev[5]*bp[1,5]*(y1[i]-mu[1]))/aux5
      y5[i] <- rnorm(1, m5, sqrt(phi5))

      m6 <- (stdev[5]*(mu[6]+bp[6,3]*y3[i])+stdev[6]*bp[5,6]*(y5[i]-mu[5]-bp[5,4]*y4[i]))/aux6
      y6[i] <- rnorm(1, m6, sqrt(phi6))
    }	
    y <- cbind(y1,y2,y3,y4,y5,y6)
    return(y[-c(1:burnin),])
  }

  n <- length(cross$pheno[,1])
  y <- matrix(0,n,6)
  for(i in 1:n){
    mu <- compute.mu(ind.geno=geno[i,],bq=bq)
    y[i,] <- as.vector(gibbs.graph.b(n=1,burnin=burnin,bp,stdev,mu=mu))
  }
  y <- data.frame(y)
  names(y) <- c("y1","y2","y3","y4","y5","y6")
  cross$pheno <- y
  return(cross)
} 


###############################
## Cyclic toy example (graph c)
generate.data.graph.c <- function(cross,burnin,bq,bp,stdev,geno)
{
  gibbs.graph.c <- function(n,burnin,bp,stdev,mu){
    aux1 <- stdev[2]+stdev[1]*bp[2,1]^2
    phi1 <- stdev[1]*stdev[2]/aux1
    aux2 <- stdev[3]*stdev[5]+stdev[2]*stdev[5]*(bp[3,2]^2)+stdev[2]*stdev[3]*(bp[5,2]^2)
    phi2 <- stdev[2]*stdev[3]*stdev[5]/aux2
    phi3 <- stdev[3]
    aux4 <- stdev[5]+stdev[4]*bp[5,4]^2
    phi4 <- stdev[4]*stdev[6]/aux4
    aux5 <- stdev[2]*stdev[6]+stdev[5]*stdev[6]*(bp[2,5]^2)+stdev[2]*stdev[5]*(bp[6,5]^2)
    phi5 <- stdev[2]*stdev[5]*stdev[6]/aux5
    phi6 <- stdev[6]
    y1 <- y2 <- y3 <- y4 <- y5 <- y6 <- rep(0,n+burnin)
    for(i in 2:(n+burnin)){
      m1 <- (stdev[2]*mu[1]+stdev[1]*bp[2,1]*(y2[i-1]-mu[2]-bp[2,5]*y5[i-1]))/aux1
      y1[i] <- rnorm(1, m1, sqrt(phi1))
      
      m2 <- (stdev[3]*stdev[5]*(mu[2]+bp[2,1]*y1[i]+bp[2,5]*y5[i-1])+stdev[2]*stdev[5]*bp[3,2]*(y3[i-1]-mu[3])+stdev[2]*stdev[3]*bp[5,2]*(y5[i-1]-mu[5]-bp[5,4]*y4[i-1]))/aux2
      y2[i] <- rnorm(1, m2, sqrt(phi2))

      m3 <- mu[3]+bp[3,2]*y2[i]
      y3[i] <- rnorm(1, m3, sqrt(phi3))

      m4 <- (stdev[5]*mu[4]+stdev[4]*bp[5,4]*(y5[i-1]-mu[5]-bp[5,2]*y2[i]))/aux4
      y4[i] <- rnorm(1, m4, sqrt(phi4))

      m5 <- (stdev[2]*stdev[6]*(mu[5]+bp[5,4]*y4[i]+bp[5,2]*y2[i])+stdev[5]*stdev[6]*bp[2,5]*(y2[i]-mu[2]-bp[2,1]*y1[i])+stdev[2]*stdev[5]*bp[6,5]*(y6[i-1]-mu[6]))/aux5
      y5[i] <- rnorm(1, m5, sqrt(phi5))

      m6 <- mu[6]+bp[6,5]*y5[i]
      y6[i] <- rnorm(1, m6, sqrt(phi6))
    }	
    y <- cbind(y1,y2,y3,y4,y5,y6)
    return(y[-c(1:burnin),])
  }

  n <- length(cross$pheno[,1])
  y <- matrix(0,n,6)
  for(i in 1:n){
    mu <- compute.mu(ind.geno=geno[i,],bq=bq)
    tmp <- gibbs.graph.c(n=1, burnin, bp, stdev, mu)
    y[i,] <- as.vector(gibbs.graph.c(n=1, burnin, bp, stdev, mu))
  }
  y <- data.frame(y)
  names(y) <- paste("y", 1:6, sep = "")
  cross$pheno <- y
  return(cross)
} 
