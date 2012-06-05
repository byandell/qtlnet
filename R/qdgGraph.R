##############################################################################
##
## $Id: codeQDG.R,v 2007/11/28 byandell Exp $
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
## Routines: graph.qdg, plot.qdg
##############################################################################

graph.qdg <- function(x, ...) igraph.qdg(x, ...)

igraph.qdg <- function(x,
                       edges = myedges, loci.list = myloci.list, ...,
                       simple = FALSE)
{

  ## Prepare parameters for plotting function.
  if(inherits(x, "qdgAlgo")){ 
    best <- which(x$Solutions$BIC == min(x$Solutions$BIC))
    pheno.output <- data.frame(x$Solutions$solutions[[best]],rep(0,nrow(x$Solutions$solutions[[best]])))
  }
  else if (inherits(x, "qdg.sem")){ 
    best <- which(x$BIC.SEM[,1] == min(x$BIC.SEM[,1]))
    pheno.output <- data.frame(x$Solutions$solutions[[best]],x$path.coeffs)
  }
  names(pheno.output) <- c(names(x$Solutions$solutions[[best]]),"path")

  dir <- (pheno.output$direction == "---->")
  node1 <- as.character(pheno.output$node1)
  node2 <- as.character(pheno.output$node2)
  if(any(!dir)) {
    tmp <- node1[!dir]
    node1[!dir] <- node2[!dir]
    node2[!dir] <- tmp
  }
  loci <- x$phenotype.names
  myedges <- data.frame(cause = factor(node1, loci), effect = factor(node2, loci),
                      prob = pchisq(log(10) * pheno.output$lod, 1))
  myloci.list <- x$marker.names

  igraph.qtlnet(x, edges, loci.list, ...)
}
