## This creates object of class igraph.
## Is there a more general format?
graph.qtlnet <- function(qtlnet.object,
                         edges = summary(qtlnet.object)$averaged.net,
                         loci.list = loci.qtlnet(qtlnet.object, chr.pos),
                         pheno.color="green", qtl.color="red",
                         include.qtl=TRUE,
                         chr.pos = FALSE,
                         vertex.color = node.color,
                         ...)
{
  if(is.null(loci.list) | !include.qtl)
    node.color <- pheno.color
  else {
    loci.data.frame <- data.frame(qtl = unlist(loci.list))
    loci.data.frame$pheno <- rep(names(loci.list), sapply(loci.list, length))

    pheno.names <- attr(edges, "node.names")
    node.names <- c(pheno.names, levels(loci.data.frame[[1]]))

    edges <- cbind.data.frame(cause = c(as.character(edges[[1]]),
                                 as.character(loci.data.frame[[1]])),
                               react = c(as.character(edges[[2]]),
                                 as.character(loci.data.frame[[2]])),
                               width = c(edges[[3]],
                                 rep(1, nrow(loci.data.frame))))
    attr(edges, "node.names") <- node.names

    node.color <- rep(qtl.color, length(node.names))
    node.color[node.names %in% pheno.names] <- pheno.color
  }

  ## Set up vertices
  node.names <- attr(edges, "node.names")
  vertices <- data.frame(name = node.names, label = node.names, color = vertex.color)

  ## Great graph object (library igraph).
  graph.data.frame(edges, TRUE, vertices = vertices)
}
