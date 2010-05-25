plot.qtlnet <- function(x, ...)
{
  gr <- graph.qtlnet(x, ...)
  tkplot(gr, ...)
  
  invisible(gr)
}
###################################################################
## This creates object of class igraph.
## Is there a more general format?
graph.qtlnet <- function(x,
                         edges = summary(x, ...)$averaged.net,
                         loci.list = loci.qtlnet(x, ...),
                         pheno.color="green", qtl.color="red",
                         vertex.color = node.color,
                         include.qtl=TRUE,
                         ...)
{
  node.names <- levels(edges[[1]])
  if(is.null(node.names))
    node.names <- unique(c(as.character(edges[[1]]), as.character(edges[[2]])))

  if(is.null(loci.list) | !include.qtl)
    node.color <- pheno.color
  else {
    loci.data.frame <- data.frame(qtl = unlist(loci.list))
    loci.data.frame$pheno <- rep(names(loci.list), sapply(loci.list, length))

    pheno.names <- node.names
    node.names <- c(pheno.names, levels(loci.data.frame[[1]]))

    edges <- cbind.data.frame(cause = c(as.character(edges[[1]]),
                                as.character(loci.data.frame[[1]])),
                              effect = c(as.character(edges[[2]]),
                                as.character(loci.data.frame[[2]])),
                              width = c(edges[[3]],
                                rep(1, nrow(loci.data.frame))))

    node.color <- rep(qtl.color, length(node.names))
    node.color[node.names %in% pheno.names] <- pheno.color
  }

  ## Set up vertices
  vertex.color <- array(vertex.color, length(node.names))
  vertices <- data.frame(name = node.names, label = node.names, color = vertex.color)

  ## Great graph object (library igraph).
  graph.data.frame(edges, TRUE, vertices = vertices)
}
