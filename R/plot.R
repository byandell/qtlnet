plot.qtlnet <- function(x, cross,
                        pheno.output = summary(x)$averaged.net,
                        marker.list = qtlnet.pheno(x, cross, ...),
                        pheno.color="green", qtl.color="red",
                        include.qtl=TRUE,
                        type = "igraph",
                        ...)
{
  if(type == "igraph") {
    igraph.qtlnet(pheno.output, make.qtl.dir(marker.list),
                  pheno.color = pheno.color, qtl.color = qtl.color, ...)
  }
  else { ## Rgraphviz
    graphviz.plot.qtlnet(x, cross, pheno.output, marker.list,
                         pheno.color, qtl.color, include.qtl, ...)
  }
}
################################################################
igraph.qtlnet <- function(pheno.output,
                          qtl.output=NULL,
                          pheno.color="green", qtl.color="red", 
                          vertex.color = node.color, ...)
{
  node.color <- rep(qtl.color, length(node.names))
  node.color[node.names %in% pheno.names] <- pheno.color

  if(is.null(qtl.output))
    output <- pheno.output
  else {
    pheno.names <- attr(pheno.output, "node.names")
    node.names <- c(pheno.names, levels(qtl.output[[1]]))

    output <- cbind.data.frame(cause = c(as.character(pheno.output[[1]]),
                                 as.character(qtl.output[[1]])),
                               react = c(as.character(pheno.output[[2]]),
                                 as.character(qtl.output[[2]])),
                               width = c(pheno.output[[3]], rep(1, nrow(qtl.output))))
    attr(output, "node.names") <- node.names
  }
  invisible(create.directed.graph.object.igraph(output, vertex.color = vertex.color, ...))
}
################################################################
make.qtl.dir <- function(x)
{
  out <- data.frame(qtl = unlist(x))
  out$pheno <- rep(names(x), sapply(x, length))
  out
}
################################################################
create.directed.graph.object.igraph <-
  function(output,
           layout=layout.lgl,
           vertex.color = rep("transparent", length(node.names)),
           ...)
{
  node.names <- attr(output, "node.names")
  vertices <- data.frame(name = node.names, label = node.names, color = vertex.color)

  ## I don't know how to force a root or to get the same layout twice!
  tkplot(graph.data.frame(output, vertices = vertices), layout = layout, ...)
  invisible(list(edges = output, vertices = vertices))
}
