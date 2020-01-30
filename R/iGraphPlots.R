plotOntology <- function(mysem1R, ruleNumber = 1)
{
  if(class(mysem1R) != "Rcpp_sem1R")
  {
    write(paste0("mysem1R must be sem1R class!"), stderr())
  }
  else
  {
    if(ruleNumber > length(mysem1R$getIgraphEdges()))
    {
      write(paste0("Max value of ruleNumber is ",length(mysem1R$getIgraphEdges())), stderr())
    }
    else
    {
      require(igraph)
      ruleIND <- ruleNumber
      rule1Edges <- mysem1R$getIgraphEdges()[[ruleIND]]
      rule1Nodes <- mysem1R$getIgraphNodes()[[ruleIND]]
      g <- graph_from_data_frame(rule1Edges, directed=TRUE, vertices=rule1Nodes)
      l <- layout.sugiyama(g, layers = rule1Nodes$nodeLevel, hgap = 5, vgap = 0.1)
      V(g)$color <- ifelse(V(g)$isRuleTerm == 1, "#FF8673", "#C4D8E2")
      plot(g,
           layout = l$layout,                   # draw graph as tree
           vertex.size = 30,                  # node size
           #vertex.color = '#C4D8E2',          # node color
           vertex.shape = "circle",
           vertex.label = gsub("(.{16,}?)\\s", "\\1\n", paste0(rule1Nodes$nodeID,"\n",rule1Nodes$nodeName)),        # node labels
           vertex.label.cex = .5,             # node label size
           vertex.label.family = "Helvetica", # node label family
           vertex.label.font = 2,             # node label type (bold)
           vertex.label.color = '#000000',    # node label size
           #edge.label = edge_labels,          # edge labels
           edge.label.cex = .7,               # edge label size
           edge.label.family = "Helvetica",   # edge label family
           edge.label.font = 2,               # edge label font type (bold)
           edge.label.color = '#000000',      # edge label color
           edge.arrow.size = .6,              # arrow size
           edge.arrow.width = .6              # arrow width
      ) 
    }
  }
}

plotHeatmap <- function(mysem1R, hypothesis, matrix, binmatrix, ruleNumber = 1)
{
  if(ruleNumber > length(hypothesis))
  {
    write(paste0("Max value of ruleNumber is ",length(mysem1R$getIgraphEdges())), stderr())
  }
  else if(class(mysem1R) != "Rcpp_sem1R")
  {
    write(paste0("mysem1R must be sem1R class!"), stderr())
  }
  else
  {
    require(pheatmap)
    if(mysem1R$ruleFormat == "col")
    {
      require(RColorBrewer)
      darkcols <- brewer.pal(8, "Dark2")
      darkcols <- darkcols[1:2]
      names(darkcols) <- c("POSITIVE", "NEGATIVE")
      anno_colors <- list(category = darkcols)
      mynodes <- mysem1R$getIgraphNodes()[[ruleNumber]]
      mynodes <- mynodes[mynodes$isRuleTerm == 1, ]
      rulename <- paste0("'",as.character(mynodes$nodeName), collapse = "' AND ")
      data2plot <- matrix[,c(hypothesis[[ruleNumber]]$coveredPOS, hypothesis[[ruleNumber]]$coveredNEG)]
      mydata_new_category  <- rep("NEGATIVE", ncol(data2plot))
      mydata_new_category[colSums(binmatrix[,c(hypothesis[[ruleNumber]]$coveredPOS, hypothesis[[ruleNumber]]$coveredNEG)]) > 0] <- "POSITIVE"
      mydata_new_anot <- data.frame(colnames(data2plot), category=mydata_new_category, row.names = 1)
      pheatmap(data2plot, annotation = mydata_new_anot, main = rulename, annotation_colors = anno_colors)
    }
  }
}
