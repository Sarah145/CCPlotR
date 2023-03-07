library(dplyr)
library(ggplot2)
library(ggraph)

#' Network Plot Function
#'
#' This function plots a network of representing the number of interactions between cell types
#' @param cc_df A dataframe with columns 'source', 'target', 'ligand', 'receptor' and 'score'. See `toy_data` for example.
#' @param colours A vector of colours for each cell type. Default is `paletteMartin()`, a colourblind-friendly palette.
#' @export
#' @import dplyr ggplot2 ggraph
#' @importFrom igraph graph_from_data_frame
#' @examples
#' cc_network(toy_data)
#' cc_network(toy_data, colours = c('hotpink', 'orange', 'cornflowerblue'))

cc_network <- function(cc_df, colours = paletteMartin()){
  graph <- graph_from_data_frame(cc_df %>% group_by(source, target) %>% tally())
  ggraph(graph, layout = 'linear', circular = T) + 
    geom_edge_loop(aes(edge_width = n, color = as.factor(from)), 
                   arrow = arrow(length = unit(5, 'mm'), angle = 20, type = 'closed'), 
                   show.legend = F, alpha = 0.8) +
    geom_edge_fan(aes(edge_width = n, color = as.factor(from)), 
                  arrow = arrow(length = unit(5, 'mm'), angle = 20, type = 'closed'), 
                  end_cap = circle(3, 'mm'), 
                  show.legend = F, alpha = 0.8) +
    geom_node_label(aes(label = name), fill = alpha('white', 0.3), size = 6, nudge_y = 0.1, fontface = 'bold') +
    scale_edge_color_manual(values = colours) +
    scale_edge_width(range = c(0.5,4)) +
    theme_void()
}
