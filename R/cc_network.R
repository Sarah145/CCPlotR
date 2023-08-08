#' Network Plot Function
#'
#' This function plots a network of representing the number of interactions between cell types
#' @param cc_df A dataframe with columns 'source', 'target', 'ligand', 'receptor' and 'score'. See `toy_data` for example.
#' @param option Either 'A' or 'B'. Option A will plot the number of interactions between pairs of cell types, option B will plot the top `n_top_ints` interactions and their scores. 
#' @param n_top_ints The number of top interactions to plot. Only required for option B. 
#' @param colours A vector of colours for each cell type. Default is `paletteMartin()`, a colourblind-friendly palette.
#' @param node_size Point size for nodes in option B.
#' @param label_size Size for labels in option B.
#' @export
#' @import dplyr ggplot2 ggraph scatterpie
#' @importFrom igraph graph_from_data_frame layout_with_kk V
#' @return Returns a plot generated with the ggplot2 package
#' @examples
#' cc_network(toy_data)
#' cc_network(toy_data, colours = c('orange', 'cornflowerblue', 'hotpink'), option = 'B')

cc_network <- function(cc_df, colours = paletteMartin(), option = 'A', n_top_ints = 20, node_size = 2.75, label_size = 4){
  target <- from <- name <- score <- ligand <- receptor <- NULL
  if(option == 'A'){
    graph <- graph_from_data_frame(cc_df %>% group_by(source, target) %>% tally())
    ggraph(graph, layout = 'linear', circular = TRUE) + 
      geom_edge_loop(aes(edge_width = n, color = as.factor(from)), end_cap = circle(20, 'pt'),
                     arrow = arrow(length = unit(5, 'mm'), angle = 20, type = 'closed'), 
                     show.legend = FALSE) +
      geom_edge_fan(aes(edge_width = n, color = as.factor(from)), 
                    arrow = arrow(length = unit(5, 'mm'), angle = 20, type = 'closed'), 
                    end_cap = circle(20, 'pt'), 
                    show.legend = FALSE) +
      geom_node_label(aes(label = name, col = name), 
                      #fill = alpha('white', 0.3), 
                      size = 6, fontface = 'bold', show.legend = FALSE, label.size = 2) +
      geom_node_label(aes(label = name), 
                      #fill = alpha('white', 0.3), 
                      size = 6, fontface = 'bold', show.legend = FALSE, label.size = 0, col = 'black') +
      scale_edge_colour_manual(values = ifelse(rep(is.null(names(colours)),length(colours)), colours, unname(colours[sort(unique(cc_df$source))]))) +
      scale_colour_manual(values = colours) +
      scale_edge_width(range = c(0.5,4)) +
      theme_void(base_size = 14) +
      coord_cartesian(clip = 'off')}
  else if(option == 'B'){
    input_df <- cc_df %>% slice_max(order_by = score, n = n_top_ints)
    graph <- graph_from_data_frame(input_df %>% select(ligand, receptor, score))
    xy <- layout_with_kk(graph)
    igraph::V(graph)$x <- xy[, 1]
    igraph::V(graph)$y <- xy[, 2]
    graph_df <- igraph::as_data_frame(graph, 'vertices')
    genes <- unique(c(input_df$ligand, input_df$receptor))
    celltypes <- unique(c(input_df$source, input_df$target))
    for(i in genes){
      for(j in celltypes){
        graph_df[i,j] <- sum((input_df$source == j & input_df$ligand == i) | (input_df$target == j & input_df$receptor == i))
      }
    }
    ggraph(graph, "manual", x = V(graph)$x, y = V(graph)$y) +
      geom_edge_loop(arrow = arrow(length = unit(2.5, 'mm'), angle = 20, type = 'closed'), 
                     show.legend = FALSE, check_overlap = TRUE, end_cap = circle(20, 'pt')) +
      geom_edge_fan(arrow = arrow(length = unit(2.5, 'mm'), angle = 20, type = 'closed'), 
                    end_cap = circle(20, 'pt'), check_overlap = TRUE,
                    show.legend = FALSE) +
      geom_scatterpie(
        cols = celltypes,
        data = graph_df,
        colour = NA,
        pie_scale = node_size
      ) + 
      geom_node_label(aes(label = name), fill = alpha('white', 0.7), size = label_size, fontface = 'bold', repel = FALSE) +
      scale_fill_manual(values = colours, name = 'Cell type') +
      coord_fixed(clip='off') +
      theme_graph(base_size = 14, base_family="sans") +
      theme(legend.position = "bottom")
  } else {print('option must be either A or B')}
}
