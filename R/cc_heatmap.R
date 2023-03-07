library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)
library(forcats)

#' Heatmap Function
#'
#' This plots a heatmap
#' @param cc_df A dataframe with columns 'source', 'target', 'ligand', 'receptor' and 'score'. See `toy_data` for example.
#' @param option Either 'A' or 'B'. Option A will plot the number of interactions between pairs of cell types, option B will plot the top `n_top_ints` interactions and their scores. 
#' @param n_top_ints The number of top interactions to plot. Only required for option B. 
#' @export
#' @import dplyr tidyr ggplot2 ggtext forcats
#' @examples
#' cc_heatmap(toy_data)
#' cc_heatmap(toy_data, option = 'B', n_top_ints = 10)

cc_heatmap <- function(cc_df, option = 'A', n_top_ints = 30){
  if(option == 'A'){
    input_df <- cc_df %>% mutate(source = factor(source), target = factor(target)) %>% 
      group_by(source, target, .drop = F) %>% tally()
    ggplot(input_df, aes(x = target, y = source, fill = n)) +
      geom_tile(col = 'black', linewidth = 0.8) +
      scale_x_discrete(expand = c(0,0), name = '\nReceiver cell type') +
      scale_y_discrete(expand = c(0,0), name = 'Sender cell type\n', limits = rev(levels(input_df$source))) +
      scale_fill_viridis_c(option = 'C', name = 'Number of interactions\n') +
      guides(fill = guide_colourbar(title.position = 'right', frame.colour = 'black', frame.linewidth = 0.4, ticks = T)) +
      theme_classic(base_size = 14) +
      theme(axis.line = element_blank(),
            axis.text = element_text(colour = 'black'),
            legend.title = element_text(angle = 270, hjust = 0.5),
            legend.key.height = unit(dev.size()[1] / 8, "inches"))}
  else if(option == 'B'){
    input_df <- cc_df %>% slice_max(order_by = score, n = n_top_ints) %>%
      mutate(lr_pair = factor(paste0(ligand, '|', receptor)), 
             cell_pair = factor(paste0(source, '&rarr;', target))) %>%
      complete(cell_pair, lr_pair) %>%
      arrange(score) %>%
      mutate(lr_pair = fct_inorder(lr_pair))
    ggplot(input_df, aes(x = cell_pair, y = lr_pair, fill = score)) +
      geom_tile(col = 'white', linewidth = 0.25) +
      scale_fill_viridis_c(option = 'C', na.value = 'black', direction = 1) +
      scale_x_discrete(expand = c(0,0)) +
      guides(fill = guide_colourbar(title.position = 'right', frame.colour = 'black', frame.linewidth = 0.4, ticks = T)) +
      labs(x = '\nSender cell type &rarr; Receiver cell type',
           y = 'Ligand|Receptor\n', fill = 'Score') +
      theme_classic(base_size = 14) +
      theme(axis.line = element_blank(),
            axis.text.y = element_text(colour = 'black'),
            axis.text.x = element_markdown(colour = 'black', angle = 90, hjust = 1, vjust = 0.5),
            axis.title.x = element_markdown(),
            legend.title = element_text(angle = 270, hjust = 0.5),
            legend.key.height = unit(dev.size()[1] / 8.8, "inches"))
  } else {print('option must be either A or B')}
}
