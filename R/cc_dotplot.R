library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)
library(forcats)

#' Dotplot Function
#'
#' This function plots a dotplot
#' @param cc_df A dataframe with columns 'source', 'target', 'ligand', 'receptor' and 'score'. See `toy_data` for example.
#' @param option Either 'A' or 'B'. Option A will plot the number of interactions between pairs of cell types, option B will plot the top `n_top_ints` interactions and their scores. 
#' @param n_top_ints The number of top interactions to plot. Only required for option B. 
#' @export
#' @import dplyr tidyr ggplot2 ggtext forcats
#' @importFrom plyr round_any
#' @examples
#' cc_dotplot(toy_data)
#' cc_dotplot(toy_data, option = 'B', n_top_ints = 10)

cc_dotplot <- function(cc_df, option = 'A', n_top_ints = 30){
  if(option == 'A'){
    input_df <- cc_df %>% mutate(source = factor(source), target = factor(target)) %>% 
      group_by(source, target, .drop = F) %>% tally()
    brks <- scales::pretty_breaks(n = 5)(c(round_any(min(input_df$n), 10, f = floor), round_any(max(input_df$n), 10, f = ceiling)))
    ggplot(input_df, 
           aes(x = target, y = source, fill = n, size = n)) +
      geom_point(pch = 21) +
      scale_x_discrete(name = '\nReceiver cell type') +
      scale_y_discrete(name = 'Sender cell type\n', limits = rev(levels(input_df$source))) +
      scale_size(range = c(2,7), limits = c(round_any(min(input_df$n), 10, f = floor), round_any(max(input_df$n), 10, f = ceiling)), breaks = rev(brks), name = 'No. of\ninteractions')  +
      scale_fill_viridis_c(option = 'B', name = 'No. of\ninteractions', limits = c(round_any(min(input_df$n), 10, f = floor), round_any(max(input_df$n), 10, f = ceiling)), breaks = brks) +
      guides(fill = 'none', size = guide_legend(override.aes = list(fill = viridis::inferno(n = length(brks), direction = -1)))) +
      theme_minimal(base_size = 14) +
      theme(axis.line = element_blank(),
            axis.text = element_text(colour = 'black'),
            panel.grid = element_line(colour = 'grey85'))
  }
  else if(option == 'B'){
    input_df <- cc_df %>% slice_max(order_by = score, n = n_top_ints) %>%
      mutate(lr_pair = factor(paste0(ligand, '|', receptor)), 
             cell_pair = factor(paste0(source, '&rarr;', target))) %>%
      arrange(score) %>%
      mutate(lr_pair = fct_inorder(lr_pair))
    brks <- scales::pretty_breaks(n = 5)(c(round_any(min(input_df$score), 1, f = floor), round_any(max(input_df$score), 1, f = ceiling)))
    ggplot(input_df, aes(x = cell_pair, y = lr_pair, size = score, fill = score)) +
      geom_point(pch = 21) +
      scale_size(range = c(2,6), limits = c(min(brks), max(brks)), breaks = rev(brks), name = 'Score')  +
      scale_fill_viridis_c(option = 'B', limits = c(min(brks), max(brks)), name = 'Score',  breaks = brks, direction = -1) +
      guides(fill = 'none', size = guide_legend(override.aes = list(fill = viridis::inferno(n = length(brks), direction = 1)))) +
      labs(x = 'Sender cell type &rarr; Receiver cell type',
           y = 'Ligand|Receptor') +
      theme_minimal(base_size = 14) +
      theme(axis.line = element_blank(),
            axis.ticks = element_line(colour = 'grey20'),
            axis.text.y = element_text(colour = 'black'),
            axis.text.x = element_markdown(colour = 'black', angle = 90, hjust = 1, vjust = 0.5),
            axis.title.x = element_markdown())
  } else {print('option must be either A or B')}
}
