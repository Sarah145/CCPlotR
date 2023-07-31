library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)
library(forcats)
library(tibble)
library(ggh4x)
library(patchwork)

#' Heatmap Function
#'
#' This plots a heatmap
#' @param cc_df A dataframe with columns 'source', 'target', 'ligand', 'receptor' and 'score'. See `toy_data` for example.
#' @param option Either 'A', 'B', 'CellPhoneDB' or 'Liana'. Option A will plot the number of interactions between pairs of cell types, option B will plot the top `n_top_ints` interactions and their scores. The 'CellPhoneDB' and 'Liana' options will generate a heatmap in the style of these popular tools.  
#' @param n_top_ints The number of top interactions to plot. Only required for option B. 
#' @export
#' @import dplyr tidyr ggplot2 ggtext forcats tibble ggh4x patchwork
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' cc_heatmap(toy_data)
#' cc_heatmap(toy_data, option = 'B', n_top_ints = 10)
#' cc_heatmap(toy_data, option = 'CellPhoneDB')

cc_heatmap <- function(cc_df, option = 'A', n_top_ints = 30){
  target <- score <- ligand <- receptor <- lr_pair <- cell_pair <- cc <- cell1 <- cell2 <- n_ints <- total <- NULL
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
  } else if(option == 'CellPhoneDB'){
    input_df <- cc_df %>% mutate(source = factor(source), target = factor(target)) %>% 
      group_by(source,target, .drop = F) %>% tally() %>% 
      mutate(cc = paste0(source, target)) %>% group_by(cc) %>% 
      mutate(cell_pair = str_c(sort(c(source, target)), collapse='|')) %>% 
      group_by(cell_pair) %>% summarise(n_ints = sum(n)) %>% ungroup() %>% 
      separate(cell_pair, into = c('cell1', 'cell2'), sep='\\|')    
    ggplot(rbind(input_df, data.frame(cell1 = input_df$cell2, cell2 = input_df$cell1, n_ints = input_df$n_ints)), 
           aes(x = cell1, y = cell2, fill = n_ints)) +
      geom_tile(col = 'white', linewidth = 1, alpha = 0.95) +
      scale_x_discrete(limits = unique(input_df$cell1), expand = c(0,0), name = NULL) +
      scale_y_discrete(limits = rev(unique(input_df$cell1)), position = 'right', expand = c(0,0), name = NULL) +
      scale_fill_gradientn(colours = c("dodgerblue4", "peachpuff", "deeppink4"), name = 'Number of\ninteractions') +
      guides(fill = guide_colourbar(label.position = 'top')) +
      theme_minimal(base_size = 14) +
      theme(axis.text = element_text(colour = 'black'),
            axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5),
            legend.key.width = unit(3, 'lines'),
            legend.position = 'bottom')
  } else if(option == 'Liana'){
    # input_mat <- as.matrix(cc_df %>% mutate(source = factor(source), target = factor(target)) %>% 
    #   group_by(source, target, .drop = F) %>% tally() %>% 
    #   pivot_wider(names_from = source, values_from = n) %>%
    #   column_to_rownames(var = 'target'))
    # liana::liana_heatmap(input_mat)
    col_func <- grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 8, name = 'Dark2')))
    
    input_df <- cc_df %>% mutate(source = factor(source), target = factor(target)) %>% 
      group_by(source, target, .drop = F) %>% tally()
    strip <- strip_themed(background_x = elem_list_rect(fill = col_func(length(unique(input_df$target)))),
                          background_y = elem_list_rect(fill = col_func(length(unique(input_df$source)))))
    p1 <- ggplot(input_df %>% group_by(source) %>% mutate(total = sum(n)), aes(x = source, y = total, fill = source)) +
      geom_col(show.legend = F) +
      scale_fill_manual(values = col_func(length(unique(input_df$source)))) +
      scale_y_continuous(expand = c(0,0)) +
      theme_classic(base_size = 12) +
      theme(axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(colour = 'black'),
            axis.title = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = margin(t=0,b=0))
    
    p2 <- ggplot(input_df %>% group_by(target) %>% mutate(total = sum(n)), aes(y = target, x = total, fill = target)) +
      geom_col(show.legend = F) +
      scale_fill_manual(values = col_func(length(unique(input_df$target)))) +
      scale_x_continuous(expand = c(0,0)) +
      theme_classic(base_size = 12) +
      theme(axis.line.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(colour = 'black', angle = 90, hjust = 1, vjust = 0.5),
            axis.title = element_blank(),
            axis.ticks.y = element_blank(),
            plot.margin = margin(l=0, r=0))
    p3 <- ggplot(input_df, aes(x = target, y = source, fill = n)) +
      geom_tile() +
      scale_x_discrete(expand = c(0,0), name = 'Receiver (Cell types)') +
      scale_y_discrete(expand = c(0,0), name = 'Sender (Cell types)') +
      scale_fill_distiller(palette = 'PuRd', direction = 1, name = 'Frequency') +
      facet_grid2(source~target, scales = 'free', switch = 'both', strip = strip) +
      theme_minimal(base_size = 14) +
      theme(panel.grid = element_blank(),
            axis.text = element_text(colour = 'black'),
            legend.title = element_text(face = 'bold'),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            strip.background = element_rect(colour = 'transparent'),
            strip.text = element_text(colour = 'transparent', margin = margin(0,0,0,0, "cm")),
            panel.spacing = unit(0, 'lines'),
            plot.margin = margin(t = 0, r = 0))
    (p1 + plot_spacer() + plot_layout(widths = c(1.7, 0.3)))/(p3 + p2 + plot_layout(widths = c(1.7, 0.3))) + plot_layout(heights = c(0.3, 1.7), guides = 'collect')
  } else {print("option must be either 'A', 'B', 'CellPhoneDB' or 'Liana'")}
}
