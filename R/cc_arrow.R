library(dplyr)
library(ggplot2)

#' Paired Arrow Plot Function
#'
#' This function plots interactions between a pair of cell types
#' @param cc_df A dataframe with columns 'source', 'target', 'ligand', 'receptor' and 'score'. See `toy_data` for example.
#' @param cell_types A vector of which two cell types to plot.
#' @param option Either 'A' or 'B'. Option A will plot the top `n_top_ints` interactions between `cell_types` and their scores. Option B will plot the top `n_top_ints` interactions between `cell_types`, their scores and the expression of the ligand/receptor genes in the sender/receiver cell types. 
#' @param n_top_ints The number of top interactions to plot.
#' @param exp_df A dataframe containing the mean expression values for each ligand/receptor in each cell type. See `toy_exp` for an example. Only required for option B.
#' @param colours A named vector of colours for each cell type. Default is `paletteMartin()`, a colourblind-friendly palette. Only used for option A.
#' @param palette Which colour palette to use to show the mean expression. Should be one of the RColorBrewer sequential palettes. Only used for option B.
#' @export
#' @import dplyr ggplot2 
#' @examples
#' cc_arrow(toy_data, cell_types = c('B', 'CD8 T'), colours = c(`B` = 'hotpink', `CD8 T` = 'orange'))
#' cc_arrow(toy_data, cell_types = c('NK', 'CD8 T'), option = 'B', exp_df = toy_exp, n_top_ints = 10, palette = 'OrRd')
#'

cc_arrow <- function(cc_df, cell_types, option = 'A', n_top_ints = 15, exp_df = NULL, colours = setNames(paletteMartin(n = 2), cell_types), palette = 'BuPu'){
  if(option == 'A'){
    input_df <- cc_df %>% filter((source == cell_types[1] & target == cell_types[2]) | (source == cell_types[2] & target == cell_types[1])) %>%
      slice_max(order_by = score, n = n_top_ints)
    
    cell1_lr <- unique(c(input_df %>% filter(source == cell_types[1]) %>% pull(ligand),
                         input_df %>% filter(target == cell_types[1]) %>% pull(receptor)))
    
    cell2_lr <- unique(c(input_df %>% filter(target == cell_types[2]) %>% pull(receptor),
                         input_df %>% filter(source == cell_types[2]) %>% pull(ligand)))
    
    ypos_1 <- setNames(1:length(cell1_lr), cell1_lr)
    ypos_2 <- setNames(1:length(cell2_lr), cell2_lr)
    
    scaler <- max(c(ypos_1, ypos_2))
    scale_func <- function(x){
      scaled_x <- (((x - min(x))/(max(x)-min(x)))*(scaler-1))+1
      return(scaled_x)
    }
    
    ypos_1 <- scale_func(ypos_1)
    ypos_2 <- scale_func(ypos_2)
    
    xpos_1 <- -((ypos_1-scaler/2)-0.5)**2
    xpos_2 <- (((ypos_2-scaler/2)-0.5)**2)+max(abs(xpos_1))*3
    
    input_df <- input_df %>% mutate(ypos1 = case_when(
      source == cell_types[1] ~ ypos_1[ligand],
      target == cell_types[1] ~ ypos_1[receptor]
    ), ypos2 = case_when(
      source == cell_types[2] ~ ypos_2[ligand],
      target == cell_types[2] ~ ypos_2[receptor]
    ), xpos1 = case_when(
      source == cell_types[1] ~ xpos_1[ligand],
      target == cell_types[1] ~ xpos_1[receptor]
    ), xpos2 = case_when(
      source == cell_types[2] ~ xpos_2[ligand],
      target == cell_types[2] ~ xpos_2[receptor]
    ))
    
    xrange <- abs(min(xpos_1)) + max(xpos_2)
    
    ggplot(input_df) +
      geom_point(aes(x = xpos1 - xrange/40, y = ypos1), col = colours[cell_types[1]], pch = 15, size = 4) +
      geom_point(aes(x = xpos2 + xrange/40, y = ypos2), col = colours[cell_types[2]], pch = 15, size = 4) +
      geom_segment(data = input_df %>% filter(source == cell_types[1]), 
                   aes(x = xpos1, xend = xpos2, y = ypos1, yend = ypos2, linewidth = score),
                   arrow = arrow(length = unit(3.5, 'mm'), type = 'closed'), show.legend = F) +
      geom_segment(data = input_df %>% filter(source == cell_types[2]), 
                   aes(x = xpos2, xend = xpos1, y = ypos2, yend = ypos1, linewidth = score),
                   arrow = arrow(length = unit(3.5, 'mm'), type = 'closed'), show.legend = F) +
      annotate('text', x = c(min(xpos_1) - xrange/8,  max(xpos_2) + xrange/8), y = scaler/2, 
               label = cell_types, size = 6, angle = c(90,270), col = colours[cell_types]) +
      annotate('text', x = xpos_1 - xrange/22, y = ypos_1, label = cell1_lr, hjust = 1) +
      annotate('text', x = xpos_2 + xrange/22, y = ypos_2, label = cell2_lr, hjust = 0) +
      scale_linewidth(range = c(0.1,1.2)) +
      scale_y_reverse() +
      scale_x_continuous(limits = c(min(xpos_1) - xrange/5, max(xpos_2) + xrange/5)) +
      theme_void(base_size = 14) +
      theme(plot.margin = margin(10,20,10,20)) +
      coord_cartesian(clip = 'off')
    } else if(option == 'B'){
      if(is.null(exp_df)){print('exp_df is required for option B')}
      
      input_df <- cc_df %>% filter((source == cell_types[1] & target == cell_types[2]) | (source == cell_types[2] & target == cell_types[1])) %>%
        slice_max(order_by = score, n = n_top_ints)
      
      cell1_lr <- unique(c(input_df %>% filter(source == cell_types[1]) %>% pull(ligand),
                           input_df %>% filter(target == cell_types[1]) %>% pull(receptor)))
      
      cell2_lr <- unique(c(input_df %>% filter(target == cell_types[2]) %>% pull(receptor),
                           input_df %>% filter(source == cell_types[2]) %>% pull(ligand)))
      
      ypos_1 <- setNames(1:length(cell1_lr), cell1_lr)
      ypos_2 <- setNames(1:length(cell2_lr), cell2_lr)
      
      scaler <- max(c(ypos_1, ypos_2))
      scale_func <- function(x){
        scaled_x <- (((x - min(x))/(max(x)-min(x)))*(scaler-1))+1
        return(scaled_x)
      }
      
      ypos_1 <- scale_func(ypos_1)
      ypos_2 <- scale_func(ypos_2)
      
      xpos_1 <- -((ypos_1-scaler/2)-0.5)**2
      xpos_2 <- (((ypos_2-scaler/2)-0.5)**2)+max(abs(xpos_1))*3
      
      gene_df <- as.data.frame(exp_df %>% mutate(cell_gene = paste0(cell_type, '|', gene))) %>%
        filter(cell_type %in% cell_types, gene %in% unique(c(cell1_lr, cell2_lr)))
      rownames(gene_df) <- gene_df$cell_gene
      
      input_df <- input_df %>% mutate(ypos1 = case_when(
        source == cell_types[1] ~ ypos_1[ligand],
        target == cell_types[1] ~ ypos_1[receptor]
      ), ypos2 = case_when(
        source == cell_types[2] ~ ypos_2[ligand],
        target == cell_types[2] ~ ypos_2[receptor]
      ), xpos1 = case_when(
        source == cell_types[1] ~ xpos_1[ligand],
        target == cell_types[1] ~ xpos_1[receptor]
      ), xpos2 = case_when(
        source == cell_types[2] ~ xpos_2[ligand],
        target == cell_types[2] ~ xpos_2[receptor]
      ), col1 = ifelse(source == cell_types[1], gene_df[paste0(cell_types[1], '|', ligand), 'mean_exp'], gene_df[paste0(cell_types[1], '|', receptor), 'mean_exp']),
      col2 = ifelse(source == cell_types[2], gene_df[paste0(cell_types[2], '|', ligand), 'mean_exp'], gene_df[paste0(cell_types[2], '|', receptor), 'mean_exp']))
      
      xrange <- abs(min(xpos_1)) + max(xpos_2)
      
      ggplot(input_df) +
        geom_point(aes(x = xpos1 - xrange/40, y = ypos1, fill = col1), col = 'black', pch = 22, size = 5) +
        geom_point(aes(x = xpos2 + xrange/40, y = ypos2, fill = col2), col = 'black', pch = 22, size = 5) +
        geom_segment(data = input_df %>% filter(source == cell_types[1]), 
                     aes(x = xpos1, xend = xpos2, y = ypos1, yend = ypos2, linewidth = score),
                     arrow = arrow(length = unit(3.5, 'mm'), type = 'closed'), show.legend = F) +
        geom_segment(data = input_df %>% filter(source == cell_types[2]), 
                     aes(x = xpos2, xend = xpos1, y = ypos2, yend = ypos1, linewidth = score),
                     arrow = arrow(length = unit(3.5, 'mm'), type = 'closed'), show.legend = F) +
        annotate('text', x = c(min(xpos_1) - xrange/8,  max(xpos_2) + xrange/8), y = scaler/2, label = cell_types, size = 6, angle = c(90,270)) +
        annotate('text', x = xpos_1 - xrange/22, y = ypos_1, label = cell1_lr, hjust = 1) +
        annotate('text', x = xpos_2 + xrange/22, y = ypos_2, label = cell2_lr, hjust = 0) +
        scale_fill_stepsn(colours = RColorBrewer::brewer.pal(8, palette), 
                          name = 'Mean expression', n.breaks = 8, show.limits = T) +
        scale_linewidth(range = c(0.1,1.2)) +
        scale_y_reverse() +
        scale_x_continuous(limits = c(min(xpos_1) - xrange/5, max(xpos_2) + xrange/5)) +
        guides(fill = guide_coloursteps(frame.colour = 'black', frame.linewidth = 0.4, ticks = F, show.limits = T)) +
        theme_void(base_size = 14) +
        theme(plot.margin = margin(10,20,10,20),
              legend.key.width = unit(dev.size()[1] / 9, "inches"),
              legend.position = 'bottom') +
        coord_cartesian(clip = 'off')
        
      } else(print('option must be either A or B'))
}
