library(dplyr)
library(ggbump)
library(ggplot2)
library(stringr)

#' Sigmoid Plot Function
#'
#' This function plots interactions using the `geom_sigmoid` function from the `ggbump` R package 
#' @param cc_df A dataframe with columns 'source', 'target', 'ligand', 'receptor' and 'score'. See `toy_data` for example.
#' @param n_top_ints The number of top interactions to plot.
#' @param colours A named vector of colours for each cell type. Default is `paletteMartin()`, a colourblind-friendly palette.
#' @export
#' @import dplyr ggplot2 ggbump stringr
#' @examples
#' cc_sigmoid(toy_data)
#' cc_sigmoid(toy_data, colours = c(`B` = 'hotpink', `CD8 T` = 'orange', `NK` = 'cornflowerblue'), n_top_ints = 25)
#'

cc_sigmoid <- function(cc_df, n_top_ints = 20, colours = paletteMartin()){
  input_df <- cc_df %>%
    slice_max(order_by = score, n = n_top_ints)
  
  source_ligands <- setNames(seq(1, length(unique(paste0(input_df$source, '|', input_df$ligand))), 1), unique((paste0(input_df$source, '|', input_df$ligand))))   
  target_receptors <- setNames(seq(1, length(unique(paste0(input_df$target, '|', input_df$receptor))), 1), unique(paste0(input_df$target, '|', input_df$receptor)))
  
  input_df <- input_df %>% mutate(y1 = source_ligands[paste0(source, '|', ligand)], y2 = target_receptors[paste0(target, '|', receptor)])
  
  ggplot(input_df, aes(y = as.factor(y1), yend = as.factor(y2), x = 0, xend = 1)) +
    geom_sigmoid() +
    annotate('text', x = c(0,1), y = as.factor(0), label = c('Sender\n', 'Reciever\n'), size = 6) +
    annotate('text', x = c(-0.5,1.5), y = as.factor(round(max(c(source_ligands, target_receptors))/2)), label = c('Ligand', 'Receptor'), angle = c(90, 270), size = 6) +
    annotate('text', x = -0.05, y = as.factor(source_ligands), label = str_extract(names(source_ligands), '[^|]+$'), hjust = 1) +
    annotate('text', x = 1.05, y = as.factor(target_receptors), label = str_extract(names(target_receptors), '[^|]+$'), hjust = 0) +
    geom_point(aes(col = source, x = 0, y = as.factor(y1)), show.legend = T, pch = 15, size = 7) +
    geom_point(aes(col = target, x = 1, y = as.factor(y2)), show.legend = F, pch = 15, size = 7) +
    scale_x_continuous(limits = c(-0.5, 1.5)) +
    scale_y_discrete(limits = as.factor(seq(max(source_ligands, target_receptors), 0, -1))) +
    scale_colour_manual(values = colours, breaks = unique(c(input_df$source, input_df$receptor)), name = 'Cell type') +
    theme_void(base_size = 14) +
    theme(legend.position = 'bottom',
          plot.margin = margin(t = 20)) +
    coord_cartesian(clip = 'off')
}


  
    