#' Dotplot Function
#'
#' This function plots a dotplot
#' @param cc_df A dataframe with columns 'source', 'target', 'ligand', 'receptor' and 'score'. See `toy_data` for example.
#' @param option Either 'A', 'B', 'CellPhoneDB' or 'Liana'. Option A will plot the number of interactions between pairs of cell types, option B will plot the top `n_top_ints` interactions and their scores. The 'CellPhoneDB' and 'Liana' options will generate a dotplot in the style of these popular tools.
#' @param n_top_ints The number of top interactions to plot. Only required for option B.
#' @export
#' @import dplyr tidyr ggplot2 ggtext forcats grDevices viridis
#' @importFrom plyr round_any
#' @importFrom scales pretty_breaks
#' @importFrom methods is
#' @return Returns a plot generated with the ggplot2 package
#' @examples
#' data(toy_data, package = 'CCPlotR')
#' cc_dotplot(toy_data)
#' cc_dotplot(toy_data, option = "B", n_top_ints = 10)
#' cc_dotplot(toy_data, option = "Liana", n_top_ints = 15)
cc_dotplot <- function(cc_df, option = "A", n_top_ints = 30) {
    stopifnot("'cc_df' must be a dataframe" = is(cc_df, "data.frame"))
    stopifnot("cc_df should contain columns named source, target, ligand, receptor and score. See `toy_data` for an example." = all(c('source', 'target', 'ligand', 'receptor', 'score') %in% colnames(cc_df)))
    stopifnot("option must be either 'A', 'B', 'CellPhoneDB' or 'Liana'" = option %in% c('A', 'B', 'CellPhoneDB', 'Liana'))
    
    target <- score <- ligand <- receptor <- lr_pair <- cell_pair <- NULL
    if (option == "A") {
        input_df <- cc_df %>%
            mutate(source = factor(source), target = factor(target)) %>%
            group_by(source, target, .drop = FALSE) %>%
            tally()
        brks <- scales::pretty_breaks(n = 5)(c(round_any(min(input_df$n), 10, f = floor), round_any(max(input_df$n), 10, f = ceiling)))
        ggplot(
            input_df,
            aes(x = target, y = source, fill = n, size = n)
        ) +
            geom_point(pch = 21) +
            scale_x_discrete(name = "\nReceiver cell type") +
            scale_y_discrete(name = "Sender cell type\n", limits = rev(levels(input_df$source))) +
            scale_size(range = c(2, 7), limits = c(round_any(min(input_df$n), 10, f = floor), round_any(max(input_df$n), 10, f = ceiling)), breaks = rev(brks), name = "No. of\ninteractions") +
            scale_fill_viridis_c(option = "B", name = "No. of\ninteractions", limits = c(round_any(min(input_df$n), 10, f = floor), round_any(max(input_df$n), 10, f = ceiling)), breaks = brks) +
            guides(fill = "none", size = guide_legend(override.aes = list(fill = inferno(n = length(brks), direction = -1)))) +
            theme_minimal(base_size = 14) +
            theme(
                axis.line = element_blank(),
                axis.text = element_text(colour = "black"),
                panel.grid = element_line(colour = "grey85")
            )
    } else if (option == "B") {
        input_df <- cc_df %>%
            slice_max(order_by = score, n = n_top_ints) %>%
            mutate(
                lr_pair = factor(paste0(ligand, "|", receptor)),
                cell_pair = factor(paste0(source, "&rarr;", target))
            ) %>%
            arrange(score) %>%
            mutate(lr_pair = fct_inorder(lr_pair))
        brks <- scales::pretty_breaks(n = 5)(c(round_any(min(input_df$score), 1, f = floor), round_any(max(input_df$score), 1, f = ceiling)))
        ggplot(input_df, aes(x = cell_pair, y = lr_pair, size = score, fill = score)) +
            geom_point(pch = 21) +
            scale_size(range = c(2, 6), limits = c(min(brks), max(brks)), breaks = rev(brks), name = "Score") +
            scale_fill_viridis_c(option = "B", limits = c(min(brks), max(brks)), name = "Score", breaks = brks, direction = -1) +
            guides(fill = "none", size = guide_legend(override.aes = list(fill = inferno(n = length(brks), direction = 1)))) +
            labs(
                x = "Sender cell type &rarr; Receiver cell type",
                y = "Ligand|Receptor"
            ) +
            theme_minimal(base_size = 14) +
            theme(
                axis.line = element_blank(),
                axis.ticks = element_line(colour = "grey20"),
                axis.text.y = element_text(colour = "black"),
                axis.text.x = element_markdown(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                axis.title.x = element_markdown()
            )
    } else if (option == "CellPhoneDB") {
        input_df <- cc_df %>%
            slice_max(order_by = score, n = n_top_ints) %>%
            mutate(
                lr_pair = factor(paste0(ligand, "_", receptor)),
                cell_pair = factor(paste0(source, "|", target))
            ) %>%
            arrange(score) %>%
            mutate(lr_pair = fct_inorder(lr_pair))
        pal <- c("#000000", "#3f008d", "#4923a9", "#ddf500", "#fa7200", "#ef000f")
        brks <- scales::pretty_breaks(n = 4)(c(round_any(min(input_df$score), 1, f = floor), round_any(max(input_df$score), 1, f = ceiling)))
        ggplot(input_df, aes(x = cell_pair, y = lr_pair, col = score, size = score)) +
            geom_point() +
            scale_colour_gradientn(colours = pal) +
            scale_size(range = c(2, 6), limits = c(min(brks), max(brks)), breaks = brks, name = "Score") +
            guides(colour = "none", size = guide_legend(override.aes = list(colour = colorRampPalette(pal)(length(brks))))) +
            theme_linedraw(base_size = 14) +
            theme(
                axis.text = element_text(colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                axis.title = element_blank(),
                panel.grid = element_blank(),
                legend.position = "bottom"
            )
    } else if (option == "Liana") {
        input_df <- cc_df %>%
            slice_max(order_by = score, n = n_top_ints) %>%
            mutate(lr_pair = factor(paste(ligand, "->", receptor))) %>%
            arrange(score) %>%
            mutate(lr_pair = fct_inorder(lr_pair))
        brks <- scales::pretty_breaks(n = 4)(c(round_any(min(input_df$score), 1, f = floor), round_any(max(input_df$score), 1, f = ceiling)))
        ggplot(input_df, aes(x = target, y = lr_pair, col = score, size = score)) +
            geom_point() +
            scale_size(range = c(2, 6), limits = c(min(brks), max(brks)), breaks = rev(brks), name = "Score") +
            scale_colour_viridis_c(option = "D", limits = c(min(brks), max(brks)), name = "Score", breaks = brks, direction = -1) +
            guides(colour = "none", size = guide_legend(override.aes = list(colour = viridis(n = length(brks), direction = 1)))) +
            labs(title = "Source", x = "Target", y = "Interactions (Ligand -> Receptor)") +
            facet_wrap(~source) +
            theme_bw(base_size = 14) +
            theme(
                strip.background = element_rect(fill = "white"),
                axis.text.x = element_text(colour = "#E69F00", face = "bold"),
                plot.title = element_text(hjust = 0.5)
            )
    } 
}
