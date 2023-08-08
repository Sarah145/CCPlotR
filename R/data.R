#' Toy data for CCPlotR
#'
#' A toy dataset of ligand-receptor interactions to demonstrate cell-cell interaction plots.
#'
#' @return ## `toy_data`
#' A data frame with 735 rows and 5 columns:
#' \describe{
#'   \item{source}{Cell type expressing the ligand}
#'   \item{target}{Cell type expressing the receptor}
#'   \item{ligand}{Ligand}
#'   \item{receptor}{Receptor}
#'   \item{score}{A score for each interaction e.g. -log10(aggregate_rank) returned by Liana}
#' }
#' @source This is a modified version of the toy dataset that comes with the Liana R package.
"toy_data"

#' Toy expression data for CCPlotR
#'
#' A dataframe showing the mean expression values for each ligand and receptor in each cell type.
#'
#' @return ## `toy_exp`
#' A data frame with 477 rows and 3 columns:
#' \describe{
#'   \item{cell_type}{Cell type}
#'   \item{gene}{Ligand/receptor gene}
#'   \item{mean_exp}{Mean (normalised) expression of lignad/receptor gene in cell type}
#' }
"toy_exp"