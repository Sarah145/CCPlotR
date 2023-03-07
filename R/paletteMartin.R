#' Discrete palette generator
#'
#' Generate a palette of up to 15 colours. The colours are from the paletteMartin palette in the colorBlindess R package.
#'
#' @param n Number of colours to return. Max = 15.
#' @export
#' @examples
#' scales::show_col(paletteMartin(n=9))

paletteMartin <- function(n = 15){
  cols <- c("#009292FF", "#490092FF", "#6DB6FFFF", "#920000FF", "#FFB6DBFF", "#FFFF6DFF", "#B66DFFFF", "#006DDBFF", "#924900FF", "#FF6DB6FF", "#004949FF", "#DB6D00FF", "#B6DBFFFF", "#24FF24FF", "#000000FF")
  return(cols[1:n])
}