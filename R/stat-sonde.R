#' @export
StatSonde <- ggplot2::ggproto("StatSonde", ggplot2::Stat,
  compute_group = function(data, scales) {
    data$y <- RadioSonde::skewty(data$y)
    data$x <- RadioSonde::skewtx(data$x, data$y)
    data
  },

  required_aes = c("x", "y")
)

#' Title
#'
#' @param mapping
#' @param data
#' @param geom
#' @param position
#' @param na.rm
#' @param show.legend
#' @param inherit.aes
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
stat_sonde <- function(mapping = NULL, data = NULL, geom = "sonde",
                       position = "identity", na.rm = FALSE, show.legend = NA,
                       inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatSonde, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}



