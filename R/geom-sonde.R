#' @export
GeomSonde <- ggplot2::ggproto("GeomSonde", ggplot2::GeomPath,
  default_aes = ggplot2::aes(colour = "blue", size = 1,
                             linetype = 1, alpha = NA)
)

#' Title
#'
#' @param mapping
#' @param data
#' @param stat
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
geom_sonde <- function(mapping = NULL, data = NULL, stat = StatSonde,
                       position = "identity", na.rm = TRUE,
                       show.legend = NA, inherit.aes = TRUE,
                       ...) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSonde,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}
