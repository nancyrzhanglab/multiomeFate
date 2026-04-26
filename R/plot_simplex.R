#' Plot a simplex (ternary) diagram
#'
#' Creates a ternary plot for three-way compositional data using ggplot2,
#' without requiring the ggtern package.
#'
#' @param df A data frame where each row is a point to plot.
#' @param x_col Character name of the column in \code{df} for the first axis (bottom-left corner).
#' @param y_col Character name of the column in \code{df} for the second axis (bottom-right corner).
#' @param z_col Character name of the column in \code{df} for the third axis (top corner).
#' @param color_col Character name of the column in \code{df} for point colors. Default \code{NULL}.
#' @param size_col Character name of the column in \code{df} for point sizes. Default \code{NULL}.
#' @param col_palette Named vector of colors passed to \code{scale_color_manual}. Default \code{NULL}.
#' @param xlab Label for the first axis (bottom-left corner). Default \code{"x"}.
#' @param ylab Label for the second axis (bottom-right corner). Default \code{"y"}.
#' @param zlab Label for the third axis (top corner). Default \code{"z"}.
#' @param title Plot title. Default \code{"title"}.
#'
#' @return A ggplot2 object.
#' @export
plot_simplex <- function(df,
                         x_col,
                         y_col,
                         z_col,
                         color_col = NULL,
                         size_col = NULL,
                         col_palette = NULL,
                         xlab = "x",
                         ylab = "y",
                         zlab = "z",
                         title = "title") {
  # Normalize to sum to 1 and convert ternary to Cartesian coordinates.
  # Convention: x -> bottom-left vertex (0,0), y -> bottom-right vertex (1,0),
  # z -> top vertex (0.5, sqrt(3)/2).
  x <- df[[x_col]]
  y <- df[[y_col]]
  z <- df[[z_col]]
  tot <- x + y + z
  x <- x / tot
  y <- y / tot
  z <- z / tot

  s3 <- sqrt(3) / 2
  plot_df <- df
  plot_df$.cx <- y + z / 2
  plot_df$.cy <- z * s3

  # Triangle boundary
  tri <- data.frame(
    .cx = c(0, 1, 0.5, 0),
    .cy = c(0, 0, s3, 0)
  )

  # Build base plot
  base_aes <- ggplot2::aes(x = .data$.cx, y = .data$.cy)
  if (!is.null(color_col)) base_aes <- utils::modifyList(base_aes, ggplot2::aes(color = .data[[color_col]]))
  if (!is.null(size_col))  base_aes <- utils::modifyList(base_aes, ggplot2::aes(size  = .data[[size_col]]))

  plot1 <- ggplot2::ggplot(plot_df, base_aes) +
    ggplot2::geom_path(data = tri,
                       mapping = ggplot2::aes(x = .data$.cx, y = .data$.cy),
                       color = "black", inherit.aes = FALSE) +
    ggplot2::geom_point() +
    ggplot2::coord_fixed(xlim = c(-0.12, 1.12), ylim = c(-0.12, s3 + 0.12)) +
    ggplot2::theme_void() +
    ggplot2::labs(title = title) +
    ggplot2::annotate("text", x = 0,   y = -0.08, label = xlab, hjust = 0.5) +
    ggplot2::annotate("text", x = 1,   y = -0.08, label = ylab, hjust = 0.5) +
    ggplot2::annotate("text", x = 0.5, y = s3 + 0.08, label = zlab, hjust = 0.5)

  if (!is.null(col_palette)) {
    plot1 <- plot1 + ggplot2::scale_color_manual(values = col_palette)
  }

  plot1
}
