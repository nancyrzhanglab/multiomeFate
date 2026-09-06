#' Cross-validation diagnostic: training and held-out curves over the lambda path
#'
#' The standard check on a \code{cyfer()} run. Two side-by-side panels, training
#' on the left and held out on the right, each showing the median across folds
#' of the unpenalized objective at every lambda, with a shaded band between the
#' outer two quantiles. A dashed vertical line on the test panel marks the
#' selected lambda.
#'
#' What a healthy curve looks like: the training objective falls monotonically
#' as lambda shrinks, while the test curve turns up again at small lambda, and
#' the dashed line sits in that interior minimum. Two failure modes are visible
#' at a glance --- a flat test curve, and a dashed line pinned at the right-hand
#' end. Both mean the fit did not move. Ties in the test curve resolve to the
#' \emph{largest} lambda, because \code{which.min()} takes the first hit and
#' \code{lambda_sequence} is decreasing, so "selected lambda equals the ceiling
#' I set" is a symptom of a broken fit rather than of under-regularization.
#'
#' The x axis is \code{lambda + 1} on a log10 scale, since the path runs down to
#' exactly 0.
#'
#' @param cv_fit_list An object of class \code{"cyfer"}, i.e. the return value of
#'   \code{\link{cyfer}}. Asserted.
#' @param axis_size Point size of the axis titles. Default \code{8}.
#' @param bool_include_lambda_title Whether to append the selected lambda to the
#'   test panel's title. Default \code{TRUE}.
#' @param fill_col Fill colour of the inter-quantile band. Default
#'   \code{"gray"}.
#' @param quantile_vec Length-3 increasing vector of quantiles across folds:
#'   band lower edge, centre line, band upper edge. Default
#'   \code{c(0.1, 0.5, 0.9)}. The centre line is also what the displayed lambda
#'   is chosen by, so \code{quantile_vec[2]} \bold{must} be \code{0.5} --- that
#'   is what \code{\link{cyfer_finalize}} selects by, and a different centre
#'   would draw a dashed line at a lambda the refit never uses. Only the two
#'   outer entries are free.
#' @param xlab X axis label, shared by both panels. Default
#'   \code{"Lambda+1 (Log10-scale tickmarks)"}.
#' @param ylab_test,ylab_train Y axis labels. Defaults name these curves
#'   "negative loglikelihood"; they are the unpenalized objective, for which
#'   lower is better, and not a literal log-likelihood.
#' @param title_size Point size of the panel titles. Default \code{10}.
#' @param title_test,title_train Panel titles. Default \code{""}.
#'
#' @returns A \code{ggplot} object: the two panels combined by
#'   \code{cowplot::plot_grid()}.
#'
#' @examples
#' \donttest{
#' data(priming_simulation)
#' cv <- cyfer(cell_features = priming_simulation$cell_features,
#'             cell_lineage = priming_simulation$cell_lineage,
#'             lineage_future_count = priming_simulation$lineage_future_count,
#'             lambda_initial = 3,
#'             lambda_sequence_length = 5,
#'             num_folds = 5,
#'             verbose = 0)
#' plot_trainTest(cv)
#' }
#' @export
plot_trainTest <- function(cv_fit_list,
                           axis_size = 8,
                           bool_include_lambda_title = TRUE,
                           fill_col = "gray",
                           quantile_vec = c(0.1, 0.5, 0.9),
                           xlab = "Lambda+1 (Log10-scale tickmarks)",
                           ylab_test = "Negative loglikelihood (Testing)",
                           ylab_train = "Negative loglikelihood (Training)",
                           title_size = 10,
                           title_test = "",
                           title_train = ""){
  stopifnot(inherits(cv_fit_list, "cyfer"))
  
  res_train <- .prepare_trainTest_data(cv_fit_list,
                                       what = "train",
                                       quantile_vec = quantile_vec) 
  res_test <- .prepare_trainTest_data(cv_fit_list,
                                      what = "test",
                                      quantile_vec = quantile_vec) 
  
  plot_list <- vector("list", 2)
  names(plot_list) <- c("train", "test")
  plot_list[["train"]] <- .plot_trainTest_helper(
    res_train$df,
    axis_size = axis_size,
    fill_col = fill_col,
    lambda_value = NULL,
    title = title_train,
    title_size = title_size,
    xlab = xlab,
    ylab = ylab_train
  )
  
  if(bool_include_lambda_title){
    title_test <- paste0(title_test, ", Lambda = ", round(res_test$lambda,3))
  }
  
  plot_list[["test"]] <- .plot_trainTest_helper(
    res_test$df,
    axis_size = axis_size,
    fill_col = fill_col,
    lambda_value = res_test$lambda,
    title = title_test,
    title_size = title_size,
    xlab = xlab,
    ylab = ylab_test
  )
  
  plot1 <- cowplot::plot_grid(plotlist = plot_list, ncol = 2)
  plot1
}

#' Reduce the per-fold CV curves to quantiles and pick the displayed lambda
#'
#' Stacks one fold's objective curve per column, takes \code{quantile_vec}
#' across folds at each lambda, and returns the result in the long form
#' \code{ggplot2} wants.
#'
#' @param cv_fit_list An object of class \code{"cyfer"}.
#' @param quantile_vec Length-3 strictly increasing vector of quantiles, whose
#'   middle entry must be \code{0.5}; asserted.
#' @param what Either \code{"train"} or \code{"test"}, selecting which curve to
#'   summarize.
#'
#' @returns A list with:
#'   \describe{
#'     \item{\code{df}}{data frame with \code{3 * length(lambda_sequence)} rows
#'       and columns \code{lambda} (already shifted by \code{+1} for the log
#'       axis), \code{value}, and \code{quantile_str}, the last taking the
#'       literal values \code{"lower"}, \code{"median"}, \code{"upper"}
#'       regardless of which outer quantiles were requested.}
#'     \item{\code{lambda}}{the unshifted lambda minimizing the median curve. At
#'       \code{what == "test"} this is what \code{cyfer_finalize()} selects.}
#'   }
#'
#' @noRd
.prepare_trainTest_data <- function(cv_fit_list,
                                    quantile_vec,
                                    what){
  stopifnot(what %in% c("train", "test"),
            length(quantile_vec) == 3,
            all(diff(quantile_vec) > 0))
  # The displayed lambda is chosen by the middle curve, and cyfer_finalize()
  # always refits at the median, so anything else would mark a lambda the refit
  # never uses.
  if(quantile_vec[2] != 0.5){
    stop("`quantile_vec[2]` must be 0.5 (got ", quantile_vec[2], "): the ",
         "marked lambda is chosen by the middle curve, and `cyfer_finalize()` ",
         "selects by the median.")
  }
  
  if(what == "train"){
    loglik_obj <- "train_loglik"
  } else {
    loglik_obj <- "test_loglik"
  }
  
  loglik_mat <- sapply(cv_fit_list, function(x){
    x[[loglik_obj]]
  })
  
  
  loglik_quantile <- apply(loglik_mat, 1, function(vec){
    stats::quantile(vec, 
                    probs = quantile_vec)
  })
  
  lambda_sequence <- cv_fit_list[[1]]$train_fit$lambda_sequence
  lambda <- lambda_sequence[which.min(loglik_quantile[2,])]
  
  df <- data.frame(lambda = rep(lambda_sequence+1, times = 3),
                   value = as.numeric(t(loglik_quantile)),
                   quantile_str = rep(c("lower", "median", "upper"), 
                                  each = length(lambda_sequence)))
  
  list(df = df,
       lambda = lambda)
}

#' Draw one panel of the cross-validation diagnostic
#'
#' The band is drawn with \code{geom_polygon()} on the lower quantile followed
#' by the reversed upper quantile, rather than \code{geom_ribbon()}, which is
#' why the data frame is reshaped rather than kept wide.
#'
#' Carries the \code{aes()} calls, and therefore the \code{@importFrom rlang
#' .data}.
#'
#' @param axis_size Point size of the axis titles.
#' @param df The \code{df} element of \code{.prepare_trainTest_data()}. Its
#'   \code{lambda} column must already carry the \code{+1} shift, since the
#'   panel uses a log10 x axis.
#' @param fill_col Fill colour of the band.
#' @param lambda_value Unshifted lambda at which to draw a dashed vertical
#'   marker; \code{NULL} draws none. Shifted by \code{+1} internally to match
#'   the axis.
#' @param title,title_size,xlab,ylab Panel title, its point size, and the axis
#'   labels.
#'
#' @returns A \code{ggplot} object.
#'
#' @importFrom rlang .data
#' @noRd
.plot_trainTest_helper <- function(axis_size,
                                   df,
                                   fill_col,
                                   lambda_value,
                                   title,
                                   title_size,
                                   xlab,
                                   ylab){
  plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$lambda, y = .data$value))
  
  # add polygon
  plot1 <- plot1 + ggplot2::geom_polygon(
    data = data.frame(x = c(df$lambda[df$quantile_str == 'lower'], rev(df$lambda[df$quantile_str == 'upper'])),
                      y = c(df$value[df$quantile_str == 'lower'], rev(df$value[df$quantile_str == 'upper'])),
                      value = "tmp"),
    ggplot2::aes(x = .data$x, 
                 y = .data$y,
                 fill = .data$value)
  )
  plot1 <- plot1 + ggplot2::scale_fill_manual(values = c(tmp = fill_col))
  
  # add lines
  plot1 <- plot1 + ggplot2::geom_point(data = subset(df, df$quantile_str == 'median'), 
                                       ggplot2::aes(x = .data$lambda, y = .data$value), 
                                       shape = 16) 
  plot1 <- plot1 + ggplot2::geom_line(data = subset(df, df$quantile_str == 'median'),
                                      ggplot2::aes(x = .data$lambda, y = .data$value),
                                      linewidth = 1)
  
  plot1 <- plot1 + ggplot2::scale_x_log10() 
  plot1 <- plot1 + ggplot2::labs(
    title = title,
    x = xlab,
    y = ylab)
  plot1 <- plot1 + Seurat::NoLegend()
  
  if(!is.null(lambda_value)){
    plot1 <- plot1 + ggplot2::geom_vline(
      xintercept = lambda_value+1,
      linetype = 2,
      color = "coral"
    )
  }
  
  plot1 <- plot1 + ggplot2::theme(
    plot.title = ggplot2::element_text(size = title_size),
    axis.title.x = ggplot2::element_text(size = axis_size),
    axis.title.y = ggplot2::element_text(size = axis_size)
  )
  
  
  plot1
}