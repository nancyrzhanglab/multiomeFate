plot_trainTest <- function(cv_fit_list,
                           axis_size = 8,
                           bool_include_lambda_title = TRUE,
                           quantile_vec = c(0.1, 0.5, 0.9),
                           xlab = "Lambda+1 (Log10-scale tickmarks)",
                           ylab_test = "Negative loglikelihood (Testing)",
                           ylab_train = "Negative loglikelihood (Training)",
                           title_size = 10,
                           title_test = "",
                           title_train = ""){
  stopifnot(inherits(cv_fit_list, "lineage_cv"))
  
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
    lambda_value = res_test$lambda,
    title = title_test,
    title_size = title_size,
    xlab = xlab,
    ylab = ylab_test
  )
  
  plot1 <- cowplot::plot_grid(plotlist = plot_list, ncol = 2)
  plot1
}

.prepare_trainTest_data <- function(cv_fit_list,
                                    quantile_vec,
                                    what){
  stopifnot(what %in% c("train", "test"),
            length(quantile_vec) == 3,
            all(diff(quantile_vec) > 0))
  
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
                   quantile = rep(c("lower", "median", "upper"), 
                                  each = length(lambda_sequence)))
  
  list(df = df,
       lambda = lambda)
}

.plot_trainTest_helper <- function(axis_size,
                                   df,
                                   lambda_value,
                                   title,
                                   title_size,
                                   xlab,
                                   ylab){
  plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = lambda, y = value))
  
  # add polygon
  plot1 <- plot1 + ggplot2::geom_polygon(
    data = data.frame(x = c(df$lambda[df$quantile == 'lower'], rev(df$lambda[df$quantile == 'upper'])),
                      y = c(df$value[df$quantile == 'lower'], rev(df$value[df$quantile == 'upper'])),
                      value = "tmp"),
    ggplot2::aes(x = x, 
                 y = y,
                 fill = value)
  )
  plot1 <- plot1 + ggplot2::scale_fill_manual(values = c(tmp = "gray"))
  
  # add lines
  plot1 <- plot1 + ggplot2::geom_point(data = subset(df, quantile == 'median'), 
                                       ggplot2::aes(x = lambda, y = value), 
                                       shape = 16) 
  plot1 <- plot1 + ggplot2::geom_line(data = subset(df, quantile == 'median'), 
                                      ggplot2::aes(x = lambda, y = value), 
                                      size = 1) 
  
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