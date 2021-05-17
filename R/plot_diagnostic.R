plot_diagnos_neigh <- function(obj, individual_plot = F, only_selected = F,
                               bool_rev = T,
                               forward_col = "darkblue", current_col = "red",
                               backward_col = "green", 
                               selected_col = "purple",
                               main_text = "", verbose = T, ...){
  stopifnot(class(obj) %in% c("chromatin_potential", "diagnos_neigh"))
  bool_chrom <- class(obj) == "chromatin_potential"
  
  if(bool_chrom){
    len <- length(obj$list_diagnos)
  } else {
    len <- length(obj)
  }
  
  if(individual_plot){
    if(bool_rev) order_vec <- rev(1:len) else order_vec <- 1:len
    for(iter in order_vec){
      if(verbose) print(iter)
      
      if(bool_chrom){
        tmp <- obj$list_diagnos[[as.character(iter)]]$recruit$postprocess$df_diag
        tmp2 <- t(tmp[,c("forward_num", "current_num", "backward_num")])
      } else {
        tmp <- obj[[as.character(iter)]]
        tmp2 <- t(tmp[,c("forward_num", "backward_num")])
      }
     
      max_val <- max(colSums(tmp2))
      tmp2 <- as.table(tmp2)
      colnames(tmp2) <- tmp$idx
      
      if(bool_chrom){
        graphics::barplot(tmp2, col=c(forward_col, current_col, backward_col), horiz = F,
              main = paste0(main_text, "Iteration ", iter), ...)
      } else {
        graphics::barplot(tmp2, col=c(forward_col, backward_col), horiz = F,
                          main = paste0(main_text, "Iteration ", iter), ...)
      }
      idx <- which(tmp$selected)
      for(i in idx){
        start <- 0.2+1.2*(i-1)
        mid <- start+0.5
        graphics::lines(rep(mid,2), c(0, max_val), col = selected_col, lty = 2, lwd = 2)
      }
    }
  } else {
    dat <- matrix(NA, nrow = len, ncol = ifelse(bool_chrom, 3, 2))
    
    for(i in 1:len){
      
      if(bool_chrom){
        tmp <- obj$list_diagnos[[as.character(i)]]$recruit$postprocess$df_diag
        tmp2 <- as.matrix(tmp[,c("forward_num", "current_num", "backward_num")])
      } else {
        tmp <- obj[[as.character(i)]]
        tmp2 <- as.matrix(tmp[,c("forward_num", "backward_num")])
      }
      
      if(only_selected){
        tmp2 <- tmp2[which(tmp$selected),]
      }
      
      tmp2 <- t(apply(tmp2, 1, function(x){x/sum(x)}))
      dat[i,] <- colMeans(tmp2)
    }
    
    dat <- as.table(t(dat))
    colnames(dat) <- 1:len
    if(bool_chrom){
      graphics::barplot(dat, col=c(forward_col, current_col, backward_col), horiz = F,
              main = main_text, ...)
    } else {
      graphics::barplot(dat, col=c(forward_col, backward_col), horiz = F,
                        main = main_text, ...)
    }
  }
  
  invisible()
}