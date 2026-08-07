#' Violin plot of fate potential within lineages, with an ANOVA summary
#'
#' The figure that shows how much of the variation in fate potential sits
#' \emph{between} lineages versus \emph{within} them --- the visual form of the
#' selection-versus-adaptation question. One violin per selected lineage, plus a
#' final \code{"All"} violin pooling every scored cell as a reference.
#'
#' Lineages are chosen by \emph{observed future size}, not by predicted score:
#' the \code{num_lineages_top} largest and \code{num_lineages_bottom} smallest
#' after filtering to lineages with at least \code{min_lineage_size} cells at
#' the current time point. If a lineage falls in both sets --- which happens
#' whenever the two counts overlap the number of surviving lineages --- it is
#' drawn once.
#'
#' The title reports two numbers. The ANOVA p-value comes from
#' \code{stats::oneway.test()}, which does \emph{not} assume equal variances
#' across lineages. The "lineage effect" percentage is the ordinary between-group
#' share of the total sum of squares, computed by \code{.anova_percentage()};
#' both are computed on the plotted lineages only, not on all cells.
#'
#' @param seurat_object A Seurat object whose \code{meta.data} supplies the
#'   lineage and time-celltype annotations. Only the metadata is used.
#' @param cell_imputed_score A named numeric vector of per-cell fate potentials,
#'   names being cell IDs --- the \code{cell_imputed_score} element of
#'   \code{\link{cyfer_finalize}}, on the log10 scale. It need not cover every
#'   cell in the object; cells absent from it are ignored, and \code{NA} entries
#'   are dropped.
#' @param assigned_lineage_variable Name of the \code{meta.data} column holding
#'   the lineage assignment.
#' @param time_celltype_variable Name of the \code{meta.data} column whose
#'   values identify time point (possibly crossed with cell type). Used only to
#'   count each lineage's cells at \code{day_later}.
#' @param day_later The value of \code{time_celltype_variable} that defines the
#'   future time point. Must occur in that column; asserted.
#' @param bool_add_future_size Whether to append each lineage's observed future
#'   size to its axis label, as \code{"lineage (n)"}. Default \code{TRUE}.
#' @param bool_anova Whether to run the ANOVA and put its result in the title.
#'   Default \code{TRUE}.
#' @param bool_mark_mean Whether to draw a crossbar at each violin's
#'   \bold{median} --- the name says mean, the code uses
#'   \code{stats::median}. Default \code{TRUE}.
#' @param bool_mark_max Whether to draw a crossbar at each violin's maximum.
#'   Default \code{FALSE}.
#' @param col Colour of those crossbars. Default \code{"#E69F00"}.
#' @param min_lineage_size Minimum number of \emph{current} time point cells for
#'   a lineage to be eligible. Default \code{2}.
#' @param num_lineages_top Number of largest-future-size lineages to draw.
#'   Default \code{10}.
#' @param num_lineages_bottom Number of smallest-future-size lineages to draw.
#'   Default \code{10}.
#' @param ylab Y axis label. Default \code{""}.
#' @param ylim Length-2 numeric passed to \code{ggplot2::ylim()}. Default
#'   \code{NA}, meaning no limits. An \code{NA} in one position leaves that end
#'   free.
#'
#' @returns A \code{ggplot} object.
#'
#' @export
plot_anova <- function(seurat_object,
                       cell_imputed_score,
                       assigned_lineage_variable,
                       time_celltype_variable,
                       day_later,
                       bool_add_future_size = TRUE,
                       bool_anova = TRUE,
                       bool_mark_mean = TRUE,
                       bool_mark_max = FALSE,
                       col = "#E69F00",
                       min_lineage_size = 2,
                       num_lineages_top = 10,
                       num_lineages_bottom = 10,
                       ylab = "",
                       ylim = NA){
  
  # grab the vector of which celltype-time each cell is
  assigned_lineage <- seurat_object@meta.data[,assigned_lineage_variable]
  names(assigned_lineage) <- Seurat::Cells(seurat_object)
  
  time_celltype <- seurat_object@meta.data[,time_celltype_variable]
  names(time_celltype) <- Seurat::Cells(seurat_object)
  stopifnot(day_later %in% time_celltype)
  
  # determine which lineages qualify to be in the plot
  lineage_vec <- assigned_lineage[names(cell_imputed_score)]
  tab_mat <- table(assigned_lineage, time_celltype)
  lineage_future_size <- tab_mat[, day_later]
  names(lineage_future_size) <- rownames(tab_mat)
  
  .plot_anova_helper(seurat_object = seurat_object,
                     cell_imputed_score = cell_imputed_score,
                     assigned_lineage_variable = assigned_lineage_variable,
                     lineage_future_size = lineage_future_size,
                     bool_add_future_size = bool_add_future_size,
                     bool_anova = bool_anova,
                     bool_mark_mean = bool_mark_mean,
                     bool_mark_max = bool_mark_max,
                     col = col,
                     min_lineage_size = min_lineage_size,
                     num_lineages_top = num_lineages_top,
                     num_lineages_bottom = num_lineages_bottom,
                     ylab = ylab,
                     ylim = ylim)
}

#########################

#' Draw the fate-potential violin plot
#'
#' Everything in \code{plot_anova()} downstream of computing each lineage's
#' observed future size. Split out so a caller who already has those sizes ---
#' from a source other than a time-celltype metadata column --- can supply them
#' directly as \code{lineage_future_size} and skip the Seurat table.
#'
#' Carries the \code{aes()} calls, and therefore the \code{@importFrom rlang
#' .data} that keeps \code{R CMD check} from flagging the column names as
#' undefined globals.
#'
#' @param lineage_future_size Named numeric vector, names being lineage IDs and
#'   values their observed cell counts at the future time point. Drives both the
#'   top/bottom selection and the \code{"(n)"} axis labels.
#' @param bool_mark_max Whether to draw a crossbar at each violin's maximum.
#'   Default \code{TRUE} here --- note this differs from \code{plot_anova()}'s
#'   \code{FALSE}.
#' @inheritParams plot_anova
#'
#' @returns A \code{ggplot} object.
#'
#' @importFrom rlang .data
#' @noRd
.plot_anova_helper <- function(seurat_object,
                               cell_imputed_score,
                               assigned_lineage_variable,
                               lineage_future_size,
                               bool_add_future_size = TRUE,
                               bool_anova = TRUE,
                               bool_mark_mean = TRUE,
                               bool_mark_max = TRUE,
                               col = "#E69F00",
                               min_lineage_size = 2,
                               num_lineages_top = 10,
                               num_lineages_bottom = 10,
                               ylab = "",
                               ylim = NA){
  stopifnot(length(names(cell_imputed_score)) == length(cell_imputed_score))
  
  if(any(is.na(cell_imputed_score))){
    cell_imputed_score <- cell_imputed_score[!is.na(cell_imputed_score)]
  }
  
  # grab the vector of which celltype-time each cell is
  assigned_lineage <- seurat_object@meta.data[,assigned_lineage_variable]
  names(assigned_lineage) <- Seurat::Cells(seurat_object)
  assigned_lineage <- assigned_lineage[names(cell_imputed_score)]
  
  # filter out lineages that too small
  lineage_vec <- assigned_lineage[names(cell_imputed_score)]
  tab_vec <- table(assigned_lineage)
  tab_vec <- tab_vec[tab_vec >= min_lineage_size] # current size needs to be big enough
  passed_lineages <- names(tab_vec)
  passed_cells_names <- names(assigned_lineage)[which(assigned_lineage %in% passed_lineages)]
  cell_imputed_score <- cell_imputed_score[passed_cells_names]
  lineage_vec <- assigned_lineage[names(cell_imputed_score)]
  lineage_future_size <- lineage_future_size[which(names(lineage_future_size) %in% passed_lineages)]
  
  # determine which lineages qualify to be in the plot
  lineage_names_ordered <- names(lineage_future_size)[order(lineage_future_size, decreasing = TRUE)]
  # Clamp both ends: asking for more lineages than exist must fall back to all
  # of them, not run off either edge of the vector.
  num_qualifying <- length(lineage_names_ordered)
  lineage_names_top <- lineage_names_ordered[
    seq_len(min(num_lineages_top, num_qualifying))]
  lineage_names_bottom <- lineage_names_ordered[
    seq(max(1, num_qualifying - num_lineages_bottom + 1), length.out =
          min(num_lineages_bottom, num_qualifying))]
  lineage_names <- unique(c(lineage_names_top, lineage_names_bottom))
  idx <- which(lineage_vec %in% lineage_names)
  
  if(bool_add_future_size){
    lineage_names_new <- sapply(lineage_names, function(lineage_name){
      paste0(lineage_name, " (", lineage_future_size[lineage_name], ")")
    })
    lineage_vec <- plyr::mapvalues(lineage_vec, from = lineage_names, to = lineage_names_new)
    lineage_names <- lineage_names_new
  }
  
  # form data frame
  df <- data.frame(lineage = lineage_vec[idx],
                   imputed_count = cell_imputed_score[idx])
  df_tmp <- df; df_tmp$lineage <- droplevels(as.factor(df_tmp$lineage))
  if(bool_anova) anova_res <- stats::oneway.test(imputed_count ~ lineage, data = df_tmp)
  df2 <- data.frame(lineage = "All",
                    imputed_count = cell_imputed_score)
  df <- rbind(df, df2)
  
  # compute percentage
  if(bool_anova){
    lineage_effect <- .anova_percentage(
      df = df_tmp,
      lineage_variable = "lineage",
      value_variable = "imputed_count"
    )
  }
  
  gray_vec <- rep("lightgray", length(unique(df$lineage)))
  names(gray_vec) <- unique(df$lineage)
  
  plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$lineage, y = .data$imputed_count))
  plot1 <- plot1 + ggplot2::geom_violin(trim = TRUE, 
                                        scale = "width", 
                                        ggplot2::aes(fill = .data$lineage))
  plot1 <- plot1 + ggplot2::scale_fill_manual(values = gray_vec) 
  plot1 <- plot1 + ggplot2::geom_jitter(shape=16, 
                                        position=ggplot2::position_jitter(0.2), alpha = 0.3, size = 0.5)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::scale_x_discrete(limits = c(lineage_names, "All"),
                                             guide = ggplot2::guide_axis(angle = 45))
  plot1 <- plot1 + ggplot2::ylab(ylab)
  
  if(!all(is.na(ylim))){
    plot1 <- plot1 + ggplot2::ylim(ylim[1], ylim[2])
  }
  
  if(bool_mark_mean) 
    plot1 <- plot1 + ggplot2::stat_summary(fun = stats::median, geom = "crossbar", 
                                           width = 0.75, color = col)
  
  if(bool_mark_max) 
    plot1 <- plot1 + ggplot2::stat_summary(fun = max, geom = "crossbar", 
                                           width = 0.75, color = col)
 
  if(bool_anova){
    plot1 <- plot1 + ggplot2::ggtitle(paste0("ANOVA -Log10(pvalue)=", round(-log10(anova_res$p.value), 2), 
                                             ", Lineage effect = ", lineage_effect, "%"))
  }
 
  plot1
}

##################

#' Between-lineage share of the total sum of squares
#'
#' The classical one-way ANOVA \eqn{R^2}: between-group sum of squares over
#' total sum of squares, as a percentage. Reported in the plot title as the
#' "lineage effect" --- the fraction of variation in fate potential explained by
#' knowing which lineage a cell belongs to, which is the quantity the
#' selection-versus-adaptation argument turns on.
#'
#' Computed directly rather than read off \code{stats::oneway.test()}, because
#' that function returns a Welch statistic and no variance decomposition.
#' Unadjusted for the number of lineages, so it rises mechanically as lineages
#' get smaller; compare across panels only at matched lineage counts.
#'
#' @param df A data frame holding one row per cell.
#' @param lineage_variable Name of the grouping column. Must be a
#'   \bold{factor}, and its \code{levels()} --- not its observed values ---
#'   define the groups, so a factor carrying unused levels will produce
#'   \code{NaN} from the empty groups. \code{.plot_anova_helper()} calls
#'   \code{droplevels()} before passing it in.
#' @param value_variable Name of the numeric column being decomposed.
#'
#' @returns A single numeric between 0 and 100, rounded to one decimal.
#'
#' @noRd
.anova_percentage <- function(df,
                              lineage_variable,
                              value_variable){
  stopifnot(is.factor(df[,lineage_variable]))
  imputed_count <- df[,value_variable]
  lineage <- df[,lineage_variable]
  
  total_std <- sum((imputed_count - mean(imputed_count))^2)
  
  within_lineage_std <- sum(sapply(levels(lineage), function(lineage_name){
    idx <- which(lineage == lineage_name)
    sum((imputed_count[idx] - mean(imputed_count[idx]))^2)
  }))
  
  across_lineage_std <- sum(sapply(levels(lineage), function(lineage_name){
    idx <- which(lineage == lineage_name)
    mean_val <- mean(imputed_count[idx])
    length(idx) * (mean_val - mean(imputed_count))^2 
  }))
  
  lineage_effect <- round(across_lineage_std/total_std*100,1)
  lineage_effect
}