#' Plot UMAP embedding (generic)
#'
#' @param obj object
#' @param ... additional parameters
#'
#' @return nothing. A plot is made
#' @export
plot_umap <- function(obj, ...) {
  UseMethod(generic = "plot_umap", object = obj)
}