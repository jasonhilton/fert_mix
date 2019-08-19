#' Plot each row or column of a matrix
#' 
#' @param sample_mat The matrix to plot
#' @param plot_dim Which dimension to plot - must be one of "columns" or "rows"
#' @param add Should the plot be added to an existing. False by default.
#' 
#' 
plot_matrix <- function(target_matrix, plot_dim =c("columns","rows"), add=F, ...){
  plot_dim <- match.arg(plot_dim)
  if (plot_dim == "rows"){
    target_matrix <- t(target_matrix)
  }
  
  iters <- 1:ncol(target_matrix)
  if (!add){
    plot(target_matrix[,1], ylim=range(target_matrix), type="l", ...)  
  }
  trash <- sapply(iters, function(j) points(target_matrix[,j], type="l", ...))
}

