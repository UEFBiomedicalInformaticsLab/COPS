#' Plot similarity matrix as a heatmap
#' 
#' Order using hierarchical clustering
#'
#' @param sim_mat similarity matrix
#' @param method hclust method
#' @param palette color distiller palette
#' @param limits bounds of similarity
#' @param palette_direction set color direction
#' @param title plot title
#'
#' @return
#' @export
plot_similarity_matrix <- function(sim_mat, 
                                   method = "average", 
                                   palette = "RdBu", 
                                   #pos_color = "#6E0700", 
                                   #neg_color = "#05006E", 
                                   #midpoint = 0, 
                                   limits = c(-1,1), 
                                   palette_direction = 1, 
                                   title = NULL) {
  # Remove NAs before hclust
  sim_mat[is.na(sim_mat)] <- 0
  hc <- hclust(as.dist(-limits[2] - sim_mat), method = method)
  dat <- reshape2::melt(sim_mat[hc$order, rev(hc$order)])
  ggplot(dat, aes(Var1, Var2, fill = value)) + geom_tile() + theme_bw() + coord_fixed() + ggtitle(title) + 
    scale_fill_distiller(palette = palette, limits = limits, direction = palette_direction) + 
    #scale_fill_gradientn(colors = c(pos_color, "#FFFFFF", neg_color), 
    #                     values = c(0, (limits[1] + midpoint) / limits[1] + limits[2], 1), 
    #                     limits = limits) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
}