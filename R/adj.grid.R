#' @title Find adjacency matrix
#' @description Find adjacency matrix
#' @param coords data frame of x, y coordinates
#' @param order order of neighborhood
#' @param spam save as sparse matrix?
#' @return Adjacency matrix
#' @import spam
#' @export
#'
#### function to compute adjacency matrix
adj.grid <- function(coords, order = 1, spam = F){
  adj <- apply(coords, 1, function(pt)
    (abs(as.numeric(pt["y"]) - as.numeric(coords[,'y'])) <= order) & (abs(as.numeric(pt["x"]) - as.numeric(coords[,"x"])) <=  order))
  diag(adj) <- F

if(spam){
  return(as.spam(adj))
}else{
	return(adj)
}
}
