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
adj.grid <- function(coords, max.dist = NULL, radius = 1, spam = F, type = "radius", sample){
  coords <- cbind(coords, adjSamp = sample)
  
  if(!(type %in% c("radius", "dist"))) stop('Invalid neighborhood type. Choose either "dist" or "radius".')
  

  if(type == "radius"){
    #### all points within radius below OR above OR left OR right
    adj <- apply(coords, 1, function(pt)
      (pt['adjSamp'] == coords[,'adjSamp']) & (abs(as.numeric(pt["y"]) - as.numeric(coords[,'y'])) <= radius) & (abs(as.numeric(pt["x"]) - as.numeric(coords[,"x"])) <= radius))
  }else if(type =="dist"){
    #### all points no more than max.dist away (euclidean)
    adj <- apply(coords, 1, function(pt)
      (pt['adjSamp'] == coords[,'adjSamp']) & (sqrt((as.numeric(pt["y"]) - as.numeric(coords[,'y']))^2  + (as.numeric(pt["x"]) - as.numeric(coords[,"x"]))^2) <= max.dist))
  }
  
#   else if(type == "compass"){
#     #### all points closest one point below, above, left, and right
#     adj <- apply(coords, 1, function(pt){
#       #point on line below that is closest
#       yloc <- as.numeric(pt["y"]) - as.numeric(coords[,'y'])[as.numeric(pt["y"]) - as.numeric(coords[,'y'])>0]
#       if(length(yloc)>0){
#         yloc <- min(yloc)
#       line <- (as.numeric(pt["y"]) - as.numeric(coords[,'y']) == yloc)
#       val <- min(abs(as.numeric(pt["x"]) - coords[line,'x']))
#         nearBelow <- line & abs(as.numeric(pt["x"]) - coords[,'x']) == val
#       }else{
#         nearBelow <- rep(F, nrow(coords))
#       }
#       
#       #point on line above that is closest
#       yloc <- as.numeric(pt["y"]) - as.numeric(coords[,'y'])[as.numeric(pt["y"]) - as.numeric(coords[,'y'])<0]
#       if(length(yloc)>0){
#         yloc <- max(yloc)
#         line <- (as.numeric(pt["y"]) - as.numeric(coords[,'y']) == yloc)
#         val <- min(abs(as.numeric(pt["x"]) - coords[line,'x']))
#         nearAbove <- line & abs(as.numeric(pt["x"]) - coords[,'x']) == val
#       }else{
#         nearAbove <- rep(F, nrow(coords))
#       }
#       
#       #point on line to left that is closest
#       line <- (as.numeric(pt["y"]) - as.numeric(coords[,'y']) == 0) & (as.numeric(pt["x"]) - coords[,'x'] > 0)
#       
#       if(any(line)){
#         val <- min(as.numeric(pt["x"]) - coords[line,'x'])
#         nearLeft <- line & abs(as.numeric(pt["x"]) - coords[,'x']) == val
#       }else{
#         nearLeft <- rep(F, nrow(coords))
#       }
#       
#       #point on line to right that is closest
#       line <- (as.numeric(pt["y"]) - as.numeric(coords[,'y']) == 0) & (as.numeric(pt["x"]) - coords[,'x'] < 0)
#     if(any(line)){
#       val <- max(as.numeric(pt["x"]) - coords[line,'x'])
#         nearRight <- line & (as.numeric(pt["x"]) - coords[,'x'] == val)
#       }else{
#         nearRight <- rep(F, nrow(coords))
#       }
#       
#       return(nearBelow | nearAbove | nearLeft | nearRight)
#       
#     })
#     
# adj <- t(adj)
#     
#   }#end compas
  
  diag(adj) <- F
  
  if(spam){
    return(as.spam(adj))
  }else{
    return(adj)
  }
}

