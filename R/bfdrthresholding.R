#' @title BFDR control helper
#' @description m
#' @param pi0 adfaf
#' @param alpha adfadf
#' @return list
#' @export
#'
### to control fdr at .05, cutoff is the smallest possible cutoff for which the corresponding BFDRhat is > alpha

BFDRhat_V <- function(pi0, cutoff2){
  sum(pi0*ifelse(pi0 <= cutoff2, 1, 0))/sum(pi0 <= cutoff2)
}
