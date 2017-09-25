
### to control fdr at .05, cutoff is the smallest possible cutoff for which the corresponding BFDRhat is > alpha

BFDRhat_V <- function(pi0, cutoff2){
  sum(pi0*ifelse(pi0 <= cutoff2, 1, 0))/sum(pi0 <= cutoff2)
}

BFDRdecision <- function(pi0, alpha = .05){
  bfdr2 <- sapply(pi0, function(x) BFDRhat_V(pi0 = pi0, cutoff2 = x)) #calculate bfdr for each potential cutoff
  cutoff <- min(pi0[which(bfdr2 > alpha)]) #best cutoff is smallest pi0 that gives bfdr > alpha
  pi0 < cutoff
}