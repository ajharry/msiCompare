#' @title bfdr
#' @description bfdr
#' @param pi1 adsfa
#' @param pi0 adsfa
#' @param cutoff adf
#' @param cutoff2 adf
#' @param pi0k adfadf
#' @return qval
#' @export
#'

##### According to Newton 2004 summary (4), estimate BFDR by


BFDRhat <- function(pi1, cutoff){
  sum((1-pi1)*ifelse(pi1 > cutoff, 1, 0))/sum(pi1 > cutoff)
}


##### Defined conversely by Ventrucci (i.e. cutoff2 = 1-cutoff):

BFDRhat_V <- function(pi0, cutoff2){
  sum(pi0*ifelse(pi0 <= cutoff2, 1, 0))/sum(pi0 <= cutoff2)
}

# bfdrv <- sapply(pi0k, function(x) BFDRhat_V(pi0k, x))
# #### Ventrucci says cutoff2 should be smallest pi0 that gives BDFRhat bigger than alpha
#
# min(pi0k[which(sapply(pi0k, function(x) BFDRhat_V(pi0k, x))>.05)])
#
# min(pvs[which(sapply(pvs, function(x) BFDRhat_V(pvs, x))>.05)])
#


####### Storey 2011, estimating FDR

FDRhat <- function(pi0, cutoff2){
  m <- length(pi0)
  Rt <- sum(pi0 <= cutoff2)
  if(Rt >0){
    m*cutoff2/Rt
  }else{
    m*cutoff2
  }
}

##pi0k is a smallest-to-largest sorted vector of posterior probabilities
# fdrs <- sapply(pi0k, function(x) FDRhat(pi0k, x))
# plot(bfdrv, fdrs)
# abline(a=0,b=1)


qvalues <- function(pi0k){
    fdrs <- sapply(pi0k, function(x) FDRhat(pi0k, x))
    sort <- order(pi0k)
    qval <- numeric(length(pi0k))
    for(i in 1:length(fdrs)){
      qval[i] <- min(fdrs[sort[which(sort == i):length(sort)]])
    }
    return(qval)
}

# q <- qvalues(pi0k)
#
#
# qq <- qvalues(pvs)
# qq
# qvs$qvalues
# plot(qq, qvs)
# abline(a=0,b=.487)

#my function differs from the qvalue package by the estimate of m_hat. i chose m.
