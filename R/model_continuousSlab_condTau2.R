#' @title fit model
#' @description fit model
#' @param msset adfdf
#' @return res
#' @import mvtnorm
#' @import lme4
#' @import spam
#' @export
#'

spatialComparison_condT2 <- function(msset,sample,conditionOfInterest,
                                     feature, nsim=5000, burnin = 2500, trace = T,
                                     piPrior = .1, seed = 1, logbase2 = F, coord = NULL,
                                     type.neighbor = "radius", radius.neighbor = 1, maxdist.neighbor = NULL,
                                     spInit = NULL){


  fitct2 <- compareMSI(msset = msset,
             conditionOfInterest = conditionOfInterest,
             feature = feature,
             nsim = nsim,
             burnin = burnin,
             trace = trace,
             piPrior = piPrior,
             seed = seed,
             logbase2= logbase2,
             coord = coord,
             type.neighbor = type.neighbor,
             radius.neighbor = radius.neighbor,
             maxdist.neighbor = maxdist.neighbor,
             spInit = spInit,
             bioRep = NULL,
             techRep = sample)

  for(i in length(fitct2)){
    names(fitct2[[i]])[which(names(fitct2[[i]]) == "sig2tec")] <- "sig2b"
    names(fitct2[[i]])[which(names(fitct2[[i]]) == "sig2tec_trace")] <- "sig2b_trace"
  }

  return(fitct2)
}#function



