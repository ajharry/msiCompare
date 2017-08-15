#' @title fit model
#' @description fit model
#' @param msset adfdf
#' @return res
#' @import mvtnorm
#' @import lme4
#' @import spam
#' @export
#'

spatialComparison_condT2nest <- function(msset,conditionOfInterest,
                                         feature, nsim=5000, burnin = 2500, trace = T,
                                         piPrior = .1, seed = 1, logbase2 = F, coord = NULL,
                                         type.neighbor = "radius", radius.neighbor = 1, maxdist.neighbor = NULL,
                                         spInit = NULL,
                                         bioRep = NULL,
                                         techRep){
  compareMSI(msset = msset,
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
             bioRep = bioRep,
             techRep = techRep)

}#function



