#' @title fit model
#' @description fit model
#' @param msset adfdf
#' @return res
#' @import mvtnorm
#' @import lme4
#' @import spam
#' @import coda
#' @export
#'

compareMSI2 <- function(msset,conditionOfInterest,
                       feature, nsim=5000, burnin = 2500, trace = T,
                       piPrior = .1, seed = 1, logbase2 = F, coord = NULL,
                       type.neighbor = "radius", radius.neighbor = 1, maxdist.neighbor = NULL,
                       spInit = NULL,
                       bioRep = NULL,
                       techRep,
                       beta0 = 0, # Prior Mean for beta, only allow intercept
                       prec0 = .01, # Prior Precision Matrix of beta (vague)  (only allow intercept)
                       precAlpha0 = .01, #Prior Precision of slab (value of condition effect if it is not zero)
                       d0=.001, g0=.001,			# Hyperpriors for tau, taubio, tautec
                       rd = .00001 # ratio of varSpike/varSlab
){
  compareMSI(msset,conditionOfInterest,
                          feature, nsim, burnin, trace ,
                          piPrior, seed, logbase2, coord ,
                          type.neighbor, radius.neighbor , maxdist.neighbor,
                          spInit,
                          bioRep ,
                          techRep,
                          beta0 ,
                          prec0,
                          precAlpha0,
                          d0, g0,
                          rd)

}#function



