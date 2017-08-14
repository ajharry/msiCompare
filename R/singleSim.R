#' @title sim single tissue
#' @description sim single tissue
#' @param tau2 asfasdf
#' @return sim
#' @import Cardinal
#' @import mvtnorm
#' @import spam
#' @export
#'
simSingle <- function(
  reps = 1,
  diff = log2(1.5),
  tau2 = 1,
  sig2 = .1,
  seed = 100,
  size1 = 15,
  center.pattern = F,
  pattern =  ifelse((expand.grid(x=1:size1, y=1:size1)$x %in% (1+floor(size1/5)):(size1-floor(size1/5)) & expand.grid(x=1:size1, y=1:size1)$y %in% (1+floor(size1/5)):(size1-floor(size1/5))), 2, 1)
  ){

  size2 <- size1^2

  diagnosis <-  ifelse(pattern == 2, "Healthy", "Disease")
  coord <- expand.grid(x=1:size1, y=1:size1)


  samples <-sampleCAR(condDiff = diff,
                      coord = expand.grid(x=1:size1, y=1:size1),
                      sig2 = sig2, tau2 = tau2,
                      rho = .9999, nrep=reps,
                      pattern = pattern,
                      save = F, randomSeed = seed, center.pattern = center.pattern)




  simSet <-  MSImageSet(spectra = matrix(samples,ncol = size2),
                        coord = expand.grid(x=1:size1, y=1:size1))

  pData(simSet)$diagnosis <- diagnosis
  pData(simSet)$sample <- factor(rep("1", size2))

  #image(simSet, feature=1)

  return(list(simSet = simSet,
              reps = reps,
              diff = diff,
              tau2 = tau2,
              sig2 = sig2,
              seed = seed,
              size1 = size1,
              size2 = size2
  ))
}
