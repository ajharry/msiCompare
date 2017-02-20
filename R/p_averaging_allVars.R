#' T-test on sample averages. A separate average is produced for each combination of tissue and condition.
#' @title t-test on averages
#' @description T-test on sample averages.
#' @param simSet An MSImageSet object.
#' @return The estimate of tissue-to-tissue variance, overall mean, estimate of condition difference, and p-value for the test.
#' @import Cardinal
#' @export
#'
p_averaging <- function(simSet){

  means <- c()
  samples <- c()
  diag <- c()

  for(s in sampleNames(simSet)){
    pixels1 <- (pData(simSet)$sample == s & pData(simSet)$diagnosis == levels(factor(pData(simSet)$diagnosis))[1])
    pixels2 <- (pData(simSet)$sample == s & pData(simSet)$diagnosis == levels(factor(pData(simSet)$diagnosis))[2])

    if(any(pixels1)){
      means <- cbind(means,rowMeans(spectra(simSet[,pixels1])))
      samples <- c(samples, s)
      diag <- c(diag, levels(factor(pData(simSet)$diagnosis))[1])
    }

    if(any(pixels2)){
      means <- cbind(means,rowMeans(spectra(simSet[,pixels2])))
      samples <- c(samples, s)
      diag <- c(diag, levels(factor(pData(simSet)$diagnosis))[2])
    }

  }

  lmfits <- apply(means,1, function(x) lm(x~diag))

  sig2b <- unlist(lapply(lmfits, function(x) summary(x)$sigma^2)) #sigma2b estimate
  pvalue <- unlist(lapply(lmfits, function(x) coef(summary(x))[2,'Pr(>|t|)'])) #pvalue
  intercept <- unlist(lapply(lmfits, function(x) coef(summary(x))[1,'Estimate'])) #int estimate
  condDiff <- unlist(lapply(lmfits, function(x) coef(summary(x))[2,'Estimate'])) #cond estimate


  return(list = list(sig2b= sig2b, pvalue = pvalue, intercept=intercept,condDiff=condDiff))
}
