#' T-test on sample averages. A separate average is produced for each combination of tissue and condition.
#' @title t-test on averages
#' @description T-test on sample averages.
#' @param msset An MSImageSet object.
#' @param setSamp sample labels for msset
#' @param setCond condition of interest labels for msset
#' @param setBioRep (optional) biological replicate labels for msset
#' @return The estimate(s) of variance, baseline, condition difference, and p-value for the test.
#' @import Cardinal
#' @import lmerTest
#' @export
#'
tissueWiseANOVA <- function(msset,
                        setSamp=pData(msset)$sample,
                        setCond = pData(msset)$diagnosis,
                        setBioRep = NULL, logbase2 = F,statusUpdates= T){

  if(logbase2){
    print("Log transforming spectra.")
    spec <- log2(spectra(msset))
  }else{
    spec <- spectra(msset)
  }

  setSamp <- factor(setSamp)
  sampNames <- levels(setSamp)
  setCond <- factor(setCond)
  condNames <- levels(setCond)

  means <- c()
  samples <- c()
  cond <- c()


  if(is.null(setBioRep)){
    if(statusUpdates) print("Getting sample means.")
    for(s in sampNames){
      pixels1 <- (setSamp == s & setCond == condNames[1])
      pixels2 <- (setSamp == s & setCond == condNames[2])

      if(any(pixels1)){
        means <- cbind(means,apply(spec[,pixels1], 1, function(sp) mean(sp[is.finite(sp)])))
        samples <- c(samples, s)
        cond <- c(cond, condNames[1])
      }

      if(any(pixels2)){
        means <- cbind(means,apply(spec[,pixels2], 1, function(sp) mean(sp[is.finite(sp)])))
        samples <- c(samples, s)
        cond <- c(cond, condNames[2])
      }

    }

    if(statusUpdates) print("Fitting ANOVA model.")
    lmfits <- apply(means,1, function(x) lm(x~cond))

    if(statusUpdates) print("Summarizing results.")
    sig2b <- unname(unlist(lapply(lmfits, function(x) summary(x)$sigma^2))) #sigma2b estimate
    pvalue <- unname(unlist(lapply(lmfits, function(x) coef(summary(x))[2,'Pr(>|t|)']))) #pvalue
    intercept <- unname(unlist(lapply(lmfits, function(x) coef(summary(x))[1,'Estimate']))) #int estimate
    condDiff <- unname(unlist(lapply(lmfits, function(x) coef(summary(x))[2,'Estimate']))) #cond estimate

    if(statusUpdates) print("Done.")
    return(list = list(sig2b= sig2b,
                       pvalue = pvalue,
                       intercept=intercept,
                       condDiff=condDiff))

  }else{ ### if there are bio replicates
    if(statusUpdates) print("Getting sample means.")
    setBioRep <- factor(setBioRep)
    bioRepNames <- levels(setBioRep)
    bRep <- c()

    for(bio in bioRepNames){
      for(s in sampNames){
        pixels1 <- (setSamp == s & setCond == condNames[1] & setBioRep == bio)
        pixels2 <- (setSamp == s & setCond == condNames[2] & setBioRep == bio)

        if(any(pixels1)){
          means <- cbind(means,apply(spec[,pixels1], 1, function(sp) mean(sp[is.finite(sp)])))
          samples <- c(samples, s)
          cond <- c(cond, condNames[1])
          bRep <- c(bRep, bio)
        }

        if(any(pixels2)){
          means <- cbind(means,apply(spec[,pixels2], 1, function(sp) mean(sp[is.finite(sp)])))
          samples <- c(samples, s)
          cond <- c(cond, condNames[2])
          bRep <- c(bRep, bio)
        }

      }
    }

    if(statusUpdates) print("Fitting ANOVA model.")
    lmmfits <- apply(means, 1, function(m) lmer(m~cond + (1|bRep)))

    if(statusUpdates) print("Summarizing results.")
    sig2tec <- unname(unlist(lapply(lmmfits, function(x) as.data.frame(VarCorr(x))[2, 'vcov']))) #sigma2btec estimate
    sig2bio <- unname(unlist(lapply(lmmfits, function(x) as.data.frame(VarCorr(x))[1, 'vcov'])) )#sigma2bio estimate
    pvalue <- unname(unlist(lapply(lmmfits, function(x) anova(x)['Pr(>F)']))) #pvalue
    intercept <- unname(unlist(lapply(lmmfits, function(x) summary(x)$coef[1,'Estimate']))) #int estimate
    condDiff <- unname(unlist(lapply(lmmfits, function(x) summary(x)$coef[2,'Estimate']))) #cond estimate

    if(statusUpdates) print("Done.")
    return(list = list(pvalue = pvalue,
                       intercept=intercept,
                       condDiff=condDiff,
                       sig2bio = sig2bio,sig2tec= sig2tec))
  }


}
