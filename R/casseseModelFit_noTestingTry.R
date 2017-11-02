#' @title model fitting cassese et al
#' @description model fitting cassese et al
#' @param msset adfaf
#' @param roiFactor adfadf
#' @param logscale asdfa
#' @param thresholds adfdf
#' @return list
#' @import spdep
#' @export
#'


cass <- function(msset, roiFactor, logscale = T, thresholds = 1:10,cutoff=.05 ){
  if(logscale){
msdata2 <- log(t(spectra(msset))) #original code
#msdata2 <- log2(t(spectra(msset))) #i typically work on base 2 scale
  }else{
    msdata2 <- t(spectra(msset))
}
xy <- as.matrix(coord(msset))
roi <- as.numeric(roiFactor)
# DATA INPUT:
# msdata: n x m matrix of MSI data; n=pixels and m=features
# xy: 2-column matrix of coordinates of pixels (1st column: x-coordinate; 2nd column: y- coordinate)
# roi: vector of integers coding the membership of each pixel to the regions of interest
# apply if necessary log-transformation

#msdata2 = msdata
#msdata2 = log(msdata)

# define significance cutoff
p_cutoff = 0.05
# calculate simple t-test p-values
ttestpvals = c()
for(i in 1:ncol(msdata2)) {
  ttestpvals = c(ttestpvals, t.test(msdata2[,i] ~ as.factor(roi))$p.value)
  }
padj_ttestpvals = p.adjust(ttestpvals, method="BH")
#table(padj_ttestpvals <= p_cutoff)

# perform CAR regression for accounting for spatial autocorrelation; load 'spdep' library


# define neighborhood thresholds (default: 1 to 10 pixel units)
#thresholds = 1
#thresholds = 1:10 #original code

# initialization of results vectors aics = c()
pvals = c()
morans_i = c()
morans_pvals = c()
aics = c()
sigma2 = c()
intercept = c()
condeffect = c()
# create distance matrix (caution: can become very big: O=f(n*n))
Euclid_distance=dist(xy)
Eucl_mat=as.matrix(Euclid_distance)

for(threshold in thresholds) {

  print(paste("threshold", threshold))
  # define neighbors
  neighbor=ifelse(Eucl_mat>threshold,0,1)
  diag(neighbor)=0
  # consider only neighbors from same ROI
  neighbor[roi==1,roi==2] = 0
  neighbor[roi==2,roi==1] = 0

  # create spdep neighbor list structure
  xy_weights = mat2listw(neighbor)

  # create for each feature a CAR model and a Moran's test
  for(i in 1:ncol(msdata2)) {
    print(paste(i, "/" , ncol(msdata2), names(msdata2)[i]))

    # CAR model
    #added try to catch singular fit
    z = try(spautolm(msdata2[,i] ~ as.numeric(roi), family = "CAR", listw = xy_weights))

    if (class(z) != "try-error"){
      z1 = summary(z)
      pvals = c(pvals, z1$Coef[8])
      aics = c(aics,AIC(z1))
      intercept = c(intercept, z1$Coef[1])
      condeffect = c(condeffect, z1$Coef[2])
      sigma2 = c(sigma2, z$fit$s2)
    }else{
      sigma2 = c(sigma2, NA)
      intercept = c(intercept, NA)
      condeffect = c(condeffect, NA)
      pvals = c(pvals, NA)
      aics = c(aics,NA)
    }

    # Moran's test #added try
    mt = try(moran.test(msdata2[,i], xy_weights))

    if (class(mt) != "try-error"){
      morans_i = c(morans_i, mt$estimate[[1]])
      morans_pvals = c(morans_pvals, mt$p.value)
    }else{
      morans_i = c(morans_i,NA)
      morans_pvals = c(morans_pvals, NA)

    }
  }
}

return(list = list(results = resCass(fit = list(aics = aics,
                morans_i = morans_i,
                morans_pvals = morans_pvals,
                pvals = pvals,
                padj_ttestpvals = padj_ttestpvals,
                sigma2 = sigma2,
                cond = condeffect,
                intercept = intercept), thresholds = thresholds, mz = mz(msset), p_cutoff = cutoff),
                rawResults = list(aics = aics,
                                              morans_i = morans_i,
                                              morans_pvals = morans_pvals,
                                              pvals = pvals,
                                              padj_ttestpvals = padj_ttestpvals,
                                              sigma2 = sigma2,
                                              cond = condeffect,
                                              intercept = intercept)))
}




resCass <- function(fit, thresholds = 1:10, p_cutoff = .05, mz){

  resultsList <- fit

  # reshape results in matrices: rows contain features and columns represent distance tresholds
  aic_matrix = t(matrix(resultsList$aics, nrow=length(thresholds), byrow=T))
  morans_matrix = t(matrix(resultsList$morans_i, nrow=length(thresholds), byrow=T))
  morans_pvals_matrix = t(matrix(resultsList$morans_pvals, nrow=length(thresholds), byrow=T))
  pvals_matrix = t(matrix(resultsList$pvals, nrow=length(thresholds), byrow=T))

  #### added to assess fit ######
  sig2_matrix = t(matrix(resultsList$sigma2, nrow=length(thresholds), byrow=T))
  cond_matrix = t(matrix(resultsList$cond, nrow=length(thresholds), byrow=T))
  int_matrix = t(matrix(resultsList$intercept, nrow=length(thresholds), byrow=T))
  ##############################

  # apply Benjamini-Hochberg correction to each distance threshold
  pvals_adj_matrix = apply(pvals_matrix,2,p.adjust, method="BH")

  # perform t-test (tests mean differences) and McNemar's test (tests number of significant features) between p-values list of different thresholds
  ttestvals = c()
  mcnvals = c()
  comparisons = c()
  maxthreshold <- NA #added for when we want threshold 1 only

  if(length(thresholds) > 1){ #added for when we want threshold 1 only
    for(i in 1:(length(thresholds)-1)) {
      ttestvals = c(ttestvals, t.test(pvals_matrix[,i], pvals_matrix[,i+1], paired=T)$p.value)

      #changed because there were times when this would throw an error due to all p-values being significant (or not)
      mcnvals = c(mcnvals, mcnemar.test(factor(pvals_matrix[,i] <= p_cutoff, levels = c("FALSE", "TRUE")), factor(pvals_matrix[,i+1] <= p_cutoff, levels = c("FALSE", "TRUE")))$p.value)
      comparisons = c(comparisons, paste0(i,"-",i+1))
    }

    disttests = data.frame(ttestvals, mcnvals, row.names=comparisons)
    disttests = disttests > p_cutoff

    # get max. threshold where both test show no significant changes
    maxthreshold = which(apply(disttests, 1, sum)>1)[1]
  }else{
    maxthreshold <- thresholds
  }
  #added because there were times when this would throw an error due to all p-values being significant (or not)
  if(is.na(maxthreshold)) maxthreshold <- 1

  # seek minimum Akaike information criterion of each model till max. treshold
  #added drop = F for situations where maxthreshold = 1
  # had to edit in cases where the model doesn't converge for the max threshold, i.e. aic = NA and minimum is integer(0)
  aic_minima = apply(aic_matrix[,1:maxthreshold,drop=F], 1, function(x){
    if(all(is.na(x))){
      NA
    }else{
      which.min(x)
    }})

  # select CAR model p-value based on minimum AIC and adjust list by Benjamini- Hochberg
  #car_pvals = p.adjust(pvals_matrix[aic_minima], method="BH") #original code doesnt seem right...
  car_pvals_unadjusted =
    sapply(1:nrow(pvals_matrix), function(x){
      pvals_matrix[x, aic_minima[x]]
    })

  car_pvals = p.adjust(
    sapply(1:nrow(pvals_matrix), function(x){
      pvals_matrix[x, aic_minima[x]]
    }), method = "BH")

  ###### added to examine fit ######
  sigma2fit =
    sapply(1:nrow(sig2_matrix), function(x){
      sig2_matrix[x, aic_minima[x]]
    })

  intfit =
    sapply(1:nrow(int_matrix), function(x){
      int_matrix[x, aic_minima[x]]
    })

  condfit =
    sapply(1:nrow(cond_matrix), function(x){
      cond_matrix[x, aic_minima[x]]
    })
  ###########################

  # sum up results in table

  res <- data.frame(mz = mz, ttest=resultsList$padj_ttestpvals,
                    AIC_minimum = aic_minima,
                    CAR_AIC_min_pvalues = car_pvals,
                    CAR_AIC_min_pvalues_unadj = car_pvals_unadjusted,
                    CAR_d_max_pvalues = pvals_adj_matrix[, maxthreshold],
                    sig2 = sigma2fit,
                    int = intfit,
                    cond = condfit
  )

  return(res)

}


