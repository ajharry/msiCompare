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


cass <- function(msset, roiFactor, logscale = T, thresholds = 1:10 ){
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

return(list = list(aics = aics,
                morans_i = morans_i,
                morans_pvals = morans_pvals,
                pvals = pvals,
                padj_ttestpvals = padj_ttestpvals,
                sigma2 = sigma2,
                cond = condeffect,
                intercept = intercept))
}

