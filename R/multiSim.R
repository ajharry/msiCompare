#' @title sim multiple images
#' @description sim multiple images
#' @param sampleVar asdfadsf
#' @return simSet
#' @import Cardinal
#' @import mvtnorm
#' @import spam
#' @export
#'
simMulti <- function(
sampleVar = .05,
reps = 1,
diff = log2(1.5),
tau2 = 1,
sig2 = .1,
seed = 100,
size1 = 7,
size2 = size1^2,
numHealthy = 2,
numDisease = 2){

  set.seed(seed)

numSample <- numHealthy + numDisease


samples <-sampleCAR(condDiff = 0,
          coord = expand.grid(x=1:size1, y=1:size1),
          pattern = rep(1, size2),
          sig2 = sig2, tau2 = tau2,
          rho = .9999, nrep=numSample*reps,
          save = F, randomSeed = seed)



sampleNames <- rep(c(paste0("H", numHealthy), paste0("D", numDisease)), each = size2)

diagnosis <- rep(c("Healthy", "Disease"), times = c(size2*numHealthy, size2*numDisease))

coord <- expand.grid(x=1:size1, y=1:size1)

for(i in 1:numSample){
assign(paste0("sample", i), MSImageSet(spectra = matrix(samples[(i*reps-(reps-1)):(i*reps),], ncol = size2),
                                       coord = expand.grid(x=1:size1, y=1:size1),
                                       pixelData = IAnnotatedDataFrame(
                                         data=cbind(coord, sample = rep(paste0("sample",i), size2)),
                                         varMetadata=data.frame(labelType=c("dim","dim","sample"))
                                         )))
  if(i == 2){
    simSet <- combine(sample1, sample2)
  }else if(i > 2){
    simSet <- combine(simSet, get(paste0("sample", i)))
  }

}

rm(list = paste0("sample", 1:numSample))

pData(simSet)$diagnosis <- diagnosis
#image(simSet, feature=1, layout = c(numSample,1))


spectra(simSet)[,simSet$diagnosis == "Healthy"] <- spectra(simSet)[,simSet$diagnosis == "Healthy"] + diff

if(sampleVar > 0){
  sampleEffect <- rnorm(numSample,mean = 0,sd = sqrt(sampleVar))
}else{
  sampleEffect <- rep(0, numSample)
}

for(i in 1:numSample){
  sam <- sampleNames(simSet)[i]
  spectra(simSet)[,simSet$sample == sam] <- spectra(simSet)[,simSet$sample == sam] + sampleEffect[i]
}

return(list(simSet = simSet,
      sampleVar = sampleVar,
       reps = reps,
       diff = diff,
       tau2 = tau2,
       sig2 = sig2,
       seed = seed,
       size1 = size1,
       size2 = size2,
       numHealthy = numHealthy,
       numDisease = numDisease))
}
