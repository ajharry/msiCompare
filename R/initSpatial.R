#' @title initSpatial
#' @description init
#' @param conditionNames condnames
#' @return init
#' @export
#'


initializeSpatial <- function(conditionNames, conditionOfInterest,
                              coord, type.neighbor, radius.neighbor, maxdist.neighbor,
                              nsl, sample){

j <- 1

m <- vector("list", length(nsl))
Wtrip <- vector("list", length(nsl))
Wbegfin <- vector("list", length(nsl))

for(l in conditionNames){
  ind_cond <- conditionOfInterest == l #pixels from condition l

  #####################################################
  ##################### Initialize W ##################
  #####################################################
  #### Create adjacency matrix for pixels from this condition
  assign(paste("W",  l, sep="_"), adj.grid(coords = coord[ind_cond,],
                                           type = type.neighbor,
                                           radius = radius.neighbor,
                                           max.dist = maxdist.neighbor,
                                           sample = sample[ind_cond])+0)
  #### number of neighbors for each pixel
  m[[j]] <- rowSums(get(paste("W",  l, sep="_")))
  names(m)[j] <- paste("m",  l, sep="_")
  ##### number of pixels from this condition and sample pair


  ###################  ###################  ###################
  ####### Stuff for CARBayes

  Wtrip[[j]] <- triplet(as.spam(get(paste("W", l, sep="_"))))
  Wtrip[[j]] <- cbind(Wtrip[[j]]$indices, Wtrip[[j]]$values)
  names(Wtrip)[j] <- paste("Wtrip", l, sep="_")

  Temp <-  array(NA, c(nsl[j], 2))
  temp <- 1


  for(i in 1:nsl[j])
  {
    Temp[i, ] <- c(temp, (temp + m[[j]][i]-1))
    temp <- temp + m[[j]][i]
  }

  Wbegfin[[j]] <- Temp
  names(Wbegfin)[j] <- paste("Wbegfin", l, sep="_")
  rm(Temp)



  #####################################################
  #####################################################
  j <- j+1

}

  return(list(m = m, Wtrip = Wtrip, Wbegfin = Wbegfin))
}
