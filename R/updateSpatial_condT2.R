#' @title update spatial vector
#' @description update spatial vector
#' @param Wtrip adfadfa
#' @return phiVec
#' @import CARBayes
#' @export
#'
updateSpatial_condT2 <- function(Wtrip,
                               Wbegfin,
                               m,
                               nsl,
                               phiVec,
                               tau2,
                               rho=1,
                               eps_m.var,
                               offset.phi,
                               tauVar.a,
                               tauVar.b,
                               proper = F,
                               sample){
  ###################  ###################  ###################
  ########  get update of spatial effects
  ###################  ###################  ###################

  sample <- factor(sample)

  ###### Debugging ######
  #phiVecold <- phiVec_m
  ###### Debugging ######


  phi <-  CARBayes:::gaussiancarupdate(Wtriplet=Wtrip,
                                       Wbegfin=Wbegfin,
                                       Wtripletsum = m,
                                       nsites= nsl,
                                       phi=phiVec, #
                                       tau2=tau2, #
                                       rho=rho ,
                                       nu2=eps_m.var,
                                       offset=offset.phi,
                                       missind = rep(1,nsl) #0 if observation missing, 1 else
  )


  ##### center within samples
  for(sLevel in levels(sample)){
    s <- sample == sLevel
    phi[s] <- phi[s] - mean(phi[s])
  }


  ###################  ###################  ###################
  ######### get update of spatial effects variance
  ###################  ###################  ###################

  temp2 <- CARBayes:::quadform(Wtrip,
                               m,
                               nrow(Wtrip),
                               nsl,
                               phi,
                               phi,
                               rho)


  tau2.posterior.scale <- temp2 + tauVar.b #.5 factor already in quadform function



  if(proper){ #proper
    tau2 <- 1 / rgamma(1, tauVar.a+nsl/2, rate = tau2.posterior.scale)

  }else{#improper
    tau2 <- 1 / rgamma(1, tauVar.a+(nsl-length(levels(sample)))/2, rate = tau2.posterior.scale)
  }

  return(list(phi = phi, tau2 = tau2))

}
