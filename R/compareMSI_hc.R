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

compareMSI_hc <- function(msset,conditionOfInterest,
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

  techRep <- factor(techRep) #factor with different levels for each tissue (like "sample" before)
  n_tec <- length(levels(techRep)) #the number of distinct tissues
  nis_tec <- sapply(levels(techRep), function(x) sum(techRep == x)) #number of pixels in each tissue


  if(n_tec == 1 | !is.null(bioRep)){
    return(compareMSI2(msset,conditionOfInterest,
                       feature, nsim, burnin, trace,
                       piPrior, seed, logbase2, coord,
                       type.neighbor, radius.neighbor, maxdist.neighbor,
                       spInit,
                       bioRep,
                       techRep,
                       beta0, # Prior Mean for beta, only allow intercept
                       prec0, # Prior Precision Matrix of beta (vague)  (only allow intercept)
                       precAlpha0, #Prior Precision of slab (value of condition effect if it is not zero)
                       d0, g0,			# Hyperpriors for tau, taubio, tautec
                       rd # ratio of varSpike/varSlab
    ))
  }


  set.seed(seed) #random seed

  if(is.null(coord)){
    coord <- coord(msset)
  }


  conditionOfInterest <- factor(conditionOfInterest) # make the condition labels a factor in case it is a character vector
  conditionNames <- levels(conditionOfInterest) #create vector of condition names
  nCond <- length(conditionNames) #how many conditions there are


  N <- nrow(coord) # how many pixels there are

  #### ONLY IF THERE ARE TWO CONDITIONS #####
  conditionVec <- ifelse(conditionOfInterest == conditionNames[1], 0, 1) #vector converts condition names from characters to numeric
  numCond2 <- sum(conditionVec == 1) #number of pixels from condition two
  numCond1 <- sum(conditionVec == 0) #number of pixels from condition 1
  ###########################################

  X <- matrix(rep(1, N), ncol = 1) #design matrix for intercept and covariates #currently set to intercept only
  X1 <- matrix(conditionVec, ncol = 1)  #design matrix without condition effect


  k <-ncol(X) #number of covariates, including intercept
  res <- list() #list that will hold results


  ####################################################################################################
  ############## Obtain neighborhood matrices for each combination of sample and condition ###########
  ####################################################################################################

  numSpatialParams <- 2 #number of spatial parameters to estimate. this will be the number unique of conditions being compared
  nsl <- c(sum(conditionVec == 0), sum(conditionVec == 1)) ##### vector of number of pixels from each condition
  names(nsl) <- conditionNames

  if(is.null(spInit)){
    print("Initializing spatial components...")
    sptime <- system.time({
      spInit <- initializeSpatial(conditionNames= conditionNames, conditionOfInterest = conditionOfInterest,
                                  coord = coord, type.neighbor = type.neighbor, radius.neighbor = radius.neighbor,
                                  maxdist.neighbor = maxdist.neighbor, nsl = nsl,
                                  sample = techRep)
    })
    print(paste0("...Initialization done in ", sptime['elapsed'], " seconds."))
  }else{
    print("Spatial components provided, no need for initialization.")
  }

  for(i in 1:length(spInit)){
    for(j in 1:length(spInit[[i]])){
      assign(names(spInit[[i]][j]), spInit[[i]][[j]])
    }
  }
  rm(spInit)


  ####################################################################################################
  ####################################################################################################
  ####################################################################################################


  ########## Get number of pixels for each sample/condition combination ########
  sampCond <- data.frame(sample = numeric(n_tec*nCond),
                         condition = numeric(n_tec*nCond),
                         numPix = numeric(n_tec*nCond))
  indsc <- 1
  for(smp in levels(techRep)){
    for(cd in unique(conditionVec)){
      sampCond$sample[indsc] <- smp
      sampCond$condition[indsc] <- cd
      sampCond$numPix[indsc] <- sum(techRep == smp & conditionVec == cd)
      indsc <- indsc + 1
    }
  }

  sampCond <- sampCond[sampCond$numPix > 0, ]

  ##### Number of tissues in each condition
  numTissueCond1 <- sum(sampCond$condition == 0)
  numTissueCond2 <- sum(sampCond$condition == 1)
  ###############################################################################


  phiVec_m <- rep(0, N) #initialize trace vector for spatial effects
  feat <- 1 #initialize feature index


  ####################################################################################################
  ##################################### Fit model feature by feature #################################
  ####################################################################################################
  minNonZero <- min(spectra(msset)[spectra(msset) != 0])/100
  for(f in feature){
    print(paste0("Feature ", f, " of ", length(feature)))
    time <- system.time({ #time the overall model fits
      y <- spectra(msset)[f,]

      if(logbase2){ #do log transformation if necessary
        y[y==0] <- minNonZero  #zeros in the image will cause problems if a log transformation is required. add a small number to the zeroes.
        y <- log2(y)
      }




      ####################################################################################################
      ################################### Initialize variables  #####################################
      ####################################################################################################

      lm <- lm(y~X+X1) ## fit linear model to get reasonable starting values
      coef <- coef(lm)[-2]
      tau<-1				# technical error precision
      eps_m.var <- 1/tau #technical error variance

      tau_tec<-1				# Random Effects precision
      beta <- coef[1:k] #initial value of intercept and covariates
      alpha <- coef[k+1] #initial value of condition effect
      b_tec <- rep(0, numTissueCond1+numTissueCond2)

      Z_tec<-as.spam(matrix(0,N,n_tec))	# Random effect design used for updating b_tec
      for(i in 1:n_tec) Z_tec[as.numeric(techRep)==i,i]<-1
      xb <-  X%*%beta
      x1a <- X1 %*% alpha
      zb_tec <-  rep(0, N)
      gamma <- 1 # initiate condition effect as nonzero
      tauVar <- rep(1,numSpatialParams) # spatial variances


      #################
      # Store Results #
      #################

      Betas<-matrix(0,nsim,k)	# Fixed Effects
      spVar<-matrix(0,nsim,numSpatialParams)
      taus<- taus_tec<-gammas <- rep(0,nsim)
      Condition <- Condition0 <- Condition1 <- rep(NA,nsim)	# Error Precision Parms

      ###############################
      # Fixed Posterior Hyperparms 	#
      #    for tau and tautech		#
      ###############################
      d<-d0+N/2
      nu_tec<-d0+n_tec/2

      ####################################################################################################
      ######################################## THE GIBBS SAMPLER  ########################################
      ####################################################################################################
      for (i in 1:nsim) { #this is an iterative method, nsim is the number of iterations


        ############# First level sample effect bij ##########
        sc2 <- 1
        for(sc in 1:nrow(sampCond)){
            ncs <- sampCond$numPix[sc]
            vb_tec <- 1/(ncs/eps_m.var + tau_tec)
            mb_tec <-vb_tec*(((beta+alpha*sampCond$condition[sc])*tau_tec) + (sum((y-phiVec_m)[conditionVec == sampCond$condition[sc] &
                                                                                             techRep == sampCond$sample[sc]])/eps_m.var))
            b_tec[sc2] <- rnorm(1,mean = mb_tec, sd = sqrt(vb_tec))

            zb_tec[conditionVec == sampCond$condition[sc] & techRep == sampCond$sample[sc]] <- b_tec[sc2]

            names(b_tec)[sc2] <- paste(sampCond$sample[sc], sampCond$condition[sc], sep = "_")
            sc2 <- sc2+1
        }
        ###################################################

        ################### Second level: baseline effect ###################
        resbeta <- sum(b_tec[sampCond$condition == 0]) #residuals for pixels in first condition only

        vbeta<- 1/(prec0+numTissueCond1*tau_tec)
        mbeta<-vbeta*(prec0*beta0 + resbeta*tau_tec)
        beta <- rnorm(n=1, mean = mbeta, sd = sqrt(vbeta))
        xb <-  X*beta
        Betas[i,]<- beta
        #####################################################################

        ################### Second level: condition effect ###################
        resa <- sum((b_tec-beta)[sampCond$condition == 1]) #residuals for pixels in second condition onlt

        # Update Condition effect
        if(gamma == 1){ #this is the estimate if the condition effect is not zero
          valph <- 1/(numTissueCond2*tau_tec + precAlpha0)
          malph <- valph*(0*precAlpha0 + resa*tau_tec)
          Condition1[i] <- alpha <- rnorm(n = 1, mean = malph, sd = sqrt(valph))
        }else{ #this is the estimate if the condition effect is zero (or very close to it)
          valph <- 1/(numTissueCond2*tau_tec + 1/rd * precAlpha0)
          malph <- valph*(0*(1/rd)*precAlpha0 + resa*tau_tec)
          Condition0[i] <- alpha <- rnorm(n = 1, mean = malph, sd = sqrt(valph))
        }
        x1a <- X1*alpha
          Condition[i] <- alpha
          #####################################################################

          ################### Third level: indicator of differential abundance ###################
          loglik_slab <- dnorm(alpha, mean = 0, sd = sqrt(1/precAlpha0), log = T)
          loglik_spike <- dnorm(alpha, mean = 0 , sd = sqrt(rd*(1/precAlpha0)), log = T)
          pi1Post <-  1/(1 + exp(loglik_spike - loglik_slab)*(1-piPrior)/piPrior )

          gamma <- rbinom(n=1, size = 1, prob = pi1Post)
          gammas[i] <-gamma
          #####################################################################


          ################### Second level: measurement error precision ###################
        g<-g0+crossprod(y-zb_tec-phiVec_m,
                        y-zb_tec-phiVec_m)/2
        taus[i]<-tau<-rgamma(1,d,g)
        eps_m.var <- 1/tau
        #####################################################################

        ################### Third level: sample-to-sample error precision ###################

          m_tec<-c(g0+sum(unique(zb_tec-x1a-xb)^2)/2)
          taus_tec[i]<-tau_tec<-rgamma(1,nu_tec,m_tec)

        offset.phi <- (y-zb_tec) / eps_m.var

        #########################################################
        ########### Update the spatial effects ##################
        #########################################################

        j <- 1
        for(l in conditionNames){
          ind_cond <- conditionOfInterest == l


          offset <- offset.phi[ind_cond]

          phiUpdate <- updateSpatial_condT2(
            Wtrip=get(paste("Wtrip", l, sep="_")),
            Wbegfin=get(paste("Wbegfin", l, sep="_")),
            m = get(paste("m", l, sep="_")),
            nsl= nsl[j],
            phiVec=phiVec_m[ind_cond], #
            tau2=tauVar[j], #
            rho=1,
            eps_m.var =eps_m.var,
            offset.phi =offset,
            tauVar.a = .001,
            tauVar.b = .001,
            sample = techRep[ind_cond]
          )



          phiVec_m[ind_cond] <- phiUpdate$phi
          spVar[i,j] <- tauVar[j] <- phiUpdate$tau2
          j <- j+1

        }


        #########################################################
        #########################################################
        #########################################################
        if (i%%1000==0 || i == 1) print(paste0("MCMC Iteration ", i, " of ", nsim))
      } #On to the next mcmc iteration

      ###########
      # Results #
      ###########

      mbeta<-apply(Betas[(burnin):nsim,, drop = F],2,mean)
      msigma.e2<-mean(1/taus[(burnin):nsim])

      msigma.b2_tec<-mean(1/taus_tec[(burnin):nsim])

      msigma.t2<-apply(spVar[(burnin):nsim,, drop = F],2,mean)
      gam <- mean(gammas[burnin:nsim])
      malpha <- mean(Condition[burnin:nsim])
      malpha1 <- mean(Condition1[burnin:nsim], na.rm = T)
      malpha0 <- mean(Condition0[burnin:nsim], na.rm = T)


        ess <- effectiveSize(mcmc(cbind(beta_trace = c(Betas)[burnin:nsim],
                          cond_trace = Condition[burnin:nsim],
                          cond1_trace = Condition1[burnin:nsim],
                          cond0_trace = Condition0[burnin:nsim],
                          sig2_trace = 1/taus[burnin:nsim],
                          sig2tec_trace = 1/taus_tec[burnin:nsim],
                          tau2_trace1 = spVar[burnin:nsim,1],
                          tau2_trace2 = spVar[burnin:nsim,2],
                          gamma_trace = gammas[burnin:nsim])))

    }) #time


   if(trace & is.null(bioRep)){
      res[[feat]] <-list(
        beta = mbeta,
        cond = malpha,
        cond0 = malpha0,
        cond1 = malpha1,
        sig2 = msigma.e2,
        sig2tec = msigma.b2_tec,
        tau2 = msigma.t2,
        gamma = gam,
        ess = ess,
        trace = mcmc(cbind(beta_trace = c(Betas),
                           cond_trace = Condition,
                           cond1_trace = Condition1,
                           cond0_trace = Condition0,
                           sig2_trace = 1/taus,
                           sig2tec_trace = 1/taus_tec,
                           tau2_trace1 = spVar[,1],
                           tau2_trace2 = spVar[,2],
                           gamma_trace = gammas)),
        time = time
      )
    }else{
      res[[feat]] <-list(
        beta = mbeta,
        cond = malpha,
        cond0 = malpha0,
        cond1 = malpha1,
        sig2 = msigma.e2,
        sig2tec = msigma.b2_tec,
        tau2 = msigma.t2,
        gamma = gam,
        ess = ess,
        time = time
      )
    }

    names(res)[feat] <- paste0("Feature",f)
    feat <- feat + 1


  } #feature


  return(res)
}#function



