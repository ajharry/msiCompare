#' @title Fit hiearchical spatial model to MSI data, with improved mixing
#' @description compareMSI_hc is used to fit a hiearchical Bayesian spatial model to MSI data using a Gibbs Sampler MCMC approach. The model is fit separately for each m/z feature. This version is distinct from compareMSI because it uses hiearchical centering to improve mixing of the MCMC. This function will be called by msiCompare if the MSImageSet is a multi tissue experiment with all tissues coming from the same donor or all from different donors.
#' @param msset an object of class "MSImageSet"
#' @param conditionOfInterest a vector or factor giving the level of the condition of interest for each pixel in msset
#' @param feature the index of the m/z features for which the model should be fit
#' @param nsim number of desired MCMC samples
#' @param burnin number of MCMC samples to discard
#' @param trace logical, should the full list of MCMC samples be returned for each variable?
#' @param piPrior prior probability of differential abundance
#' @param seed random seed
#' @param logbase2 logical, should the intensities be log transformed?
#' @param coord data fram of coordinates of the MSImageSet, with columns 'x' and 'y'
#' @param type.neighbor neighborhood type (see adj.grid)
#' @param radius.neighbor desired neighborhood radius if neighborhood type 'radius' is selected (see adj.grid)
#' @param maxdist.neighbor maximum distance for locations to be considered neighbors if neighborhood type 'max.dist' is selected (see adj.grid)
#' @param spInit optional, provide precomputed spatial information from output of intializeSpatial
#' @param bioRep optional, vector or factor giving the individual/donor to which pixel in the msset belongs
#' @param techRep vector or factor giving the tissue to which each pixel in the msset belongs
#' @param beta0 prior mean of baseline effect
#' @param prec0 prior variance of baseline effect
#' @param precAlpha0 prior mean of condition 2 effect
#' @param a0_eps shape parameter for measurment error precision hyperprior
#' @param a0_bio shape parameter for biological replicate error precision hyperprior
#' @param b0_eps scale parameter for measurment error precision hyperprior
#' @param b0_bio scale parameter for biological replicate error precision hyperprior
#' @param a0_tec shape parameter for sample to sample error precision hyperprior
#' @param b0_tec scale parameter for sample to sample error precision hyperprior
#' @param a0_sp shape parameter for spatial precision hyperprior
#' @param b0_sp scale parameter for spatial precision hyperprior
#' @param rd ratio of spike variance to slab variance for condition 2 effect
#' @return res
#' @import mvtnorm
#' @import lme4
#' @import spam
#' @import coda
#' @export
#'

compareMSI_zeros <- function(msset,conditionOfInterest,
                          feature, nsim=5000, burnin = 2500, trace = T,
                          piPrior = .1, seed = 1, logbase2 = F, coord = NULL,
                          type.neighbor = "radius", radius.neighbor = 1, maxdist.neighbor = NULL,
                          spInit = NULL,
                          bioRep = NULL,
                          techRep,
                          beta0 = 0, # Prior Mean for beta, only allow intercept
                          prec0 = .01, # Prior Precision Matrix of beta (vague)  (only allow intercept)
                          precAlpha0 = .01, #Prior Precision of slab (value of condition effect if it is not zero)
                          a0_eps=.001, b0_eps=.001,			# Hyperprior for tau (1/eps.var)
                          a0_bio=.001, b0_bio=.001,			# Hyperprior for taubio
                          a0_tec=.001, b0_tec=.001,			# Hyperprior for tautec
                          a0_sp=.001, b0_sp=.001,			# Hyperprior for tau.spatial
                          rd = .00001, # ratio of varSpike/varSlab
                          dropZeros = T #should we treat zeros as missing values?
){

  techRep <- factor(techRep) #factor with different levels for each tissue (like "sample" before)
  n_tec <- length(levels(techRep)) #the number of distinct tissues
  nis_tec <- sapply(levels(techRep), function(x) sum(techRep == x)) #number of pixels in each tissue


  if(n_tec == 1){
    return(compareMSI_zerosSingle(msset,conditionOfInterest,
                      feature, nsim, burnin, trace,
                      piPrior, seed, logbase2, coord,
                      type.neighbor, radius.neighbor, maxdist.neighbor,
                      spInit,
                      bioRep,
                      techRep,
                      beta0, # Prior Mean for beta, only allow intercept
                      prec0, # Prior Precision Matrix of beta (vague)  (only allow intercept)
                      precAlpha0, #Prior Precision of slab (value of condition effect if it is not zero)
                      a0_eps, b0_eps,			# Hyperprior for tau (1/eps.var)
                      a0_bio, b0_bio,			# Hyperprior for taubio
                      a0_tec, b0_tec,			# Hyperprior for tautec
                      a0_sp, b0_sp,			# Hyperprior for tau.spatial
                      rd, # ratio of varSpike/varSlab
                      dropZero = T
    ))
  }


  set.seed(seed) #random seed

  if(is.null(coord)){
    coord <- coord(msset)
  }


  conditionOfInterest <- factor(conditionOfInterest) # make the condition labels a factor in case it is a character vector
  conditionNames <- levels(conditionOfInterest) #create vector of condition names
  nCond <- length(conditionNames) #how many conditions there are

  if(!is.null(bioRep)){
    bioRep <- factor(bioRep) #factor with different levels for each biological unit
    n_bio <- length(levels(bioRep)) #the number of distinct biological units
    nis_bio <- sapply(levels(bioRep), function(x) sum(bioRep == x)) #number of pixels in each biological unit
  }


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


  ####### check for islands of unconnected pixels within tissues
  W <- adj.grid(coord, sample = techRep,
                type = type.neighbor,
                radius = radius.neighbor,
                max.dist = maxdist.neighbor)+0
  W.list<- mat2listw(W)
  W.nb <- W.list$neighbours
  W.islands <- n.comp.nb(W.nb)
  islands <- W.islands$comp.id
  rm(W, W.list, W.nb, W.islands)

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
for(f in feature){
    print(paste0("Feature ", f, " of ", length(feature)))
    time <- system.time({ #time the overall model fits
      y <- spectra(msset)[f,]

      if(logbase2){ #do log transformation if necessary
 y <- log2(y)
      }




      ####################################################################################################
      ################################### Initialize variables  #####################################
      ####################################################################################################

      lm <- lm(y~X+X1) ## fit linear model to get reasonable starting values
      coef <- coef(lm)[-2]
      tau<-1				# technical error precision
      eps_m.var <- 1/tau #technical error variance
      if(!is.null(bioRep)) b_bio <-rep(0,n_bio)			# Random effects bio unit
      tau_tec<- tau_bio<-1				# Random Effects precision
      beta <- coef[1:k] #initial value of intercept and covariates
      alpha <- coef[k+1] #initial value of condition effect
      b_tec <- rep(0, numTissueCond1+numTissueCond2)

      if(!is.null(bioRep)) Z_bio<-as.spam(matrix(0,N,n_bio))	# Random effect design used for updating b_bio
      Z_tec<-as.spam(matrix(0,N,n_tec))	# Random effect design used for updating b_tec
      for(i in 1:n_tec) Z_tec[as.numeric(techRep)==i,i]<-1
      if(!is.null(bioRep)){ for(i in 1:n_bio){Z_bio[as.numeric(bioRep)==i,i]<-1}}

      xb <-  X%*%beta
      x1a <- X1 %*% alpha
      zb_bio <- zb_tec <-  rep(0, N)
      gamma <- 1 # initiate condition effect as nonzero
      tauVar <- rep(1,numSpatialParams) # spatial variances


      #################
      # Store Results #
      #################

      Betas<-matrix(0,nsim,k)	# Fixed Effects
      spVar<-matrix(0,nsim,numSpatialParams)
      taus<- taus_tec<-gammas <- rep(0,nsim)
      if(!is.null(bioRep)) taus_bio <- rep(0,nsim)
      Condition <- Condition0 <- Condition1 <- rep(NA,nsim)	# Error Precision Parms

      ###############################
      # Fixed Posterior Hyperparms 	#
      #    for tau and tautech		#
      ###############################
      d<-a0_eps+N/2
      if(!is.null(bioRep)) nu_bio<-a0_bio+n_bio/2
      nu_tec<-a0_tec+n_tec/2

      ####################################################################################################
      ######################################## THE GIBBS SAMPLER  ########################################
      ####################################################################################################
      for (i in 1:nsim) { #this is an iterative method, nsim is the number of iterations


        ############# First level tissue/sample effect ##########
        sc2 <- 1
        for(sc in 1:nrow(sampCond)){
          ncs <- sampCond$numPix[sc]  #number of pixels in condition i, sample j
          vb_tec <- 1/(ncs/eps_m.var + tau_tec)
          mb_tec <-vb_tec*(((beta+alpha*sampCond$condition[sc])*tau_tec) + (sum((y-zb_bio-phiVec_m)[conditionVec == sampCond$condition[sc] &
                                                                                                      techRep == sampCond$sample[sc]])/eps_m.var))
          b_tec[sc2] <- rnorm(1,mean = mb_tec, sd = sqrt(vb_tec))

          zb_tec[conditionVec == sampCond$condition[sc] & techRep == sampCond$sample[sc]] <- b_tec[sc2]

          names(b_tec)[sc2] <- paste(sampCond$sample[sc], sampCond$condition[sc], sep = "_")
          sc2 <- sc2+1
        }
        ###################################################

        ################### First level: subject/bio effect ###################

        if(!is.null(bioRep)){
          # Update the biological Replicate effect
          vb_bio<-1/(tau_bio+nis_bio*tau)
          mb_bio<-vb_bio*(tau*t(Z_bio)%*%(y-zb_tec-phiVec_m))
          b_bio<-rnorm(n_bio,mb_bio,sqrt(vb_bio))
          zb_bio<-Z_bio%*%b_bio
        }else{
          zb_bio <- rep(0, N)
        }

        #####################################################################


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
        g<-b0_eps+crossprod(y-zb_tec-zb_bio-phiVec_m,
                        y-zb_tec-zb_bio-phiVec_m)/2
        taus[i]<-tau<-rgamma(1,d,g)
        eps_m.var <- 1/tau
        #####################################################################

        ################### ------  level: biological var error precision ###################

        if(!is.null(bioRep)){
          m_bio<-c(b0_bio+crossprod(b_bio,b_bio)/2)
          taus_bio[i]<-tau_bio<-rgamma(1,nu_bio,m_bio)
        }

        #####################################################################

        ################### Third level: sample-to-sample error precision ###################

        m_tec<-c(b0_tec+sum(unique(zb_tec-x1a-xb)^2)/2)
        taus_tec[i]<-tau_tec<-rgamma(1,nu_tec,m_tec)


        #####################################################################


        offset.phi <- (y-zb_tec-zb_bio) / eps_m.var

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
            tauVar.a = a0_sp,
            tauVar.b = b0_sp,
            sample = techRep[ind_cond],
            islands = islands[ind_cond]
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

      if(!is.null(bioRep)){
        msigma.b2_bio<-mean(1/taus_bio[(burnin):nsim])
      }else{
        msigma.b2_bio <- NA
      }

      if(n_tec >1){
        msigma.b2_tec<-mean(1/taus_tec[(burnin):nsim])
      }else{
        msigma.b2_tec<- NA
      }

      msigma.t2<-apply(spVar[(burnin):nsim,, drop = F],2,mean)
      gam <- mean(gammas[burnin:nsim])
      malpha <- mean(Condition[burnin:nsim])
      malpha1 <- mean(Condition1[burnin:nsim], na.rm = T)
      malpha0 <- mean(Condition0[burnin:nsim], na.rm = T)


      if(!is.null(bioRep)){ #both bio and tech replicates
        ess <- effectiveSize(mcmc(cbind(beta_trace = c(Betas)[burnin:nsim],
                                        cond_trace = Condition[burnin:nsim],
                                        sig2_trace = 1/taus[burnin:nsim],
                                        sig2tec_trace = 1/taus_tec[burnin:nsim],
                                        sig2bio_trace = 1/taus_bio[burnin:nsim],
                                        tau2_trace1 = spVar[burnin:nsim,1],
                                        tau2_trace2 = spVar[burnin:nsim,2],
                                        gamma_trace = gammas[burnin:nsim])))
      }else{
        if(n_tec > 1){ #only tech replicates
          ess <- effectiveSize(mcmc(cbind(beta_trace = c(Betas)[burnin:nsim],
                                          cond_trace = Condition[burnin:nsim],
                                          sig2_trace = 1/taus[burnin:nsim],
                                          sig2tec_trace = 1/taus_tec[burnin:nsim],
                                          tau2_trace1 = spVar[burnin:nsim,1],
                                          tau2_trace2 = spVar[burnin:nsim,2],
                                          gamma_trace = gammas[burnin:nsim])))
        }else{ #single sample
          ess <- effectiveSize(mcmc(cbind(beta_trace = c(Betas)[burnin:nsim],
                                          cond_trace = Condition[burnin:nsim],
                                          sig2_trace = 1/taus[burnin:nsim],
                                          tau2_trace1 = spVar[burnin:nsim,1],
                                          tau2_trace2 = spVar[burnin:nsim,2],
                                          gamma_trace = gammas[burnin:nsim])))
        }
      }


    }) #time


    if(trace){ #return trace info
      if(is.null(bioRep)){
        if(n_tec >1){ #there are no bio reps, only technical
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
        }else{ #single sample, no reps of any kind
          res[[feat]] <-list(
            beta = mbeta,
            cond = malpha,
            cond0 = malpha0,
            cond1 = malpha1,
            sig2 = msigma.e2,
            tau2 = msigma.t2,
            gamma = gam,
            ess = ess,
            trace = mcmc(cbind(beta_trace = c(Betas),
                               cond_trace = Condition,
                               cond1_trace = Condition1,
                               cond0_trace = Condition0,
                               sig2_trace = 1/taus,
                               tau2_trace1 = spVar[,1],
                               tau2_trace2 = spVar[,2],
                               gamma_trace = gammas)),
            time = time
          )
        }
      }else{ #both bio and technical reps
        res[[feat]] <-list(
          beta = mbeta,
          cond = malpha,
          cond0 = malpha0,
          cond1 = malpha1,
          sig2 = msigma.e2,
          sig2tec = msigma.b2_tec,
          sig2bio = msigma.b2_bio,
          tau2 = msigma.t2,
          gamma = gam,
          ess = ess,
          trace = mcmc(cbind(beta_trace = c(Betas),
                             cond_trace = Condition,
                             cond1_trace = Condition1,
                             cond0_trace = Condition0,
                             sig2_trace = 1/taus,
                             sig2tec_trace = 1/taus_tec,
                             sig2bio_trace = 1/taus_bio,
                             tau2_trace1 = spVar[,1],
                             tau2_trace2 = spVar[,2],
                             gamma_trace = gammas)),
          time = time
        )
      }
    }else{ #dont return trace info
      res[[feat]] <-list(
        beta = mbeta,
        cond = malpha,
        cond0 = malpha0,
        cond1 = malpha1,
        sig2 = msigma.e2,
        sig2tec = msigma.b2_tec, #will return NA
        sig2bio = msigma.b2_bio, #will return NA
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



