#' @title Fit hiearchical spatial model to MSI data
#' @description compareMSI is used to fit a hiearchical Bayesian spatial model to MSI data using a Gibbs Sampler MCMC approach. The model is fit separately for each m/z feature.
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

compareMSI <- function(msset,conditionOfInterest,
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
                       dropZeros = F
){

  set.seed(seed) #random seed

  if(is.null(coord)){
    coord <- coord(msset)
  }



  if(dropZeros){
    res <- vector("list", length(feature))
    fInd <- 1
    seeds <- sample(size = length(feature), x=1:100000)

    for(f in feature){
      msset0 <- msset[f,]

      ######## drop zero pixels ##########
      for(s in levels(factor(techRep))){
        print(paste0("For feature", f, ", dropping ",
                     sum(spectra(msset0) == 0 & techRep == s),
                     " zero pixels from tissue ", s))
      }

      techRep0 <- techRep[c(spectra(msset0) != 0)]
      bioRep0 <- bioRep[c(spectra(msset0) != 0)]
      conditionOfInterest0 <- conditionOfInterest[c(spectra(msset0) != 0)]
      coord0 <- coord[c(spectra(msset0) != 0),]
      msset0 <- msset0[,c(spectra(msset0) != 0)]



      #### drop pixels with no neighbors ####
      W <- adj.grid(coord0, sample = factor(paste0(techRep0, conditionOfInterest0)),
                    type = type.neighbor,
                    radius = radius.neighbor,
                    max.dist = maxdist.neighbor)+0
      lonely <- which(rowSums(W) == 0)

      if(length(lonely) > 0){
        rm(W)
        print(paste0("dropping ", length(lonely), " pixel(s) with no neighbors: "))
        #print(coord0[lonely,])

        coord0 <- coord0[-lonely,]
        msset0 <- msset0[,-lonely]
        techRep0 <- techRep0[-lonely]
        bioRep0 <- bioRep0[-lonely]
        conditionOfInterest0 <- conditionOfInterest0[-lonely]

       if(any(table(techRep0) == 0)) warning("At least one tissue has no non-zero, non-island pixels")
        
        print("Remaining pixels per tissue:")
        print(table(techRep0))

        W <- adj.grid(coord0, sample = techRep0,
                      type = type.neighbor,
                      radius = radius.neighbor,
                      max.dist = maxdist.neighbor)+0
      }else{
        print("All pixels have at least one neighbor.")
      }

      ### Find any disconnected islands of pixels within samples (because spatial effects must be centered within islands) ####
      ### Code thanks to CARBayes ###
      W.list<- mat2listw(W)
      W.nb <- W.list$neighbours
      W.islands <- n.comp.nb(W.nb)
      islands <- W.islands$comp.id
      n.islands <- length(unique(islands))
      rm(W, W.list, W.nb, W.islands)

      if(length(levels(factor(techRep0))) > 1){
        print("Fitting model version: drop zeros, multiple samples, hiearchical centering")
        res[[fInd]] <- compareMSI_zeros(msset0,conditionOfInterest0,
                                        feature=1, nsim, burnin, trace,
                                        piPrior, seeds[fInd], logbase2, coord0,
                                        type.neighbor, radius.neighbor, maxdist.neighbor,
                                        spInit=NULL,
                                        bioRep0,
                                        techRep0,
                                        beta0, # Prior Mean for beta, only allow intercept
                                        prec0, # Prior Precision Matrix of beta (vague)  (only allow intercept)
                                        precAlpha0, #Prior Precision of slab (value of condition effect if it is not zero)
                                        a0_eps, b0_eps,			# Hyperprior for tau (1/eps.var)
                                        a0_bio, b0_bio,			# Hyperprior for taubio
                                        a0_tec, b0_tec,			# Hyperprior for tautec
                                        a0_sp, b0_sp,			# Hyperprior for tau.spatial
                                        rd, # ratio of varSpike/varSlab
                                        dropZeros = T
        )[[1]]
        res[[fInd]]$seed <- seeds[fInd]
      }else{
        print("Fitting model version: drop zeros, single sample")
        res[[fInd]] <- compareMSI_zerosSingle(msset0,conditionOfInterest0,
                                              feature=1, nsim, burnin, trace,
                                              piPrior, seeds[fInd], logbase2, coord0,
                                              type.neighbor, radius.neighbor, maxdist.neighbor,
                                              spInit=NULL,
                                              bioRep0,
                                              techRep0,
                                              beta0, # Prior Mean for beta, only allow intercept
                                              prec0, # Prior Precision Matrix of beta (vague)  (only allow intercept)
                                              precAlpha0, #Prior Precision of slab (value of condition effect if it is not zero)
                                              a0_eps, b0_eps,			# Hyperprior for tau (1/eps.var)
                                              a0_bio, b0_bio,			# Hyperprior for taubio
                                              a0_tec, b0_tec,			# Hyperprior for tautec
                                              a0_sp, b0_sp,			# Hyperprior for tau.spatial
                                              rd, # ratio of varSpike/varSlab
                                              dropZeros = T
        )[[1]]
        res[[fInd]]$seed <- seeds[fInd]
      }
      names(res)[fInd] <- paste0("Feature",f)
      fInd <- fInd + 1
    }
    return(res)

  }else{ #don't drop zeros

    techRep <- factor(techRep) #factor with different levels for each tissue (like "sample" before)
    n_tec <- length(levels(techRep)) #the number of distinct tissues
    nis_tec <- sapply(levels(techRep), function(x) sum(techRep == x)) #number of pixels in each tissue

    if(n_tec > 1){ #do hiearchical centering for multi tissue experiments (only those without subsampling for now)
      print("Fitting model version: replace zeros with small value, multiple samples, hiearchcial centering")
      return(compareMSI_hc_sub(msset,conditionOfInterest,
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
                               rd # ratio of varSpike/varSlab
      ))
    }else{#don't do hiearchical centering if it's a single tissue experiment
      print("Fitting model version: replace zeros with small value, single sample")





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


    ####################################################################################################
    ####################################################################################################
    ####################################################################################################


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
        if(!is.null(bioRep)) b_bio <-rep(0,n_bio)			# Random effects bio unit
        b_tec <-rep(0,n_tec)			# Random effects tech unit (sample/tissue)
        tau_bio <- tau_tec<-1				# Random Effects precision
        beta <- coef[1:k] #initial value of intercept and covariates
        alpha <- coef[k+1] #initial value of condition effect
        if(!is.null(bioRep)) Z_bio<-as.spam(matrix(0,N,n_bio))	# Random effect design used for updating b_bio
        Z_tec<-as.spam(matrix(0,N,n_tec))	# Random effect design used for updating b_tec
        if(!is.null(bioRep)){ for(i in 1:n_bio){Z_bio[as.numeric(bioRep)==i,i]<-1}}
        for(i in 1:n_tec) Z_tec[as.numeric(techRep)==i,i]<-1
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
        #    for tau and taubio tautech		#
        ###############################
        d<-a0_eps+N/2
        if(!is.null(bioRep)) nu_bio<-a0_bio+n_bio/2
        nu_tec<-a0_tec+n_tec/2

        ####################################################################################################
        ######################################## THE GIBBS SAMPLER  ########################################
        ####################################################################################################
        for (i in 1:nsim) { #this is an iterative method, nsim is the number of iterations

          ######## Less General, only allow intercept, but assuring baseline constraint ####################
          resbeta <- sum((y-x1a-zb_bio-zb_tec-phiVec_m)[conditionVec == 0]) #residuals for pixels in first condition only

          vbeta<- 1/(prec0+numCond1/eps_m.var)
          mbeta<-vbeta*(prec0*beta0 + resbeta/eps_m.var)
          beta <- rnorm(n=1, mean = mbeta, sd = sqrt(vbeta))
          xb <-  X*beta
          Betas[i,]<- beta


          # Update intercept and covariates
          # vbeta<-solve(prec0+tau*crossprod(X,X))
          # mbeta<-vbeta%*%(prec0%*%beta0 + tau*crossprod(X,y-x1a-zb_bio-zb_tec-phiVec_m))
          # beta <-c(rmvnorm(1,mbeta,vbeta))
          # xb <-  X%*%beta
          # Betas[i,]<- beta


          resa <- sum((y-xb-zb_bio-zb_tec-phiVec_m)[conditionVec == 1]) #residuals for pixels in second condition onlt

          # Update Condition effect
          if(gamma == 1){ #this is the estimate if the condition effect is not zero
            valph <- 1/(numCond2/eps_m.var + precAlpha0)
            malph <- valph*resa/eps_m.var
            Condition1[i] <- alpha <- rnorm(n = 1, mean = malph, sd = sqrt(valph))
          }else{ #this is the estimate if the condition effect is zero (or very close to it)
            valph <- 1/(numCond2/eps_m.var + 1/rd * precAlpha0)
            malph <- valph*resa/eps_m.var
            Condition0[i] <- alpha <- rnorm(n = 1, mean = malph, sd = sqrt(valph))
          }

          Condition[i] <- alpha
          x1a <- X1 %*% alpha

          # update indicator of differential abundance
          loglik_slab <- dnorm(alpha, mean = 0, sd = sqrt(1/precAlpha0), log = T)
          loglik_spike <- dnorm(alpha, mean = 0 , sd = sqrt(rd/precAlpha0), log = T)
          pi1Post <-  1/(1 + exp(loglik_spike - loglik_slab)*(1-piPrior)/piPrior )


          gamma <- rbinom(n=1, size = 1, prob = pi1Post)
          gammas[i] <-gamma


          if(!is.null(bioRep)){
            # Update the biological Replicate effect
            vb_bio<-1/(tau_bio+nis_bio*tau)
            mb_bio<-vb_bio*(tau*t(Z_bio)%*%(y-x1a-xb-zb_tec-phiVec_m))
            b_bio<-rnorm(n_bio,mb_bio,sqrt(vb_bio))
            zb_bio<-Z_bio%*%b_bio
          }else{
            zb_bio <- rep(0, N)
          }

          if(n_tec > 1){
            vb_tec<-1/(tau_tec+nis_tec*tau)
            mb_tec<-vb_tec*(tau*t(Z_tec)%*%(y-x1a-xb-zb_bio-phiVec_m))
            b_tec<-rnorm(n_tec,mb_tec,sqrt(vb_tec))
            zb_tec<-Z_tec%*%b_tec
          }else{
            zb_tec <- rep(0, N)
          }




          # Update the measurment error precision

          g<-b0_eps+crossprod(y-xb-x1a-zb_bio-zb_tec-phiVec_m,
                          y-xb-x1a-zb_bio-zb_tec-phiVec_m)/2
          taus[i]<-tau<-rgamma(1,d,g)
          eps_m.var <- 1/tau

          if(n_tec > 1){
            # Update the precision of the technical replicate effect
            m_tec<-c(b0_tec+crossprod(b_tec,b_tec)/2)
            taus_tec[i]<-tau_tec<-rgamma(1,nu_tec,m_tec)
          }else{
            taus_tec[i]<-tau_tec<- NA
          }

          if(!is.null(bioRep)){
            m_bio<-c(b0_bio+crossprod(b_bio,b_bio)/2)
            taus_bio[i]<-tau_bio<-rgamma(1,nu_bio,m_bio)
          }

          offset.phi <- (y-xb-x1a-zb_bio-zb_tec) / eps_m.var

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
              islands = NULL
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


        if(trace & !is.null(bioRep)){
          ess <- effectiveSize(mcmc(cbind(beta_trace = c(Betas)[burnin:nsim],
                                          cond_trace = Condition[burnin:nsim],
                                          sig2_trace = 1/taus[burnin:nsim],
                                          sig2tec_trace = 1/taus_tec[burnin:nsim],
                                          sig2bio_trace = 1/taus_bio[burnin:nsim],
                                          tau2_trace1 = spVar[burnin:nsim,1],
                                          tau2_trace2 = spVar[burnin:nsim,2],
                                          gamma_trace = gammas[burnin:nsim])))
        }else{
          if(n_tec > 1){
            ess <- effectiveSize(mcmc(cbind(beta_trace = c(Betas)[burnin:nsim],
                                            cond_trace = Condition[burnin:nsim],
                                            sig2_trace = 1/taus[burnin:nsim],
                                            sig2tec_trace = 1/taus_tec[burnin:nsim],
                                            tau2_trace1 = spVar[burnin:nsim,1],
                                            tau2_trace2 = spVar[burnin:nsim,2],
                                            gamma_trace = gammas[burnin:nsim])))
          }else{
            ess <- effectiveSize(mcmc(cbind(beta_trace = c(Betas)[burnin:nsim],
                                            cond_trace = Condition[burnin:nsim],
                                            sig2_trace = 1/taus[burnin:nsim],
                                            tau2_trace1 = spVar[burnin:nsim,1],
                                            tau2_trace2 = spVar[burnin:nsim,2],
                                            gamma_trace = gammas[burnin:nsim])))
          }
        }

      }) #time


      if(trace & !is.null(bioRep)){
        res[[feat]] <-list(
          beta = mbeta,
          cond = malpha,
          cond0 = malpha0,
          cond1 = malpha1,
          sig2 = msigma.e2,
          sig2bio = msigma.b2_bio,
          sig2tec = msigma.b2_tec,
          tau2 = msigma.t2,
          gamma = gam,
          ess = ess,
          trace = mcmc(cbind(beta_trace = c(Betas),
                             cond_trace = Condition,
                             sig2_trace = 1/taus,
                             sig2bio_trace = 1/taus_bio,
                             sig2tec_trace = 1/taus_tec,
                             tau2_trace1 = spVar[,1],
                             tau2_trace2 = spVar[,2],
                             gamma_trace = gammas)),
          time = time
        )
      }else if(trace & is.null(bioRep)){
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
          sig2bio = msigma.b2_bio,
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
} #single or multi tissue experiment
  }#if we don't drop zeros and instead just add a small value before log transform

}#function



