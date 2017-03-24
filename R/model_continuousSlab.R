#' @title fit model
#' @description fit model
#' @param msset adfdf
#' @return res
#' @import mvtnorm
#' @import lme4
#' @import spam
#' @import Cardinal
#' @export
#'

spatialComparison <- function(msset,sample,conditionOfInterest,
                              feature, nsim=5000, burnin = 2500, trace = T,
                              piPrior = .1, seed = 1, logbase2 = F, coord = NULL,
                              type.neighbor = "radius", radius.neighbor = 1, maxdist.neighbor = NULL){
  
  if(is.null(coord)){
    coord <- coord(msset)
  }

  print("Initializing spatial components...")
  set.seed(seed) #random seed

  sample <- factor(sample) # make the sample labels a factor in case it is a character vector
  conditionOfInterest <- factor(conditionOfInterest) # make the condition labels a factor in case it is a character vector
  sampleNames <- levels(sample) #create vector of sample names
  conditionNames <- levels(conditionOfInterest) #create vector of condition names
  nCond <- length(conditionNames) #how many conditions there are

  N <- length(sample) # how many pixels there are
  nis <- sapply(sampleNames, function(x) sum(sample == x)) #how many pixels in each sample
  n <- length(sampleNames) #now many samples there are
  id <- as.numeric(sample) #convert sample names into id numbers

  #### ONLY IF THERE ARE TWO CONDITIONS #####
  conditionVec <- ifelse(conditionOfInterest == conditionNames[1], 0, 1) #vector converts condition names from characters to numeric
  numCond2 <- sum(conditionVec == 1) #number of pixels from condition two
  ###########################################

  X <- matrix(rep(1, N), ncol = 1) #design matrix for intercept and covariates #currently set to intercept only
  X1 <- matrix(conditionVec, ncol = 1)  #design matrix without condition effect


  k <-ncol(X) #number of covariates, including intercept
  res <- list() #list that will hold results


  ####################################################################################################
  ############## Obtain neighborhood matrices for each combination of sample and condition ###########
  ####################################################################################################
  condAndSample <- data.frame(index = numeric(N), name =character(N)) #dataframe that numbers the sample/condition combinations
  numSpatialParams <- 0 #number of spatial parameters to estimate. this will be the number unique of sample/condition pairs that exist in the dataset
  nsl <- c() ##### vector of number of pixels from each condition and sample pair
  j <- 1
  for(s in levels(factor(id))){
    ind_samp <- id == s # pixels from sample s
    for(l in conditionNames){
      ind_cond <- conditionOfInterest == l #pixels from condition l


      if(any(ind_samp & ind_cond)){ #check to see if sample/condition combination exists in dataset
        condAndSample[ind_samp & ind_cond,]$index <- j  #give the same index to all pixels from a given condition and sample
        levels(condAndSample$name) <- c(levels(condAndSample$name), paste(s, l,sep="_")) #add this combination of condition and sample to the list of names
        condAndSample[ind_samp & ind_cond,]$name <- paste(s, l,sep="_") #give the same name to all pixels from the same condition and sample
        numSpatialParams <-  numSpatialParams +1

        #####################################################
        ##################### Initialize W ##################
        #####################################################
        #### Create adjacency matrix for pixels from this condition and sample pair
        assign(paste("W", s, l, sep="_"), adj.grid(coords = coord[ind_cond & ind_samp,],
                                                   type = type.neighbor, 
                                                   radius = radius.neighbor, 
                                                   max.dist = maxdist.neighbor)+0)
        #### number of neighbors for each pixel
        assign(paste("m", s, l, sep="_"), rowSums(get(paste("W", s, l, sep="_"))))
        ##### number of pixels from this condition and sample pair
        nsl <- c(nsl, sum(ind_cond & ind_samp))
        names(nsl)[j] <- paste(s,l)

      ###################  ###################  ###################
      ####### Stuff for CARBayes

      assign(paste("Wtrip", s, l, sep="_"), triplet(as.spam(get(paste("W", s, l, sep="_")))))
      assign(paste("Wtrip", s, l, sep="_"), cbind(get(paste("Wtrip", s, l, sep="_"))$indices, get(paste("Wtrip", s, l, sep="_"))$values))
      Temp <-  array(NA, c(nsl[j], 2))
      temp <- 1



      for(i in 1:nsl[j])
      {
        Temp[i, ] <- c(temp, (temp + get(paste("m", s, l, sep="_"))[i]-1))
        temp <- temp + get(paste("m", s, l, sep="_"))[i]
      }

      assign(paste("Wbegfin", s, l, sep="_"), Temp)
      rm(Temp)



      #####################################################
      #####################################################
      j <- j+1
    }
  }
}

nSp <- length(nsl)

condAndSample$name<-factor(condAndSample$name)
####################################################################################################
####################################################################################################
####################################################################################################


phiVec_m <- rep(0, N) #initialize trace vector for spatial effects
feat <- 1 #initialize feature index

print("...Initialization done.")
####################################################################################################
##################################### Fit model feature by feature #################################
####################################################################################################

for(f in feature){
  print(paste0("Feature ", f, " of ", length(feature)))
  time <- system.time({ #time the overall model fits
    y <- spectra(msset)[f,]

    if(logbase2){ #do log transformation if necessary
      y[y==0] <- .001 #zeros in the image will cause problems if a log transformation is required. add a small number to the zeroes.
      y <- log2(y)
    }



    ####################################################################################################
    ################################### Set up prior distributions #####################################
    ####################################################################################################
    beta0<-rep(0,k)			# Prior Mean for beta
    prec0<-diag(.01,k)		# Prior Precision Matrix of beta (vague), independent
    precAlpha0 <- .01 #Prior Precision of slab (value of condition effect if it is not zero)
    d0<-g0<-.001			# Hyperpriors for tau, taub
    rd <- .00001 # ratio of varSpike/varSlab


    ####################################################################################################
    ################################### Initialize variables  #####################################
    ####################################################################################################

    lm <- lm(y~X+X1) ## fit linear model to get reasonable starting values
    coef <- coef(lm)[-2]
    tau<-1				# technical error precision
    eps_m.var <- 1/tau #technical error variance
    b<-rep(0,n)			# Random effects (int and slope)
    taub<-1				# Random Effects precision
    beta <- coef[1:k] #initial value of intercept and covariates
    alpha <- coef[k+1] #initial value of condition effect
    Z<-as.spam(matrix(0,N,n))	# Random effect design used for updating b
    for (i in 1:n) Z[id==i,i]<-1
    xb <-  X%*%beta
    x1a <- X1 %*% alpha
    zb <-  rep(0, N)
    gamma <- 1 # initiate condition effect as nonzero
    tauVar <- rep(1,numSpatialParams) # spatial variances


#################
# Store Results #
#################

Betas<-matrix(0,nsim,k)	# Fixed Effects
spVar<-matrix(0,nsim,nSp)
taus<-taubs<-gammas <- rep(0,nsim)
Condition <- Condition0 <- Condition1 <- rep(NA,nsim)	# Error Precision Parms

###############################
# Fixed Posterior Hyperparms 	#
#    for tau and taub		#
###############################
d<-d0+N/2
nu<-d0+n/2



####################################################################################################
######################################## THE GIBBS SAMPLER  ########################################
####################################################################################################
for (i in 1:nsim) { #this is an iterative method, nsim is the number of iterations

  # Update intercept and covariates
  vbeta<-solve(prec0+tau*crossprod(X,X))
  mbeta<-vbeta%*%(prec0%*%beta0 + tau*crossprod(X,y-x1a-zb-phiVec_m))
  beta <-c(rmvnorm(1,mbeta,vbeta))
  xb <-  X%*%beta
  Betas[i,]<- beta


  resa <- sum((y-xb-zb-phiVec_m)[conditionVec == 1]) #residuals for pixels in second condition onlt

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


  if(n > 1){
    # Update the sample-to-sample effect
    vb<-1/(taub+nis*tau)
    mb<-vb*(tau*t(Z)%*%(y-x1a-xb-phiVec_m))
    b<-rnorm(n,mb,sqrt(vb))
  }

  # Update the technical error precision
  if(n > 1){
    zb<-Z%*%b
  }else{
    zb <- rep(0, N)
  }
  g<-g0+crossprod(y-xb-x1a-zb-phiVec_m,y-xb-x1a-zb-phiVec_m)/2
  taus[i]<-tau<-rgamma(1,d,g)
  eps_m.var <- 1/tau

  if(n > 1){
    # Update the precision of the sample effect
    m<-c(g0+crossprod(b,b)/2)
    taubs[i]<-taub<-rgamma(1,nu,m)
  }else{
    taubs[i]<-taub<- NA
  }


  offset.phi <- (y-xb-x1a-zb) / eps_m.var

  #########################################################
  ########### Update the spatial effects ##################
  #########################################################

  j <- 1
  for(s in levels(factor(id))){
    ind_samp <- id == s
    for(l in conditionNames){
      ind_cond <- conditionOfInterest == l

      if(any(ind_samp & ind_cond)){

        offset <- offset.phi[ind_samp & ind_cond]

         phiUpdate <- updateSpatial(
                      Wtrip=get(paste("Wtrip", s, l, sep="_")),
                      Wbegfin=get(paste("Wbegfin", s, l, sep="_")),
                      m = get(paste("m", s, l, sep="_")),
                      nsl= nsl[j],
                      phiVec=phiVec_m[ind_samp & ind_cond], #
                      tau2=tauVar[j], #
                      rho=1,
                      eps_m.var =eps_m.var,
                      offset.phi =offset,
                      tauVar.a = .001,
                      tauVar.b = .001
        )



         phiVec_m[ind_samp & ind_cond] <- phiUpdate$phi
         spVar[i,j] <- tauVar[j] <- phiUpdate$tau2
  j <- j+1
      }
    }
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
if(n > 1){
msigma.b2<-mean(1/taubs[(burnin):nsim])
}else{
  msigma.b2<- NA
}
msigma.t2<-apply(spVar[(burnin):nsim,, drop = F],2,mean)
gam <- mean(gammas)
malpha <- mean(Condition[burnin:nsim])
malpha1 <- mean(Condition1[burnin:nsim], na.rm = T)
malpha0 <- mean(Condition0[burnin:nsim], na.rm = T)

     }) #time


if(trace){
  res[[feat]] <-list(
    beta = mbeta,
    cond = malpha,
    cond0 = malpha0,
    cond1 = malpha1,
    sig2 = msigma.e2,
    sig2b = msigma.b2,
    tau2 = msigma.t2,
    gamma = gam,
    beta_trace = Betas,
    cond_trace = Condition,
    cond1_trace = Condition1,
    cond0_trace = Condition0,
    sig2_trace = 1/taus,
    sig2b_trace = 1/taubs,
    tau2_trace = spVar,
    gamma_trace = gammas,
    time = time
  )
}else{
  res[[feat]] <-list(
    beta = mbeta,
    cond = malpha,
    cond0 = malpha0,
    cond1 = malpha1,
    sig2 = msigma.e2,
    sig2b = msigma.b2,
    tau2 = msigma.t2,
    gamma = gam,
    time = time
  )
}

names(res)[feat] <- paste0("Feature",f)
feat <- feat + 1


  } #feature


  return(res)
}#function



