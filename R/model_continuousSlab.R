#' @title fit model
#' @description fit model
#' @param msset adfdf
#' @return res
#' @import mvtnorm
#' @import lme4
#' @import spam
#' @export
#'

spatialComparison <- function(msset,sample,conditionOfInterest,
                              feature, nsim=5000, burnin = 2500, trace = T,
                              piPrior = .1, seed = 1, logbase2 = F){

  set.seed(seed)

sample <- factor(sample)
conditionOfInterest <- factor(conditionOfInterest)
sampleNames <- levels(sample)
conditionNames <- levels(conditionOfInterest)
nCond <- length(conditionNames)

N <- length(sample)
nis <- sapply(sampleNames, function(x) sum(sample == x))
n <- length(sampleNames)
id <- as.numeric(sample)
conditionVec <- ifelse(conditionOfInterest == conditionNames[1], 0, 1)
X <- matrix(rep(1, N), ncol = 1) #design matrix for intercept and nuisance variables
X1 <- matrix(conditionVec, ncol = 1)  #design matrix without condition effect
numCond2 <- sum(conditionVec == 1)

k <-ncol(X)
res <- list()

# print(paste("Cond of Interest = "))
# print(conditionOfInterest )
# print("levels = ")
# print(levels(conditionOfInterest))

condAndSample <- data.frame(index = numeric(N), name =character(N))
numSpatialParams <- 0
nsl <- c()
j <- 1
for(s in levels(factor(id))){
  ind_samp <- id == s
  for(l in conditionNames){
    ind_cond <- conditionOfInterest == l


    if(any(ind_samp & ind_cond)){
      condAndSample[ind_samp & ind_cond,]$index <- j
      levels(condAndSample$name) <- c(levels(condAndSample$name), paste(s, l,sep="_"))
      condAndSample[ind_samp & ind_cond,]$name <- paste(s, l,sep="_")
      numSpatialParams <-  numSpatialParams +1

      #####################################################
      ##################### Initialize W ##################
      #####################################################
      assign(paste("W", s, l, sep="_"), adj.grid(coord(msset[,ind_cond & ind_samp]))+0)
      assign(paste("m", s, l, sep="_"), rowSums(get(paste("W", s, l, sep="_"))))
      nsl <- c(nsl, length(get(paste("m", s, l, sep="_"))))
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
phiVec_m <- rep(0, N)


feat <- 1

   for(f in feature){
     print(paste0("Feature ", f, " of ", length(feature)))
     time <- system.time({
   y <- spectra(msset)[f,]

   if(logbase2){
     y[y==0] <- .001
     y <- log2(y)
   }


###########
# Priors  #
###########
beta0<-rep(0,k)			# Prior Mean for beta
prec0<-diag(.01,k)		# Prior Precision Matrix of beta (vague), independent
precAlpha0 <- .01 #Prior Precision of slab
d0<-g0<-.001			# Hyperpriors for tau, taub
rd <- .00001 # ratio of varSpike/varSlab


#########
# Inits #
#########
lm <- lm(y~X+X1)
coef <- coef(lm)[-2]
tau<-1				# Error precision = 1/sigma2
eps_m.var <- 1/tau
b<-rep(0,n)			# Random effects (int and slope)
taub<-1				# Random Effects precision
#beta<-rep(0,k)			#inital value of intercept
beta <- coef[1:k]
alpha <- coef[k+1]
Z<-as.spam(matrix(0,N,n))	# Random effect design used for updating b
for (i in 1:n) Z[id==i,i]<-1
xb <-  X%*%beta
x1a <- X1 %*% alpha
zb <-  rep(0, N)
gamma <- 1
tauVar <- rep(1000,nSp)


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



###################
# GIBBS SAMPLER	#
###################
for (i in 1:nsim) {



  # Update Beta
  #intercept and nuisance effects
  vbeta<-solve(prec0+tau*crossprod(X,X))
  mbeta<-vbeta%*%(prec0%*%beta0 + tau*crossprod(X,y-x1a-zb-phiVec_m))
  beta <-c(rmvnorm(1,mbeta,vbeta))
  xb <-  X%*%beta
  Betas[i,]<- beta
#print(beta)

  resa <- sum((y-xb-zb-phiVec_m)[conditionVec == 1])
  # Update Condition effect
  if(gamma == 1){ #slab
    valph <- 1/(numCond2/eps_m.var + precAlpha0)
    malph <- valph*resa/eps_m.var
    Condition1[i] <- alpha <- rnorm(n = 1, mean = malph, sd = sqrt(valph))
  }else{
    valph <- 1/(numCond2/eps_m.var + 1/rd * precAlpha0)
    malph <- valph*resa/eps_m.var
    Condition0[i] <- alpha <- rnorm(n = 1, mean = malph, sd = sqrt(valph))
  }

  Condition[i] <- alpha
#print(alpha)
  x1a <- X1 %*% alpha
#print(length(x1a))


  # update indicator of differential abundance
  loglik_slab <- dnorm(alpha, mean = 0, sd = sqrt(1/precAlpha0), log = T)
  loglik_spike <- dnorm(alpha, mean = 0 , sd = sqrt(rd/precAlpha0), log = T)
  pi1Post <-  1/(1 + exp(loglik_spike - loglik_slab)*(1-piPrior)/piPrior )


   gamma <- rbinom(n=1, size = 1, prob = pi1Post)
  gammas[i] <-gamma

  # print("xb1")
  # print(head(xb1))
  # print("xb0")
  # print(head(xb0))
  #

  # print("gamma")
  # print(gamma)
  # print("gamma = 1")
  # print(beta_g1)
  # print("gamma = 0")
  # print(beta_g0)

  if(n > 1){
  # Update b
  vb<-1/(taub+nis*tau)
  mb<-vb*(tau*t(Z)%*%(y-x1a-xb-phiVec_m))  # Apparently Crossprod doesn't work with Spammed Z if n is large (e.g., > 1500 or so)??
  b<-rnorm(n,mb,sqrt(vb))
  }

  # Update tau
  if(n > 1){
  zb<-Z%*%b
  }else{
    zb <- rep(0, N)
  }
  g<-g0+crossprod(y-xb-x1a-zb-phiVec_m,y-xb-x1a-zb-phiVec_m)/2
  taus[i]<-tau<-rgamma(1,d,g)
  eps_m.var <- 1/tau

  if(n > 1){
  # Update taub
  m<-c(g0+crossprod(b,b)/2)
  taubs[i]<-taub<-rgamma(1,nu,m)
  }else{
    taubs[i]<-taub<- NA
  }

  ####### Fix tau2 for now
  #spVar[i,] <- tauVar_m <- rep(.01, n)

  offset.phi <- (y-xb-x1a-zb) / eps_m.var


  j <- 1
  for(s in levels(factor(id))){
    ind_samp <- id == s
    for(l in conditionNames){
      ind_cond <- conditionOfInterest == l




      if(any(ind_samp & ind_cond)){

        offset <- offset.phi[ind_samp & ind_cond]

# print(phiVec_m[ind_samp & ind_cond])
# print(paste(nsl[j],tauVar[j], eps_m.var, offset[1]))
# print(get(paste("m", s, l, sep="_")))
# print(get(paste("Wtrip", s, l, sep="_")))
# print(get(paste("Wbegfin", s, l, sep="_")))

        # print(class(phiVec_m[ind_samp & ind_cond]))
        # print(paste(class(as.numeric(nsl[j])),class(tauVar[j]), class(eps_m.var), class(offset[1])))
        # print(class(get(paste("m", s, l, sep="_"))))
        # print(class(get(paste("Wtrip", s, l, sep="_"))))
        # print(class(get(paste("Wbegfin", s, l, sep="_"))))

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

#print(dim(Betas))

  if (i%%1000==0 || i == 1) print(paste0("MCMC Iteration ", i, " of ", nsim))
} #mcmc iteration

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



