#' @title sample car
#' @description sample car
#' @param tau2 safaf
#' @return sample
#' @import spam
#' @import mvtnorm
#' @export
#'
#Function to generate an nxn square of CAR noise given a spatial variance parameter tau2
# and a propriety parameter rho (optional. defaults to near ICAR), and a number of replicates
#nrep

#function generates samples and saves them to rdata file

#coord: coordinates of simulated image
#condDiff: true underlying condition difference
#pattern: factor, pattern that the two differences occur in
#sig2: measurment error variance


#source("/Users/April/Documents/Research/R Code/Model with indicator term/ICARsamplingFunctions.R") #for adj.grid
#setwd("/Users/April/Documents/Research/R Code/Model with indicator term/10x10")


#Wstd is row standardized adj. matrix
#Dw is diagonal matrix of row sums
#rho is propriety parameter, |rho| < 1
#SigmaInv = tau^2 * Dw %*% (I - diag(rho, Wstd )

sampleCAR <- function(condDiff, coord = expand.grid(x=1:10, y=1:10),
                      pattern = ifelse((coords$x %in% 3:8 & coords$y %in% 3:8), 2, 1),
                      sig2, tau2, rho = .9999, nrep=100, save = F, randomSeed = 1,
                      neighbor.type = "radius", radius = 1){

            set.seed(randomSeed)

            W <- adj.grid(coord, sample= factor(rep(1, nrow(coord))),
                          type = neighbor.type, radius = radius)
            n <- nrow(W)
            m <- rowSums(W)
            Wstd <- diag(1/m)%*%W
            Dw <- diag(m)

            phisamp <- matrix(nrow = nrep, ncol = n)
            tau_1m <- sqrt(tau2)					      # Spatial SD
            SigmaInv <-as.spam(diag(m,n) - rho * W)
            Sigma <-  tau2 * solve(SigmaInv)			# Covariance of phis

            for(i in 1:nrep){
              # Spatial Random Effects +  # Measurement error + # Condition Effect
              sp <- c(rmvnorm(1,sigma=Sigma))
              phisamp[i,] <- (sp - mean(sp)) + rnorm(n, mean = 0, sd = sqrt(sig2))  + ifelse(pattern == 2, condDiff, 0)
            }



            if(save){
              assign(paste("phi_r", rho, "_t", tau2,"_s",sig2,"_a",condDiff, sep=""), phisamp)
              save(list = paste("phi_r", rho, "_t", tau2,"_s",sig2,"_a",condDiff, sep=""),
                 file = paste("phi_r", rho, "_t", tau2,"_s",sig2,"_a",condDiff, ".rdata", sep="") )
            }

            return(phisamp)
}




