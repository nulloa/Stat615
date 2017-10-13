library("dlm")
library("dplyr")
library("reshape2")
library("ggplot2")
OGlake <- read.table("lakeSuperior.dat", col.names = c("year","precipitation_inches"))
lake <- OGlake$precipitation_inches
plot(lake, type ='o', col = "seagreen")


# using mle as starting values
buildFun <- function(x){ dlmModPoly(2, dV = exp(x[1]), dW=c(0,exp(x[2])), m0=c(30,9), C0=matrix(c(100^2,0,0,1^2), 2, 2)) }
fit <- dlmMLE(lake, parm = c(1,1), build = buildFun)
sval <- list(v=exp(fit$par[1]), w=exp(fit$par[2]))


mcmc1 <- function(y, niter, sval){
  # Summary Stats
  n <- length(y)
  
  # setup storage
  thetakeep <- matrix(data=NA, ncol=2*(n+1), nrow=niter)
  phikeep   <- matrix(data=NA, ncol=2, nrow=niter)
  
  # setup initialization
  sigma  <- sval$v
  sigmab <- sval$w
  
  for(i in 1:niter){
    if(i%%(niter/10) == 0){print(paste(100*i/niter,"% of iterations complete", sep=""))}
    # sample theta jointly via backwards sampling
    dlmM    <- dlm::dlmModPoly(2, dV = sigma, dW=c(0,sval$w))
    FiltD   <- dlm::dlmFilter(y, dlmM)
    bsample <- dlm::dlmBSample(FiltD)
    
    # sample phi 
    SSEs   <- sum((y - bsample[-1,1])^2)
    SSEb   <- sum(diff(bsample[,2])^2)
    
    sigma  <- MCMCpack::rinvgamma(1, (n/2)+1, (SSEs/2) + 1)
    sigmab <- MCMCpack::rinvgamma(1, (n/2)+1, (SSEb/2) + 1)
    
    # save samples
    thetakeep[i,]  <- c(bsample[,1],bsample[,2])
    phikeep[i,1]   <- sigma
    phikeep[i,2]   <- sigmab
  }
  
  # end of sampling mechanism
  return(list("theta"=thetakeep, "phi"=phikeep))
}



mcmc2 <- function(y, niter, sval, Sigma){
  # Summary Stats
  n <- length(y)
  
  # setup storage
  thetakeep <- matrix(data=NA, ncol=2*(n+1), nrow=niter)
  phikeep   <- matrix(data=NA, ncol=2, nrow=niter)
  
  # setup initialization
  prevphi  <- c(sval$v, sval$w)
  
  dinvgamma = function(x, a, b) dgamma(1/x,a,b)/x^2
  
  for(i in 1:niter){
    if(i%%(niter/10) == 0){print(paste(100*i/niter,"% of iterations complete", sep=""))}
    
    # sample phi 
    dphi <- function(phi, y, a, b){
      dinvgamma = function(x, a, b) dgamma(1/x,a,b)/x^2
      tmpPhi <- NULL
      dlmM      <- dlm::dlmModPoly(2, dV = phi[1], dW=c(0,phi[2]), m0=c(10,10), C0=matrix(c(10,0,0,5), 2, 2))
      FiltD     <- dlm::dlmFilter(y, dlmM)
      Rt        <- dlm::dlmSvd2var(FiltD$U.R, FiltD$D.R)
      for(t in 1:length(y)){
        prop_var <- dlmM$FF %*% Rt[[t]] %*% t(dlmM$FF) + dlmM$V
        tmpPhi[t] <- dnorm(y[t], FiltD$f[t], sqrt(prop_var), log=TRUE)
      }
      return(sum(tmpPhi) + log(dinvgamma(phi[1], a, b)) + log(dinvgamma(phi[2], a, b)))
    }
    
    propphi <- MASS::mvrnorm(n = 1, prevphi, Sigma)
    
    if(sum(propphi < 0) == 1){
      phi <- prevphi
    }else{
      u <- runif(1, 0, 1)
      # print(dphi(prevphi, y, 1, 1)/dphi(propphi, y, 1, 1))
      ifelse(u > exp(dphi(propphi, y, 1, 1)) / exp(dphi(prevphi, y, 1, 1)), phi <- prevphi, phi <- propphi)
    }
    
    # sample theta
    dlmM    <- dlm::dlmModPoly(2, dV = phi[1], dW=c(0,phi[2]))
    FiltD   <- dlm::dlmFilter(y, dlmM)
    bsample <- dlm::dlmBSample(FiltD)
    
    # save samples
    thetakeep[i,] <- c(bsample[,1],bsample[,2])
    phikeep[i,1]  <- phi[1]
    phikeep[i,2]  <- phi[2]
    prevphi       <- phi
  }
  
  # end of sampling mechanism
  return(list("theta"=thetakeep, "phi"=phikeep))
}

Sigma0 <- matrix(data=c(10,0,0,0.05), 2, 2, byrow=TRUE)
sigres <- mcmc2(lake, 1000, sval, Sigma0)
covmat <- cov(sigres$phi)

