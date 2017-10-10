library(dlm)

lake <- read.table("lakeSuperior.dat", col.names = c("year","precipitation_inches"))
plot(lake, type ='o', col = "seagreen")

lake <- lake$precipitation_inches

buildFun <- function(x){ dlmModPoly(2, dV = exp(x[1]), dW=c(0,exp(x[2])), m0=c(10,10), C0=matrix(c(10,0,0,5), 2, 2)) }
fit <- dlmMLE(lake, parm = c(1,1), build = buildFun)
(dlmLake <- buildFun(fit$par))
StructTS(x = lake, type = "trend")
FiltD   <- dlm::dlmFilter(lake, dlmLake)
bsample <- dlm::dlmBSample(FiltD)

lakeFilt <- dlmFilter(lake, dlmLake)
bsample  <- dlmBSample(lakeFilt)


mcmc1 <- function(y, niter, sval){
  
  # Summary Stats
  n <- nrow(y)
  
  # setup storage
  thetakeep <- matrix(data=NA, ncol=2*(n+1), nrow=niter)
  phikeep   <- matrix(data=NA, ncol=2, nrow=niter)
  
  for(i in 1:niter){
    # sample theta jointly via backwards sampling
    dlmM    <- dlm::dlmModPoly(2, dV = sval$v, dW=sval$w)
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
  
  #colnames(thetakeep)[1:n+1] <- paste("mu[",seq(0:n),"]",sep="")
  #colnames(thetakeep)[-c(1:n+1)] <- paste("beta[",seq(0:n),"]",sep="")
  #colnames(thetakeep)[-c(1:n+1)] <- paste("beta[",seq(0:n),"]",sep="")
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
      for(t in 1:length(y)){
        dlmM      <- dlm::dlmModPoly(2, dV = phi[1], dW=c(0,phi[2]))
        FiltD     <- dlm::dlmFilter(y[t], dlmM)
        Rt        <- dlm::dlmSvd2var(FiltD$U.R, FiltD$D.R)
        tmpPhi[t] <- dnorm(y[t], FiltD$a[1], Rt[[1]][1,1])
      }
      return(prod(tmpPhi)*dinvgamma(phi[1], a, b)*dinvgamma(phi[2], a, b))
    }
    
    propphi <- MASS::mvrnorm(n = 1, prevphi, Sigma)
    
    if(sum(propphi < 0) == 1){
      phi <- prevphi
    }else{
      u <- runif(1, 0, 1)
      ifelse(u > dphi(prevphi, y, 1, 1)/dphi(propphi, y, 1, 1), phi <- prevphi, phi <- propphi)
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

# using mles found above as starting values
Sigma <- matrix(data=c(10,0,0,0.05), 2, 2, byrow=TRUE)
system.time(res2 <- mcmc2(lake, 1000, sval, Sigma))






ggplot() + geom_histogram(aes(x = rnorm(100, 0.5, 0.0005), y = ..density../1000),bin=100)
ggplot() + geom_histogram(aes(x = rnorm(100, 0.5, 0.0005), y = ..density..),bin=100)


