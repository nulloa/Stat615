library(dlm)

lake <- read.table("lakeSuperior.dat", col.names = c("year","precipitation_inches"))
plot(lake, type ='o', col = "seagreen")

lake <- as.ts(lake$precipitation_inches)

buildFun <- function(x){ dlmModPoly(2, dV = exp(x[1]), dW=c(0,exp(x[2]))) }
fit <- dlmMLE(lake, parm = c(0,0), build = buildFun)
dlmLake <- buildFun(fit$par)

lakeFilt <- dlmFilter(lake, dlmLake)
bsample  <- dlmBSample(lakeFilt)


mcmc1 <- function(y, niter, sval){
  # setup storage
  thetakeep <- matrix(data=NA, ncol=nrow(y)+1, nrow=niter)
  phikeep   <- matrix(data=NA, ncol=nrow(y)+1, nrow=niter)
  
  for(i in 1:niter){
    # sample theta jointly via backwards sampling
    dlmM    <- dlm::dlmModPoly(2, dV = sval$v, dW=sval$w)
    FiltD   <- dlm::dlmFilter(y, dlmM)
    bsample <- dlm::dlmBSample(FiltD)
    
    # sample phi 
    dphi <- function(phi, y, a, b){
      for(t in 1:nrow(t)){
        dlmM      <- dlm::dlmModPoly(2, dV = sval$v, dW=sval$w)
        FiltD     <- dlm::dlmFilter(y[t], dlmM)
        Rt        <- dlm::dlmSvd2var(FiltD$U.C, FiltD$D.C)
        tmpPhi[t] <- dnorm(y[t], filtD$a, Rt)
      }
      return(prod(tmpPhi)*dinvgamma(phi, a, b))
    }
    phi <- 
    
    # save samples
    thetakeep[i,] <- bsample
    phikeep[i,]   <- phi
  }
  
  # end of sampling mechanism
  return(list("theta"=thetakeep, "phi"=phikeep))
}



mcmc2 <- function(y, niter, sval){
  # setup storage
  thetakeep <- matrix(data=NA, ncol=nrow(y)+1, nrow=niter)
  phikeep   <- matrix(data=NA, ncol=nrow(y)+1, nrow=niter)
  tmpPhi    <- NULL
  
  dinvgamma = function(x, a, b) dgamma(1/x,a,b)/x^2
  
  for(i in 1:niter){
    # sample phi 
    dphi <- function(phi, y, a, b){
      for(t in 1:nrow(t)){
        dlmM      <- dlm::dlmModPoly(2, dV = sval$v, dW=sval$w)
        FiltD     <- dlm::dlmFilter(y[t], dlmM)
        Rt        <- dlm::dlmSvd2var(FiltD$U.C, FiltD$D.C)
        tmpPhi[t] <- dnorm(y[t], filtD$a, Rt)
      }
      return(prod(tmpPhi)*dinvgamma(phi, a, b))
    }
    
    
    propphi <- MASS::mvrnorm(n = 1, prevphi, Sigma)
    
    u <- runif(1, 0, 1)
    ifelse(u > dphi(prevphi, y, 1, 1)/dphi(propphiy, 1, 1), phi <- prevphi, phi <- propphi)
    
    # sample theta
    dlmM    <-  dlm::dlmModPoly(2, dV = phi[1], dW=c(0,phi[2]))
    FiltD   <- dlm::dlmFilter(y, dlmM)
    bsample <- dlm::dlmBSample(FiltD)
    
    
    # save samples
    thetakeep[i,] <- bsample
    phikeep[i,]   <- phi
    prevphi <- phi
  }
  
  # end of sampling mechanism
  return(list("theta"=thetakeep, "phi"=phikeep))
}










