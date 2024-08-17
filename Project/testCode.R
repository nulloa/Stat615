##### Load Libraries #####
library("ggplot2")
library("cdcfluview")
library("plyr")
library("dplyr")
library("reshape2")
library("boot")
library("dlm")


##### Save Data #####
# fludat <- ilinet("hhs", years=2010:2013)
# save(fludat, file="flu1013.RData")


##### Setup Data #####
fludat$region <- as.factor(fludat$region)
fludat <- dplyr::filter(fludat, region=="Region 9")
names(fludat) <- tolower(names(fludat))
fludat$time <- 1:nrow(fludat)
ggplot(fludat) + geom_point(aes(x=week, y=ilitotal/`total_patients`)) + facet_wrap(~year)
ggplot(fludat) + geom_point(aes(x=time, y=ilitotal/`total_patients`))
flu <- as.numeric(unlist(fludat[,c("ilitotal")]))/as.numeric(unlist(fludat[,c("total_patients")]))
cflu <- (flu - mean(flu))
fluts <- ts(fludat$ilitotal/fludat$`total_patients`, start=c(2010, 40), frequency = 52)

##### Plot Time-Series #####
plot(decompose(fluts)) 




##### Run Seasonal Analysis #####

# Find Starting Values
buildFun <- function(x){ dlmModTrig(s=52, q = 1, dV = exp(x[1])) +
    dlmModPoly(2, dV = 0, dW = c(0, exp(x[2])), m0=rep(0,2), C0=100*diag(2)) }
fit <- dlmMLE(logit(fluts), parm = c(1,1), build = buildFun)
sval1 <- list(s=exp(fit$par[1]), b=exp(fit$par[2]))

buildFun <- function(x){ dlmModTrig(s=52, q = 1, dV = exp(x[1]), dW=exp(x[2])) +
    dlmModPoly(2, dV = 0, dW = c(0, exp(x[3])), m0=rep(0,2), C0=100*diag(2)) }
fit <- dlmMLE(logit(fluts), parm = c(1,1,1), build = buildFun)
sval1a <- list(s=exp(fit$par[1]), w=exp(fit$par[2]), b=exp(fit$par[3]))

buildFun <- function(x){ dlmModTrig(s=52, q = 2, dV = exp(x[1])) +
    dlmModPoly(2, dV = 0, dW = c(0, exp(x[2])), m0=rep(0,2), C0=100*diag(2)) }
fit <- dlmMLE(logit(fluts), parm = c(1,1), build = buildFun)
sval2 <- list(s=exp(fit$par[1]), b=exp(fit$par[2]))

buildFun <- function(x){ dlmModTrig(s=52, q = 2, dV = exp(x[1]), dW = exp(x[2])) +
    dlmModPoly(2, dV = 0, dW = c(0, exp(x[3])), m0=rep(0,2), C0=100*diag(2)) }
fit <- dlmMLE(logit(fluts), parm = c(1,1,1), build = buildFun)
sval2a <- list(s=exp(fit$par[1]), w=exp(fit$par[2]), b=exp(fit$par[3]))

# Check MLE Fit
flu.model1 = dlmModTrig(s=52, q = 1, dV = sval1$s) + 
  dlmModPoly(2, dV = 0, dW = c(0, sval1$b), m0=rep(0,2), C0=100*diag(2))
flu.model1a = dlmModTrig(s=52, q = 1, dV = sval1a$s, dW = sval1a$w) + 
  dlmModPoly(2, dV = 0, dW = c(0, sval1a$b), m0=rep(0,2), C0=100*diag(2))
flu.model2 = dlmModTrig(s=52, q = 2, dV = sval2$s) + 
  dlmModPoly(2, dV = 0, dW = c(0, sval2$b), m0=rep(0,2), C0=100*diag(2))
flu.model2a = dlmModTrig(s=52, q = 2, dV = sval2a$s, dW = sval2a$w) + 
  dlmModPoly(2, dV = 0, dW = c(0, sval2a$b), m0=rep(0,2), C0=100*diag(2))

flu.filter1 = dlmFilter(logit(fluts), flu.model1)
flu.filter1a = dlmFilter(logit(fluts), flu.model1a)
flu.filter2 = dlmFilter(logit(fluts), flu.model2)
flu.filter2a = dlmFilter(logit(fluts), flu.model2a)

bsample1 <- dlm::dlmBSample(flu.filter1)
bsample1a <- dlm::dlmBSample(flu.filter1a)
bsample2 <- dlm::dlmBSample(flu.filter2)
bsample2a <- dlm::dlmBSample(flu.filter2a)

tmpdat <- rbind(data.frame(x=seq(1, length(flu)), y=(flu), est = bsample1[-1,1]+bsample1[-1,3], model="1harm"),
                data.frame(x=seq(1, length(flu)), y=(flu), est = bsample1a[-1,1]+bsample1a[-1,3], model="1aharm"),
                data.frame(x=seq(1, length(flu)), y=(flu), est = bsample2[-1,1]+bsample2[-1,3]+bsample2[-1,5], model="2harm"))
                # data.frame(x=seq(1, length(flu)), y=(flu), est = bsample2a[-1,1]+bsample2a[-1,3]+bsample2a[-1,5], model="2aharm"))
ggplot(tmpdat) + 
  geom_line(aes(x=x, y=inv.logit(est), color=model)) +
  geom_point(aes(x=x, y=y)) + 
  labs(y="Percentage of ILI Patients", x="Time") + theme_bw()


ggplot()+geom_line(aes(x=c(1:length(flu)),y=bsample1[-1,1])) +
  geom_line(aes(x=c(1:length(flu)),y=bsample1a[-1,1]), color="red")


ggplot()+geom_line(aes(x=c(1:length(flu)),y=bsample1[-1,3])) +
  geom_line(aes(x=c(1:length(flu)),y=bsample1a[-1,3]), color="red")


# Using 1 Harmonic
mcmc1 <- function(y, niter, sval){
  # Summary Stats
  n <- length(y)
  
  # setup storage
  harmonic1a <- matrix(data=NA, ncol=(n+1), nrow=niter)
  harmonic1b <- matrix(data=NA, ncol=(n+1), nrow=niter)
  mu         <- matrix(data=NA, ncol=(n+1), nrow=niter)
  beta       <- matrix(data=NA, ncol=(n+1), nrow=niter)
  phikeep    <- matrix(data=NA, ncol=length(sval), nrow=niter)
  
  # setup initialization
  sigma  <- sval$s
  sigmab  <- sval$b
  
  for(i in 1:niter){
    if(i%%(niter/10) == 0){print(paste(100*i/niter,"% of iterations complete", sep=""))}
    # sample theta jointly via backwards sampling
    dlmM    <- dlm::dlmModTrig(om=32, q = 1, dV = sigma) + 
      dlm::dlmModPoly(2, dV = 0, dW = c(0, sigmab), m0=rep(0,2), C0=100*diag(2))
    FiltD   <- dlm::dlmFilter(y, dlmM)
    bsample <- dlm::dlmBSample(FiltD)
    
    # sample phi
    SSE  <- sum((y - (bsample[-1,1] + bsample[-1,3]))^2)
    SSEb  <- sum(diff(bsample[,4])^2)
    
    sigma <- MCMCpack::rinvgamma(1, (n/2)+2, (SSE/2) + 0.007)
    sigmab <- MCMCpack::rinvgamma(1, (n/2)+2, (SSEb/2) + 0.002)
    
    # save samples
    harmonic1a[i,]  <- c(bsample[,1])
    harmonic1b[i,]  <- c(bsample[,2])
    mu[i,]          <- c(bsample[,3])
    beta[i,]        <- c(bsample[,4])
    phikeep[i,]     <- c(sigma, sigmab)
  }
  
  # end of sampling mechanism
  return(list("harmonic1a"=harmonic1a,"harmonic1b"=harmonic1b,"mu"=mu,"beta"=beta, "phi"=phikeep))
}
time1 <- system.time(res1 <- mcmc1(logit(fluts), 10000, sval1))

# Using 2 Harmonics
mcmc2 <- function(y, niter, sval){
  # Summary Stats
  n <- length(y)
  
  # setup storage
  harmonic1a <- matrix(data=NA, ncol=(n+1), nrow=niter)
  harmonic1b <- matrix(data=NA, ncol=(n+1), nrow=niter)
  harmonic2a <- matrix(data=NA, ncol=(n+1), nrow=niter)
  harmonic2b <- matrix(data=NA, ncol=(n+1), nrow=niter)
  mu         <- matrix(data=NA, ncol=(n+1), nrow=niter)
  beta       <- matrix(data=NA, ncol=(n+1), nrow=niter)
  phikeep    <- matrix(data=NA, ncol=length(sval), nrow=niter)
  
  
  # setup initialization
  #sigmam  <- sval$vm
  sigma  <- sval$s
  sigmab  <- sval$b
  
  for(i in 1:niter){
    if(i%%(niter/10) == 0){print(paste(100*i/niter,"% of iterations complete", sep=""))}
    # sample theta jointly via backwards sampling
    dlmM    <- dlm::dlmModTrig(om=32, q = 2, dV = sigma) + 
      dlm::dlmModPoly(2, dV = 0, dW = c(0, sigmab), m0=rep(0,2), C0=100*diag(2))
    FiltD   <- dlm::dlmFilter(y, dlmM)
    bsample <- dlm::dlmBSample(FiltD)
    
    # sample phi
    SSE  <- sum((y - (bsample[-1,1] + bsample[-1,3] + bsample[-1,5]))^2)
    SSEb  <- sum(diff(bsample[,6])^2)
    
    sigma <- MCMCpack::rinvgamma(1, (n/2)+2, (SSE/2) + 0.007)
    sigmab <- MCMCpack::rinvgamma(1, (n/2)+2, (SSEb/2) + 0.002)
    
    # save samples
    harmonic1a[i,]  <- c(bsample[,1])
    harmonic1b[i,]  <- c(bsample[,2])
    harmonic2a[i,]  <- c(bsample[,3])
    harmonic2b[i,]  <- c(bsample[,4])
    mu[i,]          <- c(bsample[,5])
    beta[i,]        <- c(bsample[,6])
    phikeep[i,]     <- c(sigma, sigmab)
  }
  
  # end of sampling mechanism
  return(list("harmonic1a"=harmonic1a,"harmonic1b"=harmonic1b,
              "harmonic2a"=harmonic2a,"harmonic2b"=harmonic2b,
              "mu"=mu,"beta"=beta, "phi"=phikeep))
}
time2 <- system.time(res2 <- mcmc2(logit(fluts), 10000, sval2))





##### Diagnostics #####
# Geweke
gd <- data.frame(Harmonic1a=as.numeric(coda::geweke.diag(coda::mcmc(res1$harmonic1a))$z),
                 mu=as.numeric(coda::geweke.diag(coda::mcmc(res1$mu))$z),
                 beta=as.numeric(coda::geweke.diag(coda::mcmc(res1$beta))$z))
mgd <- melt(gd)
ggplot(mgd) + geom_histogram(aes(x=value)) + facet_wrap(~variable) + theme_bw()

# Variance plots
mphi2 <- melt(matrix(as.numeric(res2$phi[,]), ncol=2, byrow=F))
mphi2$Var2 <- factor(mphi2$Var2)
levels(mphi2$Var2) <- c("sigma", "sigma_beta")

p1 <- ggplot(mphi2) + geom_line(aes(x=Var1, y=sqrt(value))) + facet_wrap(~Var2, scales = "free") + 
  labs(x="Iterations", y="") + theme_bw()


# prior to posterior plots
d1 <- melt(data.frame(sigma=res2$phi[,1], sigma_beta=res2$phi[,2]))
dinvgamma = function(x, a, b) dgamma(1/x,a,b)/x^2
dsqrtinvgamma = function(x, a, b) dinvgamma(x^2, a, b)*2*x

p2 <- ggplot(data.frame(sigma=res2$phi[,1])) + 
  geom_histogram(aes(x = sqrt(sigma), y = ..density..), bins = 100) + 
  stat_function(fun = dsqrtinvgamma, color = "red", args = list(a = 2, b = 0.007)) +
  theme_bw() + labs(x="sigma", y="")

p3 <- ggplot(data.frame(sigma_beta=res2$phi[,2])) + 
  geom_histogram(aes(x = sqrt(sigma_beta), y = ..density..), bins = 100) + 
  stat_function(fun = dsqrtinvgamma, color = "red", args = list(a = 2, b = 0.002)) +
  theme_bw() + labs(x="sigma_beta", y="")

lay <- rbind(c(1,1),
             c(2,3))
grid.arrange(p1, p2, p3, layout_matrix = lay)


###### Results #####
n <- length(flu)
L1 <- U1 <- esttheta <- L2 <- U2 <- esttheta2 <- LT <- UT <- estT <- LT2 <- UT2 <- estT2 <- NULL
for(i in 1:n){
  L1[i]       <- quantile(inv.logit(res1$harmonic1a[,i+1] + res1$mu[,i+1]), probs=0.025)
  U1[i]       <- quantile(inv.logit(res1$harmonic1a[,i+1] + res1$mu[,i+1]), probs=0.975)
  esttheta[i] <- quantile(inv.logit(res1$harmonic1a[,i+1] + res1$mu[,i+1]), probs=0.5)
  LT[i]   <- quantile(inv.logit(res1$mu[,i+1]), probs=0.025)
  UT[i]   <- quantile(inv.logit(res1$mu[,i+1]), probs=0.975)
  estT[i] <- quantile(inv.logit(res1$mu[,i+1]), probs=0.5)
  L2[i]        <- quantile(inv.logit(res2$harmonic1a[,i+1] + res2$harmonic2a[,i+1] + res2$mu[,i+1]), probs=0.025)
  U2[i]        <- quantile(inv.logit(res2$harmonic1a[,i+1] + res2$harmonic2a[,i+1] + res2$mu[,i+1]), probs=0.975)
  esttheta2[i] <- quantile(inv.logit(res2$harmonic1a[,i+1] + res2$harmonic2a[,i+1] + res2$mu[,i+1]), probs=0.5)
  LT2[i]   <- quantile(inv.logit(res2$mu[,i+1]), probs=0.025)
  UT2[i]   <- quantile(inv.logit(res2$mu[,i+1]), probs=0.975)
  estT2[i] <- quantile(inv.logit(res2$mu[,i+1]), probs=0.5)
}

dat1 <- data.frame(L=L1, U=U1, x=seq(1,n), ogy=(flu), est=esttheta, LT=LT, UT=UT, estT=estT, num.harmonic="1")
dat2 <- data.frame(L=L2, U=U2, x=seq(1,n), ogy=(flu), est=esttheta2, LT=LT2, UT=UT2, estT=estT2, num.harmonic="2")
ggplot(rbind(dat1, dat2)) + 
  geom_segment(aes(x=x, xend=x, y=L, yend = U, color=num.harmonic)) + 
  #geom_line(aes(x=x, y=est, color=num.harmonic)) +
  #geom_segment(aes(x=x, xend=x, y=LT, yend = UT, color=num.harmonic), alpha=0.25) + 
  #geom_line(aes(x=x, y=estT, color=harmonic)) +
  #geom_point(aes(x=x, y=est, color=num.harmonic), alpha=0.5) + 
  labs(y="Percentage of ILI Patients", x="Time") + 
  theme_bw() + facet_grid(num.harmonic~.)



##### Predictions #####
true.data <- subset(get_cdc_data(2014), REGION=="Region 9")[1:4,]
true.data <- (true.data$TotILI/true.data$TotPat)

flu.model = dlmModTrig(om=52, q = 1, dV = 1000) + 
  dlmModPoly(2, dV = 0, dW = c(0, 0), m0=rep(0,2), C0=100*diag(2))
flu.model2 = dlmModTrig(om=52, q = 2, dV = 1000) + 
  dlmModPoly(2, dV = 0, dW = c(0, 2), m0=rep(0,2), C0=100*diag(2))


wk1pred1 <- wk1pred2 <- wk2pred1 <- wk2pred2 <- wk3pred1 <- wk3pred2 <- wk4pred1 <- wk4pred2 <- NULL
wk1state1 <- wk2state1 <- wk3state1 <- wk4state1 <- matrix(data=NA, ncol=4, nrow=nrow(res2$phi))
wk1state2 <- wk2state2 <- wk3state2 <- wk4state2 <- matrix(data=NA, ncol=6, nrow=nrow(res2$phi))
mse1wk1 <- mse1wk2 <- mse1wk3 <- mse1wk4 <- mse2wk1 <- mse2wk2 <- mse2wk3 <- mse2wk4 <- NULL
for(i in 1:nrow(res2$phi)){
  # Setting up theta_t for both models
  thetat1 <- c(res1$harmonic1a[i,209], res1$harmonic1b[i,209], res1$mu[i,209], res1$beta[i,209])
  thetat2 <- c(res2$harmonic1a[i,209], res2$harmonic1b[i,209], res2$harmonic2a[i,209], res2$harmonic2b[i,209], res2$mu[i,209], res2$beta[i,209])
  
  # Getting 1 wk ahead predictions
  wk1state1[i,] <- MASS::mvrnorm(1, flu.model$GG %*% thetat1, diag(c(0,0,0,res1$phi[i,2])) )
  wk1state2[i,] <- MASS::mvrnorm(1, flu.model2$GG %*% thetat2, diag(c(0,0,0,0,0,res2$phi[i,2])) )
  wk1pred1[i] <- rnorm(1, flu.model$FF %*% wk1state1[i,], res1$phi[i,1])
  wk1pred2[i] <- rnorm(1, flu.model2$FF %*% wk1state2[i,], res2$phi[i,1])
  
  # Getting 2 wk ahead predictions
  wk2state1[i,] <- MASS::mvrnorm(1, flu.model$GG %*%  wk1state1[i,], diag(c(0,0,0,res1$phi[i,2])) )
  wk2state2[i,] <- MASS::mvrnorm(1, flu.model2$GG %*% wk1state2[i,], diag(c(0,0,0,0,0,res2$phi[i,2])) )
  wk2pred1[i] <- rnorm(1, flu.model$FF %*%  wk2state1[i,], res1$phi[i,1])
  wk2pred2[i] <- rnorm(1, flu.model2$FF %*% wk2state2[i,], res2$phi[i,1])
  
  # Getting 3 wk ahead predictions
  wk3state1[i,] <- MASS::mvrnorm(1, flu.model$GG %*%  wk2state1[i,], diag(c(0,0,0,res1$phi[i,2])) )
  wk3state2[i,] <- MASS::mvrnorm(1, flu.model2$GG %*% wk2state2[i,], diag(c(0,0,0,0,0,res2$phi[i,2])) )
  wk3pred1[i] <- rnorm(1, flu.model$FF %*%  wk3state1[i,], res1$phi[i,1])
  wk3pred2[i] <- rnorm(1, flu.model2$FF %*% wk3state2[i,], res2$phi[i,1])
  
  # Getting 4 wk ahead predictions
  wk4state1[i,] <- MASS::mvrnorm(1, flu.model$GG %*%  wk3state1[i,], diag(c(0,0,0,res1$phi[i,2])) )
  wk4state2[i,] <- MASS::mvrnorm(1, flu.model2$GG %*% wk3state2[i,], diag(c(0,0,0,0,0,res2$phi[i,2])) )
  wk4pred1[i] <- rnorm(1, flu.model$FF %*%  wk4state1[i,], res1$phi[i,1])
  wk4pred2[i] <- rnorm(1, flu.model2$FF %*% wk4state2[i,], res2$phi[i,1])
}

mse1wk1 <- (mean(inv.logit(wk1pred1)) - true.data[1])^2
mse1wk2 <- (mean(inv.logit(wk2pred1)) - true.data[2])^2
mse1wk3 <- (mean(inv.logit(wk3pred1)) - true.data[3])^2
mse1wk4 <- (mean(inv.logit(wk4pred1)) - true.data[4])^2
mse2wk1 <- (mean(inv.logit(wk1pred2)) - true.data[1])^2
mse2wk2 <- (mean(inv.logit(wk2pred2)) - true.data[2])^2
mse2wk3 <- (mean(inv.logit(wk3pred2)) - true.data[3])^2
mse2wk4 <- (mean(inv.logit(wk4pred2)) - true.data[4])^2



