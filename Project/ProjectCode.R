##### Load Libraries #####
library("ggplot2")
library("cdcfluview")
library("plyr")
library("dplyr")
library("reshape2")
library("boot")
library("dlm")
source("../../ServerScripts/R/CreateCDCDF.R")


##### Setup Data #####
# OGflu <- get_cdc_data(year=2016)
# names(OGflu) <- tolower(names(OGflu))
# flu <- OGflu$totili
# ggplot(OGflu) + geom_point(aes(x=week, y=totili)) + facet_wrap(~region)

fludat <- get_flu_data(years=2010:2013)
fludat$REGION <- as.factor(fludat$REGION)
fludat <- dplyr::filter(fludat, REGION=="Region 9")
names(fludat) <- tolower(names(fludat))
fludat$time <- 1:nrow(fludat)
ggplot(fludat) + geom_point(aes(x=week, y=ilitotal/`total patients`)) + facet_wrap(~year)
ggplot(fludat) + geom_point(aes(x=time, y=ilitotal/`total patients`))
flu <- as.numeric(unlist(fludat[,c("ilitotal")]))/as.numeric(unlist(fludat[,c("total patients")]))
cflu <- (flu - mean(flu))
fluts <- ts(fludat$ilitotal/fludat$`total patients`, start=c(2010, 40), frequency = 52)

##### Plot Time-Series #####
plot(decompose(fluts)) 




##### Run Seasonal Analysis #####
flu.model = dlmModSeas(frequency=52, dV=2.315^2,dW=c(0, rep(0,50))) + dlmModPoly(2, dV = 2.315^2, dW = c(0, 2.315^(-10)))
flu.filter = dlmFilter(log(flu), flu.model)
bsample <- dlm::dlmBSample(flu.filter)



##### Find MLEs for starting values w/ no change #####
buildFun <- function(x){ dlmModSeas(frequency=52, dV = exp(x[1]), dW=rep(0,51), m0=rep(0,51), C0=100*diag(51)) +
    dlmModPoly(2, dV = 0, dW = c(0, exp(x[2])), m0=rep(0,2), C0=100*diag(2)) }
fit <- dlmMLE(logit(fluts), parm = c(1,1), build = buildFun)
sval <- list(vs=exp(fit$par[1]), w=exp(fit$par[2]))

# Check for gen fit with mle
flu.model = dlmModSeas(frequency=52, dV=sval$vs,dW=c(0, rep(0,50)), m0=rep(0,51), C0=100*diag(51)) + 
  dlmModPoly(2, dV = 0, dW = c(0, sval$w), m0=rep(0,2), C0=100*diag(2))
flu.filter = dlmFilter(log(fluts), flu.model)
bsample <- dlm::dlmBSample(flu.filter)
tmpdat <- data.frame(x=seq(1, length(flu)), y=(flu), est=bsample[-1,1], est2 = bsample[-1,1]+bsample[-1,52])
ggplot(tmpdat) + 
  geom_line(aes(x=x, y=exp(est2))) +
  geom_point(aes(x=x, y=y)) + 
  labs(y="Percentage of ILI Patients", x="Time") + theme_bw()


# Assumes no change to seasonal components
mcmc1 <- function(y, niter, sval){
  # Summary Stats
  n <- length(y)
  
  # setup storage
  thetakeep <- matrix(data=NA, ncol=3*(n+1), nrow=niter)
  phikeep   <- matrix(data=NA, ncol=length(sval), nrow=niter)
  
  # setup initialization
  #sigmam  <- sval$vm
  sigmaw  <- sval$w
  sigmas  <- sval$vs
  
  for(i in 1:niter){
    if(i%%(niter/100) == 0){print(paste(100*i/niter,"% of iterations complete", sep=""))}
    # sample theta jointly via backwards sampling
    dlmM    <- dlm::dlmModSeas(frequency=52, dV=sigmas ,dW=rep(0,51), m0=rep(0,51), C0=100*diag(51)) + 
      dlmModPoly(2, dV = 0, dW = c(0, sigmaw), m0=rep(0,2), C0=100*diag(2))
    FiltD   <- dlm::dlmFilter(y, dlmM)
    bsample <- dlm::dlmBSample(FiltD)
    
    # sample phi
    #SSEm  <- sum((y - (bsample[-1,52] + bsample[-1,53]))^2)
    SSEw  <- sum(diff(bsample[,53])^2)
    SSEs  <- sum((y - (bsample[-1,1] + bsample[-1,52]))^2)
    
    sigmas <- MCMCpack::rinvgamma(1, (n/2)+2, (SSEs/2) + 0.007)
    #sigmam <- MCMCpack::rinvgamma(1, (n/2)+1, (SSEm/2))
    sigmaw <- MCMCpack::rinvgamma(1, (n/2)+2, (SSEw/2) + 0.0008)
    
    # save samples
    thetakeep[i,]  <- c(bsample[,1], bsample[,52], bsample[,53])
    phikeep[i,]    <- c(sigmas, sigmaw)
  }
  
  # end of sampling mechanism
  return(list("theta"=thetakeep, "phi"=phikeep))
}
time <- system.time(res <- mcmc1(logit(fluts), 1000, sval))
save(res, time, file="SeasLTsamples.Rdata")
bires <- lapply(res, head, n = 1000)


ggplot() + geom_line(aes(x=1:length(res$theta[,1]), y=res$theta[,1]))
ggplot() + geom_line(aes(x=1:length(res$theta[,210]), y=res$theta[,210]))


efs <- data.frame(value=c(mcmcse::ess(bires$theta)))
efs$variable <- rep(NA, length(efs))
efs$variable[1:209] <- "theta"
efs$variable[-c(1:209)] <- "mean"
ggplot(efs) + geom_histogram(aes(x=value)) + facet_wrap(~variable) + theme_bw()


mphi <- melt(matrix(as.numeric(res$phi[,]), ncol=2, byrow=F))
mphi$Var2 <- factor(mphi$Var2)
levels(mphi$Var2) <- c("sigma_s", "sigma_w")
ggplot(mphi) + geom_line(aes(x=Var1, y=value)) + facet_wrap(~Var2, scales = "free") + 
  labs(x="Iterations", y="") + theme_bw()

gd <- data.frame(value=as.numeric(coda::geweke.diag(coda::mcmc(bires$theta))$z))
gd$variable <- rep(NA, length(gd))
gd$variable[1:209] <- "theta"
gd$variable[-c(1:209)] <- "mean"
ggplot(gd) + geom_histogram(aes(x=value)) + facet_wrap(~variable) + theme_bw()

d1 <- melt(data.frame(sigma_s=res$phi[,1], sigma_w=res$phi[,2]))
dinvgamma = function(x, a, b) dgamma(1/x,a,b)/x^2
dsqrtinvgamma = function(x, a, b) dinvgamma(x^2, a, b)*2*x

ggplot(d1) + 
  geom_histogram(aes(x = sqrt(value), y = ..density..), bins = 100) + 
  facet_grid( ~ variable, scales = "free") +
  stat_function(fun = dinvgamma, color = "red", args = list(a = 1, b = 0.007)) +
  theme_bw()


n <- length(flu)
L1 <- U1 <- esttheta <- NULL
for(i in 1:n){
  L1[i] <- quantile(exp(bires$theta[,i+1] + bires$theta[,i+210]), probs=0.025)
  U1[i] <- quantile(exp(bires$theta[,i+1] + bires$theta[,i+210]), probs=0.975)
  esttheta[i] <- quantile(exp(bires$theta[,i+1] + bires$theta[,i+210]), probs=0.5)
}
dat1 <- data.frame(L=L1, U=U1, x=seq(1,n), ogy=(flu), est=esttheta)

ggplot(dat1) + 
  geom_ribbon(aes(x=x, ymin=L, ymax = U), alpha=0.5) + 
  geom_line(aes(x=x, y=est)) +
  geom_point(aes(x=x, y=ogy)) + 
  labs(y="Percentage of ILI Patients", x="Time") + theme_bw()










##### Find MLEs for Harmonics and Linear Trend #####
buildFun <- function(x){ dlmModTrig(om=52, q = 1, dV = exp(x[1])) +
    dlmModPoly(2, dV = 0, dW = c(0, exp(x[2])), m0=rep(0,2), C0=100*diag(2)) }
fit <- dlmMLE(logit(fluts), parm = c(0.002,0.007), build = buildFun)
sval1 <- list(s=exp(fit$par[1]), b=exp(fit$par[2]))

flu.model = dlmModTrig(om=52, q = 1, dV = sval1$s) + 
  dlmModPoly(2, dV = 0, dW = c(0, sval1$b), m0=rep(0,2), C0=100*diag(2))
flu.filter = dlmFilter(logit(fluts), flu.model)
flu.forecast1 <- dlmForecast(flu.filter, nAhead=4)
bsample <- dlm::dlmBSample(flu.filter)



Rt1 <- flu.model$GG %*% dlmSvd2var(flu.filter$U.C, flu.filter$D.C)[[209]] %*% t(flu.model$GG) + flu.model$W
at1 <- flu.model$GG %*% bsample[209,]

ft1 <- flu.model$FF %*% at1
Qt1 <- flu.model$FF %*% Rt1 %*% t(flu.model$FF) + flu.model$V
rnorm(1, ft1, Qt1)

dlmForecast(flu.filter, nAhead=4)$f

tmpdat <- data.frame(x=seq(1, length(flu)), y=c(flu), est = c(bsample[-1,1]+bsample[-1,3]), harmonic="1")







buildFun <- function(x){ dlmModTrig(om=52, q = 2, dV = exp(x[1])) +
    dlmModPoly(2, dV = 0, dW = c(0, exp(x[2])), m0=rep(0,2), C0=100*diag(2)) }
fit <- dlmMLE(logit(fluts), parm = c(1,1), build = buildFun)
sval2 <- list(s=exp(fit$par[1]), b=exp(fit$par[2]))

flu.model = dlmModTrig(om=52, q = 2, dV = sval2$s) + 
  dlmModPoly(2, dV = 0, dW = c(0, sval2$b), m0=rep(0,2), C0=100*diag(2))
flu.filter = dlmFilter(logit(fluts), flu.model)
flu.forecast2 <- dlmForecast(flu.filter, nAhead=4)
bsample <- dlm::dlmBSample(flu.filter)
tmpdat2 <- data.frame(x=seq(1, length(flu)), y=c(flu), est = c(bsample[-1,1]+bsample[-1,3]+bsample[-1,5]), harmonic="2")

forecast1 <- data.frame(x=208:211, forecast=as.numeric(flu.forecast1$f), harmonic="1")
forecast2 <- data.frame(x=208:211, forecast=as.numeric(flu.forecast2$f), harmonic="2")

blank_data <- data.frame(x = c(0:(length(flu)+3)), y = seq(0,.05, length.out=length(flu)+4))

ggplot() + 
  geom_blank(data = blank_data, aes(x = x, y = y)) +
  geom_line(data=rbind(tmpdat, tmpdat2), aes(x=x, y=inv.logit(est), color=harmonic)) +
  geom_line(data=rbind(forecast1, forecast2), aes(x=x, y=inv.logit(forecast), color=harmonic)) +
  geom_point(data=rbind(tmpdat, tmpdat2), aes(x=x, y=y)) + 
  labs(y="Percentage of ILI Patients", x="Time") + theme_bw()


# Using 1 Harmonic
mcmc2 <- function(y, niter, sval){
  # Summary Stats
  n <- length(y)
  
  # setup storage
  thetakeep <- matrix(data=NA, ncol=(n+1)*3, nrow=niter)
  phikeep   <- matrix(data=NA, ncol=length(sval), nrow=niter)
  
  # setup initialization
  sigma  <- sval$s
  sigmab  <- sval$b
  
  for(i in 1:niter){
    if(i%%(niter/10) == 0){print(paste(100*i/niter,"% of iterations complete", sep=""))}
    # sample theta jointly via backwards sampling
    dlmM    <- dlm::dlmModTrig(om=52, q = 1, dV = sigma) + 
      dlm::dlmModPoly(2, dV = 0, dW = c(0, sigmab), m0=rep(0,2), C0=100*diag(2))
    FiltD   <- dlm::dlmFilter(y, dlmM)
    bsample <- dlm::dlmBSample(FiltD)
    
    # sample phi
    SSE  <- sum((y - (bsample[-1,1] + bsample[-1,3]))^2)
    SSEb  <- sum(diff(bsample[,4])^2)
    
    sigma <- MCMCpack::rinvgamma(1, (n/2)+2, (SSE/2) + 0.007)
    sigmab <- MCMCpack::rinvgamma(1, (n/2)+2, (SSEb/2) + 0.002)
    
    # save samples
    thetakeep[i,]  <- c(bsample[,1], bsample[,3], bsample[,4])
    phikeep[i,]    <- c(sigma, sigmab)
  }
  
  # end of sampling mechanism
  return(list("theta"=thetakeep, "phi"=phikeep))
}
time2 <- system.time(res2 <- mcmc2(logit(fluts), 10000, sval1))

# Using 2 Harmonics
mcmc3 <- function(y, niter, sval){
  # Summary Stats
  n <- length(y)
  
  # setup storage
  thetakeep <- matrix(data=NA, ncol=4*(n+1), nrow=niter)
  phikeep   <- matrix(data=NA, ncol=length(sval), nrow=niter)
  
  # setup initialization
  #sigmam  <- sval$vm
  sigma  <- sval$s
  sigmab  <- sval$b
  
  for(i in 1:niter){
    if(i%%(niter/10) == 0){print(paste(100*i/niter,"% of iterations complete", sep=""))}
    # sample theta jointly via backwards sampling
    dlmM    <- dlm::dlmModTrig(om=52, q = 2, dV = sigma) + 
      dlm::dlmModPoly(2, dV = 0, dW = c(0, sigmab), m0=rep(0,2), C0=100*diag(2))
    FiltD   <- dlm::dlmFilter(y, dlmM)
    bsample <- dlm::dlmBSample(FiltD)
    
    # sample phi
    SSE  <- sum((y - (bsample[-1,1] + bsample[-1,3] + bsample[-1,5]))^2)
    SSEb  <- sum(diff(bsample[,6])^2)
    
    sigma <- MCMCpack::rinvgamma(1, (n/2)+2, (SSE/2) + 0.007)
    sigmab <- MCMCpack::rinvgamma(1, (n/2)+2, (SSEb/2) + 0.002)
    
    # save samples
    thetakeep[i,]  <- c(bsample[,1], bsample[,3], bsample[,5], bsample[,6])
    phikeep[i,]    <- c(sigma, sigmab)
  }
  
  # end of sampling mechanism
  return(list("theta"=thetakeep, "phi"=phikeep))
}
time3 <- system.time(res3 <- mcmc3(logit(fluts), 10000, sval2))

save(res2, res3, file="HarmLTSamples.RData")


bires <- lapply(res, head, n = 1000)

# Traceplots
ggplot() + geom_line(aes(x=1:length(res2$theta[,1]),   y=res2$theta[,1])) + labs(y="theta")
ggplot() + geom_line(aes(x=1:length(res2$theta[,210]), y=res2$theta[,210])) + labs(y="mu")
ggplot() + geom_line(aes(x=1:length(res2$theta[,419]), y=res2$theta[,419])) + labs(y="beta")

ggplot() + geom_line(aes(x=1:length(res3$theta[,1]),   y=res3$theta[,1])) + labs(y="theta1")
ggplot() + geom_line(aes(x=1:length(res3$theta[,210]), y=res3$theta[,210])) + labs(y="theta2")
ggplot() + geom_line(aes(x=1:length(res3$theta[,419]), y=res3$theta[,419])) + labs(y="mu")
ggplot() + geom_line(aes(x=1:length(res3$theta[,420]), y=res3$theta[,420])) + labs(y="mu")

# Effective Sample Sizes
efs <- data.frame(value=c(mcmcse::ess(res2$theta)))
efs$variable <- rep(NA, length(efs))
efs$variable[1:209] <- "theta"
efs$variable[210:418] <- "mu"
efs$variable[-c(1:418)] <- "beta"
ggplot(efs) + geom_histogram(aes(x=value)) + facet_wrap(~variable) + theme_bw()

efs <- data.frame(value=c(mcmcse::ess(res3$theta)))
efs$variable <- rep(NA, length(efs))
efs$variable[1:209] <- "theta1"
efs$variable[210:418] <- "theta2"
efs$variable[419:627] <- "mu"
efs$variable[628:836] <- "beta"
ggplot(efs) + geom_histogram(aes(x=value)) + facet_wrap(~variable) + theme_bw()


# Trace plots of sigma values
mphi <- melt(matrix(as.numeric(res2$phi[,]), ncol=2, byrow=F))
mphi$Var2 <- factor(mphi$Var2)
levels(mphi$Var2) <- c("sigma", "sigma_beta")
mphi$harmonic <- "1"

mphi2 <- melt(matrix(as.numeric(res3$phi[,]), ncol=2, byrow=F))
mphi2$Var2 <- factor(mphi2$Var2)
levels(mphi2$Var2) <- c("sigma", "sigma_beta")
mphi2$harmonic <- "2"

ggplot(rbind(mphi,mphi2)) + geom_line(aes(x=Var1, y=value)) + facet_wrap(harmonic~Var2, scales = "free") + 
  labs(x="Iterations", y="") + theme_bw()

# Geweke Diagnostics
gd <- data.frame(value=as.numeric(coda::geweke.diag(coda::mcmc(res$theta))$z))
gd$variable <- rep(NA, length(gd))
gd$variable[1:209] <- "theta"
gd$variable[210:418] <- "mu"
gd$variable[-c(1:418)] <- "mean"
ggplot(gd) + geom_histogram(aes(x=value)) + facet_wrap(~variable) + theme_bw()

# prior to posterior plots
d1 <- melt(data.frame(sigma=res2$phi[,1], sigma_beta=res2$phi[,2]))
dinvgamma = function(x, a, b) dgamma(1/x,a,b)/x^2
dsqrtinvgamma = function(x, a, b) dinvgamma(x^2, a, b)*2*x

ggplot(d1) + 
  geom_histogram(aes(x = sqrt(value), y = ..density..), bins = 100) + 
  facet_grid( ~ variable, scales = "free") +
  stat_function(fun = dsqrtinvgamma, color = "red", args = list(a = 2, b = 0.007)) +
  theme_bw()

# Model Fits
n <- length(flu)
L1 <- U1 <- esttheta <- L2 <- U2 <- esttheta2 <- NULL
for(i in 1:n){
  L1[i] <- quantile(exp(res2$theta[,i+1] + res2$theta[,i+210]), probs=0.025)
  U1[i] <- quantile(exp(res2$theta[,i+1] + res2$theta[,i+210]), probs=0.975)
  esttheta[i] <- quantile(exp(res2$theta[,i+1] + res2$theta[,i+210]), probs=0.5)
  L2[i] <- quantile(exp(res3$theta[,i+1] + res3$theta[,i+210] + res3$theta[,i+419]), probs=0.025)
  U2[i] <- quantile(exp(res3$theta[,i+1] + res3$theta[,i+210] + res3$theta[,i+419]), probs=0.975)
  esttheta2[i] <- quantile(exp(res3$theta[,i+1] + res3$theta[,i+210] + res3$theta[,i+419]), probs=0.5)
}
dat1 <- data.frame(L=L1, U=U1, x=seq(1,n), ogy=(flu), est=esttheta, harmonic="1")
dat2 <- data.frame(L=L2, U=U2, x=seq(1,n), ogy=(flu), est=esttheta2, harmonic="2")
ggplot(rbind(dat1, dat2)) + 
  geom_ribbon(aes(x=x, ymin=L, ymax = U, fill=harmonic), alpha=0.5) + 
  geom_line(aes(x=x, y=est, color=harmonic)) +
  geom_point(aes(x=x, y=ogy)) + 
  labs(y="Percentage of ILI Patients", x="Time") + 
  theme_bw() + facet_grid(harmonic~.)














##### Check Number of Harmonics #####
m <- lm(flu-mean(flu)~factor(fludat$week)-1)
flu.model = dlmModSeas(52,dV=2.315^2,dW=rep(0,51))
flu.filter = dlmFilter(flu-mean(flu), flu.model)
n = length(flu) + 1 # Due to theta_0
data.frame(mean = c(-sum(flu.filter$m[n,]), rev(flu.filter$m[n,])),
           lm_mean = m$coefficients)

d = plyr::ddply(data.frame(harmonics=1:6), .(harmonics), function(x) {
  flu.model = dlmModTrig(52,x$harmonics,dV=2.315^2,dW=rep(0,51))
  flu.filter = dlmFilter(flu-mean(flu),flu.model)
  forecast = dlmForecast(flu.filter, 52)
  data.frame(month = factor(fludat$week),
             effect = forecast$f)
})


d$harmonics_f = factor(d$harmonics)
ggplot(d, aes(month, effect, color=harmonics_f, shape=harmonics_f, group=harmonics_f)) +
  geom_point() +
  geom_line() +
  theme_bw()




mod <- dlmModTrig(om=52, q = 2, dV = 1000) +
  dlmModPoly(1, dV = 1000, dW = c(0))

flu.smooth <- dlmSmooth(flu, mod)
y <- cbind(fluts,
           tcrossprod(dropFirst(flu.smooth$s[, c(1, 3, 5)]),
                      matrix(c(0, 0, 1, 1, 1, 0), nr = 2,
                             byrow = TRUE)))
colnames(y) <- c("Flu", "Level", "Periodic")
plot(y, yax.flip = TRUE, oma.multi = c(2, 0, 1,4))



