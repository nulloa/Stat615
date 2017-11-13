##### Load Libraries #####
library("dlm")
library("cdcfluview")
library("plyr")
library("dplyr")
library("reshape2")
library("ggplot2")
source("../ServerScripts/R/CreateCDCDF.R")


##### Setup Data #####
OGflu <- get_cdc_data(year=2016)
names(OGflu) <- tolower(names(OGflu))
flu <- OGflu$totili
ggplot(OGflu) + geom_point(aes(x=week, y=totili)) + facet_wrap(~region)

fludat <- get_flu_data(years=2010:2013)
fludat$REGION <- as.factor(fludat$REGION)
fludat <- dplyr::filter(fludat, REGION=="Region 6")
names(fludat) <- tolower(names(fludat))
fludat$time <- 1:nrow(fludat)
ggplot(fludat) + geom_point(aes(x=week, y=ilitotal)) + facet_wrap(~year)
ggplot(fludat) + geom_point(aes(x=time, y=ilitotal))
flu <- as.numeric(unlist(fludat[,c("ilitotal")]))
fluts <- ts(fludat$ilitotal, start=c(2010, 40), frequency = 52)

##### Plot Time-Series #####
plot(decompose(fluts)) 






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


##### Find MLEs for starting values #####
buildFun <- function(x){ dlmModPoly(2, dV = exp(x[1]), dW=c(0,exp(x[2])), m0=c(30,9), C0=matrix(c(100^2,0,0,1^2), 2, 2)) }
fit <- dlmMLE(flu, parm = c(1,1), build = buildFun)
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

time <- system.time(res <- mcmc1(flu, 9999, sval))


efs <- data.frame(value=c(mcmcse::ess(res$theta)))
efs$variable <- rep(NA, length(efs))
efs$variable[1:87] <- "mu"
efs$variable[-c(1:87)] <- "Beta"
ggplot(efs) + geom_histogram(aes(x=value)) + facet_wrap(~variable) + theme_bw()

mphi <- melt(res$phi)
mphi$Var2 <- factor(mphi$Var2)
levels(mphi$Var2) <- c("sigma", "sigma_beta")
ggplot(mphi) + geom_line(aes(x=Var1, y=value)) + facet_wrap(~Var2, scales = "free") + 
  labs(x="Iterations", y="") + theme_bw()


gd <- data.frame(value=as.numeric(coda::geweke.diag(coda::mcmc(res$theta))$z))
gd$variable <- rep(NA, length(gd))
gd$variable[1:87] <- "mu"
gd$variable[-c(1:87)] <- "Beta"
ggplot(gd) + geom_histogram(aes(x=value)) + facet_wrap(~variable) + theme_bw()


d1 <- melt(data.frame(sigma=res$phi[,1], sigma_beta=res$phi[,2]))
d <- d1

dinvgamma = function(x, a, b) dgamma(1/x,a,b)/x^2
dsqrtinvgamma = function(x, a, b) dinvgamma(x^2, a, b)*2*x

ggplot(d1) + 
  geom_histogram(aes(x = sqrt(value), y = ..density..), bins = 100) + 
  facet_grid( ~ variable, scales = "free") +
  stat_function(fun = dinvgamma, color = "red", args = list(a = 1, b = 1)) +
  theme_bw()


n <- length(lake)
L1 <- U1 <- L2 <- U2 <- NULL
for(i in 1:n){
  L1[i] <- quantile(res$theta[,i+1], probs=0.025)
  U1[i] <- quantile(res$theta[,i+1], probs=0.975)
}
dat1 <- data.frame(L=L1, U=U1, x=OGlake$year, ogy=lake)

ggplot(dat1) + 
  geom_ribbon(aes(x=x, ymin=L, ymax = U)) + 
  geom_point(aes(x=x, y=ogy)) + 
  theme_bw() + 
  labs(x="Year", y="Precipitation")

