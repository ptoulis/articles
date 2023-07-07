#' 7/5, Greek Stochastics Conference
#' Fbart examples. Panos Toulis panos.toulis@chicagobooth.edu

# Sharp test [0.17, 0.42]
rm(list=ls())

library(experimentr)
lynn$tip = lynn$tip/lynn$groupsize
library(dplyr)

NRANDS = 10000

# create dataset.
D = lynn %>% mutate(Y=tip, Z=crouch) %>% select(Y,Z)

# initialize bcf
library(bcf)
X = cbind(1, as.factor(lynn$groupsize), lynn$daytime, lynn$paid_by_credit_card)
n = nrow(X)
pihat = rep(0.5, n)
out = bcf(y=D$Y, z=D$Z, x_control=X, pihat=pihat, nburn=1000, nsim=2500)
sigma_est = median(out$sigma)

mu = colMeans(out$tau) # mean
S = cov(out$tau) + 2*sigma_est^2*diag(n) # cov matrix
# library(covglasso)  # ucomment for sparse estimation
# S = covglasso(S=S, n=nrow(S))$sigma 

# p-value finction
pval = function(tvals, tobs) {
  m = sum(tvals > tobs)
  r = length(tvals)
  p = (1/(r+1)) * (m + 1)
  2*min(p, 1-p)
}


#' FRT for H0: tau_vector = taus.
frt_tau = function(taus) {
  stopifnot(length(taus)==nrow(lynn))
  Y = D$Y
  Z = D$Z
  
  # Y = Y0 + taus * Z
  Y0 = Y - (taus*Z)
  Y1 = Y0 + taus
  
  # diff in means
  tstat = function(y1, z1){
    mean(y1[z1==1]) - mean(y1[z1==0])
  }
  
  Tobs = tstat(Y,Z)
  
  Tvals = replicate(NRANDS, {
    z_r = sample(Z)
    y_r = Y1*z_r + Y0*(1-z_r)
    tstat(y_r, z_r)
  })
  
  
  #hist(Tvals, Tobs, breaks=30)
  #abline(v=Tobs, col="red", lwd=2)
  pval(Tvals, Tobs)
}

#' Fisher sharp null for H0: tau_i = q, for all i.
#
frt = function(q) {
  
  Y = D$Y
  Z = D$Z
  
  taus = rep(q, length(Y))
  frt_tau(taus)
}

# Simple implementation of bcf BART.
#
fbart = function(ATE, nreps=500)  {
  
  # Conditional distribution 
  # H0: mean(tau) = ATE
  a = rep(1, n)/n
  # Calculate distribution of  x | a'x=0 given that x ~ N(mu, S)
  # Say N(mu_0, S_0)
  S12 = as.vector(S %*% a)
  mu2 = sum(a*mu)
  v2 = as.numeric(t(a) %*% S %*% a)
  
  mu_0 = mu + S12 * (ATE - sum(a*mu))/v2
  S_0 = S - S12 %*% t(S12) / v2
 
  # sup over the conditional distribution
  library(mvtnorm)
  x = rmvnorm(nreps, mean=mu_0, sigma = S_0)
  sup_pval = 0
  for(i in 1:nreps) {
    p = frt_tau(x[i,])
    if(p > sup_pval) {
      print(paste("> New sup pval = ", p))
      sup_pval = p
    }
  }
  
}


# Slides results
set.seed(110)
# Sharp test: 
frt(0.170) # accepted
frt(0.424) # accepted
# Interval is [0.17, 0.424]. Anything outside is rejected.

frt(0.12) ## rejected.

# fbart.
fbart(0.12) # takes time
fbart(0.48) # takes time




