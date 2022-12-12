# clear memory
rm(list = ls())
library(hdbinseg)
library(MARSS)
library(tsutils)

# Set working directory
sw_common <- read.csv("sw_2016_58vars_common.csv",header=F,sep=",",stringsAsFactors = FALSE)

# transpose data so time goes across columns
z <- t(sw_common)

## Get the number of factors using Bai-Ng
bn_factor <- get.factor.model(z, r.max = 11, ic = c("ah", "bn")[2], ic.op = 2, normalisation = TRUE)
bn_factor$r.hat

## Get the number of factors using Ah-Horenstein
ah_factor <- get.factor.model(z, r.max = 11, ic = c("ah", "bn")[1], ic.op = 2, normalisation = TRUE)
ah_factor$r.hat

sw <- read.csv("sw_2016_58vars_NAs.csv",header=F,sep=",",stringsAsFactors = FALSE)

# convert to ts object
sw.ts <- ts(sw, 
            start = c(1959, 1), 
            end = c(2014, 4), 
            frequency = 4)

## transpose data so time goes across columns
t_sw <- zscore(t(sw.ts), mean.only=F)

## get number of time series
N_ts <- dim(t_sw)[1]
N_ts
## get length of time series
TT <- dim(t_sw)[2]
TT

## Set-up the matrices required
## 'ZZ' is loadings matrix
ZZ <- matrix(c(rep(1,58),rep(0,174)),58,4)
ZZ <- matrix(list(0), 58, 4)
ZZ[,1] <- paste0("1_", seq(58))
ZZ

## RR is var-cov matrix for obs errors
RR = diag(0,58,58)
diag(RR)<-"r"

## QQ is the var-cov matrix of the state equations
QQ <- matrix(list(0),4,4)
QQ[1,1] <- "q"

## BB is the factor loadings matrix
BB <- matrix(list("b11",1,0,0,"b12",0,1,0,"b13",0,0,1,"b14",0,0,0),4,4)
BB
aa <- matrix(0,58,1) ## constants
uu <- matrix(0,4,1) ## constants
pi <- matrix(mean(t_sw, na.rm=T),4,1) # prior parameter
V <- matrix(0,4,4) # prior parameter

## Estimate the DFM
mod_list <- list(Z = ZZ, A = aa, R = RR, B = BB, 
                 U = uu, Q = QQ, x0 = pi, V0 = V)
con_list <- list(maxit = 3000, allow.degen = TRUE, conv.test.slope.tol=0.05)
init_model <- MARSS(y = t_sw[,5:216], model = mod_list, control = con_list)

# Pass initial values from max-E algorithm to do the final MLE estimation
# inits <- coef(init_model, what = "par")
# bfgs_dfm_model <- MARSS(y = t_sw, model = mod_list, inits=inits, method = "BFGS")

# Obtain factor
ss_ <- tsSmooth(init_model, type = "xtT", interval = "none") # fixed-interval smoothing on a univariate time series
states_ <- ss_$.estimate
factor1 <- states_[1:216]  # 216 is the number of observations
fac_fs <- rowSums(lagmatrix(factor1,c(0,1,2,3)))
factor1.ts <- ts(fac_fs, 
                 start = c(1959, 1), 
                 end = c(2012, 4), 
                 frequency = 4)
plot(factor1.ts)

# forecasts
fr <- predict(init_model,  newdata = list(t=217:224, y=t_sw[,217:224]), type = "ytt1")
plot(fr)
