## M_VCSM_tutorial.R
#############################################################################
## Description: A step-by-step implementation of M-VCSM and the associated  
## procedures described in "Multilevel Varying Coefficient Spatiotemporal Model".  
#############################################################################
## Functions implemented: 
## M_VCSM_simulation.R, M_VCSM_decomposition.R, M_VCSM_MCMC.R, M_VCSM_inference.R
#############################################################################
## Tutorial Outline:
## 1. Simulate hospitalization rate outcome data (M_VCSM_simulation.R)
## 2. Perform M-VCSM estimation (M_VCSM_decomposition.R, M_VCSM_MCMC.R)
## 3. Inference on multilevel risk factors and prediction of region- and facility-specific deviations (M_VCSM_inference.R)
## 4. Visualization of M-VCSM results
#############################################################################

# Install missing packages
list.of.packages <- c("MASS", "caTools", "locpol", "KernSmooth", "fANCOVA", "mgcv", "mvtnorm", "spdep", "refund")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages) 

# Load packages  
library(MASS)
library(caTools)
library(locpol)
library(KernSmooth)
library(fANCOVA)
library(mgcv)
library(mvtnorm)
library(spdep)
library(refund)


#############################################################################
# 1. Simulate hospitalization rate outcome data
#############################################################################

# NOTE: Generating one dataset with 49 regions and 4-20 facilities per region will take approximately three minutes.

# Simulate one dataset from the simulation design described in Web Appendix E with 49 regions and 4-20 facilities per region
data.G <- M_VCSM_simulation(numRegion = 2, numFacility = 1, sigma2 = 0.2)  # M_VCSM_simulation.R

# Data frame used for estimation and inference
data <- data.G[[1]]

# Adjacency matrix from the map
Adj.Mat <- data.G[[2]]

# Save the underlying true values of the generated data
data.T <- data.G[[3]]

#############################################################################
# 2. Perform M-VCSM estimation
#############################################################################

# NOTE: Performing M-VCSM estimation steps 1-4 (FPCA decomposition) with 49 regions and 4-20 facilities per region will take approximately 1 minute.

FPCAout <- M_VCSM_decomposition(data = data)  # MST_FM_decompostion.R

# NOTE: Performing M-VCSM estimation steps 5-6 (MCMC estimation) with 49 regions and 4-20 facilities per region will take approximately 20 minutes.

MCMCout <- M_VCSM_MCMC(FPCAout = FPCAout, data = data, AdjMat = Adj.Mat) # MST_FM_MCMC.R

#############################################################################
# 3. Inference on multilevel risk factors and prediction of region- and facility-specific deviations
#############################################################################

Inferenceout <- M_VCSM_inference(FPCAout = FPCAout, data = data, MCMCout = MCMCout)

#############################################################################
# 4. Visualization of MST-FM results
#############################################################################  

# Define the grid points used for the mean function and eigenfunctions
ngrid <- 24 # 2 year follow up
gridPoints <- seq(0,1,length.out = ngrid)

# True eigenfunctions and varying coefficient functions (VCFs)
# first-level eigenfunctions
psi1.1 <- function(x){
    return(sqrt(2)*sin(2*pi*x))
}
psi1.2 <- function(x){
    return(sqrt(2)*cos(2*pi*x))
}

# second-level eigenfunctions
psi2.1 <- function(x){
    return(sqrt(3)*(2*x-1))
}
psi2.2 <- function(x){
    return(sqrt(5)*(6*x^2-6*x+1))
}

# VCFs
beta0 <- function(x){
    return(1-x)
}

beta1 <- function(x){
    return(2*x)
}

beta2 <- function(x){
    return(-4*x^2+4*x-1)
}

theta1 <- function(x){
    return(2*x^2-3*x+1)
}

theta2 <- function(x){
    return(10*x^3-16*x^2+7*x-1)
}

# Plot estimates of varying coefficient functions and 95% confidence bands
par(mfrow=c(3,2))

# Plot beta0(t)
beta0CB <- Inferenceout$beta.CB[,c(1,2)]
plot(gridPoints, Inferenceout$beta[,1],'l',col="black",lwd=2,main="(a)", ylim = range(beta0CB),xlim=c(0,1),xaxs = "i",cex.main=1.8,xlab = "",ylab = "",cex.axis=1.4) 
title(xlab = "t",ylab=expression(widehat(beta)[0](t)), line=2, cex.lab=1.5)
lines(gridPoints,beta0(gridPoints),lwd=2,lty=1, col = "grey52")
lines(gridPoints,beta0CB[,1], lwd=2,lty=2)
lines(gridPoints,beta0CB[,2], lwd=2,lty=2)

# Plot beta1(t)
beta1CB <- Inferenceout$beta.CB[,c(3,4)]
plot(gridPoints,Inferenceout$beta[,2],'l',col="black",lwd=2,main="(b)", ylim = range(beta1CB),xlim=c(0,1),xaxs = "i",cex.main=1.8,xlab = "",ylab = "",cex.axis=1.4) 
title(xlab = "t",ylab=expression(widehat(beta)[1](t)), line=2, cex.lab=1.5)
lines(gridPoints,beta1(gridPoints),lwd=2,lty=1, col = "grey52")
lines(gridPoints,beta1CB[,1], lwd=2,lty=2)
lines(gridPoints,beta1CB[,2], lwd=2,lty=2)

# Plot beta2(t)
beta2CB <- Inferenceout$beta.CB[,c(5,6)]
plot(gridPoints,Inferenceout$beta[,3],'l',col="black",lwd=2,main="(c)", ylim = range(beta2CB),xlim=c(0,1),xaxs = "i",cex.main=1.8,xlab = "",ylab = "",cex.axis=1.4) 
title(xlab = "t",ylab=expression(widehat(beta)[2](t)), line=2, cex.lab=1.5)
lines(gridPoints,beta2(gridPoints),lwd=2,lty=1, col = "grey52")
lines(gridPoints,beta2CB[,1], lwd=2,lty=2)
lines(gridPoints,beta2CB[,2], lwd=2,lty=2)

# Plot theta1(t)
theta1CB <- Inferenceout$theta.CB[,c(1,2)]
plot(gridPoints,Inferenceout$theta[,1],'l',col="black",lwd=2,main="(d)", ylim = range(theta1CB),xlim=c(0,1),xaxs = "i",cex.main=1.8,xlab = "",ylab = "",cex.axis=1.4) 
title(xlab = "t",ylab=expression(widehat(theta)[1](t)), line=2, cex.lab=1.5)
lines(gridPoints,theta1(gridPoints),lwd=2,lty=1, col = "grey52")
lines(gridPoints,theta1CB[,1], lwd=2,lty=2)
lines(gridPoints,theta1CB[,2], lwd=2,lty=2)

# Plot theta2(t)
theta2CB <- Inferenceout$theta.CB[,c(3,4)]
plot(gridPoints,Inferenceout$theta[,2],'l',col="black",lwd=2,main="(e)", ylim = range(theta2CB),xlim=c(0,1),xaxs = "i",cex.main=1.8,xlab = "",ylab = "",cex.axis=1.4) 
title(xlab = "t",ylab=expression(widehat(theta)[2](t)), line=2, cex.lab=1.5)
lines(gridPoints,theta2(gridPoints),lwd=2,lty=1, col = "grey52")
lines(gridPoints,theta2CB[,1], lwd=2,lty=2)
lines(gridPoints,theta2CB[,2], lwd=2,lty=2)

## Align eigenfunctions
# Note: the estimated eigenfunctions may flip signs. Here we match the sign of the
# estimated eigenfunctions to the true eigenfunctions for plotting  
MSDE.eifun <- function(x,y){
    err1 <- trapz(gridPoints, (x-y)^2)
    err2 <- trapz(gridPoints, (x+y)^2)
    return(c(min(err1, err2),which.min(c(err1, err2))))
}

psi1.1.t <- psi1.1(gridPoints)
psi1.2.t <- psi1.2(gridPoints)
psi2.1.t <- psi2.1(gridPoints)
psi2.2.t <- psi2.2(gridPoints)

psi1.1.Est <- FPCAout$psi1[,1]
psi1.2.Est <- FPCAout$psi1[,2]
psi2.1.Est <- FPCAout$psi2[,1]
psi2.2.Est <- FPCAout$psi2[,2]


MSDEeifun1.2 <- MSDE.eifun(psi1.1.Est, psi1.1.t)
MSDEeifun2.2 <- MSDE.eifun(psi1.2.Est, psi1.2.t)

MSDEeifun1.fac <- MSDE.eifun(psi2.1.Est, psi2.1.t)
MSDEeifun2.fac <- MSDE.eifun(psi2.2.Est, psi2.2.t)

fun1.sign.2 <- MSDEeifun1.2[2] * (-2) + 3
psi1.1.Est <- psi1.1.Est * fun1.sign.2
fun2.sign.2 <- MSDEeifun2.2[2] * (-2) + 3
psi1.2.Est <- psi1.2.Est * fun2.sign.2

fun1.sign.fac <- MSDEeifun1.fac[2] * (-2) + 3
psi2.1.Est <- psi2.1.Est * fun1.sign.fac
fun2.sign.fac <- MSDEeifun2.fac[2] * (-2) + 3
psi2.2.Est <- psi2.2.Est * fun2.sign.fac

# Plot estimates of eigenfunctions
par(mfrow = c(2,2))

# First-level eigenfunctions
plot(gridPoints, psi1.1.t,"l", ylim = c(-1.4,1.4), xaxs = "i", main = "(a)",cex.main=2,xlab = "", ylab = "", lwd = 2)
title(xlab = "t",ylab=expression(paste(widehat(psi)[1]^(1),(t))), line=2, cex.lab=1.6)
lines(gridPoints, psi1.1.Est, col = "grey", lwd = 2, lty = 1)

plot(gridPoints, psi1.2.t,"l", ylim = c(-1.4,1.5), xaxs = "i", main = "(b)",cex.main=2,xlab = "", ylab = "", lwd = 2)
title(xlab = "t",ylab=expression(paste(widehat(psi)[2]^(1),(t))), line=2, cex.lab=1.6)
lines(gridPoints, psi1.2.Est, col = "grey", lwd = 2, lty = 1)

# Second-level eigenfunctions
plot(gridPoints, psi2.1.t,"l", ylim = c(-1.9,1.9), xaxs = "i", main = "(c)",cex.main=2,xlab = "", ylab = "", lwd = 2)
title(xlab = "t",ylab=expression(paste(widehat(psi)[1]^(2),(t))), line=2, cex.lab=1.6)
lines(gridPoints, psi2.1.Est, col = "grey", lwd = 2, lty = 1)

plot(gridPoints, psi2.2.t,"l", ylim = c(-1.1,2.3), xaxs = "i", main = "(d)",cex.main=2,xlab = "", ylab = "", lwd = 2)
title(xlab = "t",ylab=expression(paste(widehat(psi)[2]^(2),(t))), line=2, cex.lab=1.6)
lines(gridPoints, psi2.2.Est, col = "grey", lwd = 2, lty = 1)

par(mfrow = c(2,2))
# Region-specific deviation prediction
# True region-specific deviation functions
nFac <- aggregate(data.T$fid, by = list(data.T$rid), FUN = length)[,2] / ngrid
R.traj.T <- matrix(data.T$r.eff1 + data.T$r.eff2, nrow = ngrid)[,cumsum(nFac)]
# Predicted region-specific deviation functions
R.traj.Est <- Inferenceout$R.traj.Est
# Visualization of the region-specific trajectory predictions
matplot(R.traj.T, lty=1, type="l", ylim=range(R.traj.T, R.traj.Est), main = "True region-specific deviations")
matplot(R.traj.Est, lty=1, type="l", ylim=range(R.traj.T, R.traj.Est), main = "Estimated region-specific deviations")

# Facility-specific deviation prediction
# True facility-specific deviation functions
F.traj.T <- matrix(data.T$f.eff1 + data.T$f.eff2, nrow = ngrid)
# Predicted facility-specific deviation functions
F.traj.Est <- Inferenceout$F.traj.Est
# Visualization of the facility-specific trajectory predictions
matplot(F.traj.T, lty=1, type="l", ylim=range(F.traj.T, F.traj.Est), main = "True facility-specific deviations")
matplot(F.traj.Est, lty=1, type="l", ylim=range(F.traj.T, F.traj.Est), main = "Estimated facility-specific deviations")
