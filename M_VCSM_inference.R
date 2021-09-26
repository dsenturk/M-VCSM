M_VCSM_inference <- function(FPCAout, # output from function M_VCSM_decomposition, including the follows:
                             # Definition: n: number of regions, Ni: number of facilities in region i, T: number of time points
                             # L: number of first-level eigencomponents, M: number of second-level eigencomponents
                             # psi1: estimated region-level eigenfunctions (matrix of dimension L*T)
                             # psi2: estimated facility-level eigenfunctions (matrix of dimension M*T)
                             # sigma2: estimated measurement error variance, used as an initial value in the MCMC step (scalar)
                             # eval1: estimated region-level eigenvalues, used as initial values in the MCMC step (vector of length L)
                             # eval2: estimated facility-level eigenvalues, used as initial values in the MCMC step (vector of length M)
                             # lambda.prop: estimated omega, proportion of facility-level eigenvalues to the largest facility-level eigenvalue (vector of length M)
                             
                             data, # data.frame in long format with eight labeled columns (described below)
                             # and row length equal to the length of the vectorized observations across all 
                             # regions and facilities
                             # DATA.FRAME COLUMNS: 
                             # rid: region IDs (vector of length T*sum(Ni) )
                             # fid: facility IDs (vector of length T*sum(Ni))
                             # y: hospitalization rate data (vector of length T*sum(Ni))
                             # t: follow-up time (vector of length T*sum(Ni))                              
                             # x1: region-level covariate (vector of T*sum(Ni))
                             # x2: region-level covariate (vector of length T*sum(Ni))
                             # z1: facility-level covariate (vector of length T*sum(Ni))
                             # z2: facility-level covariate (vector of length T*sum(Ni))                             
                             
                             MCMCout # output from function M_VCSM_MCMC, including the follows:
                             # beta_mat: posterior samples of varying coefficient functions (list of 5)
                             # alpha: posterior samples of spatial variance parameter (matrix of dimension L*12000)
                             # nv: posterior samples of spatial correlation parameter (vector of length 12000)
                             # sigma2: posterior samples of measurement error variance (vector of length 12000)
                             # xi: posterior samples of region-specific PC scores (matrix of dimension n*12000)
                             # zeta: posterior samples of facility-specific PC scores (matrix of dimension sum(Ni)*12000)
                             # lambda: posterior samples of region-specific facility-level variances (matrix of dimension n*12000)
                             
                             
){
    
    #############################################################################
    ## Description: Function for obtaining inference for varying coefficient functions (VCFs), time-static parameter estimates 
    ##              and prediction of region- and facility-specific deviations
    ## Definition:  n: number of regions, Ni: number of facilities in region i, T: number of time points, 
    ##              L: number of first-level eigencomponents, M: number of second-level eigencomponents
    ## Args: see above
    ## Returns:     list()
    ##              beta: Estimated region-level risk factor effects and the intercept term (matrix of dimension 3*T) 
    ##              theta: Estimated facility-level risk factor effects (matrix of dimension 2*T)
    ##              beta.CB: Confidence band of region-level risk factor effects and the intercept term (matrix of dimension 6*T) 
    ##              theta.CB: Confidence band of facility-level risk factor effects (matrix of dimension 4*T)
    ##              R.traj.Est:  Region-specific deviation predictions (matrix of dimension n*T)
    ##              F.traj.Est: Facility-specific deviation predictions (matrix of dimension sum(Ni)*T)
    ##              lambda: Estimated region-specific facility-level variances (matrix of dimension n*M)
    ##              alpha: Estimated spatial variance parameter (vector of length L) 
    ##              sigma: Estimated measurement error variance (scalar)
    ##              nv: Estimated spatial correlation parameter (scalar)
    ############################################################################# 
    
    # Function for tuning MCMC samples 
    chain.tune <- function(chain, burnin = 2001, thin = 10){
        chain <- chain[-(1:burnin)]
        thin.index <- seq(1,length(chain), thin)
        return(chain[thin.index])
    }
    
    # Format data
    df <- data
    # Number of regions n
    nregion <- length(unique(df$rid))
    
    # Number of time points T
    ngrid <- length(unique(df$t))
    gridPoints <- unique(df$t)
    df$dt <- round(df$t * (ngrid-1)) + 1
    # Number of facilities per region Ni
    nFac <- aggregate(df$fid, by = list(df$rid), FUN = length)[,2] / ngrid
    numFac <- sum(nFac) # Total number of facilities
    
    # Extract estimated parameters from functions MST_FM_decomposition and MST_FM_MCMC
    muEst <- FPCAout$mu
    nei1 <- dim(FPCAout$psi1)[2]
    nei2 <- dim(FPCAout$psi2)[2]
    psi1.t <- FPCAout$psi1
    psi2.t <- FPCAout$psi2
    lambda.prop <- FPCAout$lambda.prop
    xi <- MCMCout$xi
    zeta <- MCMCout$zeta
    
    nbasis <- 20
    smooth_list <- smoothCon(s(gridPoints, bs = "ps", k = nbasis), data.frame(gridPoints))
    basis <- smooth_list[[1]]$X
    beta_mat <- MCMCout$beta_mat
    lambda_mat <- MCMCout$lambda
    alpha <- MCMCout$alpha
    sigma.MC <- MCMCout$sigma2
    nv <- MCMCout$nv
    
    beta0.MCs <- basis %*% t(apply(beta_mat[[1]],1,chain.tune))
    beta0Est <- rowMeans(beta0.MCs) # VCF estimates
    SD_beta <- apply(beta0.MCs, 1, sd) # sd of VCF estimates
    # 95% confidence band
    cb <- quantile(apply(abs(beta0.MCs - beta0Est)/SD_beta,2,max), .95)
    beta0.CB.lower <- beta0Est - cb*SD_beta
    beta0.CB.upper <- beta0Est + cb*SD_beta
    
    beta1.MCs <- basis %*% t(apply(beta_mat[[2]],1,chain.tune))
    beta1Est <- rowMeans(beta1.MCs) # VCF estimates
    SD_beta <- apply(beta1.MCs, 1, sd) # sd of VCF estimates
    # 95% confidence band
    cb <- quantile(apply(abs(beta1.MCs - beta1Est)/SD_beta,2,max), .95)
    beta1.CB.lower <- beta1Est - cb*SD_beta
    beta1.CB.upper <- beta1Est + cb*SD_beta
    
    beta2.MCs <- basis %*% t(apply(beta_mat[[3]],1,chain.tune))
    beta2Est <- rowMeans(beta2.MCs) # VCF estimates
    SD_beta <- apply(beta2.MCs, 1, sd) # sd of VCF estimates
    # 95% confidence band
    cb <- quantile(apply(abs(beta2.MCs - beta2Est)/SD_beta,2,max), .95)
    beta2.CB.lower <- beta2Est - cb*SD_beta
    beta2.CB.upper <- beta2Est + cb*SD_beta
    
    beta <- cbind(beta0Est, beta1Est, beta2Est)
    beta.CB <- cbind(beta0.CB.lower,beta0.CB.upper, beta1.CB.lower, beta1.CB.upper,
                     beta2.CB.lower, beta2.CB.upper)
    
    theta1.MCs <- basis %*% t(apply(beta_mat[[4]],1,chain.tune))
    theta1Est <- rowMeans(theta1.MCs) # VCF estimates
    SD_theta <- apply(theta1.MCs, 1, sd) # sd of VCF estimates
    # 95% confidence band
    cb <- quantile(apply(abs(theta1.MCs - theta1Est)/SD_theta,2,max), .95)
    theta1.CB.lower <- theta1Est - cb*SD_theta
    theta1.CB.upper <- theta1Est + cb*SD_theta
    
    theta2.MCs <- basis %*% t(apply(beta_mat[[5]],1,chain.tune))
    theta2Est <- rowMeans(theta2.MCs) # VCF estimates
    SD_theta <- apply(theta2.MCs, 1, sd) # sd of VCF estimates
    # 95% confidence band
    cb <- quantile(apply(abs(theta2.MCs - theta2Est)/SD_theta,2,max), .95)
    theta2.CB.lower <- theta2Est - cb*SD_theta
    theta2.CB.upper <- theta2Est + cb*SD_theta
    
    theta <- cbind(theta1Est, theta2Est)
    theta.CB <- cbind(theta1.CB.lower,theta1.CB.upper, theta2.CB.lower, theta2.CB.upper)
    
    LambdaEst <- colMeans(apply(lambda_mat,1,chain.tune))
    lambda <- cbind(LambdaEst, LambdaEst * lambda.prop[2])
    alphaEst <- colMeans(apply(alpha,1,chain.tune))
    sigmaEst <- mean(chain.tune(sigma.MC))
    nvEst <- mean(chain.tune(nv))
    
    # Estimated eigenscores
    xiEst <- matrix(colMeans(apply(xi,1,chain.tune)), nrow = nregion)
    zetaEst <- matrix(colMeans(apply(zeta,1,chain.tune)), nrow = numFac)
    
    # Region- and facility-specific deviation predictions
    reff.Est <- psi1.t[df$dt,1] * xiEst[df$rid,1] + psi1.t[df$dt,2] * xiEst[df$rid,2]
    feff.Est <- psi2.t[df$dt,1] * zetaEst[df$fid,1] + psi2.t[df$dt,2] * zetaEst[df$fid,2]
    # Facility-level predictions
    yEst <- beta0Est[df$dt] + df$x1 * beta1Est[df$dt] + df$x2 * beta2Est[df$dt] + 
        df$z1 * theta1Est[df$dt] + df$z2 * theta2Est[df$dt] + reff.Est + feff.Est
    
    R.traj.Est <- matrix(reff.Est, nrow = ngrid)
    R.traj.Est <- R.traj.Est[,cumsum(nFac)]
    F.traj.Est <- matrix(feff.Est, nrow = ngrid)
    
    out <- list(beta = beta, beta.CB = beta.CB, theta = theta, theta.CB = theta.CB,
                R.traj.Est = R.traj.Est, F.traj.Est = F.traj.Est, lambda = lambda, 
                alpha = alphaEst, sigma = sigmaEst, nv = nvEst)
    return(out)
}
