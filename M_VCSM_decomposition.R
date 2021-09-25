M_VCSM_decomposition <- function(data      # data.frame in long format with eight labeled columns (described below)
                                 # and row length equal to the length of the vectorized observations across all 
                                 # regions, facilities and time points (NOTE: all facilities must have the same set of time points)
                                 # DATA.FRAME COLUMNS: 
                                 # rid: region IDs (vector of length T*sum(Ni)) 
                                 # fid: facility IDs (vector of length T*sum(Ni))
                                 # y: hospitalization rate data (vector of length T*sum(Ni))
                                 # t: follow-up time (vector of length T*sum(Ni))
                                 # x1: region-level covariate (vector of T*sum(Ni))
                                 # x2: region-level covariate (vector of length T*sum(Ni))
                                 # z1: facility-level covariate (vector of length T*sum(Ni))
                                 # z2: facility-level covariate (vector of length T*sum(Ni))
){
    
    #############################################################################
    ## Description: Function for FPCA decomposition (Estimation Algorithm steps 1-4 in Section 2.2) of M-VCSM described in "Multilevel Varying Coefficient Spatiotemporal Model", 
    ##              including estimation of multilevel eigenfunctions and eigenvalues. 
    ## Definition:  n: number of regions, Ni: number of facilities in region i, T: number of time points, 
    ##              L: number of region-level eigencomponents, M: number of facility-level eigencomponents
    ## Args:        see above
    ## Returns:     list()
    ##              psi1: estimated first-level eigenfunctions (matrix of dimension L*T) 
    ##              psi2: estimated second-level eigenfunctions (matrix of dimension M*T) 
    ##              sigma2: estimated measurement error variance, used as an initial value in the MCMC step (scalar)
    ##              eval1: estimated region-level eigenvalues, used as initial values in the MCMC step (vector of length L)
    ##              eval2: estimated facility-level eigenvalues, used as initial values in the MCMC step (vector of length M)
    ##              lambda.prop: estimated omega, proportion of facility-level eigenvalues to the largest facility-level eigenvalue (vector of length M)
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
    
    # Format data
    df <- data
    # Number of regions n
    nregion <- length(unique(df$rid))
    
    # Number of time points T
    ngrid <- length(unique(df$t))
    gridPoints <- unique(df$t)
    
    # Number of facilities per region Ni
    nFac <- aggregate(df$fid, by = list(df$rid), FUN = length)[,2] / ngrid
    
    ###########################################################################
    # Implement the Estimation Algorithm steps 1-4 as described in Section 2.2
    ###########################################################################
    
    # Step 1: Fit a multilevel varying coefficient model to data under the working independence assumption
    b <- gam(y ~ s(t, m = 2) + s(t, m = 2, by = x1) + s(t, m = 2, by = x2) +
                 s(t, m = 2, by = z1) + s(t, m = 2, by = z2), data = df, method = "REML")
    # Obtain the residuals
    df$dt <- round(df$t * (ngrid-1)) + 1
    X.c <- b$residuals
    
    # Step 2: estimation of between and within facility raw covariance function
    # Calculate raw between facility covariance 
    Xmat <- matrix(X.c, nrow = ngrid)
    covmat.c <- matrix(0,ngrid,ngrid)
    Fac.Cov <- list()
    for(regionid in 1:nregion){
        X.c1 <- X.c[df$rid==regionid]
        covmat.c1 <- matrix(0,ngrid,ngrid)
        Xmat1 <- matrix(X.c1, nrow = ngrid)
        numFac.region <- nFac[regionid]
        for(i in 1:(numFac.region-2)){
            Xi <- Xmat1[,i]
            XXi <- rowSums(Xmat1[,(i+1):numFac.region])
            covmat.c1 <- covmat.c1 + tcrossprod(Xi,XXi) + tcrossprod(XXi,Xi)
        }
        covmat.c1 <- covmat.c1 + tcrossprod(Xmat1[,numFac.region-1], Xmat1[,numFac.region]) + tcrossprod(Xmat1[,numFac.region], Xmat1[,numFac.region-1])
        Fac.Cov[[regionid]] <- covmat.c1 / numFac.region / (numFac.region-1) # Between facility covariance for regions
        covmat.c <- covmat.c +   Fac.Cov[[regionid]] # Pool data from all regions
    }
    
    # Calculate raw within facility covariance
    G.WithinCov <- list()
    covmat.fac <- matrix(0,ngrid,ngrid)
    for(regionid in 1:nregion){
        X.c1 <- X.c[df$rid==regionid]
        covmat.c1 <- matrix(0,ngrid,ngrid)
        Xmat1 <- matrix(X.c1, nrow = ngrid)
        G.withinFac <- tcrossprod(Xmat1, Xmat1)
        numFac.region <- nFac[regionid]
        G.withinFac <- G.withinFac / numFac.region
        # subtract the between facility covariance from total covariance
        G.Fac <-  G.withinFac - Fac.Cov[[regionid]]
        G.WithinCov[[regionid]] <- G.Fac
        covmat.fac <- covmat.fac + G.Fac # Pool data from all regions
    }
    
    
    # Steps 3 and 4: Obtain estimators of between and within facility covariance functions and employ FPCA to estimate eigenfunctions and eigenvalues
    x0 <- rep(gridPoints,each = ngrid)
    x1 <- rep(gridPoints, ngrid)
    
    # Bivariate penalized spline smoothing of between facility covariance function
    cov.mat.c.s <- gam(as.vector(covmat.c / nregion) ~ te(x0, x1, k=10, bs = "ps"))
    cov.c.s <- matrix(cov.mat.c.s$fitted.values, nrow = ngrid)
    
    # FPCA on between facility covariance function
    eigen_temp <- eigen(cov.c.s, symmetric = TRUE)
    eigen_temp$values <- eigen_temp$values[which(eigen_temp$values > 0)]  # Obtain positive eigenvalues
    eigen_temp$vectors <- eigen_temp$vectors[, 1:length(eigen_temp$values)]  # Obtain eigenvectors associated with positive eigenvalues
    # eigenfunctions
    for(e in 1:length(eigen_temp$values)){  # Normalize the eigenvalues over the domain
        normal.factor <- trapz(gridPoints, eigen_temp$vectors[, e]^2)
        eigen_temp$vectors[, e] <- eigen_temp$vectors[, e] / sqrt(normal.factor)
        eigen_temp$values[e] <- eigen_temp$values[e] * normal.factor
    }
    L <- length(which(cumsum(eigen_temp$values) / sum(eigen_temp$values) < .90)) + 1 # Number of first-level eigen components
    
    # Estimated first-level eigenfunctions
    eifun1s.2 <- eigen_temp$vectors[,1]
    eifun2s.2 <- eigen_temp$vectors[,2] 
    
    # Estimate first-level eigenvalues(initial values for MCMC)
    eval1 <- eigen_temp$values[1:2]
    
    # Bivariate penalized spline smoothing of within facility covariance function
    # Store diagonal entries of the raw within facility covariance function
    olddiag <- diag(covmat.fac)
    # Remove diagonal entries
    diag(covmat.fac) <- NA
    cov.mat.fac.s <- gam(as.vector(covmat.fac) ~ te(x0, x1, k=10, bs = "ps"))
    grids2d <- data.frame(x0 = x0, x1 = x1)
    cov.fac.s <- matrix(predict(cov.mat.fac.s, newdata = grids2d), nrow = ngrid)
    cov.fac.s <- (cov.fac.s + t(cov.fac.s)) / 2 #  Symmetrize covariance function
    
    # FPCA on within facility covariance function
    eigen_temp <- eigen(cov.fac.s, symmetric = TRUE)
    eigen_temp$values <- eigen_temp$values[which(eigen_temp$values > 0)]  # Obtain positive eigenvalues
    eigen_temp$vectors <- eigen_temp$vectors[, 1:length(eigen_temp$values)]  # Obtain eigenvectors associated with positive eigenvalues
    # eigenfunctions
    for(e in 1:length(eigen_temp$values)){  # Normalize the eigenvalues over the domain
        normal.factor <- trapz(gridPoints, eigen_temp$vectors[, e]^2)
        eigen_temp$vectors[, e] <- eigen_temp$vectors[, e] / sqrt(normal.factor)
        eigen_temp$values[e] <- eigen_temp$values[e] * normal.factor
    }
    M <- length(which(cumsum(eigen_temp$values) / sum(eigen_temp$values) < .90)) + 1 # Number of second-level eigen components
    
    # Estimated second-level eigenfunctions
    eifun1s.fac <- eigen_temp$vectors[,1]
    eifun2s.fac <- eigen_temp$vectors[,2]
    
    # Estimated second-level eigenvalues(initial values for MCMC)
    eval2 <- eigen_temp$values[1:2] / nregion
    
    # Estimation proportion of the eigenvalues
    lambda.prop <- eigen_temp$values[1:M] / eigen_temp$values[1]
    
    #Estimated measurement error variance (initial values for MCMC)
    sigmaEst <- mean(olddiag - diag(cov.fac.s)) / nregion
    
    # Construct output
    psi1.t <- cbind(eifun1s.2, eifun2s.2)
    psi2.t <- cbind(eifun1s.fac, eifun2s.fac)
    out <- list(psi1 = psi1.t, psi2 = psi2.t, sigma2 = sigmaEst, lambda.prop = lambda.prop, 
                eval1 = eval1, eval2 = eval2)
    return(out)
}


