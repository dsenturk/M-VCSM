M_VCSM_simulation <- function(numRegion=1, # number of regions (scalar, input 1 if you want 423 regions, input 2 if you want 49 regions) 
                              numFacility=1, # number of facilities per region (scalar, input 1 if you want 4-20 facilities per region, input 2 if you want 10-30 facilities per region) 
                              sigma2=0.2 # measurement error variance
){
    
    #############################################################################
    ## Description: Function for simulating one data set under the simulation design described
    ##              in Section 4.1.
    ## Args: see above
    ## Returns: list()
    #           data, data.frame with columns c("rid", "fid", "y", "t"), for further analysis using MST-FM
    #           DATA.FRAME COLUMNS (data is stored in long format): 
    #           rid: region IDs (vector of length T*sum(Ni)) T: number of time points Ni: number of facilities in region i
    #           fid: facility IDs (vector of length T*sum(Ni))
    #           y: hospitalization rate data (vector of length T*sum(Ni))
    #           t: follow-up time (vector of length T*sum(Ni)) 
    #           x1: region-level covariate (vector of T*sum(Ni))
    #           x2: region-level covariate (vector of length T*sum(Ni))
    #           z1: facility-level covariate (vector of length T*sum(Ni))
    #           z2: facility-level covariate (vector of length T*sum(Ni))
    #           data.True, data.frame with columns c("rid", "fid", "y", "t", "r.eff1","r.eff2","f.eff1","f.eff2","f.sig1","f.sig2","r.size"),
    #           for storing the true region- and facility-specific trajectories that will be used for evaluation of multilevel trajectory predictions and visualization
    #           DATA.FRAME COLUMNS: 
    #           rid: region IDs (vector of length T*sum(Ni) )
    #           fid: facility IDs (vector of length T*sum(Ni))
    #           y: hospitalization rate data (vector of length T*sum(Ni))
    #           t: follow-up time (vector of length T*sum(Ni)) 
    #           r.eff1: first component of first-level region-specific deviation U^(1)(t) (vector of length T*sum(Ni))
    #           r.eff2: second component of first-level region-specific deviation U^(1)(t)(vector of length T*sum(Ni)) 
    #           f.eff1: first component of second-level facility-specific deviation U^(2)(t) (vector of length T*sum(Ni))
    #           f.eff2: second component of second-level facility-specific deviation U^(2)(t) (vector of length T*sum(Ni)) 
    #           f.sig1: first eigenvalue of second-level facility-specific deviation U^(2)(t) (vector of length T*sum(Ni)) 
    #           f.sig2: second eigenvalue of second-level facility-specific deviation U^(2)(t) (vector of length T*sum(Ni)) 
    #           r.size: region size 1: small, 2: medium, 3: large (vector of length T*sum(Ni))
    #           xi1: first component of first-level region-specific PC score (vector of length T*sum(Ni))
    #           xi2: second component of first-level region-specific PC score (vector of length T*sum(Ni)) 
    #           zeta1: first component of second-level facility-specific PC score (vector of length T*sum(Ni))
    #           zeta2: second component of second-level facility-specific PC score (vector of length T*sum(Ni)) 
    #           x1: region-level covariate (vector of T*sum(Ni))
    #           x2: region-level covariate (vector of length T*sum(Ni))
    #           z1: facility-level covariate (vector of length T*sum(Ni))
    #           z2: facility-level covariate (vector of length T*sum(Ni))
    ################################################################################## 
    
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
    
    
    # Define the grid points used for the eigenfunctions and varying coefficient functions (VCFs)
    ngrid <- 24 # 2 year follow up
    gridPoints <- seq(0,1,length.out = ngrid)
    
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
    
    # Load adjacency matrix from maps (need to upload)
    if(numRegion==1){
        # Merged HSA map
        load("AdjMat_423.rData")
        # Adjacency matrix
        W <- AdjMat.Merge
    } else {
        # US states map 
        load("AdjMat_49.rData")
    }
    
    # Number of neighbors for each region
    D <- rowSums(W)
    
    # Spatial correlaion parameter
    nv.True <- 0.9
    
    # Covariance matrix of first-level region-specific PC scores
    covmat <- solve(diag(D) - nv.True * W)
    
    # Total number of regions
    nregion <- length(D)
    
    # Number of facilities per region
    if(numFacility==1){
        # 4-20 facilities per region
        nFac <- sample(4:20, nregion, replace = TRUE, prob = c(rep(.33/3,3), rep(.33/4,4), rep(.34/10, 10)))
        Size.region <- (nFac<7) * 1 + (nFac>=7 & nFac<11) * 2 + (nFac>=11) * 3
    }else{
        # 10-30 facilities per region
        nFac <- sample(10:30, ncounty, replace = TRUE, prob = c(rep(.33/4,4), rep(.33/4,4), rep(.34/13, 13)))
        Size.region <- (nFac<14) * 1 + (nFac>=14 & nFac<18) * 2 + (nFac>=18) * 3
    }
    
    
    # Total number of facilities
    numFac <- sum(nFac)
    
    # eigenvalues (variance of PC scores)
    # first-level (region)
    sigma11 <- 1
    sigma12 <- .25
    
    # second-level (facility)
    lambda21 <- 2
    lambda22 <- 1
    tau <- c(.15, .1, .05)
    Fac.sig <- sample(tau, nregion, replace = TRUE)
    Fac.sig1 <- Fac.sig * lambda21
    Fac.sig2 <- Fac.sig * lambda22
    
    #######################################################
    # Construct Data set
    #######################################################
    
    # Create data.frame to store dataset
    df <- data.frame(matrix(ncol = 15, nrow = numFac * ngrid))
    colnames(df) <- c("rid","fid","y","t","r.eff1","r.eff2","f.eff1","f.eff2","f.sig1","f.sig2","r.size","xi1","xi2","zeta1","zeta2")
    df$rid <- rep(1:nregion, nFac*ngrid)
    df$fid <- rep(1:numFac, each = ngrid)
    df$t <- rep(gridPoints,numFac)
    
    # First-level region-specific deviation 
    # Generate region-specific PC scores from multivariate normal distribution
    region.eff1 <- mvrnorm(1,rep(0,nregion), covmat * sigma11)
    region.eff2 <- mvrnorm(1,rep(0,nregion), covmat * sigma12)
    df$r.eff1 <- region.eff1[df$rid] * psi1.1(df$t)
    df$r.eff2 <- region.eff2[df$rid] * psi1.2(df$t)
    df$f.sig1 <- Fac.sig1[df$rid]
    df$f.sig2 <- Fac.sig2[df$rid]
    df$r.size <- Size.region[df$rid]
    df$xi1 <- region.eff1[df$rid]
    df$xi2 <- region.eff2[df$rid]
    
    # Second-level facility-specific deviation
    # Generate facility-specific PC scores from normal distribution
    fac.eff1 <- rnorm(numFac, rep(0, numFac), rep(sqrt(Fac.sig1), nFac))
    fac.eff2 <- rnorm(numFac, rep(0, numFac), rep(sqrt(Fac.sig2), nFac))
    df$f.eff1 <- fac.eff1[df$fid] * psi2.1(df$t)
    df$f.eff2 <- fac.eff2[df$fid] * psi2.2(df$t)
    df$zeta1 <- fac.eff1[df$fid]
    df$zeta2 <- fac.eff2[df$fid]
    
    # Generate region-level covariates Xi
    region.x1 <- runif(nregion, min = 0, max = 1) 
    df$x1 <- region.x1[df$rid]
    region.x2 <- rnorm(nregion, mean = 0, sd = 1)
    df$x2 <- region.x2[df$rid]
    
    # Generate facility-level covariates Zij
    fac.z1 <- rnorm(numFac, mean = 0, sd = 1)
    df$z1 <- fac.z1[df$fid] - df$t
    fac.z2 <- rnorm(numFac, mean = 0, sd = 1)
    df$z2 <- fac.z2[df$fid] + df$t
    
    # Generate outcome
    measure.err <- rnorm(dim(df)[1], 0, sqrt(sigma2))
    df$y <- beta0(df$t) + df$x1 * beta1(df$t) + df$x2 * beta2(df$t) + 
        df$z1 * theta1(df$t) + df$z2 * theta2(df$t) + df$r.eff1 + df$r.eff2 + df$f.eff1 + df$f.eff2 + measure.err
    
    out <- list(data = df[,c(1:4,16:19)], Adj.Mat = W, data.True = df)
    return(out)
}