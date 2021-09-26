M_VCSM_MCMC <- function(FPCAout, # output from function M_VCSM_decomposition, including the follows:
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
                        AdjMat  # Adjacency matrix from the map (0-1 matrix of dimension n*n)
){
    
    #############################################################################
    ## Description: Function for MCMC estimation (Estimation Algorithm steps 5-6) of M-VCSM model described in "Multilevel Varying Coefficient Spatiotemporal Model", 
    ##              including estimation of varying coefficient functions (VCFs), spatial variance parameters, region-specific variances, measurement error variance, 
    ##              and region- and facility-specific PC scores. 
    ## Definition:  n: number of regions, Ni: number of facilities in region i, T: number of time points, 
    ##              L: number of first-level eigencomponents, M: number of second-level eigencomponents
    ## Args:        see above
    ## Returns:     list()
    ##              beta_mat: posterior samples of varying coefficient functions (list of 5)
    ##              alpha: posterior samples of spatial variance parameter (matrix of dimension L*12000)
    ##              nv: posterior samples of spatial correlation parameter (vector of length 12000)
    ##              sigma2: posterior samples of measurement error variance (vector of length 12000)
    ##              xi: posterior samples of region-specific PC scores (matrix of dimension n*12000)
    ##              zeta: posterior samples of facility-specific PC scores (matrix of dimension sum(Ni)*12000)
    ##              lambda: posterior samples of region-specific facility-level variances (matrix of dimension n*12000)
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
    
    # Define functions
    logit <- function(x){log(x/(1-x))}
    
    unlogit <- function(x){exp(x)/(1+exp(x))}
    
    repeat.row <- function(m, times){
        return(m[rep(1:nrow(m), times = times),])
    }
    
    # Format data
    df <- data
    # Number of regions n
    nregion <- length(unique(df$rid))
    
    # Number of time points T
    ngrid <- length(unique(df$t))
    gridPoints <- unique(df$t)
    
    # Number of facilities per region Ni
    nFac <- aggregate(df$fid, by = list(df$rid), FUN = length)[,2] / ngrid
    numFac <- sum(nFac) # Total number of facilities
    
    # Region ID for each facility
    F.cid <- rep(1:nregion, nFac)
    
    # Extract estimated parameters from function MST_FM_decomposition
    nei1 <- dim(FPCAout$psi1)[2]
    nei2 <- dim(FPCAout$psi2)[2]
    psi1.t <- FPCAout$psi1
    psi2.t <- FPCAout$psi2
    spsi1 <- c(sum(psi1.t[,1]^2), sum(psi1.t[,2]^2))
    spsi2 <- c(sum(psi2.t[,1]^2), sum(psi2.t[,2]^2))
    df$dt <- round(df$t * (ngrid-1)) + 1
    lambda.prop <- FPCAout$lambda.prop
    eval1 <- FPCAout$eval1
    eval2 <- FPCAout$eval2
    sigmaEst <- FPCAout$sigma2
    
    # Basis functions for VCFs
    nbeta <- 3
    ntheta <- 2
    p <- nbeta + ntheta
    nbasis <- 20
    smooth_list <- smoothCon(s(gridPoints, bs = "ps", k = nbasis), data.frame(gridPoints))
    basis <- smooth_list[[1]]$X
    penalty <- smooth_list[[1]]$S[[1]]
    rank_penalty <- smooth_list[[1]]$rank
    ## design matrix
    # Add intercept term
    X.dat <- cbind(rep(1, dim(df)[1]),df[,c("x1","x2","z1","z2")])
    names(X.dat)[1] <- "intercept" 
    X <- kronecker(matrix(1,nrow = numFac, ncol = 1), basis)
    XX <- kronecker(matrix(1,nrow = numFac, ncol = p), basis)
    Xmat <- as.matrix(XX * X.dat[, rep(1:p, each = nbasis)])
    basis.sq <- colSums(basis^2)
    
    # Initial values for MCMC iteration
    xi=rep(0, nei1*nregion)
    zeta = rep(0, nei2*numFac)
    nv=.9
    alpha = eval1
    lambda = matrix(rep(eval2, nregion), nrow = nregion, byrow = T)
    n.sample=1 
    v.accept=.4
    qv=.01
    
    # list for saving MCMC results
    mod=list(xi=xi,zeta=zeta,nv=nv,alpha=alpha, sigma = sigmaEst, seed=.Random.seed,
             v.accept=v.accept,qv=qv,n.sample=n.sample,total=n.sample)
    
    # prior for nv
    av <- 9
    bv <- 1
    # prior for alpha
    aa <- c(2,2)
    ba <- eval1
    # prior for sigma2
    as <- 2
    bs <- sigmaEst
    # prior for tau
    a <- 1
    b <- 0.0005
    # prior for lambda
    al <- 2
    bl <- eval2[1]
    
    # Number of MCMC samples
    S_inc=12000
    n.sample=mod$n.sample+S_inc
    oldT=mod$n.sample #this should always be 1
    

    
    # Create arrays/vectors for our parameters, plug in current values into the first slot
    beta_mat <- vector("list", length = p)
    tau_sq_inv_vec <- c()
    prec_y_vec <- rep(0,n.sample)
    sigmab_sq_inv_vec <- matrix(0, nbasis, n.sample)
    lambda_mat <- matrix(0, nrow = nregion, ncol = n.sample)
    # use regression for initial values of P-spline coefficients
    df.Xmat <- as.data.frame(Xmat)
    df.Xmat$y <- df$y
    LM1 <- lm(y~.-1, data = df.Xmat)
    beta <- matrix(LM1$coefficients, nrow = nbasis)
    for(j in 1:p){
        beta_mat[[j]] <- cbind(beta_mat[[j]], beta[,j])
    }
    prec_y <- 1
    tau_sq_inv <- rep(1, p)
    sigmab_sq_inv <- rep(1, nbasis)
    xi <- array(dim=c(nei1*nregion, n.sample))
    xi[,1:mod$n.sample] <- mod$xi
    zeta <- array(dim=c(nei2*numFac, n.sample))
    zeta[,1:mod$n.sample] <- mod$zeta
    nv <- array(dim=n.sample)
    nv[1:mod$n.sample] <- mod$nv
    alpha <- array(dim = c(nei1, n.sample))
    alpha[,1:mod$n.sample] <- mod$alpha
    sigma.MC <- NULL
    sigma.MC[1:mod$n.sample] <- mod$sigma
    qv <- mod$qv
    total <- mod$total
    
    # Make the CAR precision matrix and get its determinant
    W <- AdjMat
    D <- rowSums(W)
    DvW <- diag(D) - nv[oldT] * W
    detDvW <- det(DvW) 
    
    ## Fix nv 
    # nv <- rep(0.9, n.sample)
    # skip.v <- TRUE
    
    ## Estimate nv
    skip.v <- FALSE
    
    for(it in (oldT+1):n.sample){
        
        ####################### Update VCFs ##################################
        xiMat <- matrix(xi[, it-1], nrow = nregion)
        zetaMat <- matrix(zeta[, it-1], nrow = numFac)
        RE.vec <- as.vector(t(repeat.row(xiMat %*% t(psi1.t), nFac))) + as.vector(t(zetaMat %*% t(psi2.t)))
        for(j in 1:p){
            precision <- prec_y * crossprod(X*X.dat[,j]) + tau_sq_inv[j] * penalty
            Xmat.j <- Xmat[, -((j-1)*nbasis + (1:nbasis))]
            g <- prec_y * t(X*X.dat[,j]) %*% (df$y - Xmat.j %*% as.vector(beta[,-j]) - RE.vec)
            chol_precision <- t(chol(precision)) # Get lower cholesky
            w <- forwardsolve(chol_precision, g)
            mu <- backsolve(t(chol_precision), w)
            z <- rnorm(nbasis)
            v <- backsolve(t(chol_precision), z)
            beta[,j] <- mu + v
            update_a <- a + rank_penalty / 2
            update_b <- b + t(beta[,j]) %*% penalty %*% beta[,j] / 2
            tau_sq_inv[j] <- rgamma(1, update_a, update_b)
            beta_mat[[j]] <- cbind(beta_mat[[j]], beta[,j])
        }
        
        ############# Update region-specific PC score xi's ############# 
        X.c <- df$y - Xmat %*% as.vector(beta)
        xiMat <- matrix(xi[, it-1], nrow = nregion)
        zetaMat <- matrix(zeta[, it-1], nrow = numFac)
        faceff <- rowSums(repeat.row(psi2.t, numFac) * repeat.row(zetaMat, rep(ngrid, numFac)))
        XX <- X.c - faceff
        
        for(i in 1:nei1){
            if(nei1==2){
                speff <- rep(psi1.t[,-i], numFac) * rep(xiMat[,-i], nFac * ngrid)
            } else{
                speff <- rowSums(repeat.row(psi1.t[,-i], numFac) * repeat.row(xiMat[,-i], nFac * ngrid))
            }
            Y <- (XX - speff) * psi1.t[df$dt,i]
            mu1 <- aggregate(Y, by = list(df$rid), FUN = sum)[,2]
            mu1 <- mu1/sigma.MC[it-1] + nv[it-1] * W %*% xiMat[,i] / alpha[i, it-1]
            v.xi <- alpha[i,it-1] * nFac * spsi1[i] + D * sigma.MC[it-1]
            v.xi <- alpha[i,it-1] * sigma.MC[it-1] / v.xi
            mu.xi <- v.xi * mu1
            xiMat[,i] <- rnorm(nregion, mu.xi, sqrt(v.xi))
            xi[(nregion*(i-1)+1):(nregion*i), it] <- xiMat[,i]
        }
        
        ############# Update facility-specific PC score zeta's #############
        xiMat <- matrix(xi[, it], nrow = nregion)
        zetaMat <- matrix(zeta[, it-1], nrow = numFac)
        speff <- rowSums(repeat.row(psi1.t, numFac) * repeat.row(xiMat, nFac * ngrid))
        XX <- X.c - speff
        
        for(i in 1:nei2){
            if(nei2==2){
                faceff <- rep(psi2.t[,-i], numFac) * rep(zetaMat[,-i], rep(ngrid, numFac))
            } else{
                faceff <- rowSums(repeat.row(psi2.t[,-i], numFac) * repeat.row(zetaMat[,-i], rep(ngrid, numFac)))
            }
            Y <- (XX - faceff) * psi2.t[df$dt,i]
            mu1 <- aggregate(Y, by = list(df$fid), FUN = sum)[,2]
            v.zeta <- lambda[,i] * spsi2[i] + sigma.MC[it-1]
            v.zeta <- lambda[,i] * sigma.MC[it-1] / v.zeta
            v.zeta <- rep(v.zeta, nFac)
            mu.zeta <- v.zeta * mu1 / sigma.MC[it-1]
            zeta[(numFac*(i-1)+1):(numFac*i), it] <- rnorm(numFac, mu.zeta, sqrt(v.zeta))
        }
        zetaMat <- matrix(zeta[, it], nrow = numFac)
        
        ############## Update spatial variance parameter alpha's ###########
        for(i in 1:nei1){
            bap <- ba[i] + 1/2 * t(xiMat[,i]) %*% DvW %*% xiMat[,i]
            alpha[i,it] <- 1/rgamma(1,aa[i] + nregion/2, bap)
        }
        
        ############## Update region-specific variance parameter lambda's ###########
        zeta.square <- aggregate(zetaMat^2, by = list(F.cid), FUN = sum)[,-1]
        bap <- bl + rowSums(.5 * zeta.square / matrix(lambda.prop, nrow = nregion, ncol = nei2, byrow = T))
        lambda[,1] <- 1/rgamma(nregion, al + nei2*nFac/2, bap)    
        for(i in 2:nei2){
            lambda[,i] <- lambda[,1] * lambda.prop[i]
        }
        lambda_mat[,it] <- lambda[,1]
        
        ############## Update measurement error variance sigma^2 ###########
        y <- as.vector(t(repeat.row(xiMat %*% t(psi1.t), nFac))) + as.vector(t(zetaMat %*% t(psi2.t)))
        b.post <- bs + sum((X.c-y)^2)/2
        sigma.MC[it] <- 1/rgamma(1, as + numFac*ngrid/2, b.post)
        
        
        ############## Update spatial correlation parameter nv #############
        if(skip.v==FALSE){ 
            lam <- logit(nv[it-1])
            lams <- rnorm(1,lam,qv) 
            nvs <- unlogit(lams)
            DvWs <- diag(D) - nvs*W
            detDvWs <- det(DvWs)
            r1 <- nei1/2 * log(detDvWs/detDvW) 
            
            r2 <- 0
            for(i in 1:nei1){
                r2 <- r2 + t(xiMat[,i]) %*% W %*% xiMat[,i] / alpha[i,it]
            }
            
            r2 <- (nvs - nv[it-1])/2 * r2
            
            r3 <- exp(lams-lam) * ( (1+exp(lam))/(1+exp(lams)) )^2 *  # this is the jacobian of the logit transformation
                (nvs/nv[it-1])^(av-1) * ((1-nvs)/(1-nv[it-1]) )^(bv-1)
            
            r <- exp(r1+r2)*r3
            
            accept <- ifelse(r>runif(1),1,0)
            if(accept==1){
                nv[it]=nvs
                DvW=DvWs
                detDvW=detDvWs
            }else{
                nv[it] <- nv[it-1]
            }
        }
        
    }
    mod=list(beta_mat=beta_mat, xi=xi, zeta=zeta, nv=nv, alpha=alpha, 
             lambda=lambda_mat, sigma2=sigma.MC, seed=.Random.seed)
    return(mod)
}

