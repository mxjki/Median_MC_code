########################## Standard SVT
SVT = function(svdu, svdd, svdv, lambda) {
    
    nzD = sum(pmax(svdd - lambda, 0) > 1e-10)
    
    if (nzD > 1) {
        svtx = svdu[, seq(nzD)] %*% (pmax(svdd - lambda, 
            0)[seq(nzD)] * t(svdv[, seq(nzD)]))
    } else {
        svtx = pmax(svdd - lambda, 0)[1] * svdu[, 1] %*% 
            t(svdv[, 1])
    }
    
    return(svtx)
    
}

KernelK = function(x) {
    
    if (x >= -1 && x <= 1) {
        Kx = -315/64 * x^6 + 735/64 * x^4 - 525/64 * 
            x^2 + 105/64
    } else {
        Kx = 0
    }
    
    return(Kx)
}



rhotau = function(A, n, m, tau) {
    
    return(sum(A * (tau - matrix(as.numeric(A <= 0), 
        nrow = n, ncol = m)))/(n * m))
    
}


min_abs = function(Yvec, omega, Bvec, tau, mu) {
    
    U1vec = pmax(0, Yvec - Bvec - tau/(sum(omega) * 
        mu))
    
    U2vec = pmax(0, Bvec - Yvec - (1 - tau)/(sum(omega) * 
        mu))
    
    return(Yvec - U1vec + U2vec)
}

update_M = function(Aobs, omega, Lmat, Umat, tau, mu) {
    
    Mmat = Lmat - Umat
    
    Bmat = Mmat[omega == 1]
    
    Mmat[omega == 1] = min_abs(Aobs[omega == 1], omega, 
        Bmat, tau, mu)
    
    return(Mmat)
}

LADMCADMM = function(Aobs, omega, p, m, tau, mu, lambda, 
    maxit, thresh, warm.start, trace.it) {
    
    Mmat = matrix(rnorm(p * m, 0, 1), p, m)
    if (!is.null(warm.start)) {
        
        if (!all(match(c("u", "d", "v"), names(warm.start), 
            0) > 0)) 
            
        stop("warm.start does not have components u, d and v")
        
        Mmat = warm.start$u %*% (warm.start$d * t(warm.start$v))
        
    }
    
    Lmat = matrix(0, p, m)
    Umat = matrix(0, p, m)
    
    error = 1
    iter = 1
    
    while (error > thresh && iter <= maxit - 1) {
        
        Mmatkm1 = Mmat
        
        oldUmat = Umat
        
        svdProxJ = svd(Mmat + Umat, LINPACK = T)
        
        Lmat = SVT(svdProxJ$u, svdProxJ$d, svdProxJ$v, 
            lambda/mu)
        
        Mmat = update_M(Aobs, omega, Lmat, Umat, tau, 
            mu)
        
        Umat = Umat + Mmat - Lmat
        
        error = max(abs(Mmat - Mmatkm1)) + max(abs(Umat - 
            oldUmat))
        
        if (trace.it) {
            
            print(error)
            
        }
        
        iter = iter + 1
        
    }
    
    
    svdGamma = svd(Mmat, LINPACK = T)
    
    ### Final cleanup of svd
    nzJ = sum(svdGamma$d > 1e-10)
    
    return(list(u = svdGamma$u[, seq(nzJ)], d = svdGamma$d[seq(nzJ)], 
        v = svdGamma$v[, seq(nzJ)], error = error, iter = iter, 
        Lmat = Lmat, Umat = Umat))
}
