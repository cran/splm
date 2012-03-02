semsrREmod <-
function (X, y, ind, tind, n, k, t, nT, w, w2, coef0 = rep(0, 4),
    hess = FALSE, trace = trace, x.tol = 1.5e-18, rel.tol = 1e-15,
    ...)
{

    ## extensive function rewriting, Giovanni Millo 29/09/2010
    ## structure:
    ## a) specific part
    ## - set names, bounds and initial values for parms
    ## - define building blocks for likelihood and GLS as functions of parms
    ## - define likelihood
    ## b) generic part(independent from ll.c() and #parms)
    ## - fetch covariance parms from max lik
    ## - calc last GLS step
    ## - fetch betas
    ## - calc final covariances
    ## - make list of results

    # mark
    #print("uso versione 5") # fixed  'rho^2' for 'rho' in uno <- n/2 * log(1 - rho^2)

    ## set names for final parms vectors
    nam.beta <- dimnames(X)[[2]]
    nam.errcomp <- c("phi", "psi", "rho")

    ## initialize values for optimizer
    myparms0 <- coef0
    ## set bounds for optimizer
    lower.bounds <- c(1e-08, -0.999, -0.999)
    upper.bounds <- c(1e+09, 0.999, 0.999)

    ## modules for likelihood
    Vmat <- function(rho, t) {
        V1 <- matrix(ncol = t, nrow = t)
        for (i in 1:t) V1[i, ] <- rho^abs(1:t - i)
        V <- (1/(1 - rho^2)) * V1
    }
    B <- function(lambda, w) diag(1, ncol(w)) - lambda * w
    detB <- function(lambda, w) det(B(lambda, w))
    alfa2 <- function(rho) (1 + rho)/(1 - rho)
    d2 <- function(rho, t) alfa2(rho) + t - 1
    Jt <- matrix(1, ncol = t, nrow = t)
    In <- diag(1, n)
    det2 <- function(phi, rho, lambda, t) det(d2(rho, t) * (1 -
        rho)^2 * phi * In + solve(crossprod(B(lambda, w))))
    Z0 <- function(phi, rho, lambda, t) solve(d2(rho, t) * (1 -
        rho)^2 * phi * In + solve(crossprod(B(lambda, w))))
    invSigma <- function(phirholambda, n, t, w) {
        ## retrieve parms
        phi <- phirholambda[1]
        rho <- phirholambda[2]
        lambda <- phirholambda[3]
        ## psi not used: here passing 4 parms, but works anyway
        ## because psi is last one
        ## calc inverse
        invVmat <- solve(Vmat(rho, t))
        BB <- crossprod(B(lambda, w))
        invSi1 <- kronecker(invVmat, BB)
        invSi2 <- 1/(d2(rho, t) * (1 - rho)^2)
        invSi3 <- kronecker(solve(Vmat(rho, t), Jt) %*% invVmat,
            Z0(phi, rho, lambda, t) - BB)
        invSigma <- invSi1 + invSi2 * invSi3
        invSigma
    }
    ## likelihood function, both steps included
    ll.c <- function(phirholambda, y, X, n, t, w, w2, wy) {
        ## retrieve parms
        phi <- phirholambda[1]
        rho <- phirholambda[2]
        lambda <- phirholambda[3]
        ## calc inverse sigma
        sigma.1 <- invSigma(phirholambda, n, t, w)
        ## do GLS step to get e, s2e
        glsres <- GLSstep(X, y, sigma.1)
        e <- glsres[["ehat"]]
        s2e <- glsres[["sigma2"]]
        ## calc ll
        uno <- n/2 * log(1 - rho^2)
        due <- -1/2 * log(det2(phi, rho, lambda, t))
        tre <- -(n * t)/2 * log(s2e)
        quattro <- (t - 1) * log(detB(lambda, w))
        cinque <- -1/(2 * s2e) * crossprod(e, sigma.1) %*% e
        const <- -(n * t)/2 * log(2 * pi)
        ll.c <- const + uno + due + tre + quattro + cinque
        ## invert sign for minimization
        llc <- -ll.c
    }

    ## generic from here

    ## GLS step function
    GLSstep <- function(X, y, sigma.1) {
        b.hat <- solve(crossprod(X, sigma.1) %*% X,
                       crossprod(X, sigma.1) %*% y)
        ehat <- y - X %*% b.hat
        sigma2ehat <- (crossprod(ehat, sigma.1) %*% ehat)/(n * t)
        return(list(betahat=b.hat, ehat=ehat, sigma2=sigma2ehat))
    }

    ## lag y parm kept for compatibility
    wy <- NULL

    ## max likelihood
    optimum <- nlminb(start = myparms0, objective = ll.c,
                      gradient = NULL, hessian = NULL,
                      y = y, X = X, n = n, t = t, w = w, w2 = w2, wy = wy,
                      scale = 1, control = list(x.tol = x.tol,
                                 rel.tol = rel.tol, trace = trace),
                      lower = lower.bounds, upper = upper.bounds)

    ## log likelihood at optimum (notice inverted sign)
    myll <- -optimum$objective
    ## retrieve optimal parms
    myparms <- optimum$par

    ## one last GLS step at optimal vcov parms
    sigma.1 <- invSigma(myparms, n, t, w)
    beta <- GLSstep(X, y, sigma.1)

    ## final vcov(beta)
    covB <- as.numeric(beta[[3]]) *
        solve(crossprod(X, sigma.1) %*% X)

    ## final vcov(errcomp)
    covTheta <- solve(-fdHess(myparms, function(x) -ll.c(x,
        y, X, n, t, w, w2, wy))$Hessian)          # lag-specific line: wy
    covAR <- NULL
    covPRL <- covTheta

    ## final parms
    betas <- as.vector(beta[[1]])
    arcoef <- NULL
    errcomp <- myparms
    names(betas) <- nam.beta
    names(errcomp) <- nam.errcomp

    dimnames(covB) <- list(nam.beta, nam.beta)
    dimnames(covPRL) <- list(names(errcomp), names(errcomp))

    ## result
    RES <- list(betas = betas, arcoef=arcoef, errcomp = errcomp,
                covB = covB, covAR=covAR, covPRL = covPRL, ll = myll)

    return(RES)
}
