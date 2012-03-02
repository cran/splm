semmod <-
function (X, y, ind, tind, n, k, t, nT, w, w2, coef0 = rep(0, 2),
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

    ## set names for final parms vectors
    nam.beta <- dimnames(X)[[2]]
    nam.errcomp <- c("rho")

    ## initialize values for optimizer
    myparms0 <- coef0
    ## set bounds for optimizer
    lower.bounds <- c(-0.999)
    upper.bounds <- c(0.999)

    ## modules for likelihood
    B <- function(lambda, w) diag(1, ncol(w)) - lambda * w
    detB <- function(lambda, w) det(B(lambda, w))
    invSigma <- function(lambdapsi, n, t, w) {
        It <- diag(1, t)
        ## retrieve parms
        lambda <- lambdapsi[1]
        ## psi not used: here passing 4 parms, but works anyway
        ## because psi is last one
        ## calc inverse
        BB <- crossprod(B(lambda, w))
        invSigma <- kronecker(It, BB)
        invSigma
    }
    detSigma <- function(lambda, t, w) {
        detSigma <- t * log(detB(lambda, w))
        detSigma
    }

    ## likelihood function, both steps included
    ll.c <- function(lambdapsi, y, X, n, t, w, w2, wy) {
        ## retrieve parms
        lambda <- lambdapsi[1]
        ## calc inverse sigma
        sigma.1 <- invSigma(lambdapsi, n, t, w)
        ## do GLS step to get e, s2e
        glsres <- GLSstep(X, y, sigma.1)
        e <- glsres[["ehat"]]
        s2e <- glsres[["sigma2"]]
        ## calc ll
        due <- detSigma(lambda, t, w)
        tre <- -n * t/2 * log(s2e)
        quattro <- -1/(2 * s2e) * crossprod(e, sigma.1) %*% e
        const <- -(n * t)/2 * log(2 * pi)
        ll.c <- const + due + tre + quattro
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

    ## lag y not needed, here for compatibility
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
    nvcovpms <- length(nam.errcomp)
    covAR <- NULL
    covPRL <- covTheta[1:nvcovpms, 1:nvcovpms, drop=FALSE]

    ## final parms
    betas <- as.vector(beta[[1]])
    arcoef <- NULL
    errcomp <- myparms
    names(betas) <- nam.beta
    names(errcomp) <- nam.errcomp[which(nam.errcomp!="psi")]

    dimnames(covB) <- list(nam.beta, nam.beta)
    dimnames(covPRL) <- list(names(errcomp), names(errcomp))

    ## result
    RES <- list(betas = betas, arcoef=arcoef, errcomp = errcomp,
                covB = covB, covAR=covAR, covPRL = covPRL, ll = myll)

    return(RES)
}
