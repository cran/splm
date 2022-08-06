semgREmod <-
function (X, y, ind, tind, n, k, t., nT, w, w2, coef0 = rep(0, 3),
    hess = FALSE, trace = trace, x.tol = 1.5e-18, rel.tol = 1e-15,
    method="nlminb", ...)
{

    ## Trieste, 31/8/2021. From the materials to the semregmod paper
    ## (Millo, J spat Econometrics, 2022).
    
    ## general errors a la
    ## Baltagi, Egger and Pfaffermayr

    ## This is a fix of the version from 29/09/2010, had (B'B)^(-1) instead 
    ## of just B'B in the last element of invSigma(). Works now.

    ## (summary from previous version)
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
    nam.errcomp <- c("phi", "rho1", "rho2")

    ## initialize values for optimizer
    myparms0 <- coef0
    
    ## set bounds for optimizer
    lower.bounds <- c(1e-08, -0.999, -0.999)
    upper.bounds <- c(1e+09, 0.999, 0.999)
    
    ## constraints as cA %*% theta + cB >= 0
    ## equivalent to: phi>=0, -1<=(rho1, rho2)<=1
    ## NB in maxLik() optimization cannot start at the boundary of the
    ## parameter space !
    cA <- cbind(c(1, rep(0,4)),
               c(0,1,-1,rep(0,2)),
               c(rep(0,3), 1, -1))
    cB <- c(0, rep(1,4))

    ## modules for likelihood
    B <- function(lambda, w) diag(1, ncol(w)) - lambda * w
    detB <- function(lambda, w) det(B(lambda, w))
    invSigma <- function(philambda, n, t., w) {
        Jt <- matrix(1, ncol = t., nrow = t.)
        #In <- diag(1, n)
        It <- diag(1, t.)
        Jbart <- Jt/t.
        Et <- It - Jbart
        ## retrieve parms
        phi <- philambda[1]
        rho1 <- philambda[2]
        rho2 <- philambda[3]
        ## psi not used: here passing 4 parms, but works anyway
        ## because psi is last one
        ## calc inverse
        invSigma <- kronecker(Jbart,
                              solve(t.*phi*solve(crossprod(B(rho1,w)))+
                                    solve(crossprod(B(rho2,w))))) +
            kronecker(Et, crossprod(B(rho2,w))) # fixed error:
                                        # was solve(crossprod(B...

        invSigma
    }
    detSigma <- function(philambda, t., w) {
        ## used in some formulations
        Jt <- matrix(1, ncol = t., nrow = t.)
        #In <- diag(1, n)
        It <- diag(1, t.)
        Jbart <- Jt/t.
        Et <- It - Jbart
        
        phi <- philambda[1]
        rho1 <- philambda[2]
        rho2 <- philambda[3]

        sigma <- kronecker(Jbart,
                  t.*phi*solve(crossprod(B(rho1,w))) + solve(crossprod(B(rho2,w)))) +
                 kronecker(Et, solve(crossprod(B(rho2,w))))
        detSigma <- det(sigma)
        detSigma
    }
 
    ## likelihood function, both steps included
    ll.c <- function(philambda, y, X, n, t., w, w2, wy) {
        ## retrieve parms
        phi <- philambda[1]
        rho1 <- philambda[2]
        rho2 <- philambda[3]
        ## calc inverse sigma
        sigma.1 <- invSigma(philambda, n, t., w2)
        ## do GLS step to get e, s2e
        glsres <- GLSstep(X, y, sigma.1)
        e <- glsres[["ehat"]]
        s2e <- glsres[["sigma2"]]
        ## calc ll
 
        ## ex appendix to BEP, p.6, concentrated ll
        ## (notice that, from GLS step, s2e = 1/nt * crossprod(e, sigma.1) %*% e)
        ll.c <- - (n*t.)/2 * (log(2*pi) + 1) - (n*t.)/2 * log(s2e) - 1/2 *
            log(detSigma(philambda, t., w2))

        ## invert sign for minimization
        llc <- -ll.c
    }

    ## generic from here

    ## GLS step function
    GLSstep <- function(X, y, sigma.1) {
        b.hat <- solve(crossprod(X, sigma.1) %*% X,
                       crossprod(X, sigma.1) %*% y)
        ehat <- y - X %*% b.hat
        sigma2ehat <- (crossprod(ehat, sigma.1) %*% ehat)/(n * t.)
        return(list(betahat=b.hat, ehat=ehat, sigma2=sigma2ehat))
    }

    ## lag y parm kept for compatibility
    wy <- NULL

    ## max likelihood

    ## optimization

    ## adaptive scaling
    parscale <- 1/max(myparms0, 0.1)

    if(method=="nlminb") {
    
        optimum <- nlminb(start = myparms0, objective = ll.c,
                          gradient = NULL, hessian = NULL,
                          y = y, X = X, n = n, t. = t., w = w, w2 = w2, wy = wy,
                          scale = 1, control = list(x.tol = x.tol,
                                                    rel.tol = rel.tol, trace = trace),
                          lower = lower.bounds, upper = upper.bounds)

        ## log likelihood at optimum (notice inverted sign)
        myll <- -optimum$objective
        ## retrieve optimal parms
        myparms <- optimum$par
        myHessian <- fdHess(myparms,
                            function(x) -ll.c(x, y, X, n, t., w, w2, wy))$Hessian
                                        # lag-specific line: wy
        
    } else {
        
        #stop("Optim. methods other than 'nlminb' not implemented yet")

        ## initial values are not allowed to be zero
        maxout<-function(x,a) ifelse(x>a, x, a)
        myparms0 <- maxout(myparms0, 0.01)

        ## invert sign for MAXimization
        ll.c2 <- function(phirholambda, y, X, n, t., w, w2, wy) {
            -ll.c(phirholambda, y, X, n, t., w, w2, wy)
        }

        ## max likelihood
        optimum <- maxLik(logLik = ll.c2,
                          grad = NULL, hess = NULL, start=myparms0,
                          method = method,
                          parscale = parscale,
                          constraints=list(ineqA=cA, ineqB=cB),
                          y = y, X = X, n = n, t. = t., w = w, w2 = w2, wy = wy)

        ## log likelihood at optimum (notice inverted sign)
        myll <- optimum$maximum  # this one MAXimizes
        ## retrieve optimal parms and H
        myparms <- optimum$estimate
        myHessian <- optimum$hessian
        
    }
        
        
    ## one last GLS step at optimal vcov parms
    sigma.1 <- invSigma(myparms, n, t., w)
    beta <- GLSstep(X, y, sigma.1)

    ## final vcov(beta)
    covB <- as.numeric(beta[[3]]) *
        solve(crossprod(X, sigma.1) %*% X)

    ## final vcov(errcomp)
    nvcovpms <- length(nam.errcomp)
    ## error handler here for singular Hessian cases
    covTheta <- try(solve(-myHessian), silent=TRUE)
    if(inherits(covTheta, "try-error")) {
        covTheta <- matrix(NA, ncol=nvcovpms,
                           nrow=nvcovpms)
        warning("Hessian matrix is not invertible")
    }
    covAR <- NULL
    covPRL <- covTheta

    ## final parms
    betas <- as.vector(beta[[1]])
    arcoef <- NULL
    errcomp <- myparms[which(nam.errcomp!="psi")]
    names(betas) <- nam.beta
    names(errcomp) <- nam.errcomp

    dimnames(covB) <- list(nam.beta, nam.beta)
    dimnames(covPRL) <- list(names(errcomp), names(errcomp))

    ## result
    RES <- list(betas = betas, arcoef=arcoef, errcomp = errcomp,
                covB = covB, covAR=covAR, covPRL = covPRL, ll = myll)

    return(RES)
}
