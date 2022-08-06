saremgREmod <-
function (X, y, ind, tind, n, k, t., nT, w, w2, coef0 = rep(0, 4),
    hess = FALSE, trace = trace, x.tol = 1.5e-18, rel.tol = 1e-15,
    method="nlminb", ...)
{

    ## Trieste, 22/3/2022. This is the SAR extension of semgREmod2 from the
    ## materials to the semregmod paper (Millo, J spat Econometrics, 2022)
    
    ## (summary from general version)
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

    ## this version sem2gre (orig.: 14/10/2010): general errors 
    ## a la Baltagi, Egger and Pfaffermayr

    ## set names for final parms vectors
    nam.beta <- dimnames(X)[[2]]
    nam.errcomp <- c("phi", "rho1", "rho2", "lambda")

    ## initialize values for optimizer
    myparms0 <- coef0
    
    ## set bounds for optimizer
    lower.bounds <- c(1e-08, -0.999, -0.999, -0.999)  # lag-specific line (4th parm)
    upper.bounds <- c(1e+09, 0.999, 0.999, 0.999)     # lag-specific line (idem)

    ## constraints as cA %*% theta + cB >= 0
    ## equivalent to: phi>=0, -1<=(rho1, rho2, lambda)<=1
    ## NB in maxLik() optimization cannot start at the boundary of the
    ## parameter space !
    cA <- cbind(c(1, rep(0,6)),
               c(0,1,-1,rep(0,4)),
               c(rep(0,3), 1, -1, rep(0,2)),
               c(rep(0,5), 1, -1))
    cB <- c(0, rep(1,6))

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
        ## lambda not used: here passing 4 parms, but works anyway
        ## because lambda is last one
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
        lambda <- philambda[4]                       # lag-specific line        
        ## calc inverse sigma
        sigma.1 <- invSigma(philambda, n, t., w2)
        ## lag y
        Ay <- y - lambda * wy                        # lag-specific line
        ## do GLS step to get e, s2e
        glsres <- GLSstep(X, Ay, sigma.1)         # lag-specific line (Ay for y)
        e <- glsres[["ehat"]]
        s2e <- glsres[["sigma2"]]
        ## calc ll
        zero <- t.*log(detB(lambda, w))              # lag-specific line (else zero <- 0)
        ## ex appendix to BEP, p.6, concentrated ll
        ## (notice that, from GLS step, s2e = 1/nt * crossprod(e, sigma.1) %*% e)
        ll.c <- zero - (n*t.)/2 * (log(2*pi) + 1) - (n*t.)/2 * log(s2e) -  # lag-specific
            1/2 * log(detSigma(philambda, t., w2))
                                        # or: log(detOmega(philambda, t, w, s2e)) 

        ## invert sign for minimization
        llc <- -ll.c
    }

    ## generic from here
    
    ## calc. Wy (spatial lag of y)
    ## NB problems with NAs on y, see examples_sar.R in /Articoli/semregmod
    Wy <- function(y, w, tind) {                  # lag-specific line
        wy<-list()                                # lag-specific line
        for (j in 1:length(unique(tind))) {       # lag-specific line
            yT<-y[tind==unique(tind)[j]]          # lag-specific line
             wy[[j]] <- w %*% yT                  # lag-specific line
             }                                    # lag-specific line
        return(unlist(wy))                        # lag-specific line
    }                                             # lag-specific line

    ## GLS step function
    GLSstep <- function(X, y, sigma.1) {
        b.hat <- solve(crossprod(X, sigma.1) %*% X,
                       crossprod(X, sigma.1) %*% y)
        ehat <- y - X %*% b.hat
        sigma2ehat <- (crossprod(ehat, sigma.1) %*% ehat)/(n * t.)
        return(list(betahat=b.hat, ehat=ehat, sigma2=sigma2ehat))
    }

    ## lag y once for all
    wy <- Wy(y, w, tind)                          # lag-specific line

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
    sigma.1 <- invSigma(myparms, n, t., w2)
    Ay <- y - myparms[length(myparms)] * wy       # lag-specific line
    beta <- GLSstep(X, Ay, sigma.1)               # lag-specific line

    ## final vcov(beta)
    covB <- as.numeric(beta[[3]]) *
        solve(crossprod(X, sigma.1) %*% X)

    ## final vcov(errcomp)
    nvcovpms <- length(nam.errcomp) - 1            # lag-specific line: -1
    ## error handler here for singular Hessian cases
    covTheta <- try(solve(-myHessian), silent=TRUE)
    if(inherits(covTheta, "try-error")) {
        covTheta <- matrix(NA, ncol=nvcovpms+1,    # lag-specific line: +1
                           nrow=nvcovpms+1)        # lag-specific line: +1
        warning("Hessian matrix is not invertible")
    }
    covAR <- covTheta[nvcovpms+1, nvcovpms+1, drop=FALSE]  # lag-specific line
    covPRL <- covTheta[1:nvcovpms, 1:nvcovpms, drop=FALSE]

    ## final parms
    betas <- as.vector(beta[[1]])
    arcoef <- myparms[which(nam.errcomp=="lambda")]  # lag-specific line
    errcomp <- myparms[which(nam.errcomp!="lambda")]
    names(betas) <- nam.beta
    names(arcoef) <- "lambda"                        # lag-specific line
    names(errcomp) <- nam.errcomp[which(nam.errcomp!="lambda")]

    dimnames(covB) <- list(nam.beta, nam.beta)
    dimnames(covAR) <- list(names(arcoef), names(arcoef))
    dimnames(covPRL) <- list(names(errcomp), names(errcomp))

    ## result
    RES <- list(betas = betas, arcoef=arcoef, errcomp = errcomp,
                covB = covB, covAR=covAR, covPRL = covPRL, ll = myll)

    return(RES)
}
