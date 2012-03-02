`conclikpan` <- function(lambda, env){
	
	e0e0 <- get("e0e0", envir = env)
	e1e1 <- get("e1e1", envir = env)
	e0e1 <- get("e0e1", envir = env)
	NT <- get("NT", envir = env)
	T <- get("T", envir = env)
	
	Nsig <- e0e0 - 2*lambda*e0e1 + lambda*lambda*e1e1
	sigma2 <- Nsig/NT
	
	ldet <-  do_ldet(lambda, env)
	#ret <-  T * ldet - ( (NT)/2 )*log(2*pi*sigma2) - (1/(2*sigma2))*Nsig
ret <- - (NT/2)*log(Nsig)  + T * ldet  
#	ret <-  (NT/2)*log(Nsig)  - T * ldet  
#print(T*log(det))
	  if (get("verbose", envir=env)) 
        cat("rho:\t", lambda, "\tfunction value:\t", ret, 
            "\n")
	ret
}

sarpanelerror<-function (lambda, env=env) 
{

	yt<- get("yt", envir = env)
	xt<- get("xt", envir = env)
	wyt<- get("wyt", envir = env)
	wxt<- get("wxt", envir = env)
	wy<- get("wy", envir = env)
	wx<- get("wx", envir = env)

	listw<- get("listw", envir = env)
	NT<- get("NT", envir = env)
	inde<- get("inde", envir = env)
	T<- get("T", envir = env)
	
    yco <- yt - lambda * wyt
    xco <- xt - lambda * wxt
    bb<- solve(crossprod(xco),crossprod(xco, yco) )

    ehat<- yco - xco %*% bb
    SSE <- crossprod(ehat)
  ldet <- do_ldet(lambda, env)

    ret <- T*ldet - (NT/2) * log(SSE) 

if (get("verbose", envir = env)) 
        cat("lambda:", lambda, " function:", ret, " Jacobian:", ldet, " SSE:", SSE, "\n")
 ret
}



sacsarpanel<-function (coefs, env) 
{
  	 T<-get("T", envir = env)
    SSE <- sacsarpanel_sse(coefs, env)
    n <- get("n", envir = env)
    s2 <- SSE/n
    ldet1 <- do_ldet(coefs[1], env, which = 1)
    ldet2 <- do_ldet(coefs[2], env, which = 2)

            ret <- (T * ldet1 + T * ldet2 - ((n*T/2) * log(2 * pi)) - (n*T/2) * log(s2) - 
        (1/(2 * (s2))) * SSE)
        
if (get("verbose", envir = env)) cat("rho:", coefs[1], " lambda:", coefs[2], " function:", 
            ret, " Jacobian1:", ldet1, " Jacobian2:", ldet2, 
            " SSE:", SSE, "\n")
    -ret
}


sacsarpanel_sse <- function (coefs, env) 
{
    rho <- coefs[1]
    lambda <- coefs[2]
    yl <- get("yt", envir = env) - rho * get("wyt", envir = env) - 
        lambda * get("w2yt", envir = env) + rho * lambda * get("w2wyt", 
        envir = env)
    xl <- get("xt", envir = env) - lambda * get("wxt", envir = env)
    xl.q <- qr.Q(qr(xl, LAPACK = get("LAPACK", envir = env)))
    xl.q.yl <- crossprod(xl.q, yl)
    SSE <- crossprod(yl) - crossprod(xl.q.yl)
    SSE
}


spsararlm<-function(env, zero.policy = zero.policy, con = con, llprof = llprof, tol.solve= 1e-10){

xt <- get("xt", envir = env)
yt <- get("yt", envir = env)
wyt <- get("wyt", envir = env)
w2yt <- get("w2yt", envir = env)
w2wyt <- get("w2wyt", envir = env)
wxt<-get("wxt", envir = env)

wy <- get("wy", envir = env)
wx<-get("wx", envir = env)
	
NT<-get("NT", envir = env)
T<-get("T", envir = env)
listw<-get("listw", envir = env)
inde<-get("inde", envir = env)
	
    pars <- con$pars
    lower <- con$lower
    upper <- con$upper
if (!is.null(llprof)) {
        llrho <- NULL
        lllambda <- NULL
        if (length(llprof) == 1) {
            llrho <- seq(lower[1], upper[1], length.out = llprof)
            lllambda <- seq(lower[2], upper[2], length.out = llprof)
            llprof <- as.matrix(expand.grid(llrho, lllambda))
        }
        ll_prof <- numeric(nrow(llprof))
for (i in 1:nrow(llprof)) ll_prof[i] <- sacsarpanel(llprof[i,], env = env)
    }
    if (is.null(pars)) {
        if (con$npars == 4L) {
            xseq <- c(lower[1], 0, upper[1], upper[1]) * 0.8
            yseq <- c(upper[2], 0, upper[2], lower[2]) * 0.8
            mpars <- cbind(xseq, yseq)
        }
        else {
            xseq <- seq(lower[1], upper[1], (upper[1] - lower[1])/2) * 0.8
            yseq <- seq(lower[2], upper[2], (upper[2] - lower[2])/2) * 0.8
            mpars <- as.matrix(expand.grid(xseq, yseq))
        }
    }
    
    else {
        mxs <- NULL
    }

        if (is.null(pars)) {
            mxs <- apply(mpars, 1, function(pp) -nlminb(pp, sacsarpanel, env = env, control = con$opt_control, lower = lower, upper = upper)$objective)
            pars <- mpars[which.max(mxs), ]
            optres <- nlminb(pars, sacsarpanel, env = env, control = con$opt_control, lower = lower, upper = upper)
        }
        else {
            optres <- nlminb(pars, sacsarpanel, env = env, control = con$opt_control, lower = lower, upper = upper)
        }

    LL <- -optres$objective
    if (optres$convergence != 0) 
        warning(paste("convergence failure:", optres$message))
	
    rho <- optres$par[2]
    names(rho) <- "rho"
    lambda <- optres$par[1]
    names(lambda) <- "lambda"
    
    
    lm.target <- lm(I(yt - lambda * wyt - rho * w2yt + rho * lambda * 
        w2wyt) ~ I(xt - rho * wxt) - 1)

    r <- as.vector(residuals(lm.target))
    fit <- as.vector(yt - r)
    p <- lm.target$rank
    SSE <- deviance(lm.target)
    s2 <- SSE/NT
    betas <- coefficients(lm.target)
    names(betas) <- colnames(xt)  


        coefs <- c(rho, lambda, betas)
        
        fd <- fdHess(coefs, f_sacpanel_hess, env)
        mat <- fd$Hessian
		  fdHess<- solve(-(mat), tol.solve = tol.solve)
        rownames(fdHess) <- colnames(fdHess) <- c("lambda", "lambda",colnames(xt))
        

            rho.se <- fdHess[2, 2]
            lambda.se <- fdHess[1, 1]
            rest.se <- vcov(lm.target)


return<-list(coeff=betas,lambda=lambda, rho = rho, s2=s2, asyvar1 = rest.se, lambda.se = lambda.se, rho.se = rho.se)	
	}

f_sacpanel_hess <- function (coefs, env) 
{
	T<-get("T", envir = env)
	NT<-get("NT", envir = env)
    rho <- coefs[1]
    lambda <- coefs[2]
    beta <- coefs[-(1:2)]
    SSE <- sar_sac_hess_sse_panel(rho, lambda, beta, env)
    n <- NT/T
    s2 <- SSE/n
    ldet1 <- do_ldet(rho, env, which = 1)
    ldet2 <- do_ldet(lambda, env, which = 2)
    ret <- (T * ldet1 + T * ldet2 - ((n*T/2) * log(2 * pi)) - (n*T/2) * log(s2) - 
        (1/(2 * s2)) * SSE)
    if (get("verbose", envir = env)) 
        cat("rho:", rho, "lambda:", lambda, " function:", ret, 
            " Jacobian1:", ldet1, " Jacobian2:", ldet2, " SSE:", 
            SSE, "\n")
    ret
}

sar_sac_hess_sse_panel <- function (rho, lambda, beta, env) 
{
    yl <- get("yt", envir = env) - lambda * get("wyt", envir = env) - 
        rho * get("w2yt", envir = env) + rho * lambda * get("w2wyt", 
        envir = env)
    xl <- get("xt", envir = env) - rho * get("wxt", envir = env)
    res <- yl - (xl %*% beta)
    SSE <- c(crossprod(res))
    SSE
}



sperrorlm<-function(env, zero.policy = zero.policy, interval = interval){

xt <- get("xt", envir = env)
yt <- get("yt", envir = env)
wyt <- get("wyt", envir = env)
wxt<-get("wxt", envir = env)

wy <- get("wy", envir = env)
wx<-get("wx", envir = env)

con<-get("con", envir = env)
NT<-get("NT", envir = env)
T<-get("T", envir = env)
listw<-get("listw", envir = env)
inde<-get("inde", envir = env)


opt <- optimize(sarpanelerror, interval = interval, maximum = TRUE, env = env, tol = con$tol.opt)


#opt <- nlminb(0.5,sarpanelerror,lower = interval[1], upper= interval[2], env = env)
#print(opt)

        lambda <- opt$maximum
        names(lambda) <- "rho"
        LL <- opt$objective


    lm.target <- lm(I(yt - lambda * wyt) ~ I(xt - lambda * wxt) - 
        1)
    r <- as.vector(residuals(lm.target))
    p <- lm.target$rank
    s2 <- crossprod(r)/NT
    rest.se <- (summary(lm.target)$coefficients[, 2]) * sqrt((NT - p)/NT)     
    betas <- coefficients(lm.target)
    names(betas) <- colnames(xt)  

    
        tr <- function(A) sum(diag(A))
        W <- listw2dgCMatrix(listw, zero.policy = zero.policy)
        A <- solve(Diagonal(NT/T) - lambda * W)
        WA <- W %*% A
        asyvar <- matrix(0, nrow = 2 + p, ncol = 2 + p)
        asyvar[1, 1] <- NT/(2 * (s2^2))
        asyvar[2, 1] <- asyvar[1, 2] <- T*tr(WA)/s2
        asyvar[2, 2] <- T*(tr(WA %*% WA) + tr(t(WA) %*% WA))
        asyvar[3:(p + 2), 3:(p + 2)] <- 1/as.numeric(s2) * (t(xt - lambda *wxt) %*% (xt - lambda * wxt)) 
        asyv <- solve(asyvar, tol = con$tol.solve)
        rownames(asyv) <- colnames(asyv) <- c("sigma","rho", colnames(xt))
        lambda.se <- sqrt(asyv[2, 2])
        asyvar1 <- asyv[-1,-1]
        rownames(asyvar1) <- colnames(asyvar1) <- c("rho", colnames(xt))



	return<-list(coeff=betas,lambda=lambda,s2=s2, rest.se=rest.se, lambda.se=lambda.se,asyvar1=asyvar1)
}



splaglm<-function(env, zero.policy = zero.policy, interval = interval){

xt <- get("xt", envir = env)
yt <- get("yt", envir = env)
wyt <- get("wyt", envir = env)
con<-get("con", envir = env)
NT<-get("NT", envir = env)
T<-get("T", envir = env)
listw<-get("listw", envir = env)
inde<-get("inde", envir = env)

      XpX<-crossprod(xt)
		b0<-solve(XpX,crossprod(xt,yt)) ####y on X
		b1<-solve(XpX,crossprod(xt,wyt)) ####Wy on x
		e0<-yt - xt%*% b0
		e1<-wyt - xt%*% b1
		e0e0<-crossprod(e0)
		e1e1<-crossprod(e1)
		e0e1<-t(e1)%*%e0

assign("e0e0", e0e0, envir = env)		
assign("e1e1", e1e1, envir = env)		
assign("e0e1", e0e1, envir = env)		
		
 
opt <- optimize(conclikpan,  interval = interval, maximum = TRUE, env = env, tol = con$tol.opt)


#opt <- nlminb(0.02138744, conclikpan,  lower = interval[1], upper= interval[2],  env = env)

        rho <- opt$maximum
#       rho <- opt$par
        names(rho) <- "lambda"
        LL <- opt$objective
        optres <- opt

	lm.lag <- lm((yt - rho * wyt) ~ xt - 1)
	p <- lm.lag$rank
    r <- residuals(lm.lag)
    fit <- yt - r
    names(r) <- names(fit)
	betas <- coefficients(lm.lag)
	names(betas) <- colnames(xt)

	SSE <- deviance(lm.lag)
	s2 <- SSE/NT
    
        tr <- function(A) sum(diag(A))
        W <-listw2dgCMatrix(listw, zero.policy = zero.policy)
        A <- solve(diag(NT/T) - rho * W)
        WA <- W %*% A
        one  <- T*(tr(WA %*% WA) + tr(t(WA) %*% WA))

		  lag<-function(q) trash<-unlist(tapply(q,inde,function(TT) as.matrix(WA %*% TT), simplify=TRUE))		  
		  lag2<-function(q) trash<-unlist(tapply(q,inde,function(TT) as.matrix(t(WA)%*%TT), simplify=TRUE))
		  WAxt<-apply(as.matrix(xt),2,lag)
#		  print(WAxt)
        WAWAxt<-apply(WAxt,2,lag2)
   #     print(WAWAxt)
        xtWAWAxt <- crossprod(xt,WAWAxt)
        xtWAxt <- crossprod(xt,WAxt)
        xtxt <- crossprod(xt) 
        two <- 1/as.numeric(s2) * t(betas) %*% xtWAWAxt  %*% betas
		  V <- one + two
		  zero <- rbind(rep(0, length(betas)))
        col1 <- rbind(NT/(2 * (s2^2)), T*tr(WA)/s2, t(zero))
        three <- (1/as.numeric(s2)) * xtWAxt %*% betas
        col2 <- rbind(T*tr(WA)/s2, V, three )
        col3 <- rbind(zero, t(three), 1/as.numeric(s2)* xtxt)
        asyvar <- cbind(col1, col2, col3)
        asyv <- solve(asyvar, tol = con$tol.solve)
		rownames(asyv) <- colnames(asyv) <- c("sigma","lambda", colnames(xt))
        rho.se <- sqrt(asyv[2, 2])        
        rest.se <- sqrt(diag(asyv))[-c(1:2)]
        asyvar1 <- asyv[-1,-1]
        rownames(asyvar1) <- colnames(asyvar1) <- c("lambda", colnames(xt))

    	return<-list(coeff=betas,rho=rho,s2=s2, rest.se=rest.se, rho.se=rho.se,asyvar1=asyvar1)
} 