spfeml<-function(formula, data=list(), index=NULL,listw,listw2 = NULL, model=c("lag","error", "sarar"),effects=c('pooled','spfe','tpfe','sptpfe'), method="eigen",na.action=na.fail,quiet=TRUE,zero.policy = NULL, interval = NULL, tol.solve = 1e-10, control=list(), legacy = FALSE, llprof = NULL, cl = NULL, ...){

	  
        timings <- list()
        .ptime_start <- proc.time()

model<-match.arg(model)

if (model == "sarar") con <- list(tol.opt=.Machine$double.eps^0.5,fdHess=NULL, optimHess=FALSE, LAPACK = FALSE, compiled_sse=FALSE, Imult=2,cheb_q=5, MC_p=16, MC_m=30, super=FALSE, npars = 4L, lower = c(-1, -1) + .Machine$double.eps^0.5, 
        upper = c(1, 1) - .Machine$double.eps^0.5)

else con <- list(tol.opt=.Machine$double.eps^0.5,fdHess=NULL, optimHess=FALSE, compiled_sse=FALSE, Imult=2,cheb_q=5, MC_p=16, MC_m=30, super=FALSE)
	


nmsC <- names(con)
con[(namc <- names(control))] <- control
    
    if (length(noNms <- namc[!namc %in% nmsC])) 
            warning("unknown names in control: ", paste(noNms, collapse = ", "))

    if (is.null(quiet)) 
	quiet <- !get("verbose", envir = spdep:::.spdepOptions)
    stopifnot(is.logical(quiet))

	if (is.null(zero.policy))
            zero.policy <- get.ZeroPolicyOption()
        stopifnot(is.logical(zero.policy))
	  
	  
	  ## reorder data if needed
  if(!is.null(index)) {
    require(plm)
    data <- plm.data(data, index)
    }

  index <- data[,1]
  tindex <- data[,2]

	  ## record call
        ## now passed on from spml() but to be sure:
        if(is.null(cl)) cl <- match.call()

#check the model
# model<-match.arg(model)


#check the effects
effects<-match.arg(effects)


  ## check
  if(dim(data)[[1]]!=length(index)) stop("Non conformable arguments")

  ## reduce X,y to model matrix values (no NAs)
  x<-model.matrix(formula,data=data)
  y<-model.response(model.frame(formula,data=data))
  ## reduce index accordingly
  names(index)<-row.names(data)
  ind<-index[which(names(index)%in%row.names(x))]
  tind<-tindex[which(names(index)%in%row.names(x))]

  ## reorder data by cross-sections, then time
  oo<-order(tind,ind)
  x<-x[oo,]
  y<-y[oo]
  ind<-ind[oo]
  tind<-tind[oo]

  clnames<-colnames(x)
  rwnames<-rownames(x)
  #make sure that the model has no intercept if effects !=pooled
  if (effects !="pooled" && colnames(x)[1]=="(Intercept)") {
  	x<-x[,-1]
  	#cat('\n Warning: x may not contain an intercept if fixed effects are specified \n')
	}

if(is.vector(x)){
  k<-1 
  x<-matrix(x)
  colnames(x)<-clnames[-1]
  dimnames(x)[[1]]<-rwnames
  dimnames(x)[[2]]<-clnames[-1]
  }
  else   k<-dim(x)[[2]]

  


  ## det. number of groups and df
  N<-length(unique(ind))
  n<-N
  #x<-matrix(x,length(ind),k)
  ## det. max. group numerosity
  T<-max(tapply(x[,1],ind,length))
  ## det. total number of obs. (robust vs. unbalanced panels)
  NT<-length(ind)

##checks on listw
  if(is.matrix(listw)) {
    if(dim(listw)[[1]] !=N ) stop("Non conformable spatial weights")
    require(spdep)
    listw <- mat2listw(listw)
   }
  if (!inherits(listw, "listw"))
        stop("No neighbourhood list")
        
if(model == "sarar"){
	    if (is.null(listw2)) 
        listw2 <- listw
    else if (!inherits(listw2, "listw")) 
        stop("No 2nd neighbourhood list")

	}

if (is.null(con$fdHess)) con$fdHess <- method != "eigen"
        stopifnot(is.logical(con$fdHess))
	can.sim <- FALSE

if (listw$style %in% c("W", "S")) 
		can.sim <- spdep:::can.be.simmed(listw)


if(model == "sarar"){
	    if (is.null(listw2)) 
        listw2 <- listw
    else if (!inherits(listw2, "listw")) 
        stop("No 2nd neighbourhood list")

if (listw2$style %in% c("W", "S")) 
        can.sim2 <- spdep:::can.be.simmed(listw2)
	}


if(model=="sarar"){
    stopifnot(is.numeric(con$lower))
    stopifnot(length(con$lower) == 2)
    if (!is.null(con$pars)) {
        stopifnot(is.numeric(con$pars))
        stopifnot(length(con$pars) == length(con$lower))
    }
    stopifnot(is.integer(con$npars))
    stopifnot(is.numeric(con$upper))
    stopifnot(length(con$upper) == length(con$lower))
    stopifnot(is.logical(con$fdHess))
    stopifnot(is.logical(con$optimHess))
    stopifnot(is.logical(con$LAPACK))
    stopifnot(is.logical(con$super))
}

switch(model, lag = if (!quiet) cat("\n Spatial Lag Fixed Effects Model \n"),
	    error = if (!quiet) cat("\n Spatial Error Fixed Effects Model\n"),
	    sarar = if (!quiet) cat("\n Spatial SARAR Fixed Effects Model\n"),
	    stop("\nUnknown model type\n"))


  ## check whether the panel is balanced
  balanced<-N*T==NT
  if(!balanced) stop("Estimation method unavailable for unbalanced panels")


  mt<-terms(formula,data=data)
  mf<-lm(formula,data,method="model.frame")#,na.action=na.fail
  na.act<-attr(mf,'na.action')

	indic<-seq(1,T)
	inde<-as.numeric(rep(indic,each=N)) ####takes the first n observations
	indic1<-seq(1,N)
	inde1<-rep(indic1,T) ####takes observations 1,  n+1, 2n+1...
  ### generates Wy if model=='lag'
  
  
env <- new.env(parent=globalenv())
assign("y",y, envir=env)
assign("x",x, envir=env)
assign("listw",listw, envir=env)
assign("NT",NT, envir=env)
assign("T",T, envir=env)
assign("k",k, envir=env)
assign("n",n, envir=env)


wy<-unlist(tapply(y,inde, function(u) lag.listw(listw,u), simplify=TRUE))
	
#demeaning of the y and x variables depending both on model and effects

if (effects=="tpfe" | effects=="sptpfe"){
	ytms<-tapply(y,inde,mean) ####for each time period takes the mean for the cross section observations
	tpms<-function(q) tapply(q,inde,mean)
	xtms<-apply(x,2,tpms)   ###same thing for the X variable
	ytm<-rep(ytms,each=N) ##generate the NT variables
	xtm<-matrix(,NT,k)
	for (i in 1:k) xtm[,i]<-rep(xtms[,i],each=N)

if (model %in% c("lag", "sarar")) {
		wytms<-tapply(wy,inde,mean) ###same thing for Wy
		wytm<-rep(wytms,each=N)
assign("wytms",wytms, envir=env)
				
}

assign("ytms",ytms, envir=env)
assign("xtms",xtms, envir=env)
	}


if (effects=="spfe" | effects=="sptpfe"){
	ysms<-tapply(y,inde1,mean) ###for each cross-sectional unit takes the mean over the time periods
	spms<-function(q) tapply(q,inde1,mean)
	xsms<-apply(x,2,spms)
	ysm<-rep(ysms,T)
	xsm<-matrix(,NT,k)
	for (i in 1:k) xsm[,i]<-rep(xsms[,i],T)

if (model %in% c("lag", "sarar")){
			wysms<-tapply(wy,inde1,mean)
			wysm<-rep(wysms,T)
			assign("wysms",wysms, envir=env)
			}
			
assign("ysms",ysms, envir=env)
assign("xsms",xsms, envir=env)		
	}
	
	
if (effects=='pooled'){
	yt<-y  	###keep the variables with no transformation
	xt<-x
	}


if (effects=="tpfe"){ ####generate the demeaned variables for tpfe
	yt<-y-ytm
	xt<-x-xtm
						}


if(effects=="spfe"){ ####generate the demeaned variables for spfe
	yt<-y-ysm
	xt<-x-xsm

	 					}

if (effects=="sptpfe"){ ####generate the demeaned variables for both types of FE
	yt<-y - ysm - ytm + rep(mean(y),NT)
	xmm<-matrix(,NT,(k))
	for (i in 1:(k)) xmm[,i]<-rep(mean(x[,i]),NT)
	xt<-x - xsm - xtm + xmm
								}
								

	wyt<-unlist(tapply(yt,inde, function(u) lag.listw(listw,u), simplify=TRUE))

if(model=="sarar")	{
	w2yt<-unlist(tapply(yt,inde, function(u) lag.listw(listw2,u), simplify=TRUE))
	w2wyt<-unlist(tapply(wyt,inde, function(u) lag.listw(listw2,u), simplify=TRUE))
	
	}
								
								
if 	(model == "error"){
	dm<-function(A) trash<-unlist(tapply(A,inde,function(TT) lag.listw(listw,TT), simplify=TRUE))
   wxt<-apply(xt,2,dm)
   colnames(wxt)<-paste('Lag.',colnames(x), sep="")
   wx<-apply(x,2,dm)
   colnames(wx)<-paste('lag.',colnames(x), sep="")
	}

if 	(model == "sarar"){
	dm<-function(A) trash<-unlist(tapply(A,inde,function(TT) lag.listw(listw2,TT), simplify=TRUE))
   wxt<-apply(xt,2,dm)
   colnames(wxt)<-paste('Lag.',colnames(x), sep="")
   wx<-apply(x,2,dm)
   colnames(wx)<-paste('lag.',colnames(x), sep="")
	}


colnames(xt)<-dimnames(x)[[2]]

	

assign("yt",yt, envir=env)
assign("xt",xt, envir=env)
assign("wyt",wyt, envir=env)
assign("wy",wy, envir=env)


if (model %in% c("error", "sarar")){
	assign("wx",wx, envir=env)
	assign("wxt",wxt, envir=env)
if (model == "sarar")	{
	assign("w2yt",w2yt, envir=env)
	assign("w2wyt",w2wyt, envir=env)
	assign("listw2", listw2, envir = env)
	}
	} 


assign("inde",inde, envir=env)
assign("con", con, envir=env)
assign("verbose", !quiet, envir=env)
assign("can.sim", can.sim, envir=env)
assign("compiled_sse", con$compiled_sse, envir=env)
assign("similar", FALSE, envir=env)
assign("similar2", FALSE, envir = env)
assign("family", "SAR", envir = env)
assign("LAPACK", con$LAPACK, envir = env)

timings[["set_up"]] <- proc.time() - .ptime_start
.ptime_start <- proc.time()


if (!quiet) cat("Jacobian calculated using ")
	switch(method, 
		eigen = {

if (!quiet) cat("neighbourhood matrix eigenvalues\n")
			if(model != "sarar"){
                    eigen_setup(env)
                    er <- get("eig.range", envir=env)
                    if (is.null(interval)) 
                        interval <- c(er[1]+.Machine$double.eps, 
                            er[2]-.Machine$double.eps)
}
else{
        eigen_setup(env, which = 1)
        eigen_setup(env, which = 2)
	
	}

                },
	        Matrix = {
	     
	     if (listw$style %in% c("W", "S") && !can.sim) stop("Matrix method requires symmetric weights")
        if (listw$style %in% c("B", "C", "U") && !(is.symmetric.glist(listw$neighbours, 
            listw$weights))) stop("Matrix method requires symmetric weights")

if (!quiet) cat("sparse matrix Cholesky decomposition\n")

if(model=="sarar"){

        Imult <- con$Imult
        if (listw$style == "B") {
            Imult <- ceiling((2/3) * max(sapply(listw$weights, 
                sum)))
        }
        Matrix_setup(env, Imult, con$super, which = 1)
        W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
        I <- as_dsCMatrix_I(n)
        if (listw2$style %in% c("W", "S") && !can.sim2) stop("Matrix method requires symmetric weights")
        if (listw2$style %in% c("B", "C", "U") && !(is.symmetric.glist(listw2$neighbours, 
            listw2$weights))) stop("Matrix method requires symmetric weights")
        Imult <- con$Imult
        if (listw2$style == "B") {
            Imult <- ceiling((2/3) * max(sapply(listw2$weights, 
                sum)))
        }
        Matrix_setup(env, Imult, con$super, which = 2)
        W2 <- as(as_dgRMatrix_listw(listw2), "CsparseMatrix")

	
	}
else{

        Imult <- con$Imult
        if (is.null(interval)) {
            if (listw$style == "B") {
                Imult <- ceiling((2/3) * max(sapply(listw$weights, 
                  sum)))
                interval <- c(-0.5, +0.25)
            } else interval <- c(-1, 0.999)
        }
        if (is.null(con$super)) con$super <- as.logical(NA)
        Matrix_setup(env, Imult, con$super)
        W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
        I <- as_dsCMatrix_I(n)
	}

		},
	       
	       
	        spam = {
    
        if (!require(spam)) stop("spam not available")
        if (listw$style %in% c("W", "S") && !can.sim) stop("spam method requires symmetric weights")
        if (listw$style %in% c("B", "C", "U") && !(is.symmetric.glist(listw$neighbours, 
            listw$weights))) stop("spam method requires symmetric weights")
			if (!quiet) cat("sparse matrix Cholesky decomposition\n")	        	
if(model=="sarar")	{

        spam_setup(env, which = 1)
        W <- as.spam.listw(get("listw", envir = env))
        if (listw2$style %in% c("W", "S") && !can.sim2) stop("spam method requires symmetric weights")
        if (listw2$style %in% c("B", "C", "U") && !(is.symmetric.glist(listw2$neighbours, 
            listw2$weights))) stop("spam method requires symmetric weights")
        spam_setup(env, which = 2)
        W2 <- as.spam.listw(get("listw2", envir = env))
	
	}
   else{
                    spam_setup(env)
                    W <- as.spam.listw(get("listw", envir=env))
                    if (is.null(interval)) interval <- c(-1,0.999)
                    }

		},

                Chebyshev = {

        if (listw$style %in% c("W", "S") && !can.sim) stop("Chebyshev method requires symmetric weights")
        if (listw$style %in% c("B", "C", "U") && !(is.symmetric.glist(listw$neighbours, 
            listw$weights))) stop("Chebyshev method requires symmetric weights")
		    if (!quiet) cat("sparse matrix Chebyshev approximation\n")
                	
if(model=="sarar"){

        cheb_setup(env, q = con$cheb_q, which = 1)
        W <- get("W", envir = env)
        I <- as_dsCMatrix_I(n)
        if (listw2$style %in% c("W", "S") && !can.sim2) stop("Chebyshev method requires symmetric weights")
        if (listw2$style %in% c("B", "C", "U") && !(is.symmetric.glist(listw2$neighbours, 
            listw2$weights))) stop("Chebyshev method requires symmetric weights")
        cheb_setup(env, q = con$cheb_q, which = 2)
        W2 <- get("W2", envir = env)
	
	}
	
	else{                	
          cheb_setup(env, q=con$cheb_q)
          W <- get("W", envir=env)
        	 I <- as_dsCMatrix_I(n)
          if (is.null(interval)) interval <- c(-1,0.999)
                    }
                    
                },
                
                MC = {
           
          if (!listw$style %in% c("W")) stop("MC method requires row-standardised weights")
          if (!quiet) cat("sparse matrix Monte Carlo approximation\n")	
if(model=="sarar"){

        mcdet_setup(env, p = con$MC_p, m = con$MC_m, which = 1)
        W <- get("W", envir = env)
        I <- as_dsCMatrix_I(n)
        if (!listw2$style %in% c("W")) stop("MC method requires row-standardised weights")
        mcdet_setup(env, p = con$MC_p, m = con$MC_m, which = 2)
        W2 <- get("W2", envir = env)

	}
	else{                	
                    mcdet_setup(env, p=con$MC_p, m=con$MC_m)
                    W <- get("W", envir=env)
        	    I <- as_dsCMatrix_I(n)
                    if (is.null(interval)) interval <- c(-1,0.999)
                    }
                },
                
                LU = {
		    if (!quiet) cat("sparse matrix LU decomposition\n")

if(model=="sarar"){
        LU_setup(env, which = 1)
        W <- get("W", envir = env)
        I <- get("I", envir = env)
        LU_setup(env, which = 2)
        W2 <- get("W2", envir = env)

            	
            	}
            else{   
            	
            	LU_setup(env)
            W <- get("W", envir=env)
            I <- get("I", envir=env)
            if (is.null(interval)) interval <- c(-1,0.999)
            }
                },
		stop("...\nUnknown method\n"))




nm <- paste(method, "set_up", sep="_")
timings[[nm]] <- proc.time() - .ptime_start
.ptime_start <- proc.time()


	  ## uses the two estimation functions: splaglm & sperrorlm
if(model=='lag'){
    
   RES<- splaglm(env = env, zero.policy = zero.policy, interval = interval)

    res.eff<-felag(env = env, beta=RES$coeff, sige=RES$s2, effects = effects ,method =method, rho=RES$rho, legacy = legacy, zero.policy = zero.policy)    
   
    }

if (model=='error'){


  RES<- sperrorlm(env = env, zero.policy = zero.policy, interval = interval)	
    	res.eff<-feerror(env = env, beta=RES$coeff, sige=RES$s2, effects = effects ,method =method, lambda=RES$lambda, legacy = legacy)
    	
    }
    
if(model=="sarar")    {

  RES<- spsararlm(env = env, zero.policy = zero.policy, con = con, llprof = llprof, tol.solve = tol.solve)
res.eff<-felag(env = env, beta=RES$coeff, sige=RES$s2, effects = effects ,method = method, rho = RES$lambda, legacy = legacy, zero.policy = zero.policy)    	
	}

    ##calculate the R-squared
    yme <- y-mean(y)
    rsqr2 <- crossprod(yme)
    rsqr1 <- crossprod(res.eff[[1]]$res.e)
    res.R2<- 1- rsqr1/rsqr2

	#generate fixed values (from fixed_effect)
	y.hat <- res.eff[[1]]$xhat
	res <- as.numeric(res.eff[[1]]$res.e)
	N.vars<-res.eff$N.vars


	nam.rows <- dimnames(x)[[1]]
	
	
   names(y.hat) <- nam.rows
   names(res) <- nam.rows


	## make model data
   model.data <- data.frame(cbind(y,x))
   dimnames(model.data)[[1]] <- nam.rows




if (model == "lag")   spat.coef<-RES$rho
if (model == "error") spat.coef<-RES$lambda
if (model == "sarar") spat.coef <- c(RES$rho, RES$lambda)

if (is.null(RES$lambda.se) && model=="error") Coeff<-RES$coeff
else  Coeff<-c(spat.coef,RES$coeff)

type <- paste("fixed effects", model)

var<-RES$asyvar1
if(model == "sarar"){
	var<-matrix(0,(ncol(RES$asyvar1)+2),(ncol(RES$asyvar1)+2))
   var[1,1]<-	RES$lambda.se
   var[2,2]<-	RES$rho.se
   var[((2+1):ncol(var)),((2+1):ncol(var))]<-RES$asyvar1
	}

spmod <- list(coefficients=Coeff, errcomp=NULL,
                vcov = var ,spat.coef=spat.coef,
                vcov.errcomp=NULL,
                residuals=res, fitted.values=y.hat,
                sigma2=RES$s2, type=type, model=model.data,
                call=cl, logLik=RES$ll, method=method, effects=effects, 
                res.eff=res.eff)
  class(spmod) <- "splm"
  return(spmod)

	}
       
       
       
       
