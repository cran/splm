`spgm` <-
function(formula, data = list(), index = NULL, listw = NULL, listw2 = NULL, Durbin = FALSE,
         model = c("within","random"), lag = FALSE, spatial.error = TRUE,
         moments = c("initial", "weights", "fullweights"), endog = NULL, 
         instruments = NULL, lag.instruments = FALSE, verbose = FALSE, 
         method = c("w2sls", "b2sls", "g2sls", "ec2sls"), control = list(), 
         optim.method = "nlminb", pars = NULL){

## translation for uniformity
effects <- switch(match.arg(model), within = "fixed", random = "random")

# if no spatial model, there has to be at least an endogenous variable
if(is.null(listw) && is.null(endog) && !spatial.error) 
  stop("An endogenous variable should be specified")

# transofm listw into a sparse matrix
if(!is.null(listw)){
	
    if(!inherits(listw, c("listw", "Matrix", "matrix"))) stop("listw format unknown")
    if(inherits(listw, "listw"))  listw <- listw2dgCMatrix(listw)	
    if(inherits(listw, "matrix"))  listw <- Matrix(listw)

}

#check on listw2
if(!is.null(listw2) && !lag && spatial.error) stop("listw2 can be specified only with sarar")

#model sarar can have two Ws
if(lag && spatial.error){

if(is.null(listw2)) {
	
  twow <- FALSE		
	listw2 <- listw

	}
else{

    twow <- TRUE
    if(!inherits(listw2,c("listw", "Matrix", "matrix"))) stop("listw2 format unknown")
    if(inherits(listw2,"listw"))  listw2<-listw2dgCMatrix(listw2)	
    if(inherits(listw2,"matrix"))  listw2<-Matrix(listw2)
    
}	

}





if(length(model) != 1) model <- "within" 		
if((model == "within") && ((attr(terms(formula), "intercept")) == 1 )) 
  formula <- as.formula(paste(attr(terms(formula), "variables")[1+attr(terms(formula),"response")], paste(attr(terms(formula),"term.labels"), collapse="+"), sep="~"))

	


cl<-match.call()
if(!spatial.error){
	
	results<-ivsplm(formula = formula, effects = effects, 
	                data=data, index = index, endog = endog, 
	                instruments = instruments, method = method, 
	                lag = lag, listw = listw, 
	                lag.instruments = lag.instruments, Durbin = Durbin)
	
	results$nfimpacts <- "lag_gm"
	}



else{
	
	
if(!lag) results <- sperrorgm(formula = formula, data = data, index = index, 
                              listw = listw, moments = moments, endog = endog, 
                              instruments = instruments, verbose = verbose, 
                              effects = effects, control = control, 
                              lag.instruments = lag.instruments, 
                              optim.method = optim.method, pars = pars, Durbin = Durbin)
#, initial.GMerror = initial.GMerror

else {
  results <- spsarargm(formula = formula, data = data, index = index, 
                          listw = listw, listw2 = listw2,  
                          moments = moments, lag = lag, endog = endog, 
                          instruments = instruments, verbose = verbose, 
                          effects = effects, control = control, 
                          lag.instruments = lag.instruments, 
                          optim.method = optim.method, pars = pars, twow = twow, Durbin = Durbin)
  results$nfimpacts <- "sarar_gm"

}
}
#results$lag <- lag
#results$error <- error
results$call <- cl
results$ef.sph<- effects
results$legacy <- c(lag, spatial.error)
results$endog <- endog
results$est.meth <- "GM"
results$Durbin <- Durbin
class(results) <- c("splm_GM", "splm")
results

}




panel.transformations <- function(arg1, indicator, type = c("within","between","both")){

	type<-match.arg(type)
	num <- length(indicator)
	den <- length(unique(indicator))
	tmp <- num/den  
	# print(tmp)
 	if(!is.matrix(arg1)) arg1 <- as.matrix(arg1)
 	
spms <- function(q) tapply(q,indicator,mean)

 switch(type, between = {

	arg1_b<-matrix(apply(arg1,2,spms), den, ncol(arg1))
	out<-arg1_b
		
		}, 
		within ={			

	arg1_b<-matrix(apply(arg1,2,spms), den, ncol(arg1))
	 # print(arg1_b)
	arg1_bnt<- apply(arg1_b, 2, function(bb) rep(bb,  tmp))             
	# print(arg1_bnt)
	arg1_w	<- arg1 - arg1_bnt
	out<-arg1_w
		
			}, 			
			both = {
				
	arg1_b<-matrix(apply(arg1,2,spms), den, ncol(arg1))
	arg1_bnt<- apply(arg1_b, 2, function(bb) rep(bb, tmp ))  
	arg1_w	<- arg1 - arg1_bnt
	out<-list(arg1_w, arg1_b)			
				
				}		
		)	

out


	}


spgm.tsls <- function(y, yend, X, Zinst, Hinst= NULL, instr = FALSE){

	yend<-as.matrix(yend)
if (instr) H <- Hinst	
else    H <- cbind(X, Zinst)
   Z <- cbind(yend, X)
   Znames<- colnames(Z) 
	 yendp<-matrix(,nrow(yend), ncol(yend))
for(i in 1:ncol(yend))	yendp[,i] <- fitted.values(lm(yend[,i] ~ H-1  ))
    Zp <- cbind(yendp,X)
    model.fit<-lm(y~Zp-1)
    biv <- coefficients(model.fit)
readout<- which(is.na(biv))
if(any(is.na(biv)))  yp <- as.matrix(Z[,-which(is.na(biv))]) %*% biv[-which(is.na(biv))]
else yp <- Z %*% biv



        ehat <- y-yp
        sse <- c(crossprod(ehat, ehat))
        df <- model.fit$df.residual
        s2 <- sse/df
        
if(any(is.na(biv)))   Zp<-as.matrix(Zp)[,-which(is.na(biv))]  


ZpZi<-solve(crossprod(Zp))   

         varb<-ZpZi *  s2



names(biv)<-Znames
#print(length(biv))
#print(Znames)
#print(dim(varb))
colnames(varb) <- rownames(varb) <- Znames[-which(is.na(biv))] 

if(any(is.na(biv))) biv<-biv[-which(is.na(biv))] 
else	biv <- biv

#print(length(biv))



        result <- list(coefficients = biv, vcov = varb, sse = sse, 
            residuals = as.numeric(ehat), df = df, Zp = Zp, readout = readout)
    
    result
}



sperrorgm <- function(formula, data = list(), 
                    index = NULL, listw = NULL, 
                    moments = c("initial","weights","fullweights"), 
                    endog = NULL, instruments = NULL, verbose = FALSE, 
                    effects = c("fixed","random"), control = list(), 
                    lag.instruments = lag.instruments, 
                    optim.method = optim.method, pars = pars, Durbin){

effects<-match.arg(effects)
moments<-match.arg(moments)
indes<-index

 #if(!is.null(index)) {
    #require(plm)
 #   data <- plm.data(data, index)
 #   }

 ## substituted module based on plm.data (deprecated)
 ## with one based on pdata.frame
  if(!("pdata.frame" %in% class(data))) {
    data <- pdata.frame(data, index)
  }   

  index <- data[,1]
  tindex <- data[,2]

  names(index)<-row.names(data)
  ind <-index[which(names(index)%in%row.names(data))]
  tind<-tindex[which(names(index)%in%row.names(data))]
  spord <- order(tind, ind)
  data <-  data[spord,]


  cl <- match.call()
  if(dim(data)[[1]]!=length(index)) stop("Non conformable arguments")


    mt <- terms(formula, data = data)
    mf <- lm(formula, data, na.action = na.fail, method = "model.frame")

    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
    N <-length(unique(ind))
    k <-dim(x)[[2]]
    t <-max(tapply(x[,1],ind,length))
    NT<-length(ind)

    
  
  balanced<-N*t==NT
if(!balanced) stop("Estimation method unavailable for unbalanced panels")

I_T <- Diagonal(t)
Ws <- kronecker(I_T, listw)


if(!is.null(endog)){
	endog <- as.matrix(lm(endog, data, na.action = na.fail, method = "model.frame"))
if(is.null(instruments)) stop("No instruments specified  for the additional endogenous variable")
else instruments <- as.matrix(lm(instruments, data, na.action = na.fail, method = "model.frame"))	
	}

if(lag.instruments){
  winst <- Ws %*% instruments
  wwinst <- Ws %*% winst
  instruments <- cbind(instruments, winst, wwinst)
}

if(inherits(Durbin, "formula")) xdur <- as.matrix(lm(Durbin, data, na.action = na.fail, method="model.frame"))
else xdur <- NULL


indic <- rep(1:N,t)

#generate an X for w2sls 
#X <- x



switch(effects, 

	fixed = {

               if(is.null(endog)){
                        if(isTRUE(Durbin)  | inherits(Durbin, "formula")){
    if(inherits(Durbin, "formula")){
      
      xdur <- as.matrix(lm(Durbin, data, na.action = na.fail, method="model.frame"))
      colnmx <- colnames(x)
      colnameswx <- paste("lag_", colnames(xdur), sep="")
      wx <- Ws %*% xdur
      x <- as.matrix(cbind(x, wx))
      colnames(x) <- c(colnmx, colnameswx)
     
    }
    else{
      
      colnmx <- colnames(x)
      
      if(colnmx[1] == "(Intercept)"){
        
        wx <- Ws %*% x[,-1]
        colnameswx <- paste("lag_", colnames(x)[-1], sep = "")
        
      }
      else{
        
        wx <- Ws %*% x
        colnameswx <- paste("lag_", colnames(x), sep = "")
        
      }
      
      x <- as.matrix(cbind(x, wx))
      colnames(x) <- c(colnmx, colnameswx)
      
    }
  }   
                                Xwithin  <- panel.transformations(x, indic, type= "within")
                                Xwithin  <- as.matrix(Xwithin)
                                del      <- which(diag(var(as.matrix(Xwithin)))==0)
                                ywithin  <- panel.transformations(y, indic, type= "within")	
                                result   <- lm(ywithin ~ Xwithin[,-del] -1)
                                
}	

                  else 	{
                    result <- ivplm.w2sls(Y = y, X = x, H = instruments, endog = endog, 
                                              twow = FALSE, lag = FALSE, listw = Ws,  listw2 = NULL,
                                              lag.instruments = lag.instruments,  t = t, N = N, NT = NT, 
                                              Durbin = Durbin, xdur = xdur)
                    
                    del <- result$del
                    x   <- result$xdu
                  
                  }
	  
	  res <- as.matrix(residuals(result))

Gg<-fswithin(Ws, res, N, t)

if(is.null(pars)) {
	
    wres <- as.matrix(Ws %*% res)
    r.init <- solve(crossprod(res),crossprod(res,wres))
if(is.null(endog))	v.init <- crossprod(res)/NT	
else    	        v.init <- result$sigmav
	pars <- c(r.init, v.init)	
}

if (optim.method == "nlminb") estim1 <- nlminb(pars, arg, v = Gg,
                                               verbose = verbose, control = control, 
                                               lower=c(-0.999,0), upper=c(0.999, Inf))

else estim1 <- optim(pars, arg, v = Gg, verbose = verbose, control = control, 
                     method = optim.method)

	 finrho=estim1$par[1]
	 finsigmaV=estim1$par[2]

   wy <- as.matrix(Ws %*% y)
   yt <- y-finrho*wy
   xl<- as.matrix(Ws %*%  x)
   #print(head(xl))
   xt <- x-finrho*xl

#print(head(xt))
	yf<-panel.transformations(yt,indic, type= "within")
  xf<-panel.transformations(xt,indic, type= "within")
	xf<-xf[,-del]
	xf<-as.matrix(xf)
	colnames(xf) <- colnames(x)[-del]
	wxf <- as.matrix(Ws %*% xf)

if (is.null(endog)){

result <- lm(as.matrix(yf)~as.matrix(xf)-1)
betaGLS <- coefficients(result)[which(!is.na(coefficients(result)))]
vcov <- vcov(result)[which(!is.na(coefficients(result))), which(!is.na(coefficients(result)))]


	names(betaGLS) <- colnames(xf)[which(!is.na(coefficients(result)))]
  errcomp<-rbind(finrho,finsigmaV)
  nam.errcomp <- c("rho","sigma^2_v")
  rownames(errcomp) <- nam.errcomp
  colnames(errcomp)<-"Estimate"
  model.data <- data.frame(cbind(y,x[,-1]))

  type <- "Spatial fixed effects error model (GM estimation)"
    spmod <- list(coefficients= betaGLS, errcomp=errcomp,
                vcov=vcov, vcov.errcomp=NULL,
                residuals=residuals(result), fitted.values=(y-as.vector(residuals(result))),
                sigma2=crossprod(residuals(result))/result$df.residual, 
                type=type, rho=errcomp, model=model.data, logLik=NULL)
	
	}

else{
	
   endogl <- as.matrix(Ws %*% endog)
   endogt <- endog - finrho* endogl
   endogf <- panel.transformations(endogt,indic, type= "within")

  instwithin <- result$Hwithin
  instwithin<-cbind(xf, instwithin)
  
  model.fit <- spgm.tsls(yf, endogf, xf, Hinst = instwithin, instr  = T )  
  betaGLS <- model.fit$coefficients 


  Zp <- model.fit$Zp
  Z  <- cbind(endog,x[,-del])
  fv <- as.matrix(Z)[, qr(Z)$pivot[seq_len(qr(Z)$rank)]] %*% as.matrix(betaGLS)
  egls <- y - fv
  SGLS <- sum(egls^2)/(N-1)
  xfpxfNT <- (1/NT)*crossprod(Zp)/finsigmaV
  PSI <- solve(xfpxfNT)
  covbeta <- PSI/NT

  nam.beta <- names(model.fit$coefficients)
  names(betaGLS) <- nam.beta
  errcomp<-rbind(finrho,finsigmaV)
  nam.errcomp <- c("rho","sigma^2_v")
  rownames(errcomp) <- nam.errcomp
  colnames(errcomp)<-"Estimate"
   model.data <- data.frame(cbind(y, as.matrix(x[,-1])))

  type <- "Spatial fixed effects error model with additional endogenous variables (GM estimation)"
    spmod <- list(coefficients=betaGLS, errcomp=  errcomp,
                vcov=covbeta, vcov.errcomp=NULL,
                residuals=as.vector(egls), fitted.values=fv,
                sigma2=SGLS, type=type, rho=errcomp, model=model.data, sigmav = errcomp[2] , logLik=NULL)
  
	}

		},
		
	random = {
		
if(is.null(endog))	{
  
  if(isTRUE(Durbin)  | inherits(Durbin, "formula")){
    if(inherits(Durbin, "formula")){
      
      xdur <- as.matrix(lm(Durbin, data, na.action = na.fail, method="model.frame"))
      colnmx <- colnames(x)
      colnameswx <- paste("lag_", colnames(xdur), sep="")
      wx <- Ws %*% xdur
      x <- as.matrix(cbind(x, wx))
      colnames(x) <- c(colnmx, colnameswx)
      
    }
    else{
      
      colnmx <- colnames(x)
      
      if(colnmx[1] == "(Intercept)"){
        
        wx <- Ws %*% x[,-1]
        colnameswx <- paste("lag_", colnames(x)[-1], sep = "")
        
      }
      else{
        
        wx <- Ws %*% x
        colnameswx <- paste("lag_", colnames(x), sep = "")
        
      }
      
      x <- as.matrix(cbind(x, wx))
      colnames(x) <- c(colnmx, colnameswx)
      
    }
  }   
  
result<-lm(y~x-1) 
#print(coefficients((result)))
res<-as.matrix(residuals(result))
Gg<-fs(Ws,res,N,t)

## parameter initial values 
 if(is.null(pars)) {
    wres <- as.matrix(Ws %*% res)
    r.init <- solve(crossprod(res),crossprod(res,wres))
	v.init <- crossprod(res)/NT	
	pars <- c(r.init, v.init)	
}


 if (optim.method == "nlminb") estim1 <- nlminb(pars, arg, v = Gg, 
                                                verbose = verbose, 
                                                control = control, 
                                                lower=c(-0.999,0), upper=c(0.999,Inf))
else estim1 <- optim(pars, arg, v = Gg, verbose = verbose, 
                     control = control, method = optim.method)

urub<-res- estim1$par[1]*Gg$ub
Q1urQ1ub<-Gg$Q1u - estim1$par[1]*Gg$Q1ub
S1 <- crossprod(urub, Q1urQ1ub)/N

switch(moments, 
	  
	  initial = {
	  	
	finrho=estim1$par[1]
	finsigmaV=estim1$par[2]
	finsigma1=S1
	
    }, 
    
    weights = {
    	
  	   Ggw<-pw(bigG=Gg$bigG, smallg=Gg$smallg, Q1u=Gg$Q1u,Q1ub=Gg$Q1ub,Q1ubb=Gg$Q1ubb, u=res, ub=Gg$ub,ubb=Gg$ubb,N=N, TR=Gg$TR)
      pars2<-c(estim1$par[1],estim1$par[2],S1)

 if (optim.method == "nlminb") estim2 <- nlminb(pars2, arg1, v = Ggw,t=t,ss=estim1$par[2] ,SS=S1, verbose = verbose, control = control, lower=c(-0.999,0,0), upper=c(0.999,Inf,Inf))
 else      estim2 <- optim(pars2, arg1, v = Ggw,t=t,ss=estim1$par[2] ,SS=S1, verbose = verbose, control = control, method = optim.method)

	finrho=estim2$par[1]
	finsigmaV=estim2$par[2]
	finsigma1=estim2$par[3]

    }, 
    
    fullweights = {

	   Ggw<-pw(bigG=Gg$bigG, smallg=Gg$smallg, Q1u=Gg$Q1u,Q1ub=Gg$Q1ub,Q1ubb=Gg$Q1ubb, u=res, ub=Gg$ub,ubb=Gg$ubb,N=N, TR=Gg$TR)
      weights<-tw(listw, N)
      pars2<-c(estim1$par[1],estim1$par[2],S1)

      if (optim.method == "nlminb") estim3 <-nlminb(pars2, arg2, v = Ggw, t=t, 
                                               ss=estim1$par[2] ,SS=S1, TW=weights$TW, 
                                               verbose = verbose, control = control, 
                                               lower=c(-0.999,0,0), upper=c(0.999,Inf,Inf))

      else   estim3 <-optim(pars2, arg2, v = Ggw, t=t, 
                      ss=estim1$par[2] ,SS=S1, TW=weights$TW, 
                      verbose = verbose, control = control, 
                      method = optim.method)
      
	finrho=estim3$par[1]
	finsigmaV=estim3$par[2]
	finsigma1=estim3$par[3]


	
		}, 
		stop("...\nUnknown method\n"))
	
	
	}
else{
 
result1<-ivplm.w2sls(Y = y,X =x, H = instruments, endog = endog, twow = FALSE, 
                     lag = FALSE, listw = Ws,  listw2 = NULL, lag.instruments = lag.instruments, 
                     t, N, NT, Durbin = Durbin, xdur = xdur)
result2<-ivplm.b2sls(Y = y,X =x, H = instruments, endog = endog,  twow = FALSE,
                     lag = FALSE, listw = Ws, listw2 = NULL, lag.instruments = lag.instruments,
                     t, N, NT, Durbin = Durbin, xdur = xdur)



xw   <- result1$xdu[, qr(result1$xdu)$pivot[seq_len(qr(result1$xdu)$rank)]]
xb   <- result2$xdu[, qr(result2$xdu)$pivot[seq_len(qr(result2$xdu)$rank)]]
x <- cbind(xw, xb)
x <- x[,unique(colnames(x))]

res1<-as.matrix(as.numeric(residuals(result1)))
res2<-as.matrix(as.numeric(residuals(result2)))

Gg<-fswithin(Ws,res1,N,t)

if(is.null(pars)) {
    wres <- as.matrix(Ws %*% res1)
    r.init <- solve(crossprod(res1),crossprod(res1,wres))
	v.init <- result1$sigmav	
	pars <- c(r.init, v.init)	
}


if (optim.method == "nlminb")  estim1 <- nlminb(pars, arg, v = Gg, 
                                                verbose = verbose, 
                                                control = control, 
                                                lower=c(-0.999, 0), 
                                                upper=c(0.999,Inf))

else estim1 <- optim(pars, arg, v = Gg, 
                     verbose = verbose, 
                     control = control, 
                     method = optim.method)


Wres2 <- as.matrix(listw %*% res2)
urhoWu<-res2 - estim1$par[1] * Wres2
finsigma1<-crossprod(urhoWu)/N


switch(moments, 
	  
	  initial = {
	  	
	finrho=estim1$par[1]
	finsigmaV=estim1$par[2]
	finsigma1=finsigma1
	
    }, 
    
    weights = {
    	
    	Ggw<-pwbetween(bigG=Gg$bigG, smallg=Gg$smallg, 
    	               u=res2, N=N, t=t, TR=Gg$TR, listw = listw)
      pars2<-c(estim1$par[1],estim1$par[2],finsigma1)

 if (optim.method == "nlminb")  estim2 <- nlminb(pars2, arg1, v = Ggw, t=t, 
                                                 ss= estim1$par[2], SS=finsigma1, 
                                                 verbose = verbose, 
                                                 control = control, 
                                                 lower=c(-0.999,0,0), 
                                                 upper=c(0.999,Inf,Inf))

else    estim2 <- optim(pars2, arg1, v = Ggw, t=t, 
                        ss= estim1$par[2], SS=finsigma1 , 
                        verbose = verbose, control = control, 
                        method = optim.method)
      
	finrho=estim2$par[1]
	finsigmaV=estim2$par[2]
	finsigma1=estim2$par[3]

    }, 
    
    fullweights = {

	   Ggw<-pwbetween(bigG = Gg$bigG, smallg = Gg$smallg, 
	                  u = res2, N = N, t = t, 
	                  TR = Gg$TR, listw = listw)
      weights <- tw(listw, N)
      pars2 <- c(estim1$par[1],estim1$par[2],finsigma1)

 if (optim.method == "nlminb") estim3 <-nlminb(pars2, arg2, v = Ggw, 
                                               t = t, ss = estim1$par[2], 
                                               SS=finsigma1, TW = weights$TW, 
                                               verbose = verbose, 
                                               control = control, 
                                               lower=c(-0.999,0,0), 
                                               upper=c(0.999,Inf,Inf))
 
else    estim3 <-optim(pars2, arg2, v = Ggw, t = t, 
                       ss = estim1$par[2], SS = finsigma1, 
                       TW = weights$TW, verbose = verbose, 
                       control = control, method = optim.method)

	finrho=estim3$par[1]
	finsigmaV=estim3$par[2]
	finsigma1=estim3$par[3]


	
		}, 
		stop("...\nUnknown method\n"))
	}

theta<- 1-(sqrt(finsigmaV)/sqrt(finsigma1))	
wy <- as.matrix(Ws %*% y)
yt <- y-finrho*wy
xl<- as.matrix(Ws %*% x)
xt <- x-finrho*xl
#[, qr(Z)$pivot[seq_len(qr(Z)$rank)]]
#if(is.null(endog)) xt <- xt[,which(!is.na(coefficients(result)))]
xt <- xt[, qr(xt)$pivot[seq_len(qr(xt)$rank)]]

ytmt<-tapply(yt, indic, mean)
ytNT<-rep(ytmt, t)
yf<-(yt - as.numeric(theta)*ytNT)

dm1<- function(A) rep(unlist(tapply(A, indic, mean, simplify=TRUE)), t)
xtNT<-apply(xt,2,dm1)
xf<- as.matrix(xt - as.numeric(theta)*xtNT)
wxf <- as.matrix(Ws %*% xf)


if (is.null(endog)){
	
  xfpxf<-crossprod(xf)
  xfpxfi<-solve(xfpxf)
  betaGLS<-xfpxfi%*%crossprod(xf,yf) 
  betaGLS<- as.numeric(betaGLS)

  fv<-as.vector(x[, qr(xt)$pivot[seq_len(qr(xt)$rank)]] %*% betaGLS)
  egls<-y - fv
  SGLS<-sum(egls^2)/(N-1)
  xfpxfNT<-(1/NT)*xfpxf/finsigmaV
  PSI<-solve(xfpxfNT)
  covbeta<-PSI/NT

  errcomp<-rbind(finrho,finsigmaV,finsigma1,theta)
  nam.beta <- names(betaGLS)
  nam.errcomp <- c("rho","sigma^2_v",'sigma^2_1',"theta")
  names(betaGLS) <- nam.beta
  rownames(errcomp) <- nam.errcomp
  colnames(errcomp)<-"Estimate"
model.data <- data.frame(cbind(y,as.matrix(x)))
sigma2 <- SGLS
  type <- "Spatial random effects error model (GM estimation)"
    spmod <- list(coefficients=betaGLS, errcomp=errcomp,
                vcov=covbeta, vcov.errcomp=NULL,
                residuals=as.vector(egls), fitted.values=fv,
                sigma2=sigma2,type=type, rho=errcomp, model=model.data,
                call=cl, logLik=NULL, coy=yt, cox=xt, rhs=k)
  

	}

else{
	
	endogl<- as.matrix(Ws %*% endog)
   endogt <- endog - finrho*endogl

endogt<-as.matrix(endogt)
endogtNT<-apply(endogt,2,dm1)
endogf<-(endogt - as.numeric(theta)*endogtNT)

# instt<-panel.transformations(instruments,indic, type= "both")
# instbetween<-instt[[2]]
# instwithin<-instt[[1]]
instbetween <- result2$Hbetween
instwithin <- result1$Hwithin

instbetweennt<-matrix(,NT,ncol(instbetween))
for (i in 1:ncol(instbetween)) instbetweennt[,i]<-rep(instbetween[,i], t)

 Zf<-cbind(endog,x[, qr(x)$pivot[seq_len(qr(x)$rank)]])
 # Hins<-cbind(xf,wxf, instwithin,instbetweennt)
 Hins<-cbind(xf, instwithin,instbetweennt)
 Hins <- as.matrix(Hins[, qr(Hins)$pivot[seq_len(qr(Hins)$rank)]])
 
  model.fit <- spgm.tsls(yf,as.matrix(endogf),xf, Hinst = Hins, instr  = T)  
  betaGLS <- model.fit$coefficients 



  Zp<-model.fit$Zp
  fv<-as.matrix(Zf) %*% as.matrix(betaGLS)
  egls<-y - fv
  SGLS<-sum(egls^2)/(N-1)
  xfpxfNT<-(1/NT)*crossprod(Zp)/finsigmaV
  PSI<-solve(xfpxfNT)
  covbeta<-PSI/NT

  errcomp<-rbind(finrho,finsigmaV,finsigma1,theta)
  nam.beta <- colnames(Zf)
  names(betaGLS) <- nam.beta
  nam.errcomp <- c("rho","sigma^2_v",'sigma^2_1',"theta")
  rownames(errcomp) <- nam.errcomp
  colnames(errcomp)<-"Estimate"
model.data <- data.frame(cbind(y,as.matrix(x)))
sigma2 <- SGLS
  type <- "Spatial random effects error model with additional endogenous variables (GM estimation)"
    spmod <- list(coefficients=betaGLS, errcomp=errcomp,
                vcov=covbeta, vcov.errcomp=NULL,
                residuals=as.vector(egls), fitted.values=fv,
                sigma2=sigma2,type=type, rho=errcomp, model=model.data,
                call=cl, logLik=NULL, coy=yt, cox=xt, rhs=k)

}

	
    }, 

    
    stop("...\nUnknown method\n"))
	
	
return(spmod)	
	
	}




spsarargm<-function(formula, data = list(), 
                    index = NULL, listw, listw2 = NULL, 
                    moments = c("initial", "weights", "fullweights"), 
                    lag= FALSE, endog = NULL, instruments = NULL, 
                    verbose = FALSE, effects = c("fixed","random"), 
                    control = list(), lag.instruments = lag.instruments, 
                    optim.method = optim.method, pars = pars, twow, Durbin){


effects<-match.arg(effects)
moments<-match.arg(moments)
indes<-index

 #if(!is.null(index)) {
    #require(plm)
 #   data <- plm.data(data, index)
 #   }

 ## substituted module based on plm.data (deprecated)
 ## with one based on pdata.frame
  if(!("pdata.frame" %in% class(data))) {
    data <- pdata.frame(data, index)
  } 
  
  index <- data[,1]
  tindex <- data[,2]

  names(index)<-row.names(data)
  ind <-index[which(names(index)%in%row.names(data))]
  tind<-tindex[which(names(index)%in%row.names(data))]
   spord <- order(tind, ind)
   data <-  data[spord,]


  cl <- match.call()
  if(dim(data)[[1]]!=length(index)) stop("Non conformable arguments")


    mt <- terms(formula, data = data)
    mf <- lm(formula, data, na.action = na.fail, method = "model.frame")

    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
    namesx <- colnames(x)
    N<-length(unique(ind))
    k<-dim(x)[[2]]
    t <-max(tapply(x[,1],ind,length))
    NT <-length(ind)
    indic <- rep(1:N,t)

  balanced<-N*t == NT
if(!balanced) stop("Estimation method unavailable for unbalanced panels")


I_T <- Diagonal(t)
Ws <- kronecker(I_T, listw)

if(twow) Ws2 <- kronecker(I_T, listw2)
else Ws2 <-Ws

listw2nn <- Ws2[1:N, 1:N]

if(!is.null(endog)){
	endog <- as.matrix(lm(endog, data, na.action = na.fail, method = "model.frame"))
if(is.null(instruments)) stop("No instruments specified  for the additional endogenous variables")
else instruments <- as.matrix(lm(instruments, data, na.action = na.fail, method = "model.frame"))	
	}

if(lag.instruments){
  winst <- Ws %*% instruments
  wwinst <- Ws %*% winst
  instruments <- cbind(instruments, winst, wwinst)
}

if(inherits(Durbin, "formula")) xdur <- as.matrix(lm(Durbin, data, na.action = na.fail, method="model.frame"))
else xdur <- NULL

##transform X	
# transx<-panel.transformations(x, indic, type= "both")
# Xbetween<-transx[[2]]
# Xwithin<-transx[[1]]
# Xbetween <- as.matrix(Xbetween)
# Xwithin <- as.matrix(Xwithin)
# del<-which(diag(var(as.matrix(Xwithin)))==0)
# if (namesx[1] == "(Intercept)") Xbetween <- Xbetween[,-1]
# delb<-which(diag(var(as.matrix(Xbetween)))==0)
# if(length(delb)==0) Xbetween<-Xbetween
# else Xbetween<-Xbetween[,-delb]
# Xbetween<-cbind(1,Xbetween)
# deltot<-c(del,delb)
# Xbetweennt<-matrix(,NT, ncol(Xbetween))
# for (i in 1:ncol(Xbetween)) Xbetweennt[,i]<-rep(Xbetween[,i], t)

switch(effects, 

	fixed = {
		
if(is.null(endog)) result<-ivplm.w2sls(Y = y,X = x, 
                                       lag = TRUE, listw = Ws, 
                                       listw2 = Ws2, twow = twow, 
                                       lag.instruments = lag.instruments, 
                                       t = t, N = N, NT = NT, Durbin = Durbin, xdur = xdur)

else 	result<-ivplm.w2sls(Y = y, X = x, 
                          H = instruments, endog = endog, 
                          lag = TRUE, listw = Ws, 
                          listw2 = Ws2, twow = twow, 
                          lag.instruments = lag.instruments, 
                          t = t, N = N, NT = NT, Durbin = Durbin, xdur =xdur)

del <- result$del
x   <- result$xdu

# print(result$coefficients)

res <- as.matrix(as.numeric(residuals(result)))
Hwithin <- cbind(result$Xwithin, result$Hwithin)


Gg<-fswithin(Ws2,res,N,t)

if(is.null(pars)) {
	
    wres <- as.matrix(Ws2 %*% res)
    r.init <- solve(crossprod(res),crossprod(res,wres))
	v.init <- crossprod(res)/NT	
	pars <- c(r.init, v.init)	
}




if (optim.method == "nlminb")  estim1 <- nlminb(pars, arg, v = Gg, verbose = verbose,  lower=c(-0.999,0), upper=c(0.9999,Inf))

else estim1 <- optim(pars, arg, v = Gg, verbose = verbose, control = control, method = optim.method)

finrho<-estim1$par[1]
finsigmaV<-estim1$par[2]
# print(c(finrho,finsigmaV))

   wy <- as.matrix(Ws2 %*% y)
   yt <- y-finrho*wy
   xl <- as.matrix(Ws2 %*%  x)
   xt <- x-finrho*xl
   
	wwy <- as.matrix(Ws2 %*% wy)
	wyt<-wy-finrho*wwy


	yf<-panel.transformations(yt,indic, type= "within")
   xf<-panel.transformations(xt,indic, type= "within")
   wyf<-panel.transformations(wyt,indic, type= "within")
	xf<-xf[,-del]
	xf<-as.matrix(xf)
	#colnames(xf)<- namesx
	# wxf <- as.matrix(Ws %*% xf)

	#print(colnames(xf))
if (is.null(endog)){
	

model.fit <- spgm.tsls(yf, wyf, xf, Hinst =  Hwithin, instr = TRUE)  	
	  betaGLS <- model.fit$coefficients 
	  
	  names(betaGLS)[1]<-"lambda"
  Zp <- model.fit$Zp
  Z <- cbind(wy,x[,-del])
  Z <- Z[, qr(Z)$pivot[seq_len(qr(Z)$rank)]]
  fv<- Z %*% betaGLS
  egls<-y - fv
  SGLS<-sum(egls^2)/(N-1)
  xfpxfNT<-(1/NT)*crossprod(Zp)/finsigmaV
  PSI<-solve(xfpxfNT)
  covbeta<-PSI/NT



  errcomp<-rbind(finrho,finsigmaV)
  nam.errcomp <- c("rho","sigma^2_v")
  rownames(errcomp) <- nam.errcomp
  colnames(errcomp)<-"Estimate"
model.data <- data.frame(cbind(y,as.matrix(x)))
#print(betaGLS)
  type <- "Spatial fixed effects  SARAR model (GM estimation)"
    spmod <- list(coefficients=betaGLS, errcomp=errcomp,
                vcov=covbeta, vcov.errcomp=NULL,
                residuals=as.numeric(egls), fitted.values=fv,
                sigma2=SGLS,type=type, rho=errcomp, model=model.data, logLik=NULL, coy=yt, cox=xt, rhs=k, type = type)

	}

else{
	
	endogl<- as.matrix(Ws2 %*% endog)
   endogt <- endog - finrho*endogl

	endogf<-panel.transformations(endogt,indic, type= "within")

  yend<-cbind(endogf, wyf)
  
model.fit <- spgm.tsls(yf, yend, xf, Hinst = Hwithin, instr = TRUE)  	
	  betaGLS <- model.fit$coefficients 
	  names(betaGLS)[ncol(endog)+1]<-"lambda"
#print(betaGLS)

  Zp<-model.fit$Zp
  Z<-cbind(endog, wy, x[,-del])
  fv<-as.matrix(Z)[, qr(Z)$pivot[seq_len(qr(Z)$rank)]] %*% betaGLS
  egls<-y - fv
  SGLS<-sum(egls^2)/(N-1)
  xfpxfNT<-(1/NT)*crossprod(Zp)/finsigmaV
  PSI<-solve(xfpxfNT)
  covbeta<-PSI/NT
  
  errcomp<-rbind(finrho,finsigmaV)
  nam.errcomp <- c("rho","sigma^2_v")
  rownames(errcomp) <- nam.errcomp
  colnames(errcomp)<-"Estimate"
model.data <- data.frame(cbind(y,as.matrix(x)))

  type <- "Spatial fixed effects SARAR model with additional endogenous variables (GM estimation)"
    spmod <- list(coefficients=betaGLS, errcomp=  errcomp,
                vcov=  covbeta, vcov.errcomp=NULL,
                residuals=as.numeric(egls), fitted.values=fv,
                sigma2=SGLS,type=type, rho=errcomp, model=model.data, logLik=NULL, coy=yt, cox=xt, rhs=k, type = type)

	}

		},
		
	random = {

if(is.null(endog)){

	result1<-ivplm.w2sls(Y = y,X =x, lag = TRUE, 
	                     listw = Ws, listw2 = Ws2, 
	                     twow = twow, lag.instruments = lag.instruments, 
	                     t = t, N = N, NT = NT, Durbin = Durbin, xdur =xdur)
	result2<-ivplm.b2sls(Y = y,X =x, lag = TRUE, listw = Ws, 
	                     listw2 = Ws2, twow = twow, 
	                     lag.instruments = lag.instruments, 
	                     t = t, N = N, NT = NT, Durbin = Durbin, xdur = xdur)
	
	xw   <- result1$xdu[, qr(result1$xdu)$pivot[seq_len(qr(result1$xdu)$rank)]]
	xb   <- result2$xdu[, qr(result2$xdu)$pivot[seq_len(qr(result2$xdu)$rank)]]
	x <- cbind(xw, xb)
	x <- x[,unique(colnames(x))]
	
	}
	
	else{
		result1<-ivplm.w2sls(Y = y,X =x, H = instruments, 
		                     endog = endog, lag = TRUE, listw = Ws, 
		                     listw2 = Ws2, twow = twow, 
		                     lag.instruments = lag.instruments, 
		                     t = t, N = N, NT = NT, Durbin = Durbin, xdur = xdur)
		result2<-ivplm.b2sls(Y = y,X =x, H = instruments, 
		                     endog = endog, lag = TRUE, listw = Ws, 
		                     listw2 = Ws2, twow = twow, 
		                     lag.instruments = lag.instruments, 
		                     t = t, N = N, NT = NT, Durbin = Durbin, xdur = xdur)

		xw   <- result1$xdu[, qr(result1$xdu)$pivot[seq_len(qr(result1$xdu)$rank)]]
		xb   <- result2$xdu[, qr(result2$xdu)$pivot[seq_len(qr(result2$xdu)$rank)]]
		x <- cbind(xw, xb)
		x <- x[,unique(colnames(x))]
	}

Hwithin <- result1$Hwithin
Hbetween <- result2$Hbetween
Hbetweennt<-matrix(,NT, ncol(Hbetween))
for (i in 1:ncol(Hbetween)) Hbetweennt[,i]<-rep(Hbetween[,i], t)

res1<-as.matrix(as.numeric(residuals(result1)))
res2<-as.matrix(as.numeric(residuals(result2)))

Gg<-fswithin(Ws2, res1, N, t)

if(is.null(pars)) {
    wres <- as.matrix(Ws2 %*% res1)
    r.init <- solve(crossprod(res1),crossprod(res1,wres))
	v.init <- crossprod(res1)/NT	
	pars <- c(r.init, v.init)	
}


if (optim.method == "nlminb")  estim1 <- nlminb(pars, arg, v = Gg, verbose = verbose, control = control, lower=c(-0.999,0), upper=c(0.999,Inf))

else estim1 <- optim(pars, arg, v = Gg, verbose = verbose, method = optim.method)

Wres2 <- as.matrix(listw2nn %*% res2)
urhoWu<-res2 - estim1$par[1] * Wres2
finsigma1<-crossprod(urhoWu)/N


switch(moments, 
	  
	  initial = {
	  	
	finrho=estim1$par[1]
	finsigmaV=estim1$par[2]
	finsigma1=finsigma1
	
    }, 
    
    weights = {
    	
    	Ggw<-pwbetween(bigG=Gg$bigG, smallg=Gg$smallg, u=res2, N=N, t=t, TR=Gg$TR, listw=listw2nn)
      pars2<-c(estim1$par[1],estim1$par[2],finsigma1)

if (optim.method == "nlminb")  estim2 <- nlminb(pars2, arg1, v = Ggw,t=t, ss= estim1$par[2], SS=finsigma1 ,verbose = verbose, control = control, lower=c(-0.999,0,0), upper=c(0.999,Inf,Inf))

else      estim2 <- optim(pars2, arg1, v = Ggw, t=t, ss= estim1$par[2], SS=finsigma1 ,verbose = verbose, method = optim.method)
      
      
	finrho=estim2$par[1]
	finsigmaV=estim2$par[2]
	finsigma1=estim2$par[3]

    }, 
    
    fullweights = {

	   Ggw<-pwbetween(bigG=Gg$bigG, smallg=Gg$smallg, u=res2, N=N,t=t, TR=Gg$TR, listw = listw2nn)
      weights<-tw(W=listw2nn, N)
      pars2<-c(estim1$par[1],estim1$par[2],finsigma1)
     # 
if(optim.method == "nlminb")   estim3 <-nlminb(pars2, arg2, v = Ggw, t = t, ss= estim1$par[2], SS=finsigma1 ,TW = weights$TW, verbose = verbose, control = control, lower=c(-0.999,0,0), upper=c(0.999,Inf,Inf))
else estim3 <-optim(pars2, arg2, v = Ggw, t = t, ss= estim1$par[2], SS=finsigma1 ,TW = weights$TW, verbose = verbose, method = optim.method)

	finrho=estim3$par[1]
	finsigmaV=estim3$par[2]
	finsigma1=estim3$par[3]


	
		}, 
		stop("...\nUnknown method\n"))



 theta<- 1-(sqrt(finsigmaV)/sqrt(finsigma1))
 
	wy<- as.matrix(Ws2 %*% y)
   yt <- y-finrho*wy
   xl<- as.matrix(Ws2 %*% x)
   xt <- x-finrho*xl
   wwy<- as.matrix(Ws2 %*% wy)
	wyt<-wy - finrho*wwy

	
	

  ytmt<-tapply(yt, indic, mean)
  ytNT<-rep(ytmt, t)
  yf<-(yt - as.numeric(theta)*ytNT)
  
  dm1<- function(A) rep(unlist(tapply(A, indic, mean, simplify=TRUE)), t)
  xtNT<-apply(xt,2,dm1)
  xf<-(xt - as.numeric(theta)*xtNT)
  
  wytmt<-tapply(wyt, indic, mean)
  wytNT<-rep(wytmt, t)  
  wyf<-wyt - as.numeric(theta)*wytNT

Hgls <- cbind(1, result1$Xwithin,result2$Xbetweennt, Hwithin, Hbetweennt)	
Hgls <- as.matrix(Hgls[, qr(Hgls)$pivot[seq_len(qr(Hgls)$rank)]])

if (is.null(endog)){
model.fit <- spgm.tsls(yf, as.matrix(wyf), as.matrix(xf),  as.matrix(Hgls))  	
	  betaGLS <- model.fit$coefficients 
names(betaGLS)[1]<-"lambda"
#print(betaGLS)

  Zp<-model.fit$Zp
  Z<-cbind(wy,x)
  fv<-as.matrix(Z)[, qr(Z)$pivot[seq_len(qr(Z)$rank)]] %*% as.matrix(betaGLS)
  egls<-y - fv
  SGLS<-sum(egls^2)/(N-1)
  xfpxfNT<-(1/NT)*crossprod(Zp)/finsigmaV
  PSI<-solve(xfpxfNT)
  covbeta<-PSI/NT
	

  errcomp<-rbind(finrho,finsigmaV,finsigma1,theta)
  nam.errcomp <- c("rho","sigma^2_v",'sigma^2_1',"theta")
  rownames(errcomp) <- nam.errcomp
  colnames(errcomp)<-"Estimate"
model.data <- data.frame(cbind(y,as.matrix(x)))

  type <- "Spatial random effects SARAR model (GM estimation)"
    spmod <- list(coefficients=betaGLS, errcomp=errcomp,
                vcov=covbeta, vcov.errcomp=NULL,
                residuals=as.numeric(egls), fitted.values=fv,
                sigma2=SGLS,type=type, rho=errcomp, model=model.data, logLik=NULL, coy=yt, cox=xt, rhs=k, type = type)

	}

else{
	
  endogl<-as.matrix(Ws2 %*% endog)
  endogt <- endog - finrho*endogl

  endogtNT<-apply(as.matrix(endogt), 2, dm1)
  endogf<-(endogt - as.numeric(theta)*endogtNT)


  yend<-cbind(endogf, wyf)  
  
model.fit <- spgm.tsls(yf, as.matrix(yend), as.matrix(xf),  as.matrix(Hgls))  	
	  betaGLS <- model.fit$coefficients 
names(betaGLS)[ncol(as.matrix(endog)) + 1 ]<-"lambda"
#print(betaGLS)

  Zp<-model.fit$Zp
  Z<-cbind(endog, wy,x)
  fv<-as.matrix(Z)[, qr(Z)$pivot[seq_len(qr(Z)$rank)]] %*% as.matrix(betaGLS)
  egls<-y - fv
  SGLS<-sum(egls^2)/(N-1)
  xfpxfNT<-(1/NT)*crossprod(Zp)/finsigmaV
  PSI<-solve(xfpxfNT)
  covbeta<-PSI/NT
  

  errcomp<-rbind(finrho,finsigmaV,finsigma1,theta)
  nam.errcomp <- c("rho","sigma^2_v",'sigma^2_1',"theta")
  rownames(errcomp) <- nam.errcomp
  colnames(errcomp)<-"Estimate"
model.data <- data.frame(cbind(y,as.matrix(x)))

  type <- "Spatial random effects SARAR model with additional endogenous variables (GM estimation)"
    spmod <- list(coefficients=betaGLS, errcomp=  errcomp,
                vcov=  covbeta, vcov.errcomp=NULL,
                residuals=as.numeric(egls), fitted.values=fv,
                sigma2=SGLS,type=type, rho=errcomp, model=model.data,
                logLik=NULL, coy=yt, cox=xt, rhs=k, type = type, 
                nfimpacts = "sarar_gm")
  
	}
		}, 
		stop("...\nUnknown method\n"))	
	
	
  return(spmod)
  

	}
	
	





	
	
