`spgm` <-
function(formula, data=list(), index=NULL, listw,
         model=c("within","random"), lag=FALSE, spatial.error=FALSE,
         moments = c("initial", "weights", "fullweights"), endog = NULL, instruments= NULL, verbose = FALSE, method = c("w2sls", "b2sls", "g2sls", "ec2sls"), control = list()){

## translation for uniformity
effects <- switch(match.arg(model), within="fixed", random="random")

# source('tss.R')
# source('sumres.R')
# source('utilities_GM.R')
# source('listw2dgCMatrix.R')
# source('summary.splm.R')
# source('print.splm.R')
# source('print.summary.splm.R')
#
#
#
	# source("ivsplm.R")
	# source("ivplm.w2sls.R")
	# source("ivplm.b2sls.R")	
	# source("ivplm.g2sls.R")
	# source("ivplm.ec2sls.R")
#	source("sperrorgm2.R")		
##	source("spsarargm3.R")		##working with the Mutl and Pfaff. procedure
#	source("spsarargm4.R")	#works fine with a modified procedure
#source("spsarargm5.R")	


if(model == "within" && attr(terms(formula), "intercept") == 0 ) formula <- as.formula(paste(attr(terms(formula),"variables")[1+attr(terms(formula),"response")], paste(attr(terms(formula),"term.labels"), collapse="+"), sep="~"))

	
	


cl<-match.call()
if(!spatial.error){
	
	results<-ivsplm(formula = formula, effects = effects, data=data, index = index, endog = endog, instruments = instruments, method = method, lag = lag, listw = listw)
	
	}



else{
	
	
if(!lag) results <- sperrorgm(formula = formula, data = data, index = index, listw = listw, moments = moments, endog = endog, instruments = instruments, verbose = verbose, effects = effects, control = control)
#, initial.GMerror = initial.GMerror
else results <- spsarargm(formula = formula, data = data, index = index, listw = listw, moments = moments, lag = lag, endog = endog, instruments = instruments, verbose = verbose, effects = effects, control = control)

	}

results$call <- cl
results$ef.sph<- effects
results$legacy <- c(lag, spatial.error)
results

}



lag.listwpanel<-function(arg1, arg2,indicator){
	
	if(!is.matrix(arg2)) arg2<-as.matrix(arg2)

`inter` <- function(arg4) hold <- unlist(tapply(arg4, indicator, function(arg3) lag.listw(arg1,arg3)))

out<-apply(arg2, 2,inter)
	}
	
	


panel.transformations<-function(arg1, indicator, type=c("within","between","both")){

	type<-match.arg(type)
 	if(!is.matrix(arg1)) arg1<-as.matrix(arg1)


spms<-function(q) tapply(q,indicator,mean)

switch(type, between = {

	arg1_b<-matrix(apply(arg1,2,spms), length(unique(indicator)), ncol(arg1))
	out<-arg1_b
		
		}, 
		within ={			

	arg1_b<-matrix(apply(arg1,2,spms), length(unique(indicator)), ncol(arg1))
	arg1_bnt<- apply(arg1_b, 2, function(bb) rep(bb, each = length(indicator)/ length(unique(indicator)) ))             
	arg1_w	<- arg1 - arg1_bnt
	out<-arg1_w
		
			}, 			
			both = {
				
	arg1_b<-matrix(apply(arg1,2,spms), length(unique(indicator)), ncol(arg1))
	arg1_bnt<- apply(arg1_b, 2, function(bb) rep(bb, each = length(indicator)/ length(unique(indicator)) ))  
	arg1_w	<- arg1 - arg1_bnt
	out<-list(arg1_w, arg1_b)			
				
				}		
		)	

out


	}


spgm.tsls<-function(y, yend, X, Zinst, Hinst= NULL, instr = FALSE){
	yend<-as.matrix(yend)

if (instr) H<-Hinst	
else    H <- cbind(X, Zinst)

   Z <- cbind(yend, X)
  Znames<- colnames(Z) 
	 yendp<-matrix(,nrow(yend), ncol(yend))
for(i in 1:ncol(yend))    yendp[,i] <- fitted.values(lm(yend[,i] ~ H -1 ))
    Zp <- cbind(yendp,X)
    #print(Zp[1:5,])
    model.fit<-lm(y~Zp-1)
    biv <- coefficients(model.fit)
#print(biv)

if(any(is.na(biv)))  yp <- Z[,-which(is.na(biv))] %*% biv[-which(is.na(biv))]
else yp <- Z %*% biv



        ehat <- y-yp
        sse <- c(crossprod(ehat, ehat))
        df <- model.fit$df.residual
        s2 <- sse/df
        
if(any(is.na(biv)))   Zp<-Zp[,-which(is.na(biv))]  


ZpZi<-solve(crossprod(Zp))   

         varb<-ZpZi *  s2



names(biv)<-Znames

if(any(is.na(biv))) biv<-biv[-which(is.na(biv))] 
else	biv<-biv



        result <- list(coefficients = biv, var = varb, sse = sse, 
            residuals = as.numeric(ehat), df = df, Zp = Zp)
    
    result
}




sperrorgm<-function(formula, data = list(), index = NULL, listw , moments = c("initial", "weights", "fullweights"), endog = NULL, instruments = NULL, verbose = FALSE, effects = c("fixed","random"), control = list() ){

effects<-match.arg(effects)
moments<-match.arg(moments)
	indes<-index
if(!is.null(index)) {
    require(plm)
    data <- plm.data(data, index)
    }
    
  index <- data[,1]
  tindex <- data[,2]

  cl <- match.call()
  if(dim(data)[[1]]!=length(index)) stop("Non conformable arguments")


    mt <- terms(formula, data = data)
    mf <- lm(formula, data, na.action = na.fail, method = "model.frame")

    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)


   
  names(index)<-row.names(data)
  ind<-index[which(names(index)%in%row.names(x))]
  tind<-tindex[which(names(index)%in%row.names(x))]
  N<-length(unique(ind))
  k<-dim(x)[[2]]
  T<-max(tapply(x[,1],ind,length))
  NT<-length(ind)
  indt<-rep(seq(1,N), T)
  tindt<-rep(seq(1,T), each = N)
  tord<-order(indt,tindt)
  
  balanced<-N*T==NT
if(!balanced) stop("Estimation method unavailable for unbalanced panels")

if (is.matrix(listw)){
        if(dim(listw)[[1]] != N) 
            stop("Non conformable spatial weights")
        require(spdep)
        listw <- mat2listw(listw)
    }
    else{
if(length(listw$neighbours)!=N	)stop("Non conformable spatial weights")
    	}

if (!inherits(listw, c("listw", "matrix"))) 
        stop("listw should be either a matrix of an object of class listw")


if(!is.null(endog)){
	if(is.character(endog)){
	xend<- match(endog, colnames(data))
	endog <- data[,xend]
}
kend <- ncol(as.matrix(endog)) 

if(is.null(instruments)) stop("No instruments specified")
if(is.character(instruments)){
	inst <- match(instruments, colnames(data))
	instruments <- data[,inst]
	}
kinst <- ncol(instruments) 

if(kinst < kend) stop("The model is not identified: Not enough instruments specified")
	}

require(plm)	

transx<-panel.transformations(x, ind, type= "both")
Xbetween<-transx[[2]]
Xwithin<-transx[[1]]



del<-which(diag(var(as.matrix(Xwithin)))==0)

if (colnames(x)[1] == "(Intercept)") Xbetween <- Xbetween[,-1]
delb<-as.numeric(which(diag(var(Xbetween))==0))
if(length(delb)==0) Xbetween<-Xbetween
else Xbetween<-Xbetween[,-delb]
Xbetween<-cbind(1,Xbetween)



WX <- lag.listwpanel(listw, x, tind)
WX<-WX[tord,]
WWX <- lag.listwpanel(listw, WX, tind)
WWX<-WWX[tord,]
WX<-WX[,-del]
WWX<-WWX[,-del]
HX<-cbind(WX, WWX)

#

WXwithin <- lag.listwpanel(listw, Xwithin, tind)
WXwithin<-WXwithin[tord,]
WWXwithin <- lag.listwpanel(listw, WXwithin, tind)
WWXwithin<-WWXwithin[tord,]
WXwithin<-WXwithin[,-del]
WWXwithin<-WWXwithin[,-del]


spms<-function(q) tapply(q,ind,mean)
Xbetweennt<-matrix(,NT,ncol(Xbetween))
for (i in 1:ncol(Xbetween)) Xbetweennt[,i]<-rep(Xbetween[,i],each=T)
if (colnames(x)[1] == "(Intercept)") Xbetweennt <- Xbetweennt[,-1]

WXbetween <- lag.listwpanel(listw, Xbetweennt, tind)
WXbetween <- lag.listwpanel(listw, Xbetweennt, tind)
if(length(delb)==0) WXbetween<-WXbetween
else WXbetween<-WXbetween[,-delb]
WXbetween<-WXbetween[tord,]
WWXbetween <- lag.listwpanel(listw, WXbetween, tind)
if(length(delb)==0) WWXbetween<-WWXbetween
else WWXbetween<-WWXbetween[,-delb]
WWXbetween<-WWXbetween[tord,]


Hwithin<-cbind(Xwithin[,-del],WXwithin, WWXwithin)
Hbetween<-cbind(1, Xbetweennt,WXbetween, WWXbetween)
Hgls<-cbind(1, Hwithin, Hbetween[,-1])


	switch(effects, 

	fixed = {

formula1<- formula(paste(formula[2],"~",formula[3],-1))
if(is.null(endog))	 result<-plm(formula = formula1, data = data, effect = "individual", model = "within") 
else 	result<-ivplm.w2sls(Y = y,X =x, H = instruments, endog = endog, ind = ind, tind = tind, lag = FALSE, listw = listw)
#print(result$coefficients)

residuals<-as.matrix(as.numeric(residuals(result)))
oo<-order(tind,ind)
res<-residuals[oo,]

Gg<-fswithin(listw,res,N,T)
pars<-c(0.2, 0.5)
#estim1 <- nlminb(pars, arg, v = Gg, verbose = verbose, control = control, lower=c(-0.999,0), upper=c(0.999,Inf))

estim1 <- optim(pars, arg, v = Gg, verbose = verbose)

	finrho=estim1$par[1]
	finsigmaV=estim1$par[2]
#print(c(finrho,finsigmaV))	

	wy<-lag.listwpanel(listw, y, tind)
	wy<-wy[tord]
   yt <- y-finrho*wy
   xl<- lag.listwpanel(listw, x, tind)
   xl<-xl[tord,]
   xt <- x-finrho*xl


	yf<-panel.transformations(yt,ind, type= "within")
    xf<-panel.transformations(xt,ind, type= "within")
	xf<-xf[,-del]
	xf<-as.matrix(xf)
	colnames(xf)<-colnames(x)[-del]

if (is.null(endog)){

result<-lm(as.matrix(yf)~as.matrix(xf)-1)
vcov<-vcov(result)
betaGLS<-coefficients(result)

	names(betaGLS)<-colnames(xf)
  errcomp<-rbind(finrho,finsigmaV)
  nam.errcomp <- c("rho","sigma^2_v")
    rownames(errcomp) <- nam.errcomp
  colnames(errcomp)<-"Estimate"
   model.data <- data.frame(cbind(y,x[,-1]))

  type <- "fixed effects GM"
    spmod <- list(coefficients= betaGLS, errcomp=errcomp,
                vcov=vcov, vcov.errcomp=NULL,
                residuals=residuals(result), fitted.values=(y-as.vector(residuals(result))),
                sigma2=crossprod(residuals(result))/result$df.residual, type=type, rho=errcomp, model=model.data, logLik=NULL)
  class(spmod) <- "splm"
  return(spmod)

	
	}

else{
	
	endogl<-lag.listwpanel(listw, endog, tind)
	endogl<-endogl[tord,]
   endogt <- endog - finrho*endogl

endogf<-panel.transformations(endogt,ind, type= "within")
instwithin<-panel.transformations(instruments,ind, type= "within")

  instwithin<-cbind(xf, instwithin)
#   instwithin<-cbind(Hwithin, instwithin)


  model.fit <- spgm.tsls(yf,endogf, xf, Hinst = instwithin, instr  = TRUE )  
  betaGLS <- model.fit$coefficients 

#print(betaGLS)

  Zp<-model.fit$Zp
  Z<-cbind(endog,x[,-del])
  fv<-as.matrix(Z) %*% as.matrix(betaGLS)
  egls<-y - fv
  SGLS<-sum(egls^2)/(N-1)
  xfpxfNT<-(1/NT)*crossprod(Zp)/finsigmaV
  PSI<-solve(xfpxfNT)
  covbeta<-PSI/NT

  nam.beta <- c(colnames(endog), colnames(x)[-del])
  names(betaGLS) <- nam.beta
  errcomp<-rbind(finrho,finsigmaV)
  nam.errcomp <- c("rho","sigma^2_v")
   model.data <- data.frame(cbind(y,x[,-1]))

  type <- "fixed effects GM"
    spmod <- list(coefficients=betaGLS, errcomp=  errcomp,
                vcov=covbeta, vcov.errcomp=NULL,
                residuals=as.vector(egls), fitted.values=fv,
                sigma2=SGLS, type=type, rho=errcomp[1], model=model.data, logLik=NULL)
  class(spmod) <- "splm"
  return(spmod)
	
	
	}
 return(spmod)
		},
		
	random = {
		
if(is.null(endog))	{

result<-lm(y~x-1) 

residuals<-as.matrix(as.numeric(residuals(result)))
oo<-order(tind,ind)
res<-residuals[oo,]

Gg<-fs(listw,res,N,T)
pars<-c(0,0)
#estim1 <- nlminb(pars, arg, v = Gg, verbose = verbose, control = control, lower=c(-0.999,0), upper=c(0.999,Inf))
estim1 <- optim(pars, arg, v = Gg, verbose = verbose)

urub<-res- estim1$par[1]*Gg$ub
Q1urQ1ub<-Gg$Q1u - estim1$par[1]*Gg$Q1ub
S1<- crossprod(urub, Q1urQ1ub)/N

switch(moments, 
	  
	  initial = {
	  	
	finrho=estim1$par[1]
	finsigmaV=estim1$par[2]
	finsigma1=S1
	
    }, 
    
    weights = {
    	
  	   Ggw<-pw(bigG=Gg$bigG, smallg=Gg$smallg, Q1u=Gg$Q1u,Q1ub=Gg$Q1ub,Q1ubb=Gg$Q1ubb, u=res, ub=Gg$ub,ubb=Gg$ubb,N=N, TR=Gg$TR)
      pars2<-c(estim1$par[1],estim1$par[2],S1)
#      estim2 <- nlminb(pars2, arg1, v = Ggw,T=T,ss=estim1$par[2] ,SS=S1, verbose = verbose, control = control, lower=c(-0.999,0,0), upper=c(0.999,Inf,Inf))
#
      estim2 <- optim(pars2, arg1, v = Ggw,T=T,ss=estim1$par[2] ,SS=S1, verbose = verbose)

	finrho=estim2$par[1]
	finsigmaV=estim2$par[2]
	finsigma1=estim2$par[3]

    }, 
    
    fullweights = {

	   Ggw<-pw(bigG=Gg$bigG, smallg=Gg$smallg, Q1u=Gg$Q1u,Q1ub=Gg$Q1ub,Q1ubb=Gg$Q1ubb, u=res, ub=Gg$ub,ubb=Gg$ubb,N=N, TR=Gg$TR)
      weights<-tw(W=listw,N)
      pars2<-c(estim1$par[1],estim1$par[2],S1)
#      estim3 <-nlminb(pars2, arg2, v = Ggw, T=T,ss=estim1$par[2] ,SS=S1, TW=weights$TW, verbose = verbose, control = control, lower=c(-0.999,0,0), upper=c(0.999,Inf,Inf))
#
      estim3 <-optim(pars2, arg2, v = Ggw, T=T,ss=estim1$par[2] ,SS=S1, TW=weights$TW, verbose = verbose)
      
	finrho=estim3$par[1]
	finsigmaV=estim3$par[2]
	finsigma1=estim3$par[3]


	
		}, 
		stop("...\nUnknown method\n"))
	
	
	}
else{
	
		result1<-ivplm.w2sls(Y = y,X =x, H = instruments, endog = endog,  ind = ind, tind = tind, lag = FALSE, listw = listw)

		result2<-ivplm.b2sls(Y = y,X =x, H = instruments, endog = endog,  ind = ind, tind = tind, lag = FALSE, listw = listw)
#print(result2$coefficients)

residuals1<-as.matrix(as.numeric(residuals(result1)))
oo<-order(tind,ind)
res1<-residuals1[oo]
res2<-as.matrix(as.numeric(residuals(result2)))

Gg<-fswithin(listw,res1,N,T)

pars<-c(0.2, result1$sigmav)
#estim1 <- nlminb(pars, arg, v = Gg, verbose = verbose, control = control, lower=c(-0.999, 0), upper=c(0.999,Inf))

estim1 <- optim(pars, arg, v = Gg, verbose = verbose)


Wres2 <-lag.listw(listw, res2)
urhoWu<-res2 - estim1$par[1] * Wres2
finsigma1<-crossprod(urhoWu)/N


switch(moments, 
	  
	  initial = {
	  	
	finrho=estim1$par[1]
	finsigmaV=estim1$par[2]
	finsigma1=finsigma1
	
    }, 
    
    weights = {
    	
    	Ggw<-pwbetween(bigG=Gg$bigG, smallg=Gg$smallg, u=res2, N=N, T=T, TR=Gg$TR, listw=listw)
      pars2<-c(estim1$par[1],estim1$par[2],finsigma1)
#      estim2 <- nlminb(pars2, arg1, v = Ggw,T=T, ss= estim1$par[2], SS=finsigma1 , verbose = verbose, control = control, lower=c(-0.999,0,0), upper=c(0.999,Inf,Inf))

      estim2 <- optim(pars2, arg1, v = Ggw,T=T, ss= estim1$par[2], SS=finsigma1 , verbose = verbose)
      
	finrho=estim2$par[1]
	finsigmaV=estim2$par[2]
	finsigma1=estim2$par[3]

    }, 
    
    fullweights = {

	   Ggw<-pwbetween(bigG=Gg$bigG, smallg=Gg$smallg, u=res2, N=N,T=T, TR=Gg$TR, listw = listw)
      weights<-tw(W=listw,N)
      pars2<-c(estim1$par[1],estim1$par[2],finsigma1)
#      estim3 <-nlminb(pars2, arg2, v = Ggw, T = T,ss= estim1$par[2], SS=finsigma1 , TW = weights$TW, verbose = verbose, control = control, lower=c(-0.999,0,0), upper=c(0.999,Inf,Inf))
 
       estim3 <-optim(pars2, arg2, v = Ggw, T = T,ss= estim1$par[2], SS=finsigma1 , TW = weights$TW, verbose = verbose)

	finrho=estim3$par[1]
	finsigmaV=estim3$par[2]
	finsigma1=estim3$par[3]


	
		}, 
		stop("...\nUnknown method\n"))
	}

theta<- 1-(sqrt(finsigmaV)/sqrt(finsigma1))	
wy <- lag.listwpanel(listw, y, tind)
wy<-wy[tord]
yt <- y-finrho*wy
xl<- lag.listwpanel(listw, x, tind)
xl<-xl[tord,]
xt <- x-finrho*xl

 
ytmt<-tapply(yt, ind, mean)
ytNT<-rep(ytmt, each = T)
yf<-(yt - theta*ytNT)

dm1<- function(A) rep(unlist(tapply(A, ind, mean, simplify=TRUE)), each= T)
xtNT<-apply(xt,2,dm1)
xf<-(xt - as.numeric(theta)*xtNT)


if (is.null(endog)){
	
  xfpxf<-crossprod(xf)
  xfpxfi<-solve(xfpxf)
  betaGLS<-xfpxfi%*%crossprod(xf,yf) 
  betaGLS<- as.numeric(betaGLS)
#print(betaGLS)  
  fv<-as.vector(x %*% betaGLS)
  egls<-y - fv
  SGLS<-sum(egls^2)/(N-1)
  xfpxfNT<-(1/NT)*xfpxf/finsigmaV
  PSI<-solve(xfpxfNT)
  covbeta<-PSI/NT

  errcomp<-rbind(finrho,finsigmaV,finsigma1,theta)
  nam.beta <- dimnames(x)[[2]]
  nam.errcomp <- c("rho","sigma^2_v",'sigma^2_1',"theta")
  names(betaGLS) <- nam.beta
  rownames(errcomp) <- nam.errcomp
  colnames(errcomp)<-"Estimate"
model.data <- data.frame(cbind(y,x))
sigma2 <- SGLS
  type <- "random effects GM"
    spmod <- list(coefficients=betaGLS, errcomp=errcomp,
                vcov=covbeta, vcov.errcomp=NULL,
                residuals=as.vector(egls), fitted.values=fv,
                sigma2=sigma2,type=type, rho=errcomp, model=model.data,
                call=cl, logLik=NULL, coy=yt, cox=xt, rhs=k)
  class(spmod) <- "splm"
  return(spmod)

	}

else{
	
	endogl<-lag.listwpanel(listw, endog, tind)
	endogl<-endogl[tord,]
   endogt <- endog - finrho*endogl

endogt<-as.matrix(endogt)
endogtNT<-apply(endogt,2,dm1)
endogf<-(endogt - as.numeric(theta)*endogtNT)

instt<-panel.transformations(instruments,ind, type= "both")
instbetween<-instt[[2]]
instwithin<-instt[[1]]

instbetweennt<-matrix(,NT,ncol(instbetween))
for (i in 1:ncol(instbetween)) instbetweennt[,i]<-rep(instbetween[,i],each=T)

 Zf<-cbind(endog,x)

#Hins<-cbind(xf,Xwithin,Xbetweennt,instwithin,instbetweennt)
#Hgls
Hins<-cbind(xf,instwithin,instbetweennt)
#Hins<-cbind(Hgls, instwithin, instbetweennt)

  model.fit <- spgm.tsls(yf,endogf,xf, Hinst = Hins, instr  = TRUE )  
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
model.data <- data.frame(cbind(y,x))
sigma2 <- SGLS
  type <- "random effects GM"
    spmod <- list(coefficients=betaGLS, errcomp=errcomp,
                vcov=covbeta, vcov.errcomp=NULL,
                residuals=as.vector(egls), fitted.values=fv,
                sigma2=sigma2,type=type, rho=errcomp, model=model.data,
                call=cl, logLik=NULL, coy=yt, cox=xt, rhs=k)
  class(spmod) <- "splm"
  return(spmod)

}

	
    }, 

    
    stop("...\nUnknown method\n"))
	
	
return(spmod)	
	
	}





spsarargm<-function(formula, data = list(), index = NULL, listw , moments = c("initial", "weights", "fullweights"), lag= FALSE, endog = NULL, instruments = NULL, verbose = FALSE, effects = c("fixed","random"), control = list()){

effects<-match.arg(effects)
moments<-match.arg(moments)


	
if(!is.null(index)) {
    require(plm)
    data <- plm.data(data, index)
    }
  index <- data[,1]
  tindex <- data[,2]
  cl <- match.call()
if(dim(data)[[1]]!=length(index)) stop("Non conformable arguments")
  
    mt <- terms(formula, data = data)
    mf <- lm(formula, data, na.action = na.fail, method = "model.frame")

    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
    
  names(index)<-row.names(data)
  ind<-index[which(names(index)%in%row.names(x))]
  tind<-tindex[which(names(index)%in%row.names(x))]
  N<-length(unique(ind))
  k<-dim(x)[[2]]
  T<-max(tapply(x[,1],ind,length))
  NT<-length(ind)
  indt<-rep(seq(1,N), T)
  tindt<-rep(seq(1,T), each = N)
  tord<-order(indt,tindt)

  balanced<-N*T==NT
if(!balanced) stop("Estimation method unavailable for unbalanced panels")

if (is.matrix(listw)){
        if(dim(listw)[[1]] != N) 
            stop("Non conformable spatial weights")
        require(spdep)
        listw <- mat2listw(listw)
    }
    else{
if(length(listw$neighbours)!=N	)stop("Non conformable spatial weights")
			}

if (!inherits(listw, c("listw", "matrix"))) 
        stop("listw should be either a matrix of an object of class listw")



if(!is.null(endog)){
	if(is.character(endog)){
	xend<-match(endog, colnames(data))
	endog <- data[,xend]
	endog<-as.matrix(endog)
if(is.null(instruments)) warnings("Spatially lagged exogenous variables will be used as instruments")
kend <- ncol(as.matrix(endog))  + 1
	}
	
if(is.character(instruments)){
	inst<-match(instruments, colnames(data))
	instruments <- data[,inst]
	instruments<-as.matrix(instruments)
kinst <- ncol(instruments)  + ifelse(colnames(x)[1]=="(Intercept)", 1, 0) + ncol(x[,-ifelse(colnames(x)[1]=="(Intercept)", 1, 0)])*3
	}

else kinst <-  ifelse(colnames(x)[1]=="(Intercept)", 1, 0) + ncol(x[,-ifelse(colnames(x)[1]=="(Intercept)", 1, 0)])*3

if(kinst < kend) stop("The model is not identified: Not enough instruments specified")
	}



require(plm)

##transform X	
transx<-panel.transformations(x, ind, type= "both")
Xbetween<-transx[[2]]
Xwithin<-transx[[1]]
del<-which(diag(var(as.matrix(Xwithin)))==0)


if (colnames(x)[1] == "(Intercept)") Xbetween <- Xbetween[,-1]
delb<-as.numeric(which(diag(var(Xbetween))==0))
if(length(delb)==0) Xbetween<-Xbetween
else Xbetween<-Xbetween[,-delb]
Xbetween<-cbind(1,Xbetween)
deltot<-c(del,delb)

WX <- lag.listwpanel(listw, x, tind)
WX<-WX[tord,]
WWX <- lag.listwpanel(listw, WX, tind)
WWX<-WWX[tord,]
WX<-WX[,-del]
WWX<-WWX[,-del]
HX<-cbind(WX, WWX)


WXwithin <- lag.listwpanel(listw, Xwithin, tind)
WXwithin<-WXwithin[tord,]
WWXwithin <- lag.listwpanel(listw, WXwithin, tind)
WWXwithin<-WWXwithin[tord,]
WXwithin<-WXwithin[,-del]
WWXwithin<-WWXwithin[,-del]


spms<-function(q) tapply(q,ind,mean)
Xbetweennt<-matrix(,NT,ncol(Xbetween))
for (i in 1:ncol(Xbetween)) Xbetweennt[,i]<-rep(Xbetween[,i],each=T)
if (colnames(x)[1] == "(Intercept)") Xbetweennt <- Xbetweennt[,-1]

WXbetween <- lag.listwpanel(listw, Xbetweennt, tind)
WXbetween <- lag.listwpanel(listw, Xbetweennt, tind)
if(length(delb)==0) WXbetween<-WXbetween
else WXbetween<-WXbetween[,-delb]
WXbetween<-WXbetween[tord,]
WWXbetween <- lag.listwpanel(listw, WXbetween, tind)
if(length(delb)==0) WWXbetween<-WWXbetween
else WWXbetween<-WWXbetween[,-delb]
WWXbetween<-WWXbetween[tord,]


Hwithin<-cbind(Xwithin[,-del],WXwithin, WWXwithin)
Hbetween<-cbind(1, Xbetweennt,WXbetween, WWXbetween)
Hgls<-cbind(1, Hwithin, Hbetween[,-1])



switch(effects, 

	fixed = {
		
if(is.null(endog)) result<-ivplm.w2sls(Y = y,X =x, endog = NULL, ind = ind, tind = tind, lag = TRUE, listw = listw)
else 	result<-ivplm.w2sls(Y = y,X =x, H = instruments, endog = endog, ind = ind, tind = tind, lag = TRUE, listw = listw)
#plot(density(result$residuals))
#print(result$coefficients)


residuals<-as.matrix(as.numeric(residuals(result)))
oo<-order(tind,ind)
res<-residuals[oo,]

#print(result$sigmav)

Gg<-fswithin(listw,res,N,T)
pars<-c(0, 0)
estim1 <- nlminb(pars, arg, v = Gg, verbose = verbose,  lower=c(-0.999,0), upper=c(0.9999,Inf))
#estim1 <- optim(pars, arg, v = Gg, verbose = verbose)
finrho<-estim1$par[1]
finsigmaV<-estim1$par[2]
#print(c(finrho,finsigmaV))

	wy<-lag.listwpanel(listw, y, tind)
	wy<-wy[tord]
   yt <- y-finrho*wy
   xl<- lag.listwpanel(listw, x, tind)
   xl<-xl[tord,]
   xt <- x-finrho*xl
	wwy<-lag.listwpanel(listw, wy, tind)
	wwy<-wwy[tord]
	wyt<-wy-finrho*wwy


	yf<-panel.transformations(yt,ind, type= "within")
   xf<-panel.transformations(xt,ind, type= "within")
   wyf<-panel.transformations(wyt,ind, type= "within")
	colnames(xf)<-colnames(x)


	
if (is.null(endog)){
	
	

model.fit <- spgm.tsls(yf, wyf, xf,  Hwithin)  	
	  betaGLS <- model.fit$coefficients 
	  names(betaGLS)[1]<-"lambda"

  Zp<-model.fit$Zp
  Z<-cbind(wy,x[,-del])
  fv<-as.matrix(Z) %*% betaGLS
  egls<-y - fv
  SGLS<-sum(egls^2)/(N-1)
  xfpxfNT<-(1/NT)*crossprod(Zp)/finsigmaV
  PSI<-solve(xfpxfNT)
  covbeta<-PSI/NT



  errcomp<-rbind(finrho,finsigmaV)
  nam.errcomp <- c("rho","sigma^2_v")
  rownames(errcomp) <- nam.errcomp
  colnames(errcomp)<-"Estimate"
model.data <- data.frame(cbind(y,x))

  type <- "fixed effects GM"
    spmod <- list(coefficients=betaGLS, errcomp=errcomp,
                vcov=covbeta, vcov.errcomp=NULL,
                residuals=as.numeric(egls), fitted.values=fv,
                sigma2=SGLS,type=type, rho=errcomp, model=model.data, logLik=NULL, coy=yt, cox=xt, rhs=k)
  class(spmod) <- "splm"
  return(spmod)

	}

else{
	
	endogl<-lag.listwpanel(listw, endog, tind)
	endogl<-endogl[tord,]
   endogt <- endog - finrho*endogl

	endogf<-panel.transformations(endogt,ind, type= "within")

instwithin<-panel.transformations(instruments,ind, type= "within")

instwithin<-cbind(Hwithin,instwithin)

  yend<-cbind(endogf, wyf)
  
model.fit <- spgm.tsls(yf, yend, xf, Hinst = instwithin, instr = TRUE)  	
	  betaGLS <- model.fit$coefficients 
	  names(betaGLS)[ncol(endog)+1]<-"lambda"
#print(betaGLS)

  Zp<-model.fit$Zp
  Z<-cbind(endog, wy, x[,-del])
#  print(colnames(endog))
#  print(colnames(x[,-del]))
#  print(colnames(Z))
#  print(betaGLS)
  fv<-as.matrix(Z) %*% betaGLS
  egls<-y - fv
  SGLS<-sum(egls^2)/(N-1)
  xfpxfNT<-(1/NT)*crossprod(Zp)/finsigmaV
  PSI<-solve(xfpxfNT)
  covbeta<-PSI/NT
  
  errcomp<-rbind(finrho,finsigmaV)
  nam.errcomp <- c("rho","sigma^2_v")
  rownames(errcomp) <- nam.errcomp
  colnames(errcomp)<-"Estimate"
model.data <- data.frame(cbind(y,x))

  type <- "fixed effects GM"
    spmod <- list(coefficients=betaGLS, errcomp=  errcomp,
                vcov=  covbeta, vcov.errcomp=NULL,
                residuals=as.numeric(egls), fitted.values=fv,
                sigma2=SGLS,type=type, rho=errcomp, model=model.data, logLik=NULL, coy=yt, cox=xt, rhs=k)
  class(spmod) <- "splm"
  return(spmod)
	
	
	}

		},
		
	random = {

if(is.null(endog)){
	result1<-ivplm.w2sls(Y = y,X =x, H=NULL, endog = NULL, ind = ind, tind = tind, lag = TRUE, listw = listw)
	result2<-ivplm.b2sls(Y = y,X =x, H=NULL, endog = NULL, ind = ind, tind = tind, lag = TRUE, listw = listw)

#		fsig1<-result$sigma1
#		fsigv<-result$sigmav
	}
	
	else{
		result1<-ivplm.w2sls(Y = y,X =x, H = instruments, endog = endog,  ind = ind, tind = tind, lag = TRUE, listw = listw)
		result2<-ivplm.b2sls(Y = y,X =x, H = instruments, endog = endog,  ind = ind, tind = tind, lag = TRUE, listw = listw)

#		fsig1<-result$sigma1
#		fsigv<-result$sigmav		
		}



residuals1<-as.matrix(as.numeric(residuals(result1)))
oo<-order(tind,ind)
res1<-residuals1[oo,]
res2<-as.matrix(as.numeric(residuals(result2)))
Gg<-fswithin(listw,res1,N,T)

pars<-c(0.2,0.5)
#estim1 <- nlminb(pars, arg, v = Gg, verbose = verbose, control = control, lower=c(-0.999,0), upper=c(0.999,Inf))
estim1 <- optim(pars, arg, v = Gg, verbose = verbose)

Wres2 <-lag.listw(listw, res2)
urhoWu<-res2 - estim1$par[1] * Wres2
finsigma1<-crossprod(urhoWu)/N


switch(moments, 
	  
	  initial = {
	  	
	finrho=estim1$par[1]
	finsigmaV=estim1$par[2]
	finsigma1=finsigma1
	
    }, 
    
    weights = {
    	
    	Ggw<-pwbetween(bigG=Gg$bigG, smallg=Gg$smallg, u=res2, N=N, T=T, TR=Gg$TR, listw=listw)
      pars2<-c(estim1$par[1],estim1$par[2],finsigma1)
#      estim2 <- nlminb(pars2, arg1, v = Ggw,T=T, ss= estim1$par[2], SS=finsigma1 ,verbose = verbose, control = control, lower=c(-0.999,0,0), upper=c(0.999,Inf,Inf))
      estim2 <- optim(pars2, arg1, v = Ggw,T=T, ss= estim1$par[2], SS=finsigma1 ,verbose = verbose)
      
      
	finrho=estim2$par[1]
	finsigmaV=estim2$par[2]
	finsigma1=estim2$par[3]

    }, 
    
    fullweights = {

	   Ggw<-pwbetween(bigG=Gg$bigG, smallg=Gg$smallg, u=res2, N=N,T=T, TR=Gg$TR, listw = listw)
      weights<-tw(W=listw,N)
      pars2<-c(estim1$par[1],estim1$par[2],finsigma1)
     # estim3 <-nlminb(pars2, arg2, v = Ggw, T = T, ss= estim1$par[2], SS=finsigma1 ,TW = weights$TW, verbose = verbose, control = control, lower=c(-0.999,0,0), upper=c(0.999,Inf,Inf))

estim3 <-optim(pars2, arg2, v = Ggw, T = T, ss= estim1$par[2], SS=finsigma1 ,TW = weights$TW, verbose = verbose)

	finrho=estim3$par[1]
	finsigmaV=estim3$par[2]
	finsigma1=estim3$par[3]


	
		}, 
		stop("...\nUnknown method\n"))



 theta<- 1-(sqrt(finsigmaV)/sqrt(finsigma1))
 
	wy<-lag.listwpanel(listw, y, tind)
	wy<-wy[tord]
   yt <- y-finrho*wy
   xl<- lag.listwpanel(listw, x, tind)
   xl<-xl[tord,]
   xt <- x-finrho*xl
   wwy<-lag.listwpanel(listw, wy, tind)
   wwy<-wwy[tord]
	wyt<-wy - finrho*wwy

	
	

  ytmt<-tapply(yt, ind, mean)
  ytNT<-rep(ytmt,each = T)
  yf<-(yt - theta*ytNT)
  dm1<- function(A) rep(unlist(tapply(A, ind, mean, simplify=TRUE)), each= T)
  xtNT<-apply(xt,2,dm1)
  xf<-(xt - as.numeric(theta)*xtNT)
  wytmt<-tapply(wyt, ind, mean)
  wytNT<-rep(wytmt,each = T)  
  wyf<-wyt - as.numeric(theta)*wytNT


if (is.null(endog)){
	
	
model.fit <- spgm.tsls(yf, wyf, xf,  Hgls)  	
	  betaGLS <- model.fit$coefficients 
names(betaGLS)[1]<-"lambda"
#print(betaGLS)

  Zp<-model.fit$Zp
  Z<-cbind(wy,x)
  fv<-as.matrix(Z) %*% as.matrix(betaGLS)
  egls<-y - fv
  SGLS<-sum(egls^2)/(N-1)
  xfpxfNT<-(1/NT)*crossprod(Zp)/finsigmaV
  PSI<-solve(xfpxfNT)
  covbeta<-PSI/NT
	

  errcomp<-rbind(finrho,finsigmaV,finsigma1,theta)
  nam.errcomp <- c("rho","sigma^2_v",'sigma^2_1',"theta")
  rownames(errcomp) <- nam.errcomp
  colnames(errcomp)<-"Estimate"
model.data <- data.frame(cbind(y,x))

  type <- "random effects GM"
    spmod <- list(coefficients=betaGLS, errcomp=errcomp,
                vcov=covbeta, vcov.errcomp=NULL,
                residuals=as.numeric(egls), fitted.values=fv,
                sigma2=SGLS,type=type, rho=errcomp, model=model.data, logLik=NULL, coy=yt, cox=xt, rhs=k)
  class(spmod) <- "splm"
  return(spmod)

	}

else{
	
  endogl<-lag.listwpanel(listw, endog, tind)
  endogl<-endogl[tord,]
  endogt <- endog - finrho*endogl

  endogtNT<-apply(as.matrix(endogt), 2, dm1)
  endogf<-(endogt - as.numeric(theta)*endogtNT)



transinst<-panel.transformations(instruments,ind, type= "both")
instbetween<-transinst[[2]]
instwithin<-transinst[[1]]
instbetweennt<-matrix(,NT,ncol(instbetween))
for (i in 1:ncol(instbetween)) instbetweennt[,i]<-rep(instbetween[,i],each=T)

  yend<-cbind(endogf, wyf)  
  Hgls<-cbind(Hgls,instwithin,instbetweennt)  
  
model.fit <- spgm.tsls(yf, yend, xf,  Hgls)  	
	  betaGLS <- model.fit$coefficients 
names(betaGLS)[ncol(as.matrix(endog)) + 1 ]<-"lambda"
#print(betaGLS)

  Zp<-model.fit$Zp
  Z<-cbind(endog, wy,x)
  fv<-as.matrix(Z) %*% as.matrix(betaGLS)
  egls<-y - fv
  SGLS<-sum(egls^2)/(N-1)
  xfpxfNT<-(1/NT)*crossprod(Zp)/finsigmaV
  PSI<-solve(xfpxfNT)
  covbeta<-PSI/NT
  

  errcomp<-rbind(finrho,finsigmaV,finsigma1,theta)
  nam.errcomp <- c("rho","sigma^2_v",'sigma^2_1',"theta")
  rownames(errcomp) <- nam.errcomp
  colnames(errcomp)<-"Estimate"
model.data <- data.frame(cbind(y,x))

  type <- "random effects GM"
    spmod <- list(coefficients=betaGLS, errcomp=  errcomp,
                vcov=  covbeta, vcov.errcomp=NULL,
                residuals=as.numeric(egls), fitted.values=fv,
                sigma2=SGLS,type=type, rho=errcomp, model=model.data, logLik=NULL, coy=yt, cox=xt, rhs=k)
  class(spmod) <- "splm"
  return(spmod)

	}
		}, 
		stop("...\nUnknown method\n"))	
	
	
  return(spmod)
  

	}
	
	
	
	
	
	
