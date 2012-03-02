
ivsplm<-function(formula,data=list(), index=NULL, endog = NULL, instruments= NULL, method = c("w2sls", "b2sls", "g2sls", "ec2sls"), lag = FALSE, listw, effects = NULL){

if(length(method) !=1 && effects == "fixed") method <- "w2sls" 	
if(length(method) !=1 && effects == "random") method <- "ec2sls" 	
		
  if(!is.null(index)) {
    require(plm)
    data <- plm.data(data, index)
    }
  
  index <- data[,1]
  tindex <- data[,2]

  ## record call
  cl <- match.call()

  ## check
  if(dim(data)[[1]]!=length(index)) stop("Non conformable arguments")
  
    mt <- terms(formula, data = data)
    mf <- lm(formula, data, na.action = na.fail, method = "model.frame")

    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
    
  names(index)<-row.names(data)
  ind<-index[which(names(index)%in%row.names(x))]
  tind<-tindex[which(names(index)%in%row.names(x))]
  
#if (lag)  {
#	 oo <- order(tind, ind)
#    xs <- x[oo, ]
#    ys <- y[oo]
#    inds <- ind[oo]
#    tinds <- tind[oo]
#
#	}
#  
  N<-length(unique(ind))
  k<-dim(x)[[2]]
  T<-max(tapply(x[,1],ind,length))
  NT<-length(ind)


  balanced<-N*T==NT
if(!balanced) stop("Estimation method unavailable for unbalanced panels")



if(!lag){
if(is.null(endog)) stop("No engogenous variables specified")
if(is.character(endog)){
	xend<-match(endog,colnames(data))
	endog <- data[,xend]
	}
kend <- ncol(endog) 
# print(kend)
	
if(is.null(instruments)) stop("No instruments specified")
if(is.character(instruments)){
	inst<-match(instruments,colnames(data)) 
	instruments <- data[,inst]
	}
kinst <- ncol(instruments) 

if(kinst < kend) stop("The model is not identified: Not enough instruments specified")
}




if(lag){
	
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

if(is.null(endog) && lag == FALSE) stop("No engogenous variables specified")

if(is.character(endog)){
	xend<- match(endog,colnames(data))  
	endog <- data[,xend]
if(is.null(instruments)) warnings("Spatially lagged exogenous variables will be used as instruments")
kend <- ncol(endog)  + 1
	}
else kend <-  1


kint<-ifelse(colnames(x)[1]=="(Intercept)", 1, 0)
#print(kint)
if(is.character(instruments)){
	inst<- match(instruments,colnames(data))
	instruments <- data[,inst]
kinst <- ncol(instruments)  + kint + ncol(as.matrix(x[,-kint]))*3
	}
	
else kinst <-  kint + ncol(as.matrix(x[,-kint]))*3
#print(kinst)
if(kinst < kend) stop("The model is not identified: Not enough instruments specified")
	
	}




switch(method, 
w2sls = {
	result <- ivplm.w2sls(Y = y,X =x, H = instruments, endog = endog, ind = ind, tind = tind, lag = lag, listw = listw)
	},
b2sls = {
	result <- ivplm.b2sls(Y = y,X =x, H = instruments, endog = endog, ind = ind, tind = tind, lag = lag, listw = listw)
	},
ec2sls = {
	result <- ivplm.ec2sls(Y = y,X =x, H = instruments, endog = endog, ind = ind, tind = tind, lag = lag, listw = listw)
	},
g2sls = {
	result <-ivplm.g2sls(Y = y,X =x, H = instruments, endog = endog, ind = ind, tind = tind, lag = lag, listw = listw)
	},
stop("...\nUnknown method\n"))


    result$zero.policy <- FALSE
    result$robust <- FALSE
    result$legacy <- FALSE
    result$listw_style <- FALSE
    result$call <- match.call()

class(result) <- "stsls"
result
}










 