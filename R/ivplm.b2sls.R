###### between 2sls

ivplm.b2sls <- function(Y,X,H = NULL, endog = NULL, 
                        twow, lag = FALSE, listw = NULL, listw2 = NULL,
                        lag.instruments = NULL, t = t, N = N, NT = NT, Durbin = FALSE, xdur){

indic <- rep(1:N,t)
xdu <- X

##transform y	
ybetween<-panel.transformations(Y,indic, type= "between")
ndim <- length(ybetween)
listwnn <- listw[1:ndim, 1:ndim]


if(isTRUE(Durbin)  | inherits(Durbin, "formula")){
  
  if(inherits(Durbin, "formula")){
    
    colnmx <- colnames(X)
    colnamesbx <- paste("lag_", colnames(xdur), sep="")
    wx <- listw %*% xdur
    X <- cbind(X, wx)
    Xbetween <- panel.transformations(X, indic, type= "between")
    colnames(Xbetween) <- c(colnmx, colnamesbx)
    
    if (colnames(Xbetween)[1] == "(Intercept)") Xbetween <- Xbetween[,-1]
    delb <- as.numeric(which(diag(var(Xbetween)) == 0))
    if(length(delb) == 0) Xbetween <- Xbetween
    else Xbetween <- Xbetween[,-delb]
    
    if (colnmx[1] == "(Intercept)") Xbetween <- cbind(1, Xbetween)
    colnames(Xbetween)[1] <- "(Intercept)"
    xdu <- cbind(xdu, wx)
    colnames(xdu) <- c(colnmx, colnamesbx)
    
  }
  else{
    
    colnmx <- colnames(X)
    
    if(colnmx[1] == "(Intercept)"){
    
        wx <- listw %*% X[,-1]
        colnameswx <- paste("lag_", colnames(X)[-1], sep = "")
        xdu <- cbind(xdu, wx)
        colnames(xdu) <- c(colnmx, colnameswx)
        
    }
    else{
    
      wx <- listw %*% X
      colnameswx <- paste("lag_", colnames(X), sep = "")
      xdu <- cbind(xdu, wx)
      colnames(xdu) <- c(colnmx, colnameswx)
      
    }
    
    X <- cbind(X, wx)
    
    Xbetween <- panel.transformations(X, indic, type = "between")
    colnames(Xbetween) <- c(colnmx, colnameswx) 
    
    if (colnames(Xbetween)[1] == "(Intercept)") Xbetween <- Xbetween[,-1]
    delb <- as.numeric(which(diag(var(Xbetween)) == 0))
    if(length(delb) == 0) Xbetween <- Xbetween
    else Xbetween <- Xbetween[,-delb]
    #print(colxbeet)
    if (colnmx[1] == "(Intercept)") Xbetween <- cbind(1, Xbetween)
    colnames(Xbetween)[1] <- "(Intercept)"
    
  }
}
else{
  
  colnmx <- colnames(X)
  Xbetween <- panel.transformations(X, indic, type = "between")
  colnames(Xbetween) <- colnmx
  
  if (colnames(Xbetween)[1] == "(Intercept)") Xbetween <- Xbetween[,-1]
  delb <- as.numeric(which(diag(var(Xbetween)) == 0))
  if(length(delb) == 0) Xbetween <- Xbetween
  else Xbetween <- Xbetween[,-delb]
  
  if (colnmx[1] == "(Intercept)") Xbetween <- cbind(1, Xbetween)
  colnames(Xbetween)[1] <- "(Intercept)"
#  xdu <- cbind(xdu, wx)
#  colnames(xdu) <- c(colnmx, colnameswx)
  
}



if(!lag){
##transform the instruments H and the endogenous variable
	Hbetween<-panel.transformations(H,indic, type= "between")

#transorm the endogenous 
	endogbetween<-panel.transformations(endog,indic, type= "between")
   colnames(endogbetween)<-colnames(endog)

res <-spgm.tsls(sqrt(t)*as.matrix(ybetween), sqrt(t)*endogbetween, sqrt(t)*Xbetween, sqrt(t)*as.matrix(Hbetween) )
res$Hbetween <- Hbetween
res$xdu <- xdu
res$type <- "b2sls model without spatial lag"

}

else{
	
	wybetween <- listwnn %*% as.matrix(ybetween)
	wybetween <- as.matrix(wybetween)
  colnames(wybetween) <- ("lambda")
	
	if(is.null(endog)){
		#no external instruments 
		            
	  if(twow){
		            	
     listw2nn     <- listw2[1:ndim, 1:ndim]       	
	   WXbetween    <- listwnn %*%  Xbetween
     WWXbetween   <- listwnn %*% WXbetween
	   W2Xbetween   <- listw2nn %*%  Xbetween
     W2WXbetween  <- listw2nn %*% WXbetween
     W2WWXbetween <- listw2nn %*% WWXbetween
            	
 	Hbetween <- cbind(as.matrix(WXbetween), as.matrix(WWXbetween), as.matrix(W2Xbetween), as.matrix(W2WXbetween), as.matrix(W2WWXbetween))            	
    
            }
else{            
	
            WXbetween  <- as.matrix(listwnn %*% Xbetween)
            WWXbetween <- as.matrix(listwnn %*% WXbetween)
        
  Hbetween <- cbind(WXbetween, WWXbetween)        
 	
 	}


res<-spgm.tsls(sqrt(t)*as.matrix(ybetween), sqrt(t)*as.matrix(wybetween), sqrt(t)*Xbetween, sqrt(t)*as.matrix(Hbetween) )
res$Hbetween <- Hbetween
res$xdu <- xdu
res$type <- "Spatial b2sls model"

		}
		
else{
	
		Hbetween <- panel.transformations(H,indic, type= "between") 
			

            if(twow){
    listw2nn <- listw2[1:ndim, 1:ndim]        	
	  WXbetween <- listwnn %*%  Xbetween
    WWXbetween <- listwnn %*% WXbetween
	  W2Xbetween <- listw2nn %*%  Xbetween
    W2WXbetween <- listw2nn %*% WXbetween
    W2WWXbetween <- listw2nn %*% WWXbetween
            	
 	Hbetween <-cbind(Hbetween, as.matrix(WXbetween), as.matrix(WWXbetween), as.matrix(W2Xbetween), as.matrix(W2WXbetween), as.matrix(W2WWXbetween))            	
    
            }
else{            
	
	WXbetween <- listwnn %*%  Xbetween
  WWXbetween <- listwnn %*% WXbetween
 	Hbetween <-cbind(Hbetween, as.matrix(WXbetween), as.matrix(WWXbetween))
 	
 	}


	##transform the endogenous variables endog
	endogbetween<-panel.transformations(endog,indic, type= "between")
	endogbetween<-cbind(endogbetween, wybetween)

colnames(endogbetween)<-c(colnames(endog), "lambda")


	
res<-spgm.tsls(sqrt(t)*as.matrix(ybetween), sqrt(t)*endogbetween, sqrt(t)*Xbetween, sqrt(t)*as.matrix(Hbetween) )
res$Hbetween <- Hbetween
res$xdu <- xdu
res$type <- "Spatial b2sls model with additional endogenous variables"
	}			
	}	

res
}




