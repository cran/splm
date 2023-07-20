###### within 2sls
ivplm.w2sls <- function(Y,X,H = NULL, endog = NULL, 
                        twow, lag = FALSE, listw = NULL,  listw2 = NULL, 
                        lag.instruments = NULL, t, N, NT, Durbin = FALSE, xdur){
indic <- rep(1:N,t)

xdu <- X
#print(head(xdu))
##transform y and X	

ywithin <- panel.transformations(Y, indic, type = "within")

#transorm x incuding durbin
if(isTRUE(Durbin)  | inherits(Durbin, "formula")){
  
  if(inherits(Durbin, "formula")){
    
    colnmx <- colnames(X)
    colnameswx <- paste("lag_", colnames(xdur), sep="")
    wx <- listw %*% xdur
    colnames(wx) <- colnameswx
    X <- cbind(X, wx)
    Xwithin <- panel.transformations(X, indic, type = "within")
    colnames(Xwithin) <- c(colnmx, colnameswx)
    del <- which(diag(var(Xwithin)) == 0)
    Xwithin <- Xwithin[,-del]
    xdu <- cbind(xdu, wx)
    colnames(xdu) <- c(colnmx, colnameswx)
    #print(xdu)
  }
  else{
    
   colnmx <- colnames(X)
   wx <- listw %*% X
   colnameswx <- paste("lag_", colnames(X), sep="")
   X <- cbind(X, wx)
   Xwithin <- panel.transformations(X, indic, type = "within")
   colnames(Xwithin) <- c(colnmx, colnameswx)
   del <- which(diag(var(Xwithin)) == 0)
   Xwithin <- Xwithin[,-del]
   xdu <- cbind(xdu, wx)
   colnames(xdu) <- c(colnmx, colnameswx)
  }
}
else{
  xdu <- xdu
  colnmx <- colnames(X)
  Xwithin <- panel.transformations(X, indic, type = "within")
  colnames(Xwithin) <- colnmx
  del <- which(diag(var(Xwithin)) == 0)
  Xwithin <- Xwithin[,-del]
  colnameswx <- NULL
  colnames(xdu) <- c(colnmx)
}
	#print(head(xdu))
#xdu <- X[, -del]
#colnames(xdu) <- c(colnmx, colnameswx)[-del]

if(!lag){
##transform the instruments H

Hwithin <- panel.transformations(H, indic, type= "within")

##transform the endogenous variables endog

endogwithin <-panel.transformations(endog, indic, type= "within")
colnames(endogwithin) <- colnames(endog)



res<-spgm.tsls(as.matrix(ywithin), as.matrix(endogwithin), Xwithin, as.matrix(Hwithin))
#print(names(res$coefficients))
#print(colnames(Xwithin))
varb<-res$var *res$df /((N * (t -1)) - ncol(as.matrix(Xwithin)) - ncol(endogwithin)) 
res$var<-varb
sigma2v1<- res$sse/ ((N * (t -1)) - ncol(as.matrix(Xwithin)) - ncol(endogwithin)) 
res$sigmav<- sigma2v1	
res$Hwithin <- Hwithin
res$Xwithin <- Xwithin
res$del <- del
#print(names(res$coefficients))
#print(head(xdu[,which(colnames(xdu) %in% names(res$coefficients))]))
res$xdu <- xdu
#print(head(res$xdu))
res$type <- "w2sls model without spatial lag"

	}
	
else{

   wywithin <- listw %*% as.matrix(ywithin)
   wywithin <- as.matrix(wywithin)
   colnames(wywithin)<-"lambda"
   

   
if(is.null(endog)){
#no external instruments
            
      if(twow){

   	WXwithin <- listw %*%  Xwithin
    WWXwithin <- listw %*% WXwithin
	  W2Xwithin <- listw2 %*%  Xwithin
    W2WXwithin <- listw2 %*% WXwithin
    W2WWXwithin <- listw2 %*% WWXwithin

 	Hwithin <-cbind(as.matrix(WXwithin), as.matrix(WWXwithin), as.matrix(W2Xwithin), as.matrix(W2WXwithin), as.matrix(W2WWXwithin))

            }
else{

	WXwithin <- listw %*%  Xwithin
  WWXwithin <- listw %*% WXwithin
 	Hwithin <-cbind(as.matrix(WXwithin), as.matrix(WWXwithin))

 	}


res<-spgm.tsls(ywithin, wywithin, Xwithin, Hwithin)
varb<-res$var *res$df / ((N * (t -1)) - ncol(as.matrix(Xwithin)) - 1) 
res$var<-varb
sigma2v1<- res$sse / ((N * (t -1)) - ncol(as.matrix(Xwithin)) - 1) 
res$sigmav <- sigma2v1
res$del <- del
res$xdu <- xdu
#print(head(res$xdu))
res$Hwithin <- Hwithin
res$Xwithin <- Xwithin
res$type <- "Spatial w2sls model"
		}
		
else{
  
			#add the external instruments
			Hwithin <-panel.transformations(H, indic, type= "within")

            if(twow){
            	
	  WXwithin <- listw %*%  Xwithin
    WWXwithin <- listw %*% WXwithin
	  W2Xwithin <- listw2 %*%  Xwithin
    W2WXwithin <- listw2 %*% WXwithin
    W2WWXwithin <- listw2 %*% WWXwithin
            	
 	  Hwithin <-cbind(Hwithin, as.matrix(WXwithin), as.matrix(WWXwithin), as.matrix(W2Xwithin), as.matrix(W2WXwithin), as.matrix(W2WWXwithin))            	
    
            }
else{            
	
	  WXwithin <- listw %*%  Xwithin
    WWXwithin <- listw %*% WXwithin
 	  Hwithin <-cbind(Hwithin, as.matrix(WXwithin), as.matrix(WWXwithin))
 	
 	}

endogwithin <- panel.transformations(endog, indic, type= "within")

endogwithin <-cbind(endogwithin, wywithin)
colnames(endogwithin)<-c(colnames(endog), "lambda")
# colnames(Xwithin)<-colnames(X)[-del]

 

res <- spgm.tsls(ywithin, endogwithin, Xwithin, Hwithin)

res$xdu <- c(endog,  xdu[,which(colnames(xdu) %in% names(res$coefficients))])
varb<-res$var *res$df / ((N * (t -1)) - ncol(as.matrix(Xwithin)) - ncol(endogwithin)) 
res$var<-varb
sigma2v1<- res$sse / ((N * (t -1)) - ncol(as.matrix(Xwithin)) - ncol(endogwithin)) 
res$sigmav <- sigma2v1
res$Hwithin <- Hwithin
res$Xwithin <- Xwithin
res$del <- del 
res$xdu <- xdu
#print(head(res$xdu))
res$type <- "Spatial w2sls model with additional endogenous variables"
	}		
	
}	
res
}




