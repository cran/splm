
# # # # # ec 2sls

ivplm.ec2sls <- function(Y,X,H = NULL, endog = NULL, 
                         twow, lag = FALSE, listw = NULL,  listw2 = NULL, 
                         lag.instruments = NULL, t, N, NT, Durbin = FALSE, xdur){

indic <- rep(1:N,t)
listwnn <- listw[1:N, 1:N]
if(twow)  listw2nn <- listw2[1:N,1:N]

##transform y	
transy<-panel.transformations(Y,indic, type= "both")
ybetween<-transy[[2]]
ywithin<-transy[[1]]
ybetweennt<- rep(ybetween, t)	


#check if the model is durbin and include the spatial lags

if(isTRUE(Durbin)  | inherits(Durbin, "formula")){
  
  if(inherits(Durbin, "formula")){
    
    colnmx <- colnames(X)
    colnameswx <- paste("lag_", colnames(xdur), sep="")
    wx <- listw %*% xdur
    X <- cbind(X, wx)
    
    transx   <- panel.transformations(X,indic, type= "both")
    
    Xbetween <- transx[[2]]
    Xwithin  <- transx[[1]]
    
    
    colnames(Xwithin) <- c(colnmx, colnameswx)
    colnames(Xbetween) <- c(colnmx, colnameswx)
    
    Xbetweennt<-matrix(,NT, ncol(Xbetween))
    for (i in 1:ncol(Xbetween)) Xbetweennt[,i]<-rep(Xbetween[,i], t)
    del <- which(diag(var(Xwithin))==0)
    colnames(Xbetweennt)<- c(colnmx, colnameswx)
    
    
  }
  else{
    
    colnmx <- colnames(X)
    
    if(colnames(X)[1] == "(Intercept)"){
      wx <- listw %*% X[,-1]  
      colnameswx <- paste("lag_", colnames(X)[-1], sep="")
    }   
    else {
      wx <- listw %*% X
      colnameswx <- paste("lag_", colnames(X), sep="")
    }
    
    
    X <- cbind(X, wx)
    
    transx   <- panel.transformations(X,indic, type= "both")
    
    Xbetween <- transx[[2]]
    Xwithin  <- transx[[1]]
    
    colnames(Xwithin)  <- c(colnmx, colnameswx)
    colnames(Xbetween) <- c(colnmx, colnameswx)
    Xbetweennt<-matrix(,NT, ncol(Xbetween))
    for (i in 1:ncol(Xbetween)) Xbetweennt[,i]<-rep(Xbetween[,i], t)
    del <- which(diag(var(Xwithin))==0)
    colnames(Xbetweennt) <- c(colnmx, colnameswx)
    
    
  }
}
else{
  
  transx<-panel.transformations(X,indic, type= "both")
  Xbetween<-transx[[2]]
  Xwithin<-transx[[1]]
  colnames(Xwithin)<-colnames(X)
  colnames(Xbetween)<-colnames(X)
  Xbetweennt<-matrix(,NT, ncol(Xbetween))
  for (i in 1:ncol(Xbetween)) Xbetweennt[,i]<-rep(Xbetween[,i], t)
  del <- which(diag(var(Xwithin))==0)
  colnames(Xbetweennt) <- colnames(X)
  
  
  
}

if(!lag){
##transform the instruments H
transH     <- panel.transformations(H,indic, type= "both")
Hbetween   <- transH[[2]]
Hwithin    <- transH[[1]]
Hbetweennt <- matrix(,NT, ncol(Hbetween))
for (i in 1:ncol(Hbetween)) Hbetweennt[,i] <- rep(Hbetween[,i],t)


##transform the endogenous variables endog
transendog                <- panel.transformations(endog,indic, type= "both")
endogbetween              <- transendog[[2]]
endogwithin               <- transendog[[1]]
endogbetweennt            <- matrix(,NT, ncol(endogbetween))
for (i in 1:ncol(endogbetween)) endogbetweennt[,i]<-rep(endogbetween[,i],t)
colnames(endogbetweennt)  <- colnames(endog)
colnames(endogwithin)     <- colnames(endog)

#W2SLS
resw<-spgm.tsls(as.matrix(ywithin), endogwithin, Xwithin, Hwithin )

sigma2v1<-resw$sse / ((N * (t -1)) - ncol(as.matrix(Xwithin[,-del])) - ncol(endogwithin)) 

#B2SLS
resb<-spgm.tsls(sqrt(t)*as.matrix(ybetween), sqrt(t)*as.matrix(endogbetween), sqrt(t)*Xbetween, sqrt(t)*as.matrix(Hbetween) )

sigma21<-resb$sse /  resb$df

ystar<-ywithin/sqrt(sigma2v1) + ybetweennt/sqrt(sigma21)
xstar<-Xwithin/sqrt(sigma2v1) + Xbetweennt/sqrt(sigma21)
endogstar<-endogwithin/sqrt(sigma2v1) + endogbetweennt/sqrt(sigma21)

Hins <- cbind(Xwithin,Xbetweennt,Hwithin,Hbetweennt)
res  <- spgm.tsls(ystar, endogstar, xstar, Hinst = Hins, instr = TRUE)

res$sigma1<-sigma21
res$sigmav<-sigma2v1
res$type <- "ec2sls model without spatial lag"
}


else{
	
     wy        <- listw %*%  Y
     wywithin  <- listw %*% ywithin
     wywithin  <- as.matrix(wywithin)
     colnames(wywithin)<-"lambda"
  	 wybetween <- listwnn %*% as.matrix(ybetween)
     colnames(wybetween) <- ("lambda")
           
           #WXwithin <- as.matrix(listw %*% Xwithin)
           #WWXwithin <- as.matrix(listw %*%  WXwithin)

            #WXbetween <- as.matrix(listwnn %*% Xbetween)
            #WWXbetween <- as.matrix(listwnn %*% WXbetween)

if(is.null(endog)){

  if(twow){
    
    WXwithin <- listw %*%  Xwithin
    WWXwithin <- listw %*% WXwithin
    W2Xwithin <- listw2 %*%  Xwithin
    W2WXwithin <- listw2 %*% WXwithin
    W2WWXwithin <- listw2 %*% WWXwithin
    
    WXbetween <- listwnn %*%  Xbetween
    WWXbetween <- listwnn %*% WXbetween
    W2Xbetween <- listw2nn %*%  Xbetween
    W2WXbetween <- listw2nn %*% WXbetween
    W2WWXbetween <- listw2nn %*% WWXbetween
    
    Hwithin <-cbind(as.matrix(WXwithin), as.matrix(WWXwithin), as.matrix(W2Xwithin), as.matrix(W2WXwithin), as.matrix(W2WWXwithin))
    Hbetween <-cbind(as.matrix(WXbetween), as.matrix(WWXbetween), as.matrix(W2Xbetween), as.matrix(W2WXbetween), as.matrix(W2WWXbetween))
    
  }
  else{
    
    WXwithin <- listw %*%  Xwithin
    WWXwithin <- listw %*% WXwithin
    Hwithin <-cbind(as.matrix(WXwithin), as.matrix(WWXwithin))
    
    WXbetween <- listwnn %*%  Xbetween
    WWXbetween <- listwnn %*% WXbetween
    Hbetween <-cbind(as.matrix(WXbetween), as.matrix(WWXbetween))
    
  }
  

  
       
resw<-spgm.tsls(ywithin, wywithin, Xwithin, Hwithin)
sigma2v1<- resw$sse / ((N * (t -1)) - ncol(as.matrix(Xwithin[,-del])) - 1) 


resb<-spgm.tsls(sqrt(t)*as.matrix(ybetween), sqrt(t)*as.matrix(wybetween), sqrt(t)*Xbetween, sqrt(t)*as.matrix(Hbetween) )
sigma21<-resb$sse /  resb$df

ystar<-ywithin/sqrt(sigma2v1) + ybetweennt/sqrt(sigma21)
xstar<-Xwithin/sqrt(sigma2v1) + Xbetweennt/sqrt(sigma21)
endogstar<-wywithin/sqrt(sigma2v1) + rep(as.matrix(wybetween), t)/sqrt(as.numeric(sigma21))
endogstar<-as.matrix(endogstar)
colnames(endogstar)<-"lambda"


Hbetweennt<-matrix(,NT, ncol(Hbetween))
for (i in 1:ncol(Hbetween)) Hbetweennt[,i]<-rep(Hbetween[,i],t)

A <- cbind(1, Xwithin, Xbetweennt, Hwithin, Hbetweennt)

res <- spgm.tsls(ystar, endogstar, xstar, Hinst = A, instr = TRUE)
res$sigma1 <- sigma21
res$sigmav <- sigma2v1
res$type <- "Spatial ec2sls model"
}	

else{

transH <- panel.transformations(H,indic, type= "both")
Hbetween <- transH[[2]]
Hwithin<-transH[[1]]
Hbetweennt<-matrix(,NT, ncol(Hbetween))
for (i in 1:ncol(Hbetween)) Hbetweennt[,i]<-rep(Hbetween[,i], t)


#Hwithin<-cbind(Hwithin, WXwithin, WWXwithin)

if(twow){
  
  WXwithin <- listw %*%  Xwithin
  WWXwithin <- listw %*% WXwithin
  W2Xwithin <- listw2 %*%  Xwithin
  W2WXwithin <- listw2 %*% WXwithin
  W2WWXwithin <- listw2 %*% WWXwithin
  
  WXbetween <- listwnn %*%  Xbetween
  WWXbetween <- listwnn %*% WXbetween
  W2Xbetween <- listw2nn %*%  Xbetween
  W2WXbetween <- listw2nn %*% WXbetween
  W2WWXbetween <- listw2nn %*% WWXbetween
  
  Hwithin <-cbind(Hwithin, as.matrix(WXwithin), as.matrix(WWXwithin), as.matrix(W2Xwithin), as.matrix(W2WXwithin), as.matrix(W2WWXwithin))
  Hbetween <-cbind(Hbetween, as.matrix(WXbetween), as.matrix(WWXbetween), as.matrix(W2Xbetween), as.matrix(W2WXbetween), as.matrix(W2WWXbetween))
  
}
else{
  
  WXwithin <- listw %*%  Xwithin
  WWXwithin <- listw %*% WXwithin
  Hwithin <-cbind(Hwithin, as.matrix(WXwithin), as.matrix(WWXwithin))
  
  WXbetween <- listwnn %*%  Xbetween
  WWXbetween <- listwnn %*% WXbetween
  Hbetween <-cbind(Hbetween, as.matrix(WXbetween), as.matrix(WWXbetween))
  
}





transendog<-panel.transformations(endog,indic, type= "both")
endogbetween<-transendog[[2]]
endogwithin<-transendog[[1]]

# endogbetweennt<-matrix(,NT, ncol(endogbetween))
# for (i in 1:ncol(endogbetween)) endogbetweennt[,i]<-rep(endogbetween[,i], T)

endogwithin<-cbind(endogwithin, wywithin)
colnames(endogwithin) <- c(colnames(endog), "lambda")
    
resw<-spgm.tsls(as.matrix(ywithin), as.matrix(endogwithin), Xwithin, Hwithin )
sigma2v1<-resw$sse / ((N * (t -1)) - ncol(as.matrix(Xwithin[,-del])) - ncol(endogwithin)) 


Hbetween <- cbind(Hbetween, WXbetween, WWXbetween)
endogbetween <- cbind(endogbetween, as.matrix(wybetween))
colnames(endogbetween) <- c(colnames(endog), "lambda")
endogbetweennt<-matrix(,NT, ncol(endogbetween))
for (i in 1:ncol(endogbetween)) endogbetweennt[,i]<-rep(endogbetween[,i], t)


resb<-spgm.tsls(sqrt(t)*as.matrix(ybetween),  sqrt(t)*as.matrix(endogbetween), sqrt(t)*Xbetween, sqrt(t)*as.matrix(Hbetween))
sigma21<-resb$sse / resb$df


ystar<-ywithin/sqrt(sigma2v1) + ybetweennt/sqrt(sigma21)
xstar<-Xwithin/sqrt(sigma2v1) + Xbetweennt/sqrt(sigma21)
# print(dim(endogwithin))
# print(dim(as.matrix(endogbetween)))
endogstar<-endogwithin/sqrt(sigma2v1) + as.matrix(endogbetweennt)/sqrt(as.numeric(sigma21))

# print(sigma2v1)
# print(sigma21)

Hbetweennt<-matrix(,NT, ncol(Hbetween))
for (i in 1:ncol(Hbetween)) Hbetweennt[,i]<-rep(Hbetween[,i], t)

A<-cbind(1,Xwithin[,-del],Xbetweennt[,-del], Hwithin, Hbetweennt)

res <- spgm.tsls(ystar, endogstar, xstar, Hinst = A, instr = TRUE)
# print(res$coefficients)

res$sigma1<- sigma21
res$sigma1<- sigma2v1
res$type <- "Spatial ec2sls model with additional endogenous variables"

	}	
	
	}

res 

}
 


