ivplm.g2sls<-function(Y,X,H,endog, ind, tind, lag=FALSE, listw){

   N<-length(unique(ind))
  T<-length(unique(tind))
  indt<-rep(seq(1,N), T)
  tindt<-rep(seq(1,T), each = N)
  tord<-order(indt,tindt)

##transform y	
transy<-panel.transformations(Y,ind, type= "both")
ybetween<-transy[[2]]
ywithin<-transy[[1]]
ybetweennt<- rep(ybetween, each=length(unique(tind)))	

##transform X	
transx<-panel.transformations(X,ind, type= "both")
Xbetween<-transx[[2]]
Xwithin<-transx[[1]]
colnames(Xwithin)<-colnames(X)
colnames(Xbetween)<-colnames(X)
Xbetweennt<-matrix(,nrow(Xwithin), ncol(Xbetween))
for (i in 1:ncol(Xbetween)) Xbetweennt[,i]<-rep(Xbetween[,i],each= length(unique(tind)))
del<- which(diag(var(Xwithin))==0)
colnames(Xbetweennt)<-colnames(X)

if(!lag){

##transform the instruments H
transH<-panel.transformations(H,ind, type= "both")
Hbetween<-transH[[2]]
Hwithin<-transH[[1]]
Hbetweennt<-matrix(,nrow(Hwithin), ncol(Hbetween))
for (i in 1:ncol(Hbetween)) Hbetweennt[,i]<-rep(Hbetween[,i],each= length(unique(tind)))

##transform the endogenous variables endog
transendog<-panel.transformations(endog,ind, type= "both")
endogbetween<-transendog[[2]]
endogwithin<-transendog[[1]]
endogbetweennt<-matrix(,length(ind), ncol(endogbetween))
for (i in 1:ncol(endogbetween)) endogbetweennt[,i]<-rep(endogbetween[,i],each= length(unique(tind)))
colnames(endogbetweennt)<-colnames(endog)
colnames(endogwithin)<-colnames(endog)

#W2SLS
resw<-spgm.tsls(as.matrix(ywithin), endogwithin, Xwithin, Hwithin )

sigma2v1<-resw$sse / ((length(unique(ind)) * (length(unique(tind)) -1)) - ncol(as.matrix(Xwithin[,-del])) - ncol(endogwithin)) 

#B2SLS
resb<-spgm.tsls(sqrt(length(unique(tind)))*as.matrix(ybetween), sqrt(length(unique(tind)))*as.matrix(endogbetween), sqrt(length(unique(tind)))*Xbetween, sqrt(length(unique(tind)))*as.matrix(Hbetween) )

sigma21<-resb$sse /  resb$df

ystar<-ywithin/sqrt(sigma2v1) + ybetweennt/sqrt(sigma21)
xstar<-Xwithin/sqrt(sigma2v1) + Xbetweennt/sqrt(sigma21)
endogstar<-endogwithin/sqrt(sigma2v1) + endogbetweennt/sqrt(sigma21)


Hstar<-Hwithin/sqrt(sigma2v1) + Hbetweennt/sqrt(sigma21)

res<-spgm.tsls(ystar, endogstar, xstar, Hstar )
res$sigma1<-sigma21
res$sigmav<-sigma2v1

}

else{
     wy<-lag.listwpanel(listw, Y, tind)
     wy<-wy[tord]

	  wywithin <- lag.listwpanel(listw, ywithin, tind)
	  wywithin <-wywithin[tord]
     wywithin <- as.matrix(wywithin)
     colnames(wywithin)<-"lambda"
  	  wybetween <- lag.listw(listw, as.matrix(ybetween))
     colnames(wybetween) <- ("lambda")


        WXwithin <- matrix(nrow = nrow(Xwithin), ncol = ncol(Xwithin))
        WWXwithin <- matrix(nrow = nrow(Xwithin), ncol = ncol(Xwithin))

for (i in 1:ncol(Xwithin)) {
           wx<- lag.listwpanel(listw, Xwithin[,i], tind)
           wwx<- lag.listwpanel(listw, wx, tind)

            if (any(is.na(wx))) 
                stop("NAs in lagged independent variable")
            WXwithin[, i] <- wx
            WWXwithin[, i] <- wwx
        }

WXwithin <- WXwithin[tord,] 
WWXwithin <- WWXwithin[tord,]



        WXbetween <- matrix(nrow = nrow(Xbetween), ncol = ncol(Xbetween))
        WWXbetween <- matrix(nrow = nrow(Xbetween), ncol = ncol(WXbetween))
for (i in 1:ncol(Xbetween)) {
            wx <- lag.listw(listw,Xbetween[,i])
            wwx <- lag.listw(listw,wx)
            if (any(is.na(wx))) 
                stop("NAs in lagged independent variable")
            WXbetween[, i] <- wx
            WWXbetween[, i] <- wwx
        }

	
	if(is.null(endog)){

Hwithin<-cbind(WXwithin, WWXwithin)    
    
resw<-spgm.tsls(ywithin, wywithin, Xwithin, Hwithin)

sigma2v1<- resw$sse / ((length(unique(ind)) * (length(unique(tind)) -1)) - ncol(as.matrix(Xwithin[,-del])) - 1) 

        
Hbetween<-cbind(WXbetween, WWXbetween)        

resb<-spgm.tsls(sqrt(length(unique(tind)))*as.matrix(ybetween), sqrt(length(unique(tind)))*as.matrix(wybetween), sqrt(length(unique(tind)))*Xbetween, sqrt(length(unique(tind)))*as.matrix(Hbetween) )
sigma21<-resb$sse /  resb$df


ystar<-ywithin/sqrt(sigma2v1) + ybetweennt/sqrt(sigma21)
xstar<-Xwithin/sqrt(sigma2v1) + Xbetweennt/sqrt(sigma21)
endogstar<-wywithin/sqrt(sigma2v1) + rep(wybetween, each=length(unique(tind)))/sqrt(sigma21)
endogstar<-as.matrix(endogstar)
colnames(endogstar)<-"lambda"

Hbetweennt<-matrix(,length(ind), ncol(Hbetween))
for (i in 1:ncol(Hbetween)) Hbetweennt[,i]<-rep(Hbetween[,i],each= length(unique(tind)))

Hstar<-Hwithin/sqrt(sigma2v1) + Hbetweennt/sqrt(sigma21)

res <- spgm.tsls(ystar, endogstar, xstar, Hstar)

res$sigma1 <- sigma21
res$sigmav <- sigma2v1

}


else{

transH<-panel.transformations(H,ind, type= "both")
Hbetween<-transH[[2]]
Hwithin<-transH[[1]]
Hbetweennt<-matrix(,length(ind), ncol(Hbetween))
for (i in 1:ncol(Hbetween)) Hbetweennt[,i]<-rep(Hbetween[,i],each= length(unique(tind)))

Hwithin<-cbind(Hwithin, WXwithin, WWXwithin)


transendog<-panel.transformations(endog,ind, type= "both")
endogbetween<-transendog[[2]]
endogwithin<-transendog[[1]]
endogbetweennt<-matrix(,length(ind), ncol(endogbetween))
for (i in 1:ncol(endogbetween)) endogbetweennt[,i]<-rep(endogbetween[,i],each= length(unique(tind)))

endogwithin<-cbind(endogwithin, wywithin)
    
resw<-spgm.tsls(as.matrix(ywithin), as.matrix(endogwithin), Xwithin, Hwithin )
sigma2v1<-resw$sse / ((length(unique(ind)) * (length(unique(tind)) -1)) - ncol(as.matrix(Xwithin[,-del])) - ncol(endogwithin)) 


Hbetween<-cbind(Hbetween, WXbetween, WWXbetween)
endogbetween<-cbind(endogbetween, wybetween)

resb<-spgm.tsls(sqrt(length(unique(tind)))*as.matrix(ybetween), sqrt(length(unique(tind)))*as.matrix(endogbetween), sqrt(length(unique(tind)))*Xbetween, sqrt(length(unique(tind)))*as.matrix(Hbetween))
sigma21<-resb$sse / resb$df

ystar<-ywithin/sqrt(sigma2v1) + ybetweennt/sqrt(sigma21)
xstar<-Xwithin/sqrt(sigma2v1) + Xbetweennt/sqrt(sigma21)
endogstar<-endogwithin/sqrt(sigma2v1) + rep(endogbetween, each=length(unique(tind)))/sqrt(sigma21)

Hbetweennt<-matrix(,length(ind), ncol(Hbetween))
for (i in 1:ncol(Hbetween)) Hbetweennt[,i]<-rep(Hbetween[,i],each= length(unique(tind)))

Hstar<- Hwithin/sqrt(sigma2v1) + Hbetweennt/sqrt(sigma21)

res <- spgm.tsls(ystar, endogstar, xstar, Hstar)

res$sigma1<- sigma21
res$sigma1<- sigma2v1
	}
}


res
}


