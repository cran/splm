ivplm.w2sls<-function(Y,X,H,endog, ind, tind, lag=FALSE, listw){

  N<-length(unique(ind))
  T<-length(unique(tind))
  indt<-rep(seq(1,N), T)
  tindt<-rep(seq(1,T), each = N)
  tord<-order(indt,tindt)


##transform y	
transy<-panel.transformations(Y,ind, type= "both")
ybetween<-transy[[2]]
ywithin<-transy[[1]]
	

##transform X	
transx<-panel.transformations(X,ind, type= "both")
Xbetween<-transx[[2]]
Xwithin<-transx[[1]]
colnames(Xwithin)<-colnames(X)
colnames(Xbetween)<-colnames(X)
del<- which(diag(var(Xwithin))==0)



if(!lag){

##transform the instruments H
transH<-panel.transformations(H,ind, type= "both")
Hbetween<-transH[[2]]
Hwithin<-transH[[1]]

##transform the endogenous variables endog
transendog<-panel.transformations(endog,ind, type= "both")
endogbetween<-transendog[[2]]
endogwithin<-transendog[[1]]
colnames(endogwithin)<-colnames(endog)
colnames(endogbetween)<-colnames(endog)



res<-spgm.tsls(as.matrix(ywithin), as.matrix(endogwithin), Xwithin, as.matrix(Hwithin) )
#print(res$coefficients)
varb<-res$var *res$df /((length(unique(ind)) * (length(unique(tind)) -1)) - ncol(as.matrix(Xwithin[,-del])) - ncol(endogwithin)) 

res$var<-varb

sigma2v1<- res$sse / ((length(unique(ind)) * (length(unique(tind)) -1)) - ncol(as.matrix(Xwithin[,-del])) - ncol(endogwithin)) 

res$sigmav<-	sigma2v1	


res
	}
	
else{
	
	wywithin <- lag.listwpanel(listw, ywithin, tind)
	wywithin <-wywithin[tord]
   wywithin <- as.matrix(wywithin)
   colnames(wywithin)<-"lambda"

	
if(is.null(endog)){



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
        
Hwithin<-cbind(WXwithin, WWXwithin)        
Hwithin<-Hwithin[tord,]

    
res<-spgm.tsls(ywithin, wywithin, Xwithin, Hwithin )


varb<-res$var *res$df / ((length(unique(ind)) * (length(unique(tind)) -1)) - ncol(as.matrix(Xwithin[,-del])) - 1) 
res$var<-varb

sigma2v1<- res$sse / ((length(unique(ind)) * (length(unique(tind)) -1)) - ncol(as.matrix(Xwithin[,-del])) - 1) 
res$sigmav<-	sigma2v1

		}
		
else{
	

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
	
WXwithin<-WXwithin[tord,]
WWXwithin<-WWXwithin[tord,]

##transform the instruments H
transH<-panel.transformations(H,ind, type= "both")
Hbetween<-transH[[2]]
Hwithin<-transH[[1]]


Hwithin<-cbind(Hwithin, WXwithin, WWXwithin)

##transform the endogenous variables endog
transendog<-panel.transformations(endog,ind, type= "both")
endogbetween<-transendog[[2]]
endogwithin<-transendog[[1]]

colnames(endogbetween)<-colnames(endog)


endogwithin<-cbind(endogwithin, wywithin)

if(is.null(colnames(endog))) colnames(endogwithin)<-c(rep("endog", (ncol(endogwithin)-1)), "lambda")
else 	colnames(endogwithin)<-c(colnames(endog), "lambda")

colnames(Xwithin)<-colnames(X)


res<-spgm.tsls(ywithin, endogwithin, Xwithin, Hwithin)


varb<-res$var *res$df / ((length(unique(ind)) * (length(unique(tind)) -1)) - ncol(as.matrix(Xwithin[,-del])) - ncol(endogwithin)) 
res$var<-varb

sigma2v1<- res$sse / ((length(unique(ind)) * (length(unique(tind)) -1)) - ncol(as.matrix(Xwithin[,-del])) - ncol(endogwithin)) 
res$sigmav<-	sigma2v1

	}		
	
	}	


res

}
