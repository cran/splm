ivplm.b2sls<-function(Y,X,H,endog, ind, tind, lag=FALSE, listw){

##transform y	
ybetween<-panel.transformations(Y,ind, type= "between")

Xbetween<-panel.transformations(X,ind, type= "between")
colnames(Xbetween)<-colnames(X)

if (colnames(Xbetween)[1] == "(Intercept)") Xbetween<-Xbetween[,-1]
delb<-as.numeric(which(diag(var(Xbetween))==0))
if(length(delb)==0) Xbetween<-Xbetween
else Xbetween<-Xbetween[,-delb]

if (colnames(X)[1] == "(Intercept)") Xbetween<-cbind(1,Xbetween)
colnames(Xbetween)[1]<-"(Intercept)"



if(!lag){
##transform the instruments H
	Hbetween<-panel.transformations(H,ind, type= "between")
	endogbetween<-panel.transformations(endog,ind, type= "between")
   colnames(endogbetween)<-colnames(endog)
##tsls
res<-spgm.tsls(sqrt(length(unique(tind)))*as.matrix(ybetween), sqrt(length(unique(tind)))*endogbetween, sqrt(length(unique(tind)))*Xbetween, sqrt(length(unique(tind)))*as.matrix(Hbetween) )
#print(res$coefficients)
}


else{
	
	wybetween <- lag.listw(listw, as.matrix(ybetween))
   colnames(wybetween) <- ("lambda")

	
	if(is.null(endog)){

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
        
Hbetween<-cbind(WXbetween, WWXbetween)        

res<-spgm.tsls(sqrt(length(unique(tind)))*as.matrix(ybetween), sqrt(length(unique(tind)))*as.matrix(wybetween), sqrt(length(unique(tind)))*Xbetween, sqrt(length(unique(tind)))*as.matrix(Hbetween) )


		}
		
else{
	
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
        
	##transform the instruments H
	Hbetween<-panel.transformations(H,ind, type= "between") 

Hbetween<-cbind(Hbetween, WXbetween, WWXbetween)


##transform the endogenous variables endog
	endogbetween<-panel.transformations(endog,ind, type= "between")

	endogbetween<-cbind(endogbetween, wybetween)

if(is.null(colnames(endog))) colnames(endogbetween)<-c(rep("endog", (ncol(endogbetween)-1)), "lambda")
else 	colnames(endogbetween)<-c(colnames(endog), "lambda")


	
res<-spgm.tsls(sqrt(length(unique(tind)))*as.matrix(ybetween), sqrt(length(unique(tind)))*endogbetween, sqrt(length(unique(tind)))*Xbetween, sqrt(length(unique(tind)))*as.matrix(Hbetween) )

	

	}		
	
	}	

res
}



