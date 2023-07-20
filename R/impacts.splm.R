
impacts <- function(obj, ...){
  UseMethod("impacts", obj)
}

impacts.splm_ML <- function(obj, listw = NULL,
                         time = NULL, ...,
                         tr = NULL, R = 200,
                         type = "mult",
                         empirical = FALSE, Q = NULL){
  
  if(is.null(listw) && is.null(tr)) stop("either listw or tr should be provided")
  
  
  if(!is.null(listw) ){
    if(listw$style != "W") stop("Only row-standardised weights supported")
    if(is.null(time) && is.null(tr)) stop("time periods should be provided")
  }
  
  
  if(is.null(tr)){
    
    sparse.W <- listw2dgCMatrix(listw)
    s.lws <- kronecker(Diagonal(time) , sparse.W)
    tr <- trW(s.lws, type= type)
    
  }
  
  if(is.na(match(obj$type, c("fixed effects lag","fixed effects sarar",
                             "random effects ML", "fixed effects GM","lag GM",
                             "fixed effects GM")))) stop("object type not recognized")
  
  if(obj$type == "fixed effects lag"){
    
    class(obj)<- "Gmsar"
    obj$type <- "SARAR"
    obj$data <- as.vector(obj$model)
    obj$s2 <- obj$sigma2
    obj$secstep_var <- obj$vcov
    imp <- spatialreg::impacts(obj, tr=tr, R=R, ...)
    
  }
  
  if(obj$type == "fixed effects sarar"){
    
    class(obj)<- "Gmsar"
    obj$type <- "SARAR"
    rho <- obj$coefficients[2]
    obj$coefficients <- obj$coefficients[-2]
    obj$data <- as.vector(obj$model)
    obj$s2 <- obj$sigma2
    obj$secstep_var <- obj$vcov[-2,-2]
    imp <- spatialreg::impacts(obj, tr=tr, R=R,...)
    
  }
  
  if(obj$type == "fixed effects error") stop("Impacts Estimates are not available for Error Model")
  
  if(obj$type == "random effects ML")	{
    
    if(!is.null(obj$arcoef)) {
      class(obj)<- "Gmsar"
      obj$type <- "SARAR"
      
      obj$coefficients <- c(obj$arcoef, obj$coefficients)
      obj$data <- as.vector(obj$model)
      obj$s2 <- obj$sigma2
      obj$secstep_var <- matrix(0,nrow(obj$vcov)+1,nrow(obj$vcov)+1)
      obj$secstep_var[1,1] <- obj$vcov.arcoef
      obj$secstep_var[(2:(nrow(obj$vcov)+1)),(2:(nrow(obj$vcov)+1))] <- obj$vcov
      imp <- spatialreg::impacts(obj, tr=tr, R=R, ...)
    }
    else stop("Impacts Estimates are not available for Error Model")
    
  }
  
  
  if(obj$type == "fixed effects GM"){
    
    if(is.null(obj$endog)) {
      obj$secstep_var <- vcov(obj)
      class(obj)<- "Gmsar"
      obj$type <- "SARAR"
      obj$data <- as.vector(obj$model)
      obj$s2 <- obj$sigma2
      
      imp <- impacts(obj, tr=tr, R=R, ...)
      
      
    }
    
    else stop("No impacts estimates when endogenous variables are present in the system")
    
  }
  
  if(obj$type == "lag GM")			{
    
    if(is.null(obj$endog)) {
      
      class(obj)<- "Gmsar"
      obj$type <- "SARAR"
      obj$secstep_var <- obj$var
      obj$data <- as.vector(obj$model)
      obj$s2 <- obj$sigma2
      
      imp <- impacts(obj, tr=tr, R=R, ...)
      
      
    }
    
    else stop("No impacts estimates when endogenous variables are present in the system")
    
    
    
  }
  
  
  if(obj$type == "random effects GM")			{
    
    if(is.null(obj$endog)) {
      
      class(obj)<- "Gmsar"
      obj$type <- "SARAR"
      obj$secstep_var <- obj$vcov
      obj$data <- as.vector(obj$model)
      obj$s2 <- obj$sigma2
      
      imp <- impacts(obj, tr=tr, R=R, ...)
      
      
    }
    
    else stop("No impacts estimates when endogenous variables are present in the system")
    
    
    
  }
  
  
  
  
  
  return(imp)
  
}


impacts.splm_GM <- function(obj, ..., tr=NULL, 
                            R=NULL, listw=NULL,
                            type = "mult",
                            time = NULL,
                            evalues=NULL, tol=1e-6, 
                            empirical=FALSE, Q=NULL, 
                            KPformula = FALSE, prt = TRUE) {
  
  object <- obj
  
  
  if(is.null(listw) && is.null(tr)) stop("either listw or tr should be provided")
  if (!is.null(object$endo)) stop("impacts for model 
                                  with additional endogenous 
                                  variables not yet available in splm")
  
  if(!is.null(listw) ){	
    if(!inherits(listw,c("listw", "Matrix", "matrix"))) stop("listw format unknown")
    if(inherits(listw,"listw"))  Ws <- listw2dgCMatrix(listw)	
    if(inherits(listw,"matrix"))  Ws <- Matrix(listw)	
    if(inherits(listw,"Matrix"))  Ws <- listw	
    
    if(all(rowSums(Ws) != rep(1,Ws@Dim[1]))) stop("Only row-standardised weights supported")
    if(is.null(time) && is.null(tr)) stop("time periods should be provided")
  }
  
  
  
  if(is.null(tr)){
    
    s.lws <- kronecker(Diagonal(time), Ws)
    tr <- trW(s.lws, type= type)
    
  }
  
  
  if(object$nfimpacts == "lag_gm"){
    
    if(isTRUE(KPformula)){
      
      stop("KP formula not yet implemented")
      # if(is.null(tr)){
      #   if(is.null(evalues)) evalues <- eigen(object$listw)$values
      # }
      # if(isTRUE(object$Durbin) | inherits(object$Durbin, "formula")) vc_impacts_formula_lag_mixed(object, evalues, tr, prt)
      # else vc_impacts_formula_lag(object, evalues, tr, prt)
    }
    else{ 
      
      if(isTRUE(object$Durbin) | inherits(object$Durbin, "formula")){
        
        type <- "mixed"
        interval <- c(-1,0.999)
        
        coefs <- drop(object$coefficients)
       # print(coefs)
        p2 <- length(coefs)
        lambda <- coefs[1]
        #print(lambda)
        beta <- coefs[2:p2]
        #print(beta)
        p <- length(beta)
        p1 <- p + 1
        names(beta) <- names(object$coefficients)[2:p2]
        #print(names(beta))
        icept <- grep("(Intercept)", names(beta))
        iicept <- length(icept) > 0L
        n <- length(object$residuals)
        
        irho <- 1
        drop2beta <-  1
        
        if((lambda > interval[2] ) | (lambda < interval[1])) warning("Value of the spatial parameter outside of parameter space")    
        
            
          dvars <- NULL
          if (iicept) b1 <- beta[-icept]
          else  b1 <- beta
          p <- length(b1)
          #print(names(beta))
          dn <- grep("lag_", names(beta)) #which of the names of beta has "lag_"
          #print(dn)
          dc <- beta[dn] # betas that are lagged
          #print(dc)
          beta1 <- beta[-dn] # all the betas (and)including the lagged) 
          #print(beta1)
          xb <- beta1[which(names(beta1) %in% stringr::str_remove(names(dc),"lag_") )]
          #print(order(names(xb)))
          xb <- xb[order(names(xb))]
          l_xb <- length(xb)
          if(length(which(names(beta1) %in% stringr::str_remove(names(dc),"lag_") ))>=1) xo <- beta1[ -which(names(beta1) %in% stringr::str_remove(names(dc),"lag_") )]
          else xo <- beta1
          
          
          gamma <- dc[which( stringr::str_remove(names(dc),"lag_")  %in% names(beta1))]
          gamma <- gamma[order(names(gamma))]
          l_gamma <- length(gamma)
          
          if(length(which(stringr::str_remove(names(dc),"lag_") %in% names(beta1)))>=1) don <- dc[-which( stringr::str_remove(names(dc),"lag_")  %in% names(beta1))]
          else don <- dc
          l_don <- length(don)
          
          
          
          if (iicept) {
            xo <- xo[-icept]
            l_xo <- length(xo)  
            if(l_xb != 0){
              b2 <- c(xo, xb, rep(0, l_don), rep(0, l_xo), gamma, don)
              bnames <- c(names(xo), sort(names(xb)), names(don)) 
            } 
            else{
              b2 <- c(xo, rep(0, l_don), rep(0, l_xo),  don)
              bnames <- c(names(xo), names(don)) 
            }
            P <- matrix(b2, ncol = 2)
          } 
          else {
            xo <- xo
            l_xo <- length(xo)
            if(l_xb != 0){
              b2 <- c(xo, xb, rep(0, l_don), rep(0, l_xo), gamma, don)
              bnames <- c(names(xo), sort(names(xb)), names(don)) 
            } 
            else{
              b2 <- c(xo, rep(0, l_don), rep(0, l_xo),  don)
              bnames <- c(names(xo), names(don)) 
            }
            P <- matrix(b2, ncol = 2)
          }
          #print(P)
          
          if(!is.null(R)){
            
            #first build the zero fills
            if (l_xo != 0 && l_don != 0 && l_xb !=0) {
              zero_fill <- list(1: l_xo, (l_xo+1):(l_xo+ l_xb),  l_don, 
                                l_xo, (l_xo + l_xb+1):  (l_xo + l_xb+l_gamma),
                                (l_xo + l_xb+l_gamma+1):(l_xo + l_xb+l_gamma+l_don))
              l_zero_fill <- 6
            }
            
            if(l_xo != 0 && l_don == 0 && l_xb != 0) {
              zero_fill <- list(1: l_xo, (l_xo+1):(l_xo+ l_xb),  
                                l_xo, (l_xo + l_xb+1):  (l_xo + l_xb+l_gamma))
              
              l_zero_fill <- 5
            }
            if(l_xo == 0 && l_don != 0 && l_xb !=0)  {
              zero_fill <- list(1: l_xb,  l_don,(l_xb+1):  (l_xb+l_gamma), (l_xb+l_gamma+1):( l_xb+l_gamma+l_don))
              l_zero_fill <- 4
            }
            
            if(l_xo != 0 && l_don != 0 && l_xb ==0)  {
              zero_fill <- list(1: l_xo, l_don, 
                                l_xo, (l_xo+1):(l_xo +l_don))
              l_zero_fill <- 3
            }
            
            # if(l_xo == 0 && l_don != 0 && l_xb ==0){
            #   zero_fill <- list(l_don, 1: l_don)
            #   l_zero_fill <- 2   
            #   
            # }
            
            # print(class(zero_fill))
            attr(zero_fill, "l_zero_fill") <- l_zero_fill
            #   print(l_zero_fill)
            
            if(iicept) {
              mu <- c(lambda, beta[1], xo, xb, gamma, don)
              m0 <- c(xo, xb, gamma, don)
            }
            else{ 
              mu <- c(lambda, xo, xb, gamma, don)
              m0 <- c(xo, xb, gamma, don)
              }
            
           # print(mu)
            #print(object$vcov)
            Sigma <- object$vcov[match(names(mu), rownames(object$vcov)),match(names(mu), rownames(object$vcov)) ]
            
            #print(Sigma)
          }
          ####CHECK THE 
          #  if (!requireNamespace("spatialreg", quietly=TRUE))
          #  stop("install spatialreg")
          # res <- spatialreg::intImpacts(rho=lambda, beta=beta, P=P, n=n, mu=mu,
          #                               Sigma=Sigma, irho=irho, drop2beta=drop2beta, bnames=bnames,
          #                               interval=interval, type = type, tr=tr, R=R, listw=listw, evalues=evalues,
          #                               tol=tol, empirical=empirical, Q=Q, icept=icept, iicept=iicept, p=p, 
          #                               zero_fill = zero_fill)
          
          res <- spatialreg::intImpacts(rho=lambda, beta=beta, P=P, n=n, mu=mu,
                                        Sigma=Sigma, irho=irho, drop2beta=drop2beta, bnames=bnames,
                                        interval=interval, type = type, tr=tr, R=R, listw=NULL,
                                        evalues=NULL,
                                        tol=tol, empirical=empirical, Q=Q, icept=icept, iicept=iicept, p=p, 
                                        zero_fill = zero_fill)
        }
    
      
      else{
        
        type <- "lag" 
        if (is.null(object$interval)) interval <- c(-1,0.999)
        
        coefs <- drop(object$coefficients)
        p2 <- length(coefs)
        lambda <- coefs[1]
        beta <- coefs[2:p2]
        names(beta) <- names(object$coefficients)[2:p2]
        colnames(object$vcov) <- rownames(object$vcov) <- names(object$coefficients)
        
        p <- length(beta)
        p1 <- p + 1
        
        if((lambda > interval[2] ) | (lambda < interval[1])) warning("Value of the spatial parameter outside of parameter space")    
        
        icept <- grep("(Intercept)", names(beta))
        iicept <- length(icept) > 0L
        zero_fill <- NULL
        dvars <- NULL
        
        if (iicept) {
          P <- matrix(beta[-icept], ncol=1)
          bnames <- names(beta[-icept])
        } else {
          P <- matrix(beta, ncol=1)
          bnames <- names(beta)
        }
        
        n <- length(object$residuals)
        mu <- NULL
        Sigma <- NULL
        
        if (!is.null(R)) {
          
          mu <- c(lambda, beta)
          Sigma <- object$vcov[match(names(mu), rownames(object$vcov)),match(names(mu), rownames(object$vcov))]
          
        }      

        irho <- length(beta) + 1
        drop2beta <- length(beta) + 1
        
        #  if (!requireNamespace("spatialreg", quietly=TRUE))
        #  stop("install spatialreg")
        # res <- spatialreg::intImpacts(rho=lambda, beta=beta, P=P, n=n, mu=mu,
        #                               Sigma=Sigma, irho=irho, drop2beta=drop2beta, bnames=bnames,
        #                               interval=interval, type = type, tr=tr, R=R, listw=listw, evalues=evalues,
        #                               tol=tol, empirical=empirical, Q=Q, icept=icept, iicept=iicept, p=p)
        # 
    
        res <- spatialreg::intImpacts(rho=lambda, beta=beta, P=P, n=n, mu=mu,
                                      Sigma=Sigma, irho=irho, drop2beta=drop2beta, bnames=bnames,
                                      interval=interval, type = type, tr=tr, R=R, listw=NULL,
                                      evalues=NULL,
                                      tol=tol, empirical=empirical, Q=Q, icept=icept, iicept=iicept, p=p, 
                                      zero_fill = zero_fill)    
      }
      
      # if (!requireNamespace("spatialreg", quietly=TRUE))
      #   stop("install spatialreg")
    
      
      attr(res, "iClass") <- class(object)
      
      res
    }
  }
  
  else stop("impacts for sarar not yet implemented")
  
  
}



