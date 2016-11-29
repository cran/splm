spml <- function(formula, data, index=NULL, listw, listw2=listw, na.action,
                 model=c("within","random","pooling"),
                 effect=c("individual","time","twoways"),
                 lag=FALSE, spatial.error=c("b","kkp","none"),
                 ...) {

  ## wrapper function for all ML models

  ## record call
  cl <- match.call()

  ## check class(listw)
  checklw <- function(x) {
    
    if(!("listw" %in% class(x))) {
        x <- x
      if("matrix" %in% class(x)) {
        #require(spdep)
        x <- mat2listw(x)
      } 
      else {
        stop("'listw' has to be either a 'listw' or a 'matrix' object")
      }}
      # }
    return(x)
  }

  listw <- checklw(listw)
  listw2 <- checklw(listw2)

  ## dimensions check is moved downstream

  ##added by gpiras on November 25, 2015 for consistency with the test bsk
  ## removed by the_sculler on Jan 8 2016 because bsktest() never calls spml()

  #if(model == 'pooling' && spatial.error == 'b' && lag ==FALSE){
  #
  #	res <- spfeml(formula=formula, data=data, index=index,
  #                listw=listw, listw2=listw2, na.action,
  #                model = 'error', effects = "pooling",
  #                cl=cl, ...)
  #}
  #else{
  switch(match.arg(model), within={
  
    if(lag) {
        model <- switch(match.arg(spatial.error),
                        b="sarar",
                        kkp="sarar",
                        none="lag")
    } else {
        model <- switch(match.arg(spatial.error),
                        b="error",
                        kkp="error",
                        none="plm")
        if(model == "plm") stop("No spatial component, use plm instead")
        ## put call to plm() here, fetch results
        ## and suitably transform them for compliance
    }
    effects <- switch(match.arg(effect), individual="spfe",
                      time="tpfe", twoways="sptpfe")

    res <- spfeml(formula=formula, data=data, index=index,
                  listw=listw, listw2=listw2, na.action,
                  model=model, effects=effects,
                  cl=cl, ...)
  }, random={
    switch(match.arg(effect),
           time={stop("time random effects not implemented")},
           twoways={stop("twoway random effects not implemented")},
           individual={
             errors <- switch(match.arg(spatial.error),
                              b="semre", kkp="sem2re", none="re")})
    res <- spreml(formula=formula, data=data, index=index,
                  w=listw2mat(listw), w2=listw2mat(listw2),
                  lag=lag, errors=errors, cl=cl, ...)
  }, pooling={
           errors <- switch(match.arg(spatial.error),
                              b="sem", kkp="sem", none="ols")
    res <- spreml(formula=formula, data=data, index=index,
                  w=listw2mat(listw), w2=listw2mat(listw2),
                  lag=lag, errors=errors, cl=cl, ...)
         })

   #}
  return(res)
}

    
