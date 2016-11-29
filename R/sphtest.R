sphtest <- function (x, ...)
{
    UseMethod("sphtest")
}

sphtest.formula <- function (x, data, index = NULL, listw,
                             spatial.model = c("lag", "error", "sarar"),
                             method = c("ML", "GM"), errors = c("KKP", "BSK"),...) {
    ## performs a Hausman test of a RE model with spatial lag or error
    ## against FE "alternative" with same spatial specification

    switch(match.arg(spatial.model),
    lag = {
    	lag = TRUE
    	spatial.error = FALSE
    	},
    error = {
    	lag = FALSE
    	spatial.error = TRUE
    	},
    sarar = {
    	lag = TRUE
    	spatial.error = TRUE
    	})

    errors <- match.arg(errors)

    x0 <- update(x, .~.-1)

    method <- switch(match.arg(method),
                     ML = {
                         ## adapt argument
                         spatial.error <- if(spatial.error) {
                             spatial.error <- if(errors=="BSK") "b" else "kkp"
                         } else {
                             spatial.error <- "none"
                         }
                         femod <- spml(x, data = data, index = index, listw = listw, lag = lag,
                                       spatial.error = spatial.error, model = "within")
                         remod <- spml(x, data = data, index = index, listw = listw, lag = lag,
                                       spatial.error = spatial.error, model = "random")
                         },
                     GM = {
                         femod <- spgm(x, data = data, index = index, listw = listw, lag = lag,
                                       spatial.error = spatial.error, model = "within", moments = "fullweights")
                         remod <- spgm(x, data = data, index = index, listw = listw, lag = lag,
                                       spatial.error = spatial.error, model = "random", moments = "fullweights")
                         },
                     stop("\n Unknown method"))
    
    return(sphtest(femod, remod, ...))
    }

sphtest.splm <- function (x, x2, ...){
    ## check whether the models have been estimated by GM (different slots...)
    is.gm <- !is.null(x$ef.sph)

    if(is.gm) {
        ## check that the models have the same specification but different effects
        if (!all.equal(x$legacy, x2$legacy)) stop("The models are different")
        if(x$ef.sph == x2$ef.sph) stop("Effects should be different")

        ran <- match("random", c(x$ef.sph, x2$ef.sph))
        if(ran == 1){
            xwith <- x2
            xbetw <- x
        }
        if(ran == 2){
	    xwith <- x
	    xbetw <- x2
        }
    
      ## test on coefficients (excluding SAR)
      ## model order is irrelevant
  
      tc <- match(names(coef(xwith)), names(coef(xbetw)) )

      coef.wi <- coef(xwith)
      coef.re <- coef(xbetw)[tc]
      vcov.wi <- xwith$vcov
      vcov.re <- xbetw$vcov[tc,tc]

    } else {

        ## then they are ML

        ## determine which is FE
        if(is.null(dimnames(x$vcov))) {
            xwith <- x
            xbetw <- x2
        } else {
            xwith <- x2
            xbetw <- x
        }
            

        tc <- intersect(names(coef(xwith)), names(coef(xbetw)))
        ## fix because vcov for FE is not named. Aaaargh!
        wtc <- match(tc, names(coef(xwith)))
        
        coef.wi <- coef(xwith)[wtc]
        coef.re <- coef(xbetw)[tc]
        vcov.wi <- xwith$vcov[wtc,wtc]
        vcov.re <- xbetw$vcov[tc,tc]
        
    }


    
    dbeta <- coef.wi - coef.re
    df <- length(dbeta)
    dvcov <- vcov.re - vcov.wi
    stat <- abs(t(dbeta) %*% solve(dvcov) %*% dbeta)
    pval <- pchisq(stat, df = df, lower.tail = FALSE)
    names(stat) <- "chisq"
    parameter <- df
    names(parameter) <- "df"
    data.name <- paste(deparse(x$call$formula))
    alternative <- "one model is inconsistent"
    res <- list(statistic = stat, p.value = pval, parameter = parameter,
                method = "Hausman test for spatial models",
                data.name = data.name, alternative = alternative)
    class(res) <- "htest"
    return(res)
}
