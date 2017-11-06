slmtest<-function (x, ...)
 {
    UseMethod("slmtest")
 }

slmtest.plm <- function(x, listw,
                        test=c("lme", "lml",
                               "rlme", "rlml"),
                        ...)
{
    ## plm method for slmtest
    return(slmtestres(mod = x, listw = listw,
                      test = match.arg(test)))   
}

slmtest.formula <- function(formula, data, listw,
                            model="pooling",
                            test=c("lme", "lml",
                                   "rlme", "rlml"),
                            index=NULL, ...) {
    ## estimate pooled model
    ## (notice you get "standard" panel ordering!
    ## i.e., id is slow and time is fast)
    mod <- plm(formula=formula, data=data,
               model=model, index=index, ...)
    return(slmtestres(mod = mod, listw = listw,
                      test = match.arg(test)))
    }
                      

slmtestres <- function(mod, listw, test) {

    ## Computing engine for slmtest; expects a 'plm' object
    N <- pdim(mod)$nT$n
    t. <- pdim(mod)$nT$T
    NT <- pdim(mod)$nT$N

    X <- model.matrix(mod)
    y <- as.numeric(pmodel.response(mod))
    hatY <- X %*% coef(mod)

    ## check w or listw and in case transform to matrix
    w <- listw
    if (!is.matrix(w)) {
        if ("listw" %in% class(w)) {
            w <- listw2mat(w)
        }
        else {
            stop("listw has to be either a 'matrix' or a 'listw' object")
        }
    }
    
    ## instead of reordering y, X, e... the 'spatial panel'
    ## way, we just swap I and W in
    ## making the bigW: the rest is all full panel vectors.
    ## Optimization is here done by using sparse matrix
    ## methods.
    W <- as.spam(w)
    bigW <- kronecker(W, diag.spam(1, t.)) 
    tr <- function(x) sum(diag(x))
    Tw <- tr(W%*%W + crossprod(W))
    M <- diag.spam(1, NT) - X %*% solve(crossprod(X)) %*% t(X)
    Whaty <- bigW %*% hatY

    ## extract residuals as a vector (no pseries features)
    e <- as.numeric(mod$residuals)

    sigma2 <- crossprod(e)/NT

    J <- (crossprod(Whaty, M) %*% Whaty) / sigma2 + t.*Tw

    switch(test, lml = {
        statistic <- ((crossprod(e, bigW) %*% y)/sigma2)^2 / J
        descr <- "lag"
        rob <- ""
    }, lme = {
        statistic <- ((crossprod(e, bigW) %*% e)/sigma2)^2 /
            (t. * Tw)
        descr <- "error"
        rob <- ""
    }, rlml = {
        nume <- ((crossprod(e, bigW) %*% y)/sigma2)-
            ((crossprod(e, bigW) %*% e)/sigma2)
        deno <- J - t.*Tw
        statistic <- nume^2/deno
        descr <- "lag"
        rob <- " sub spatial error"
    }, rlme = {
        nume <- ((crossprod(e, bigW) %*% e)/sigma2)-
            t.*Tw/J * ((crossprod(e, bigW) %*% y)/sigma2)
        deno <- t.*Tw * (1 - t.*Tw/J)
        statistic <- nume^2/deno
        descr <- "error"
        rob <- " sub spatial lag"
    })

    names(statistic) <- "LM"
    df <- 1
    names(df) <- "df"
    model <- mod$args$model
    pre <- if(rob == "") "" else "Locally robust "
    transf.type <- if(model == "pooling") "" else {
         paste("(", model, " transformation)", sep="")}     
    alternative = paste("spatial", descr, "dependence")
    form <- paste(deparse(substitute(formula)),
                  transf.type)
    
    p.value <- pchisq(statistic, df=df, lower.tail=F)

    RVAL <- list(statistic = statistic, parameter = df,
                 method = paste(pre, "LM test for spatial ",
                                descr, " dependence",
                                rob, sep=""),
        alternative = alternative, p.value = p.value, 
        data.name = form)
    class(RVAL) <- "htest"
    return(RVAL)
}

## check
## lm.LMtests(lm(fm, Produc), listw=mat2listw(kronecker(diag(1, 17), usaww)))
## lm.LMtests(lm(fm, Produc), listw=mat2listw(kronecker(diag(1, 17), usaww)), test="LMlag")
## lm.LMtests(lm(fm, Produc), listw=mat2listw(kronecker(diag(1, 17), usaww)), test="RLMerr")
## lm.LMtests(lm(fm, Produc), listw=mat2listw(kronecker(diag(1, 17), usaww)), test="RLMlag")
