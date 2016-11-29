`print.summary.splm` <- function(x, digits=max(3, getOption("digits") - 2),
                                 width=getOption("width"), ...) {

    ## manage model description
    if(grepl("random", x$type)) {
            ## "random" models by spreml() have a more complicated description
            m.des <- x$type.des
        } else {
            m.des <- paste("Spatial panel", x$type,"model\n")
        }
    cat(paste(m.des, "\n"))

    ## print call
    cat("\nCall:\n")
    print(x$call)

    ## print residual des
    cat("\nResiduals:\n")
    save.digits <- unlist(options(digits=digits))
    on.exit(options(digits=save.digits))
    print(sumres(x))

    ## if model is of 'random' type ex spreml():
    if(grepl("random", x$type)) {

        ## print error components' table for 'random' models
        if(!is.null(x$ErrCompTable)) {
            cat("\nError variance parameters:\n")
            printCoefmat(x$ErrCompTable, digits=digits, signif.legend=FALSE)
        }

        ## print spatial lag coefficient for 'random' models
        if(!is.null(x$ARCoefTable)) {
            cat("\nSpatial autoregressive coefficient:\n")
            printCoefmat(x$ARCoefTable, digits=digits, signif.legend=FALSE)
        }

        ## print betas
        cat("\nCoefficients:\n")
        printCoefmat(x$CoefTable,  digits=digits)
        cat("\n")
        
    } else {

        ## then it is of 'fixed' type ex spfeml()

        ## print spatial lag coefficient (is this condition ever true??)
        if(is.numeric(x$lambda)) {
            cat("\nEstimated spatial coefficient, variance components and theta:\n")
            print(x$lambda)
        }
        
        ## print error components' table for 'random' models
        if("rho" %in% dimnames(x$CoefTable)[[1]]) {
            cat("\nSpatial error parameter:\n")
            printCoefmat(x$CoefTable["rho", , drop=FALSE], digits=digits, signif.legend=FALSE)
        }

        ## print spatial lag coefficient for 'random' models
        if("lambda" %in% dimnames(x$CoefTable)[[1]]) {
            cat("\nSpatial autoregressive coefficient:\n")
            printCoefmat(x$CoefTable["lambda", , drop=FALSE], digits=digits, signif.legend=FALSE)
        }

        ## print betas (w/o spatial coefs)
        cat("\nCoefficients:\n")
        spat.nam <- dimnames(x$CoefTable)[[1]] %in% c("rho","lambda")
        printCoefmat(x$CoefTable[!spat.nam, , drop=FALSE], digits=digits)
        cat("\n")

    }
        
    invisible(x)
}

