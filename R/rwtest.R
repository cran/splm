
## this version 2 (10 for mcwtest): pass on pseries of residuals,
## then reshape and use fast cor() to calculate the rho matrix
## (much speedier); from changes to pcdtest in plm_1.6-2
## moreover, parallelized simulation step.
## Giovanni Millo, 17/11/2016


rwtest<-function (x, ...)
 {
    UseMethod("rwtest")
 }


rwtest.formula <- function (x, data, w, index = NULL,
                            model = NULL,
                            replications=99, seed=NULL,
                            order=1, mc=1,
                            test = c("rho", "cd", "sclm"),
                            alternative=c("twosided",
                                          "onesided",
                                          "symmetric"), ...) {

    ## much like pcdtest.formula
    mymod <- plm(x, data, index = index,
                 model = "pooling", ...)
    if (is.null(model) & min(pdim(mymod)$Tint$Ti) <
        length(mymod$coefficients) + 1) {
        warning("Insufficient number of observations in time to estimate heterogeneous model: using within residuals",
            call. = FALSE)
        model <- "within"
    }
    
    ind0 <- attr(model.frame(mymod), "index")
    tind <- as.numeric(ind0[[2]])
    ind <- as.numeric(ind0[[1]])
    
    if (is.null(model)) {
        ## estimate individual regressions one by one
        X <- model.matrix(mymod)
        y <- model.response(model.frame(mymod))
        unind <- unique(ind)
        n <- length(unind)
        ti.res <- vector("list", n)
        ind.res <- vector("list", n)
        tind.res <- vector("list", n)
        for (i in 1:n) {
            tX <- X[ind == unind[i], , drop = FALSE]
            ty <- y[ind == unind[i]]
            res.i <- lm.fit(tX, ty)$resid
            ti.res[[i]] <- res.i
            names(ti.res[[i]]) <- tind[ind == unind[i]]
            ind.res[[i]] <- rep(i, length(res.i))
            tind.res[[i]] <- tind[ind == unind[i]]
        }
        ## make pseries of (all) residuals
        resdata <- data.frame(ee = unlist(ti.res),
                              ind = unlist(ind.res),
                              tind = unlist(tind.res))
        pee <- pdata.frame(resdata, index = c("ind", "tind"))
        tres <- pee$ee
    }
    else {
        mymod <- plm(x, data, index = index, model = model, ...)
        tres <- resid(mymod)
        unind <- unique(ind)
        n <- length(unind)
        t <- min(pdim(mymod)$Tint$Ti)
        nT <- length(ind)
        k <- length(mymod$coefficients)
        }

    return(mcwres(tres=tres, n=n, w=w,
                  form=paste(deparse(substitute(formula))),
                  test=match.arg(test),
                  alternative=match.arg(alternative),
                  replications=replications,
                  seed=seed, order=order, mc=mc))

}

rwtest.pseries <- function(x, w, replications=99,
                           seed=NULL, order=1, mc=1,
                           test = c("rho", "cd", "sclm"),
                           alternative=c("twosided",
                                         "onesided",
                                         "symmetric"), ...) {
  ## get indices
  index <- attr(x, "index")
  tind <- as.numeric(index[[2]])
  ind <- as.numeric(index[[1]])
  n <- length(unique(ind))

  return(mcwres(tres=x, n=n, w=w,
                form=paste(deparse(substitute(formula))),
                test=match.arg(test),
                alternative=match.arg(alternative),
                replications=replications, seed=seed,
                order=order, mc=mc))

}

    
rwtest.panelmodel <- function(x, w, replications=99,
                              seed=NULL, order=1, mc=1,
                              test = c("rho", "cd", "sclm"),
                              alternative=c("twosided",
                                            "onesided",
                                           "symmetric"), ...) {

  ## this is taken from the last piece of pcdtest.formula,
  ## after estimating relevant model

  ## fetch residuals
  tres <- resid(x)

  ## get indices
  index <- attr(model.frame(x), "index")
  tind <- as.numeric(index[[2]])
  ind <- as.numeric(index[[1]])

  ## det. number of groups
  unind<-unique(ind)
  n<-length(unind)

  return(mcwres(tres=tres, n=n, w=w,
                form=paste(deparse(substitute(formula))),
                test=match.arg(test),
                alternative=match.arg(alternative),
                replications=replications, seed=seed,
                order=order, mc=mc))

}

mcwres<-function(tres, n, w, form, test, alternative,
                 replications, seed, order, mc) {
    ## as in pcdtest;
    ## now tres is the pseries of model residuals

    ## calc matrix of all possible pairwise corr.
    ## coeffs. (200x speedup from using cor())
    wideres <- t(preshape(tres, na.rm=FALSE))
    rho <- cor(wideres, use="pairwise.complete.obs")
    
    ## find length of intersecting pairs
    ## fast method, times down 200x
    data.res <- data.frame(time=attr(tres, "index")[[2]],
                           indiv=attr(tres, "index")[[1]])
    ## tabulate which obs in time for each ind are !na
    presence.tab <- table(data.res)
    ## calculate t.ij
    t.ij <- crossprod(presence.tab)
    
    ## input check
    if (!is.null(w)) {
        dims.w <- dim(w)
        if(dims.w[1] != n || dims.w[2] != n)
            stop(paste0("matrix 'w' describing proximity of individuals has wrong dimensions: ",
                        "should be ", n, " x ", n,
                        " (no. of individuals) but is ",
                        dims.w[1], " x ", dims.w[2]))
    }
   
  ## 'til here as in panelCDtest. Now for the selection of
  ## rho's based on (p-th order) contiguity
    
  ## make proximity matrix up to the p-th lag
  ## TODO: outsource this to wlag()

  if(order>1) {
      ## if it is a matrix make it an nb object
      if(is.matrix(w)) w<-mat2listw(w)
      ## make lagged proximity matrix
      w<-wlag(w$neighbours,maxlag=order)
  } else {
      ## if order=1, just make w<-listw2mat(w) if it is
      ## a listw obj.
      if("listw"%in%class(w)) {
          w<-listw2mat(w)
      }
  }

    ## until here all the necessary rho_ij's have been
    ## calculated;
    ## the relevant W matrix has been determined;
    ## now for the randomization of W and collection of the
    ## distr. of CD(p)s

    ##### begin: random W snippet (from MCWtest.R) #####
    ##### (until end) ##################################

    mcwdist<-rep(NA,replications+1)
    
    ## calc. ''true'' CD(p) statistic:

    ## make (binary) selector matrix based on the contiguity
    ## matrix w and extracting elements corresponding to ones
    ## in the lower triangle excluding the diagonal

    ## transform in logicals (0=FALSE, else=TRUE: no need to
    ## worry about row-std. matrices)
    selector.mat<-matrix(as.logical(w),ncol=n)
    ## set upper tri and diagonal to false
    selector.mat[upper.tri(selector.mat,diag=TRUE)]<-FALSE

    ## number of elements in selector.mat
    elem.num<-sum(selector.mat)

    ## set components not dependent on w once for all,
    ## depending on test type

    switch(test,
           cd = {
               ## Pesaran statistic for (p-th order) local
               ## cross-sectional dependence
               rad.n<-sqrt(1/elem.num)
               t.rhos<-sqrt(t.ij)*rho
               testname<-"CD test"
           },
           sclm = {
               ## Scaled LM statistic for (p-th order) local
               ## cross-sectional dependence
               rad.n<-sqrt(1/(2*elem.num))
               t.rhos<-t.ij*rho^2-1
               testname<-"NLM test"
           },
           rho = {
               rad.n <- 1/elem.num
               t.rhos <- rho
               testname <- "Average correlation coefficient"
           }
           )

    cdpstat <- rad.n*sum(t.rhos[selector.mat])

    ## calc. ''randomized'' CD(p) stats
    ## reuse t.rhos, rad.n, independent from w's

    if(!is.null(seed)) set.seed(seed)
    rwfun <- function(x) rad.n*sum(t.rhos[brandW(w)])
    lfun <- if(mc==1) {
                lapply}
            else {
                    function(x, FUN, ...) {
                        mclapply(x, FUN, mc.cores=mc, ...)
                    }
                }
    mcwlist <- lfun(seq_along(1:replications), FUN=rwfun)
    mcwdist <- c(cdpstat, unlist(mcwlist))
    
    ## how many random draws are more extreme?
    switch(alternative,
           twosided={quasi.pval <- 2 * min(sum(mcwdist < cdpstat)/(replications+1), sum(mcwdist >= cdpstat)/(replications+1))
           },
           onesided={
               quasi.pval <- sum(mcwdist >= cdpstat)/(replications+1)
           },
           symmetric={
               quasi.pval <- sum(abs(mcwdist) >= abs(cdpstat))/(replications+1)
           })
    myres <- quasi.pval
    
    ## setup htest object
    RVAL <- list(statistic = NULL, parameter = NULL,
                 method = paste("Randomized W test for spatial correlation of order", order),  #, " \n\n Based on:", testname), 
                 alternative = alternative,
                 p.value = quasi.pval,
                 data.name = form)
    class(RVAL) <- "htest"
    return(RVAL)
}

## utilities

wlag<-function(x,maxlag) {
    ## define (incremental!) lag function for w
    ## returns the nb list of all neighbours up to order=maxlag
    ## standalone version under wlag.R

    n<-length(x)
    mynb<-nblag(x,maxlag=maxlag)
    
    mytot<-vector("list",n)
    
    for(i in 1:n) {
        mytot[[i]]<-mynb[[1]][[i]]
        for(j in 2:maxlag) mytot[[i]] <- c(mytot[[i]],
                                           mynb[[j]][[i]])
        ## reorder
        mytot[[i]]<-mytot[[i]][order(mytot[[i]])]
    }
    lagmat<-matrix(0,ncol=n,nrow=n)
    for(i in 1:n) lagmat[i,mytot[[i]]]<-1
    return(lagmat)
}

brandW<-function(w) {
  ## randomize a contiguity matrix, boolean output matrix
  ## argument: a contiguity matrix (doesn't matter whether
  ## standardized or not)

  n<-dim(w)[[1]]
  brW<-matrix(FALSE,ncol=n,nrow=n)

  ## calc. filling rate
  ## (sum(ones)/(total cells-diagonal) must be
  ## integer and even, by construction)

  ## make everything binary:
  booleW<-as.logical(w)
  ## calc. (diagonal is always excluded)
  fillW<-sum(booleW)/2

  ## fill in randomly the lower part of matrix with
  ## fillW 'TRUE's
  smplW<-sample(n*(n-1)/2,fillW)   # w/o replacement
  brW[lower.tri(brW)][smplW]<-TRUE

  ## don't need to mirror upper triangle any more
  ## (only lower tri is used in CD(p) test
  return(brW)
 }

preshape <- function(x, na.rm=TRUE, ...) {
    ## to be substituted by importing plm:::preshape
    
    ## reshapes pseries, e.g. of residuals from a
    ## panelmodel, in wide form
    inames <- names(attr(x, "index"))
    mres <- reshape(cbind(as.vector(x), attr(x, "index")),
                    direction="wide",
                    timevar=inames[2], idvar=inames[1])
    ## drop ind in first column
    mres <- mres[,-1]
    ## reorder columns (may be scrambled depending on first
    ## available obs in unbalanced panels)
    mres <- mres[, order(dimnames(mres)[[2]])]
    ## if requested, drop columns (time periods) with NAs
    if(na.rm) {
        rmc <- which(is.na(apply(mres, 2, sum)))
        if(sum(rmc)>0) mres <- mres[,-rmc]
    }
    return(mres)
}

