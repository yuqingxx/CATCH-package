##data 
##x a list of arrays
##y vector



# extract fortran outputs and format it into sparse matries
formatoutputt <- function(fit, maxit, pmax, dfmax, nvars, vnames, nk) {
    nalam <- fit$nalam
    ntheta <- fit$ntheta[seq(nalam)]
    nthetamax <- max(ntheta)
    lam <- fit$alam[seq(nalam)]
  	obj <- fit$obj[seq(nalam)]
    stepnames <- paste("s", seq(nalam) - 1, sep = "")
    resnames <- paste("delta", seq(nk), sep = "")
    
    errmsg <- errt(fit$jerr, maxit, pmax,dfmax)  ### error messages from fortran
    switch(paste(errmsg$n), `1` = stop(errmsg$msg, call. = FALSE), `-1` = cat(errmsg$msg))
    
    dd <- c(nvars, nk)
    df <- rep(0, nalam)
    if (nthetamax > 0) {
        ja <- fit$itheta[seq(nthetamax)]
        oja <- order(ja)
        ja <- rep(ja[oja], nk)
        itheta <- cumsum(c(1, rep(nthetamax, nk)))
        pos <- rep(1:nalam, each = nk * pmax)
        theta <- split(fit$theta[seq(nk * pmax * nalam)], pos)
        for (l in 1:nalam) {
            theta[[l]] <- matrix(theta[[l]], pmax, nk, byrow = TRUE)[seq(nthetamax), 
                , drop = FALSE]
            df[l] <- sum(rowSums(abs(theta[[l]])) != 0)
            theta[[l]] <- new("dgCMatrix", Dim = dd, Dimnames = list(vnames, 
                resnames), x = as.vector(theta[[l]][oja, ]), p = as.integer(itheta - 
                1), i = as.integer(ja - 1))
        }
    } else {
        theta <- list()
        for (l in 1:nalam) {
            theta[[l]] <- zeromat(nvars, nk, vnames, resnames)
        }
        df <- rep(0, nalam)
    }
    list(beta = theta, df = df, dim = dd, lambda = lam, obj = obj)
}


# generate sigma, delta and mu from x, y
prept <- function(x, y) {
    	# data setup
    	nclass <- as.integer(length(unique(y)))
    	prior <- rep(0, nclass)
    	for (k in 1:nclass) {
     	   prior[k] <- mean(y == k)
    	}
	    ldim=length(dim(x[[1]])) #m dimensional
    	p=matrix(0,ncol=ldim,nrow=1)
    	nvars=1
    	for (i in 1:ldim){
    		p[i]=as.integer(dim(x[[1]])[i])
    		nvars=nvars*p[i]
    	}
     	nobs <- length(x)
    	nres <- length(y)
      if (nres != nobs) 
         stop("x and y have different number of observations")
	
    	dim=dim(x[[1]])
    	num=matrix(0,nrow=k,ncol=1)
    	# estimate mean
    	es_M=array(list(),k)
    	for (i in 1:k){
    		es_M[[i]]=array(0,dim=dim)
    	}
      nk=nclass-1
    	for (i in 1:nobs){
	    	es_M[[y[i]]]=es_M[[y[i]]]+x[[i]]	
    		num[[y[i]]]=num[[y[i]]]+1
    	}
    	for (i in 1:k){
    		es_M[[i]]=es_M[[i]]/num[[i]]
    	}

	#estimate sigma
    	n=nobs
    	matx=matrix(list(),n,1)
    	es_sigma=matrix(list(),ldim,1)
    	sum_sigma=matrix(list(),ldim,1)
      matmu=matrix(list(),k,1)
    	tr=1
    	for (j in 1:ldim){
		    for (i in 1:k){
			    matmu[[i]]=mat(es_M[[i]],j)
		    }
		    sum_sigma[[j]]=matrix(0,nrow=dim[j],ncol=dim[j])
		    for (i in 1:n){
			    matx[[i]]=mat(x[[i]],j)
			    sum_sigma[[j]]=sum_sigma[[j]]+(matx[[i]]-matmu[[y[i]]])%*%t(matx[[i]]-matmu[[y[i]]])
		    }
		    sum_sigma[[j]]=sum_sigma[[j]]/(n-k)
		    if (j<ldim){
			    es_sigma[[j]]=sum_sigma[[j]]/sum_sigma[[j]][1,1]
			    tr=tr*sum(diag(es_sigma[[j]]))
		    }else
	    	{
		    	es_sigma[[j]]=sum_sigma[[j]]/tr
		    }
  	  }
      delta=matrix(rep(0),nrow=nclass-1,ncol=prod(dim))
    	for (i in 1:(nclass-1)){
		    delta[i,]=as.matrix(as.vector(es_M[[i+1]])-as.vector(es_M[[1]]))
    	}
      
      maxd=max(dim)    
      sigma=matrix(0,nrow=ldim,ncol=(maxd^2))
      for (i in 1:ldim){
        t<-matrix(0,nrow=maxd,ncol=maxd)
        t[1:dim[i],1:dim[i]]=es_sigma[[i]]
        sigma[i,]=as.vector(t)
      }    
	    outlist <- list(sigma=sigma, delta = delta, mu =es_M, prior = prior)
    	outlist
}

errt <- function(n, maxit, pmax,dfmax) {
    if (n == 0) 
        msg <- ""
    if (n > 0) {
        # fatal error
        if (n < 7777) 
            msg <- "Memory allocation error; contact package maintainer"
        if (n == 10000) 
            msg <- "All penalty factors are <= 0"
        n <- 1
        msg <- paste("in the fortran code -", msg)
    }
    if (n < 0) {
        # non fatal error
        if (n > -10000) 
            msg <- paste("Convergence for ", -n, "th lambda value not reached after maxit=", 
                maxit, " iterations; solutions for larger lambdas returned.\n", 
                sep = "")
        if (n < -10000) 
            msg <- paste("Number of nonzero coefficients along the path exceeds pmax=", 
                pmax, " at ", -n - 10000, "th lambda value; solutions for larger lambdas returned.\n", 
                sep = "")
        if (n < -20000) 
            msg <- paste("Number of nonzero coefficients along the path exceeds dfmax=", 
                dfmax, " at ", -n - 20000, "th lambda value; solutions for larger lambdas returned.\n", 
                sep = "")
        n <- -1
    }
    list(n = n, msg = msg)
}

zeromat <- function(nvars, nalam, vnames, stepnames) {
    ca <- rep(0, nalam)
    ia <- seq(nalam + 1)
    ja <- rep(1, nalam)
    dd <- c(nvars, nalam)
    new("dgCMatrix", Dim = dd, Dimnames = list(vnames, stepnames), x = as.vector(ca), 
        p = as.integer(ia - 1), i = as.integer(ja - 1))
}

lamfix <- function(lam) {
    llam <- log(lam)
    lam[1] <- exp(2 * llam[2] - llam[3])
    lam
}
