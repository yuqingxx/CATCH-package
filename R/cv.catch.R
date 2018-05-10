cv.catch<-function(x, z=NULL,y,nfolds=5, lambda=NULL, lambda.opt="min",...){
	y<-drop(y)
	n<-length(x)
	dimen=dim(x[[1]])
	nvars=prod(dimen)
	K<-length(unique(y))
	prior<-rep(0,K)
	for (i in 1:K){
		prior[i]<-mean(y==i)
	}
	### Fit the model once to get dimensions etc of output
	if (!is.null(z)){
	  z = as.matrix(z)
	  q = dim(z)[2]
		obj <- adjten(x,z,y)
	  tmp <- tsda(obj$xres, y, lambda=lambda,...)
	  }else
	  {
	    tmp <- tsda(x,y,lambda=lambda,...)
	  }
	lambda<- tmp$lambda
	### Now Fit the nfold models and store them
	foldid <- sample(rep(seq(nfolds), length = n))

 	if (nfolds < 3) 
          stop("nfolds must be bigger than 3; nfolds=5 recommended")
      if (nfolds > n) 
          stop("The number of folds should be smaller than the sample size.")
	residmat<-matrix(NA,nfolds,length(lambda))
	for (i in seq(nfolds)){
		which<- foldid==i
		
		if (is.null(z) == TRUE){
		  fit <- tsda(x[!which, drop=FALSE], y[!which],lambda=lambda,...)
		  preds <- predict.tsda(fit, x[which, drop=FALSE])
		
		}else{
	  	xtr=x[!which]
		  xte=x[which]
		  ztr=z[!which,]
		  zte=z[which,]
	  	obj<-adjten(xtr,ztr,y[!which],xte,zte)
		  fit<-tsda(obj$xres,y[!which],lambda=lambda,...)
		  preds<-predict.catch(fit,obj$testxres,ztr,zte,obj$alpha,...)			
		}
		nlami <- length(fit$lambda)
		residmat[i,seq(nlami)] <- colMeans(y[which]!=preds)
		
	}
	residmat[is.na(residmat)] <- 1
	cvm <- colMeans(residmat, na.rm = TRUE)
	cvsd <- sqrt(colMeans(scale(residmat, cvm, FALSE)^2, na.rm = TRUE)/(nfolds - 1))
      if (lambda.opt == "min") {
        lambda.min <- min(lambda[which(cvm == min(cvm, na.rm = TRUE))])
      } else {
        lambda.min <- max(lambda[which(cvm == min(cvm, na.rm = TRUE))])
      }
	id.min<-which.min(cvm)
  if (cvsd[id.min]==0){
    lambda.1se=lambda.min
    }else{
  	lambda.1se <- max(lambda[cvm < min(cvm, na.rm = TRUE) + cvsd[id.min]], na.rm = TRUE)
    }
	obj <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, 
        lambda.min = lambda.min, lambda.1se = lambda.1se, catch.fit = tmp)
    	class(obj) <- "cv.catch"
    	obj
}
