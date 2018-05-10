#predictadj.msdat<-function(object,newx,z,ztest,gamma,...){
predict.catch<-function(object,newx,z=NULL,ztest=NULL,gamma=NULL,...){

  if (is.null(gamma)){
    pred=predict.tsda(object,newx)
  }else{
    thetatm <- object$beta
	  theta = array(list(),length(thetatm))
	  mu<- object$mu
	  prior <- object$prior
	  nclass <- length(prior)	
	  dimen <- dim(newx[[1]])
	  nvars <- prod(dimen)
	  nlambda <- length(theta)
	  gamma= as.matrix(gamma)
	  q = dim(gamma)[1]
	  z=as.matrix(z)
	  ztest=as.matrix(ztest)
	  for (i in 1:nlambda){
		  theta[[i]] = matrix(0,nrow=dim(thetatm[[i]])[1]+q,ncol=nclass-1)
		  for (j in 1:nclass-1){
			  theta[[i]][1:nvars,j] = matrix(thetatm[[i]][,j],ncol=1)
			  for (qq in 1:q){
				  theta[[i]][nvars+qq,j] = gamma[qq,j]
			  }
		  }
	  }
	  mubar = matrix(list(),nclass-1,1)
	  for (i in 1:(nclass-1)){
		  mubar[[i]] = (mu[[i+1]]+mu[[1]])/2
	  }
	  n <- length(newx)
	  nn <- length(object$x)
	  x.train <- object$x
	  vecx.train = matrix(0,ncol=nn,nrow=nvars+q)
	  vecnewx = matrix(0,ncol=n,nrow=nvars+q)
	  for (i in 1:nn){
		  vecx.train[1:nvars,i] <- matrix(x.train[[i]],ncol=1)
		  for (qq in 1:q){
			  vecx.train[nvars+q,i] = z[i,q]
		  }
	  }
	  vecx.train = t(vecx.train)
	  for (i in 1:(length(newx))){
		  vecnewx[1:nvars,i] <- matrix(newx[[i]],ncol=1)
		  for (qq in 1:q){
			  vecnewx[nvars+q,i] = ztest[i,q]
		  }
	  }
	  vecnewx = t(vecnewx)
	  y.train <- object$y
	  pred <- matrix(0,n,nlambda)
	  pred[1] <- which.max(prior)
	  for (i in 1:nlambda){
		  nz <- sum(theta[[i]][,1]!=0)	
		  if (nz == 0){
		  	pred[,i] <- which.max(prior)
		  }else{	
			  xfit <- vecx.train %*% theta[[i]][,1:(min(nclass-1,nz)),drop=FALSE]
			  xfit.sd <- matrix(0,nclass,ncol(xfit))
			  for (j in 1:nclass){
			  	xfit.sd[j,] <- apply(xfit[y.train==j,,drop=FALSE],2,sd)	
			  }
			  xfit.sd <- apply(xfit.sd,2,min)
			  if (min(xfit.sd)<1e-4){pred[,i]<-which.max(prior)}else{
			  	l <- lda(xfit, y.train)
			  	pred[,i] <- predict(l,vecnewx%*%theta[[i]][,1:(min(nclass-1,nz))])$class
			  }
		  }
	  }
  }
	pred
}