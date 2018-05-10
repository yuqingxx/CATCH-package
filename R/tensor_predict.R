predict.tsda<-function(object,newx,...){
	beta<-object$beta
	mu<-object$mu
	prior<-object$prior
	nclass<-length(prior)		
	dimen<-dim(newx[[1]])
	nvars<-prod(dimen)
	mubar=matrix(list(),nclass-1,1)
	for (i in 1:(nclass-1)){
		mubar[[i]] = (mu[[i+1]]+mu[[1]])/2
	}
	n<-length(newx)
	nn<-length(object$x)
	x.train<-object$x
	vecx.train = matrix(0,ncol=nn,nrow=nvars)
	vecnewx = matrix(0,ncol=n,nrow=nvars)
	for (i in 1:nn){
		vecx.train[,i]<-matrix(x.train[[i]],ncol=1)	
	}
	vecx.train = t(vecx.train)
	for (i in 1:(length(newx))){
		vecnewx[,i]<-matrix(newx[[i]],ncol=1)	
	}
	vecnewx = t(vecnewx)
	y.train<-object$y
	nlambda<-length(beta)
	pred<-matrix(0,n,nlambda)
	pred[1]<-which.max(prior)
	for (i in 1:nlambda){
		nz<-sum(beta[[i]][,1]!=0)	
		if (nz == 0){
			pred[,i]<-which.max(prior)
		}else{	
			xfit<-vecx.train %*% beta[[i]][,1:(min(nclass-1,nz)),drop=FALSE]
			xfit.sd<-matrix(0,nclass,ncol(xfit))
			for (j in 1:nclass){
				xfit.sd[j,]<-apply(xfit[y.train==j,,drop=FALSE],2,sd)			
			}
			xfit.sd<-apply(xfit.sd,2,min)
			if (min(xfit.sd)<1e-4){pred[,i]<-which.max(prior)}else{
				l<-lda(xfit, y.train)
				pred[,i]<-predict(l,vecnewx%*%beta[[i]][,1:(min(nclass-1,nz))])$class
			}
		}
	}
	pred
}