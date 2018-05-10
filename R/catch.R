catch <- function(x, z=NULL, y, testx=NULL,testz=NULL, nlambda = 100, 
                  lambda.factor = ifelse((nobs - nclass) <= nvars, 0.2, 1E-03), 
                  lambda = NULL, dfmax = nobs, pmax = min(dfmax * 2 + 20, nvars), 
                  pf = rep(1, nvars), eps = 1e-04, maxit = 1e+05, sml = 1e-06, verbose = FALSE, perturb = NULL){
 # dyn.load('tensor.dll')
  pred = NULL
  nobs=length(y)
  nvars=prod(dim(x[[1]]))
  nclass=length(unique(as.factor(y)))
  if (is.null(testx)){
    if (is.null(z)){
      objt <- tsda(x,y,nlambda,lambda.factor,lambda,dfmax,pmax,pf,eps,maxit,sml,verbose,perturb)
    }else{
      obj <- adjten(x,z,y)
      objt <- tsda(obj$xres,y,nlambda,lambda.factor,lambda,dfmax,pmax,pf,eps,maxit,sml,verbose,perturb)
    }
  }
  else{  
   if (is.null(z) && is.null(testz)){
     objt <- tsda(x,y,nlambda,lambda.factor,lambda,dfmax,pmax,pf,eps,maxit,sml,verbose,perturb)
     pred <- predict.tsda(objt,testx)
   }else{
      if (is.null(z)){
         print('Covariates for training data are missing.')
         return()
      }
      if (is.null(testz)){
         print('Covariates for testing data are missing.')
        return()
      }
      obj <- adjten(x,z,y,testx,testz)
      objt <- tsda(obj$xres,y,nlambda,lambda.factor,lambda,dfmax,pmax,pf,eps,maxit,sml,verbose,perturb)
      pred <- predict.catch(objt,obj$testxres,z,testz,obj$gamma)
    }
  }
  outlist=c(objt,list(pred=pred))
  class(outlist) <- c("catch")
  return(outlist)
  
}

                                                                