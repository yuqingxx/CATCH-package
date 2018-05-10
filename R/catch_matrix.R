catch_matrix <- function(x, z=NULL, y, testx=NULL,testz=NULL,...){
  #input x is matrix * size
  pred = NULL
  nobs=length(y)
  nvars=prod(dim(x))/nobs
  nclass=length(unique(as.factor(y)))
  
  ntrain=dim(x)[3]
  newx=array(list(),ntrain)
  for (i in 1:ntrain){
    newx[[i]]=x[,,i]
  }
  newtestx=NULL
  if (!is.null(testx)){
    ntest=dim(testx)[3]
    newtestx=array(list(),ntest)
    for (i in 1:ntest){
      newtestx[[i]]=testx[,,i]
    }
  }
  catch(newx,z,y,newtestx,testz,...)
  
}

