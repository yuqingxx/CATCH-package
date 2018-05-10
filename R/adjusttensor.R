adjten <- function(x,z,y,testx=NULL,testz=NULL,is.centered=FALSE){
  

  n = length(x)
  dimen = dim(x[[1]])
  nvars = prod(dimen)
  vec_x = matrix(0,nvars,n)
  for (i in 1:n){
    vec_x[,i] = matrix(x[[i]],nvars,1)
  }
  vec_x = t(vec_x)
  
  if (!is.null(testx)){
    tesize=length(testx)
    vec_testx = matrix(0,nvars,tesize)
    for (i in 1:tesize){
      vec_testx[,i]=matrix(testx[[i]],nvars,1)
    }
    vec_testx = t(vec_testx)
  }
  z=as.matrix(z)
  cz = z
  n=dim(z)[1]
  q = dim(z)[2]
  cvecx = vec_x
  nclass <- as.integer(length(unique(y)))
  if (is.centered==FALSE){
    for (i in 1:nclass){
      if (q>1) {cz[y == i,] = sweep(z[y == i,],2,colMeans(z[y == i,]))}
      else {cz[y==i]=z[y==i]-mean(z[y==i])}
      cvecx[y == i,] = sweep(vec_x[y == i,],2,colMeans(vec_x[y == i,]))
    }
  }
  c = solve(t(cz) %*% cz) %*% t(cz)
  c = c %*% cvecx

  vec_xres = vec_x-z%*%c

  xres = array(list(),n)
  for (i in 1:n){
    xres[[i]] = array(vec_xres[i,],dim=dimen)
  }
 
  
  if (!is.null(testx)){
    vec_testxres = vec_testx-testz %*% c
    testxres = array(list(),tesize)
    for (i in 1:tesize){
      testxres[[i]] = array(vec_testxres[i,],dim=dimen)
    } 
  }else{
    testxres=NULL
  }
  
  muz = matrix(0,nrow=q,ncol=(nclass-1))
  for (i in 2:nclass){
    if (q>1){muz[,(i-1)] = apply(z[y==i,],2,mean)-apply(z[y==1,],2,mean)}
    else {muz[i-1]=mean(z[y==i])-mean(z[y==1])}
  }
  gamma = solve(cov(z)) %*% muz
  
  outlist=list(gamma=gamma,xres=xres,testxres=testxres)
  outlist
}