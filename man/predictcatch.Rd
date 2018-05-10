\name{predict.catch}
\alias{predict.catch}
\title{
Predict categorical responses.
}
\description{
Predict categorical responses on new data given the fitted model input.
}
\usage{
\method{predict}{catch}(object, newx, z = NULL, ztest = NULL, gamma = NULL,...)
}
\arguments{
  \item{object}{Input \code{catchobj} class object as fitted model.}
  \item{newx}{Input adjusted testing tensor or matrix list. Each element of the list is a tensor. The tensor should of the same dimension as training data.}
  \item{z}{Input training covariates matrix. \code{z} can be omitted if there is no covariate. }
  \item{ztest}{Input testing covariates matrix. \code{ztest} can be omitted if there is no covariate.}
  \item{gamma}{Coefficients of covariates obtained from \code{\link{adjten}}. \code{gamma} is \code{NULL} if there is no covariate.}
  \item{\dots}{Other arguments that can be passed to \code{predict}.}
}
\details{
The function fits LDA model on selected discriminant vectors. Call \code{predict} or \code{predict.catch} to perform predictions.

There are two ways to make predictions. One way is to directly predict at the same time as fitting model by \code{\link{catch}} since \code{\link{predict.catch}} has already been embedded in \code{\link{catch}} and it will predicts response when testing data is provided. The other way is to first use \code{\link{adjten}} to adjuste tensor and \code{\link{catch}} to fit model. \code{\link{predict.catch}} will take the input adjusted tensor list \code{newx}, covariate coefficient \code{gamma} from \code{\link{adjten}} and the fitted model from \code{\link{catch}} to perform prediction. The prediction is identical to providing \code{\link{catch}} testing data.
}
\value{
\item{pred}{Predicted response of \code{newx} for each \code{lambda} in model \code{object}.}
}
\references{
Pan, Y., Mai, Q., and Zhang, X. (2018) \emph{Covariate-Adjusted Tensor Classification in High-Dimensions, Manuscript.}
}
\author{
Yuqing Pan, Qing Mai, Xin Zhang
}
\seealso{
\code{\link{catch}}, \code{\link{adjten}}
}
\examples{
#generate training data
n <- 20
p <- 4
k <- 2
nvars <- p*p*p
x <- array(list(),n)
vec_x <- matrix(rnorm(n*nvars),nrow=n,ncol=nvars)
vec_x[1:10,] <- vec_x[1:10,]+2
z <- matrix(rnorm(n*2),nrow=n,ncol=2)
z[1:10,] <- z[1:10,]+0.5
y <- c(rep(1,10),rep(2,10))
for (i in 1:n){
  x[[i]] <- array(vec_x[i,],dim=c(p,p,p))
}

#generate testing data
newx <- array(list(),n)
vec_newx <- matrix(rnorm(n*nvars),nrow=n,ncol=nvars)
vec_newx[1:10,] <- vec_newx[1:10,]+2
newz <- matrix(rnorm(n*2),nrow=n,ncol=2)
newz[1:10,] <- newz[1:10,]+0.5
for (i in 1:n){
  newx[[i]] <- array(vec_newx[i,],dim=c(p,p,p))
}

#Make adjustment and fit model
obj <- adjten(x, z, y, newx, newz)
fit <- catch(x, z, y)
#Predict
pred <- predict.catch(fit, obj$testxres, z, newz, obj$gamma)

#The adjusting, fitting model and predicting step can also be completed
#by one command.
pred <- catch(x, z, y, newx, newz)$pred
}
