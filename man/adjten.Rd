\name{adjten}
\alias{adjten}
\title{Adjust tensor for covariates.
}
\description{
Adjusts tensor with respect to covariates to achieve a more accurate performance. Tensor depends on the covariates through a linear regression model. The function returns the coefficients of covariates in regression and adjusted tensor list for further classifier modeling. It estimates coefficients based on training data, and then adjusts training tensor. When testing data is provided, the function will automatically adjust testing data by learned coefficients as well.
}
\usage{
adjten(x, z, y, testx = NULL, testz = NULL, is.centered = FALSE)
}
\arguments{
  \item{x}{Input tensor or matrix list of length \eqn{N}{N}, where \eqn{N}{N} is the number of observations. Each element of the list is a tensor or matrix. The order of tensor can be any integer not less than 2.}
  \item{z}{Input covariate matrix of dimension \eqn{N \times q}{N*q}, where \eqn{q<N}{q<N}. Each row of \code{z} is an observation.}
  \item{y}{Class label vector of dimention \eqn{N\times 1}{N*1}. For \code{K} class problems, \code{y} takes values in \eqn{\{1,\cdots,\code{K}\}}.}
  \item{testx}{Input testing tensor or matrix list. Each element of the list is a test case. When \code{testx} is not provided, the function will only adjust training data.}
  \item{testz}{Input testing covariate matrix with each row being an observation.}
  \item{is.centered}{Indicates whether the input tensor and covariates have already been centered by their within class mean or not. If \code{is.centered} is FALSE, the function \code{adjten} will center data by class. If \code{is.centered} is TRUE, the function will skip the centering step.}
}
\details{
The model CATCH assumes the linear relationship bewteen covariates and tensor as
\deqn{\mathbf{X}=\boldsymbol{\mu}_k+\boldsymbol{\alpha}\overline{\times}_{M+1}\mathbf{Z}+\mathbf{E},} where \eqn{\boldsymbol{\alpha}} is the matrix of estimated coefficient of covariates.
The function removes the effects of covariates on response variable through tensor and obtain \eqn{\mathbf{X}-\boldsymbol{\alpha}\overline{\times}_{M+1}\mathbf{Z}} as adjusted tensor to fit tensor discriminant analysis model.

In estimating \eqn{\boldsymbol{\alpha}}, which is the \code{alpha} in the package, \code{\link{adjten}} first centers both tensor and covariates within their individual classes, then performs tensor response regression which regresses \eqn{{\mathbf{X}}} on \eqn{{\mathbf{Z}}}. 
}
\value{
  \item{gamma}{The estimated coefficients of covariates to plug in classifier. \code{gamma} is the \eqn{\boldsymbol{\gamma}_k} defined function \code{\link{catch}} of dimension \eqn{q\times (K-1)}{q*(K-1)}, where \code{q} is the size of covariates and \code{K} is the number of classes.}
  \item{xres}{Adjusted training tensor list \eqn{\mathbf{X}-\boldsymbol{\alpha}\overline{\times}_{M+1}\mathbf{Z}} after adjusting for covariates. The effect of the covariate is removed.}
  \item{testxres}{Adjusted testing tensor list \eqn{\mathbf{X}-\boldsymbol{\alpha}\overline{\times}_{M+1}\mathbf{Z}} after adjusting for covariates. The effect of the covariate is removed.}
}
\references{
Pan, Y., Mai, Q., and Zhang, X. (2018) \emph{Covariate-Adjusted Tensor Classification in High-Dimensions, Manuscript.}
}
\author{
Yuqing Pan, Qing Mai, Xin Zhang
}

\seealso{
\code{\link{catch}}
}
\examples{
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
obj <- adjten(x, z, y)
}


