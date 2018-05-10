\name{cv.catch}
\alias{cv.catch}
\title{
Cross-validation for CATCH
}
\description{
Performs k-fold cross validation for CATCH and returns the best tuning parameter \eqn{\lambda} in the user-specified or automatically generated choices.
}
\usage{
cv.catch(x, z = NULL, y, nfolds = 5, lambda = NULL, 
lambda.opt = "min",...)
}
\arguments{
  \item{x}{Input tensor or matrix list of length \eqn{N}{N}, where \eqn{N}{N} is the number of observations. Each element of the list is a tensor or matrix. The order of tensor can be any number and not limited to three.}
  \item{z}{Input covariate matrix of dimension \eqn{N \times q}{N*q}, where \eqn{q<N}{q<N}. \code{z} can be omitted if covariate is absent. }
  \item{y}{Class label. For \code{K} class problems, \code{y} takes values in \eqn{\{1,\cdots,\code{K}\}}.}
  \item{nfolds}{Number of folds. Default value is \code{5}.}
  \item{lambda}{User-specified \code{lambda} sequence for cross validation. If not specified, the algorithm will generate a sequence of \code{lambda}s based on all data and cross validate on the sequence.}
  \item{lambda.opt}{The optimal criteria when multiple elements in \code{lambda} return the same minimum classification error. "\code{min}" will return the smallest \code{lambda} with minimum cross validation error. "\code{max}" will return the largest \code{lambda} with the minimum cross validation error.}
  \item{\dots}{Other arguments that can be passed to \code{\link{catch}}.}
}
\details{
The function \code{\link{cv.catch}} runs function \code{\link{catch}} \code{nfolds+1} times. The first one fits model on all data. If \code{lambda} is specified, it will check if all \code{lambda} satisfies the constraints of \code{dfmax} and \code{pmax} in \code{\link{catch}}. If not, a \code{lambda} sequence will be generated according to \code{lambda.factor} in \code{\link{catch}}. Then the rest \code{nfolds} many replicates will fit model on \code{nfolds-1} many folds data and predict on the omitted fold, repectively. Return the \code{lambda} with minimum average cross validation error and the largest \code{lambda} within one standard error of the minimum.
}
\value{
  \item{lambda}{The actual \code{lambda} sequence used. The user specified sequence or automatically generated sequence could be truncated by constraints on \code{dfmax} and \code{pmax}.}
  \item{cvm}{The mean of cross validation errors for each \code{lambda}.}
  \item{cvsd}{The standard error of cross validaiton errors for each \code{lambda}.}
  \item{lambda.min}{The \code{lambda} with minimum cross validation error. If \code{lambda.opt} is \code{min}, then returns the smallest \code{lambda} with minimum cross validation error. If \code{lambda.opt} is \code{max}, then returns the largest\code{lambda} with minimum cross validation error.}
  \item{lambda.1se}{The largest \code{lambda} with cross validation error within one standard error of the minimum.}
  \item{catch.fit}{The fitted \code{catchobj} object.}
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
vec_x <- matrix(rnorm(n*nvars), nrow=n, ncol=nvars)
vec_x[1:10,] <- vec_x[1:10,]+2
z <- matrix(rnorm(n*2),nrow=n,ncol=2)
z[1:10,] <- z[1:10,]+0.5
y <- c(rep(1,10),rep(2,10))
for (i in 1:n){
  x[[i]] <- array(vec_x[i,], dim=c(p,p,p))
}
objcv <- cv.catch(x, z, y=y)
}

