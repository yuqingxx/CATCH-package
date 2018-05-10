\name{catch_matrix}
\alias{catch_matrix}

\title{
Fit a CATCH model for matrix and predict categorical response.
}
\description{
Fits a classifier for matrix data. \code{catch_matrix} is a special case of \code{\link{catch}} when each observation \eqn{\mathbf{X}_i} is a matrix. Different from \code{\link{catch}} takes list as input, data need to be formed in an array to call the function (see arguments). The function will perform prediction as well.
}
\usage{
catch_matrix(x, z = NULL, y, testx = NULL, testz = NULL, \dots)
}
\arguments{

  \item{x}{Input matrix array. The array should be organized with dimension \eqn{p_1 \times p_2 \times N}{p1*p2*N, where \eqn{N} is the number of observations}.
}
  \item{z}{Input covariate matrix of dimension \eqn{N \times q}{N*q}, where \eqn{q<N}{q<N}. \code{z} can be omitted if covariate is absent. }
  \item{y}{Class label. For \code{K} class problems, \code{y} takes values in \eqn{\{1,\cdots,\code{K}\}}.}
  \item{testx}{Input testing matrix array. When \code{testx} is not provided, the function will only fit model. When \code{testx} is provided, the function will predict response on \code{testx} as well.}
  \item{testz}{Input testing covariate matrix. Can be omitted if there is no covariate.}
  \item{\dots}{Other arguments that can be passed to \code{\link{catch}}.}

}
\details{
The function fits a matrix classifier as a special case of \code{\link{catch}}. The fitted model and predictions should be identical to \code{\link{catch}} when matrix data is provided. Input matrix should be organized as three-way array where sample size is the last dimension. If the matrix is organized in a list, users can either reorganize it or use \code{\link{catch}} directly to fit model, which takes a matrix or tensor list as input and has the same output as \code{\link{catch_matrix}}.
}

\value{
\item{beta}{Output variable coefficients for each \code{lambda}. \code{beta} is a list of length being the number of \code{lambda}s. Each element of \code{beta} is a matrix of dimension \eqn{(p_1\times p_2)\times (nclass-1)}{(p_1\times p_2)*(nclass-1)}.}
\item{df}{The number of nonzero variables for each value in sequence \code{lambda}.}
\item{dim}{Dimension of coefficient array.}
\item{lambda}{The actual \code{lambda} sequence used. The user specified sequence or automatically generated sequence could be truncated by constraints on \code{dfmax} and \code{pmax}.}
\item{obj}{Objective function value for each value in sequence \code{lambda}.}
\item{x}{The matrix list after adjustment in training data. If covariate is absent, this is the original input matrix.}
\item{y}{Class label in training data.}
\item{npasses}{Total number of iterations.}
\item{jerr}{Error flag.}
\item{sigma}{Estimated covariance matrix on each mode. \code{sigma} is a list with the \code{i}th element being covariance matrix on \code{i}th mode.}
\item{delta}{Estimated delta matrix \eqn{(vec(\widehat{\boldsymbol{\mu}}_2-\widehat{\boldsymbol{\mu}}_1),\cdots,vec(\widehat{\boldsymbol{\mu}}_K-\widehat{\boldsymbol{\mu}}_1))}.}
\item{mu}{Estimated mean array.}
\item{prior}{Prior proportion of observations in each class.}
\item{call}{The call that produces this object.}
\item{pred}{Predicted categorical response for each value in sequence \code{lambda} when \code{testx} is provided.}
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
#without prediction
n <- 20
p <- 4
k <- 2
nvars <- p*p
x=array(rnorm(n*nvars),dim=c(p,p,n))
x[,,11:20]=x[,,11:20]+0.3
z <- matrix(rnorm(n*2), nrow=n, ncol=2)
z[1:10,] <- z[1:10,]+0.5
y <- c(rep(1,10),rep(2,10))
obj <- catch_matrix(x,z,y=y)
}

