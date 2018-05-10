\name{csa}
\alias{csa}
\docType{data}
\title{
Colorimetric sensor array (CSA) data}
\description{
A dataset collected from a series of CSA experiments to identify volatile chemical toxicants (VCT). Chemical dyes were exposed to VCT under different concentration conditions and colors of dyes were recorded to identify the class of VCT. There are two concentration conditions PEL (permissible exposure level) and IDLH (immediately dangerous to life of health). 

}
\usage{
data(csa)
}

\format{
  Two lists, \emph{PEL} and \emph{IDLH}, and a numeric vector \emph{y}. Each list contains 147 matrices of dimension \eqn{36 \times 3}{36*3}. 
  \describe{
    \item{\code{PEL}}{A list of matrices containing the observations after exposure at PEL.}
    \item{\code{IDLH}}{A list of matrices containing the observations after exposure at IDLH level.}
    \item{\code{y}}{Class label ranging from 1 to 21.}
  
  }
}


\details{
This dataset is provided in the Supplementary matrial of Zhong (2015). In each concentration case, there are 147 observations and 21 classes. We reorganize the data into a list to be directly called by \code{catch}. For matrices in the list, each row represents a dye and the three columns correspond to red, green and blue.
}
\source{
Wenxuan Zhong and Kenneth S. Suslick (2015). "Matrix discriminant analysis with application to colorimetric sensor array data" \emph{Technometrics} \bold{57}(4), 524--534.
}
\keyword{datasets}
