\name{SCA}
\alias{SCA}

\title{SparseCA}

\usage{
SCA(X0,
    para.pmd,
    para.uv,
    Sparse.type,
    decomposition.type,
    Criterion.type,
    index.Col.Nonzero, #
    cK,
    R)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X0}{
%%     ~~Describe \code{res.pmd} here~~
}
  \item{para.pmd}{
%%     ~~Describe \code{Yn1} here~~
}
  \item{para.uv}{
%%     ~~Describe \code{lambda0} here~~
}
  \item{Sparse.type}{
%%     ~~Describe \code{l} here~~
}
  \item{decomposition.type}{
%%     ~~Describe \code{r} here~~
}
  \item{Criterion.type}{
%%     ~~Describe \code{r} here~~
}
  \item{index.Col.Nonzero}{
%%     ~~Describe \code{r} here~~
}
  \item{cK}{
%%     ~~Describe \code{r} here~~
}
  \item{R}{
%%     ~~Describe \code{r} here~~
}
}
\description{
Sparse correspondence analysis (CA) for large contingency tables. Sparse CA remedies to the difficulty of interpreting CA results when the size of the table is large. Since CA is a double weighted PCA or a weighted generalized SVD, we adapt here known sparse versions of these methods with specific developments to obtain orthogonal solutions and a search method for the optimal sparsity level. We distinguish two cases where sparsity is asked for both rows and columns, or only for one set.
}
\references{
Ruiping Liu, Ndeye Niang Keita, Gilbert Saporta, Huiwen Wang. (2022). Sparse Correspondence Analysis for Contingency Tables. Advances in Data Analysis and Classification. https://doi.org/10.1007/s11634-022-00531-5
}
\author{
Ruiping Liu
}
\examples{
load("~/Data/ColorSound.rda")
SCA(X0 = ColorSound,
    para.pmd = 1, para.uv=cbind(seq(1,sqrt(10),length=10), seq(1,sqrt(9),length=10)),
    Sparse.type='tau',
    decomposition.type='BiOPD',
    Criterion.type='IS',
    index.Col.Nonzero=NULL,
    cK=NULL, R=2)
}

