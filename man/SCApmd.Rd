\name{SCApmd}
\alias{SCApmd}

\title{Return the svd results of SCA by BiOPD}
\usage{
SCApmd(X,
        para.pmd,
        para.uv,
        v0,
        cK,
        Sparse.type,
        decomposition.type,
        Criterion.type)
}
\arguments{
  \item{X}{
%%     ~~Describe \code{res.pmd} here~~
}
  \item{para.pmd}{
%%     ~~Describe \code{Yn1} here~~
}
  \item{para.uv}{
%%     ~~Describe \code{lambda0} here~~
}
  \item{v0}{
%%     ~~Describe \code{r} here~~
}
  \item{cK}{
%%     ~~Describe \code{r} here~~
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
}
\details{
Function used in function 'SCA'
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{res.pmd.dim1}{res.pmd.dim1}
  \item{temp.Criterion.dim1}{temp.Criterion.dim1}
  \item{res.pmd.dim2}{temp.Criterion.dim1$res.pmd.dim2}
  \item{res.BiOPD.dim2}{temp.Criterion.dim1$res.BiOPD.dim2}
  \item{Criterion.BiOPD}{Criterion.BiOPD.dim2}
  \item{PEV}{PEV}
}
\description{
  Return the svd results of SCA by BiOPD
}
\note{
can be further optimized
}
\references{
Ruiping Liu, Ndeye Niang Keita, Gilbert Saporta, Huiwen Wang. (2022). Sparse Correspondence Analysis for Contingency Tables. Advances in Data Analysis and Classification. https://doi.org/10.1007/s11634-022-00531-5
}
\examples{
load("~/Data/ColorSound.rda")
X0 <- ColorSound
library(PMA, FactoMineR)
X <- X0/sum(X0);
c <- as.matrix(apply(X,2,sum)) # col sum
r <- as.matrix(apply(X,1,sum)) # row sum
Dp <- diag(as.numeric(c)); Dn <- diag(as.numeric(r))
Yn1 <- sqrt(solve(Dn)) %*% (X- r %*% t(c)) %*% sqrt(solve(Dp)) ;
Yn1 <- as.matrix(Yn1)
res.pmd0 <- PMD(Yn1, type="standard",sumabs=1,
                        sumabsu=NULL, sumabsv=NULL, lambda=NULL,
                        niter=20, K=1, v=NULL, trace=FALSE,
                        center=FALSE, chrom=NULL, rnames=rownames(X0), cnames=colnames(X0),
                        upos=FALSE, uneg=FALSE, vpos=FALSE, vneg=FALSE);
v0 <- res.pmd0$v.init;
v0 = as.vector(v0);
SCApmd(X,
       para.pmd = seq(1/sqrt(min(dim(X))), 1, len=10),
       para.uv = cbind(seq(1,sqrt(10),length=10), seq(1,sqrt(9),length=10)),
       v0,
       cK=NULL,
       Sparse.type ='tau',
       decomposition.type ='BiOPD',
       Criterion.type ='IS')
}
