\name{ISpmd}
\alias{ISpmd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return the best parameter based on IS.pmd Criterion
%%  ~~function to do ... ~~
}
\description{
Return the best parameter based on IS.pmd Criterion
}
\usage{
ISpmd(res.pmd, Yn1, lambda0, l, r)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res.pmd}{
%%     ~~Describe \code{res.pmd} here~~
}
  \item{Yn1}{
%%     ~~Describe \code{Yn1} here~~
}
  \item{lambda0}{
%%     ~~Describe \code{lambda0} here~~
}
  \item{l}{
%%     ~~Describe \code{l} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
}
\author{
Ruiping Liu
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (res.pmd, Yn1, lambda0, l, r)
{
    ISpmduv <- matrix(NA, l, 4)
    num.zero.pmduv <- matrix(NA, l, 2)
    Va.u <- matrix(NA, l, 1)
    Va.v <- matrix(NA, l, 1)
    Vs <- rep(NA, l)
    Vo <- lambda0
    for (i in 1:l) {
        A.v = Yn1 \%*\% res.pmd[[i]]$v
        A.u = t(Yn1) \%*\% res.pmd[[i]]$u
        qrresult.v <- qr(A.v)
        qrresult.u <- qr(A.u)
        Qreslt.u <- qr.Q(qrresult.u)
        Qreslt.v <- qr.Q(qrresult.v)
        if (r == 1) {
            Vs[i] <- (res.pmd[[i]]$d[1])^2
            Va.v[i] <- as.numeric(crossprod(A.v[, 1:r])) * as.numeric(cor(Qreslt.v[,
                1:r], A.v[, 1:r]))^2
            Va.u[i] <- as.numeric(crossprod(A.u[, 1:r])) * as.numeric(cor(Qreslt.u[,
                1:r], A.u[, 1:r]))^2
        }
        else {
            Vs[i] <- sum((res.pmd[[i]]$d[1:r])^2)
            Va.v[i] <- as.numeric(t(as.matrix(diag(crossprod(A.v[,
                1:r])))) \%*\% as.matrix(diag(cor(Qreslt.v[, 1:r],
                A.v[, 1:r])))^2)
            Va.u[i] <- as.numeric(t(as.matrix(diag(crossprod(A.u[,
                1:r])))) \%*\% as.matrix(diag(cor(Qreslt.u[, 1:r],
                A.u[, 1:r])))^2)
        }
        num.zero.pmduv[i, 1] <- length(which(abs(res.pmd[[i]]$u[,
            1:r]) < 1e-06))
        num.zero.pmduv[i, 2] <- length(which(abs(res.pmd[[i]]$v[,
            1:r]) < 1e-06))
        ISpmduv[i, 1] <- (Vs[i] * Va.u[i]/(Vo^2)) * (num.zero.pmduv[i,
            1]/(n * r))
        ISpmduv[i, 2] <- Vs[i] * Va.v[i]/(Vo^2) * (num.zero.pmduv[i,
            2]/(p * r))
        ISpmduv[i, 3] <- sum(ISpmduv[i, 1:2])
    }
    ISpmduv[, 4] <- ISpmduv[, 3]
    i0.IS = which.max(ISpmduv[, 3])
    list(ISpmduv = ISpmduv, Vs = Vs, Va.u = Va.u, Va.v = Va.v,
        num.zero.pmduv = num.zero.pmduv, i0.IS = i0.IS)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
