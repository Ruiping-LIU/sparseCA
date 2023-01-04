#' Title Return the best parameter based on IS.pmd Criterion
#'
#' @param res.pmd for Case.tau (sumabs) & Case.grid(sumabsu, sumabsv)
#' @param lambda0
#' @param l length(para.pmd)
#' @param r Default=1 (for Dim by Dim), r be the requied dimention
#'
#' @return Best parameter based on IS.pmd Criterion
#' @export
#'
ISpmd <- function(res.pmd,
                   Yn1,
                   lambda0,
                   l,
                   r
){
  n <- dim(Yn1)[1]; p <- dim(Yn1)[2]
  IS.pmduv <- matrix(NA, l, 4)
  num.zero.pmduv <- matrix(NA, l, 2) #U & V
  Va.u <- matrix(NA, l, 1)
  Va.v <- matrix(NA, l, 1)
  Vs <- rep(NA, l)
  Vo <-
    lambda0 #Dim by Dim
  for (i in 1:l) {
    # QR decomposition for Components A
    A.v = Yn1 %*% res.pmd[[i]]$v
    A.u = t(Yn1) %*% res.pmd[[i]]$u

    qrresult.v <- qr(A.v)  # QR-decomposition
    qrresult.u <- qr(A.u)

    Qreslt.u <- qr.Q(qrresult.u)  #
    Qreslt.v <- qr.Q(qrresult.v)
    # Rreslt <- qr.R(qrresult)  #
    # Xreslt <- qr.X(qrresult)  #
    if (r == 1) {
      Vs[i] <- (res.pmd[[i]]$d[1])^2
      Va.v[i] <-
        as.numeric(crossprod(A.v[, 1:r])) * as.numeric(cor(Qreslt.v[, 1:r], A.v[, 1:r])) ^ 2
      Va.u[i] <-
        as.numeric(crossprod(A.u[, 1:r]))* as.numeric(cor(Qreslt.u[, 1:r], A.u[, 1:r])) ^ 2
    }
    else {
      Vs[i] <- sum((res.pmd[[i]]$d[1:r]) ^ 2)
      Va.v[i] <-
        as.numeric(t(as.matrix(diag(crossprod(
          A.v[, 1:r]
        )))) %*% as.matrix(diag(cor(Qreslt.v[, 1:r], A.v[, 1:r]))) ^ 2)
      Va.u[i] <-
        as.numeric(t(as.matrix(diag(crossprod(
          A.u[, 1:r]
        )))) %*% as.matrix(diag(cor(Qreslt.u[, 1:r], A.u[, 1:r]))) ^ 2)
    }
    ##
    num.zero.pmduv[i, 1] <-
      length(which(abs(res.pmd[[i]]$u[, 1:r]) < 1e-6)) #for U
    num.zero.pmduv[i, 2] <-
      length(which(abs(res.pmd[[i]]$v[, 1:r]) < 1e-6)) #for V
    IS.pmduv[i, 1] <-
      # (sum((res.pmd[[i]]$d[1:r])^2)/sum(re$eig[1:r,1]))^2*(num.zero.pmduv[i,1]/(n*r)) #for U
      # (sum((res.pmd[[i]]$d[1:r])^2)/sum(re$eig[1:r,1]))*(num.zero.pmduv[i,1]/(n*r))
      (Vs[i] * Va.u[i] / (Vo ^ 2)) * (num.zero.pmduv[i, 1] / (n * r))

    IS.pmduv[i, 2] <-
      # (sum((res.pmd[[i]]$d[1:r])^2)/sum(re$eig[1:r,1]))^2*(num.zero.pmduv[i,2]/(p*r)) #for V
      # (sum((res.pmd[[i]]$d[1:r])^2)/sum(re$eig[1:r,1]))*(num.zero.pmduv[i,2]/(p*r))
      Vs[i] * Va.v[i] / (Vo ^ 2) * (num.zero.pmduv[i, 2] / (p * r))
    IS.pmduv[i, 3] <-
      sum(IS.pmduv[i, 1:2])
  }
  IS.pmduv[, 4] <- IS.pmduv[, 3]
  # for(i in 1:l){
  #   if (n-num.zero.pmduv[i,1] <2 | p-num.zero.pmduv[i,2] <2){
  #     IS.pmduv[i, 3] = 0
  #   }
  # }
  i0.IS = which.max(IS.pmduv[, 3])
  list(IS.pmduv = IS.pmduv, Vs = Vs, Va.u = Va.u, Va.v = Va.v,
       num.zero.pmduv = num.zero.pmduv, i0.IS = i0.IS)
}
