#' Title Return the svd results for SCA by BiOPD
#'
#' @param X data matrix (Frequency data)
#' @param para.pmd parameter range for type 'tau'
#' @param para.uv parameter range for type 'grid'
#' @param v0 the initial v
#' @param cK used when Criterion=='Nonzeros'
#' @param Sparse.type 'tau' or 'grid'
#' @param decomposition.type 'BiOPD'
#' @param Criterion.type default=='IS'
#'
#' @return BiOPD results
#' @export
#'
#' @examples SCApmd(X,para.pmd,para.uv,v0,cK,Sparse.type,decomposition.type,Criterion.type)

SCApmd <- function(
  X, para.pmd, para.uv,
  v0,
  cK,
  Sparse.type,
  decomposition.type,
  Criterion.type){
  # -- Obtain matrix Yn1: sqrt(solve(Dn))%*%(X-r%*%t(c))%*%sqrt(solve(Dp)) # with the 'centered' X
  re.CA <- CA(X);
  n <- dim(X)[1]; p <- dim(X)[2] # dim of data matrix
  Fd <- prop.table(X) # X0/sum(X0); the frequency data matrix;
  #
  ObtainYn <- function(X){
    X <- X/sum(X)
    c <- as.matrix(apply(X,2,sum)) # col sum
    r <- as.matrix(apply(X,1,sum)) # row sum
    Dp <- diag(as.numeric(c)); Dn <- diag(as.numeric(r))
    Y <- X- r%*% t(c)                              # 'centered' frequency data table, to delete the eigenvalue=1
    Yn1 <- sqrt(solve(Dn)) %*%Y%*% sqrt(solve(Dp)) # with the 'centered' X
    Yn <- sqrt(solve(Dn)) %*%X%*% sqrt(solve(Dp))
    Yn1 <- as.matrix(Yn1)
    list(
      r=r,
      c=c,
      Dn=Dn,
      Dp=Dp,
      Yn=Yn,
      Yn1=Yn1)
  }
  tempYn <- ObtainYn(Fd);  Dn <- tempYn$Dn; Dp <- tempYn$Dp
  # r <- ObtainYn(Fd)$r; c <- ObtainYn(Fd)$c
  Yn1 <- tempYn$Yn1;
  if (Sparse.type == 'tau'){
    # Dim1
    res.pmd.dim1 <- NULL;
    res.pmd.sumabs.dim1 <- NULL;
    for (i in 1:length(para.pmd)){
      res.pmd.sumabs.dim1[[i]] <- PMD(Yn1, type="standard",
                                      sumabs=para.pmd[i],
                                      sumabsu=NULL, sumabsv=NULL, lambda=NULL,
                                      niter=20, K=1, #min(dim(Yn1)),
                                      v=v0, trace=F, #
                                      center=FALSE, chrom=NULL, rnames=NULL, cnames=NULL,
                                      upos=FALSE,
                                      uneg=FALSE, vpos=FALSE, vneg=FALSE)
    }
    res.pmd.dim1 = res.pmd.sumabs.dim1
  }
  else if (Sparse.type == 'grid'){ # (sumabsu, sumabsv)
    luv = dim(para.uv)[1];
    # Dim1
    res.pmd.sumabsuv.dim1 <- NULL; #len=lu*lv
    for (i in 1:luv){
      res.pmd.sumabsuv.dim1 <- c(res.pmd.sumabsuv.dim1,
        list(PMD(Yn1, type="standard",sumabs=NULL,
                 sumabsu=
                   # para.pmdu[iu],
                   para.uv[i, 1],
                 sumabsv=
                   # para.pmdv[iv],
                   para.uv[i, 2],
                 lambda=NULL,
                 niter=20, K=1,
                 v=v0, # NULL,
                 trace=F,
                 center=F, chrom=NULL, rnames=NULL, cnames=NULL, upos=FALSE,
                 uneg=FALSE, vpos=FALSE, vneg=FALSE))
      )
      res.pmd.dim1 <- res.pmd.sumabsuv.dim1
    }
    res.pmd.dim1 <- res.pmd.sumabsuv.dim1
  }
  # Dim2
  temp.Criterion.dim1 <-
    Decomposition.Criterion.Dim1(Fd = X,
                                 Yn1,
                                 v0,
                                 Dn, Dp,
                                 re = re.CA,
                                 para.pmd,
                                 para.uv,
                                 cK,
                                 Sparse.type,
                                 Criterion.type,
                                 decomposition.type,
                                 res.pmd.dim1);
  i0.dim1 <- temp.Criterion.dim1$i0.dim1
  # Criterion.Dim2
  Criterion.BiOPD.dim2 <-
    Decomposition.Criterion.Dim2(Fd = X,
                                 Yn1, v0,
                                 Dn, Dp,
                                 re = re.CA,
                                 cK,
                                 Sparse.type,
                                 Criterion.type,
                                 decomposition.type,
                                 temp.Decomposition.dim1 = temp.Criterion.dim1);
  i0.dim2 <- Criterion.BiOPD.dim2$i0.dim2;
  # i0.dim2.union <- Criterion.BiOPD.dim2$i0.dim2.union
  #PEV
  PEV <- NULL;
  PEV$PEV0 <- c(res.pmd.dim1[[i0.dim1]]$d^2,
                temp.Criterion.dim1$res.BiOPD.dim2[[i0.dim2]]$d^2)/sum(re.CA$eig[,1]);
  # PEV$PEV.union <- c(res.pmd.dim1[[i0.dim1]]$d^2,
  #                    temp.Criterion.dim1$res.BiOPD.dim2[[i0.dim2.union]]$d^2)/sum(re$eig[,1])
  list(
    res.pmd.dim1 = res.pmd.dim1,
    temp.Criterion.dim1 = temp.Criterion.dim1,
    res.pmd.dim2 = temp.Criterion.dim1$res.pmd.dim2,
    res.BiOPD.dim2 = temp.Criterion.dim1$res.BiOPD.dim2,
    Criterion.BiOPD = Criterion.BiOPD.dim2,
    PEV = PEV
  )
}
