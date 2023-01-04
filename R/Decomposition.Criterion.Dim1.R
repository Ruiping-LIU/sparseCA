#' Title Decomposition for Dim1
#'
#' @param Fd Frequency data
#' @param Dn diag matrix of RowSum
#' @param Dp diag matrix of ColumnSum
#' @param re results of CA
#' @param para.pmd parameter range for 'tau' type
#' @param para.uv parameter range for 'grid' type
#' @param cK cumulative Nonzeros. For Criterion=='Nonzeros'
#' @param Sparse.type tau or grid
#' @param Criterion.type default=IS
#' @param decomposition.type 'BiOPD'
#' @param res.pmd.dim1 pmd results for Dim1
#'
#' @return Results by Decomposition for Dim1
#' @export
#'
#' @examples Decomposition.Criterion.Dim1()
Decomposition.Criterion.Dim1 <- function(Fd,
                                         Yn1,
                                         v0,
                                         Dn, Dp,
                                         re,
                                         para.pmd,
                                         para.uv,
                                         cK,
                                         Sparse.type,
                                         Criterion.type,
                                         decomposition.type,
                                         res.pmd.dim1 # Eg. res.pmd.BothSide.dim1,or res.pmd.OneSide.dim1
){
  n <- dim(Fd)[1]; p <- dim(Fd)[2]
  r <- as.matrix(diag(Dn)); c<-as.matrix(diag(Dp))
  l = length(res.pmd.dim1); #length(para.pmd), length(para.pmduv[,1])
  Nonzero.pmduv <- matrix(NA, l, 2);
  k <- cK[1];
  #
  res.pmd.dim2 <- vector("list", l); #list for Dim2
  res.BiOPD.dim2 <- vector("list", l)
  #
  if(Sparse.type == 'tau'){# "BothSide"
    # res.pmd.dim1 = res.pmd.BothSide.dim1
    #Criterion
    if(Criterion.type == 'IS'){
      tm.IS.dim1 <- ISpmd(res.pmd = res.pmd.dim1,
                           Yn1,
                           lambda0 = re$eig[1,1],
                           l=l, r=1)
      i0.dim1 <- tm.IS.dim1$i0.IS;
      tm.Criterion.dim1 = tm.IS.dim1
    }
    #
    #deflation
    if(decomposition.type == 'PMD'){
      tm.dataMat.dim2 <-
        Yn1 - res.pmd.dim1[[i0.dim1]]$u %*%
        t(res.pmd.dim1[[i0.dim1]]$d *
            res.pmd.dim1[[i0.dim1]]$v) #
      #
      for (i in 1:l){
        # parameter: sumabs 50, 100 parameters
        res.pmd.dim2[[i]] <- PMD(tm.dataMat.dim2,
                                 type="standard",sumabs=para.pmd[i],
                                 sumabsu=NULL, sumabsv=NULL, lambda=NULL,
                                 niter=20, K=1, v = v0, trace=F, #v0[,1], v0[,1:K]
                                 center=F, chrom=NULL, rnames=NULL, cnames=NULL, upos=FALSE,
                                 uneg=FALSE, vpos=FALSE, vneg=FALSE)
        #
        res.BiOPD.dim2[[i]] <- NULL
      }
    }
    else if (decomposition.type == 'mPMD'){
      tm.dataMat.dim2 <-
        (diag(n)-res.pmd.dim1[[i0.dim1]]$u %*%
           t(res.pmd.dim1[[i0.dim1]]$u)) %*% Yn1 %*%
        (diag(p)-res.pmd.dim1[[i0.dim1]]$v %*%
           t(res.pmd.dim1[[i0.dim1]]$v))
      #
      for (i in 1:l){
        # parameter: sumabs 50, 100 parameters
        res.pmd.dim2[[i]] <- PMD(tm.dataMat.dim2,
                                 type="standard",sumabs=para.pmd[i],
                                 sumabsu=NULL, sumabsv=NULL, lambda=NULL,
                                 niter=20, K=1, v = v0, trace=F, #v0[,1], v0[,1:K]
                                 center=F, chrom=NULL, rnames=NULL, cnames=NULL, upos=FALSE,
                                 uneg=FALSE, vpos=FALSE, vneg=FALSE)
        #
      }
    }
    else if (decomposition.type == 'BiOPD'){
      q1 <- res.pmd.dim1[[i0.dim1]]$v
      r1 <- res.pmd.dim1[[i0.dim1]]$u
      B1 <- diag(p)-q1%*%t(q1)
      C1 <- diag(n)-r1%*%t(r1)
      #deflation
      tm.dataMat.dim2 <-
        (diag(n)-res.pmd.dim1[[i0.dim1]]$u %*%
           t(res.pmd.dim1[[i0.dim1]]$u)) %*% Yn1 %*%
        (diag(p)-res.pmd.dim1[[i0.dim1]]$v %*%
           t(res.pmd.dim1[[i0.dim1]]$v))
      #
      for (i in 1:l){
        # parameter: sumabs 50, 100 parameters
        res.pmd.dim2[[i]] <- PMD(tm.dataMat.dim2,
                                 type="standard",sumabs=para.pmd[i],
                                 sumabsu=NULL, sumabsv=NULL, lambda=NULL,
                                 niter=20, K=1, v = v0, trace=F, #v0[,1], v0[,1:K]
                                 center=F, chrom=NULL, rnames=NULL, cnames=NULL, upos=FALSE,
                                 uneg=FALSE, vpos=FALSE, vneg=FALSE)
        #
        # OPD (Orthogonalize)
        q2 <- B1%*%res.pmd.dim2[[i]]$v;
        r2 <- C1%*%res.pmd.dim2[[i]]$u;
        #
        a <- NULL;
        #
        a$v <- q2/as.numeric(sqrt(crossprod(q2)));
        a$u <- r2/as.numeric(sqrt(crossprod(r2)));
        a$d <- res.pmd.dim2[[i]]$d;
        res.BiOPD.dim2[[i]] <- a
        #
        # rm(q2, r2, a)
        rm(q2);rm(r2); rm(a)
      }
    }
    #
  }
  else if(Sparse.type == 'grid'){
    # res.pmd.dim1 = res.pmd.OneSide.dim1
    #Criterion
    if(Criterion.type == 'IS'){
      tm.IS.dim1 <- ISpmd(res.pmd = res.pmd.dim1,
                           lambda0 = re$eig[1,1],
                           l=l, r=1)
      i0.dim1 <- tm.IS.dim1$i0.IS;
      tm.Criterion.dim1 = tm.IS.dim1
    }
    else if(Criterion.type == 'Nonzeros'){
      for(i in 1:l){
        Nonzero.pmduv[i, 1] <-  sum(abs(res.pmd.dim1[[i]]$u[, 1]) != 0)#length(which(abs(res.pmd.dim1[[i]]$u[, 1]) > 1e-6)); #for U
        Nonzero.pmduv[i, 2] <-  sum(abs(res.pmd.dim1[[i]]$v[, 1]) != 0)#length(which(abs(res.pmd.dim1[[i]]$v[, 1]) > 1e-6)); #for V
      }
      i0.dim1 <- which.min(abs(Nonzero.pmduv[,2]-k));
      tm.Criterion.dim1 <- Nonzero.pmduv
    }
    #
    #deflation
    if(decomposition.type == 'PMD'){
      tm.dataMat.dim2 <-
        Yn1 - res.pmd.dim1[[i0.dim1]]$u %*%
        t(res.pmd.dim1[[i0.dim1]]$d *
            res.pmd.dim1[[i0.dim1]]$v) #
    }
    else if (decomposition.type == 'mPMD'){
      tm.dataMat.dim2 <-
        (diag(n)-res.pmd.dim1[[i0.dim1]]$u %*%
           t(res.pmd.dim1[[i0.dim1]]$u))%*%Yn1%*%
        (diag(p)-res.pmd.dim1[[i0.dim1]]$v %*%
           t(res.pmd.dim1[[i0.dim1]]$v)) #
      #
      for (i in 1:l){
        # parameter: sumabs 50, 100 parameters
        res.pmd.dim2[[i]] <- PMD(tm.dataMat.dim2,
                                 type="standard",sumabs=NULL,
                                 sumabsu=para.uv[i,1],
                                 sumabsv=para.uv[i,2],
                                 lambda=NULL,
                                 niter=20, K=1, v=v0, trace=F, #v0[,1], v0[,1:K]
                                 center=F, chrom=NULL, rnames=NULL, cnames=NULL, upos=FALSE,
                                 uneg=FALSE, vpos=FALSE, vneg=FALSE)
      }
    }
    else if (decomposition.type == 'BiOPD'){
      #deflation
      tm.dataMat.dim2 <-
        (diag(n)-res.pmd.dim1[[i0.dim1]]$u %*%
           t(res.pmd.dim1[[i0.dim1]]$u))%*%Yn1%*%
        (diag(p)-res.pmd.dim1[[i0.dim1]]$v %*%
           t(res.pmd.dim1[[i0.dim1]]$v))
      #
    }
    q1 <- res.pmd.dim1[[i0.dim1]]$v
    r1 <- res.pmd.dim1[[i0.dim1]]$u
    B1 <- diag(p) - q1%*% t(q1)
    C1 <- diag(n) - r1%*% t(r1)
    #res.pmd.dim2   #
    for (i in 1:l){
      res.pmd.dim2[[i]] <- PMD(tm.dataMat.dim2,
                               type="standard",sumabs=NULL,
                               sumabsu=para.uv[i,1],
                               sumabsv=para.uv[i,2],
                               lambda=NULL,
                               niter=20, K=1, v=v0, trace=F, #v0[,1], v0[,1:K]
                               center=F, chrom=NULL, rnames=NULL, cnames=NULL, upos=FALSE,
                               uneg=FALSE, vpos=FALSE, vneg=FALSE)
      #
      # OPD (Orthogonalize)
      q2 <- B1%*%res.pmd.dim2[[i]]$v;
      r2 <- C1%*%res.pmd.dim2[[i]]$u;
      #
      a <- NULL;
      #
      a$v <- q2/as.numeric(sqrt(crossprod(q2)));
      a$u <- r2/as.numeric(sqrt(crossprod(r2)));
      a$d <- res.pmd.dim2[[i]]$d #*as.numeric(sqrt(crossprod(q2)))*as.numeric(sqrt(crossprod(r2)))
      res.BiOPD.dim2[[i]] <- a
      #
      # rm(q2, r2, a)
      eval(call("rm", q2))
      eval(call("rm", r2))
      eval(call("rm", a))
    }
  }
  #
  list(
    i0.dim1 = i0.dim1,
    tm.Criterion.dim1 = tm.Criterion.dim1,
    u = res.pmd.dim1[[i0.dim1]]$u,
    v = res.pmd.dim1[[i0.dim1]]$v,
    d = res.pmd.dim1[[i0.dim1]]$d,
    tm.dataMat.dim2 = tm.dataMat.dim2,
    res.pmd.dim2 = res.pmd.dim2,
    res.BiOPD.dim2 = res.BiOPD.dim2
  )
}
