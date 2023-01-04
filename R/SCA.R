#' Title Return the coordinates for SCA
#'
#' @param X0 data matrix (Contingency table)
#' @param para.pmd parameters for 'tau'
#' @param para.uv parameters for 'grid'
#' @param Sparse.type 'tau' or 'grid'
#' @param Criterion.type default=='IS'
#' @param decomposition.type 'BiOPD'
#' @param index.Col.Nonzero  == NULL or 1 #
#' @param cK cumulative Nonzeros. For Criterion=='Nonzeros'
#' @param R dimensions for output
#'
#' @return the coordinates(a,b,f,g) and (eig,PEV) for SCA
#' @export
#'
#' @examples SCA()
#' #' load 'Data/ColorSound.rda'
SCA <- # Dim by Dim
  function(X0,
           para.pmd,
           para.uv,
           Sparse.type,
           decomposition.type,
           Criterion.type,
           index.Col.Nonzero, #
           cK,
           R){
    #
    # -- Obtain matrix Yn1: sqrt(solve(Dn))%*%(X-r%*%t(c))%*%sqrt(solve(Dp)) # with the 'centered' X
    #
    re.CA <- CA(X0);
    n <- dim(X0)[1]; p <- dim(X0)[2] # dim of data matrix
    Fd <- prop.table(X0) # X0/sum(X0); the frequency data matrix;
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
    # -- Return the weighted Variance
    weightVar <- function(vec, W){
      return(t(vec)%*%W%*%vec)
    }
    tempYn <- ObtainYn(Fd)
    Dn <- tempYn$Dn; Dp <- tempYn$Dp
    # r <- ObtainYn(Fd)$r; c <- ObtainYn(Fd)$c
    Yn1 <- tempYn$Yn1;
    X <- Fd;
    Y=Yn1;
    res.pmd0 <- PMD(Yn1, type="standard",sumabs=1,  # using $Yn1$, already substract the eigenvector for $\lambda = 1$
                        sumabsu=NULL, sumabsv=NULL, lambda=NULL,
                        niter=20, K=1, v=NULL, trace=FALSE,
                        center=FALSE, chrom=NULL, rnames=rownames(X), cnames=colnames(X),
                        upos=FALSE, uneg=FALSE, vpos=FALSE, vneg=FALSE);
    v0 <- res.pmd0$v.init;
    #
    SCA.out <-
      SCApmd(
        X =Fd,
        para.pmd,
        para.uv,
        v0 = v0,
        cK=cK,
        Sparse.type,
        decomposition.type,
        Criterion.type)
    #
    tm.res.pmd.Dim <- list(SCA.out$res.pmd.dim1[[SCA.out$temp.Criterion.dim1$i0.dim1]],
                           SCA.out$res.BiOPD.dim2[[SCA.out$Criterion.BiOPD$i0.dim2]]);
    K <- min(n, p)-1;  m <- length(index.Col.Nonzero)
    r <- as.matrix(diag(Dn)); c<-as.matrix(diag(Dp))
    # Score
    # (u, d, v) by SVD by Bi-OPD
    tm.res.pmd.Dim.u <- NULL; tm.res.pmd.Dim.v <- NULL; tm.res.pmd.Dim.d <- NULL;
    for(j in 1:R){
      tm.res.pmd.Dim.u <- cbind(tm.res.pmd.Dim.u, tm.res.pmd.Dim[[j]]$u);
      tm.res.pmd.Dim.v <- cbind(tm.res.pmd.Dim.v, tm.res.pmd.Dim[[j]]$v)
      tm.res.pmd.Dim.d <- c(tm.res.pmd.Dim.d, tm.res.pmd.Dim[[j]]$d)
    }
    # check signs (using the obtained No.l u1, u2, compared with a1,a2 resulted by CA. )
    a0 <-
      # cbind(-re.CA$row$coord[, 1], re.CA$row$coord[, 2]); #for Speech
      cbind(re.CA$row$coord[, 1], re.CA$row$coord[, 2]);    #for Toy data
    b0 <-
      # cbind(-re.CA$col$coord[,1], re.CA$col$coord[,2])
      cbind(re.CA$col$coord[,1], re.CA$col$coord[,2])
    #
    # for (j in 1:R) {
    #   if (sum(Dn %*% tm.res.pmd.Dim.u[, j] * a0[, j]) < 0){
    #     tm.res.pmd.Dim.u[, j] <- -1 * tm.res.pmd.Dim.u[, j]
    #   }
    #   if(sum(Dp %*% tm.res.pmd.Dim.v[, j] * b0[, j]) < 0){
    #     tm.res.pmd.Dim.v[, j] <- -1 * tm.res.pmd.Dim.v[, j]
    #   }
    # }
    #Output
    # Coordinates by 'row & column' Score
    rScore.pmd <-
      solve(Dn) %*% as.matrix(Fd-r%*%t(c)) %*% sqrt(solve(Dp))  %*% tm.res.pmd.Dim.v;      #
    cScore.pmd <-
      solve(Dp) %*% t(Fd-r%*%t(c)) %*% sqrt(solve(Dn)) %*% tm.res.pmd.Dim.u;      #
    rownames(rScore.pmd) <- rownames(Fd); rownames(cScore.pmd) <- colnames(Fd)
    #
    #
    #index.Weights; index.Score; (based on the obtained weights & score)
    {
      # index.Col <- paste0('index.Col.Weights.dim', 1:R);
      # index.Row <- paste0('index.Row.Weights.dim', 1:R);
      # for(j in 1:R){
      #   assign(index.Col[j], which(tm.res.pmd.Dim[[j]]$v != 0));
      #   assign(index.Row[j], which(tm.res.pmd.Dim[[j]]$u != 0))
      # }
      #
      index.Col.Weights.dim1 <- which(tm.res.pmd.Dim[[1]]$v != 0);
      index.Col.Weights.dim2 <- which(tm.res.pmd.Dim[[2]]$v != 0);
      index.Row.Weights.dim1 <- which(tm.res.pmd.Dim[[1]]$u != 0);
      index.Row.Weights.dim2 <- which(tm.res.pmd.Dim[[2]]$u != 0);
      index.col.Weights <- union(index.Col.Weights.dim1, index.Col.Weights.dim2)
      index.row.Weights <- union(index.Row.Weights.dim1, index.Row.Weights.dim2)
      #index.Score
      # index.Col.Score <- paste0('index.Col.Score.dim', 1:R);
      # index.Row.Score <- paste0('index.Row.Score.dim', 1:R);
      # for(j in 1:R){
      #   assign(index.Row.Score[j], which(rScore.pmd[,j] != 0));
      #   assign(index.Col.Score[j], which(cScore.pmd[,j] != 0))
      # }
      index.Col.Score.dim1 <- which(cScore.pmd[,1] != 0)
      index.Col.Score.dim2 <- which(cScore.pmd[,2] != 0)
      index.Row.Score.dim1 <- which(rScore.pmd[,1] != 0)
      index.Row.Score.dim2 <- which(rScore.pmd[,2] != 0)
      index.col.Score <- union(index.Col.Score.dim1, index.Col.Score.dim2)
      index.row.Score <- union(index.Row.Score.dim1, index.Row.Score.dim2)
      #
      index.col <- NULL; index.row <- NULL;
      index.col$Score <- index.col.Score; index.col$Weights <- index.col.Weights;
      index.row$Score <- index.row.Score; index.row$Weights <- index.row.Weights
    }
    #
    # Coordinates by SVD sparse weights
    U <- sqrt(solve(Dn)) %*% tm.res.pmd.Dim.u %*% diag(tm.res.pmd.Dim.d[1:R])#==Dn^{-1/2}%*%u*sqrt(lambda
    V <- sqrt(solve(Dp)) %*% tm.res.pmd.Dim.v %*% diag(tm.res.pmd.Dim.d[1:R])
    rownames(U) <- rownames(Fd);   rownames(V) <- colnames(Fd)
    #
    # pseudo Eigenvalue
    eig1 <- NULL;  eig2 <- NULL;
    for(j in 1:R){
      eig1 <- c(eig1, tm.res.pmd.Dim[[j]]$d^2)    #based on SVD
      eig2 <- c(eig2, sum(Dn%*% rScore.pmd[,j]^2))#based on row contribution
      # diag(crossprod(Y%*%tm.res.pmd.Dim$v)) #before standardizd
    }
    # PEV1 <- rep(0, R);  PEV2 <- NULL;
    PEV1 = eig1/sum(re.CA$eig[,1]);
    PEV2 = eig2/sum(re.CA$eig[,1])
    #
    # rScore and cScore, (Original Score)
    F.. <- rScore.pmd;   G.. <- cScore.pmd
    #
    temp.F <- NULL; temp.G <- NULL;
    rFs <- matrix(rep(0,2*n),nrow=n,ncol=2, byrow = T);
    cGs <- matrix(rep(0,2*p),nrow=p,ncol=2, byrow = T);
    #
    if(length(index.Col.Nonzero)==0){ #for normal case
      F. <- rScore.pmd
      G. <- cScore.pmd
      #
      for(j in 1:R){
        temp.F <- c(temp.F, tm.res.pmd.Dim[[j]]$d / sqrt(weightVar(rScore.pmd, Dn)[j, j]))
        temp.G <- c(temp.G, tm.res.pmd.Dim[[j]]$d / sqrt(weightVar(cScore.pmd, Dp)[j, j]))
      }
      a <- rScore.pmd %*% diag(temp.F);
      b <- cScore.pmd %*% diag(temp.G)
      #
      for(j in 1:R){
        # assign(rFs[j], Dn%*% F.[,j]^2)  #r*F1^2, # (sum(rF1sqr),sum(rF2sqr)) = (lambda1, lambda2)
        # assign(cGs[j], diag(weightVar(G, Dp.update))[j])  # Dp%*%G1^2; cG2sqr <- Dp%*%G2^2; cG3sqr <- Dp%*%G3^2
        rFs[,j] <- Dn %*% a[,j]^2
        cGs[,j] <- Dp %*% b[,j]^2 #for normal case
      }
    }
    else{
      # update Dp, s.t. sum(Dp.update)==1; G (value in [ setdiff(1:p, index.Col.Nonzero) ] be zero)
      index.Col.Nonzero <-
        index.col$Weights #based on Nonzeros in weights; # index.col$Score
      # index.nonzero(cScore.pmd[,1:2]) #based on Nonzeros
      #
      # G <- D.Select(p, index.Col.Nonzero)%*% cScore.pmd  #based on index.Col.Nonzero
      # temp <- Dp.Update(Dp, index.Col.Nonzero); Dp.update <- temp$Dp.update  #sum(Dp.update)=1
      # Solve.Dp.update <- temp$solve.Dp.update  #Inverse for Dp.update
      #
      for(j in 1:R){
        temp.F <- c(temp.F, tm.res.pmd.Dim[[j]]$d / sqrt(weightVar(rScore.pmd, Dn)[j, j]))
        temp.G <- c(temp.G, tm.res.pmd.Dim[[j]]$d / sqrt(weightVar(cScore.pmd, Dp)[j, j]))
      }
      a <- rScore.pmd %*% diag(temp.F);
      b <- cScore.pmd %*% diag(temp.G)
      #
      for(j in 1:R){
        rFs[,j] <- #r*F1^2, # (sum(rF1sqr),sum(rF2sqr)) = (lambda1, lambda2)
          Dn %*% a[,j]^2
        cGs[,j] <- # Dp%*%G1^2; cG2sqr <- Dp%*%G2^2
          # Dp.update%*%G.[,j]^2
          Dp %*% b[,j]^2
      }
    }
    #
    F.Ctr <- NULL; G.Ctr <- NULL; U.Ctr <- NULL; V.Ctr <- NULL;
    for(j in 1:R){
      #   F.Ctr <- cbind(F.Ctr, rFs[,j]/sum(rFs[,j]));
      #   G.Ctr <- cbind(G.Ctr, cGs[,j]/sum(cGs[,j]));
      U.Ctr <-
        tm.res.pmd.Dim.u^2
      # cbind(U.Ctr, Dn%*%U[,j]^2/tm.res.pmd.Dim[[j]]$d^2);
      V.Ctr <-
        tm.res.pmd.Dim.v^2
      # cbind(V.Ctr, Dp%*%V[,j]^2/tm.res.pmd.Dim[[j]]$d^2)
    }
    F.Ctr <- rFs%*%diag(1/apply(rFs,2,sum));
    G.Ctr <- cGs%*%diag(1/apply(cGs,2,sum))
    # F.Ctr <- round(1000*F.Ctr, 0); G.Ctr <- round(1000*G.Ctr, 0)#eig0[1]
    #
    # ## row Score and Column Score; after [ Normlized ] !!!!!!!
    # FF <- c(paste0('F',1:R)); GG <- c(paste0('G',1:R));
    # for(j in 1:R){
    #   #variance==1
    #   assign(FF[j], sqrt(n)* sqrt(Dn)%*%F.[,j]/tm.res.pmd.Dim[[j]]$d);
    #   assign(GG[j], sqrt(m)* sqrt(Dp.update)%*%G.[,j]/tm.res.pmd.Dim[[j]]$d);
    #   #variance!=1, crossprod(F)==lambda
    #   # assign(FF[j],  sqrt(Dn)%*%F.[,j]);
    #   # assign(GG[j],  sqrt(Dp.update)%*%G.[,j]);
    #   #
    #   # check signs (using the obtained No.l u1, u2, compared with a1,a2 resulted by CA. )
    #   # if (sum(Dn %*% U[,j] * re.CA$row$coord[, 1]) < 0  &&
    #   #     sum(Dp %*% V[,j] * re.CA$col$coord[, 1]) < 0)
    #   # {
    #   #   U[,j] <- -1 * U[,j]
    #   #   V[,j] <- -1 * V[,j]
    #   # }
    # }
    # #
    # rScore.Normd <- cbind(F1, F2); cScore.Normd <- cbind(G1, G2)
    #
    # F1.cos2 <- 1000*round(F1^2/diag(F.%*%t(F.)),3); F2.cos2 <- 1000*round(F2^2/diag(F.%*%t(F.)),3)
    # G1.cos2 <- 1000*round(G1^2/diag(G.%*%t(G.)),3); G2.cos2 <- 1000*round(G2^2/diag(G.%*%t(G.)),3)
    #
    # rd2 <- round(diag(Dn%*%F.%*%t(F.)),3); cd2 <- round(diag(Dp%*%G.%*%t(G.)),3) #sum=0.746    #
    # # PEV for Row and Column
    # rPEV <- c(sum(rF1sqr), sum(rF2sqr))/sum(rd2); cPEV <- c(sum(cG1sqr), sum(cG2sqr))/sum(cd2)
    #
    F. <- a; G. <- b;
    ##Cos2row & Cos2col
    CosRow <-
      # diag(1/(diag(F.%*%t(F.))))%*%F.^2
      diag(1/sqrt(diag(crossprod(t(F.)))))%*%F.
    # (diag(1/diag(re$row$coord[,1:2]%*%t(re$row$coord[,1:2])))%*%re$row$coord[,1:2]^2)
    # prop <- res.pmd$prop.var.explained
    # #
    row.names(F.) <- rownames(Fd); row.names(G.) <- colnames(Fd)
    row.names(rScore.pmd) <- rownames(Fd);
    # row.names(F.Ctr) <- rownames(Fd); row.names(G.Ctr) <- colnames(Fd)
    # rownames(U.Ctr) <- rownames(Fd); rownames(V.Ctr) <- colnames(Fd)
    # row.names(CosRow) <- rownames(Fd)
    # prop <- res.pmd$prop.var.explained
    # #
    row.names(F.) <- rownames(Fd); row.names(G.) <- colnames(Fd)
    row.names(rScore.pmd) <- rownames(Fd);
    ##
    list(
      a = F., b = G., # weighted.Var = d^2
      #Weights
      f = U, g = V,
      eig1 = eig1, eig2 = eig2,
      PEV1 = PEV1, PEV2 = PEV2,
      rFsqr = rFs, cGsqr = cGs,
      F.Ctr = F.Ctr, G.Ctr = G.Ctr,
      U.Ctr = U.Ctr, V.Ctr = V.Ctr,
      # rd2 = rd2, cd2 =cd2,
      # prop = prop
      # rPEV = rPEV, cPEV = cPEV,
      index.col = index.col,
      index.row = index.row,
      # SVD
      u = tm.res.pmd.Dim.u,
      v = tm.res.pmd.Dim.v,
      d = tm.res.pmd.Dim.d
    )#
  }
