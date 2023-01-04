#' Title Decomposition for Dim2
#'
#' @param Fd Frequency data
#' @param Dn diag matrix of RowSum
#' @param Dp diag matrix of ColumnSum
#' @param results of CA
#' @param cK cumulative Nonzeros. For Criterion=='Nonzeros'
#' @param Sparse.type tau or grid
#' @param Criterion.type default=IS
#' @param decomposition.type 'BiOPD'
#' @param temp.Decomposition.dim1
#'
#' @return Results by Decomposition for Dim2
#' @export
#'
#' @examples Decomposition.Criterion.Dim2
Decomposition.Criterion.Dim2 <- function(Fd,
                                         Yn1, v0,
                                         Dn, Dp,
                                         re,
                                         cK,
                                         Sparse.type,
                                         Criterion.type,
                                         decomposition.type,
                                         temp.Decomposition.dim1
){
  n <- dim(Fd)[1]; p <- dim(Fd)[2]
  # r <- as.matrix( diag(unlist(Dn) )); c<-as.matrix(diag(Dp))
  l = length(temp.Decomposition.dim1$res.BiOPD.dim2);
  #
  i0.dim1 <- temp.Decomposition.dim1$i0.dim1;
  Nonzero.pmduv <- matrix(NA, l, 2);
  Nonzero.v.union <- rep(NA, l); #union NonzeroIndex for the First 2 dim.
  #
  k <- cK[2];
  #
  if(decomposition.type == 'PMD'){
    res.pmd.dim2 <- temp.Decomposition.dim1$res.pmd.dim2
  }
  else if(decomposition.type == 'mPMD'){
    res.pmd.dim2 <- temp.Decomposition.dim1$res.pmd.dim2
  }
  else if (decomposition.type == 'BiOPD'){
    res.pmd.dim2 <- temp.Decomposition.dim1$res.BiOPD.dim2 ####
  }
  if(Sparse.type == 'tau'){ #"BothSide"
    # res.pmd.dim1 = res.pmd.BothSide.dim1
    #Criterion
    if(Criterion.type == 'IS'){
      tm.IS.dim2 <- ISpmd(res.pmd =
                             res.pmd.dim2,
                           Yn1,
                           lambda0 = re$eig[2,1],
                           l=l, r=1)
      i0.dim2 <- tm.IS.dim2$i0.IS;
      tm.Criterion.dim2 = tm.IS.dim2
    }
  }
  else if(Sparse.type == 'grid'){ #"Joint"
    #Criterion
    if(Criterion.type == 'IS'){
      tm.IS.dim2 <- ISpmd(res.pmd = res.pmd.dim2,
                           lambda0 = re$eig[2,1],
                           l=l, r=1)
      i0.dim2 <- tm.IS.dim2$i0.IS;
      tm.Criterion.dim2 = tm.IS.dim2
    }
    else if(Criterion.type == 'Nonzeros'){
      # k <- 50;
      tm.Criterion.dim2 <- NULL;
      for(i in 1:l){
        Nonzero.pmduv[i, 1] <-  sum(abs(res.pmd.dim2[[i]]$u[, 1]) != 0)#length(which(abs(res.pmd.dim2[[i]]$u[, 1]) > 1e-6)); #for U
        Nonzero.pmduv[i, 2] <-  sum(abs(res.pmd.dim2[[i]]$v[, 1]) != 0)#length(which(abs(res.pmd.dim2[[i]]$v[, 1]) > 0)); #for V
        Nonzero.v.union[i] <-  length(union(
          which( abs(as.numeric(temp.Decomposition.dim1$res.pmd.dim1[[i0.dim1]]$v) ) > 0),
          which( abs( res.pmd.dim2[[i]]$v) > 0)
          ))
      }
      i0.dim2 <-
        which.min(abs(Nonzero.pmduv[,2]- k));
      i0.dim2.union <-
        which.min(abs(Nonzero.v.union- k));
      tm.Criterion.dim2$Nonzero.pmduv <- Nonzero.pmduv;
      tm.Criterion.dim2$Nonzero.v.union <- Nonzero.v.union
    }
  }
  #
  list(
    i0.dim1 = temp.Decomposition.dim1$i0.dim1,
    i0.dim2 = i0.dim2,
    # i0.dim2.union = i0.dim2.union,
    tm.Criterion.dim1 = temp.Decomposition.dim1$tm.Criterion.dim1,
    tm.Criterion.dim2 = tm.Criterion.dim2
  )
}
