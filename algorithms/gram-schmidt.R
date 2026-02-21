## Goal of Gram-Schmidt Process is to output an orthonormal basis given some 
## arbitrary basis.

## input is a basis in the form of columns of a matrix A
gSchmidt <- function(A){
  if (qr(A)$rank == ncol(A)){ ##not an accurate condition!!
    
    n <- ncol(A)
    
    ## init orthonormal basis
    W <- as.matrix(A[,1] / norm(A[,1], "2"))
    
    ## induction
    for (i in 2:n){
      U <- W
      v <- A[,i]
  
      projs <- rep(0,nrow(A))
      
      for (j in 1:ncol(U)){
        u <- U[,j]
        projs <- projs + (as.numeric(crossprod(v,u))) * u
      }
      
      w <- v - projs
      w <- w/(norm(w, type = "2"))
      W <- cbind(W,w)
    }
    
    W <- as.data.frame(W)
    colnames(W) <- paste0("u", seq_len(ncol(W)))
    
    return(W)
  }
  else{
    print("Error: list provided is not a basis")
  }
}

## Test matrix #1
A <- matrix(c(1,0,1,-1,2,0,0,0,3), ncol = 3)
out <-gSchmidt(A)
out


## Test matrix #2
B <- matrix(c(3,0,1,-2,1,0), ncol = 2)
out <- gSchmidt(B)
out

## Verification of Orthonormal Basis
Out <- matrix(0, 2, 2)
for (i in 1:2){
  for (j in 1:2){
    Out[i,j] <- crossprod(out[,i], out[,j])
  }
}
## Rounding to see kronecker-delta
round(Out, 12)



