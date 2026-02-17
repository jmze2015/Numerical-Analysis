install.packages("matlib")
library(matlib)

## swap row i with row j
elem1 <- function(A, i, j){
  e <- diag(nrow(A))
  c <- rowswap(e, from = i, to = j)
  return(c)
}

## scale row i with scalar s
elem2 <- function(A, i, s){
  e <- diag(nrow(A))
  e[i,] <- s*e[i,]
  return(e)
}

## to row i add s*row j
elem3 <- function(A, i, s, j){
  e <- diag(nrow(A))
  e[i,] <- e[i,] + s * e[j,]
  return(e)
}

## REF algorithm
ref <- function(A){
  e_list <- list()
  inv_list <- list()
  N <- nrow(A)
  M <- ncol(A)
  B <- A
  
  for (j in 1:min(N,M)) {
    
    ## find pivot in column j
    idx <- NA
    for (i in j:N) {
      if (B[i, j] != 0) {
        idx <- i
        break
      }
    }
    
    ## if entire column below is zero, move on
    if (is.na(idx)) next
    
    ## swap pivot row into position j
    E1 <- elem1(B, j, idx)
    B <- E1 %*% B
    e_list[[length(e_list) + 1]] <- E1
    inv_list[[length(inv_list) + 1]] <- E1
    
    ## eliminate entries below pivot
    if (j < N) {
      for (i in (j + 1):N) {
        
        if (B[i, j] == 0) next
        
        s <- - B[i, j] / B[j, j]
        E3 <- elem3(B, i, s, j)
        inv_E3 <- elem3(B, i, -s, j)
        
        B <- E3 %*% B
        e_list[[length(e_list) + 1]] <- E3
        inv_list[[length(inv_list) + 1]] <- inv_E3
      }
    }
  }
  
  L <- Reduce(`%*%`, inv_list)
  
  return(list(
    b = B,
    actions = e_list,
    inv_actions = inv_list,
    l = L
  ))
}


LU_Decomp <- function(A){
  L <- ref(A)$l
  U <- ref(A)$b
  return(
    list(
      l = L, u = U
    )
  )
}

LDU_Decomp <- function(A){
  L <- ref(A)$l
  U_til <- ref(A)$b
  
  d <- c()
  d_inv <- c()
  for(i in 1:ncol(U_til)){
    d <- c(d, U_til[i,i])
    d_inv <- c(d_inv, 1/(U_til[i,i]))
  }
  D <- diag(d)
  D_inv <- diag(d_inv)
  U <- D_inv %*% U_til
  
  return(list(
    l = L, d = D, u = U, a = L %*% D %*% U
  ))
}

Cholesky_Factorization <- function(A){
  if (identical (A, t(A))){
    if(det(A) != 0){
      L <- LDU_Decomp(A)$l
      D <- LDU_Decomp(A)$d
      U <- LDU_Decomp(A)$u
      
      D_bar <- sqrt(D)
      Chol <- L %*% D_bar ## cholesky factorization
      return(Chol)
    }
    else("Matrix input is singular")
  }
  else(print("Matrix input is not symmetric"))
}

C = matrix(c(1,2,3,4,2,5,7,9,3,7,11,14,4,9,14,19), ncol = 4,
           byrow = TRUE)

Cholesky_Factorization(C) 
t(Cholesky_Factorization(C))















# LU_Decomp(A)
# LDU_Decomp(A)
# A <- matrix(c(1,2,3,4,5,6,7,13,9), ncol = 3)
# S <- matrix(c(1,2,3,2,5,7,3,7,13), ncol=3)
# U <- ref(A)$b
# L <- ref(A)$l
# Cholesky_Factorization(A)
# Cholesky_Factorization(S) %*% t(Cholesky_Factorization(S))



# elem1(A,1,2) # swap rows 1 and 2
# elem2(A,3,10) # scale row 3 by 10
# elem3(A,2,-3,1) # add to row 2 (-3)*row 1




