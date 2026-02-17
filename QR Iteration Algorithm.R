## QR Iteration Algorithm

## Turns a square matrix into an upper Hessenberg matrix
Hessenberg <- function(A){
  ## check for square matrix
  if(nrow(A) == ncol(A)){
    N <- nrow(A)
    H <- A
    for(i in 1:(N-2)){
      x <- as.matrix(H[(i+1):N, i])
      e1 <- c(1, rep(0, length(x)-1))
      
      q <- x + sign(x[1])*norm(x,"2") * e1
      q <- q / (norm(q, type="2"))
      Psmall <- diag(length(q)) - 2*q %*% t(q)
      
      P <- diag(N)
      P[(i+1):N, (i+1):N] <- Psmall
      
      
      H <- P %*% H %*% P
    }
    rtol <- 1e-12
    H[abs(H) < rtol] <- 0
    return(H)
  }
  
  else{
    print("Matrix input is not square")
  }
}

## Runs one QR step on an arbitrary matrix
QR_step <- function(H){
  N <- ncol(H)
  
  for (i in 1:(N-1)){
    
    v <- H[c(i,i+1),i]
    r <- norm(v, "2")
    if (r == 0) next
    
    c <- v[1] / r
    s <- v[2] / r
    
    BR <- diag(N)
    BR[c(i, i+1), c(i, i+1)] <- matrix(c(c, s, -s, c), nrow = 2)
  
    H <- t(BR) %*% H %*% BR
  }
  return(H)
}

## Takes arb. matrix -> Hessenberg -> Runs many QR steps
QR_algorithm <- function(A, maxiter = 500, rtol = 1e-12){
  H <- Hessenberg(A)
  for (k in 1:maxiter){
    H <- QR_step(H)
  }
  tol <- rtol
  H[abs(H) < tol] <- 0
  return(H)
}

## test matrix from HW
test_matrix <- matrix(c(1,2,3,4,2,4,9,7,3,9,1,1,4,7,1,10), ncol = 4)

# part a : H^{(0)}
H <- Hessenberg(test_matrix)
H
# part b : H^{(1)}
H1 <- round(QR_step(H), 10)
H1
# part c : H^{*}
round(QR_algorithm(test_matrix, 2000),10)

# part d : Verification of eigenvalues
eigen(test_matrix)$val








# test_matrix_2 <- matrix(c(2,1,1,1,2,1,1,1,2), ncol = 3)
# QR_algorithm(test_matrix_2)




