## Jacobi Method Algorithm ##


jacobi_algorithm <- function(A, b, x_init, tol = 1e-6, max_steps = 100){
  ## A must be square otherwise D^{-1} is impossible.
  if (ncol(A) == nrow(A)){
    n <- ncol(A)
    ## diagonal matrix
    D <- diag(diag(A))
    ## lower triangular matrix
    L <- A
    L[upper.tri(L)] <- 0
    diag(L) <- 0
    ## upper triangular matrix
    U <- A
    U[lower.tri(U)] <- 0
    diag(U) <- 0
    
    ## iteration
    invD <- solve(D)
    x0 <- x_init
    xseq <- matrix(x0, nrow = 2)
    for (i in 1: max_steps){
      x_new <- -1 * (invD %*% (L + U) %*% x0) + invD %*% b
      xseq <- cbind(xseq, x_new)
      
      if (norm(x_new - x0, type = "2") < tol){ break}
      x0 <- x_new
    }
    
    index <- ncol(xseq)
    
    return(list(
      seq = xseq, x = xseq[,index]
    ))
  }
  
  else{
    print("Error: Input matrix is not square")
  }
}

## Inputs from problem 1 
A <- matrix(c(2,1,1,3), ncol = 2)
b <- matrix(c(1,2), nrow = 2) 
x0 <- matrix(c(1,1), nrow = 2)

## Calling Jacobi Algorithm
jacobi_algorithm(A, b, x0, tol = 1e-10)











