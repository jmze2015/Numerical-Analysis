## Power Method Algorithm ##

PM <- function(A, x_init, tol = 1e-6, max_steps = 500){
  x_old <- x_init
  
  i <- 1
  while(i <= max_steps){
    x_new <- (A %*% x_old) / norm(A %*% x_old, type="2")
    i <- i + 1
    if(norm(x_new - x_old, type="2") < tol){
      break
    }
    x_old <- x_new
  }
  return(list(
    vec = x_new, val = (A%*%x_new)[1]/(x_new[1])
  ))
}

## Given SPD Matrix
A <- matrix(c(2,1,1,1,2,1,1,1,2), byrow = TRUE, ncol = 3)
## Given initial vector
x0 <- matrix(c(1,-1,2), nrow = 3)
## Calling Power Method Algorithm
PM(A,x0)
