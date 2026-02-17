## Higher Dimensional Newton's Method ##

NM <- function(var_list, f_list, x_init, max_step, tol = 1e-6) {
  ## list of gradients (Jacobian rows)
  grad_list <- lapply(f_list, Deriv, x = var_list)
  
  ## initialize an empty matrix
  M <- matrix(NA_real_, nrow = length(f_list), ncol = length(var_list))
  
  ## starting point
  x_old <- x_init
  step <- 1
  
  ## store iterates
  results <- c(x_old)
  
  ## main loop
  while (step <= max_step){
    
    ## Inverse Jacobian at a point
    for (i in 1:length(var_list)){
      M[i,] <- eval(grad_list[[i]], list(x=x_old[1], y=x_old[2]))
    }
    inv_M <- solve(M)
    
    ## F eval at x_n
    f_x <- eval(f_list[1], list(x=x_old[1], y=x_old[2]))
    f_y <- eval(f_list[2], list(x=x_old[1], y=x_old[2]))
    f <- as.matrix(c(f_x, f_y), nrow = 2)
    
    ## NM update
    x_new <- x_old - inv_M %*% f
    
    ##record results
    results <- rbind(results, t(x_new))
  
    diff <- norm(x_new - x_old, type="2")
    if (diff < tol){
      break
    }
    
    ##loop back process
    x_old <- x_new
    step <- step + 1
  }
  
  colnames(results) <- var_list
  return(results)
}

## list of independent variables
v <- c("x","y")
## list of f_{i} for F: |R^{n} \to |R^{n}
f <- c(expression(x^2 + y^2 -25), expression(x^2 -y -2))
## x_{0}
init <- c(1,4)


## output results
print(NM(v,f,init,10,1e-8), digits =12)


