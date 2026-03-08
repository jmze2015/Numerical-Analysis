## Steepest Descent Algorithm ##
library(rlang)
g <- expression(cos(x+y)+sin(x)+cos(y))

## Gradient at a point function (inputs: expression, output: eval at (x,y))

Grad <- function(f, x, y){
  grad_f <- Deriv(f, c("x", "y"))
  x <- x
  y <- y
  return(eval(grad_f))
}

## 3D Steepest Descent Algorithm (2D expression input)

SDM <- function(f, x_init, y_init, max_steps=100, precision=100){
  
  x_old <- x_init
  y_old <- y_init
  
  
  results <- matrix(c(x_old, y_old), ncol = 2, byrow = TRUE)
  
  #main loop
  i <- 1
  while (i<= max_steps){
    
    grad <- Grad(f, x_old, y_old)
    
    ## stopping criteria
    if (sqrt(sum(grad^2)) < 1e-6) break
    
    #create random positive scaling
    a_list <- rweibull(precision, 1.5, 1)
    
    #eval along descent line
    x_test <- x_old - a_list * grad[1]
    y_test <- y_old - a_list * grad[2]
    
    
    vals <- eval(f, env(x = x_test, y = y_test))
    
    ## pick minimizing step size
    a <- a_list[which.min(vals)]
    
    ## update
    x_new <- x_old - a * grad[1]
    y_new <- y_old - a * grad[2]
    
    results <- rbind(results, c(x_new, y_new))
    
    x_old <- x_new
    y_old <- y_new
    
    i <- i + 1
  }
  
  colnames(results) <- c("x","y")
  results
}

## Run Steepest Descent x_0 = (0,1).
print(SDM(g, 0, 1), digits = 10)













