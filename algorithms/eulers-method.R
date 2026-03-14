## Euler's Method 


## Current Euler Method has some limitations
## Need to add:
## 1) general function input
## 2) more robust update loop

euler_method <- function(h, t0, y0){
  
  y_old <- y0
  t_old <- t0
  
  l <- seq(0,5, by = h)
  
  y_seq <- c(y0)
  
  ## loop here
  for (i in l){
    y_new <- y_old + h * (-15*y_old)
    y_seq <- c(y_seq, y_new)
    y_old <- y_new
  }
  
  l <- c(l, l[length(l)]+h)

  plot(l,y_seq, type = "b", xlab="t", ylab = "y(t)",
       main =paste("Euler Method Approximation h=",h), pch = 19)
  curve(-exp(-15*x), add = TRUE, col="red")
  legend("bottomright", legend = c("Approximate Solution", "Analytic Solution"),
         col = c("black", "red"), pch = c(19,10))
  
  
}

euler_method(0.25,0,-1)
euler_method(2/15,0,-1)
euler_method(1/100,0,-1)




