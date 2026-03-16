## RK4 Code##

## inputs: h, y, yo, a, b
RK4 <- function(h, y, y0, a, b){
  
  y_expr <- substitute(y)
  
  if (is.symbol(y_expr)) {
    y_expr <- eval.parent(y_expr)
  }
  
  ## num eval functions
  y_prime <- function(f, t, y){
    eval(f, list(t = t, y = y))
  }
  
  ## sequences for t and y
  t_seq <- seq(a, b, by = h)
  n <- length(t_seq)
  y_seq <- numeric(n)
  y_seq[1] <- y0
  
  ## main update
  for(i in 1:(n-1)){
    ## using RK4 definitions
    F1 <- h*y_prime(y_expr, t_seq[i], y_seq[i])
    F2 <- h*y_prime(y_expr, t_seq[i]+(h/2), y_seq[i]+(1/2)*F1)
    F3 <- h*y_prime(y_expr, t_seq[i]+(h/2), y_seq[i]+(1/2)*F2)
    F4 <- h*y_prime(y_expr, t_seq[i]+h, y_seq[i]+F3)
    ## main explicit update
    y_seq[i + 1] <- y_seq[i] + (1/6)*(F1 + 2* F2 + 2 *F3 + F4)
  }

  ## plotting approximation
  plot(t_seq, y_seq,
       type="b",
       pch=19,
       xlab="t",
       ylab="y(t)",
       main=paste("RK4 Approximation of, y'=",as.expression(y_expr),
                  " with mesh size h =",h))
}

RK4(0.01, (y+t)^2, -1, 0, 1)
curve(tan(x-pi/4)-x, add = TRUE, col="red")
legend("bottomright", legend = c("Approximate Solution", "Analytic Solution"),
       col = c("black", "red"), pch = c(19,10))

