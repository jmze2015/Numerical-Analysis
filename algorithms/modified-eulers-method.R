## Modified Euler's Method

## give expression in terms of y and t
mod_euler_method <- function(y, h, t0, y0, a, b){
  
  y_expr <- substitute(y)
  
  if (is.symbol(y_expr)) {
    y_expr <- eval.parent(y_expr)
  }
  
  ## num evalutaion function
  EVAL <- function(f, t, y){
    eval(f, list(t = t, y = y))
  }
  
  t_f <- seq(t0, b, by = h)
  y_f <- numeric(length(t_f))
  y_f[1] <- y0
  
  for(i in seq_len(length(t_f) - 1)){
    y_f[i + 1] <- y_f[i] + (h/2) * EVAL(y_expr, t_f[i], y_f[i]) + (h/2) * EVAL(y_expr, t_f[i+1], y_f[i] + h*EVAL(y_expr, t_f[i], y_f[i])) 
  }
  
  t_b <- seq(t0, a, by = -h)
  y_b <- numeric(length(t_b))
  y_b[1] <- y0
  
  for(i in seq_len(length(t_b) - 1)){
    y_b[i + 1] <- y_b[i] - h * EVAL(y_expr, t_b[i], y_b[i])
  }
  
  t_all <- c(rev(t_b), t_f[-1])
  y_all <- c(rev(y_b), y_f[-1])
  
  plot(t_all, y_all,
       type="b",
       pch=19,
       xlab="t",
       ylab="y(t)",
       main=paste("Modified Euler Method Approximation of, y'=",as.expression(y_expr),
                  " with mesh size h =",h))
}


mod_euler_method((y+t)^2, 0.01, 0, -1, 0,1)
curve(tan(x-pi/4)-x, add = TRUE, col="red")
legend("bottomright", legend = c("Approximate Solution", "Analytic Solution"),
       col = c("black", "red"), pch = c(19,10))

