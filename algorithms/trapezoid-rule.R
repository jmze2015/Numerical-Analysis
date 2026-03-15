## Trapezoid Rule Algorithm ##

options(digits = 16)

## inputs: function f, interval [a,b], and number of subintervals N

trap_int <- function(f, a, b, N){
  if (b < a) {
    stop("Improper interval was input")
  }
  
  ## main interval
  A <- round(a, 10)
  B <- round(b, 10)
  
  ## all subintervals
  all_ints <- seq(A, B, length.out = N+1)

  ## trapezoid rule
  sum <- f(all_ints[1]) + f(all_ints[N+1])
  
  for (i in 2:(N)){
    sum <- sum + 2*f(all_ints[i])
  }
  
  sum <- sum*(b-a)/(2*N)
  ## output
  return(sum)
}


## Results for Class

for (i in seq(10,50,10)){
  print(trap_int(cos, 0, pi/2, i))
}


## Adding Simpson's Rule as well / pretty easy addition

simpson_int <- function(f, a, b, N){
  if (b < a) {
    stop("Improper interval was input")
  }
  
  ## main interval
  A <- round(a, 10)
  B <- round(b, 10)
  
  ## all subintervals
  all_ints <- seq(A, B, length.out = N+1)
  
  ## Simpson's rule
  sum <- f(all_ints[1]) + f(all_ints[N+1])
  
  for (i in 2:(N)){
    sum <- sum + 2*f(all_ints[i])
  }
  
  for (i in 1:N){
    sum <- sum + 4*f(0.5*(all_ints[i]+all_ints[i+1]))
  }
  
  sum <- ((b-a)/(6*N))*sum
  ## output
  return(sum)
}

## Adding a Comparison of the Methods

Y1 <- c()
Y2 <- c()

for (i in (3:40)){
  Y1 <- c(Y1, trap_int(cos, 0, pi/2, i))
  Y2 <- c(Y2, simpson_int(cos, 0, pi/2, i))
}

plot(c(x,x),c(Y1,Y2), col = c(rep("red",38),rep("blue", 38)), type = "b",
     pch = 19, xlab = 'Number of subintervals', ylab = "Area under the curve",
     main = "Method Comparison")

legend("bottomright", legend = c("Trapezoid Rule", "Simpson's Rule"), col = c("red", "blue"), pch = 19)





