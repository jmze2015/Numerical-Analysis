GNM <- function(var_list, f_list, x_init, max_step = 50, tol = 1e-6) {
  n <- length(var_list)
  if (length(f_list) != n) stop("Number of variables must equal number of equations.")
  
  grad_list <- lapply(f_list, Deriv, x = var_list)
  
  x_old <- as.numeric(x_init)
  results <- matrix(x_old, nrow = 1)
  
  for (step in seq_len(max_step)) {
    env <- setNames(as.list(x_old), var_list)
    
    ## Jacobian (robust build)
    J <- matrix(0, n, n, dimnames = list(NULL, var_list))
    for (i in seq_len(n)) {
      raw <- eval(grad_list[[i]], env)
      gi <- as.numeric(raw)
      nm <- names(raw)
      
      if (!is.null(nm) && length(nm) == length(gi)) {
        J[i, nm] <- gi
      } else if (length(gi) == n) {
        J[i, ] <- gi
      } else {
        cat("\nstep =", step, "x =", paste(x_old, collapse=","), "\n")
        print(raw)
        stop("Gradient returned with unexpected length/names.")
      }
    }
    
    ## F(x)
    Fval <- numeric(n)
    for (i in seq_len(n)) Fval[i] <- as.numeric(eval(f_list[[i]], env))
    
    if (any(!is.finite(J)) || any(!is.finite(Fval))) {
      cat("\nstep =", step, "x =", paste(x_old, collapse=","), "\n")
      print(J); print(Fval)
      stop("Non-finite J or F.")
    }
    
    ## Solve for delta with rank check + fallback
    qrJ <- qr(J, tol = 1e-12)
    if (qrJ$rank < n) {
      s <- svd(J)
      tol_s <- 1e-12 * max(s$d)
      d_inv <- ifelse(s$d > tol_s, 1 / s$d, 0)
      J_pinv <- s$v %*% (d_inv * t(s$u))
      delta <- as.vector(J_pinv %*% Fval)
    } else {
      delta <- as.vector(qr.coef(qrJ, Fval))
    }
    
    x_new <- x_old - delta
    results <- rbind(results, x_new)
    
    if (sqrt(sum(delta^2)) < tol || sqrt(sum(Fval^2)) < tol) break
    
    x_old <- as.numeric(x_new)
  }
  
  colnames(results) <- var_list
  results
}




# v <- c("x", "y")
# f <- list(
#   expression(x^2 + y^2 - 25),
#   expression(x^2 - y - 2)
# )
# init <- c(1, 4)

v <- c("x", "y", "z")
f <- list(
  expression(x^2 + y + z - 3),
  expression(x^2 + y^2 + z^2 - 5),
  expression(x - y)
)
init <- c(0, 0, 2)

GNM(v, f, init)





