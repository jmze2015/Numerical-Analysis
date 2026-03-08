## Lagrangian Interpolation ##

## install.packages("caracas")
library(caracas)

## Two lists for inputs: (x_{i}, f(x_{i}))_{i=0}^{n}

## n doesn't really seem to be a parameter??
Lnk <- function(k, nodes){
  ## correction to preserve theoretical indexing
  K <- k  
  ## save the k-th list element
  xk <- nodes[K]
  ## omit k-th element and save list
  xlist <- nodes[-K]
  
  
  expr_list <- list()
  
  for (i in 1:length(xlist)){
    
    g <- substitute((x - a) / (b - a), list(a = xlist[i], b = xk))
    expr_list[[i]] <- g
  }
  
  Reduce(function(a, b) call("*", a, b), expr_list)
}


lagrangian_interpolation <- function(nodex, nodef){
  if (length(nodex) != length(nodef)){
    stop("Error: Mismatch in node list lengths")
  }
  
  N <- length(nodex)
  sym_list <- vector("list", N)
    
  for (k in seq_len(N)){
    
    sym_exp <- as_sym(deparse(Lnk(k, nodex)))
    sym_list[[k]] <- nodef[k] * sym_exp
  }
  
  sum_sym <- Reduce(`+`, sym_list)
  return(expand(sum_sym))
}

x_list <- c(0,0.5,1,3)
f_list <- c(0,1,0,9)

lagrangian_interpolation(x_list, f_list)








