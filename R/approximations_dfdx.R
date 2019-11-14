fda_dfdx <- function(i=1, x, f=function(y){y}){
  dx <- x[i+1]-x[i]
  result <- f(x[i+1])-f(x[i])
  result <- result/dx
  return(result)
}

bda_dfdx <- function(i=2, x, f=function(y){y}){
  dx <- x[i]-x[i-1]
  result <- f(x[i])-f(x[i-1])
  result <- result/dx
  return(result)
}

cda_dfdx <- function(i=2, x, f=function(y){y}){
  dx <- x[i]-x[i-1]
  result <- f(x[i+1])-f(x[i-1])
  result <- result/(2.0*dx)
  return(result)
}
