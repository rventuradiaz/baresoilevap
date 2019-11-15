#' Forward difference approximation of first derivative
#'
#' @param i : position on discretised interval
#' @param x : mesh
#' @param f : function
#'
#' @return
#' @export
#'
#' @examples
fda_dfdx <- function(i=1, x, f=function(y){y}){
  dx <- x[i+1]-x[i]
  result <- f(x[i+1])-f(x[i])
  result <- result/dx
  return(result)
}

#' Backward difference approximation of first derivative
#'
#' @param i : position on discretised interval
#' @param x : mesh
#' @param f : function
#'
#' @return
#' @export
#'
#' @examples
bda_dfdx <- function(i=2, x, f=function(y){y}){
  dx <- x[i]-x[i-1]
  result <- f(x[i])-f(x[i-1])
  result <- result/dx
  return(result)
}

#' Central difference approximation of first derivative
#'
#' @param i : position on discretised interval
#' @param x : mesh
#' @param f : function
#'
#' @return
#' @export
#'
#' @examples
cda_dfdx <- function(i=2, x, f=function(y){y}){
  dx <- x[i]-x[i-1]
  result <- f(x[i+1])-f(x[i-1])
  result <- result/(2.0*dx)
  return(result)
}

#' Forward difference approximation of second order derivative
#'
#' @param i : position on discretised interval
#' @param x : mesh
#' @param f : function
#'
#' @return
#' @export
#'
#' @examples
fda_d2fdx2 <- function(i=1, x, f=function(y){y}){
  dx <- x[i+1]-x[i]
  result <- f(x[i+2])-2.0*f(x[i+1])+f(x[i])
  result <- result/(dx*dx)
  return(result)
}

#' Backward difference approximation of second order derivative
#'
#' @param i : position on discretised interval
#' @param x : mesh
#' @param f : function
#'
#' @return
#' @export
#'
#' @examples
bda_d2fdx2 <- function(i=3, x, f=function(y){y}){
  dx <- x[i]-x[i-1]
  result <- f(x[i])-2.0*f(x[i-1])+f(x[i-2])
  result <- result/(dx*dx)
  return(result)
}

#' Central difference approximation of second order derivative
#'
#' @param i : position on discretised interval
#' @param x : mesh
#' @param f : function
#'
#' @return
#' @export
#'
#' @examples
cda_d2fdx2 <- function(i=2, x, f=function(y){y}){
  dx <- x[i]-x[i-1]
  result <- f(x[i+1])-2.0*f(x[i])+f(x[i-1])
  result <- result/(dx*dx)
  return(result)
}
