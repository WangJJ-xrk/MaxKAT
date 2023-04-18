# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#' calcKerMat.gau calculate similarity matrix via Gaussian kernel function
#'
#' @param X the data to calculate its similarity matrix, with each row representing
#'  a subject
#' @param r the scale parameter used in Gaussian kernel
#'
#' @return the centered similarity matrix
#' @export
#'
#' @examples library(MASS)
#' n = 100
#' k = 10
#' x = mvrnorm(n, rep(0,k), diag(k))
#' calcKerMat.gau(x)
calcKerMat.gau <- function(X, r=NULL){
  dis.X = dist(X)
  if(is.null(r) == TRUE){
    r = median(dis.X)^2
  }
  mat.X = as.matrix(dis.X)
  n = dim(mat.X)[1]
  K.X = exp(-mat.X^2/(r))
  H = diag(rep(1,n)) - matrix(1/n,n,n)
  cK.X = H %*% K.X %*% H
  return(cK.X)
}


#' calcKerMat.poly calculate similarity matrix via polynomial kernel function
#'
#' @param X the data to calculate its similarity matrix, with each row representing
#'  a subject
#' @param r the offset value uesed in kernel function
#' @param d the order of the kernel function
#'
#' @return the centered similarity matrix
#' @export
#'
#' @examples library(MASS)
#' n = 100
#' k = 10
#' x = mvrnorm(n, rep(0,k), diag(k))
#' calcKerMat.poly(x)
calcKerMat.poly <- function(X,r = 0,d = 5){
  n = dim(X)[1]
  H = diag(rep(1,n)) - matrix(1/n,n,n)
  K.X = (X %*% t(X) + r)^d
  cK.X = H %*% K.X %*% H
  return(cK.X)
}


#' modify.normalization data normalization
#'
#' @param x the data to be normalized
#' @return the normalized data
#' @export
#'
#' @examples library(MASS)
#' n = 100
#' x = rnorm(n, 0, 1)
#' modify.normalization(x)
modify.normalization <- function (x){
  x <-as.vector(x)
  dex <- which(!is.na(x))

  row.sum <- sum(x[dex])
  row.valid <- length(dex)

  row.mean <- row.sum/row.valid
  p <- (row.sum+0.5)/(1+row.valid)

  y <- rep(0,length(x))
  y[dex] <- (x[dex]-row.mean)/sqrt(p*(1-p))

  y
}



#' geneSimi calculate similarity matrix for genetic data
#'
#' @param x genetic data, with each column representing a subject
#'
#' @return a similarity matrix
#' @export
#'
#' @examples library(MASS)
#' n = 100
#' k = 10
#' x = mvrnorm(n, rep(0,k), diag(k))
#' geneSimi(t(x))
geneSimi <- function(x){
  n = dim(x)[1]

  k = dim(x)[2]

  x <- x/2
  if( k>1 )
  {
    z <- apply(x, 1, modify.normalization)
  }
  else
  {
    z <- modify.normalization(x)
    z <- matrix(z, ncol=1)
  }

  S0 = z %*% t(z)

  S = S0/k

  return(S)

}

#' calcKer.IBS calculate similarity matrix for genetic data
#'
#' @param X genetic data, with each row representing a subject
#'
#' @return a similarity matrix
#' @export
#'
#' @examples library(MASS)
#' n = 100
#' k = 10
#' X = mvrnorm(n, rep(0,k), diag(k))
#' calcKer.IBS(X)
calcKer.IBS <- function(X){
  n = dim(X)[1]
  p = dim(X)[2]
  mat0 = matrix(NA, n, n)
  for(i in 1:n){
    x1 = X[i,]
    for(j in i:n){
      x2 = X[j,]
      y = x1-x2
      y1 = as.vector(unlist(y))
      z = mean(2-abs(y), na.rm = TRUE)/2
      mat0[i,j] = z
      mat0[j,i] = z
      print(c(i,j))
    }
  }
  return(mat0)
}


