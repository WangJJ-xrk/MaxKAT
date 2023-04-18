#' calc.skew calculates the skew for a vector
#'
#' @param x the vector whose skew is to be calcualted
#'
#' @return the skew
#' @export
calc.skew <- function(x){
  m0 = mean(x)
  v0 = var(x)
  r = mean((x-m0)^3)/(v0^(3/2))
  r
}


#' MaxKAT A maximum kernel-based association test to detect the
#' pleiotropic genetic effects on multiple phenotypes
#'
#' @param E the similarity matrix for genetic data
#' @param G the similarity matrix for phenotypic data
#' @param kVec the index set for E
#' @param mVec the index set for G
#' @param numPer the number for permutation
#' @param ifGEV the method used for significance calculation, must be one
#' "NULL","TRUE" and "Both". If ifGEV = NULL, the p value based on permutation
#' is returned. If ifGEV = TRUE, the p value based on GEV approximation
#' is returned. If ifGEV = "Both", both p values based on permutation and GVE
#' approximation is returned.
#' @importFrom extRemes pevd
#' @importFrom stats dist
#' @importFrom stats median
#' @importFrom stats var
#' @return the p value for MaxKAT
#' @export
#'
#' @examples
#' library(MASS)
#' library(ismev)
#' library(extRemes)
#' n = 100
#' p = 10
#' k = 10
#' betax = matrix(rnorm(k*p,0,1),k,p)
#' x = mvrnorm(n, rep(0,k), diag(k))
#' y = x%*%betax + mvrnorm(n, rep(0,p), diag(p))
#' G = calcKerMat.gau(x)
#' E = calcKerMat.gau(y)
#' kVec = c(1/2,3/4,1,2,3)
#' mVec = c(1/2,3/4,1,2,3)
#' MaxKAT(E,G,kVec,mVec)
MaxKAT <- function(E,G,kVec,mVec,numPer = 1000, ifGEV = TRUE){
  n = dim(E)[1]
  k1 = length(kVec)
  k2 = length(mVec)

  eigE = eigen(E)
  valE0 = eigE$values
  kE = sum(valE0>0.0001)
  valE = valE0[1:kE]
  vecE = eigE$vectors[,1:kE]

  eigG = eigen(G)
  valG0 = eigG$values
  kG = sum(valG0>0.0001)
  valG = valG0[1:kG]
  vecG = eigG$vectors[,1:kG]

  staMat = matrix(NA, k1, k2)

  for(i in 1:k1){
    valE1 = valE^kVec[i]
    E1 = vecE %*% diag(valE1) %*% t(vecE)
    for(j in 1:k2){
      valG1 = valG^mVec[j]
      G1 = vecG %*% diag(valG1) %*% t(vecG)
      sta0 = sum(E1*G1)/(n^(kVec[i]+mVec[j]-1))
      mat0 = outer(valE1, valG1, "*")/(n^(kVec[i]+mVec[j]))
      mu0 = sum(mat0)
      var0 = 2*sum(mat0^2)
      staMat[i,j] = (sta0-mu0)/sqrt(var0)
    }
  }
  max0 = max(staMat)

  maxVec = rep(NA, numPer)
  for(i1 in 1:numPer){
    ind0 = sample(n,n)
    staMatPer = matrix(NA, k1, k2)
    vecE.per = vecE[ind0,]

    for(ii in 1:k1){
      valE1 = valE^kVec[ii]
      E1 = vecE.per %*% diag(valE1) %*% t(vecE.per)
      for(jj in 1:k2){
        valG1 = valG^mVec[jj]
        G1 = vecG %*% diag(valG1) %*% t(vecG)
        sta0 = sum(E1*G1)/(n^(kVec[ii]+mVec[jj]-1))
        mat0 = outer(valE1, valG1, "*")/(n^(kVec[ii]+mVec[jj]))
        mu0 = sum(mat0)
        var0 = 2*sum(mat0^2)
        staMatPer[ii,jj] = (sta0-mu0)/sqrt(var0)
      }
    }
    maxVec[i1] = max(staMatPer)
  }


  if(ifGEV == "Both"){
    TS = maxVec
    r=calc.skew(TS)
    xi = seq(from=-0.3, to=2, by=0.001)[seq(from=-0.2, to=2, by=0.001)!=0]
    r1 = xi/(abs(xi))*(-gamma(1+3*xi)+3*gamma(1+xi)*gamma(1+2*xi)-2*gamma(1+xi)^3)/((gamma(1+2*xi)-gamma(1+xi)^2)^(3/2))
    xiS = xi[which.min(abs(r-r1))]
    alphaS=sqrt(xiS^2*var(TS)/(gamma(1+2*xiS)-gamma(1+xiS)^2))
    betaS = mean(TS)-alphaS/xiS*(1-gamma(1+xiS))
    est0 = c(betaS, alphaS, -xiS)

    pval1 = 1 - pevd(max0, loc = est0[1], scale = est0[2], shape = est0[3], type = "GEV")

    pval0 = sum(max0<maxVec)/numPer

    return(c(pval0,pval1))
  }else if(ifGEV == "TRUE"){
    TS = maxVec
    r=calc.skew(TS)
    xi = seq(from=-0.3, to=2, by=0.001)[seq(from=-0.2, to=2, by=0.001)!=0]
    r1 = xi/(abs(xi))*(-gamma(1+3*xi)+3*gamma(1+xi)*gamma(1+2*xi)-2*gamma(1+xi)^3)/((gamma(1+2*xi)-gamma(1+xi)^2)^(3/2))
    xiS = xi[which.min(abs(r-r1))]
    alphaS=sqrt(xiS^2*var(TS)/(gamma(1+2*xiS)-gamma(1+xiS)^2))
    betaS = mean(TS)-alphaS/xiS*(1-gamma(1+xiS))
    est0 = c(betaS, alphaS, -xiS)

    pval1 = 1 - pevd(max0, loc = est0[1], scale = est0[2], shape = est0[3], type = "GEV")

    return(c(pval1))
  }else{
    pval0 = sum(max0<maxVec)/numPer
    return(pval0)
  }
}

#' KAT Kernel based association test
#'
#' @param E the similarity matrix for genetic data
#' @param G the similarity matrix for phenotypic data
#' @param B1 permutation time
#'
#' @return the p value
#' @export
#'
#' @examples
#' library(MASS)
#' n = 100
#' p = 10
#' k = 10
#' betax = matrix(rnorm(k*p,0,1),k,p)
#' x = mvrnorm(n, rep(0,k), diag(k))
#' y = x%*%betax + mvrnorm(n, rep(0,p), diag(p))
#' G = calcKerMat.gau(x)
#' E = calcKerMat.gau(y)
#' KAT(E,G,100)
KAT <- function(E, G, B1 = 1000){
  n = dim(E)[1]
  sta0 = sum(E*G)
  staVec = rep(NA, B1)
  for(i in 1:B1){
    ind1 = sample(n,n)
    E1 = E[ind1, ind1]
    staVec[i] = sum(E1*G)
  }
  pval = sum(sta0<staVec)/B1
  return(pval)
}
