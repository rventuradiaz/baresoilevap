
#' Partial first derivative in backward difference approximation
#'
#' @param j
#' @param Y
#' @param X
#'
#' @return
#' @export
#'
#' @examples
firstderivative_backwardiff <- function (j=4
                                         , Y=matrix(rep(0,4),nrow=4, ncol=1)
                                         , X=matrix(rep(0,4),nrow=4, ncol=1)){
  nrow <- dim(X)[[1]]

  jminus1 <- if(j>1){j-1}else {1}
  paddingj <- matrix(rep(0,nrow-j), nrow = nrow-j, ncol = 1)
  paddingjminus1 <- matrix(rep(0,nrow-jminus1), nrow = nrow-jminus1, ncol = 1)

  xj <- rbind(paddingj ,X[1:j,,drop = FALSE])
  xjminus1 <- rbind(paddingjminus1 ,X[1:jminus1,,drop = FALSE])

  yj <- rbind(paddingj ,Y[1:j,,drop = FALSE])
  yjminus1 <- rbind(paddingjminus1 ,Y[1:jminus1,,drop = FALSE])

  dyj <- yj-yjminus1
  dxj <- xj-xjminus1

  dydx <- dyj
  if(abs(dyj[nrow])>0){
    dydx <- mapply(function(x,y)x/y, dydx , dxj)
    return(dydx)
  }else{
      return(matrix(rep(0.0,nrow,nrow = nrow, ncol = 1)))
    }

  return(dydx)


}

#' Partial second derivative in backward difference approximation
#'
#' @param j : index
#' @param Y : Dependent variable
#' @param X : Independent variable
#'
#' @return
#' @export
#'
#' @examples
secondderivative_backwarddiff <- function(j=4
                                          , Y=matrix(rep(0,4),nrow=4, ncol=1)
                                          , X=matrix(rep(0,4),nrow=4, ncol=1))
  {
    nrow <- dim(X)[[1]]

    jminus1 <- if(j>1){j-1}else {1}
    jminus2 <- if(j>2){j-2}else{1}
    jminus3 <- if(j>3){j-3}else{1}

    paddingj <- matrix(rep(0,nrow-j), nrow = nrow-j, ncol = 1)
    paddingjminus1 <- matrix(rep(0,nrow-jminus1), nrow = nrow-jminus1, ncol = 1)
    paddingjminus2 <- matrix(rep(0,nrow-jminus2), nrow = nrow-jminus2, ncol = 1)
    paddingjminus3 <- matrix(rep(0,nrow-jminus3), nrow = nrow-jminus3, ncol = 1)


    xj <- rbind(paddingj ,X[1:j,,drop = FALSE])
    xjminus1 <- rbind(paddingjminus1 ,X[1:jminus1,,drop = FALSE])

    yj <- rbind(paddingj ,Y[1:j,,drop = FALSE])
    yjminus1 <- rbind(paddingjminus1 ,Y[1:jminus1,,drop = FALSE])
    yjminus2 <- rbind(paddingjminus2 ,Y[1:jminus2,,drop = FALSE])
    yjminus3 <- rbind(paddingjminus3 ,Y[1:jminus3,,drop = FALSE])

    dxj <- xj-xjminus1

    dy2dx2 <- yj-2*yjminus1+yjminus2
    if (dy2dx2[nrow] == 0 ) {return(matrix(rep(0.0,nrow),nrow = nrow, ncol = 1))} else{
      dy2dx2 <- mapply(function(x,y)x/y, yj-2*yjminus1+yjminus2, dxj)
      dy2dx2 <- mapply(function(x,y)x/y, dy2dx2, dxj)
      return(dy2dx2)  # Devolver el último valor del vector: si

    }
}


# Flux in liquid phase

#' Product of coefficient of diffusion of volumetric humidity in liquid phase
#'
#' @param i
#' @param j
#' @param X
#' @param tim
#' @param thetasat
#'
#' @return
#' @export
#'
#' @examples
prod_Dthetaliq_d2thetadz <- function (i,j, X=matrix(NA,nrow = 4, ncol = 4), tim, thetasat = 0.1){
  nrow <- dim(X)[[1]]
  paddingzj <- matrix(rep(0,nrow-j), nrow = nrow-j, ncol = 1)
  jminus1 <- if(j>1){j-1}else {1}
  jminus2 <- if(j>2){j-2}else{1}
  jminus3 <- if(j>3){j-3}else{1}
  paddingzminus1 <- matrix(rep(0,nrow-jminus1), nrow = nrow-jminus1, ncol = 1)
  paddingzminus2 <- matrix(rep(0,nrow-jminus2), nrow = nrow-jminus2, ncol = 1)
  paddingzminus3 <- matrix(rep(0,nrow-jminus3), nrow = nrow-jminus3, ncol = 1)
  z <- X[,1,drop = FALSE]
  zj <- rbind(paddingzj ,z[1:j,,drop = FALSE])
  zminus1 <- rbind(paddingzminus1 ,z[1:jminus1,,drop = FALSE])
  zminus2 <- rbind(paddingzminus2 ,z[1:jminus2,,drop = FALSE])
  dzj <- zj-zminus1
  dzjminus1 <- zminus1-zminus2

  # volumetric water content (theta)

  thetai <- X[,i, drop = FALSE]
  thetaij <- rbind(paddingzj, thetai[1:j, ,drop = FALSE])
  thetaiminus1 <- X[,if(i>1){i-1}else{1}, drop = FALSE]
  thetaijminus1 <- rbind(paddingzminus1 ,thetai[1:jminus1,,drop = FALSE])
  thetaijminus2 <- rbind(paddingzminus2 ,thetai[1:jminus2,,drop = FALSE])
  thetaijminus3 <- rbind(paddingzminus3 ,thetai[1:jminus3,,drop = FALSE])
  thetaiminus1j <- rbind(paddingzj, thetaiminus1[1:j, ,drop = FALSE])
  diftheta <- thetaij-thetaijminus1

  # D(theta,liq)
  dthetaliqij <- apply (thetaij, MARGIN = 1, FUN = baresoilevap::Dthetaliq)
  dthetaliqijminus1 <- apply (thetaijminus1, MARGIN = 1, FUN = baresoilevap::Dthetaliq)

  # AkBk <- Akm1Bk + AkBkm1 - Akm1Bkm1 : linearization using Newton formula

  Akm1Bk <- mapply(function(x,y)x*y, dthetaliqijminus1, secondderivative_backwarddiff(j,thetai, z))
  AkBkm1 <- mapply(function(x,y)x*y, dthetaliqij, secondderivative_backwarddiff(j-1,thetai, z))
  Akm1Bkm1 <- mapply(function(x,y)x*y, dthetaliqijminus1, secondderivative_backwarddiff(j-1,thetai, z))

  AkBk <- Akm1Bk+AkBkm1-Akm1Bkm1
  return (AkBk)
}


#' Derivative of mixed terms theta and Dthetaliq
#'
#' @param i
#' @param j
#' @param X
#' @param tim
#'
#' @return
#' @export
#'
#' @examples
prod_dthetadz_dDtheliqdz <- function (i, j,X=matrix(NA,nrow = 4, ncol = 4),tim){
  nrow <- dim(X)[[1]]
  paddingzj <- matrix(rep(0,nrow-j), nrow = nrow-j, ncol = 1)
  jminus1 <- if(j>1){j-1}else {1}
  jminus2 <- if(j>2){j-2}else{1}
  jminus3 <- if(j>3){j-3}else{1}
  paddingzminus1 <- matrix(rep(0,nrow-jminus1), nrow = nrow-jminus1, ncol = 1)
  paddingzminus2 <- matrix(rep(0,nrow-jminus2), nrow = nrow-jminus2, ncol = 1)
  paddingzminus3 <- matrix(rep(0,nrow-jminus3), nrow = nrow-jminus3, ncol = 1)
  z <- X[,1,drop = FALSE]
  zj <- rbind(paddingzj ,z[1:j,,drop = FALSE])
  zminus1 <- rbind(paddingzminus1 ,z[1:jminus1,,drop = FALSE])
  zminus2 <- rbind(paddingzminus2 ,z[1:jminus2,,drop = FALSE])
  dzj <- zj-zminus1
  dzjminus1 <- zminus1-zminus2

  # volumetric water content (theta)
  thetasat <- 0.1
  thetai <- X[,i, drop = FALSE]
  thetaij <- rbind(paddingzj, thetai[1:j, ,drop = FALSE])
  thetaiminus1 <- X[,if(i>1){ i-1}else{1}, drop = FALSE]
  thetaijminus1 <- rbind(paddingzminus1 ,thetai[1:jminus1,,drop = FALSE])
  thetaijminus2 <- rbind(paddingzminus2 ,thetai[1:jminus2,,drop = FALSE])
  thetaijminus3 <- rbind(paddingzminus3 ,thetai[1:jminus3,,drop = FALSE])
  thetaiminus1j <- rbind(paddingzj, thetaiminus1[1:j, ,drop = FALSE])
  diftheta <- thetaij-thetaijminus1

  # D(theta,liq)
  dthetaliqij <- apply (thetaij, MARGIN = 1, FUN = baresoilevap::Dthetaliq)
  dthetaliqijminus1 <- apply (thetaijminus1, MARGIN = 1, FUN = baresoilevap::Dthetaliq)
  dthetaliqijminus2 <- apply (thetaijminus2, MARGIN = 1, FUN = baresoilevap::Dthetaliq)

  if (j < 3){
    # use lagging formula and theta i, j-1 = theta i,j-2
    dthetadzdDdz <- matrix(rep(0.0,nrow), nrow = nrow, ncol = 1)
    result <- dthetadzdDdz
  } else {
    # use newton formula

    Akm1Bk <- mapply(function(x,y){x/y},(thetaijminus1-thetaijminus2),dzjminus1)
    Akm1Bk <- mapply(function(x,y){x*y},Akm1Bk,(dthetaliqij - dthetaliqijminus1))
    Akm1Bk <- mapply(function(x,y){x/y},Akm1Bk,dzjminus1)
    dthetadzdDdz <- Akm1Bk
    AkBkm1 <- mapply(function(x,y){x/y},(thetaij-thetaijminus1),dzj)
    AkBkm1 <- mapply(function(x,y){x*y},AkBkm1,(dthetaliqijminus1 - dthetaliqijminus2))
    AkBkm1 <- mapply(function(x,y){x/y},AkBkm1,dzjminus1)
    dthetadzdDdz <- dthetadzdDdz + AkBkm1
    Akm1Bkm1 <- mapply(function(x,y){x/y},(thetaijminus1-thetaijminus2),dzjminus1)
    Akm1Bkm1 <- mapply(function(x,y){x*y},Akm1Bkm1,(dthetaliqijminus1 - dthetaliqijminus2))
    Akm1Bkm1 <- mapply(function(x,y){x/y},Akm1Bkm1,dzjminus1)
    dthetadzdDdz <- dthetadzdDdz - Akm1Bkm1
    result <- dthetadzdDdz
  }
  return(result)
}

#' Calculation of second partial derivative of temperature multipied by diffusion coefficient of water
#'
#' @param i
#' @param j
#' @param X
#' @param Y
#' @param tim
#'
#' @return
#' @export
#'
#' @examples
prod_dtliq_d2Tdz2 <- function(i, j,X=matrix(NA,nrow = 4, ncol = 4), Y = matrix(NA,nrow = 4, ncol = 4),tim){
  nrow <- dim(X)[[1]]
  paddingzj <- matrix(rep(0,nrow-j), nrow = nrow-j, ncol = 1)
  jminus1 <- if(j>1){j-1}else {1}
  jminus2 <- if(j>2){j-2}else{1}
  jminus3 <- if(j>3){j-3}else{1}
  paddingzminus1 <- matrix(rep(0,nrow-jminus1), nrow = nrow-jminus1, ncol = 1)
  paddingzminus2 <- matrix(rep(0,nrow-jminus2), nrow = nrow-jminus2, ncol = 1)
  paddingzminus3 <- matrix(rep(0,nrow-jminus3), nrow = nrow-jminus3, ncol = 1)
  z <- X[,1,drop = FALSE]
  zj <- rbind(paddingzj ,z[1:j,,drop = FALSE])
  zminus1 <- rbind(paddingzminus1 ,z[1:jminus1,,drop = FALSE])
  zminus2 <- rbind(paddingzminus2 ,z[1:jminus2,,drop = FALSE])
  dzj <- zj-zminus1
  dzjminus1 <- zminus1-zminus2

  # volumetric water content (theta)
  thetasat <- 0.1
  thetai <- X[,i, drop = FALSE]
  thetaij <- rbind(paddingzj, thetai[1:j, ,drop = FALSE])
  thetaiminus1 <- X[,if(i>1){i-1}else{1}, drop = FALSE]
  thetaijminus1 <- rbind(paddingzminus1 ,thetai[1:jminus1,,drop = FALSE])
  thetaijminus2 <- rbind(paddingzminus2 ,thetai[1:jminus2,,drop = FALSE])
  thetaijminus3 <- rbind(paddingzminus3 ,thetai[1:jminus3,,drop = FALSE])
  thetaiminus1j <- rbind(paddingzj, thetaiminus1[1:j, ,drop = FALSE])
  diftheta <- thetaij-thetaijminus1

  # temperature
  ti <- Y[,i, drop = FALSE]
  tij <- rbind(paddingzj, ti[1:j, ,drop = FALSE])
  timinus1 <- Y[,if(i>1) {i-1}else{1}, drop = FALSE]
  tijminus1 <- rbind(paddingzminus1 ,ti[1:jminus1,,drop = FALSE])
  timinus1j <- rbind(paddingzj, timinus1[1:j, ,drop = FALSE])
  tijminus2 <- rbind(paddingzminus2 ,ti[1:jminus2,,drop = FALSE])
  tijminus3 <- rbind(paddingzminus3 ,ti[1:jminus3,,drop = FALSE])

  # Diffusion coefficient for liquid water
  dtliqij <- mapply (baresoilevap::Dtliq, thetaij, tij)
  dtliqijminus1 <- mapply (baresoilevap::Dtliq, thetaijminus1, tijminus1)

  if (j<4) {
    # Using lagging formula
    dtliqd2tdz2 <- mapply(function(x, y){x * y}, dtliqijminus1, (tij - 2 * tijminus1 + tijminus2))
    dtliqd2tdz2 <- mapply(function(x, y){x / y}, dtliqd2tdz2 , dzj)
    dtliqd2tdz2 <- mapply(function(x, y){x / y}, dtliqd2tdz2 , dzj)
    dtliqd2tdz2 <- mapply(function(x,y)if(x==0.0){0.0}else{y}, (tij - 2 * tijminus1 + tijminus2),dtliqd2tdz2)
    result <- dtliqd2tdz2

  } else {
    # Using Newton formula
    Akm1Bk <- mapply(function(x, y){x * y}, dtliqijminus1, (tij - 2 * tijminus1 + tijminus2))
    Akm1Bk <- mapply(function(x, y){x / y}, Akm1Bk , dzj)
    Akm1Bk <- mapply(function(x, y){x / y}, Akm1Bk , dzj)
    dtliqd2tdz2 <- Akm1Bk
    AkBkm1 <- mapply(function(x, y){x * y}, dtliqij, (tijminus1 - 2 * tijminus2 + tijminus3))
    AkBkm1 <- mapply(function(x, y){x / y}, AkBkm1 , dzjminus1)
    AkBkm1 <- mapply(function(x, y){x / y}, AkBkm1 , dzjminus1)
    dtliqd2tdz2 <- dtliqd2tdz2 + AkBkm1
    Akm1Bkm1 <- mapply(function(x, y){x * y}, dtliqijminus1, (tijminus1 - 2 * tijminus2 + tijminus3))
    Akm1Bkm1 <- mapply(function(x, y){x / y}, Akm1Bkm1 , dzjminus1)
    Akm1Bkm1 <- mapply(function(x, y){x / y}, Akm1Bkm1 , dzjminus1)
    dtliqd2tdz2 <- dtliqd2tdz2 - Akm1Bkm1
    result <- dtliqd2tdz2
  }
  return(result)
}


#' Calculation of product of partial derivative of temperatura to z and partial derivative of liquid diffusion coefficient to z
#'
#' @param i
#' @param j
#' @param X
#' @param Y
#' @param tim
#'
#' @return
#' @export
#'
#' @examples
prod_dtdz_dDtlidz <- function (i, j,X=matrix(NA,nrow = 4, ncol = 4), Y = matrix(NA,nrow = 4, ncol = 4),tim){
  nrow <- dim(X)[[1]]
  paddingzj <- matrix(rep(0,nrow-j), nrow = nrow-j, ncol = 1)
  jminus1 <- if(j>1){j-1}else {1}
  jminus2 <- if(j>2){j-2}else{1}
  jminus3 <- if(j>3){j-3}else{1}
  paddingzminus1 <- matrix(rep(0,nrow-jminus1), nrow = nrow-jminus1, ncol = 1)
  paddingzminus2 <- matrix(rep(0,nrow-jminus2), nrow = nrow-jminus2, ncol = 1)
  paddingzminus3 <- matrix(rep(0,nrow-jminus3), nrow = nrow-jminus3, ncol = 1)
  z <- X[,1,drop = FALSE]
  zj <- rbind(paddingzj ,z[1:j,,drop = FALSE])
  zminus1 <- rbind(paddingzminus1 ,z[1:jminus1,,drop = FALSE])
  zminus2 <- rbind(paddingzminus2 ,z[1:jminus2,,drop = FALSE])
  dzj <- zj-zminus1
  dzjminus1 <- zminus1-zminus2

  # volumetric water content (theta)
  thetasat <- 0.1
  thetai <- X[,i, drop = FALSE]
  thetaij <- rbind(paddingzj, thetai[1:j, ,drop = FALSE])
  thetaiminus1 <- X[,if(i>1){i-1}else{1}, drop = FALSE]
  thetaijminus1 <- rbind(paddingzminus1 ,thetai[1:jminus1,,drop = FALSE])
  thetaijminus2 <- rbind(paddingzminus2 ,thetai[1:jminus2,,drop = FALSE])
  thetaijminus3 <- rbind(paddingzminus3 ,thetai[1:jminus3,,drop = FALSE])
  thetaiminus1j <- rbind(paddingzj, thetaiminus1[1:j, ,drop = FALSE])
  diftheta <- thetaij-thetaijminus1

  # temperature
  ti <- Y[,i, drop = FALSE]
  tij <- rbind(paddingzj, ti[1:j, ,drop = FALSE])
  timinus1 <- Y[,if(i>1){i-1}else{1}, drop = FALSE]
  tijminus1 <- rbind(paddingzminus1 ,ti[1:jminus1,,drop = FALSE])
  timinus1j <- rbind(paddingzj, timinus1[1:j, ,drop = FALSE])
  tijminus2 <- rbind(paddingzminus2 ,ti[1:jminus2,,drop = FALSE])
  tijminus3 <- rbind(paddingzminus3 ,ti[1:jminus3,,drop = FALSE])

  # Diffusion coefficient for liquid water
  dtliqij <- mapply (baresoilevap::Dtliq, thetaij, tij)
  dtliqijminus1 <- mapply (baresoilevap::Dtliq, thetaijminus1, tijminus1)
  dtliqijminus2 <- mapply (baresoilevap::Dtliq, thetaijminus2, tijminus2)

  if (j< 3) {
    # use lagging formula
    Akm1Bk <- mapply(function(x, y){x/y}, tijminus1-tijminus2, dzjminus1)
    Akm1Bk <- mapply(function(x, y){x*y}, Akm1Bk, dtliqij-dtliqijminus1)
    Akm1Bk <- mapply(function(x, y){x/y}, Akm1Bk, dzj)
    result <- Akm1Bk
  } else {
    Akm1Bk <- mapply(function(x, y){x/y}, tijminus1-tijminus2, dzjminus1)
    Akm1Bk <- mapply(function(x, y){x*y}, Akm1Bk, dtliqij-dtliqijminus1)
    Akm1Bk <- mapply(function(x, y){x/y}, Akm1Bk, dzj)
    dtdzddtliqdz <- Akm1Bk
    AkBkm1 <- mapply(function(x, y){x/y}, tij-tijminus1, dzj)
    AkBkm1 <- mapply(function(x, y){x*y}, AkBkm1, dtliqijminus1-dtliqijminus2)
    AkBkm1 <- mapply(function(x, y){x/y}, AkBkm1, dzjminus1)
    dtdzddtliqdz <- dtdzddtliqdz + AkBkm1
    Akm1Bkm1 <- mapply(function(x, y){x/y}, tijminus1-tijminus2, dzjminus1)
    Akm1Bkm1 <- mapply(function(x, y){x*y}, Akm1Bkm1, dtliqijminus1-dtliqijminus2)
    Akm1Bkm1 <- mapply(function(x, y){x/y}, Akm1Bkm1, dzjminus1)
    dtdzddtliqdz <- dtdzddtliqdz - Akm1Bkm1
    result <- dtdzddtliqdz
  }
  return(result)
}

#' Calculation of partial derivative of K to z
#'
#' @param i
#' @param j
#' @param X
#' @param tim
#'
#' @return
#' @export
#'
#' @examples
dK_dz <- function (i, j,X=matrix(NA,nrow = 4, ncol = 4),tim){
  nrow <- dim(X)[[1]]
  paddingzj <- matrix(rep(0,nrow-j), nrow = nrow-j, ncol = 1)
  jminus1 <- if(j>1){j-1}else {1}
  jminus2 <- if(j>2){j-2}else{1}
  jminus3 <- if(j>3){j-3}else{1}
  paddingzminus1 <- matrix(rep(0,nrow-jminus1), nrow = nrow-jminus1, ncol = 1)
  z <- X[,1,drop = FALSE]
  zj <- rbind(paddingzj ,z[1:j,,drop = FALSE])
  zminus1 <- rbind(paddingzminus1 ,z[1:jminus1,,drop = FALSE])
  dzj <- zj-zminus1

  # volumetric water content (theta)
  thetasat <- 0.1
  thetai <- X[,i, drop = FALSE]
  thetaij <- rbind(paddingzj, thetai[1:j, ,drop = FALSE])
  thetaiminus1 <- X[,if(i>1){i-1}else{1}, drop = FALSE]
  thetaijminus1 <- rbind(paddingzminus1 ,thetai[1:jminus1,,drop = FALSE])
  diftheta <- thetaij-thetaijminus1

  # Calculation of K
  kthetaij <- apply(thetaij, MARGIN = 1, FUN = baresoilevap::K_theta)
  kthetaijminus1 <- apply(thetaijminus1, MARGIN = 1, FUN = baresoilevap::K_theta)

  dkdz <- mapply(function(x,y){x/y}, kthetaij-kthetaijminus1, dzj)

  result <- dkdz
  return(result)
}

#' Calculation of derivative of flux of liquid water (Q theta,liq)
#'
#' @param i : time node number
#' @param j : spatial node number
#' @param X : mesh for theta: row) z ; column) time ; value: volumetric water content (theta)
#' @param tim : time nodes
#'
#' @return
#' @export
#'
#' @examples
dQthetliq_dz <- function (i
                          ,j
                          , rhowater=1.00
                          , X=matrix(NA,nrow = 4, ncol = 4)
                          , Y=matrix(NA,nrow = 4, ncol = 4)
                          , Z=matrix(NA,nrow = 4, ncol = 4),tim, thetasat){
  nrow <- dim(X)[[1]]
  paddingzj <- matrix(rep(0,nrow-j), nrow = nrow-j, ncol = 1)
  jminus1 <- if(j>1){j-1}else {1}
  jminus2 <- if(j>2){j-2}else{1}
  jminus3 <- if(j>3){j-3}else{1}
  paddingzminus1 <- matrix(rep(0,nrow-jminus1), nrow = nrow-jminus1, ncol = 1)
  paddingzminus2 <- matrix(rep(0,nrow-jminus2), nrow = nrow-jminus2, ncol = 1)
  paddingzminus3 <- matrix(rep(0,nrow-jminus3), nrow = nrow-jminus3, ncol = 1)

  # length
  z <- X[,1,drop = FALSE]
  zj <- rbind(paddingzj ,z[1:j,,drop = FALSE])
  zminus1 <- rbind(paddingzminus1 ,z[1:jminus1,,drop = FALSE])
  zminus2 <- rbind(paddingzminus2 ,z[1:jminus2,,drop = FALSE])
  dzj <- zj-zminus1
  dzjminus1 <- zminus1-zminus2

  # Derivative of liquid flux
  dQthetaliqdz <- prod_Dthetaliq_d2thetadz(i,j,X,tim,thetasat)
  dQthetaliqdz <- dQthetaliqdz + prod_dthetadz_dDtheliqdz(i,j,X,tim)
  dQthetaliqdz <- dQthetaliqdz + prod_dtliq_d2Tdz2(i,j,X,Y,tim)
  dQthetaliqdz <- dQthetaliqdz + prod_dtdz_dDtlidz(i,j,X,Y,tim)
  dQthetaliqdz <- dQthetaliqdz +dK_dz(i,j,X,tim)
  dQthetaliqdz <- - (rhowater)*(dQthetaliqdz)
  return(dQthetaliqdz)
}





# Flux in vapor phase




#' Product of density and molecular diffusion coefficient
#'
#' @param i
#' @param j
#' @param Y
#' @param tim
#' @param pij
#'
#' @return
#' @export
#'
#' @examples
prod_rho_datm <- function(i, j, Y = matrix(NA,nrow = 4, ncol = 4),tim, pij) {

  nrow <- dim(Y)[[1]]
  paddingzj <- matrix(rep(0,nrow-j), nrow = nrow-j, ncol = 1)
  jminus1 <- if(j>1){j-1}else {1}
  jminus2 <- if(j>2){j-2}else{1}
  jminus3 <- if(j>3){j-3}else{1}
  paddingzminus1 <- matrix(rep(0,nrow-jminus1), nrow = nrow-jminus1, ncol = 1)
  paddingzminus2 <- matrix(rep(0,nrow-jminus2), nrow = nrow-jminus2, ncol = 1)
  paddingzminus3 <- matrix(rep(0,nrow-jminus3), nrow = nrow-jminus3, ncol = 1)
  z <- Y[,1,drop = FALSE]
  zj <- rbind(paddingzj ,z[1:j,,drop = FALSE])
  zminus1 <- rbind(paddingzminus1 ,z[1:jminus1,,drop = FALSE])
  zminus2 <- rbind(paddingzminus2 ,z[1:jminus2,,drop = FALSE])
  dzj <- zj-zminus1
  dzjminus1 <- zminus1-zminus2

  # temperature
  ti <- Y[,i, drop = FALSE]
  tij <- rbind(paddingzj, ti[1:j, ,drop = FALSE])
  timinus1 <- Y[,if(i>1){i-1}else{1}, drop = FALSE]
  tijminus1 <- rbind(paddingzminus1 ,ti[1:jminus1,,drop = FALSE])
  timinus1j <- rbind(paddingzj, timinus1[1:j, ,drop = FALSE])
  tijminus2 <- rbind(paddingzminus2 ,ti[1:jminus2,,drop = FALSE])
  tijminus3 <- rbind(paddingzminus3 ,ti[1:jminus3,,drop = FALSE])

  # air density
  rhoairij <- mapply(baresoilevap::rho_air, tij, pij )
  rhoairijminus1 <- mapply(baresoilevap::rho_air, tijminus1, pij )
  rhoairiminus1j <- mapply(baresoilevap::rho_air, timinus1j, pij )

  # Molecular coefficient of diffusion Dtam
  datmij <- mapply(baresoilevap::Datm, tij, pij )
  datmijminus1 <- mapply(baresoilevap::Datm, tijminus1, pij )

  Akm1Bk <- mapply(function(x, y){x*y}, rhoairijminus1, datmij)
  AkBkm1 <- mapply(function(x, y){x*y}, rhoairij, datmijminus1)
  Akm1Bkm1 <- mapply(function(x, y){x*y}, rhoairijminus1, datmijminus1)
  prodhrodatm <- Akm1Bk+ AkBkm1-Akm1Bkm1
  return(prodhrodatm)
}

#' product of density and porosity factor
#'
#' @param i
#' @param j
#' @param X
#' @param Y
#' @param tim
#' @param pij
#' @param thetasat
#'
#' @return
#' @export
#'
#' @examples
prod_ftheta_rho <- function (i, j,X=matrix(NA,nrow = 4, ncol = 4), Y = matrix(NA,nrow = 4, ncol = 4),tim, pij, thetasat){
  nrow <- dim(X)[[1]]
  paddingzj <- matrix(rep(0,nrow-j), nrow = nrow-j, ncol = 1)
  jminus1 <- if(j>1){j-1}else {1}
  jminus2 <- if(j>2){j-2}else{1}
  jminus3 <- if(j>3){j-3}else{1}
  paddingzminus1 <- matrix(rep(0,nrow-jminus1), nrow = nrow-jminus1, ncol = 1)
  z <- X[,1,drop = FALSE]
  zj <- rbind(paddingzj ,z[1:j,,drop = FALSE])
  zminus1 <- rbind(paddingzminus1 ,z[1:jminus1,,drop = FALSE])
  dzj <- zj-zminus1

  # volumetric water content (theta)
  thetasat <- 0.1
  thetai <- X[,i, drop = FALSE]
  thetaij <- rbind(paddingzj, thetai[1:j, ,drop = FALSE])
  thetaiminus1 <- X[,if(i>1){i-1}else{1}, drop = FALSE]
  thetaijminus1 <- rbind(paddingzminus1 ,thetai[1:jminus1,,drop = FALSE])
  diftheta <- thetaij-thetaijminus1

  # temperature
  ti <- Y[,i, drop = FALSE]
  tij <- rbind(paddingzj, ti[1:j, ,drop = FALSE])
  timinus1 <- Y[,if(i>1){i-1}else{1}, drop = FALSE]
  tijminus1 <- rbind(paddingzminus1 ,ti[1:jminus1,,drop = FALSE])
  timinus1j <- rbind(paddingzj, timinus1[1:j, ,drop = FALSE])
  tijminus2 <- rbind(paddingzminus2 ,ti[1:jminus2,,drop = FALSE])
  tijminus3 <- rbind(paddingzminus3 ,ti[1:jminus3,,drop = FALSE])

  # air density
  rhoairij <- mapply(baresoilevap::rho_air, tij, pij )
  rhoairijminus1 <- mapply(baresoilevap::rho_air, tijminus1, pij )
  rhoairiminus1j <- mapply(baresoilevap::rho_air, timinus1j, pij )

  # Calculation of f(theta); porosity and tortuosity factor
  fthetaij <- mapply(baresoilevap::porosity_factor, thetaij, thetasat)
  fthetaijminus1 <- mapply(baresoilevap::porosity_factor, thetaijminus1, thetasat)

  Akm1Bk <- mapply(function(x, y){x*y}, rhoairijminus1, fthetaij)
  AkBkm1 <- mapply(function(x, y){x*y}, rhoairij, fthetaijminus1)
  Akm1Bkm1 <- mapply(function(x, y){x*y}, rhoairijminus1, fthetaijminus1)
  prodhrodatm <- Akm1Bk+ AkBkm1-Akm1Bkm1
  return(prodhrodatm)



}


#' Product of porosity factor and molecular diffusion coefficient
#'
#' @param i
#' @param j
#' @param X
#' @param Y
#' @param tim
#' @param pij
#' @param thetasat
#'
#' @return
#' @export
#'
#' @examples
prod_ftheta_datm <- function(i, j,X=matrix(NA,nrow = 4, ncol = 4), Y = matrix(NA,nrow = 4, ncol = 4),tim, pij, thetasat){

  nrow <- dim(X)[[1]]
  paddingzj <- matrix(rep(0,nrow-j), nrow = nrow-j, ncol = 1)
  jminus1 <- if(j>1){j-1}else {1}
  jminus2 <- if(j>2){j-2}else{1}
  jminus3 <- if(j>3){j-3}else{1}
  paddingzminus1 <- matrix(rep(0,nrow-jminus1), nrow = nrow-jminus1, ncol = 1)
  z <- X[,1,drop = FALSE]
  zj <- rbind(paddingzj ,z[1:j,,drop = FALSE])
  zminus1 <- rbind(paddingzminus1 ,z[1:jminus1,,drop = FALSE])
  dzj <- zj-zminus1

  # volumetric water content (theta)
  thetasat <- 0.1
  thetai <- X[,i, drop = FALSE]
  thetaij <- rbind(paddingzj, thetai[1:j, ,drop = FALSE])
  thetaiminus1 <- X[,if(i>1){i-1}else{1}, drop = FALSE]
  thetaijminus1 <- rbind(paddingzminus1 ,thetai[1:jminus1,,drop = FALSE])
  diftheta <- thetaij-thetaijminus1

# temperature
  ti <- Y[,i, drop = FALSE]
  tij <- rbind(paddingzj, ti[1:j, ,drop = FALSE])
  timinus1 <- Y[,if(i>1){i-1}else{1}, drop = FALSE]
  tijminus1 <- rbind(paddingzminus1 ,ti[1:jminus1,,drop = FALSE])
  timinus1j <- rbind(paddingzj, timinus1[1:j, ,drop = FALSE])
  tijminus2 <- rbind(paddingzminus2 ,ti[1:jminus2,,drop = FALSE])
  tijminus3 <- rbind(paddingzminus3 ,ti[1:jminus3,,drop = FALSE])

  # Calculation of f(theta); porosity and tortuosity factor
  fthetaij <- mapply(baresoilevap::porosity_factor, thetaij, thetasat)
  fthetaijminus1 <- mapply(baresoilevap::porosity_factor, thetaijminus1, thetasat)

    # Molecular coefficient of diffusion Dtam
  datmij <- mapply(baresoilevap::Datm, tij, pij )
  datmijminus1 <- mapply(baresoilevap::Datm, tijminus1, pij )

  Akm1Bk <- mapply(function(x, y){x*y}, fthetaijminus1, datmij)
  AkBkm1 <- mapply(function(x, y){x*y}, fthetaij, datmijminus1)
  Akm1Bkm1 <- mapply(function(x, y){x*y}, fthetaijminus1, datmijminus1)
  prodfthetadatm <- Akm1Bk+ AkBkm1-Akm1Bkm1
  return(prodfthetadatm)


}

#' Calculation of product of density, molecular diffusion coefficient, and porosity factor
#'
#' @param i
#' @param j
#' @param X
#' @param Y
#' @param tim
#' @param pij
#' @param thetasat
#'
#' @return
#' @export
#'
#' @examples
prod_rho_datm_ftheta <- function (i, j,X=matrix(NA,nrow = 4, ncol = 4), Y = matrix(NA,nrow = 4, ncol = 4),tim, pij, thetasat){

  nrow <- dim(X)[[1]]
  paddingzj <- matrix(rep(0,nrow-j), nrow = nrow-j, ncol = 1)
  jminus1 <- j-1
  paddingzminus1 <- matrix(rep(0,nrow-jminus1), nrow = nrow-jminus1, ncol = 1)
  z <- X[,1,drop = FALSE]
  zj <- rbind(paddingzj ,z[1:j,,drop = FALSE])
  zminus1 <- rbind(paddingzminus1 ,z[1:jminus1,,drop = FALSE])
  dzj <- zj-zminus1

  # volumetric water content (theta)
  thetasat <- 0.1
  thetai <- X[,i, drop = FALSE]
  thetaij <- rbind(paddingzj, thetai[1:j, ,drop = FALSE])
  thetaiminus1 <- X[,if(i>1){i-1}else{1}, drop = FALSE]
  thetaijminus1 <- rbind(paddingzminus1 ,thetai[1:jminus1,,drop = FALSE])
  diftheta <- thetaij-thetaijminus1

  # Calculation of f(theta); porosity and tortuosity factor
  fthetaij <- mapply(baresoilevap::porosity_factor, thetaij, thetasat)
  fthetaijminus1 <- mapply(baresoilevap::porosity_factor, thetaijminus1, thetasat)

  Akm1Bk <- mappy(function(x,y){x*y}, baresoilevap::prod_rho_datm(i,j-1,Y,tim,pij),fthetaij )
  AkBkm1 <- mapply(function(x,y){x*y}, baresoilevap::prod_rho_datm(i,j,Y,tim,pij),fthetaijminus1)
  Akm1Bkm1 <- mapply(function(x,y){x*y}, baresoilevap::prod_rho_datm(i,j-1,Y,tim,pij),fthetaijminus1)

  prodrhodatmftheta <- Akm1Bk + AkBkm1 - Akm1Bkm1

  return (prodrhodatmftheta)

}

#' Calculation of product of density , molecular diffusion coefficient, and derivative of porosity factor to z
#'
#' @param i
#' @param j
#' @param X
#' @param Y
#' @param tim
#' @param pij
#' @param thetasat
#'
#' @return
#' @export
#'
#' @examples
prod_rhodatm_dfthetadz <- function (i
                                    , j
                                    ,X=matrix(NA,nrow = 4, ncol = 4)
                                    , Y = matrix(NA,nrow = 4, ncol = 4)
                                    ,tim
                                    , pij
                                    , thetasat)
  {
  nrow <- dim(X)[[1]]
  paddingzj <- matrix(rep(0,nrow-j), nrow = nrow-j, ncol = 1)
  jminus1 <- if(j>1){j-1}else {1}
  jminus2 <- if(j>2){j-2}else{1}
  jminus3 <- if(j>3){j-3}else{1}
  paddingzminus1 <- matrix(rep(0,nrow-jminus1), nrow = nrow-jminus1, ncol = 1)
  paddingzminus2 <- matrix(rep(0,nrow-jminus2), nrow = nrow-jminus2, ncol = 1)
  paddingzminus3 <- matrix(rep(0,nrow-jminus3), nrow = nrow-jminus3, ncol = 1)
  z <- X[,1,drop = FALSE]
  zj <- rbind(paddingzj ,z[1:j,,drop = FALSE])
  zminus1 <- rbind(paddingzminus1 ,z[1:jminus1,,drop = FALSE])
  zminus2 <- rbind(paddingzminus2 ,z[1:jminus2,,drop = FALSE])
  dzj <- zj-zminus1
  dzjminus1 <- zminus1-zminus2

  # volumetric water content (theta)
  thetasat <- 0.1
  thetai <- X[,i, drop = FALSE]
  thetaij <- rbind(paddingzj, thetai[1:j, ,drop = FALSE])
  thetaiminus1 <- X[,if(i>1){i-1}else{1}, drop = FALSE]
  thetaijminus1 <- rbind(paddingzminus1 ,thetai[1:jminus1,,drop = FALSE])
  thetaijminus2 <- rbind(paddingzminus2 ,thetai[1:jminus2,,drop = FALSE])
  diftheta <- thetaij-thetaijminus1

  # Calculation of f(theta); porosity and tortuosity factor
  fthetaij <- mapply(baresoilevap::porosity_factor, thetaij, thetasat)
  fthetaijminus1 <- mapply(baresoilevap::porosity_factor, thetaijminus1, thetasat)
  fthetaijminus2 <- mapply(baresoilevap::porosity_factor, thetaijminus2, thetasat)

  Akm1Bk <- mapply(function(x,y){x*y}, baresoilevap::prod_rho_datm(i,j-1,Y = matrix(NA,nrow = 4, ncol = 4),tim, pij) ,(fthetaij-fthetaijminus1) )
  Akm1Bk <- mapply(function(x, y){x/y}, Akm1Bk, dzj)

  AkBkm1 <- mapply(function(x,y){x*y}, baresoilevap::prod_rho_datm(i,j,Y = matrix(NA,nrow = 4, ncol = 4),tim, pij) ,(fthetaijminus1-fthetaijminus2))
  AkBkm1 <- mapply(function(x,y){x*y}, AkBkm1 , dzjminus1)

  Akm1Bkm1 <- mapply(function(x,y){x*y}, baresoilevap::prod_rho_datm(i,j-1,Y = matrix(NA,nrow = 4, ncol = 4),tim, pij) ,(fthetaijminus1-fthetaijminus2))
  Akm1Bkm1 <- mapply(function(x,y){x*y}, AkmBkm1 , dzjminus1)

  if (j > 2 ) {prodrhodatmdfthetadz <- Akm1Bk + AkBkm1 - Akm1Bkm1 } else { prodrhodatmdfthetadz <- Akm1Bk}

  return(prodrhodatmdfthetadz)

  }

#' Product of tortuosity factor and density by derivative of molecular diffusion coefficient
#'
#' @param i
#' @param j
#' @param X
#' @param Y
#' @param tim
#' @param pij
#' @param thetasat
#'
#' @return
#' @export
#'
#' @examples
prod_fthetarho_dDatmdz <- function (i
                                    , j
                                    ,X=matrix(NA,nrow = 4, ncol = 4)
                                    , Y = matrix(NA,nrow = 4, ncol = 4)
                                    ,tim
                                    , pij
                                    , thetasat
                                    ){
  nrow <- dim(X)[[1]]
  paddingzj <- matrix(rep(0,nrow-j), nrow = nrow-j, ncol = 1)
  jminus1 <- if(j>1){j-1}else {1}
  jminus2 <- if(j>2){j-2}else{1}
  jminus3 <- if(j>3){j-3}else{1}
  paddingzminus1 <- matrix(rep(0,nrow-jminus1), nrow = nrow-jminus1, ncol = 1)
  paddingzminus2 <- matrix(rep(0,nrow-jminus2), nrow = nrow-jminus2, ncol = 1)
  paddingzminus3 <- matrix(rep(0,nrow-jminus3), nrow = nrow-jminus3, ncol = 1)
  z <- X[,1,drop = FALSE]
  zj <- rbind(paddingzj ,z[1:j,,drop = FALSE])
  zminus1 <- rbind(paddingzminus1 ,z[1:jminus1,,drop = FALSE])
  zminus2 <- rbind(paddingzminus2 ,z[1:jminus2,,drop = FALSE])
  dzj <- zj-zminus1
  dzjminus1 <- zminus1-zminus2

  # temperature
  ti <- Y[,i, drop = FALSE]
  tij <- rbind(paddingzj, ti[1:j, ,drop = FALSE])
  timinus1 <- Y[,if(i>1){i-1}else{1}, drop = FALSE]
  tijminus1 <- rbind(paddingzminus1 ,ti[1:jminus1,,drop = FALSE])
  timinus1j <- rbind(paddingzj, timinus1[1:j, ,drop = FALSE])
  tijminus2 <- rbind(paddingzminus2 ,ti[1:jminus2,,drop = FALSE])
  tijminus3 <- rbind(paddingzminus3 ,ti[1:jminus3,,drop = FALSE])

  # Molecular coefficient of diffusion Dtam
  datmij <- mapply(baresoilevap::Datm, tij, pij )
  datmijminus1 <-  mapply(baresoilevap::Datm, tijminus1, pij )
  datmijminus2 <-  mapply(baresoilevap::Datm, tijminus2, pij )

  Akm1Bk <- mapply(function(x,y){x*y}
                   ,baresoilevap::prod_ftheta_rho(i, jminus1,X, Y ,tim, pij, thetasat)
                   , datmij-datmijminus1)
                   Akm1Bk <- mapply(function(x,y){x*y}, Akm1Bk, dzj)
                   AkBkm1 <- if (j > 2) {mapply (function(x,y){x*y}
                                                 , baresoilevap::prod_ftheta_rho(i, j,X, Y ,tim, pij, thetasat)
                                                 , datmijminus1-datmijminus2)} else {0}
                   AkBkm1 <- mapply(function(x,y){x/y}, AkBkm1, dzjminus1)
                   Akm1Bkm1 <- if (j >2) {
                     mapply (function(x,y){x*y}
                                                  , baresoilevap::prod_ftheta_rho(i, jminus1,X, Y ,tim, pij, thetasat)
                                                  , datmijminus1-datmijminus2
                                                  )
                     } else {0}
                   Akm1Bkm1 <- mapply(function(x,y){x*y}, Akm1Bkm1, dzjminus1)
                   result <- Akm1Bk+AkBkm1-Akm1Bkm1
                   return(result)
}


#' product of tortuosity factor, molecular diffusion coefficient, and partial derivative of air density relative to z
#'
#' @param i
#' @param j
#' @param X
#' @param Y
#' @param tim
#' @param pij
#' @param thetasat
#'
#' @return
#' @export
#'
#' @examples
prod_fthetaDatm_drhodz <- function(i
                                   , j
                                   ,X=matrix(NA,nrow = 4, ncol = 4)
                                   , Y = matrix(NA,nrow = 4, ncol = 4)
                                   ,tim
                                   , pij
                                   , thetasat
)
  {
  nrow <- dim(X)[[1]]
  paddingzj <- matrix(rep(0,nrow-j), nrow = nrow-j, ncol = 1)
  jminus1 <- if(j>1){j-1}else {1}
  jminus2 <- if(j>2){j-2}else{1}
  jminus3 <- if(j>3){j-3}else{1}
  paddingzminus1 <- matrix(rep(0,nrow-jminus1), nrow = nrow-jminus1, ncol = 1)
  paddingzminus2 <- matrix(rep(0,nrow-jminus2), nrow = nrow-jminus2, ncol = 1)
  paddingzminus3 <- matrix(rep(0,nrow-jminus3), nrow = nrow-jminus3, ncol = 1)
  z <- X[,1,drop = FALSE]
  zj <- rbind(paddingzj ,z[1:j,,drop = FALSE])
  zminus1 <- rbind(paddingzminus1 ,z[1:jminus1,,drop = FALSE])
  zminus2 <- rbind(paddingzminus2 ,z[1:jminus2,,drop = FALSE])
  dzj <- zj-zminus1
  dzjminus1 <- zminus1-zminus2

  # temperature
  ti <- Y[,i, drop = FALSE]
  tij <- rbind(paddingzj, ti[1:j, ,drop = FALSE])
  timinus1 <- Y[,if(i>1){i-1}else{1}, drop = FALSE]
  tijminus1 <- rbind(paddingzminus1 ,ti[1:jminus1,,drop = FALSE])
  timinus1j <- rbind(paddingzj, timinus1[1:j, ,drop = FALSE])
  tijminus2 <- rbind(paddingzminus2 ,ti[1:jminus2,,drop = FALSE])
  tijminus3 <- rbind(paddingzminus3 ,ti[1:jminus3,,drop = FALSE])

  # Molecular coefficient of diffusion Dtam
  datmij <- mapply(baresoilevap::Datm, tij, pij )
  datmijminus1 <-  mapply(baresoilevap::Datm, tijminus1, pij )
  datmijminus2 <-  mapply(baresoilevap::Datm, tijminus2, pij )

  # AkBk <- Akm1Bk + AkBkm1 - Akm1Bkm1
  AkBkm1 <- if (j > 2) {
    mapply(function(x,y) {x*y}
           ,baresoilevap::prod_ftheta_datm(i, j,X, Y ,tim, pij, thetasat)
           , baresoilevap::rho_air(tijminus1,pij)-baresoilevap::rho_air(tijminus2,pij))
  } else {0}
  AkBkm1 <- mapply(function(x,y) {x/y},AkBkm1, dzjminus1)
  Akm1Bk <- mapply(function(x,y){x*y}
                   ,baresoilevap::prod_ftheta_datm(i, jminus1,X, Y ,tim, pij, thetasat)
                   ,baresoilevap::rho_air(tij,pij)-baresoilevap::rho_air(tijminus1,pij))
  Akm1Bk <- mapply(function(x,y){x/y},Akm1Bk, dzj)
  Akm1Bkm1 <- if (j >2) {
    mapply(function(x,y) {x*y},AkBkm1, baresoilevap::rho_air(tijminus1,pij)-baresoilevap::rho_air(tijminus2,pij))
  } else {0}
  Akm1Bkm1 <- mapply(function(x,y){x/y}, Akm1Bkm1, dzjminus1)
  Akm1Bkm1 <- mapply(function(x,y){x*y}, Akm1Bkm1
                     , baresoilevap::prod_ftheta_datm(i, jminus1,X, Y ,tim, pij, thetasat))
  result <- Akm1Bk+AkBkm1-Akm1Bkm1
  return(result)

}


#' Title
#'
#' @param i
#' @param j
#' @param X
#' @param Y
#' @param Z
#' @param tim
#' @param pij
#' @param thetasat
#'
#' @return
#' @export
#'
#' @examples
prod_rho_datm_ftheta_dqdz <- function(i
                                              , j
                                              ,X=matrix(NA,nrow = 4, ncol = 4)
                                              , Y = matrix(NA,nrow = 4, ncol = 4)
                                              , Z = matrix(NA,nrow = 4, ncol = 4)
                                              ,tim
                                              , pij
                                              , thetasat){
  nrow <- dim(X)[[1]]
  paddingzj <- matrix(rep(0,nrow-j), nrow = nrow-j, ncol = 1)
  jminus1 <- if(j>1){j-1}else {1}
  jminus2 <- if(j>2){j-2}else{1}
  jminus3 <- if(j>3){j-3}else{1}
  paddingzminus1 <- matrix(rep(0,nrow-jminus1), nrow = nrow-jminus1, ncol = 1)
  paddingzminus2 <- matrix(rep(0,nrow-jminus2), nrow = nrow-jminus2, ncol = 1)
  paddingzminus3 <- matrix(rep(0,nrow-jminus3), nrow = nrow-jminus3, ncol = 1)
  z <- X[,1,drop = FALSE]
  zj <- rbind(paddingzj ,z[1:j,,drop = FALSE])
  zminus1 <- rbind(paddingzminus1 ,z[1:jminus1,,drop = FALSE])
  zminus2 <- rbind(paddingzminus2 ,z[1:jminus2,,drop = FALSE])
  dzj <- zj-zminus1
  dzjminus1 <- zminus1-zminus2

  # temperature
  ti <- Y[,i, drop = FALSE]
  tij <- rbind(paddingzj, ti[1:j, ,drop = FALSE])
  timinus1 <- Y[,if(i>1){i-1}else{1}, drop = FALSE]
  tijminus1 <- rbind(paddingzminus1 ,ti[1:jminus1,,drop = FALSE])
  timinus1j <- rbind(paddingzj, timinus1[1:j, ,drop = FALSE])
  tijminus2 <- rbind(paddingzminus2 ,ti[1:jminus2,,drop = FALSE])
  tijminus3 <- rbind(paddingzminus3 ,ti[1:jminus3,,drop = FALSE])

  # Molecular coefficient of diffusion Dtam
  datmij <- mapply(baresoilevap::Datm, tij, pij )
  datmijminus1 <-  mapply(baresoilevap::Datm, tijminus1, pij )
  datmijminus2 <-  mapply(baresoilevap::Datm, tijminus2, pij )

  # specific humidity
  qi <- Z[,i, drop = FALSE]
  qij <- rbind(paddingzj,qi[1:j, ,drop = FALSE])
  qijminus1 <- rbind(paddingzminus1,qi[1:jminus1, ,drop = FALSE])
  qijminus2 <- rbind(paddingzminus2,qi[1:jminus2, ,drop = FALSE])
  qijminus3 <- rbind(paddingzminus3,qi[1:jminus3, ,drop = FALSE])

    # AkBk <- Akm1Bk + AkBkm1 - Akm1Bkm1
  Akm1Bk <- if(j>2){
    mapply(function(x,y){x*y},baresoilevap::prod_rho_datm_ftheta(i,j-1,X,Y,tim,pij,thetasat),qij-2*qijminus1+qijminus2)
  } else {0}
  Akm1Bk <- mapply(function(x,y){x/y},Akm1Bk,dzj)
  Akm1Bk <- mapply(function(x,y){x/y},Akm1Bk,dzj)
  AkBkm1 <- if(j>3){
    mapply(function(x,y){x*y},baresoilevap::prod_rho_datm_ftheta(i,j,X,Y,tim,pij,thetasat),qijminus1-2*qijminus2+qijminus3)
  } else {0}
  AkBkm1 <- mapply(function(x,y){x/y},AkBkm1,dzjminus1)
  AkBkm1 <- mapply(function(x,y){x/y},AkBkm1,dzjminus1)
  Akm1Bkm1 <- if(j>3){
    mapply(function(x,y){x*y},baresoilevap::prod_rho_datm_ftheta(i,j-1,X,Y,tim,pij,thetasat),qijminus1-2*qijminus2+qijminus3)
  } else {0}
  Akm1Bkm1 <- mapply(function(x,y){x/y},Akm1Bkm1,dzjminus1)
  Akm1Bkm1 <- mapply(function(x,y){x/y},Akm1Bkm1,dzjminus1)
  result <- Akm1Bk + AkBkm1 - Akm1Bkm1
  return(result)
}

#' Calculation of coefficient of first order partial derivative of specific humidity (q) with length (z)
#'
#' @param i
#' @param j
#' @param X
#' @param Y
#' @param Z
#' @param tim
#' @param pij
#' @param thetasat
#'
#' @return
#' @export
#'
#' @examples
coeff_dqdz <- function(i
                        , j
                        ,X=matrix(NA,nrow = 4, ncol = 4)
                        , Y = matrix(NA,nrow = 4, ncol = 4)
                        , Z = matrix(NA,nrow = 4, ncol = 4)
                        ,tim
                        , pij
                        , thetasat){
  nrow <- dim(X)[[1]]
  paddingzj <- matrix(rep(0,nrow-j), nrow = nrow-j, ncol = 1)
  jminus1 <- if(j>1){j-1}else {1}
  jminus2 <- if(j>2){j-2}else{1}
  jminus3 <- if(j>3){j-3}else{1}
  paddingzminus1 <- matrix(rep(0,nrow-jminus1), nrow = nrow-jminus1, ncol = 1)
  paddingzminus2 <- matrix(rep(0,nrow-jminus2), nrow = nrow-jminus2, ncol = 1)
  paddingzminus3 <- matrix(rep(0,nrow-jminus3), nrow = nrow-jminus3, ncol = 1)
  z <- X[,1,drop = FALSE]
  zj <- rbind(paddingzj ,z[1:j,,drop = FALSE])
  zminus1 <- rbind(paddingzminus1 ,z[1:jminus1,,drop = FALSE])
  zminus2 <- rbind(paddingzminus2 ,z[1:jminus2,,drop = FALSE])
  dzj <- zj-zminus1
  dzjminus1 <- zminus1-zminus2

  # temperature
  ti <- Y[,i, drop = FALSE]
  tij <- rbind(paddingzj, ti[1:j, ,drop = FALSE])
  timinus1 <- Y[,if(i>1){i-1}else{1}, drop = FALSE]
  tijminus1 <- rbind(paddingzminus1 ,ti[1:jminus1,,drop = FALSE])
  timinus1j <- rbind(paddingzj, timinus1[1:j, ,drop = FALSE])
  tijminus2 <- rbind(paddingzminus2 ,ti[1:jminus2,,drop = FALSE])
  tijminus3 <- rbind(paddingzminus3 ,ti[1:jminus3,,drop = FALSE])

  # Density of air
  rhoi <- matrix(apply(ti, MARGIN=2, FUN = baresoilevap::rho_air),nrow = nrow, ncol = 1)

  # specific humidity
  qi <- Z[,i, drop = FALSE]
  qij <- rbind(paddingzj,qi[1:j, ,drop = FALSE])
  qijminus1 <- rbind(paddingzminus1,qi[1:jminus1, ,drop = FALSE])
  qijminus2 <- rbind(paddingzminus2,qi[1:jminus2, ,drop = FALSE])
  qijminus3 <- rbind(paddingzminus3,qi[1:jminus3, ,drop = FALSE])

  # AkBk <- Akm1Bk + AkBkm1 - Akm1Bkm1 : linearization using Newton formula
  # firstderivative_backwardiff <- function (j=4
  #                                          , Y=matrix(rep(0,4),nrow=4, ncol=1)
  #                                          , X=matrix(rep(0,4),nrow=4, ncol=1)){
  Akm1Bk <- mapply(function(x,y)x*y, baresoilevap::prod_rho_datm(i,j-1,Y,tim,pij),baresoilevap::firstderivative_backwardiff(j-1,rhoi,z) )
  Akm1Bk <- mapply(function(x,y)x+y, Akm1Bk + baresoilevap::prod_fthetarho_dDatmdz(i,j-1,X, Y, tim, pij, thetasat))
  Akm1Bk <- mapply(function(x,y)x+y, Akm1Bk + baresoilevap::prod_rhodatm_dfthetadz(i,j-1,X,Y,tim,pij,thetasat))
  Akm1Bk <- mapply(function(x,y)x*y, Akm1Bk,baresoilevap::firstderivative_backwardiff(j,qi, zi))

  AkBkm1 <- mapply(function(x,y)x*y, baresoilevap::prod_rho_datm(i,j,Y,tim,pij),baresoilevap::firstderivative_backwardiff(j,rhoi,z) )
  AkBkm1 <- mapply(function(x,y)x+y, AkBkm1 + baresoilevap::prod_fthetarho_dDatmdz(i,j,X, Y, tim, pij, thetasat))
  AkBkm1 <- mapply(function(x,y)x+y, AkBkm1 + baresoilevap::prod_rhodatm_dfthetadz(i,j,X,Y,tim,pij,thetasat))
  AkBkm1 <- mapply(function(x,y)x*y, AkBkm1,baresoilevap::firstderivative_backwardiff(j-1,qi, zi))

  Akm1Bkm1 <- mapply(function(x,y)x*y, baresoilevap::prod_rho_datm(i,j,Y,tim,pij),baresoilevap::firstderivative_backwardiff(j,rhoi,z) )
  Akm1Bkm1 <- mapply(function(x,y)x+y, Akm1Bkm1 + baresoilevap::prod_fthetarho_dDatmdz(i,j,X, Y, tim, pij, thetasat))
  Akm1Bkm1 <- mapply(function(x,y)x+y, Akm1Bkm1 + baresoilevap::prod_rhodatm_dfthetadz(i,j,X,Y,tim,pij,thetasat))
  Akm1Bkm1 <- mapply(function(x,y)x*y, Akm1Bkm1,baresoilevap::firstderivative_backwardiff(j-1,qi, zi))

  K <- Akm1Bk +AkBkm1 - Akm1Bkm1

  return(K)



}

#' Coefficient of second order derivative of specific humidity (q) to length (z)
#'
#' @param i
#' @param j
#' @param X
#' @param Y
#' @param Z
#' @param tim
#' @param pij
#' @param thetasat
#'
#' @return
#' @export
#'
#' @examples
coeff_dq2dz2 <- function(i
                         , j
                         ,X=matrix(NA,nrow = 4, ncol = 4)
                         , Y = matrix(NA,nrow = 4, ncol = 4)
                         , Z = matrix(NA,nrow = 4, ncol = 4)
                         ,tim
                         , pij
                         , thetasat){
  nrow <- dim(X)[[1]]
  paddingzj <- matrix(rep(0,nrow-j), nrow = nrow-j, ncol = 1)
  jminus1 <- if(j>1){j-1}else {1}
  jminus2 <- if(j>2){j-2}else{1}
  jminus3 <- if(j>3){j-3}else{1}
  paddingzminus1 <- matrix(rep(0,nrow-jminus1), nrow = nrow-jminus1, ncol = 1)
  paddingzminus2 <- matrix(rep(0,nrow-jminus2), nrow = nrow-jminus2, ncol = 1)
  paddingzminus3 <- matrix(rep(0,nrow-jminus3), nrow = nrow-jminus3, ncol = 1)
  z <- X[,1,drop = FALSE]
  zj <- rbind(paddingzj ,z[1:j,,drop = FALSE])
  zminus1 <- rbind(paddingzminus1 ,z[1:jminus1,,drop = FALSE])
  zminus2 <- rbind(paddingzminus2 ,z[1:jminus2,,drop = FALSE])
  dzj <- zj-zminus1
  dzjminus1 <- zminus1-zminus2

  # temperature
  ti <- Y[,i, drop = FALSE]
  tij <- rbind(paddingzj, ti[1:j, ,drop = FALSE])
  timinus1 <- Y[,if(i>1){i-1}else{1}, drop = FALSE]
  tijminus1 <- rbind(paddingzminus1 ,ti[1:jminus1,,drop = FALSE])
  timinus1j <- rbind(paddingzj, timinus1[1:j, ,drop = FALSE])
  tijminus2 <- rbind(paddingzminus2 ,ti[1:jminus2,,drop = FALSE])
  tijminus3 <- rbind(paddingzminus3 ,ti[1:jminus3,,drop = FALSE])

  # Density of air
  rhoi <- matrix(apply(ti, MARGIN=2, FUN = baresoilevap::rho_air),nrow = nrow, ncol = 1)

  # specific humidity
  qi <- Z[,i, drop = FALSE]
  qij <- rbind(paddingzj,qi[1:j, ,drop = FALSE])
  qijminus1 <- rbind(paddingzminus1,qi[1:jminus1, ,drop = FALSE])
  qijminus2 <- rbind(paddingzminus2,qi[1:jminus2, ,drop = FALSE])
  qijminus3 <- rbind(paddingzminus3,qi[1:jminus3, ,drop = FALSE])

  # AkBk <- Akm1Bk + AkBkm1 - Akm1Bkm1 : linearization using Newton formula

  Akm1Bk <- mapply(function(x,y)x*y, baresoilevap::prod_rho_datm_ftheta(i,j-1,Y,tim,pij, thetasat),baresoilevap::secondderivative_backwarddiff(j,rhoi,z) )

  AkBkm1 <- mapply(function(x,y)x*y, baresoilevap::prod_rho_datm_ftheta(i,j,Y,tim,pij, thetasat),baresoilevap::secondderivative_backwarddiff(j-1,rhoi,z) )

  Akm1Bkm1 <- mapply(function(x,y)x*y, baresoilevap::prod_rho_datm_ftheta(i,j-1,Y,tim,pij, thetasat),baresoilevap::secondderivative_backwarddiff(j-1,rhoi,z) )

  K <- Akm1Bk + AkBkm1 - Akm1Bkm1

  return (K)

}

dQthetvap_dz <- function (i
                          , j
                          ,X=matrix(NA,nrow = 4, ncol = 4)
                          , Y = matrix(NA,nrow = 4, ncol = 4)
                          , Z = matrix(NA,nrow = 4, ncol = 4)
                          ,tim
                          , pij
                          , thetasat){

  dQthetadz <- -(baresoilevap::coeff_dq2dz2(i,j,X,Y,Z,tim,pij,thetasat))+(baresoilevap::coeff_dqdz(i,j,X,Y,Z,tim,pij, thetasat))

  return(dQthetadz)

}
