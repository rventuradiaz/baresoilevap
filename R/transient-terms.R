#' Title
#' Derivative od mass water vapor per unit volume of soil
#'
#' @param i : time node number
#' @param X : mesh for theta: row) z ; column) time ; value: volumetric water content (theta)
#' @param Y :mesh for temperature: row) z ; column) time ; value: temperature (t)
#' @param Z :mesh for specific humidity: row) z ; column) time ; value: specific humidity (q)
#' @param tim : time nodes
#' @param j : spatial node number
#' @param pij : pressure (in KPa)
#'
#' @return
#' @export
#'
#' @examples
dmasswatervaporperunitvolumesoil_dt <- function(i,j, X=matrix(NA,nrow = 4, ncol = 4), Y=matrix(NA,nrow = 4, ncol = 4), Z = matrix(NA,nrow = 4, ncol = 4),tim, pij){

  nrow <- dim(X)[[1]]
  paddingzj <- matrix(rep(0,nrow-j), nrow = nrow-j, ncol = 1)
  jminus1 <- j-1
  paddingzminus1 <- matrix(rep(0,nrow-jminus1), nrow = nrow-jminus1, ncol = 1)

  # volumetric water content (theta)
  thetasat <- 0.1
  thetai <- X[,i, drop = FALSE]
  thetaij <- rbind(paddingzj, thetai[1:j, ,drop = FALSE])
  thetaiminus1 <- X[,i-1, drop = FALSE]
  thetaijminus1 <- rbind(paddingzminus1 ,thetai[1:jminus1,,drop = FALSE])
  thetaiminus1j <- rbind(paddingzj, thetaiminus1[1:j, ,drop = FALSE])
  diftheta <- thetaij-thetaijminus1
  if (i >= 3) {
    thetaiminus2 <- X[,i-2, drop = FALSE]
    thetaiminus2j <- rbind(paddingzj, thetaiminus2[1:j, ,drop = FALSE])
    thetaijminus2 <- rbind(paddingzminus1 ,thetai[1:jminus1,,drop = FALSE])
  }

  # temperature
  ti <- Y[,i, drop = FALSE]
  tij <- rbind(paddingzj, ti[1:j, ,drop = FALSE])
  timinus1 <- Y[,i-1, drop = FALSE]
  timinus1j <- rbind(paddingzj, timinus1[1:j, ,drop = FALSE])

  if (i >= 3) {
    timinus2 <- Y[,i-2, drop = FALSE]
    timinus2j <- rbind(paddingzj, timinus2[1:j, ,drop = FALSE])
  }

  # specific humidity
  qi <- Z[,i, drop = FALSE]
  qiminus1 <- Z[,i-1, drop = FALSE]
  qij <- rbind(paddingzj, qi[1:j, ,drop = FALSE])
  if (i >= 3) {
    qiminus2 <- Z[,i-2, drop = FALSE]
    qiminus2j <- rbind(paddingzj, qiminus2[1:j, ,drop = FALSE])
  }
  qiminus1j <- rbind(paddingzj, qiminus1[1:j, ,drop = FALSE])
  difq <- qij - qiminus1j

  # air density
  rhoairij <- mapply(baresoilevap::rho_air, tij, pij )
  rhoairiminus1j <- mapply(baresoilevap::rho_air, timinus1j, pij )
  if (i >= 3) {rhoairiminus2j <- mapply(baresoilevap::rho_air, timinus2j, pij )}
  difrhoair <- rhoairij - rhoairiminus1j

  # calculate mass water vapor per volume unit of soil i,j
  thetasatv <- matrix(data = rep(thetasat,nrow), nrow = nrow, ncol = 1)
  masswatvapunitvolsoilij <-  thetasatv - thetaij
  masswatvapunitvolsoiliminus1j <- thetasatv - thetaiminus1j
  if(i >= 3) {masswatvapunitvolsoiliminus2j <- thetasatv - thetaiminus2j}


  # Calculate time derivatives
  dtimi <- tim[i]-tim[i-1]
  dtimiminus1 <- tim[i-1]-tim[i-2]

    # volumetric content "theta" time derivative
    dthetadtim <- diftheta/dtimi


  if (i <3){

    dmassvapvolliqdt1 <- mapply(function(x,y){x/y}, (thetaiminus1j-thetaij), dtimi)
    dmassvapvolliqdt1 <- mapply(function(x,y){x*y},dmassvapvolliqdt1,rhoairiminus1j)
    dmassvapvolliqdt1 <- mapply(function(x,y){x*y},dmassvapvolliqdt1,qiminus1j)
    dmassvapvolliqdt2 <- mapply(function(x,y){x*y},masswatvapunitvolsoiliminus1j,qiminus1j)
    dmassvapvolliqdt2 <- mapply(function(x,y){x*y},dmassvapvolliqdt2,difrhoair/dtimi)
    dmassvapvolliqdt3 <- mapply(function(x,y){x*y},masswatvapunitvolsoiliminus1j,rhoairiminus1j)
    dmassvapvolliqdt3 <- mapply(function(x,y){x*y},masswatvapunitvolsoiliminus1j,difq/dtimi)
    dmassvapvolliqdt <- dmassvapvolliqdt1 + dmassvapvolliqdt2 + dmassvapvolliqdt3
    result <- dmassvapvolliqdt
  } else {
    dmassvapvolliqdt1 <- mapply(function(x,y){x/y}, (thetaiminus1j-thetaij), dtimi)
    dmassvapvolliqdt1 <- mapply(function(x,y){x*y},dmassvapvolliqdt1,rhoairiminus1j)
    dmassvapvolliqdt1 <- mapply(function(x,y){x*y},dmassvapvolliqdt1,qiminus1j)
    dmassvapvolliqdt2 <- mapply(function(x,y){x/y}, (thetaiminus2j-thetaiminus1j), dtimiminus1)
    dmassvapvolliqdt2 <- mapply(function(x,y){x*y},dmassvapvolliqdt2,rhoairij)
    dmassvapvolliqdt2 <- mapply(function(x,y){x*y},dmassvapvolliqdt2,qij)
    dmassvapvolliqdt3 <- mapply(function(x,y){x/y}, (thetaiminus2j-thetaiminus1j), dtimiminus1)
    dmassvapvolliqdt3 <- mapply(function(x,y){x*y},dmassvapvolliqdt3,rhoairiminus1j)
    dmassvapvolliqdt3 <- mapply(function(x,y){x*y},dmassvapvolliqdt3,qiminus1j)
    dmassvapvolliqdt <- dmassvapvolliqdt1+dmassvapvolliqdt2-dmassvapvolliqdt3
    dmassvapvolliqdt1 <- mapply(function(x,y){x/y}, (rhoairij-rhoairiminus1j), dtimi)
    dmassvapvolliqdt1 <- mapply(function(x,y){x*y},dmassvapvolliqdt1,qiminus1j)
    dmassvapvolliqdt1 <- mapply(function(x,y){x*y},dmassvapvolliqdt1,(thetasatv-thetaiminus1j))
    dmassvapvolliqdt2 <- mapply(function(x,y){x/y}, (rhoairiminus1j-rhoairiminus2j), dtimiminus1)
    dmassvapvolliqdt2 <- mapply(function(x,y){x*y},dmassvapvolliqdt2,(thetasatv-thetaij))
    dmassvapvolliqdt2 <- mapply(function(x,y){x*y},dmassvapvolliqdt2,qij)
    dmassvapvolliqdt3 <- mapply(function(x,y){x/y}, (rhoairiminus1j-rhoairiminus2j), dtimiminus1)
    dmassvapvolliqdt3 <- mapply(function(x,y){x*y},dmassvapvolliqdt3,(thetasatv-thetaiminus1j))
    dmassvapvolliqdt3 <- mapply(function(x,y){x*y},dmassvapvolliqdt3,qiminus1j)
    dmassvapvolliqdt <- dmassvapvolliqdt + dmassvapvolliqdt1 + dmassvapvolliqdt2 - dmassvapvolliqdt3
    dmassvapvolliqdt1 <- mapply(function(x,y){x/y}, (qij-qiminus1j), dtimi)
    dmassvapvolliqdt1 <- mapply(function(x,y){x*y},dmassvapvolliqdt1,rhoairiminus1j)
    dmassvapvolliqdt1 <- mapply(function(x,y){x*y},dmassvapvolliqdt1,(thetasatv-thetaiminus1j))
    dmassvapvolliqdt2 <- mapply(function(x,y){x/y}, (qiminus1j-qiminus2j), dtimiminus1)
    dmassvapvolliqdt2 <- mapply(function(x,y){x*y},dmassvapvolliqdt2,(thetasatv-thetaij))
    dmassvapvolliqdt2 <- mapply(function(x,y){x*y},dmassvapvolliqdt2,rhoairij)
    dmassvapvolliqdt3 <- mapply(function(x,y){x/y}, (qiminus1j-qiminus2j), dtimiminus1)
    dmassvapvolliqdt3 <- mapply(function(x,y){x*y},dmassvapvolliqdt3,(thetasatv-thetaiminus1j))
    dmassvapvolliqdt3 <- mapply(function(x,y){x*y},dmassvapvolliqdt3,rhoairiminus1j)
    dmassvapvolliqdt <- dmassvapvolliqdt + dmassvapvolliqdt1 + dmassvapvolliqdt2 - dmassvapvolliqdt3
    result <- dmassvapvolliqdt
  }
return(result)
}
