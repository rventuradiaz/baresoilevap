## Preprocessing module



# require(CHNOSZ)
# require(marelac)
# require(humidity)
# Parameters
theta_sat_ <- 0.40 # saturated volumetric soil water content
c_theta_ <- 0.12
psi_a1_ <- -4.4*10^4
psi_a2_ <- -15*10^4
psi_sat_ <- -0.05
a1_ <- 54
a2_ <- 900
b_ <- 6.0
c_ <- 15
C_soil <- 1.26e06 #[Jm^-03K^-01]
C_water <- 4.20e06 #[Jm^-03K^-01]
K_sat_ <- 3.5*10^(-5)

# Albedo
a <- 0.3 # albedo


# Air calorific capacity [J*Kg^-1*K^-1]
c_p = 1.007e03

#solar constant [W/m**2]
s_c  <- 750

# soil emissivipty
epsilon <- 0.1 #check narita sand's

# effective sky temperature [K]
 t_sky <- 285 #warm cloudy conditios / 230 cold, clear sky

# long wave radiation
sigma_q <- 5.67e-08 # W*m**-2*K**-4

# Water density
phrow <- 1e03 #[Kg/m^3]


# Soil surface temperature [K]
T_s <- 30+273.15 #[K]

# Soil temperature [K]
t <- 30+273.15 #[K]

# Bulk transfer coefficient for water vapor
C_E <-3.0e-03
C_H <- C_E

# Wind velocity
u <- 1 # [m/s]

# Atmospheric pressure
P_atm <- 101325.0

# Initial volumetric soil water content
theta <- 0.10

# Soil water potential (psi)
#' Title
#'
#' @param x:  theta
#'
#' @return
#' @export
#'
#' @examples
psi <- function(x){
  -0.05*(x/0.4)**(-6)
}

# Relation among soil water potential and temperature d(psi)/d(T)
dpsidt <- function(x){
  gamma <- -1.48e-05*x[2]+2.26e-03
  y <- psi(x[1])*gamma
  y
}


# Resitance of water vapor transport
#' Title
#'
#' @param x: vector of variables tetha (x[1]), T (x[2]), and q (x[3])
#'
#' @return
#' @export
#'
#' @examples
Fprima <- function(x){
  0.04*exp(-200.0*x[1]^2)+0.003*exp(-100.0*x[1]^2)
}

#' Title
#'
#' @param theta
#'
#' @return
#' @export
#'
#' @examples
Rtheta <- function(theta){
  0.02*baresoilevap::Fprima(theta)
}
# Thermodynamic properties

# Relative humidity (h)
#' Title
#'
#' @param x theta
#' @param y: temperature in K
#'
#' @return z
#' @export
#'
#' @examples
  relhumidity <- function(x,y) {
    # Acceleration of gravity (g)
    g <- 9.80665 # [m2]*[s^(-2)]

    # R constant for water vapor (RW) [J/kg K]
    Rw <- 461.52

    z <- exp(baresoilevap::psi(x)*g/(Rw*y))
    return(z)
  }


# Saturation specific humidity at the soil surface temperature
#' Title
#'
#' @param x: T
#' @param P
#'
#' @return
#' @export
#'
#' @examples
  q_Ts <- function(x, P=101300){
    # require(humidity)
    e <- humidity::WVP1(x)
    y <- humidity::SH(e, P) #[kg]/[Kg]
    return(y)
  }

# Thermal conductivity (lambda)
#' Title
#'
#' @param x: vector of variables theta (x[1]), T (x[2]), and q (x[3])
#'
#' @return
#' @export
#'
#' @examples
  lambda <- function (x){
    y <- 0.251+0.5*x[1]^(1/3)
    return(y)
  }

# Volumetric heat capacity (Cvh)
#' Title
#'
#' @param theta_sat
#' @param C_soil
#' @param C_water
#' @param x: vector of variables theta (x[1]), T (x[2]), and q (x[3])
#'
#' @return
#' @export
#'
#' @examples
  Cvh <- function(theta_sat=0.1, C_soil = 1, C_water = 0.1, x){
    (1-theta_sat)*C_soil + x[1]*C_water
  }

# Specific humidity of air (qs) in Pa [kg*m^-1*s^-2]
#' Title
#'
#' @param x: T
#'
#' @return
#' @export
#'
#' @examples
  wvp_antoine<- function(x){
    A<-18.3036
    B<-3816.44
    C<-(-46.13)
    p <- exp(A-B/(C+x)) # In mmHg
    p <- p * 133.322387415 #in Pa
    return(p)
  }
  # e <- wvp_antoine(t)
  # qs <- SH(e, p=P_atm)
  # qs

# Air density : t_air in Kelvin, p in Pa
#' Title
#'
#' @param t_air
#' @param p
#'
#' @return
#' @export
#'
#' @examples
rho_air <- function(t_air = 298.73,p = 101300 ){
  # require(marelac)
  t <- t_air - 273.15 # t_air in degrees celsius
  p <- p /100000.0 # p in bar
  marelac::air_density(t, p)
}


# Water Density
#' Title
#'
#' @param x: vector of variables theta (x[1]), T (x[2]), and q (x[3])
#'
#' @return
#' @export
#'
#' @examples
rho_w <- function(x){
  # require(CHNOSZ)
  CHNOSZ::data(thermo) #Initialize package with thermodynamics data
  y <- water(property = "rho", T = x[2]) # P = Psat
  return(y)
}

# Hydraulic Conductivity (K (theta))
#' @param x
K_theta <- function(x) {
  y <- 3.5e-05*(x/0.4)**15
  return(y)
}


# Heat flux Qh boundary condition
#' Title
#'
#' @param j node index
#' @param a albedo
#' @param sc
#' @param L
#' @param sigma_q
#' @param x: vector of variables theta (x[1]), T (x[2]), and q (x[3])
#' @param cp
#' @param t_air
#' @param p
#' @param c_h
#' @param u
#'
#' @return
#' @export
#'
#' @examples
Qh <- function(j=0, a=0.3,sc = 750, L=0, sigma_q=5e-08, x, cp = 1.007, t_air=298.73, p=101300, c_h = 10e-03, u=1) {
  q_h <- if(j==0){
    -((1-a)*sc+L-sigma_q*x[2]^4-cp*rho_air(t_air , p)*c_h*u*(x[2]-t_air))}
  else{0}
  return(q_h)
}

# Upward water flux boundary condition
#' Title
#'
#' @param j node index
#' @param theta
#'
#' @return
#' @export
#'
#' @examples
Q_theta_liq <- function (j=0,theta = 0.1){
  q_theta_liq <- if(j==0){
    0
  }else{
    baresoilevap::K_theta(theta)
  }
  return(q_theta_liq)
}

# Upward vapor flux boundary condition
#' Title
#'
#' @param j node index
#' @param t_air air temperature in K
#' @param p_air pressure in Pa
#' @param c_e calorific capacity of air
#' @param u wind speed
#' @param x: Volumetric content of water theta
#' @param y: Temperature
#' @param p
#'
#' @return
#' @export
#'
#' @examples
Q_theta_vap <- function(x, y, j=0, t_air = 298.73, c_e = 10e-03, u = 1, p=101300){
  if(j==0){
    e <- wvp_antoine(t_air)
    qs <- humidity::SH(e, p)
    qthetavap <- baresoilevap::rho_air(t_air, p )*c_e*u*(baresoilevap::relhumidity(x, y)*q_Ts(y,p)-qs)
  }else {
    0
  }
  return(qthetavap)
}

# Molecular diffusion coefficient for water vapor
#' Title
#'
#' @param t_air
#' @param p
#'
#' @return
#' @export
#'
#' @examples
Datm <- function(t_air, p){
  # t in K, p in atm
  # coefficients for pair consisting of H2O and a nonpolar gas
  # Dab in m^2*s^(-01)
  a<- 3.640e-04
  b<- 2.334
  # pressure at critical point (pc)
  pca <-  218.3 # a: water vapor
  pcb <- 37.2 # b: air
  # Temperatur at critical point (tc)
  tca <- 647.4
  tcb <- 132.5
  # molecular mass
  Ma <- 18.06
  Mb <- 28.97
  Dab <- a*( t_air / (sqrt(tca*tcb)))^2*(pca*pcb)^(1/3)*(tca*tcb)^(5/12)*sqrt(1/Ma+1/Mb)/p
  return (Dab * 1e-02) #' In m^2*s^-1
}

#' Porosity and tortuosity factor f(theta)
#'
#' @param x
#' @param thetasat
#'
#' @return
#' @export
#'
#' @examples
porosity_factor <- function ( x, thetasat) {
  porosityfactor <- 0.66*(thetasat-x)
  return(porosityfactor)
}

# source("mesh_generation_class.R")

# mat_len_nodes <- matrix(c(c(0.3),c(40)), nrow = 2, ncol = 1, byrow = TRUE)
# inparam <- list(di = 1
#                 , len_nod = mat_len_nodes
#                 , geom_dist = xyz_dist
#                 , time_dist = t_dist
#                 ,sigma = 0.5
# )
#
# mesh_c <- mesh(inparam)
# mesh1 <- print.mesh(mesh_c)

# Boundary condition
#' Title
#'
#' @param mesh
#' @param nbc
#'
#' @return
#' @export
#'
#' @examples
bc <- function(mesh, nbc){
  require(humidity)
  nodes <- mesh_c$nodes
  ix <- sum1toi(nodes,nbc-1)
  bc <- vector(length = ix, mode = "numeric")
  mat_bc <- matrix(bc, nrow = nodes, ncol = nbc, byrow = TRUE)
  bc <- c(0,Q_theta_vap(),Qh())
  bc <- append(bc, c(rep(0,ix-6)))
  bc <- append(bc, c(Q_theta_liq(),0, 0))
  for (i in (1:nbc)){
    mat_bc[1:nodes,i] <- bc[seq(i,ix,nbc)]
  }
  return(mat_bc)
}

#' Title
#'
#' @param x temperature (in Kelvin degrees)
#'
#' @return
#' @export
#'
#' @examples
gammat <- function(x){
  y <- -1.48e-05*x+2.26e-03
  return(y)

}

#' Title
#'
#' @param x : theta
#' @param y : temperature (in K)
#'
#' @return z: Diffusion coefficient
#' @export
#'
#' @examples
Dtliq <- function(x, y){
  z <- baresoilevap::K_theta(x)
  z <- z*baresoilevap::psi(x)
  z <- z*baresoilevap::gammat(y)
  return(z)
}


#' Title
#'
#' @param x
#'
#' @return d: diffusion coefficient in liquid
#' @export
#'
#' @examples
Dthetaliq <- function (x)  {
  d <- 0.0400543212890625*x^8
  return(d)
}

# Rate of vaporisation from the water in small pores in the soil

Esoil <- function(t=298.0 # in Kelving degrees
                  , theta # dimensionless
                  , q  # in
) {
  esoil <- baresoilevap::rho_air(t)
  esoil <- esoil * baresoilevap::Datm(t,p=P_atm)
  esoil <- esoil * (baresoilevap::relhumidity(theta,t)*baresoilevap::q_Ts(t)-q)
  esoil <- esoil / baresoilevap::Rtheta(theta)
  return(esoil)
}
