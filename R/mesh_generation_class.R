# Generate mesh 1D problem

# class "cfd_mesh", storage of mesh
## matrix with 1:3 coordinates plus boundary condition (Dirichlet or Neumann)
## Methods: 1) dimensioning ; 2) feeding boundary conditions at initial time
## inparam is a list with dimensionality di (1D,2D, or 3D), nodes, node distribution function in space xyf_dist, node distribution in time t_dist
# sigma <- 0.5
# l <- c(0.3,0.0,0.0)
# nxyz <- c(40,0,0)
# len_nod <- rbind(l,nxyz)
# di <- 1

xyz_dist <- function(len_nod=matrix(data = c(c(rep(0.1,3)), c(rep(10,3))),nrow =2,ncol=3, byrow = TRUE) ){
  len <- len_nod[1,]
  innod <- len_nod[2,]
  rtrn <- vector(mode="numeric", length = ncol(len_nod))
  for (i in 1:length(rtrn)){
    rtrn[i]<-len[i]/(innod[i]-1)
  }
  rtrn
}
# Test
# lnod <- matrix(c(c(0.3,0,0), c(40,0,0)), nrow = 2, ncol = 3)
# xyz_dist(len_nod = lnod)

#' Title
#'
#' @param xyz_d
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
t_dist <- function(xyz_d=vector(mode = "numeric",length = 3), sigma=0.5){
  rtrn <- vector(mode = "numeric", length = 1)
  for (i in 1:length(rtrn)){
    rtrn[i] <- sigma * max(xyz_d)
  }
  return(rtrn)
}

# Test
# t_dist(xyz_dist(len_nod = lnod), 0.5)

inparam <- list(di=1
                ,len_nod=matrix(data = c(c(rep(0.1,3)), c(rep(10,3))),nrow =2,ncol=3, byrow = TRUE)
                ,geom_dist = function(len_nod=matrix(data=0.0, nrow= 2,ncol=3) ){
                  len <- len_nod[1,]
                  innod <- len_nod[2,]
                  rtrn <- vector(mode="numeric", length = 3)
                  for (i in 1:length(rtrn)){
                    rtrn[i]<-len[i]/(innod[i]-1)
                  }
                  rtrn
                }
                ,time_dist = function(xyz_d=vector(mode = "numeric",length = 3), sigma=0.5){
                  rtrn <- vector(mode = "numeric", length = 1)
                  for (i in 1:length(rtrn)){
                    rtrn[i] <- sigma * max(xyz_d)
                  }
                  rtrn
                }
                ,sigma= 0.5)

sum1toi <- function(i,j) return(sum(i)*(j+1)) #i is mesh nodes vector, j is dimension

mesh <- function(inparam=list(di=1
                              ,len_nod=matrix(data = c(c(rep(0.1,3)), c(rep(10,3))),nrow =2,ncol=3, byrow = TRUE)
                              ,geom_dist = function(len_nod=matrix(data=0.0, nrow= 2,ncol=3) ){
                                len <- len_nod[1,]
                                innod <- len_nod[2,]
                                rtrn <- vector(mode="numeric", length = 3)
                                for (i in 1:length(rtrn)){
                                  rtrn[i]<-len[i]/(innod[i]-1)
                                }
                                rtrn
                              }
                              ,time_dist = function(xyz_d=vector(mode = "numeric",length = 3), sigma=0.5){
                                rtrn <- vector(mode = "numeric", length = 1)
                                for (i in 1:length(rtrn)){
                                  rtrn[i] <- sigma * max(xyz_d)
                                }
                                rtrn
                              }
                              ,sigma= 0.5)){
  # debug(inparam[[4]])
  rtrn <- list()
  class(rtrn) <- "cdf_mesh"
  spa_di <- inparam[[1]]
  lnods <- (inparam[[2]])
  n <- sum1toi(max(lnods),spa_di)
  rtrn$mat <- vector(length= n)
  rtrn$ix <- (1:n)
  rtrn$dimensions <- spa_di
  rtrn$nodes <- lnods[2,]                               #! number of nodes
  rtrn$n <- sum1toi(max(lnods[2,]),spa_di)
  fun_spa_dist <- inparam[[3]]
  fun_t_dist <- inparam[[4]]
  sigma <- inparam[[5]]
  t_i <- 0
  for (i in seq(1,n,spa_di+1)){
    t_i <- t_i + fun_t_dist(fun_spa_dist(lnods),sigma)
    rtrn$mat[i]<- t_i
  }
  xyz_step <- 0
  for (j in 1:spa_di){
    for (i in seq(j+1,n,spa_di+1)){
      xyz_step <- xyz_step +fun_spa_dist(lnods)[j]
      rtrn$mat[i]<- xyz_step
    }
    # The mesh does not have boundary conditions
    # for (i in seq(spa_di+2,n,spa_di+2)){
    #   xyz_step <- xyz_step +fun_spa_dist(lnods)[j]
    #   theta <- 0.1/xyz_step
    #   rtrn$mat[i]<- fun_bc(theta=theta)
    # }

  }
  rtrn
}


print.mesh <- function(mesh){
  m_n <- length(mesh$mat)
  geomdim <- mesh$dimension
  n_col <- mesh$dimension+1
  nodes <- max(mesh$nodes)
  mat_vector <- mesh$mat
  n <- mesh$n
  full_mesh <- matrix(nrow=nodes, ncol = n_col)
  for (i in 1:(geomdim+1)) {
    full_mesh[1:nodes,i] <- mat_vector[seq(i,m_n,n_col)]
  }
  return(full_mesh)
}

