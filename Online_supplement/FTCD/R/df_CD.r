#File: df_CD.r
#Article: Goodness-of-fit tests for compound distributions with applications in insurance
#Package: testCD
#Author: Pierre-O Goffard
#Approximation of the compound distribution function via Panjer algorithm


#' df_X
#' Distribution function of the aggregated losses approximated via Panjer algorithm
#' @param X_par Parameter of the compound distribution
#' @param N_dist Claim frequency distribution c('pois', 'geom')
#' @param U_dist Claim sizes distribution c('exp', 'gamma', 'weibull', 'invgauss')
#' @param max_X Truncation above the threshold max_X
#' @param par_method_proxy Width of the discretization of the claim sizes distribution
#' when Panjer algorithm is used, order of truncation when the truncation method is used
#' @param method_proxy Aproximation method for the compound distribution DF c('recursive', 'truncation')
#'
#' @return The distribution function of X as a function
#' @export
#'
#' @examples
df_X <- function(X_par, N_dist, U_dist, max_X = NULL, par_method_proxy, method_proxy){

  if(method_proxy=='truncation'){
    if(N_dist == 'pois'){
      if(U_dist=='exp'){

      F_X <- function(x){
        dpois(0, X_par[2]) +
          sum(dpois(1:par_method_proxy, X_par[2]) *
                pgamma(x, shape = 1:par_method_proxy, scale = X_par[4]))}

      }else if(U_dist=='gamma'){

        F_X <- function(x){
          dpois(0, X_par[2]) +
            sum(dpois(1:par_method_proxy, X_par[2]) *
                  pgamma(x, shape = (1:par_method_proxy) * X_par[3], scale = X_par[4]))}

      }else if(U_dist=='invgauss'){

        F_X <- function(x){
          dpois(0, X_par[2]) +
            sum(dpois(1:par_method_proxy, X_par[2]) *
                  actuar::pinvgauss(x, mean = (1:par_method_proxy) * X_par[3],
                                    dispersion = X_par[4] / (1:par_method_proxy)^2))}

      }
    }else if(N_dist == 'geom'){
      if(U_dist=='exp'){
        F_X <- function(x){
          dgeom(0, X_par[1]) +
            sum(dgeom(1:par_method_proxy, X_par[1]) *
                  pgamma(x, shape = 1:par_method_proxy, scale = X_par[4]))}

      }else if(U_dist=='gamma'){

        F_X <- function(x){
          dgeom(0, X_par[1]) +
            sum(dgeom(1:par_method_proxy, X_par[1]) *
                  pgamma(x, shape = (1:50) * X_par[3], scale = X_par[4]))}

      }else if(U_dist=='invgauss'){

        F_X <- function(x){
          dgeom(0, X_par[1]) +
            sum(dgeom(1:par_method_proxy, X_par[1]) *
                  actuar::pinvgauss(x, mean = (1:par_method_proxy) * X_par[3],
                                    dispersion = X_par[4] / (1:par_method_proxy)^2))}

      }
    }
  }else if(method_proxy=='recursive'){
    if(N_dist == "geom"){

      if(U_dist == "exp"){

        f_U <- actuar::discretize(pexp(x, 1 / X_par[4]), method = "lower", from = 0, to = max_X, step = h)


      }else if(U_dist == "gamma"){

        f_U <- actuar::discretize(pgamma(x, shape = X_par[3], scale = X_par[4]),
                                  method = "rounding", from = 0, to = max_X, step = h)

      }else if(U_dist == "weibull"){

        f_U <- actuar::discretize(pweibull(x, shape = X_par[3], scale = X_par[4]),
                                  method = "rounding", from = 0, to = max_X, step = h)
      }else if(U_dist == "invgauss"){

        f_U <- actuar::discretize(actuar::pinvgauss(x, mean = X_par[3], dis = X_par[4]),
                                  method = "rounding", from = 0, to = max_X, step = h)
      }

      F_X <- actuar::aggregateDist(method = "recursive", model.freq = "geometric", x.scale = h,
                                   model.sev = f_U, prob = 1 - X_par[1])
    }else if(N_dist == "pois"){

      if(U_dist == "exp"){

        f_U <- actuar::discretize(pexp(x, 1 / X_par[4]), method = "lower", from = 0, to = max_X, step = h)

      }else if(U_dist == "gamma"){

        f_U <- actuar::discretize(pgamma(x, shape = X_par[3], scale = X_par[4]),
                                  method = "rounding", from = 0, to = max_X, step = h)

      }else if(U_dist == "weibull"){

        f_U <- actuar::discretize(pweibull(x, shape = X_par[3], scale = X_par[4]),
                                  method = "rounding", from = 0, to = max_X, step = h)

      }else if(U_dist == "invgauss"){

        f_U <- actuar::discretize(actuar::pinvgauss(x, mean = X_par[3], dis = X_par[4]),
                                  method = "rounding", from = 0, to = max_X, step = h)

      }

      F_X <- actuar::aggregateDist(method = "recursive", model.freq = "poisson", x.scale = h,
                                   model.sev = f_U, lambda = X_par[2])
  }

  }
  return(F_X)

}
