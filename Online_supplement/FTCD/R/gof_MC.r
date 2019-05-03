#File: gof_MC_Gamma.r
#Article: Goodness-of-fit tests for compound distributions with applications in insurance
#Package: testCD
#Author: Pierre-O Goffard
#Function that perform all the tests and also produce the desired outputs for the Monte Carlo experiment



#' p_value
#' Compute the p-value of a Compound distribution goodness-of-fit test using
#' the Warp speed method
#' @param n Size of the sample
#' @param N_dist_alt Claim frequency distribution of the input sample
#' @param U_dist_alt Claim sizes distribution of the input sample
#' @param N_dist_H0 Claim frequency distribution under H0
#' @param U_dist_H0 Claim sizes distribution under H0
#' @param MC Number of Monte Carlo iteration
#' @param method_infer Type of MM estimation c('partial-MME', 'MME', )
#' @param method_gof Type of GOF test c('KS', 'CvM' 'LT_ODE', 'LT_L2')
#' @param par_method Parameter of the method discretization width for the distribution
#' function gof test and parameter of the weight function for the Laplace transform based
#' test c(h, beta)
#'
#' @return the p value of the test
#' @export
#'
#' @examples
p_value <- function(n, N_dist_alt, U_dist_alt, par_alt_dist, N_dist_H0, U_dist_H0, MC,
                    method_infer, method_gof, par_method){

  TS <- matrix(0, MC, 2)

  if(method_gof == 'KS'){

    for(k in 1:MC){

      X_data <- simulate_X(n, N_dist_alt, U_dist_alt, par_alt_dist)
      X_par <- MME_CD(X_data, N_dist_H0, U_dist_H0, method_infer)

      while((X_par[1]<0|X_par[1]>1) | is.na(X_par[1]) | X_par[3]<0 |
            is.na(X_par[3]) | X_par[4]<0 | is.na(X_par[4])){

        X_data <- simulate_X(n, N_dist_alt, U_dist_alt, par_alt_dist)
        X_par <- MME_CD(X_data, N_dist_H0, U_dist_H0, method_infer)

      }



      TS[k, ] <- ks_gof(X_data , N_dist_H0, U_dist_H0, par_method,
                        'truncation', 1, method_infer)
    }


    }else if(method_gof == 'CvM'){

      for(k in 1:MC){

        X_data <- simulate_X(n, N_dist_alt, U_dist_alt, par_alt_dist)
        X_par <- MME_CD(X_data, N_dist_H0, U_dist_H0, method_infer)

        while((X_par[1]<0|X_par[1]>1) | is.na(X_par[1]) | X_par[3]<0 |
              is.na(X_par[3]) | X_par[4]<0 | is.na(X_par[4])){

          X_data <- simulate_X(n, N_dist_alt, U_dist_alt, par_alt_dist)
          X_par <- MME_CD(X_data, N_dist_H0, U_dist_H0, method_infer)

        }


        TS[k, ] <- cvm_gof(X_data , N_dist_H0, U_dist_H0, par_method,
                          'truncation', 1, method_infer)
      }


    }else if(method_gof == 'LT_L2'){

      for(k in 1:MC){

        X_data <- simulate_X(n, N_dist_alt, U_dist_alt, par_alt_dist)
        X_par <- MME_CD(X_data, N_dist_H0, U_dist_H0, method_infer)

        while((X_par[1]<0|X_par[1]>1) | is.na(X_par[1]) | X_par[3]<0 |
              is.na(X_par[3]) | X_par[4]<0 | is.na(X_par[4])){

          X_data <- simulate_X(n, N_dist_alt, U_dist_alt, par_alt_dist)
          X_par <- MME_CD(X_data, N_dist_H0, U_dist_H0, method_infer)

        }


        TS[k, ] <- lt_L2_gof(X_data, N_dist_H0,
                             U_dist_H0, par_method, 1, method_infer, TRUE)
      }


    }else if(method_gof == 'LT_ODE'){
      for(k in 1:MC){

        X_data <- simulate_X(n, N_dist_alt, U_dist_alt, par_alt_dist)
        X_par <- MME_CD(X_data, N_dist_H0, U_dist_H0, method_infer)

        while((X_par[1]<0|X_par[1]>1) | is.na(X_par[1]) | X_par[3]<0 |
              is.na(X_par[3]) | X_par[4]<0 | is.na(X_par[4])){

          X_data <- simulate_X(n, N_dist_alt, U_dist_alt, par_alt_dist)
          X_par <- MME_CD(X_data, N_dist_H0, U_dist_H0, method_infer)

        }


        TS[k, ] <- lt_ODE_gof(X_data, N_dist_H0,
                           U_dist_H0, par_method, 1, method_infer, TRUE)
      }
    }
  return(mean(TS[is.na(TS[ , 1])==FALSE ,1] > quantile(TS[, 2], 0.95, na.rm = TRUE)))
#  return(TS)

}
