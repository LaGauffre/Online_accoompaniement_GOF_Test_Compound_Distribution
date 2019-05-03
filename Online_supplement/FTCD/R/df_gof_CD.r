#File: df_gof_CD.r
#Article: Goodness-of-fit tests for compound distributions with applications in insurance
#Package: testCD
#Author: Pierre-O Goffard
#DF based goodness-of-fit-test for compound distributions

#' ks_gof
#' Perform a goodness-of-fit test based on the Kolmogorov-Smirnov distance
#' @param X_data Sample of aggregated losses
#' @param N_dist Claim frequency distribution
#' @param U_dist Claim sizes distribution
#' @param par_method_proxy Width of the discretization of the claim sizes distribution
#' when Panjer algorithm is used, order of truncation when the truncation method is used
#' @param method_proxy_df Approximation method for the compound distribution DF c('recursive', truncation)
#' @param B Number of bootstrap loop
#' @param method Type of MM estimation c('MME', 'partial-MME')
#'
#' @return The value of the test statistic and the critical value of the test
#' @export
#'
#' @examples
ks_gof <- function(X_data, N_dist, U_dist, par_method_proxy,
                   method_proxy_df, B, method = 'partial-MME'){

  # Inference of the parameters of the compound distribution
  ## If no zeros use the full-MME method otherwise always the other one!
  if(sum(X_data==0)==0){

    X_par <- MME_CD(X_data, N_dist, U_dist, 'MME')
    # Case when the data generates estimates that are negative
    while((X_par[1]<0|X_par[1]>1) | is.na(X_par[1]) | X_par[3]<0 | is.na(X_par[3]) | X_par[4]<0 | is.na(X_par[4])){

      X_boot <- sample(n, X_data, replace = TRUE)
      X_par <- MME_CD(X_boot, N_dist, U_dist, "MME")

    }

  }else{

    X_par <- MME_CD(X_data, N_dist, U_dist, method)

  }
  # DF of the compound distribution
  F_X <- df_X(X_par, N_dist, U_dist, max(X_data), par_method_proxy, method_proxy_df)

  n <- length(X_data)

  X_plus <- sort(X_data[X_data>0])

  n0 <- n - length(X_plus)

  # KS distance
  if(n0 < n){
    KS <- max(
      c(
        (n0 + 0:(n - n0)) / n - sapply(c(0, X_plus), function(x) F_X(x)),
        c(sapply(X_plus, function(x) F_X(x)), 1) - (n0 + 0:(n - n0)) / n)
    )
  }else{
    KS <- NA
  }

  # Parametric bootstrap
  ## Simulation of the bootstrap sample
  KS_boot <- {}

  for(k in 1:B){

    X_boot <- simulate_X(n, N_dist, U_dist, X_par)

    # Inference of the parameters of the compound distribution
    if(sum(X_boot==0)==0){

      X_par_boot <- MME_CD(X_boot, N_dist, U_dist, "MME")

      while((X_par_boot[1]<0|X_par_boot[1]>1) | is.na(X_par_boot[1]) | X_par_boot[3]<0 |
            is.na(X_par_boot[3])| X_par_boot[4]<0 | is.na(X_par_boot[4])){

        X_boot <- simulate_X(n, N_dist, U_dist, X_par)
        X_par_boot <- MME_CD(X_boot, N_dist, U_dist, "MME")

      }

    }else{

      X_par_boot <- MME_CD(X_boot, N_dist, U_dist, method)
      while((X_par_boot[1]<0|X_par_boot[1]>1) | is.na(X_par_boot[1]) | X_par_boot[3]<0 |
            is.na(X_par_boot[3])| X_par_boot[4]<0 | is.na(X_par_boot[4])){

        X_boot <- simulate_X(n, N_dist, U_dist, X_par)
        X_par_boot <- MME_CD(X_boot, N_dist, U_dist, method)

      }

    }

    # Cdf of the compound distribution

    F_X_boot <- df_X(X_par_boot, N_dist, U_dist, max(X_boot), par_method_proxy,
                     method_proxy_df)
    X_boot_plus <- sort(X_boot[X_boot>0])
    n0_boot <- n - length(X_boot_plus)
    # KS distance
    if(n0_boot < n){
      KS_boot[k] <- max(
        c(
          (n0_boot + 0:(n - n0_boot)) / n - sapply(c(0, X_boot_plus), function(x) F_X_boot(x)),
          c(sapply(X_boot_plus, function(x) F_X_boot(x)), 1) - (n0_boot + 0:(n - n0_boot)) / n)
      )
    }else{
      KS_boot[k] <- NA
    }
  }
  return(c(KS, quantile(KS_boot, 0.95, na.rm = TRUE)))
}

# X_data = itamtplcost_monthly$total_claim_size
# X_par = MME_CD(X_data, "pois", "gamma", "partial-MME")
# N_dist = "pois"
# U_dist = "gamma"
# par_method_proxy <- 35
# method_proxy_df <- 'truncation'
# B = 20
# method <- 'partial-MME'
# ks_gof(X_data, N_dist, U_dist, par_method_proxy,
#           method_proxy_df, B, method = 'partial-MME')


#' CvM_gof
#' Goodness-of-fit-test based on the Cramer-von Mises criterium
#' @param X_data Samle of aggregated losses
#' @param N_dist Claim frequency distribution under H0
#' @param U_dist Claim sizes distribution under H0
#' @param par_method_proxy Width of the discretization of the claim sizes distribution
#' when Panjer algorithm is used, order of truncation when the truncation method is used
#' @param method_proxy_df Appproximation method for the df of the compound distribution
#' @param B Number of bootstrap loops
#' @param method MM estimation method c('MME', 'partial-MME')
#'
#' @return Value of the test statistic and critical value of the test
#' @export
#'
#' @examples
cvm_gof <- function(X_data, N_dist, U_dist, par_method_proxy, method_proxy_df, B, method = 'partial-MME'){

  # Inference of the parameters of the compound distribution
  ## If no zeros use the full-MME method otherwise always the other one!
  if(sum(X_data==0)==0){

    X_par <- MME_CD(X_data, N_dist, U_dist, 'MME')
    # Case when the data generates estimates that are negative
    while((X_par[1]<0|X_par[1]>1) | is.na(X_par[1]) | X_par[3]<0 | is.na(X_par[3]) | X_par[4]<0 | is.na(X_par[4])){

      X_boot <- sample(n, X_data, replace = TRUE)
      X_par <- MME_CD(X_boot, N_dist, U_dist, "MME")

    }

  }else{

    X_par <- MME_CD(X_data, N_dist, U_dist, method)


  }

  # DF of the compound distribution
  F_X <- df_X(X_par, N_dist, U_dist, max(X_data), par_method_proxy, method_proxy_df)

  n <- length(X_data)

  X_plus <- sort(X_data[X_data>0])

  n0 <- n - length(X_plus)

  # Cramer-von-Mises criterion
  if(n0 < n){
    CvM <- n * (F_X(0) - n0 / n)^(2) * F_X(0) +
      (1 - F_X(0)) / (n - n0) * (1 / 12 / (n - n0) +
                                   sum(((sapply(X_plus, function(x) F_X(x)) - F_X(0)) / (1 - F_X(0))
                                        -  (2 * 1:(n - n0) - 1) / 2 / (n - n0))^(2)))

  }else{
    CvM <- NA
  }

  # Parametric bootstrap
  ## Simulation of the bootstrap sample
  CvM_boot <- {}
  for(k in 1:B){

    X_boot <- simulate_X(n, N_dist, U_dist, X_par)

    # Inference of the parameters of the compound distribution
    if(sum(X_boot==0)==0){

      X_par_boot <- MME_CD(X_boot, N_dist, U_dist, "MME")

      while((X_par_boot[1]<0|X_par_boot[1]>1) | is.na(X_par_boot[1]) | X_par_boot[3]<0 | is.na(X_par_boot[3]) | X_par_boot[4]<0 | is.na(X_par_boot[4])){

        X_boot <- simulate_X(n, N_dist, U_dist, X_par)
        X_par_boot <- MME_CD(X_boot, N_dist, U_dist, "MME")

      }

    }else{

      X_par_boot <- MME_CD(X_boot, N_dist, U_dist, method)
      while((X_par_boot[1]<0|X_par_boot[1]>1) | is.na(X_par_boot[1]) | X_par_boot[3]<0 | is.na(X_par_boot[3]) | X_par_boot[4]<0 | is.na(X_par_boot[4])){

        X_boot <- simulate_X(n, N_dist, U_dist, X_par)
        X_par_boot <- MME_CD(X_boot, N_dist, U_dist, method)

      }

    }
    # Cdf of the compound distribution
    F_X_boot <- df_X(X_par_boot, N_dist, U_dist, max(X_boot), par_method_proxy, method_proxy_df)
    X_boot_plus <- sort(X_boot[X_boot>0])
    n0_boot <- n - length(X_boot_plus)
    # Cramer-von-Mises criterion
    if(n0_boot < n){
      CvM_boot[k] <- n * (F_X_boot(0) - n0_boot / n)^(2) * F_X_boot(0) +
        (1 - F_X_boot(0)) / (n - n0_boot) * (1 / 12 / (n - n0_boot) +
                                               sum(((sapply(X_boot_plus, function(x) F_X_boot(x)) - F_X_boot(0)) / (1 - F_X_boot(0))
                                                    -  (2 * 1:(n - n0_boot) - 1) / 2 / (n - n0_boot))^(2)))
    }else{
      CvM_boot[k] <- NA
    }

  }
  return(c(CvM, quantile(CvM_boot, 0.95, na.rm = TRUE)))
}


# X_data = itamtplcost_monthly$total_claim_size
# X_par = MME_CD(X_data, "geom", "exp", "MME")
# N_dist = "geom"
# U_dist = "exp"
# beta= 10^(-1)
# scaled = TRUE
# X_data <- rep(0, 100)
# par_method_proxy <- 35
# method_proxy_df <- 'truncation'
# cvm_gof(X_data, N_dist, U_dist, par_method_proxy,
#           method_proxy_df, 2000, method = 'partial-MME')
