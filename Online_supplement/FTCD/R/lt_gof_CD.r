#File: lt_gof_CD.r
#Article: Goodness-of-fit tests for compound distributions with applications in insurance
#Package: testCD
#Author: Pierre-O Goffard
#LT based goodness-of-fit-test for compound distributions


#' LT_gamma
#' Laplace transform of the gamma distribution
#'
#' @param t point at which the Laplace transform is evaluated
#' @param r Shape parameter of the gamma distribution
#' @param theta Scale parameter of the gamma distribution
#'
#' @return The value of the Laplace transform at t
#' @export
#'
#' @examples
LT_gamma <- function(t, r, theta){
  return((1 + theta * t)^(-r))
}


#' LT_invgauss
#' Laplace transform of the inverse gaussian distribution
#' @param t point at which the Laplace transform is evaluated
#' @param mu mean parameter of the inverse gaussian distribution
#' @param phi dispersion parameter of the inverse gaussian
#'
#' @return Value of the Laplace transform at t
#' @export
#'
#' @examples
LT_invgauss <- function(t, mu, phi){
  return(exp(1 / mu / phi * (1 - sqrt(1 + 2 * mu^(2) * t * phi))))
}

#' dLT_gamma
#' First derivative of the Laplace transform of the gamma distribution
#'
#' @param t Point at which the derivative of the Laplace transform is evaluated
#' @param r Shape parameter of the gamma distribution
#' @param theta Scale paameter of the gamma distribution
#'
#' @return the value of the first derivative of the Laplace transform of a gamma distribuytion at t
#' @export
#'
#' @examples
dLT_gamma <- function(t, r, theta){
  return(- r * theta * (1 + theta * t)^(-r - 1))
}

#' dLT_invgauss
#' First derivative of the Laplace transform of the inverse Gaussian distribution
#' @param t point at which the first derivative of the Laplace transform of the inverse
#' Gaussian distribution is evaluated
#' @param mu Mean parameter of the inverse Gaussian distribution
#' @param phi Dispersion parameter of the inverse Gausian distribution
#'
#' @return The value of first derivative of the Laplace transform of the inverse
#' Gaussian distribution at t
#' @export
#'
#' @examples
dLT_invgauss <- function(t, mu, phi){
  return(- mu / sqrt(1 + 2 * mu^2 * phi * t) * LT_invgauss(t, mu, phi))
}

#' K_1_gamma
#' Function for the LT based goodness of fit test for the KF-gamma distribution
#'
#' @param x Point at which the function is evaluated
#' @param r Shape paremeter of the gamma distribution
#' @param theta Scale parameter of the gamma distribution
#' @param a Parameter of the claim frequency distribution
#'
#' @return The value of the function K_1 at x
#' @export
#'
#' @examples
K_1_gamma <- function(x, r, theta, a){

  A <- - exp(x / theta) / r * theta^(r) / x^(r+2) * expint::gammainc(r + 2, x / theta)

  B <- - (x + theta) / x^(2) /r / theta

  return(A - a * B)
}

#' K_2_gamma
#' Function to evaluate the test statistic in the LT based GOF test with a distance defined as an ODE for the
#'  KF-Gamma model
#'
#' @param x Point at which the function is evaluated
#' @param r Shape parameter of the gamma distribution
#' @param theta Scale parameter of the gamma distribution
#' @param a Parameter of the claim frequency distribution
#'
#' @return the value of the function at x
#' @export
#'
#' @examples
K_2_gamma <- function(x, r, theta, a){
  C <- exp(x / theta) / r^(2) * theta^(2 * r) / x^(2 * r + 3) * expint::gammainc(2 * r + 3, x / theta)
  D <- exp(x / theta) / r^(2) * theta^(r) / x^(r + 3) * expint::gammainc(r + 3, x / theta)
  E <- (x^(2) + 2 * theta * x + 2 * theta^(2)) / r^(2) / theta^(2) / x^(3)
  return(C - 2 * a * D + a^(2) * E)
}

#' K_1_invgauss
#' Function for the evaluation of test statistic for the LT based test using an ODE
#' as a distance when testing a KF(a,b)-invgauss(mu, phi) model
#'
#' @param x point at which the function K_1 is evaluated
#' @param mu Mean parameter of the inverse gaussian distribution governing the claim sizes under
#' H0
#' @param phi Dispersion parameter of the inverse gaussian distribution governing the
#' claim sizes under H0
#' @param a Parameter of the claim frequency distribution
#'
#' @return the value of the function at any given point
#' @export
#'
#' @examples
K_1_invgauss <- function(x, mu, phi, a){
  c <- (x - mu) / mu / sqrt(phi*x)
  d <- mu^2 * phi / x
  A <- - exp(1 / 2 / d - 1 / mu / phi + 1 / 2 / phi / x)/ mu^(2) / sqrt(x * phi) *
    (
      d * (c * exp(-c^2 / 2) + sqrt(2 * pi) * (1 - pnorm(c))) +
        2 * mu^(2) * sqrt(phi) / x / sqrt(x) * exp(- c^(2) / 2) +
        mu^(2) / x^(2) * sqrt(2 * pi) * (1  - pnorm(c))
      )
  B <- - sqrt(phi) / x/ sqrt(x) * exp(1 / 2/ d) *
    (
      1 / sqrt(d) * exp( -1 / 2 / d) + sqrt(2 * pi) * (1 - pnorm(1 / sqrt(d)))
      )

  return(A - a * B)
}

#' K_2_invgauss
#' Function for the evaluation of the test statistic associated to the LT based goodness of fit test for the
#' KF(a,b)-Inverse Gaussian model
#'
#' @param x Point at which the function is evaluated
#' @param mu Mean parameter of the inverse gaussian distribution
#' @param phi Dispersion parameter of the inverse gaussian distribution
#' @param a Parameter of the claim frequency distribution
#'
#' @return The value of the funtion at x
#' @export
#'
#' @examples
K_2_invgauss <- function(x, mu, phi, a){
  d <- mu^2 * phi / x
  e <- exp(- 2 / mu / phi + 1 / 2 / d + 2 / x/ phi) / sqrt(phi * x) / x^(3) * 2^(3)
  f <- (x - 2 * mu) / mu / sqrt(phi * x)
  A <- e * (
    sqrt(2 * pi) * (1 - pnorm(f)) +
      3 * sqrt(phi * x) / 2 * exp(- f^(2) / 2) +
      3 * phi * x / 4 * (f * exp(- f^(2)/2) + sqrt(2 * pi)* (1 - pnorm(f))) +
      (sqrt(phi * x) / 2)^(3) * (f^(2)+2)* exp(-f^(2)/2)
  )

  g <- exp(- 1 / mu / phi + 1 / 2 / d + 1 / x / phi / 2) / sqrt(phi * x) / x^(3)
  h <- (x - mu) / mu / sqrt(phi * x)
  B <- g * (
    sqrt(2 * pi) * (1 - pnorm(h)) +
      3 * sqrt(phi * x) * exp(- h^(2) / 2) +
      3 * phi * x * (h * exp(- h^(2)/2) + sqrt(2 * pi)* (1 - pnorm(h))) +
      (sqrt(phi * x))^(3) * (h^(2)+2)* exp(-h^(2)/2)
  )

   C <- 1 / mu^(2) / x + 2 * phi / x^(2)

  return(A - 2 * a * B + a^(2) * C)
}



#' Compute_T_n
#' Compute the test statistics for the KF-Gamma model based on an ODE
#' @param X_data Sample of aggregated losses
#' @param X_par Parameter of the compound distribution
#' @param N_dist Distribution if the claim sizes c('pois', 'geom')
#' @param U_dist Distribution if the claim sizes c('exp', 'gamma', invgauss)
#' @param beta Parameter in the weight function
#' @param scaled Indicates wether or not the data is being scaled
#'
#' @return The value of the test statistic
#' @export
#'
#' @examples
Compute_T_n <- function(X_data, X_par, N_dist, U_dist, beta, scaled){

  n <- length(X_data)
  if(N_dist == 'pois'){
    a <- 0
    b <- X_par[2]
  }else if(N_dist == 'geom'){
    a <- 1- X_par[1]
    b <- 0
  }

  if(U_dist == 'gamma' | U_dist == 'exp'){
    if(scaled){

      Y <- X_data / X_par[4]
      theta <- 1

    }else{

      Y <- X_data
      theta <- X_par[4]

    }

    Y_sum <- matrix(Y, n, n) + t(matrix(Y, n, n))
    Y_sum_vec <- as.vector(Y_sum)
    ID <- unique(Y_sum_vec)
    temp <- ID + beta

    r <- X_par[3]

    K_2 <-  K_2_gamma(temp, r, theta, a)

    K_1 <- K_1_gamma(temp, r, theta, a)
  }
  else if(U_dist == 'invgauss'){
    if(scaled){

      Y <-  X_data / X_par[3]
      mu <- 1
      phi <- X_par[3] * X_par[4]

    }else{

      Y <-  X_data
      mu <- X_par[3]
      phi <- X_par[4]

    }

    Y_sum <- matrix(Y, n, n) + t(matrix(Y, n, n))
    Y_sum_vec <- as.vector(Y_sum)
    ID <- unique(Y_sum_vec)
    temp <- ID + beta

    K_1 <- K_1_invgauss(temp, mu, phi, a)
    K_2 <- K_2_invgauss(temp, mu, phi, a)


  }


  K1_vec <- K_1[match(x = Y_sum_vec, ID)]
  K2_vec <- K_2[match(x = Y_sum_vec, ID)]
  K_2_mat <- matrix(K2_vec, n, n)
  K_1_mat <- matrix(K1_vec, n, n)

  Y_prod <- Y %*% t(Y)
  Y_prod_1 <- Y %*% t(rep(1, n))

  T_n <- sum(Y_prod * K_2_mat) / n  +
    2 * (a + b) * (Y %*% colSums(K_1_mat)) / n +
    (a + b)^(2) * sum(1 / (Y_sum + beta)) / n

  return(T_n)
}


#' lt_ODE_gof
#' Compute the test statistic and the critical value of the test via a parametric bootstrap routine
#' when the distance is based on an ODE
#'
#' @param X_data Sample of aggregated losses
#' @param N_dist Claim frequency distribution c('pois', 'geom')
#' @param U_dist Claim sizes distribution c('exp', 'gamma', 'invgauss')
#' @param beta Parameter of the exponential weight function
#' @param B Number of loops in the parametric bootstrap routine
#' @param method Method of moment estimation method c('MME', 'partial-MME')
#' @param scaled Indicates wether or not the data is being scaled
#'
#' @return The value of the test statistic and the critical value at a 95% confidence level
#' @export
#'
#' @examples
lt_ODE_gof <- function(X_data, N_dist, U_dist, beta, B, method = 'partial-MME', scaled){
  n <- length(X_data)
  # Inference of the parameters of the compound distribution
  if(sum(X_data==0)==0){

    X_par <- MME_CD(X_data, N_dist, U_dist, "MME")

    # Case when the data generates estimates that are negative
    while((X_par[1]<0|X_par[1]>1) | is.na(X_par[1]) | X_par[3]<0 | is.na(X_par[3]) | X_par[4]<0 | is.na(X_par[4])){

      X_boot <- sample(n, X_data, replace = TRUE)
      X_par <- MME_CD(X_boot, N_dist, U_dist, "MME")

    }

    }else{

      X_par <- MME_CD(X_data, N_dist, U_dist, method)

    }

  # Test statistic on the data
  T_n <- Compute_T_n(X_data, X_par, N_dist, U_dist, beta, scaled)

  # Parametric bootstrap for loop
  T_n_boot <- {}
  for(k in 1:B){

    if(scaled){

      if(U_dist == 'exp' | U_dist == 'gamma'){

        X_boot <- simulate_X(n, N_dist, U_dist, c(X_par[1:3],1))

      }else if(U_dist == 'invgauss'){

        X_boot <- simulate_X(n, N_dist, U_dist, c(X_par[1:2], 1, X_par[3] * X_par[4]))

      }

    }else{

      X_boot <- simulate_X(n, N_dist, U_dist, X_par)

    }

    if(sum(X_boot==0)==0){

      X_par_boot <- MME_CD(X_boot, N_dist, U_dist, "MME")

      while((X_par_boot[1]<0|X_par_boot[1]>1) | is.na(X_par_boot[1]) | X_par_boot[3]<0 | is.na(X_par_boot[3]) | X_par_boot[4]<0 | is.na(X_par_boot[4])){

        X_boot <- simulate_X(n, N_dist, U_dist, X_par)
        X_par_boot <- MME_CD(X_boot, N_dist, U_dist, "MME")

      }

      }else{

        X_par_boot <- MME_CD(X_boot, N_dist, U_dist, method)
        while((X_par_boot[1]<0|X_par_boot[1]>1) | is.na(X_par_boot[1]) | X_par_boot[3]<0 | is.na(X_par_boot[3]) |
              X_par_boot[4]<0 | is.na(X_par_boot[4])){
          X_boot <- simulate_X(n, N_dist, U_dist, X_par)
          X_par_boot <- MME_CD(X_boot, N_dist, U_dist, method)
          }

        }

    T_n_boot[k] <- Compute_T_n(X_boot, X_par_boot, N_dist, U_dist, beta, scaled)

  }
  return(
    c(T_n, quantile(T_n_boot, 0.95, na.rm = TRUE))
  )
}


# X_data = itamtplcost_small_monthly$total_claim_size
# X_par = MME_CD(X_data, "pois", "gamma", "partial-MME")
# X_data = simulate_X(length(X_data), 'pois', 'gamma', c(X_par[1:3], 1))
# N_dist = "pois"
# U_dist = "gamma"
# beta= 0.1
# B =20
# method = 'partial-MME'
# scaled = T
# lt_ODE_gof(X_data, 'pois', 'gamma', beta, B, method = 'partial-MME', scaled)

#' Compute_S_n_pois
#' Compute the test statistic of the LT based GOF when the distance is
#' the integrated squared difference between the empirical Laplace transform
#' and the theoretical Laplace transform
#'
#' @param X_data Sample of aggregated losses
#' @param X_par Parameter of the compound distribution
#' @param N_dist Claim frequency distribution c('pois', 'geom')
#' @param U_dist Claim sizes distribution c('exp', 'gamma', 'invgauss')
#' @param beta Parameter of the weight function
#' @param scaled boolean indicating wether a scale is applied to the data or not
#'
#' @return The value of the test statistic
#' @export
#'
#' @examples
Compute_S_n <- function(X_data, X_par, N_dist, U_dist, beta, scaled){
  n <- length(X_data)
  if(N_dist == 'pois'){
    if(U_dist == 'gamma' | U_dist == 'exp'){
      if(scaled){

        r <- X_par[3]
        theta <- 1
        lam <- X_par[2]
        Y <- X_data / X_par[4]

      }else{

        r <- X_par[3]
        theta <- X_par[4]
        lam <- X_par[2]
        Y <- X_data

      }


      C <- sum(sapply(Y,
                      function(x) integrate(
                        function(t) exp(lam * ((1 + theta * t)^(-r) - 1)) * exp(-(x + beta) * t)
                        , lower = 0, upper = Inf, stop.on.error = FALSE)$value))
      D <- integrate(
        function(t) exp(2 * lam * ((1 + theta * t)^(-r) - 1)) * exp(- beta * t),
        lower = 0, upper = Inf, stop.on.error = FALSE)$value


    }else if(U_dist == 'invgauss'){
      if(scaled){

        Y <-  X_data / X_par[3]
        lam <- X_par[2]
        mu <- 1
        phi <- X_par[3] * X_par[4]

      }else{

        mu <- X_par[3]
        phi <- X_par[4]
        lam <- X_par[2]
        Y <- X_data

      }


      C <- sum(sapply(Y,
                      function(x) integrate(
                        function(t) exp(lam * (LT_invgauss(t, mu, phi) - 1)) * exp(-(x + beta) * t)
                        , lower = 0, upper = Inf, stop.on.error = FALSE)$value))
      D <- integrate(
        function(t) exp(2 * lam * (LT_invgauss(t, mu, phi)-1)) * exp(- beta * t),
        lower = 0, upper = Inf, stop.on.error = FALSE)$value
      }
    }else if(N_dist == 'geom'){

      if(U_dist == 'gamma' | U_dist == 'exp'){
        if(scaled){

          Y <-  X_data / X_par[3]
          p <- X_par[1]
          theta <- 1
          r <- X_par[3]

        }else{
          Y <- X_data
          r <- X_par[3]
          theta <- X_par[4]
          p <- X_par[1]

      }

        C <- sum(sapply(Y,
                        function(x) integrate(
                          function(t) p / (1 - (1 - p) * (1 + theta * t)^(-r))
                          * exp(-(x + beta) * t),
                          lower = 0, upper = Inf, stop.on.error = FALSE)$value))
        D <- integrate(
          function(t) p^2 / (1 - (1 - p) * (1 + theta * t)^(-r))^2 * exp(- beta * t),
          lower = 0, upper = Inf, stop.on.error = FALSE)$value


      }else if(U_dist == 'invgauss'){
        if(scaled){

          Y <-  X_data / X_par[3]
          p <- X_par[1]
          mu <- 1
          phi <- X_par[3] * X_par[4]

        }else{

          Y <- X_data
          mu <- X_par[3]
          phi <- X_par[4]
          p <- X_par[1]
        }


        C <- sum(sapply(Y,
                        function(x) integrate(
                          function(t) p / (1 - (1 - p) * LT_invgauss(t, mu, phi))
                          * exp(-(x + beta) * t)
                          , lower = 0, upper = Inf, stop.on.error = FALSE)$value))
        D <- integrate(
          function(t) p^2 / (1 - (1 - p) * LT_invgauss(t, mu, phi))^2 * exp(- beta * t),
          lower = 0, upper = Inf, stop.on.error = FALSE)$value
        }
      }
  X_sum <- matrix(Y, n, n) + t(matrix(Y, n, n))
  S_n <- sum(1 / (X_sum + beta)) / n - 2 * C + n * D
  return(S_n)

}

#' lt_L2_gof_pois
#' Compute the test statistic and critical value via a parametric bootstrap routine
#' for the LT based GOF test with a distance defined as an integrated square differrence
#' @param X_data Sample of aggregated losses
#' @param N_dist Claim sizes distribution c('pois', 'geom')
#' @param U_dist Claim sizes distribution c('exp', 'gamma', 'invgauss')
#' @param beta Parameter in the weight function
#' @param B Number of loops in the parametric bootstrap
#' @param method Method of moment estimation c('MME', 'partial-MME')
#' @param scaled boolean to indicate wether the data is scaled
#'
#' @return the test statistic and the critrical value at a 95% confidence level
#' @export
#'
#' @examples
lt_L2_gof <- function(X_data, N_dist, U_dist, beta, B, method = 'partial-MME', scaled){

  n <- length(X_data)
  # Inference of the parameters of the compound distribution
  if(sum(X_data==0)==0){

    X_par <- MME_CD(X_data, N_dist, U_dist, "MME")
    # Case when the data generates estimates that are negative
    while((X_par[1]<0|X_par[1]>1) | is.na(X_par[1]) | X_par[3]<0 | is.na(X_par[3]) | X_par[4]<0 | is.na(X_par[4])){

      X_boot <- sample(n, X_data, replace = TRUE)
      X_par <- MME_CD(X_boot, N_dist, U_dist, "MME")

    }

  }else{
    X_par <- MME_CD(X_data, N_dist, U_dist, method)

  }
  # Test statistic on the data

  S_n <- Compute_S_n(X_data, X_par, N_dist, U_dist, beta, scaled)
  # Parametric bootstrap for loop
  S_n_boot <- {}
  for(k in 1:B){

    X_boot <- simulate_X(n, N_dist, U_dist, X_par)

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
    #T_n_fast(X_boot, params_X_boot, loss_dist, beta)
    S_n_boot[k] <- Compute_S_n(X_boot, X_par_boot, N_dist, U_dist, beta, scaled)
  }
  return(
    c(S_n, quantile(S_n_boot, 0.95, na.rm = TRUE))
  )
}

# X_data = itamtplcost_small_monthly$total_claim_size
# N_dist = 'poisson'
# U_dist = 'invgauss'
# beta = 10^(-3)
# B = 200
# method = 'partial-MME'
# scaled = T
# lt_L2_gof(itamtplcost_monthly$total_claim_size, 'pois', U_dist, beta, B, method = 'partial-MME', scaled)
# X_data = simulate_X(50, "geom", "exp", c(0.9, 0, 1, 1))
# X_par = MME_CD(X_data, "geom", "exp", "MME")
# N_dist = "geom"
# U_dist = "exp"
# beta= 10^(-1)
# scaled = TRUE
# X_data <- rep(0, 100)
# lt_L2_gof(X_data, N_dist, U_dist, beta, 1, method = 'partial-MME', scaled)





