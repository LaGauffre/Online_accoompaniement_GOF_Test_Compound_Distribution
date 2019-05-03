#File: Simulate_CD.r
#Article: Goodness-of-fit tests for compound distributions with applications in insurance
#Package: testCD
#Author: Pierre-O Goffard
#Functions to generate sample of aggregated claim sizes



#' sum_U
#' Sum N iid observations drawn from the distribution U_dist
#'
#' @param N Number of claim sizes
#' @param U_dist Distribution of the claim sizes c('exp', 'gamma', 'weibull', invgauss)
#' @param U_par Parameter of the distribution of the claim sizes
#'
#' @return The sum of N claim sizes
#' @export
#'
#' @examples
sum_U <- function(N, U_dist, U_par){

  if(U_dist == "exp"){

    res <- sum(rgamma(N, shape = 1, scale = U_par[2]))

  } else if(U_dist == "gamma"){

    res <- sum(rgamma(N, shape = U_par[1], scale = U_par[2]))

  } else if(U_dist == "weibull"){

    res <- sum(rweibull(N, shape = U_par[1], scale = U_par[2]))

  } else if(U_dist == "lnorm"){

    res <- sum(rlnorm(N, meanlog = U_par[1], sdlog = U_par[2]))

  } else if(U_dist == "invgauss"){

    res <- sum(actuar::rinvgauss(N, mean = U_par[1], dis = U_par[2]))

  }

  return(res)

}


#' simulate_X
#' Generate a sample of aggregated claim sizes
#'
#' @param n Size of the sample
#' @param N_dist Claim frequency distribution
#' c('pois', 'geom', 'dunif', 'nbinom', 'zmpois', 'mpois')
#' @param U_dist Claim sizes distribution c('exp', 'gamma', 'lnorm', 'weibull', 'invgauss')
#' @param X_par Parameters of the aggregated claim sizes distribution
#'
#' @return Sample of size n of the aggregated claim sizes
#' @export
#'
#' @examples
simulate_X <- function(n, N_dist, U_dist, X_par){

  if(N_dist == 'geom'){

    res <- sapply(rgeom(n, X_par[1]), function(claim_number) sum_U(claim_number, U_dist, X_par[3:4]))

  }else if(N_dist == 'pois'){

    res <- sapply(rpois(n, X_par[2]), function(claim_number) sum_U(claim_number, U_dist, X_par[3:4]))

  } else if(N_dist == 'dunif'){

    res <- sapply(rdunif(n, X_par[2], a = X_par[1]),
                  function(claim_number) sum_U(claim_number, U_dist, X_par[3:4]))
  }else if(N_dist == 'zmpois'){

    res <- sapply(actuar::rzmpois(n, X_par[1], X_par[2]),
                  function(claim_number) sum_U(claim_number, U_dist, X_par[3:4]))

  } else if(N_dist == 'nbinom'){

    res <- sapply(rnbinom(n,
                          size = X_par[2],
                          prob =  X_par[1]),
                  function(claim_number) sum_U(claim_number, U_dist, X_par[3:4]))

  } else if(N_dist == 'mpois'){

    bernoulli <- rbinom(n, size = 1, prob = X_par[1])
    res <- sapply((1-bernoulli)*rpois(n, 1) + bernoulli*rpois(n, X_par[2]),
                  function(claim_number) sum_U(claim_number, U_dist, X_par[3:4]))

  }
  return(res)
}



