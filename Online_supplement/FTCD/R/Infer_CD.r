#File: Infer_CD.r
#Article: Goodness-of-fit tests for compound distributions with applications in insurance
#Package: testCD
#Author: Pierre-O Goffard
#Method of moments::moment estimations of compound distributions


#' MME_CD
#' Infer the parameter of a compound distribution
#' @param X_data Sample of aggregated losses
#' @param N_dist Claim frequency distribution c('pois', 'geom')
#' @param U_dist Claim sizes distribution c('exp', 'gamma', weibull', 'invgauss')
#' @param method MME of partial-MME (usimng the number of zeros in the data)
#'
#' @return the parameters of the compound distribution
#' @export
#'
#' @examples
MME_CD <- function(X_data, N_dist, U_dist, method = 'partial-MME'){

  p0 <- sum(X_data==0) / length(X_data)
  if(N_dist == "pois"){

    if(U_dist == "gamma"){

      if(method == 'partial-MME'){

        a <- 0
        b <- - log(p0)
        m_U <- c(
          (1 - a)/(a + b) * moments::moment(X_data, order = 1, central = FALSE),
          (1 - a)/(a + b) * moments::moment(X_data, order = 2, central = FALSE) -
            (2 * a + b) * (1-a) / (a + b) ^ (2) * moments::moment(X_data, order = 1, central = FALSE)^(2)
        )
        U_par <- c(
          m_U[1]^(2) / (m_U[2] - m_U[1]^(2)),
          (m_U[2] - m_U[1]^(2)) / m_U[1]
        )
        res <- c(a, b, U_par)

      }else if(method == 'MME'){

        a <- 0
        r <- (2 * var(X_data)^(2) - moments::moment(X_data, order = 3, central = TRUE) * mean(X_data)) /
          (moments::moment(X_data, order = 3, central = TRUE) * mean(X_data) - var(X_data)^2)
        theta <- var(X_data) / mean(X_data) / (r + 1)
        b <-  mean(X_data) / theta / r
        U_par = c(r, theta)
        res <- c(a, b, U_par)

      }

    }else if(U_dist == "exp"){

      if(method == 'partial-MME'){

        a <- 0
        b <- - log(p0)
        E_U <- (1 - a)/(a + b) * moments::moment(X_data, order = 1, central = FALSE)
        U_par <-  c(1, E_U)

        res <- c(a, b, U_par)

      }else if(method == 'MME'){

        a <- 0
        b <- 2 * mean(X_data)^(2) / var(X_data)
        theta <- var(X_data) / 2 / mean(X_data)
        U_par <- c(1, theta)
        res <- c(a, b, U_par)

      }
    }else if(U_dist == "weibull"){

      if(method == 'partial-MME'){

        a <- 0
        b <- - log(p0)
        m_U <- c(
          (1 - a)/(a + b) * moments::moment(X_data, order = 1, central = FALSE),
          (1 - a)/(a + b) * moments::moment(X_data, order = 2, central = FALSE) -
            (2 * a + b) * (1-a) / (a + b) ^ (2) * moments::moment(X_data, order = 1, central = FALSE)^(2)
        )
        shape <- uniroot(function(x) gamma(1 + 2 / x) / gamma(1 + 1 / x)^(2)
                         - m_U[2]  / m_U[1]^(2), c(0.05, 500))$root
        scale <- m_U[1] / gamma(1 + 1 / shape)
        U_par <- c(shape, scale)
        res <- c(a, b, U_par)
      }
    }else if(U_dist == "invgauss"){

      if(method == 'partial-MME'){

        a <- 0
        b <- - log(p0)
        m_U <- c(
          (1 - a)/(a + b) * moments::moment(X_data, order = 1, central = FALSE),
          (1 - a)/(a + b) * moments::moment(X_data, order = 2, central = FALSE) -
            (2 * a + b) * (1-a) / (a + b) ^ (2) * moments::moment(X_data, order = 1, central = FALSE)^(2)
        )
        mu = m_U[1]
        phi = (m_U[2]-m_U[1]^(2)) / m_U[1]^(3)
        U_par <- c(mu, phi)
        res <- c(a, b, U_par)
      }
    }
    }else if(N_dist == "geom"){

    if(U_dist == "gamma"){

      a <- 1 - p0
      b <- 0
      m_U <- c(
        (1 - a)/(a + b) * moments::moment(X_data, order = 1, central = FALSE),
        (1 - a)/(a + b) * moments::moment(X_data, order = 2, central = FALSE) -
          (2 * a + b) * (1-a) / (a + b) ^ (2) * moments::moment(X_data, order = 1, central = FALSE)^(2)
      )
      U_par <- c(
        m_U[1]^(2) / (m_U[2] - m_U[1]^(2)),
        (m_U[2] - m_U[1]^(2)) / m_U[1]
      )
      res <- c(p0, b, U_par)

      }else if(U_dist == "exp"){

        if(method == 'partial-MME'){

          a <- 1 - p0
          b <- 0
          E_U <- (1 - a)/(a + b) * moments::moment(X_data, order = 1, central = FALSE)
          U_par <-  c(1, E_U)
          res <- c(p0, b, U_par)

      }else if(method == 'MME'){

        theta <- (var(X_data) - mean(X_data)^(2)) / 2 / mean(X_data)
        a <-  mean(X_data) / (theta + mean(X_data))
        b <- 0
        U_par <- c(1, theta)
        res <- c(1 - a, b, U_par)

        }
      }else if(U_dist == "weibull"){

        if(method == 'partial-MME'){

          a <- 1 - p0
          b <- 0
          m_U <- c(
            (1 - a)/(a + b) * moments::moment(X_data, order = 1, central = FALSE),
            (1 - a)/(a + b) * moments::moment(X_data, order = 2, central = FALSE) -
              (2 * a + b) * (1-a) / (a + b) ^ (2) * moments::moment(X_data, order = 1, central = FALSE)^(2)
          )

          shape <- uniroot(function(x) gamma(1 + 2 / x) / gamma(1 + 1 / x)^(2) - m_U[2]  / m_U[1]^(2), c(0.05, 2000))$root
          scale <- m_U[1] / gamma(1 + 1 / shape)
          U_par <- c(shape, scale)
          res <- c(p0, b, U_par)
        }
      }else if(U_dist == "invgauss"){

        if(method == 'partial-MME'){

          a <- 1 - p0
          b <- 0
          m_U <- c(
            (1 - a)/(a + b) * moments::moment(X_data, order = 1, central = FALSE),
            (1 - a)/(a + b) * moments::moment(X_data, order = 2, central = FALSE) -
              (2 * a + b) * (1-a) / (a + b) ^ (2) * moments::moment(X_data, order = 1, central = FALSE)^(2)
          )
          mu = m_U[1]
          phi = (m_U[2]-m_U[1]^(2)) / m_U[1]^(3)
          U_par <- c(mu, phi)
          res <- c(p0, b, U_par)
          }
      }
    }
  return(res)
}


#' Poisson_gamma_to_Poisson_exp
#' Parameters of a Compound Poisson exponential model obtained by matching th emoment of the compound-Gamma model
#' @param lam Intensity of the Poisson distribution
#' @param r shape of the gamma distribution
#' @param theta scale parameter of the gamma distribution
#'
#' @return The parameters of the Poisson exponential model
#' @export
#'
#' @examples
Poisson_gamma_to_Poisson_exp <- function(lam, r, theta){
  theta_pois_exp <- (r + 1) * theta / 2
  lam_pois_exp <- 2 * lam * r / (r + 1)

  return(c(0, lam_pois_exp, 1, theta_pois_exp))
}


