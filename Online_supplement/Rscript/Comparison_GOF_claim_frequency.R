####################################################################################
# power of the GOF test when testing for Poisson-exp when it is ZIP or mixed Poisson
###################################################################################
library(FTCD)
p_vec <- seq(0.1, 0.9, 0.1)
nb_p <- length(p_vec)
dist_U_vec <- c('exp')
dist_N_vec <- c('zmpois', 'mpois')
MC <- 10000

##################################################
# Comparison of the GOF procedure based on the CDF
##################################################
methods_GOF <- c("CvM", "KS")
methods_infer <- c('MME', 'partial-MME')
truncation_order <- 50
res_CDF <- NULL
for(method_infer in methods_infer){
  print(method_infer)
  for(method_GOF in methods_GOF){
    print(method_GOF)
    for(dist_U in dist_U_vec){
      print(paste("The claim sizes distribution is", dist_U))
      for(dist_N in dist_N_vec){
        print(paste("The claim frequency distribution under H0 is", dist_N))
        power_CDF <- data.frame(
          U_alt = rep( dist_U,nb_p),
          param_N = p_vec,
          method_inference = rep(method_infer, nb_p),
          U_H0 = rep(dist_U, nb_p),
          N_H0 = rep('pois', nb_p),
          N_alt = rep(dist_N,nb_p),
          method_GOF = rep(method_GOF,nb_p),
          beta = rep(truncation_order, nb_p),
          power = sapply(p_vec, function(p) p_value(100, dist_N, dist_U, c(p, 5, 1, 1), 'pois', 'exp', MC, method_infer, method_GOF, truncation_order)))
        res_CDF <- rbind(res_CDF, power_CDF)
      }
    }
  }
}


###############################################
# Power of the LT test based on the L2 distance
###############################################
beta_vec <- c(10^(-3), 10^(-2))
methods_GOF <- c("LT_L2")
methods_infer <- c('MME', 'partial-MME')
res_LT_L2 <- NULL
for(beta in beta_vec){
  print(beta)
  for(method_infer in methods_infer){
    print(method_infer)
    for(method_GOF in methods_GOF){
      print(method_GOF)
      for(dist_U in dist_U_vec){
        print(paste("The claim sizes distribution is", dist_U))
        for(dist_N in dist_N_vec){
          print(paste("The claim frequency distribution under H0 is", dist_N))
          power_LT_L2 <- data.frame(
            U_alt = rep( dist_U,nb_p),
            param_N = p_vec,
            method_inference = rep(method_infer, nb_p),
            U_H0 = rep(dist_U,nb_p),
            N_H0 = rep('pois',nb_p),
            N_alt = rep(dist_N,nb_p),
            method_GOF = rep(method_GOF,nb_p),
            beta = rep(beta, nb_p),
            power = sapply(p_vec, function(p) p_value(100, dist_N, dist_U, c(p, 5, 1, 1), 'pois', 'exp', MC, method_infer, method_GOF, beta)))
          res_LT_L2 <- rbind(res_LT_L2, power_LT_L2)
        }
      }
    }
  }
}


###############################################
# Power of the LT test based on the ODE distance
###############################################
beta_vec <- c(10^(-1), 1)
methods_GOF <- c("LT_ODE")
methods_infer <- c('MME', 'partial-MME')
res_LT_ODE <- NULL
for(beta in beta_vec){
  print(beta)
  for(method_infer in methods_infer){
    print(method_infer)
    for(method_GOF in methods_GOF){
      print(method_GOF)
      for(dist_U in dist_U_vec){
        print(paste("The claim sizes distribution is", dist_U))
        for(dist_N in dist_N_vec){
          print(paste("The claim frequency distribution under H0 is", dist_N))
          power_LT_ODE <- data.frame(
            U_alt = rep(dist_U,nb_p),
            param_N = p_vec,
            method_inference = rep(method_infer, nb_p),
            U_H0 = rep(dist_U,nb_p),
            N_H0 = rep('pois',nb_p),
            N_alt = rep(dist_N,nb_p),
            method_GOF = rep(method_GOF,nb_p),
            beta = rep(beta, nb_p),
            power = sapply(p_vec, function(p) p_value(100, dist_N, dist_U, c(p, 5, 1, 1), 'pois', 'exp', MC, method_infer, method_GOF, beta)))
          res_LT_ODE <- rbind(res_LT_ODE, power_LT_ODE)
        }
      }
    }
  }
}

# write_csv(rbind(res_CDF, res_LT_L2, res_LT_ODE),
#           path = "Comparison_power_pois.txt")

Comparison_power_claim_frequency <- read_csv("Comparison_power_Claim_Frequency.txt")



