
##################################################
# Alternative claim sizes distribution = Gamma,
# Inverse Gaussian and gamma
##################################################
library(FTCD)
shape_vec <- seq(0.2, 2, 0.2)
nb_shape = length(shape_vec)
dist_U_vec <- c('gamma', 'weibull', 'invgauss')
dist_U_H0_vec <- c('exp', 'gamma', 'invgauss')
MC <- 10000

##################################################
# Comparison of the GOF procedure based on the CDF
##################################################
methods <- c("CvM", "KS")
truncation_order <- 35
res_CDF <- NULL
for(method in methods){
  print(method)
  for(dist_U in dist_U_vec){
    print(paste("The alternative distribution is", dist_U))
    for(dist_U_H0 in dist_U_H0_vec){
      print(paste("The distribution under H0 is", dist_U_H0))
      eval(parse(text = paste0(
"power_CDF <- data.frame(
U_alt = rep('", dist_U,"',",nb_shape,"),
shape_U =shape_vec,
U_H0 = rep('", dist_U_H0,"',",nb_shape,"),
method_GOF = rep('",method,"',",nb_shape,"),
beta = rep(", truncation_order,",", nb_shape,"),
power = sapply(shape_vec, function(shape)
p_value(100, 'pois', '",dist_U,"',c(0, 1, shape, 1), 'pois', '",dist_U_H0,"',", MC,",'partial-MME','",method,"',", truncation_order,")))")))
      res_CDF <- rbind(res_CDF, power_CDF)
    }
  }
}

###############################################
# Power of the LT test based on the L2 distance
###############################################
beta_vec <- c(10^(-5),10^(-4), 10^(-3), 10^(-2), 0.1, 1, 10)
method <- "LT_L2"
res_LT_L2 <- NULL
for(beta in beta_vec){
  print(beta)
  for(dist_U in dist_U_vec){
    print(paste("The alternative distribution is", dist_U))
    for(dist_U_H0 in dist_U_H0_vec){
      print(paste("The distribution under H0 is", dist_U_H0))
      eval(parse(text = paste0(
"power_LT_L2 <- data.frame(
U_alt = rep('", dist_U,"',",nb_shape,"),
shape_U =shape_vec,
U_H0 = rep('", dist_U_H0,"',",nb_shape,"),
method_GOF = rep('",method,"',",nb_shape,"),
beta = rep(",beta,",",nb_shape,"),
power = sapply(shape_vec, function(shape)
p_value(100, 'pois', '",dist_U,"',c(0, 1, shape, 1), 'pois', '",dist_U_H0,"',", MC,",'partial-MME','",method,"',",beta,")))")))
      res_LT_L2 <- rbind(res_LT_L2, power_LT_L2)
    }
  }
}

###############################################
# Power of the LT test based on the ODE distance
###############################################
beta_vec <- c(10^(-3), 10^(-2), 10^(-1), 1, 10, 100)
method <- "LT_ODE"
res_LT_ODE <- NULL
for(beta in beta_vec){
  print(beta)
  for(dist_U in dist_U_vec){
    print(paste("The alternative distribution is", dist_U))
    for(dist_U_H0 in dist_U_H0_vec){
      print(paste("The distribution under H0 is", dist_U))
      eval(parse(text = paste0(
        "power_LT_ODE <- data.frame(
U_alt = rep('", dist_U,"',",nb_shape,"),
shape_U =shape_vec,
U_H0 = rep('", dist_U_H0,"',",nb_shape,"),
method_GOF = rep('",method,"',",nb_shape,"),
beta = rep(",beta,",",nb_shape,"),
power = sapply(shape_vec, function(shape)
p_value(100, 'pois', '",dist_U,"',c(0, 1, shape, 1), 'pois', '",dist_U_H0,"',", MC,",'partial-MME','",method,"',",beta,")))")))
      res_LT_ODE <- rbind(res_LT_ODE, power_LT_ODE)
    }
  }
}



##################################################
# Alternative claim sizes distribution = lognormal
##################################################
shape_vec <- seq(0.2, 2, 0.2)
nb_shape = length(shape_vec)
dist_U_vec <- c('lnorm')
dist_U_H0_vec <- c('exp', 'gamma', 'invgauss')
MC <- 10000

##################################################
# Comparison of the GOF procedure based on the CDF
##################################################
methods <- c("CvM", "KS")
truncation_order <- 35
res_CDF <- NULL
for(method in methods){
  print(method)
  for(dist_U in dist_U_vec){
    print(paste("The alternative distribution is", dist_U))
    for(dist_U_H0 in dist_U_H0_vec){
      print(paste("The distribution under H0 is", dist_U_H0))
      eval(parse(text = paste0(
        "power_CDF <- data.frame(
        U_alt = rep('", dist_U,"',",nb_shape,"),
        shape_U =shape_vec,
        U_H0 = rep('", dist_U_H0,"',",nb_shape,"),
        method_GOF = rep('",method,"',",nb_shape,"),
        beta = rep(", truncation_order,",", nb_shape,"),
        power = sapply(shape_vec, function(shape)
        p_value(100, 'pois', '",dist_U,"',c(0, 1, 0, shape), 'pois', '",dist_U_H0,"',", MC,",'partial-MME','",method,"',", truncation_order,")))")))
      res_CDF_lnorm <- rbind(res_CDF, power_CDF)
    }
  }
  }

###############################################
# Power of the LT test based on the L2 distance
###############################################
beta_vec <- c(10^(-4), 10^(-3), 10^(-2))
method <- "LT_L2"
res_LT_L2 <- NULL
for(beta in beta_vec){
  print(beta)
  for(dist_U in dist_U_vec){
    print(paste("The alternative distribution is", dist_U))
    for(dist_U_H0 in dist_U_H0_vec){
      print(paste("The distribution under H0 is", dist_U_H0))
      eval(parse(text = paste0(
        "power_LT_L2 <- data.frame(
        U_alt = rep('", dist_U,"',",nb_shape,"),
        shape_U =shape_vec,
        U_H0 = rep('", dist_U_H0,"',",nb_shape,"),
        method_GOF = rep('",method,"',",nb_shape,"),
        beta = rep(",beta,",",nb_shape,"),
        power = sapply(shape_vec, function(shape)
        p_value(100, 'pois', '",dist_U,"',c(0, 1, 0, shape), 'pois', '",dist_U_H0,"',", MC,",'partial-MME','",method,"',",beta,")))")))
      res_LT_L2_lnorm <- rbind(res_LT_L2, power_LT_L2)
    }
  }
  }

###############################################
# Power of the LT test based on the ODE distance
###############################################
beta_vec <- c( 10^(-2), 10^(-1), 1, 10)
method <- "LT_ODE"
res_LT_ODE <- NULL
for(beta in beta_vec){
  print(beta)
  for(dist_U in dist_U_vec){
    print(paste("The alternative distribution is", dist_U))
    for(dist_U_H0 in dist_U_H0_vec){
      print(paste("The distribution under H0 is", dist_U))
      eval(parse(text = paste0(
        "power_LT_ODE <- data.frame(
        U_alt = rep('", dist_U,"',",nb_shape,"),
        shape_U =shape_vec,
        U_H0 = rep('", dist_U_H0,"',",nb_shape,"),
        method_GOF = rep('",method,"',",nb_shape,"),
        beta = rep(",beta,",",nb_shape,"),
        power = sapply(shape_vec, function(shape)
        p_value(100, 'pois', '",dist_U,"',c(0, 1, 0, shape), 'pois', '",dist_U_H0,"',", MC,",'partial-MME','",method,"',",beta,")))")))
      res_LT_ODE_lnorm <- rbind(res_LT_ODE, power_LT_ODE)
    }
  }
}

write_csv(rbind(res_CDF, res_LT_L2, res_LT_ODE, res_CDF_lnorm, res_LT_L2_lnorm, res_LT_ODE_lnorm),
          path = "Comparison_power_claim_sizes.txt")

power <- read_csv("Comparison_power_claim_sizes.txt")

