library(FTCD)
inference_methods <- c("MME", "partial-MME")
methods <- c("LT_L2", "LT_ODE")
rvec <- c(0.5, 0.75, 1, 2, 4)
MC <- 10000
beta_vec <- 10^seq(-13 , 3, 1)
res <- NULL

for(method in methods){
  print(method)
  for(r in rvec){
    print(r)
    for(inference_method in inference_methods){
      print(inference_method)
      eval(
        parse(
          text = paste0("power_",
                        r,
                        "_",
                        method,
                        " <- data.frame(
                        r = rep(",r,", length(beta_vec)),
                        LT_gof_method = rep('", method,"', length(beta_vec)),
                        Inference_method = rep('", inference_method,"', length(beta_vec)),
                        beta = beta_vec,
                        power = sapply(beta_vec, function(beta) p_value(100, 'pois', 'gamma',c(0, 1, ",
                        r,
                        ", 1), 'pois', 'exp', MC,'",
                        inference_method,
                        "', '",
                        method,
                        "', beta)))")
        )
      )
      res <- rbind(res, eval(parse(text = paste0("power_",r,"_",method))))
    }

  }
}

write_csv(res, path = "~/Goffard-Jammalamadaka-Meintanis/Simulation_DATA_Export/Calibration.txt")

#Calibration <- read.csv("~/Goffard-Jammalamadaka-Meintanis/Simulation_DATA_Export/Calibration.txt")



