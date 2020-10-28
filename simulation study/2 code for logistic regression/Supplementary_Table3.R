###############################################################################
##############  Supplementary Table 3: Model selection frequencies ############
###############################################################################
library(xlsx)


source("function_beta.R")
load("sim_results_beta.RData")
load("setup.RData")

effect <- setup[1, 3:ncol(setup)] != 0


# compute MSF

corr.mod <- function(object = results_boot, method = "bwAIC_betas") {
  
  # extract all betas obtained with resamples
  mat <-
    lapply(object, function(x)
      cbind(x[[method]][x[[method]][, "bootstrap"] != 0, names(effect)],
            "sim" = x[[method]][x[[method]][, "bootstrap"] !=
                                  0, "simulation"]))
  
  # extract all betas obtained in a simulated data set ('true betas')
  mat_tr <-
    lapply(object, function(x)
      x[[method]][x[[method]][, "bootstrap"] == 0, names(effect)])
  
  # how many resampled models include the correct variables?
  correct <-
    sapply(mat, function(y)
      sum(apply(y[,-ncol(y)], 1, function(x)
        identical(x != 0, effect))))

  
  # how many models with the simulated data sets include the correct variables?  
  correct_tr <-
    sapply(mat_tr, function(y)
      sum(apply(y, 1, function(x)
        identical(x != 0, effect))))
  

  
  # MSF in each simulation run for each simulation szenario
  correct2 <-
    sapply(mat, function(y)
      aggregate(
        apply(y[, -ncol(y)], 1, function(x)
          identical(x != 0, effect)),
        by = list(y[, "sim"]),
        FUN = mean
      )[, 2])
  
  # MSF estimand for each simulation szenario
  correct2_tr <- correct_tr/ nrow(mat_tr[[1]])
  
  

  # compute the RMSE
  rmse_correct <-
    apply((100 * (correct2 - matrix(
      rep(correct2_tr, each = 1000), nrow = 1000
    ))) ^ 2, 2, function(x)
      sqrt(mean(x)))
  
  
  # merge results
  res <- round(correct / nrow(mat[[1]]) * 100, 2)
  res_tr <- round(correct_tr / nrow(mat_tr[[1]]) * 100, 2)
  rmse_correct <-  round(rmse_correct, 2)
  
  return(cbind(res, res_tr, rmse_correct))
}


# run the function for each subsampling approach and the bootstrap
sub_bwAIC_mat50 <- corr.mod(results_sub50, "bwAIC_betas")[1:3, ]
sub_bwAIC_mat63 <- corr.mod(results_sub63, "bwAIC_betas")[1:3, ]
sub_bwAIC_mat80 <- corr.mod(results_sub80, "bwAIC_betas")[1:3, ]
boot_bwAIC_mat <- corr.mod(results_boot, "bwAIC_betas")[1:3, ]



# combine results
tab_msf <- function(boot, sub50, sub63, sub80) {
  msf <- rbind(
    "Estimand" = apply(rbind(sub50[, 2], sub63[, 2], sub80[, 2], boot[, 2]), 2, mean),
    "S5" = sub50[, 1],
    "S5 RMSE" = sub50[, 3],
    "S6" = sub63[, 1],
    "S6 RMSE" = sub63[, 3],
    "S8" = sub80[, 1],
    "S8 RMSE" = sub80[, 3],
    "B" = boot[, 1],
    "B RMSE" = boot[, 3]
  )
  
  colnames(msf) <- paste0("N=", c(300, 750, 1000))
  return(msf)
  
}


# save results as .xlsx


x="bwAIC"
  
write.xlsx(
    tab_msf(boot=get(paste0("boot_",x,"_mat")), 
            sub50=get(paste0("sub_",x,"_mat50")),
            sub63=get(paste0("sub_",x,"_mat63")), 
            sub80=get(paste0("sub_",x,"_mat80"))), 
    file = paste0("MSF_",x,".xlsx"))



