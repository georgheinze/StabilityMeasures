##################################################################################
####################  functions to calculate VIF, RCB and RMSDR  #################
##################################################################################

library(data.table)



# compute mean of coefficients obtained with the resamples
boot_beta_mean <- function(matrix, beta) {
  do.call(cbind, lapply(apply(matrix[bootstrap != 0, mget(names(beta))], 2,
                              function(x)
                                aggregate(x, by = list(matrix[bootstrap != 0, simulation]),
                                          function(z)
                                            mean(z, na.rm = T))),
                        function(y)
                          y[, 2]))
}


# compute variable inclusion frequency (VIF)
bif <- function(matrix, beta) {
  do.call(cbind, lapply(apply(matrix[bootstrap != 0, mget(names(beta))], 2,
                              function(x)
                                aggregate(x, by = list(matrix[bootstrap != 0, simulation]),
                                          function(z)
                                            sum(z != 0, na.rm = T) / length(z[!is.na(z)]))),
                        function(y)
                          y[, 2]))
}


# compute relative conditional bias (RCB)
rcb <- function(matrix, fu, beta) {
  tf <- fu$bootstrap == 0 & fu$se == 0
  (boot_beta_mean(matrix, beta) /
      (bif(matrix, beta) * fu[tf, mget(names(beta))]) - 1) * 100
}

# compute root mean squared difference (RMSD)
rmsd <- function(matrix, fu, nboot, beta) {
  tf <- fu$bootstrap == 0 & fu$se == 0
  beta_hats <- fu[tf, mget(names(beta))]
  
  beta_hats_large <-
    beta_hats[rep(seq_len(nrow(beta_hats)), each = nboot),]
  
  betas_diff <-
    matrix[bootstrap != 0, mget(names(beta))] - beta_hats_large
  
  msd <-
    aggregate(betas_diff, by = list(matrix[bootstrap != 0, simulation]),
              function(x)
                mean(x ^ 2, na.rm = T))
  rmsd <- sqrt(msd[,-1])
  
  return(rmsd)
}


# compute RMSD ratio
rmsdratio <- function(rmsd, fu, beta) {
  sapply(rmsd, as.numeric) / sapply(fu[se == 1, mget(names(beta))], as.numeric)
}



# compute the approximated estimand for RCB
rcb_population <- function(matrix, beta, nsim) {
  tf <- matrix$bootstrap == 0
  if (sum(colnames(matrix) == "se") > 0)
    tf <- tf & matrix$se == 0
  
  beta_large <-
    matrix(rep(beta, each = nsim), nrow = nsim)
  
  betas_diff <- (matrix[tf, mget(names(beta))] - beta_large)
  
  betas_zeros <-
    apply(matrix[tf, mget(names(beta))], 2, function(x)
      x != 0)
  
  betas_sum <- betas_diff_sum <- NULL
  for (i in 1:length(beta)) {
    ni <- names(beta)[i]
    betas_diff_sum[i] <-
      sum(betas_diff[, get(ni)][betas_zeros[, ni]], na.rm = T)
    betas_sum[i] <- sum(beta_large[, i][betas_zeros[, ni]], na.rm = T)
  }
  
  rcb_pop <- betas_diff_sum / betas_sum * 100
  rcb_pop[is.infinite(rcb_pop) | is.nan(rcb_pop)] <- NA
  return(rcb_pop)
  
}



# compute the approximated estimand of RMSD ratio
rmsdratio_population <- function(matrix, fu, beta, nsim) {
  beta_large <-
    matrix(rep(beta, each = nsim), nrow = nsim)
  
  tf <- matrix$bootstrap == 0
  if (!is.null(matrix$se))
    tf <- tf & matrix$se == 0
  tf_fu <- fu$bootstrap == 0 & fu$se == 0
  
  nom <- denom <- NULL
  
  for (i in 1:length(beta)) {
    ni <- names(beta)[i]
    nom[i] <-
      sum((matrix[tf, mget(names(beta))][, get(ni)] - beta_large[, i]) ^ 2, na.rm =
            T)
    denom[i] <-
      sum((fu[tf_fu, mget(names(beta))][, get(ni)] - beta_large[, i]) ^ 2, na.rm =
            T)
  }
  
  rmsdratio <- sqrt(nom / denom)
  
  return(rmsdratio)
}


# calculate these measures per estimation/selection procedure
calc.measures <- function(matrix, fu, beta, nsim, nboot) {
  res_bif <- bif(matrix, beta)
  res_rcb <- rcb(matrix, fu, beta)
  res_rmsd <- rmsd(matrix, fu, nboot, beta)
  res_rmsdratio <- rmsdratio(res_rmsd, fu, beta)
  res_rcb_pop <- rcb_population(matrix, beta, nsim)
  res_rmsdratio_pop <- rmsdratio_population(matrix, fu, beta, nsim)
  
  return(
    list(
      bif = res_bif,
      rcb = res_rcb,
      rmsd = res_rmsd,
      rmsdratio = res_rmsdratio,
      rcb_pop = res_rcb_pop,
      rmsdratio_pop = res_rmsdratio_pop
    )
  )
}




# calculate these measures for all simulation szenarios
calc.measures.all <- function(results) {
  nsim <- results$nsim
  nboot <- results$nboot
  beta <- results$beta
  
  
  fu <- data.table(results$fu_betas)
  tr <- data.table(results$tr_betas)
  bwAIC <- data.table(results$bwAIC_betas)
  bw5 <- data.table(results$bw5_betas)
  la <- data.table(results$la_betas)
  laAIC <- data.table(results$laAIC_betas)
  
  res_fu <- calc.measures(fu, fu, beta, nsim, nboot)
  res_tr <- calc.measures(tr, fu, beta, nsim, nboot)
  res_bwAIC <- calc.measures(bwAIC, fu, beta, nsim, nboot)
  res_bw5 <- calc.measures(bw5, fu, beta, nsim, nboot)
  res_la <- calc.measures(la, fu, beta, nsim, nboot)
  res_laAIC <- calc.measures(laAIC, fu, beta, nsim, nboot)
  
  return(
    list(
      fu = res_fu,
      tr = res_tr,
      bwAIC = res_bwAIC,
      bw5 = res_bw5,
      la = res_la,
      laAIC = res_laAIC
    )
  )
  
}
