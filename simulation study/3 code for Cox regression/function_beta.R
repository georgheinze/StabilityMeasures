####################################################################################
##############  function to calculate betas for simulations study  ################
####################################################################################

library(mvtnorm)
library(glmnet)
library(parallel)
library(abe)
library(survival)


simulate <-
  function(setup,
           nsim = 40,
           nboot = 10,
           sampling = "bootstrap") {
    # extract information from 'setup'
    
    n <- setup["n"]
    seed <- setup["seed"]
    beta <- setup[3:length(setup)]
    
    global <-
      as.formula(paste("Surv(y, cens) ~ ", paste(names(beta), collapse = "+")))
    true <-
      as.formula(paste("Surv(y, cens) ~ ", paste(names(beta[beta != 0]), collapse =
                                                   "+")))
    
    
    # set parameters for parallelization
    
    w <- detectCores() - 1
    npar <- floor(nsim / w)
    npar_last <- nsim %% w
    npar_cl <- rep(npar, w)
    if (npar_last != 0)
      npar_cl[1:npar_last] <- npar + 1
    
    
    set.seed(seed)
    
    
    cl <- makeCluster(w, type = "PSOCK")
    clusterExport(
      cl,
      c("w", "npar",  "true", "global", "beta",
        "seed", "n", "nsim", "nboot", "sampling"),
      envir = environment()
    )
    clusterEvalQ(cl, library("glmnet"))
    clusterEvalQ(cl, library("mvtnorm"))
    clusterEvalQ(cl, library("survival"))
    clusterEvalQ(cl, library("abe"))
    
    res_boot <- parLapply(cl, npar_cl, function(nsimx) {
      # sometimes you have to use the absolute path here
      source("function_setup.R")
      source("function_glmnetAIC_cox.R")
      
    
      # predefine objects
      
      lacoef <- matrix(0, nsimx, length(beta))
      colnames(lacoef) <- names(beta)
      fuse <-
        trse <-
        bwAICcoef <- bw5coef <- fucoef <- trcoef <- laAICcoef <- lacoef
      
      data_sim <- list()
      
      aux_mat <- matrix(
        0,
        ncol = length(beta),
        nrow = nboot,
        dimnames = list(NULL, names(beta))
      )
      
      boot_bwAIC_est <-
        boot_bw5_est <-
        boot_la_est <-
        boot_laAIC_est <-
        boot_fu_est <-
        boot_tr_est <- replicate(nsimx, aux_mat, simplify = F)
      
      for (isim in 1:nsimx) {
        
        #simulate data with the simulation setup
        x <- pretrans.x(gen.x.from.z(gen.z(n)))
        pred <- partial.predictor(x)
        ystar <- apply(pred, 1, sum)
        dat <- x.to.data(x = x, y = ystar)
        linpred <- as.matrix(dat[, 2:ncol(dat)]) %*% beta - 3.5
        y <- (-log(runif(length(linpred))) / (1 / 5 * exp(linpred))) ** (1 /
                                                                           3)
        cens <- rep(1, length = length(linpred))
        tau <- 3.85
        futime <- 0.001 + runif(length(linpred)) * tau
        cens[y > futime] <- 0
        y[cens == 0] <- futime[cens == 0]
        dat$y <- y
        dat <- cbind(dat, cens)
        dat <<- dat
        
        data_sim[[isim]] <- dat
        
        
        # fit models and obtain betas
        
        ## global model
        fu <-
          coef(fufit <- coxph(
            formula = global,
            data = dat,
            x = T,
            y = T
          ))
        fucoef[isim,] <- fu
        fuse[isim,] <-
          summary(fufit)$coefficients[, "se(coef)"][names(fu)]
        
        ## true model
        tr <- coef(trfit <- coxph(formula = true, data = dat))
        trcoef[isim,] <- tr[names(beta)]
        trcoef[is.na(trcoef)] <- 0
        colnames(trcoef) <- names(beta)
        trse[isim,] <-
          summary(trfit)$coefficients[, "se(coef)"][names(beta)]
        trse[is.na(trse)] <- 0
        colnames(trse) <- names(beta)
        
        ## model with BE(AIC)
        bwAICfit <-
          try(abe(
            fufit,
            data = dat,
            tau = Inf,
            criterion = "AIC",
            verbose = F
          ))
        if (class(bwAICfit) != "try-error") {
          bwAICcoef[isim, names(bwAIC)] <- bwAIC <- coef(bwAICfit)
        } else{
          bwAICcoef[isim,] <- bwAIC <- rep(NA, time = length(beta))
        }
        
        
        ## model with BE(0.05)
        bw5fit <-
          try(abe(
            coxph(global, data = dat, x = T),
            data = dat,
            tau = Inf,
            criterion = "alpha",
            alpha = 0.05,
            verbose = F
          ))
        if (class(bw5fit) != "try-error") {
          bw5coef[isim, names(bw5)] <- bw5 <- coef(bw5fit)
        } else{
          bw5fit[isim,] <- bw5 <- rep(NA, time = length(beta))
        }
        
        
        ## model with Lasso
        surv_ob <- cbind(dat$y, dat$cens)
        colnames(surv_ob) <- c("time", "status")
        
        la <-
          (coef(
            lafit <-
              cv.glmnet(
                y = surv_ob,
                x = as.matrix(dat[, names(beta)]),
                family = "cox"
              ),
            s = "lambda.min"
          ))
        lacoef[isim,] <- la[names(beta), ]
        
        ## model with Lasso(AIC)
        ### we do not report this in the article or in the supplement
        laAIC <-
          laAICfit <-
          AIC.glmnet(y = surv_ob,
                     x = as.matrix(dat[, names(beta)]),
                     family = "cox")
        laAICcoef[isim,] <- laAICfit$beta.opt
        
        
        # do subsampling or bootstrapping and obtain betas with every resample
        for (i in 1:nboot) {
          if (sampling == "bootstrap")
            data_id <- sample(1:dim(dat)[1], nrow(dat), replace = T)
          if (sampling == "subsampling")
            data_id <- sample(1:dim(dat)[1], nrow(dat) * 0.632, replace = F)
          if (sampling == "subsampling50")
            data_id <- sample(1:dim(dat)[1], nrow(dat) * 0.5, replace = F)
          if (sampling == "subsampling80")
            data_id <- sample(1:dim(dat)[1], nrow(dat) * 0.8, replace = F)
          
          
          boot_data <- dat[data_id,]
          boot_data <<- boot_data
          
          ## global model
          boot_fu_mod <- coxph(
            formula = global,
            data = boot_data,
            singular.ok = T,
            x = T,
            y = T
          )
          boot_fu_est[[isim]][i, names(boot_fu_mod$coefficients)] <-
            boot_fu_mod$coefficients
          
          ## true model
          boot_tr_mod <- coxph(formula = true, data = boot_data)
          boot_tr_est[[isim]][i, names(boot_tr_mod$coefficients)] <-
            boot_tr_mod$coefficients
          
          ## model with BE(AIC)
          boot_bwAIC_mod <- try(abe(
            boot_fu_mod,
            data = boot_data,
            tau = Inf,
            criterion = "AIC",
            verbose = F
          ))
          if (class(boot_bwAIC_mod) != "try-error") {
            boot_bwAIC_est[[isim]][i, names(boot_bwAIC_mod$coefficients)] <-
              boot_bwAIC_mod$coefficients
          }
          
          ## model with BE(0.05)
          boot_bw5_mod <- try(abe(
            boot_fu_mod,
            data = boot_data,
            verbose = F,
            tau = Inf,
            criterion = "alpha",
            alpha = 0.05
          ))
          
          if (class(boot_bw5_mod) != "try-error") {
            boot_bw5_est[[isim]][i, names(boot_bw5_mod$coefficients)] <-
              boot_bw5_mod$coefficients
          }
          
          
          ## model with Lasso
          boot_surv_ob <- cbind(dat[data_id, "y"], dat[data_id, "cens"])
          colnames(boot_surv_ob) <- c("time", "status")
          
          boot_la_est[[isim]][i, ] <-
            as.vector(coef(
              cv.glmnet(
                y = boot_surv_ob,
                x = as.matrix(dat[data_id, names(beta)]),
                family = "cox"
              ),
              s = "lambda.min"
            ))
          
          ## model with Lasso(AIC)
          boot_laAIC_est[[isim]][i, ] <-
            as.vector(AIC.glmnet(
              y = boot_surv_ob,
              x = as.matrix(dat[data_id, names(beta)]),
              family = "cox"
            )$beta.opt)
          print(paste("isim=", isim, ", nboot=", i))
        }
        
        
        print(paste("isim=", isim))
        
        
      }
      
      
      return(
        list(
          data = data_sim,
          fucoef = fucoef,
          fuse = fuse,
          trcoef = trcoef,
          trse = trse,
          bwAICcoef = bwAICcoef,
          bw5coef = bw5coef,
          lacoef = lacoef,
          laAICcoef = laAICcoef,
          boot_fu_coef = boot_fu_est,
          boot_tr_coef = boot_tr_est,
          boot_bwAIC_coef = boot_bwAIC_est,
          boot_bw5_coef = boot_bw5_est,
          boot_la_coef = boot_la_est,
          boot_laAIC_coef = boot_laAIC_est
        )
      )
    })
    
    
    stopCluster(cl)
    
 
    
    # merge the results of all parallelized simulations   
    
    la_betas <- matrix(0, nsim * (nboot + 1), 2 + length(beta))
    colnames(la_betas) <- c("simulation", "bootstrap",  names(beta))
    la_betas[, "simulation"] <- rep(1:nsim, each = (nboot + 1))
    la_betas[, "bootstrap"] <- rep(c(0:nboot), times = nsim)
    bwAIC_betas <-
      bw5_betas <- tr_betas <- fu_betas <- laAIC_betas <- la_betas
    
    fu_betas <- matrix(0, nsim * (nboot + 2), 3 + length(beta))
    colnames(fu_betas) <-
      c("simulation", "bootstrap",  "se", names(beta))
    fu_betas[, "simulation"] <- rep(1:nsim, each = (nboot + 2))
    fu_betas[, "bootstrap"] <- rep(c(0, 0:nboot), times = nsim)
    fu_betas[, "se"] <- rep(c(0, 1, rep(0, times = nboot)), times = nsim)
    
    
    data_sim <- list()
    
    fu_col <- 4:ncol(fu_betas)
    fu_row <- cumsum(c(1, rep(nboot + 2, times = nsim - 1)))
    fuse_row <- cumsum(c(2, rep(nboot + 2, times = nsim - 1)))
    fuboot_row <- cumsum(c(3, rep(nboot + 2, times = nsim - 1)))
    
    all_col <- 3:ncol(la_betas)
    all_row <- cumsum(c(1, rep(nboot + 1, times = nsim - 1)))
    allboot_row <- cumsum(c(2, rep(nboot + 1, times = nsim - 1)))
    
    for (clust in 1:w) {
      if (clust == 1)
        start <- 1
      else
        start <-  cumsum(npar_cl[1:clust])[clust - 1] + 1
      end <- cumsum(npar_cl[1:clust])[clust]
      range <- start:end
      
      data_sim[range] <- res_boot[[clust]]$data
      
      
      fu_boot_rows <-
        c(apply(cbind(fuboot_row[range], (
          fuboot_row[range] + nboot - 1
        )),
        1, function(x)
          c(x[1]:x[2])))
      fu_betas[fu_row[range], fu_col] <- res_boot[[clust]]$fucoef
      fu_betas[fuse_row[range], fu_col] <- res_boot[[clust]]$fuse
      fu_betas[fu_boot_rows , fu_col] <-
        do.call(rbind, res_boot[[clust]]$boot_fu_coef)
      
      
      all_boot_rows <-
        c(apply(cbind(
          allboot_row[range], (allboot_row[range] + nboot - 1)
        ),
        1, function(x)
          c(x[1]:x[2])))
      tr_betas[all_row[range], all_col] <- res_boot[[clust]]$trcoef
      tr_betas[all_boot_rows, all_col] <-
        do.call(rbind, res_boot[[clust]]$boot_tr_coef)
      bwAIC_betas[all_row[range], all_col] <-
        res_boot[[clust]]$bwAICcoef
      bwAIC_betas[all_boot_rows, all_col] <-
        do.call(rbind, res_boot[[clust]]$boot_bwAIC_coef)
      bw5_betas[all_row[range], all_col] <- res_boot[[clust]]$bw5coef
      bw5_betas[all_boot_rows, all_col] <-
        do.call(rbind, res_boot[[clust]]$boot_bw5_coef)
      la_betas[all_row[range], all_col] <- res_boot[[clust]]$lacoef
      la_betas[all_boot_rows, all_col] <-
        do.call(rbind, res_boot[[clust]]$boot_la_coef)
      laAIC_betas[all_row[range], all_col] <-
        res_boot[[clust]]$laAICcoef
      laAIC_betas[all_boot_rows, all_col] <-
        do.call(rbind, res_boot[[clust]]$boot_laAIC_coef)
      
    }
    
    
    
    
    return(
      list(
        n = n,
        nsim = nsim,
        nboot = nboot,
        seed = seed,
        beta = beta,
        fu_betas = fu_betas,
        tr_betas = tr_betas,
        bwAIC_betas = bwAIC_betas,
        bw5_betas = bw5_betas,
        la_betas = la_betas,
        laAIC_betas = laAIC_betas
      )
    )
    
  }



#load("setup_surv_20180312.RData")
#setup[,"seed"] <- c(11:17)

#test <- simulate(setup=setup[1,], nsim =24, nboot = 10, sampling="subsampling")
#nsim=100
#nboot=100
