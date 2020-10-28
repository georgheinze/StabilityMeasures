###############################################################################
################## Modelling with synthesized data set #######################
###############################################################################

# Note: The original data set was synthesized. The results of the synthesized data set obtained with this code will be different to those in the article.

library(foreach)
library(doParallel)
library(doRNG)
library(data.table)
library(survival)
library(shrink)
library(mgcv)


# load synthesized data
dataset <- read.csv(file = "dataset_synthesized.csv")


# design variables
effects <- c(
  "age10",
  "cholratio_log",
  "triglycerides10",
  "sbp10",
  "bpmedication",
  "smoke",
  "diabetes",
  "bmi",
  "physicalactivity",
  "waist",
  "uprotein",
  "uglucose"
)


form <- as.formula(paste0("Surv(time_cvd, cvd) ~ ",
                          paste(effects, collapse = " + ")))



surv_ob <- cbind(dataset$time_cvd, dataset$cvd)
colnames(surv_ob) <- c("time", "status")


set.seed(298043)

# fit a global model
FULL_fit <- coxph(form, data = dataset, x = T, y = T)
FULL <- coef(FULL_fit)
FULL_se <- summary(FULL_fit)$coefficients[, "se(coef)"]

# fit a model with BE(AIC)
BE_AIC_fit <- step(FULL_fit,
                   k = 2,
                   direction = "backward",
                   trace = F)
BE_AIC <- coef(BE_AIC_fit)[names(FULL)]
BE_AIC_se <-
  summary(step(
    FULL_fit,
    k = 2,
    direction = "backward",
    trace = F
  ))$coefficients[, "se(coef)"][names(FULL)]


# predefine objects
boot_data <-
  data_id_boot <- data_id_sub <- boot_data <- sub_data <- NULL
boot_FULL <-  boot_BE_AIC <- boot_FULL_se <-
sub_FULL <-  sub_BE_AIC <- sub_FULL_se <- list()


cl <- makeCluster(detectCores())
clusterExport(cl, ls()[!grepl("cl", ls())], envir = environment())
clusterEvalQ(cl, library(data.table))
clusterEvalQ(cl, library(rms))
registerDoParallel(cl)


set.seed(123656)

#parallelize calculation with bootstrap and subsampling
bootnum <- 1000

res <- foreach(i = 1:bootnum) %dorng% {
  data_id_boot <-
    sample(1:nrow(dataset), nrow(dataset), replace = T)
  boot_data <- dataset[data_id_boot,]
  boot_data <<- boot_data
  
  data_id_sub  <- sample(1:nrow(dataset), nrow(dataset) * 0.5)
  sub_data  <- dataset[data_id_sub,]
  sub_data <<- sub_data
  
  
  boot_FULL_fit <- coxph(form, data = boot_data, x = T, y = T)
  boot_FULL <- coef(boot_FULL_fit)
  boot_FULL_se <- summary(boot_FULL_fit)$coefficients[, "se(coef)"]
  
  boot_BE_AIC <- coef(step(
    boot_FULL_fit,
    k = 2,
    direction = "backward",
    trace = F
  ))[names(coef(boot_FULL_fit))]
  
  
  
  sub_FULL_fit <- coxph(form, data = sub_data, x = T, y = T)
  sub_FULL <- coef(sub_FULL_fit)
  sub_FULL_se <- summary(sub_FULL_fit)$coefficients[, "se(coef)"]
  
  sub_BE_AIC <- coef(step(
    sub_FULL_fit,
    k = 2,
    direction = "backward",
    trace = F
  ))[names(coef(sub_FULL_fit))]
  
  
  list(boot_FULL,
       boot_FULL_se,
       boot_BE_AIC,
       sub_FULL,
       sub_FULL_se,
       sub_BE_AIC)
  
}

stopCluster(cl)



names_res <- c("boot_FULL",
               "boot_FULL_se",
               "boot_BE_AIC",
               "sub_FULL",
               "sub_FULL_se",
               "sub_BE_AIC")

for (k in 1:length(names_res)) {
  assign(names_res[k], t(sapply(res, function(x)
    x[[k]])))
}



#save.image("modelling.RData")

##### compute stablity measures ########

source("function_measures.R")

boot_BE_AIC[is.na(boot_BE_AIC)] <- 0
sub_BE_AIC[is.na(sub_BE_AIC)] <- 0

# pairwise shrinkage factors
BE_shrink <- shrink(BE_AIC_fit)


# Bootstrap
summ_measures_BE_boot <-
  matrix_summary(boot_BE_AIC, FULL, FULL_se, BE_AIC, BE_AIC_se)

# Subsampling
summ_measuresBE <-
  matrix_summary(sub_BE_AIC, FULL, FULL_se, BE_AIC, BE_AIC_se)


# Combination
## VIF, MSF with Subsampling
## RCB & RMSDR with Bootstrap

ordering <- rownames(summ_measuresBE)

# calculate SD for RCB; SE for shrinkage factors (SF) using the delta-method
rcbs <-
  (boot_BE_AIC - matrix(rep(FULL, bootnum), nrow = bootnum, byrow = T)) /
  matrix(rep(FULL, bootnum), nrow = bootnum, byrow = T) * (boot_BE_AIC != 0)

# RCB for every resample
rcbs <- ifelse(rcbs == 0, NA, rcbs)
colnames(rcbs) <- names(FULL)

# mean RCB
rcb_mean <- apply(rcbs, 2, function(x)
  mean(x, na.rm = T))

# Standard deviation for RCB
rcb_sd <- apply(rcbs, 2, function(x)
  sd(x, na.rm = T))

# Delta Method
# Tanner 1993, Tools for Statistical Inference, p. 17
# Variance of f(x) = f'(x)^2*B
## B is the variance of RCB (RCB=x and SF=f(x))
## f(x)=1/(1+x) ## formula for SF
## f'(x)=-1/(1+x)^(-2)
## using the formula gives:
sf_se <- sqrt(1 / (1 + rcb_mean) ^ 4 * rcb_sd ^ 2)


# Table 4 in the article
round(
  cbind(
    summ_measuresBE[, 1:5],
    summ_measures_BE_boot[ordering, 6:7],
    "SF" = 1 / (1 + summ_measures_BE_boot[ordering, "Relative conditional bias (%)"] /
                  100),
    "SF SE" = sf_se[ordering],
    "PWSF" = BE_shrink$ShrinkageFactors[ordering],
    "PWSF SE" = sqrt(diag(BE_shrink$ShrinkageFactorsVCOV))[ordering]
  ),
  3
)




# Model selection frequency
mat_01 <- apply(sub_BE_AIC, 2, function(x)
  1 * (x != 0))
colnames(mat_01) <- names(FULL)

sel_mod <- 1 * (!is.na(BE_AIC))
names(sel_mod) <- names(FULL)

msf <-
  sum(apply(mat_01, 1, function(x)
    identical(x, sel_mod))) / nrow(mat_01)

msf * 100



other_msf <- uniquecombs(mat_01)
# model selection frequencies of the 5 most selected combinations
(table(attr(other_msf, "index")) / 
    nrow(mat_01))[order(table(attr(other_msf, "index")), decreasing = T)][1:5]
identical(other_msf[10, ], sel_mod) # the final selected model is not the most frequently occurring model in the subsampling approach
# most frequently occurring combinations
cbind(sel_mod,
      other_msf[10, ],
      other_msf[1, ],
      other_msf[5, ],
      other_msf[4, ],
      other_msf[11, ])


# pairwise inclusion frequencies
# Supplementary Table 5


vif_sub <- summ_measuresBE[, 3]
effects_ord  <- names(vif_sub)[order(vif_sub, decreasing = T)]
mat_01 <- mat_01[, effects_ord]

boot_pairfreq <-
  matrix(
    100,
    ncol = ncol(mat_01),
    nrow = ncol(mat_01),
    dimnames = list(colnames(mat_01), colnames(mat_01))
  )

pval <- 0.01
combis <- combn(colnames(mat_01), 2)
expect_pairfreq <- NULL

for (i in 1:dim(combis)[2]) {
  boot_pairfreq[combis[1, i], combis[2, i]] <-
    sum(apply(mat_01[, combis[, i]], 1, sum) == 2) / bootnum * 100
  expect_pairfreq[i] <-
    vif_sub[grepl(combis[1, i], effects_ord)][1] *
    vif_sub[grepl(combis[2, i], effects_ord)][1] / 100
  boot_pairfreq[combis[2, i], combis[1, i]] <-
    ifelse(
      is(suppressWarnings(try(chisq.test(mat_01[, combis[1, i]],
                                         mat_01[, combis[2, i]]),
                              silent = T)), "try-error"), NA,
      ifelse(
        suppressWarnings(chisq.test(mat_01[, combis[1, i]],
                                    mat_01[, combis[2, i]])$p.value) > pval,
        "",
        ifelse(as.numeric(boot_pairfreq[combis[1, i], combis[2, i]]) <
                 expect_pairfreq[i], "-", "+")
      )
    )
}

diag(boot_pairfreq) <- vif_sub[effects_ord]


boot_pairfreq

