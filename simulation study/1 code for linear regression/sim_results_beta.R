#######################################################################
##############  calculate betas for simulations study  ################
#######################################################################


load("setup.RData")

source("function_beta.R")

results_sub50 <- list()
results_sub63 <- list()
results_sub80 <- list()
results_boot <- list()

# run all simulation szenarios for each subsampling approach and the bootstrap
## It takes weeks to run the code on a single PC. We ran it on the Vienna Scientific Cluster.

for(i in 1:6){
  results_sub50[[i]] <- simulate(setup=setup[i,], nsim=1000, nboot=1000, sampling="subsampling50")  
  results_sub63[[i]] <- simulate(setup=setup[i,], nsim=1000, nboot=1000, sampling="subsampling")
  results_sub80[[i]] <- simulate(setup=setup[i,], nsim=1000, nboot=1000, sampling="subsampling80")
  results_boot[[i]] <- simulate(setup=setup[i,], nsim=1000, nboot=1000, sampling="bootstrap")
}


# save the results

save(results_sub50, results_sub63, results_sub80, results_boot, file="sim_results_beta.RData")


