#######################################################################
####################  calculate RCB and RMSDR  ######################
#######################################################################

source("function_measures.R")
load("sim_results_beta.RData")


meas_sub50 <- lapply(results_sub50, calc.measures.all)
meas_sub63 <- lapply(results_sub63, calc.measures.all)
meas_sub80 <- lapply(results_sub80, calc.measures.all)
meas_boot <- lapply(results_boot, calc.measures.all)

save(meas_sub50, meas_sub63, meas_sub80, meas_boot, file="sim_measures.RData")


