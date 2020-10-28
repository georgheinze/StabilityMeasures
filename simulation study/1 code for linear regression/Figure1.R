###############################################################################################
##############  Figure 1: VIF for S0.5, S0.632, S0.8 and the bootstrap ########################
###############################################################################################

library(ggplot2)
library(gridExtra)
library(grid)
library(reshape)


source("function_beta.R")
load("sim_results_beta.RData")
load("setup.RData")

effect <- setup[1, 3:ncol(setup)] != 0


############# estimates ##################


# compute estimates and estimand 

corr.beta <- function(object = results_boot, method = "la_betas") {
  
  # extract all betas obtained with resamples
  mat <-
    lapply(object, function(x)
      x[[method]][x[[method]][, "bootstrap"] != 0, names(effect)])
  
  # extract all betas obtained in a simulated data set ('true betas')
  mat_tr <-
    lapply(object, function(x)
      x[[method]][x[[method]][, "bootstrap"] == 0, names(effect)])
  
  # vIF estimates
  varselfreq <-
    t(sapply(mat, function(y)
      apply(y, 2, function(x)
        mean(x != 0)))) * 100
  
  # approximated estimand
  varselfreq_tr <-
    t(sapply(mat_tr, function(y)
      apply(y, 2, function(x)
        mean(x != 0)))) * 100
  
  rownames(varselfreq) <-
    paste("n=", lapply(object, function(x)
      x$n), sep = "")
  rownames(varselfreq_tr) <-
    paste("n=", lapply(object, function(x)
      x$n), ", true", sep = "")
  
  num <- length(lapply(object, function(x)
    x$n))
  ord <- NULL
  for (i in 1:num)
    ord <- c(ord, i, i + num)
  
  res <- rbind(varselfreq, varselfreq_tr)[ord, ]
  
  return(round(res, 1))
}

# plot function for estimates and estimand
plot.corr.beta <- function(mat, main = "", ylab = "") {
  plmat <- mat[!grepl(",", rownames(mat)), ]
  plot(
    plmat[, 1],
    type = "l",
    ylim = c(0, 100),
    col = col19[1],
    xaxt = "n",
    ylab = ylab,
    xlab = "",
    main = main
  )
  axis(1, at = 1:nrow(plmat), labels = rownames(plmat))
  for (i in 2:ncol(plmat)) {
    lines(plmat[, i], col = col19[i],
          lty = ifelse(effect[i] == T, 1, 2))
  }
}


# compute estimates and estimand
sub_bwAIC_mat50 <- corr.beta(results_sub50, "bwAIC_betas")
sub_bwAIC_mat63 <- corr.beta(results_sub63, "bwAIC_betas")
sub_bwAIC_mat80 <- corr.beta(results_sub80, "bwAIC_betas")
boot_bwAIC_mat <- corr.beta(results_boot, "bwAIC_betas")

# estimand
bwAIC_tr <- sub_bwAIC_mat50[seq(2, nrow(sub_bwAIC_mat50), 2),]

# take only the estimates in the table
bwAIC_sub50 <-
  sub_bwAIC_mat50[seq(1, nrow(sub_bwAIC_mat50) - 1, 2),][1:3, -1]
bwAIC_sub63 <-
  sub_bwAIC_mat63[seq(1, nrow(sub_bwAIC_mat63) - 1, 2),][1:3, -1]
bwAIC_sub80 <-
  sub_bwAIC_mat80[seq(1, nrow(sub_bwAIC_mat80) - 1, 2),][1:3, -1]
bwAIC_boot <-
  boot_bwAIC_mat[seq(1, nrow(boot_bwAIC_mat) - 1, 2),][1:3, -1]


# preparation for ggplot

rownames(bwAIC_tr) <-
  rownames(bwAIC_sub50) <- rownames(bwAIC_sub63) <-
  rownames(bwAIC_sub80) <-
  rownames(bwAIC_boot) <- c("N = 150", "N = 300", "N = 750")

bwAIC <-
  list(bwAIC_tr, bwAIC_sub50, bwAIC_sub63, bwAIC_sub80, bwAIC_boot)



bwAIC_melt <- melt(bwAIC)
colnames(bwAIC_melt) <- c("N", "Variable", "VIF", "Sampling")
bwAIC_melt$Sampling <- factor(bwAIC_melt$Sampling,
                              labels = c("Estimand", "50%", "63.2%", "80%", "bootstrap"))


bwAIC_estimates_melt <- bwAIC_melt


############# RMSE ##################


# calulate the RMSE
rmse.beta <- function(object = results_sub50, method = "la_betas") {
  
  # extract all betas obtained with resamples
  mat_tr <-
    lapply(object, function(x)
      x[[method]][x[[method]][, "bootstrap"] == 0, names(effect)])
  
  # calculate the VIF for each simulated data set
  varselfreq <-
    lapply(object, function(x)
      aggregate(
        x[[method]][x[[method]][, "bootstrap"] != 0, names(effect)],
        by = list(x[[method]][x[[method]][, "bootstrap"] !=
                                0, "simulation"]),
        FUN = function(z)
          mean(z != 0)
      )[, -1] * 100)
  
  # approximated estimand
  varselfreq_tr <-
    t(sapply(mat_tr, function(y)
      apply(y, 2, function(x)
        mean(x != 0)))) * 100
  
  # calculation of RMSE
  VIF_rmse <- matrix(NA, nrow = length(object), ncol = length(effect))
  for (i in 1:length(object)) {
    VIF_rmse[i, ] <- apply(varselfreq[[i]] - matrix(
      rep(varselfreq_tr[i, ], times = 1000),
      nrow = 1000,
      byrow = T
    ),
    2,
    function(x)
      sqrt(1 / 1000 * sum(x ^ 2)))
  }
  
  rownames(VIF_rmse) <-
    paste("N = ", sapply(object, function(x)
      x$n), sep = "")
  colnames(VIF_rmse) <- names(effect)
  
  
  res <- VIF_rmse
  return(round(res, 3))
}


# calculate the RMSE for the omission/selection strategy (here:benchmark)
rmse.bench <- function(object = results_sub50, method = "bwAIC_betas") {
  
  # extract all betas obtained with the simulated data sets
  mat_tr <-
    lapply(object, function(x)
      x[[method]][x[[method]][, "bootstrap"] == 0, names(effect)])
  
  # set VIF to 0 or 1 depending on omission/selection
  varselfreq <- lapply(mat_tr, function(x)
    (x != 0) * 100)
  
  # approximated estimand 
  varselfreq_tr <-
    t(sapply(varselfreq, function(y)
      apply(y, 2, mean)))
  
  # RMSE for this strategy
  VIF_rmse <- matrix(NA, nrow = length(object), ncol = length(effect))
  for (i in 1:length(object)) {
    VIF_rmse[i, ] <- apply(varselfreq[[i]] - matrix(
      rep(varselfreq_tr[i, ], times = 1000),
      nrow = 1000,
      byrow = T
    ),
    2,
    function(x)
      sqrt(1 / 1000 * sum(x ^ 2)))
  }
  
  rownames(VIF_rmse) <-
    paste("N = ", sapply(object, function(x)
      x$n), sep = "")
  colnames(VIF_rmse) <- names(effect)
  
  
  res <- VIF_rmse
  return(round(res, 3))
}


# preparation for ggplot

sub_bwAIC_mat50 <- rmse.beta(results_sub50, "bwAIC_betas")[, -1]
sub_bwAIC_mat63 <- rmse.beta(results_sub63, "bwAIC_betas")[, -1]
sub_bwAIC_mat80 <- rmse.beta(results_sub80, "bwAIC_betas")[, -1]
boot_bwAIC_mat <- rmse.beta(results_boot, "bwAIC_betas")[, -1]
bench_bwAIC_mat <- rmse.bench(results_sub50, "bwAIC_betas")[, -1]

bwAIC_melt <-
  melt(
    list(
      sub_bwAIC_mat50,
      sub_bwAIC_mat63,
      sub_bwAIC_mat80,
      boot_bwAIC_mat,
      bench_bwAIC_mat
    )
  )




# plot functions for VIF and the RMSE

plot_varsel <- function(mat_melt, n = c(150, 300, 750)) {
  # rank by R2
  mat_melt$Variable <- factor(
    mat_melt$Variable,
    levels =
      c("x11","x5", "x6", "x3", "x8", "x10", "x4a", "x1",  "x2",
        "x4b", "x7", "x9a", "x9b", "x12", "x13",  "x14", "x15")
  )
  mat_melt$Variable <-
    as.numeric(mat_melt$Variable)
  
  bwAIC_plot <- ggplot(data = mat_melt,
                       aes(
                         y = VIF,
                         x = Variable,
                         shape = Sampling,
                         color = Sampling
                       )) +
    geom_hline(
      yintercept = 15.7,
      col = "grey55",
      size = 0.7,
      linetype = "dashed"
    ) +
    scale_color_manual(values = c("black", "blue", "purple", "red", "orange")) +
    scale_shape_manual(values = c(
      "Estimand" = "X",
      "50%" = "5",
      "63.2%" = "6",
      "80%" = "8",
      "bootstrap" = "B"
    )) +
    geom_vline(xintercept = 8.5,
               col = "grey55",
               size = 1) +
    
    geom_point(size = 4)  +
    geom_line() +
    scale_x_continuous(
      name = "",
      breaks = c(4.5, 13),
      labels = c("predictors", "non-predictors")
    ) +
    ylab("VIF (x100)") +
    theme_bw() +
    theme(
      text = element_text(size = 22),
      panel.grid.major.x = element_blank(),
      # remove vertical grid lines
      panel.grid.minor.x = element_blank(),
      # remove vertical grid lines
      legend.position = "none",
      # remove legend
      axis.ticks.x = element_blank()
    ) + # remove x-axis ticks
    ylim(0, 100) +
    facet_grid(N ~ .)
  
  return(bwAIC_plot)
  
}




plot_rmse <- function(mat_melt, n = c(150, 300, 750)) {
  # rank by R2
  mat_melt$Variable <- factor(
    mat_melt$Variable,
    levels =
      c("x11","x5", "x6", "x3", "x8", "x10", "x4a", "x1",  "x2",
        "x4b", "x7", "x9a", "x9b", "x12", "x13",  "x14", "x15")
  )
  mat_melt$Variable <-
    as.numeric(mat_melt$Variable)
  
  mat_melt <- mat_melt[mat_melt$Samplesize %in% paste("N =", n), ]
  
  bwAIC_plot <- ggplot(data = mat_melt,
                       aes(
                         y = RMSEofVIF,
                         x = Variable,
                         shape = Sampling,
                         color = Sampling
                       )) +
    scale_color_manual(values = c("blue", "purple", "red", "orange", "black")) +
    scale_shape_manual(values = c(
      "50%" = "5",
      "63.2%" = "6",
      "80%" = "8",
      "bootstrap" = "B",
      "benchmark" = "O"
    )) +
    geom_vline(xintercept = 8.5,
               col = "grey55",
               size = 1) +
    geom_point(size = 4)  +
    geom_line() +
    scale_x_continuous(
      name = "",
      breaks = c(4.5, 13),
      labels = c("predictors", "non-predictors")
    ) +
    ylab("RMSE of VIF (x100)") +
    theme_bw() +
    theme(
      text = element_text(size = 22),
      panel.grid.major.x = element_blank(),
      # remove vertical grid lines
      panel.grid.minor.x = element_blank(),
      # remove vertical grid lines
      legend.position = "none",
      # remove legend
      axis.ticks.x = element_blank()
    ) + # remove x-axis ticks
    ylim(0, 50) +
    facet_grid(Samplesize ~ .)
  
  return(bwAIC_plot)
  
}


# preparation for ggplot

colnames(bwAIC_melt) <-
  c("Samplesize", "Variable", "RMSEofVIF", "Sampling")
bwAIC_melt$Sampling <- factor(bwAIC_melt$Sampling,
                              labels = c("50%", "63.2%", "80%", "bootstrap", "benchmark"))


colnames(bwAIC_estimates_melt) <-
  c("N", "Variable", "VIF", "Sampling")


# combine both plots

plots <- grid.arrange(
  plot_varsel(mat_melt = bwAIC_estimates_melt) +
    theme(plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm")),
  plot_rmse(mat_melt = bwAIC_melt) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm")),
  nrow = 1
)



pdf("VIF.pdf",
    height = 10,
    width = 16,
    onefile = T)

grid.arrange(plots)

dev.off()
