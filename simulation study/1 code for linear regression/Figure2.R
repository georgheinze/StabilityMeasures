###############################################################################
##############  Figure 2: RCB & RMSDR for the bootstrap ########################
###############################################################################
library(car)

load("setup.RData")
load("sim_measures.RData")


# function to calculate median deviation and median absolute deviation
calc.dif <- function(object, method, measure, true) {
  median_m <- sapply(object, function(y)
    apply(y[[method]][[measure]], 2, function(x)
      median(x, na.rm = T)))
  
  # approximated estimand
  pop_value <-
    sapply(object, function(y)
      simplify2array(y[[method]][[true]]))
  
  # difference between estimate and estimand
  dif <-
    lapply(object, function(y)
      y[[method]][[measure]] - matrix(rep(simplify2array(y[[method]][[true]]), each =
                                            1000), nrow = 1000))
  # median deviation
  md <-
    sapply(dif, function(x)
      apply(x, 2, function(x)
        median(x, na.rm = T)))
  
  # median absolute deviation
  mad <-
    sapply(dif, function(x)
      apply(x, 2, function(y)
        median(abs(y), na.rm = T)))
  
  md[is.infinite(md)] <- NA
  mad[is.infinite(mad)] <- NA
  
  rownames(median_m) <-
    rownames(mad) <- rownames(md) <- colnames(setup[, 3:ncol(setup)])
  
  return(list(
    "median" = median_m,
    "md" = md,
    "mad" = mad,
    "pop_value" = pop_value
  ))
}


# determine colors for lines by their multiple correlation
colfunc <- colorRampPalette(c("blue", "red"))
colfunc(5)

col_cat <- cut(multi_corr, seq(0, 0.75, 0.15))
cols <- colfunc(5)[col_cat]
effect <- setup[1, 3:ncol(setup)] != 0


lwds_cat <- cut(R2[names(effect)[-1]] * 100, c(0.00001, 0.5, 2, 3, 6))
lwds_cat <- ifelse(is.na(lwds_cat), 0, lwds_cat)
lwds  <- recode(lwds_cat, "0=1;1=1;2=1.5; 3=2;4=2.5")


# plot line graph
plot.res <-
  function(object = boot,
           method = "bw5",
           measure = "rmsdratio",
           type = "median",
           ylim = NULL,
           line_at = NULL,
           xlab = "",
           ylab = "",
           main = "") {
    item <- object[[method]][[measure]][[type]]
    if (is.null(ylim))
      ylim <- c(min(item, na.rm = T), max(item, na.rm = T))
    plot(
      item[2, 1:ncol(item)],
      type = "l",
      xaxt = 'n',
      xlab = xlab,
      ylab = ylab,
      ylim = ylim,
      font.lab = 2,
      main = main,
      col = cols[1],
      lwd = lwds[1],
      cex.main = 2,
      cex.lab = 2,
      cex.axis = 2
    )
    axis(
      1,
      at = 1:ncol(item),
      labels = c("150", "300", "750", "1000", "5000", "10000"),
      cex.axis = 1.8
    )
    if (!is.null(line_at))
      abline(h = line_at, col = "grey65")
    
    for (i in 3:nrow(item)) {
      lines(
        item[i, 1:ncol(item)],
        col = cols[i - 1],
        lty = ifelse(effect[i] == T, 1, 2),
        lwd = lwds[i - 1]
      )
    }
  }




# bootstrap results

boot <- list(bwAIC = list())

boot$bwAIC$rmsdratio <-
  calc.dif(
    object = meas_boot,
    method = "bwAIC",
    measure = "rmsdratio",
    true = "rmsdratio_pop"
  )
boot$bwAIC$rcb <-
  calc.dif(
    object = meas_boot,
    method = "bwAIC",
    measure = "rcb",
    true = "rcb_pop"
  )



# produce a pdf for Figure 2

pdf("rcb_rmsdratio.pdf", height = 8.5, width = 12.5)

par(mfrow = c(2, 2), mar = c(5, 7, 2, 2))
j = "bwAIC"


plot.res(
  object = boot,
  method = j,
  measure = "rcb",
  type = "median",
  ylim = c(-10, 150),
  line_at = 0,
  xlab = "",
  main = "",
  ylab = "Median of \n estimated RCB (x100)"
)

plot.res(
  object = boot,
  method = j,
  measure = "rmsdratio",
  type = "median",
  ylim = c(0.9, 1.2),
  line_at = 1,
  xlab = "",
  main = "",
  ylab = "Median of \n estimated RMSD ratio"
)

plot.res(
  object = boot,
  method = j,
  measure = "rcb",
  type = "md",
  ylim = c(-15, 10),
  line_at = 0,
  xlab = "N",
  main = "",
  ylab = "Median deviation \n of estimated RCB (x100)"
)

plot.res(
  object = boot,
  method = j,
  measure = "rmsdratio",
  type = "md",
  ylim = c(-0.05, 0.3),
  line_at = 0,
  xlab = "N",
  ylab = "Median deviation \n of estimated RMSD ratio",
  main = ""
)



dev.off()
