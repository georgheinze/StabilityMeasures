###############################################################################
##############  Supplementary Material S6-S7: RCB & RMSDR ######################
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
plot(rep(1, 5),
     col = colfunc(5),
     pch = 19,
     cex = 3)

col_cat <- cut(multi_corr, seq(0, 0.75, 0.15))

cols <- colfunc(5)[col_cat]


effect <- setup[1, 3:ncol(setup)] != 0


lwds_cat <- cut(R2[names(effect)[-1]] * 100, c(0.00001, 0.5, 2, 3, 6))
lwds_cat <- ifelse(is.na(lwds_cat), 0, lwds_cat)
lwds  <- recode(lwds_cat, "0=1;1=1;2=1.5; 3=2;4=2.5")


# plot function line graph
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
    #yaxislabel <- paste(type, "of", measure, "with", method)
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
      cex.lab = 2.2,
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



# MD and MAD with bootstrap

boot <-
  list(
    fu = list(),
    bw5 = list(),
    bwAIC = list(),
    la = list()
  )


boot$fu$rmsdratio <-
  calc.dif(
    object = meas_boot,
    method = "fu",
    measure = "rmsdratio",
    true = "rmsdratio_pop"
  )
boot$fu$rcb <-
  calc.dif(
    object = meas_boot,
    method = "fu",
    measure = "rcb",
    true = "rcb_pop"
  )


boot$bw5$rmsdratio <-
  calc.dif(
    object = meas_boot,
    method = "bw5",
    measure = "rmsdratio",
    true = "rmsdratio_pop"
  )
boot$bw5$rcb <-
  calc.dif(
    object = meas_boot,
    method = "bw5",
    measure = "rcb",
    true = "rcb_pop"
  )

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

boot$la$rmsdratio <-
  calc.dif(
    object = meas_boot,
    method = "la",
    measure = "rmsdratio",
    true = "rmsdratio_pop"
  )
boot$la$rcb <-
  calc.dif(
    object = meas_boot,
    method = "la",
    measure = "rcb",
    true = "rcb_pop"
  )





# MD and MAD with subsampling with m=N*0.632

sub <-
  list(
    fu = list(),
    bw5 = list(),
    bwAIC = list(),
    la = list()
  )

# correction of RMSD ratio
for (i in 1:6) {
  for (j in 1:6)
    meas_sub63[[i]][[j]]$rmsdratio <-
      meas_sub63[[i]][[j]]$rmsdratio / sqrt((n[i] - floor(n[i] * 0.632)) / floor(n[i] *
                                                                                 0.632))
}


sub$fu$rmsdratio <-
  calc.dif(
    object = meas_sub63,
    method = "fu",
    measure = "rmsdratio",
    true = "rmsdratio_pop"
  )
sub$fu$rcb <-
  calc.dif(
    object = meas_sub63,
    method = "fu",
    measure = "rcb",
    true = "rcb_pop"
  )


sub$bw5$rmsdratio <-
  calc.dif(
    object = meas_sub63,
    method = "bw5",
    measure = "rmsdratio",
    true = "rmsdratio_pop"
  )
sub$bw5$rcb <-
  calc.dif(
    object = meas_sub63,
    method = "bw5",
    measure = "rcb",
    true = "rcb_pop"
  )

sub$bwAIC$rmsdratio <-
  calc.dif(
    object = meas_sub63,
    method = "bwAIC",
    measure = "rmsdratio",
    true = "rmsdratio_pop"
  )
sub$bwAIC$rcb <-
  calc.dif(
    object = meas_sub63,
    method = "bwAIC",
    measure = "rcb",
    true = "rcb_pop"
  )

sub$la$rmsdratio <-
  calc.dif(
    object = meas_sub63,
    method = "la",
    measure = "rmsdratio",
    true = "rmsdratio_pop"
  )
sub$la$rcb <-
  calc.dif(
    object = meas_sub63,
    method = "la",
    measure = "rcb",
    true = "rcb_pop"
  )





# MD and MAD with subsampling with m=N*0.5

sub50 <-
  list(
    fu = list(),
    bw5 = list(),
    bwAIC = list(),
    la = list(),
    laAIC = list()
  )


#correction of RMSDR
for (i in 1:6) {
  for (j in 1:6)
    meas_sub50[[i]][[j]]$rmsdratio <-
      meas_sub50[[i]][[j]]$rmsdratio / sqrt((n[i] - floor(n[i] * 0.5)) / floor(n[i] *
                                                                                 0.5))
}

sub50$fu$rmsdratio <-
  calc.dif(
    object = meas_sub50,
    method = "fu",
    measure = "rmsdratio",
    true = "rmsdratio_pop"
  )
sub50$fu$rcb <-
  calc.dif(
    object = meas_sub50,
    method = "fu",
    measure = "rcb",
    true = "rcb_pop"
  )



sub50$bw5$rmsdratio <-
  calc.dif(
    object = meas_sub50,
    method = "bw5",
    measure = "rmsdratio",
    true = "rmsdratio_pop"
  )
sub50$bw5$rcb <-
  calc.dif(
    object = meas_sub50,
    method = "bw5",
    measure = "rcb",
    true = "rcb_pop"
  )

sub50$bwAIC$rmsdratio <-
  calc.dif(
    object = meas_sub50,
    method = "bwAIC",
    measure = "rmsdratio",
    true = "rmsdratio_pop"
  )
sub50$bwAIC$rcb <-
  calc.dif(
    object = meas_sub50,
    method = "bwAIC",
    measure = "rcb",
    true = "rcb_pop"
  )

sub50$la$rmsdratio <-
  calc.dif(
    object = meas_sub50,
    method = "la",
    measure = "rmsdratio",
    true = "rmsdratio_pop"
  )
sub50$la$rcb <-
  calc.dif(
    object = meas_sub50,
    method = "la",
    measure = "rcb",
    true = "rcb_pop"
  )





# MD and MAD with subsampling with m=N*0.8

sub80 <-
  list(
    fu = list(),
    bw5 = list(),
    bwAIC = list(),
    la = list()
  )


# correction of RMSD ratio
for (i in 1:6) {
  for (j in 1:6)
    meas_sub80[[i]][[j]]$rmsdratio <-
      meas_sub80[[i]][[j]]$rmsdratio / sqrt((n[i] - floor(n[i] * 0.8)) / floor(n[i] *
                                                                                 0.8))
}


sub80$fu$rmsdratio <-
  calc.dif(
    object = meas_sub80,
    method = "fu",
    measure = "rmsdratio",
    true = "rmsdratio_pop"
  )
sub80$fu$rcb <-
  calc.dif(
    object = meas_sub80,
    method = "fu",
    measure = "rcb",
    true = "rcb_pop"
  )

sub80$bw5$rmsdratio <-
  calc.dif(
    object = meas_sub80,
    method = "bw5",
    measure = "rmsdratio",
    true = "rmsdratio_pop"
  )
sub80$bw5$rcb <-
  calc.dif(
    object = meas_sub80,
    method = "bw5",
    measure = "rcb",
    true = "rcb_pop"
  )

sub80$bwAIC$rmsdratio <-
  calc.dif(
    object = meas_sub80,
    method = "bwAIC",
    measure = "rmsdratio",
    true = "rmsdratio_pop"
  )
sub80$bwAIC$rcb <-
  calc.dif(
    object = meas_sub80,
    method = "bwAIC",
    measure = "rcb",
    true = "rcb_pop"
  )

sub80$la$rmsdratio <-
  calc.dif(
    object = meas_sub80,
    method = "la",
    measure = "rmsdratio",
    true = "rmsdratio_pop"
  )
sub80$la$rcb <-
  calc.dif(
    object = meas_sub80,
    method = "la",
    measure = "rcb",
    true = "rcb_pop"
  )








# plot RCB for BE(AIC), BE(0.05) and Lasso
# suppplementary figure 5B

pdf("rcb_all.pdf", height = 10, width = 15)

par(mfrow = c(3, 3), mar = c(5, 7, 2, 2))

j = "bwAIC"
plot.res(
  object = boot,
  method = j,
  measure = "rcb",
  type = "median",
  ylim = c(-50, 200),
  line_at = 0,
  xlab = "",
  main = "BE(AIC)",
  ylab = "Median \n of estimated RCB (x100)"
)
j = "bw5"
plot.res(
  object = boot,
  method = j,
  measure = "rcb",
  type = "median",
  ylim = c(-50, 200),
  line_at = 0,
  xlab = "",
  main = "BE(0.05)",
  ylab = ""
)
j = "la"
plot.res(
  object = boot,
  method = j,
  measure = "rcb",
  type = "median",
  ylim = c(-50, 200),
  line_at = 0,
  xlab = "",
  main = "Lasso",
  ylab = ""
)


j = "bwAIC"
plot.res(
  object = boot,
  method = j,
  measure = "rcb",
  type = "md",
  ylim = c(-25, 10),
  line_at = 0,
  xlab = "",
  main = "",
  ylab = "Median deviation \n of estimated RCB (x100)"
)
j = "bw5"
plot.res(
  object = boot,
  method = j,
  measure = "rcb",
  type = "md",
  ylim = c(-25, 10),
  line_at = 0,
  xlab = "",
  main = "",
  ylab = ""
)
j = "la"
plot.res(
  object = boot,
  method = j,
  measure = "rcb",
  type = "md",
  ylim = c(-25, 10),
  line_at = 0,
  xlab = "",
  main = "",
  ylab = ""
)



j = "bwAIC"
plot.res(
  object = boot,
  method = j,
  measure = "rcb",
  type = "mad",
  ylim = c(0, 70),
  line_at = 0,
  xlab = "N",
  main = "",
  ylab = "Median absolute deviation \n of estimated RCB (x100)"
)
j = "bw5"
plot.res(
  object = boot,
  method = j,
  measure = "rcb",
  type = "mad",
  ylim = c(0, 70),
  line_at = 0,
  xlab = "N",
  main = "",
  ylab = ""
)
j = "la"
plot.res(
  object = boot,
  method = j,
  measure = "rcb",
  type = "mad",
  ylim = c(0, 70),
  line_at = 0,
  xlab = "N",
  main = "",
  ylab = ""
)

dev.off()




# plot RMSDR for BE(AIC), BE(0.05) and Lasso
# suppplementary figure 6B

pdf("rmsdratio_all.pdf", height = 10, width = 15)
par(mfrow = c(3, 3), mar = c(5, 7, 2, 2))


j = "bwAIC"
plot.res(
  object = boot,
  method = j,
  measure = "rmsdratio",
  type = "median",
  ylim = c(0.8, 1.4),
  line_at = 1,
  xlab = "",
  main = "BE(AIC)",
  ylab = "Median \n of estimated RMSD ratio"
)


j = "bw5"
plot.res(
  object = boot,
  method = j,
  measure = "rmsdratio",
  type = "median",
  ylim = c(0.8, 1.4),
  line_at = 1,
  xlab = "",
  main = "BE(0.05)",
  ylab = ""
)

j = "la"
plot.res(
  object = boot,
  method = j,
  measure = "rmsdratio",
  type = "median",
  ylim = c(0.8, 1.4),
  line_at = 1,
  xlab = "",
  main = "Lasso",
  ylab = ""
)


j = "bwAIC"
plot.res(
  object = boot,
  method = j,
  measure = "rmsdratio",
  type = "md",
  ylim = c(-0.2, 0.5),
  line_at = 0,
  xlab = "",
  main = "",
  ylab = "Median deviation \n of estimated RMSD ratio"
)


j = "bw5"
plot.res(
  object = boot,
  method = j,
  measure = "rmsdratio",
  type = "md",
  ylim = c(-0.2, 0.5),
  line_at = 0,
  xlab = "",
  ylab = "",
  main = ""
)

j = "la"
plot.res(
  object = boot,
  method = j,
  measure = "rmsdratio",
  type = "md",
  ylim = c(-0.2, 0.5),
  line_at = 0,
  xlab = "",
  ylab = "",
  main = ""
)



j = "bwAIC"
plot.res(
  object = boot,
  method = j,
  measure = "rmsdratio",
  type = "mad",
  ylim = c(0, 0.5),
  line_at = 0,
  xlab = "N",
  main = "",
  ylab = "Median absolute deviation \n of estimated RMSD ratio"
)


j = "bw5"
plot.res(
  object = boot,
  method = j,
  measure = "rmsdratio",
  type = "mad",
  ylim = c(0, 0.5),
  line_at = 0,
  xlab = "N",
  ylab = "",
  main = ""
)

j = "la"
plot.res(
  object = boot,
  method = j,
  measure = "rmsdratio",
  type = "mad",
  ylim = c(0, 0.5),
  line_at = 0,
  xlab = "N",
  ylab = "",
  main = ""
)

dev.off()





# function for plotting RCB and RMSDR for S0.5, S0.632, S0.8 and the bootstrap
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
    item <- object[[method]][[measure]][[type]][, 1:3]
    #yaxislabel <- paste(type, "of", measure, "with", method)
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
      cex.lab = 2.2,
      cex.axis = 2
    )
    axis(
      1,
      at = 1:ncol(item),
      labels = c("150", "300", "750"),
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




# plot RCB for S0.5, S0.632, S0.8 and the bootstrap
# suppplementary figure 5A

pdf("rcb_m_bwAIC.pdf", height = 10, width = 18)

par(mfrow = c(3, 4), mar = c(5, 7, 2, 2))

j = "bwAIC"
plot.res(
  object = sub50,
  method = j,
  measure = "rcb",
  type = "median",
  ylim = c(-50, 200),
  line_at = 0,
  xlab = "",
  main = expression('S'[0.5]),
  ylab = "Median \n of estimated RCB (x100)"
)

plot.res(
  object = sub,
  method = j,
  measure = "rcb",
  type = "median",
  ylim = c(-50, 200),
  line_at = 0,
  xlab = "",
  main = expression('S'[0.632]),
  ylab = ""
)

plot.res(
  object = sub80,
  method = j,
  measure = "rcb",
  type = "median",
  ylim = c(-50, 200),
  line_at = 0,
  xlab = "",
  main = expression('S'[0.8]),
  ylab = ""
)

plot.res(
  object = boot,
  method = j,
  measure = "rcb",
  type = "median",
  ylim = c(-50, 200),
  line_at = 0,
  xlab = "",
  main = "B",
  ylab = ""
)


ylim = c(-40, 20)

plot.res(
  object = sub50,
  method = j,
  measure = "rcb",
  type = "md",
  ylim = ylim,
  line_at = 0,
  xlab = "",
  main = "",
  ylab = "Median deviation \n of estimated RCB (x100)"
)

plot.res(
  object = sub,
  method = j,
  measure = "rcb",
  type = "md",
  ylim = ylim,
  line_at = 0,
  xlab = "",
  main = "",
  ylab = ""
)

plot.res(
  object = sub80,
  method = j,
  measure = "rcb",
  type = "md",
  ylim = ylim,
  line_at = 0,
  xlab = "",
  main = "",
  ylab = ""
)

plot.res(
  object = boot,
  method = j,
  measure = "rcb",
  type = "md",
  ylim = ylim,
  line_at = 0,
  xlab = "",
  main = "",
  ylab = ""
)



ylim = c(0, 80)

plot.res(
  object = sub50,
  method = j,
  measure = "rcb",
  type = "mad",
  ylim = ylim,
  line_at = 0,
  xlab = "N",
  main = "",
  ylab = "Median absolute deviation \n of estimated RCB (x100)"
)

plot.res(
  object = sub,
  method = j,
  measure = "rcb",
  type = "mad",
  ylim = ylim,
  line_at = 0,
  xlab = "N",
  main = "",
  ylab = ""
)

plot.res(
  object = sub80,
  method = j,
  measure = "rcb",
  type = "mad",
  ylim = ylim,
  line_at = 0,
  xlab = "N",
  main = "",
  ylab = ""
)

plot.res(
  object = boot,
  method = j,
  measure = "rcb",
  type = "mad",
  ylim = ylim,
  line_at = 0,
  xlab = "N",
  main = "",
  ylab = ""
)

dev.off()




# plot RMSDR for S0.5, S0.632, S0.8 and the bootstrap
# suppplementary figure 6A

pdf("rmsdratio_m_full.pdf", height = 10, width = 18)
par(mfrow = c(3, 4), mar = c(5, 7, 2, 2))
j = "fu"

plot.res(
  object = sub50,
  method = j,
  measure = "rmsdratio",
  type = "median",
  ylim = c(0.95, 1.3),
  line_at = 1,
  xlab = "",
  ylab = "Median of \n estimated RMSD ratio",
  main = expression('S'[0.5])
)
plot.res(
  object = sub,
  method = j,
  measure = "rmsdratio",
  type = "median",
  ylim = c(0.95, 1.3),
  line_at = 1,
  xlab = "",
  ylab = "",
  main = expression('S'[0.632])
)
plot.res(
  object = sub80,
  method = j,
  measure = "rmsdratio",
  type = "median",
  ylim = c(0.95, 1.3),
  line_at = 1,
  xlab = "",
  ylab = "",
  main = expression('S'[0.8])
)
plot.res(
  object = boot,
  method = j,
  measure = "rmsdratio",
  type = "median",
  ylim = c(0.95, 1.3),
  line_at = 1,
  xlab = "",
  ylab = "",
  main = "B"
)


ylim = c(0, 0.3)

plot.res(
  object = sub50,
  method = j,
  measure = "rmsdratio",
  type = "md",
  ylim = ylim,
  line_at = 0,
  xlab = "",
  main = "",
  ylab = "Median deviation \n of estimated RMSD ratio"
)

plot.res(
  object = sub,
  method = j,
  measure = "rmsdratio",
  type = "md",
  ylim = ylim,
  line_at = 0,
  xlab = "",
  main = "",
  ylab = ""
)

plot.res(
  object = sub80,
  method = j,
  measure = "rmsdratio",
  type = "md",
  ylim = ylim,
  line_at = 0,
  xlab = "",
  main = "",
  ylab = ""
)

plot.res(
  object = boot,
  method = j,
  measure = "rmsdratio",
  type = "md",
  ylim = ylim,
  line_at = 0,
  xlab = "",
  main = "",
  ylab = ""
)



ylim = c(0, 0.3)

plot.res(
  object = sub50,
  method = j,
  measure = "rmsdratio",
  type = "mad",
  ylim = ylim,
  line_at = 0,
  xlab = "N",
  main = "",
  ylab = "Median absolute deviation \n of estimated RMSD ratio"
)

plot.res(
  object = sub,
  method = j,
  measure = "rmsdratio",
  type = "mad",
  ylim = ylim,
  line_at = 0,
  xlab = "N",
  main = "",
  ylab = ""
)

plot.res(
  object = sub80,
  method = j,
  measure = "rmsdratio",
  type = "mad",
  ylim = ylim,
  line_at = 0,
  xlab = "N",
  main = "",
  ylab = ""
)

plot.res(
  object = boot,
  method = j,
  measure = "rmsdratio",
  type = "mad",
  ylim = ylim,
  line_at = 0,
  xlab = "N",
  main = "",
  ylab = ""
)


dev.off()

