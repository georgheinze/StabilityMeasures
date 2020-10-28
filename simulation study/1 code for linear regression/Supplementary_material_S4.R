
###############################################################################################
##############  Supplementary Material S4: VIF  ###############################################
###############################################################################################

library(reshape)
library(ggplot2)
library(gridExtra)
library(grid)

source("function_beta.R")
load("sim_results_beta.RData")
load("setup.RData")

effect <- setup[1, 3:ncol(setup)] != 0


# function to calculate VIF and approximated estimands
corr.beta <- function(object=results_boot, method="la_betas"){
  
  mat <- lapply(object, function(x) x[[method]][x[[method]][,"bootstrap"]!=0,names(effect)])
  mat_tr <- lapply(object, function(x) x[[method]][x[[method]][,"bootstrap"]==0,names(effect)]) 
  
  varselfreq <- t(sapply(mat, function(y) apply(y,2, function(x) mean(x!=0) )))*100
  varselfreq_tr <- t(sapply(mat_tr, function(y) apply(y,2, function(x) mean(x!=0) )))*100
  
  rownames(varselfreq) <- paste("n=",lapply(object, function(x) x$n), sep="")
  rownames(varselfreq_tr) <- paste("n=",lapply(object, function(x) x$n),", true", sep="")
  
  num<-length(lapply(object, function(x) x$n))
  ord <- NULL
  for(i in 1:num) ord <- c(ord,i, i+num)
  
  res <- rbind(varselfreq, varselfreq_tr)[ord,]
  
  return(round(res,1))
}


# function for plotting estimates and estimand
plot.corr.beta <- function(mat, main="", ylab=""){
  plmat <- mat[!grepl(",",rownames(mat)),]
  plot(plmat[,1], type="l", ylim=c(0,100), col=col19[1],
       xaxt="n", ylab=ylab,
       xlab="", main=main)
  axis(1, at=1:nrow(plmat), labels=rownames(plmat))
  for(i in 2:ncol(plmat)){
    lines(plmat[,i], col=col19[i],
          lty=ifelse(effect[i]==T, 1,2))
  }
}


# calculate the measures for all selection procedures
sub_la_mat50 <- corr.beta(results_sub50, "la_betas")
sub_bw5_mat50 <- corr.beta(results_sub50, "bw5_betas")
sub_bwAIC_mat50 <- corr.beta(results_sub50, "bwAIC_betas")

sub_la_mat63 <- corr.beta(results_sub63, "la_betas")
sub_bw5_mat63 <- corr.beta(results_sub63, "bw5_betas")
sub_bwAIC_mat63 <- corr.beta(results_sub63, "bwAIC_betas")

sub_la_mat80 <- corr.beta(results_sub80, "la_betas")
sub_bw5_mat80 <- corr.beta(results_sub80, "bw5_betas")
sub_bwAIC_mat80 <- corr.beta(results_sub80, "bwAIC_betas")

boot_la_mat <- corr.beta(results_boot, "la_betas")
boot_bw5_mat <- corr.beta(results_boot, "bw5_betas")
boot_bwAIC_mat <- corr.beta(results_boot, "bwAIC_betas")



# plot function for line graph
plot_varsel <- function(mat_melt, n=c(150, 300, 750) ){
  
  
  # rank by R2
  mat_melt$Variable <- factor(mat_melt$Variable, levels=
                                c("x11", "x5", "x6", "x3",  "x8", "x10","x4a", "x1",
                                  "x2", "x4b", "x7", "x9a", "x9b", "x12", "x13", "x14", "x15"))
  mat_melt$Variable <- as.numeric(mat_melt$Variable)                         
  
  bwAIC_plot<- ggplot(data=mat_melt,
                      aes(y=VIF, 
                          x=Variable,
                          shape=Sampling,
                          color=Sampling)) + 
    scale_color_manual(values=c("black", "blue", "purple", "red", "orange")) +
    scale_shape_manual(values=c("Estimand"="X","50%"="5", 
                                "63.2%"="6", "80%"="8", "bootstrap"="B"))+
    geom_vline(xintercept=8.5,
               col="grey55", size=1)+
    geom_point(size=4)  + 
    geom_line()+
    scale_x_continuous(name="", breaks=c(4.5,13),
                       labels=c("predictors", "non-predictors")) + 
    ylab("VIF (x100)") + 
    theme_bw()+
    theme(text = element_text(size=22),
          panel.grid.major.x = element_blank(), # remove vertical grid lines
          panel.grid.minor.x = element_blank(), # remove vertical grid lines
          legend.position="none", # remove legend
          axis.ticks.x = element_blank()) + # remove x-axis ticks 
    ylim(0,100) +
    facet_grid(N ~ VSmethod)
  
  print(bwAIC_plot)
  
}




# preparation for ggplot

## BW(AIC)
bwAIC_tr <- sub_bwAIC_mat50[seq(2, nrow(sub_bwAIC_mat50),2), ]
  
  
bwAIC_sub50 <- sub_bwAIC_mat50[seq(1, nrow(sub_bwAIC_mat50)-1,2), ][-c(4:6),-1]
bwAIC_sub63 <- sub_bwAIC_mat63[seq(1, nrow(sub_bwAIC_mat63)-1,2), ][-c(4:6),-1]
bwAIC_sub80 <- sub_bwAIC_mat80[seq(1, nrow(sub_bwAIC_mat80)-1,2), ][-c(4:6),-1]
bwAIC_boot <- boot_bwAIC_mat[seq(1, nrow(boot_bwAIC_mat)-1,2), ][-c(4:6),-1]

rownames(bwAIC_tr) <- rownames(bwAIC_sub50) <- rownames(bwAIC_sub63) <- 
  rownames(bwAIC_sub80) <- rownames(bwAIC_boot) <- c( "N = 150", "N = 300", "N = 750")

bwAIC<- list(bwAIC_tr, bwAIC_sub50, bwAIC_sub63, bwAIC_sub80, bwAIC_boot)



bwAIC_melt <- melt(bwAIC)
colnames(bwAIC_melt) <- c( "N", "Variable","VIF", "Sampling")
bwAIC_melt$Sampling <- factor(bwAIC_melt$Sampling, 
                              labels= c("Estimand", "50%", "63.2%", "80%", "bootstrap"))



## BW(0.05)
bw5_tr <- sub_bw5_mat50[seq(2, nrow(sub_bw5_mat50),2), ]


bw5_sub50 <- sub_bw5_mat50[seq(1, nrow(sub_bw5_mat50)-1,2), ][-c(4:6),-1]
bw5_sub63 <- sub_bw5_mat63[seq(1, nrow(sub_bw5_mat63)-1,2), ][-c(4:6),-1]
bw5_sub80 <- sub_bw5_mat80[seq(1, nrow(sub_bw5_mat80)-1,2), ][-c(4:6),-1]
bw5_boot <- boot_bw5_mat[seq(1, nrow(boot_bw5_mat)-1,2), ][-c(4:6),-1]

rownames(bw5_tr) <- rownames(bw5_sub50) <- rownames(bw5_sub63) <- 
  rownames(bw5_sub80) <- rownames(bw5_boot) <- c( "N = 150", "N = 300", "N = 750")

bw5<- list(bw5_tr, bw5_sub50, bw5_sub63, bw5_sub80, bw5_boot)



bw5_melt <- melt(bw5)
colnames(bw5_melt) <- c( "N", "Variable","VIF", "Sampling")
bw5_melt$Sampling <- factor(bw5_melt$Sampling, 
                              labels= c("Estimand", "50%", "63.2%", "80%", "bootstrap"))







## Lasso

la_tr <- sub_la_mat50[seq(2, nrow(sub_la_mat50),2), ][,-1]

la_sub50 <- sub_la_mat50[seq(1, nrow(sub_la_mat50)-1,2), ][-c(4:6),-1]
la_sub63 <- sub_la_mat63[seq(1, nrow(sub_la_mat63)-1,2), ][-c(4:6),-1]
la_sub80 <- sub_la_mat80[seq(1, nrow(sub_la_mat80)-1,2), ][-c(4:6),-1]
la_boot <- boot_la_mat[seq(1, nrow(boot_la_mat)-1,2), ][-c(4:6),-1]

rownames(la_tr) <- rownames(la_sub50) <- rownames(la_sub63) <- 
  rownames(la_sub80) <- rownames(la_boot) <- c( "N = 150", "N = 300", "N = 750")

la<- list(la_tr, la_sub50, la_sub63, la_sub80, la_boot)



la_melt <- melt(la)
colnames(la_melt) <- c( "N", "Variable","VIF",  "Sampling")
la_melt$Sampling <- factor(la_melt$Sampling, 
                              labels= c("Estimand", "50%", "63.2%", "80%", "bootstrap"))




mat_melt <- rbind(
  cbind(bwAIC_melt, "VSmethod"= rep("BE(AIC)", nrow(bwAIC_melt))),
  cbind(bw5_melt, "VSmethod"= rep("BE(0.05)", nrow(bw5_melt))),
  cbind(la_melt, "VSmethod"= rep("Lasso", nrow(la_melt))))



#save plot 
pdf("varsel_all.pdf", height=12, width=18, onefile=T)

plot_varsel(mat_melt=mat_melt)

dev.off()


bwAIC_estimates_melt <- bwAIC_melt





rmse.beta <- function(object=results_sub50, method="la_betas"){
  
  mat_tr <- lapply(object, function(x) x[[method]][x[[method]][,"bootstrap"]==0,names(effect)]) 
  
  varselfreq <-lapply(object, function(x) aggregate(x[[method]][x[[method]][,"bootstrap"]!=0,names(effect)],
                                       by=list(x[[method]][x[[method]][,"bootstrap"]!=0, "simulation"]), 
                                            FUN= function(z) mean(z!=0))[,-1]*100)  
  
  varselfreq_tr <- t(sapply(mat_tr, function(y) apply(y,2, function(x) mean(x!=0) )))*100
  
  VIF_rmse <- matrix(NA, nrow=length(object), ncol=length(effect))
  for(i in 1:length(object)){
    VIF_rmse[i,] <- apply(
        varselfreq[[i]] - matrix(rep(varselfreq_tr[i,], times=1000), nrow=1000, byrow=T),
        2,
        function(x) sqrt(1/1000*sum(x^2)))
  }

  rownames(VIF_rmse) <- paste("N = ",sapply(object, function(x) x$n), sep="")
  colnames(VIF_rmse) <- names(effect)
  
  res <- VIF_rmse
  return(round(res,3))
}



rmse.bench <- function(object=results_sub50, method="bwAIC_betas"){
  
  mat_tr <- lapply(object, function(x) x[[method]][x[[method]][,"bootstrap"]==0,names(effect)]) 
  
  varselfreq <- lapply(mat_tr, function(x) (x != 0) * 100) 
  
  
  varselfreq_tr <- t(sapply(varselfreq, function(y) apply(y,2, mean)))
  
  VIF_rmse <- matrix(NA, nrow=length(object), ncol=length(effect))
  for(i in 1:length(object)){
    VIF_rmse[i,] <- apply(
      varselfreq[[i]] - matrix(rep(varselfreq_tr[i,], times=1000), nrow=1000, byrow=T),
      2,
      function(x) sqrt(1/1000*sum(x^2)))
  }
  
  rownames(VIF_rmse) <- paste("N = ",sapply(object, function(x) x$n), sep="")
  colnames(VIF_rmse) <- names(effect)

  
  res <- VIF_rmse
  return(round(res,3))
}



sub_la_mat50 <- rmse.beta(results_sub50, "la_betas")[,-1]
sub_bw5_mat50 <- rmse.beta(results_sub50, "bw5_betas")[,-1]
sub_bwAIC_mat50 <- rmse.beta(results_sub50, "bwAIC_betas")[,-1]

sub_la_mat63 <- rmse.beta(results_sub63, "la_betas")[,-1]
sub_bw5_mat63 <- rmse.beta(results_sub63, "bw5_betas")[,-1]
sub_bwAIC_mat63 <- rmse.beta(results_sub63, "bwAIC_betas")[,-1]

sub_la_mat80 <- rmse.beta(results_sub80, "la_betas")[,-1]
sub_bw5_mat80 <- rmse.beta(results_sub80, "bw5_betas")[,-1]
sub_bwAIC_mat80 <- rmse.beta(results_sub80, "bwAIC_betas")[,-1]


boot_la_mat <- rmse.beta(results_boot, "la_betas")[,-1]
boot_bw5_mat <- rmse.beta(results_boot, "bw5_betas")[,-1]
boot_bwAIC_mat <- rmse.beta(results_boot, "bwAIC_betas")[,-1]

bench_la_mat <- rmse.bench(results_sub50, "la_betas")[,-1]
bench_bw5_mat <- rmse.bench(results_sub50, "bw5_betas")[,-1]
bench_bwAIC_mat <- rmse.bench(results_sub50, "bwAIC_betas")[,-1]



bwAIC_melt <- melt(list(sub_bwAIC_mat50, sub_bwAIC_mat63, sub_bwAIC_mat80, 
                        boot_bwAIC_mat, bench_bwAIC_mat))
bw5_melt <- melt(list(sub_bw5_mat50, sub_bw5_mat63, sub_bw5_mat80, 
                        boot_bw5_mat, bench_bw5_mat))
la_melt <- melt(list(sub_la_mat50, sub_la_mat63, sub_la_mat80, 
                     boot_la_mat, bench_la_mat))

mat_melt <- cbind(rbind(bwAIC_melt,bw5_melt, la_melt), 
                  "Method" = rep(c("BW(AIC)","BW(0.05)", "Lasso"), each=nrow(la_melt)))

colnames(mat_melt) <- c("Samplesize", "Variable", "RMSEofVIF", "Sampling", "Method")

mat_melt$Method  <- relevel(mat_melt$Method, "BW(AIC)")

mat_melt$Sampling <- factor(mat_melt$Sampling, 
                              labels= c("50%", "63.2%", "80%", "bootstrap", "benchmark"))


# plot function
plot_rmse <- function(mat_melt, n=c(150, 300, 750)){

  
  # rank by R2
  mat_melt$Variable <- factor(mat_melt$Variable, levels=
                                c("x11", "x5", "x6", "x3",  "x8", "x10","x4a", "x1",
                                  "x2", "x4b", "x7", "x9a", "x9b", "x12", "x13", "x14", "x15"))
  mat_melt$Variable <- as.numeric(mat_melt$Variable)                         
  
  mat_melt <- mat_melt[mat_melt$Samplesize %in% paste("N =", n),]
  
  bwAIC_plot <- ggplot(data=mat_melt,
                      aes(y=RMSEofVIF, 
                          x=Variable,
                          #group=Samplesize,
                          shape=Sampling,
                          color=Sampling)) + 
    scale_color_manual(values=c( "blue", "purple", "red", "orange", "black")) +
    scale_shape_manual(values=c("50%"="5", 
                                "63.2%"="6", "80%"="8",
                                "bootstrap"="B",
                                "benchmark"="O")) +
  
    geom_vline(xintercept=8.5,
               col="grey55", size=1)+
    geom_point(size=4)  +
    geom_line()+

    scale_x_continuous(name="", breaks=c(4.5,13),
                       labels=c("predictors", "non-predictors")) + 
    ylab("RMSE of VIF (x100)") + 
    theme_bw()+
    theme(text = element_text(size=22),
          panel.grid.major.x = element_blank(), # remove vertical grid lines
          panel.grid.minor.x = element_blank(), # remove vertical grid lines
          legend.position="none", # remove legend
          axis.ticks.x = element_blank()) + # remove x-axis ticks 
    ylim(0,50)  + 
    facet_grid( Samplesize ~ Method)
  
  print(bwAIC_plot)
  
}


pdf("varsel_rmse_benchmark.pdf", height=12, width=18, onefile=T)

plot_rmse(mat_melt=mat_melt)

dev.off()


