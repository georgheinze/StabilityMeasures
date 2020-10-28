###############################################################################
##############  Table 2: Coefficients for the simulation setup  ################
############# + Supplementary Material S2: Histograms of variables #############
###############################################################################


load("setup.RData")


# order betas by their partial R2

names(R2)[order(R2, decreasing=T)]

ran <- c(
  "x11",
  "x5",
  "x6",
  "x3",
  "x8",
  "x10",
  "x4a",
  "x1",
  "x2",
  "x4b",
  "x7",
  "x9a",
  "x9b",
  "x12",
  "x13",
  "x14",
  "x15"
)

# Table 2

cbind(
  "coefficients"=beta[ran], 
  "Sd coefficients"= (beta[-1] * x.sd[-1])[ran],
  "partial R2"= R2[ran],
  "multiple R2"=multi_corr[ran]
)


# Supplementary Material S2
## Plot of histograms

set.seed(8643)
take <- sample(1:50000, 1000)

dat_cor<- dat[take,ran]
colnames(dat_cor) <- paste0("X", 1:17)


pdf("histograms.pdf",height=10, width=15)
par(mfrow=c(4,5), mar=c(5,5,2,2))
for(i in 1:17) hist(dat_cor[,i], xlab=paste0("x",i), main="", probability=T,
                    ylab=ifelse(i%in% c(1,6,11,16), "Density", ""),
                    cex.lab=2)
dev.off()



