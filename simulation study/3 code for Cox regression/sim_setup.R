#####################################################################
##############  genereate and save simulation setup ################
#####################################################################


# note: the code for this simulation study ran on the Vienna Scientific Cluster. Results may be different if the code is executed on a single PC.


library(mvtnorm)
library(rms)
library(data.table)


# generate a large data set 
set.seed(8586216)

source("function_setup.R")

z<-gen.z(n=50000)
x<-gen.x.from.z(z=z)

pred<-partial.predictor(x)
ystar<-apply(pred,1,sum)

x.t<-pretrans.x(x)

round(cor(x.t)*100,0)


dat<-x.to.data(x=x.t, y=ystar)

x.sd<-apply(dat,2,sd)
scope<-colnames(dat)[-1]

formr<-paste(scope,collapse="+")


leastfalse<-lm(ystar~x1+x3+x4a+x5+x6+x8+x10+x11, data=dat)

# true beta
beta<-(1:18)*0
names(beta)<-c("(Intercept)",colnames(dat)[2:ncol(dat)])
beta[names(coef(leastfalse))]<-coef(leastfalse)


# generate y

beta<- beta[-1]
linpred<-as.matrix(dat[,2:ncol(dat)]) %*% beta - 3.5


y <-(-log(runif(length(linpred)))/(1/5*exp(linpred)))**(1/3)
cens <- rep(1, length=length(linpred))
tau <- 3.85
futime <- 0.001+runif(length(linpred))*tau
cens[y > futime] <- 0
1-mean(cens)

y[cens==0] <- futime[cens==0] 
plot(survfit(Surv(y, cens)~1))


dat$y<-y

dat <- cbind(dat, cens)


global<-as.formula(paste("Surv(y,cens)~",formr))



# determine R2

source("function_schemper.R")

rand <- sample(1:50000, 50000)
mod <- coxph(global, data=dat[rand,])
lp <- predict(mod, dat[rand,], "lp")
full.r2 <- inaccuracy_surv(model2= mod, lp2=lp, data=dat[rand,], time_name="y", event_name="cens")[[2]]
full.r2


drop.r2<-function(scope=scope, y=y, data=dat, dr="x1"){
  which(scope!=dr)
  fr<-as.formula(paste("Surv(y,cens)~",paste(scope[which(scope!=dr)],collapse="+")))
  fr.fit <- coxph(fr, data=dat[rand,])
  r2<- inaccuracy_surv(model2= fr.fit, lp2=predict(fr.fit, data=dat[rand,], type="lp"), 
                       data=dat[rand,], time_name="y", event_name="cens")[[2]]
  return(full.r2-r2)
}

drop.r2(scope=scope, y=y, dr="x1")
drop.r2(scope=scope, y=y, dr="x3")
drop.r2(scope=scope, y=y, dr="x4a")
drop.r2(scope=scope, y=y, dr="x5")
drop.r2(scope=scope, y=y, dr="x6")
drop.r2(scope=scope, y=y, dr="x8")
drop.r2(scope=scope, y=y, dr="x10")
drop.r2(scope=scope, y=y, dr="x11")
drop.r2(scope=scope, y=y, dr="x15")



# save matrix with different simulation szenarios

n_szen <- 5
n <- c(300, 750, 1000, 5000, 10000)

setup <- matrix(cbind(n, 1:n_szen, 
                  matrix(beta, nrow = n_szen, ncol = length(beta), byrow = T)), 
                nrow = n_szen, ncol = 2 + length(beta),
                dimnames = list(1:n_szen, c("n", "seed", names(beta))))




save.image("setup.RData")

