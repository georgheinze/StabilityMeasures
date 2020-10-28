#####################################################################
##############  genereate and save simulation setup ################
#####################################################################


# note: the code for this simulation study ran on the Vienna Scientific Cluster. Results may be different if the code is executed on a single PC.


library(mvtnorm)
library(rsq)


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
beta[1]<- -2.9
linpred<-cbind(1,as.matrix(dat[,2:ncol(dat)])) %*% beta

ypred <- 1/(1+exp(-linpred))

hist(ypred)

y<- apply(ypred,1, function(x) sample(c(0,1), 1, prob=c(1-x,x)))
dat$y<-y

global<-as.formula(paste("y~",formr))





# determine Nagelkerke R2
full.r2 <- rsq.n(glm(formula=global, data=dat))

drop.r2<-function(scope=scope, y=y, data=dat, dr="x1"){
  which(scope!=dr)
  fr<-as.formula(paste("y~",paste(scope[which(scope!=dr)],collapse="+")))
  r2<-rsq.n(glm(formula=fr, data=dat))
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


# proportion of events
prop.table(table(dat$y))



# save matrix with different simulation szenarios

n_szen <- 5
n <- c(300, 750, 1000, 5000, 10000)

setup <- matrix(cbind(n, 1:n_szen, 
                      matrix(beta, nrow = n_szen, ncol = length(beta), byrow = T)), 
                nrow = n_szen, ncol = 2 + length(beta),
                dimnames = list(1:n_szen, c("n", "seed", names(beta))))





save.image("setup.RData")

