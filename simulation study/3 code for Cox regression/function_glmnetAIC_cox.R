### utility functions for glmnet objects


# fast implementation using Akaike Information Criterion (AIC)
AIC.glmnet<-function(x, y, alpha=1, family, ...){
  require(glmnet)
  fit<-glmnet(x=x, y=y, alpha=alpha, family=family, ...)
  dev<-(1-fit$dev.ratio)*fit$nulldev
  index.opt<-(1:length(fit$lambda))[(dev+2*fit$df)==min(dev+2*fit$df)]
  df.opt<-fit$df[index.opt]
  lambda.opt<-fit$lambda[index.opt]
  beta.opt<- fit$beta[,index.opt]
  #names(beta.opt)[1] <- "(Intercept)"
  aic<-dev+2*fit$df
  aic.opt<-dev[index.opt]+2*fit$df[index.opt]
  res<-list(glmnet.fit=fit, dev=dev, index.opt=index.opt, df.opt=df.opt, lambda.opt=lambda.opt, lambda=fit$lambda, nzero=fit$df, beta.opt=beta.opt, aic=aic, aic.opt=aic.opt,
  selected=(beta.opt!=0))
  attr(res, "class")<-"AIC.glmnet"
  return(res)
  }
