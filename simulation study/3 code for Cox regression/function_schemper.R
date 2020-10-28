inaccuracy_surv<-function(model2, model1=NULL, lp2, lp1=NULL, data, time_name, event_name){
  
  ### for Cox Models
  # model2 = model of interest
  # model1 = submodel of model2 - for comparison how much is gained by adding additional variables (model2)
  # data = data matrix
  # time_name = columnname of time-to-event from data matrix
  # event_name = columnname of event 0/1 from data matrix
  
  
  KM<-survfit(Surv(get(time_name),get(event_name))~1, data=data)
  op<-KM$surv
  G<-survfit(Surv(get(time_name),abs(1-get(event_name)))~1, data=data)
  times<- KM$time[which(KM$n.event>0)]
  times_ind <- which(KM$n.event>0)
  
  M0<-2*op*(1-op)
  
  D0<-1/sum(1/G$surv[times_ind]*KM$n.event[times_ind]) * 
      sum(1/G$surv[times_ind]*KM$n.event[times_ind]*M0[times_ind])
  
  
  M<-NULL
  bhaz_m2 <- basehaz(model2)
  for(i in 1:length(times)){
    t=times[i]
    pp= 1-exp(-bhaz_m2[which(KM$time==t),1])^exp(lp2)
    M[i]<-2*1/length(pp)*sum(pp*(1-pp))  
  }
  
  D2<-1/sum(1/G$surv[times_ind]*KM$n.event[times_ind]) * 
      sum(1/G$surv[times_ind]*KM$n.event[times_ind]*M)
  
  inacc<-round(c("D0"=D0,"D2"=D2),5)
  
  D1<-D0
  if(!is.null(model1)){
    
    M<-NULL
    bhaz_m1 <- basehaz(model1)
    for(i in 1:length(times)){
      t=times[i]
      pp= 1-exp(-bhaz_m1[which(KM$time==t),1])^exp(lp1)
      M[i]<-2*1/length(pp)*sum(pp*(1-pp))  
    }
    
    D1<-1/sum(1/G$surv[times_ind]*KM$n.event[times_ind]) * 
        sum(1/G$surv[times_ind]*KM$n.event[times_ind]*M)
    
    inacc<-round(c("D0"=D0,"D1"=D1,"D2"=D2),5)
  }
  
  expvar<-round((D1-D2)/D0,5)
  
  return(list("inaccuracy"=inacc, "explained_variation"=expvar))
  
}
