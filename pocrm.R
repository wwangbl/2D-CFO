library(pocrm)


pocrm.imp <- function(alpha,prior.o,theta,y,combos){
    
  data<-as.matrix(table(combos,y))
  ##############
  if(length(data[1,])==1){
    data <- cbind(data,rep(0,length(data[,1])))
    colnames(data) <- c(0,1)
  }
  ##############
  level<-as.numeric(row.names(data))
  nontox<-as.numeric(data[,1])
  tox<-as.numeric(data[,2])
  
  apred<-rep(0,nrow(alpha))
  lik<-rep(0,nrow(alpha))
  for(k in 1:nrow(alpha)){			
    ll<-function(a){
      la<-0 #value of the log-likelihood
      for(i in level){
        index<-match(i,level)
        la<-la+tox[index]*a*log(alpha[k,][i])+nontox[index]*log((1-alpha[k,][i]**a))
      }
      la
    }
    apred[k]<-optimize(f=ll,interval=c(0,500),maximum=T)$maximum
    lik[k]<-ll(apred[k])
  }
  pord<-(exp(lik)*prior.o)/sum(exp(lik)*prior.o)
  #library("nnet") not necessary because listed as package dependency
  ord<-which.is.max(pord)
  ahat<-apred[ord]
  rpred<-alpha[ord,]**ahat
  next.lev<-which.is.max(-(abs(rpred-theta)))
  out<-list(ord.prob=round(pord,3),order.est=ord,a.est=round(ahat,3),ptox.est=round(rpred,3),dose.rec=next.lev)
}


pocrm.sim <- function(r,alpha,prior.o,x0,n,theta,seed){
  set.seed(seed)
  y <- combos <- c()
  correct <- npercent <- ntox <- ptoxic <- 0

  # preliminary stage
  # for (x in x0) {
  #   combos <- c(combos, x)
  # 
  #   p<-runif(1)
  #   if (p<=r[x]){
  #     y <- c(y,1)
  #   } else {
  #     y <- c(y,0)
  #   }
  # 
  # 
  #   if(sum(y)>0){
  #     break
  #   }
  # }
  
  # without pre stage
  combos <- c(combos, rep(1,3))
  for (i in 1:3) {
    p<-runif(1)
    if (p<=r[1]){
      y <- c(y,1)
    } else {
      y <- c(y,0)
    }
  }
  

  while ((length(y)+3)<=n) {
    idx <- pocrm.imp(alpha,prior.o,theta,y,combos)$dose.rec
   # message(idx)
    
    combos <- c(combos, rep(idx,3))
    for (i in 1:3) {
      p<-runif(1)
      if (p<=r[idx]){
        y <- c(y,1)
      } else {
        y <- c(y,0)
      }
    }
  }
  
  ## patients left
  if(length(y)<n){
    idx <- pocrm.imp(alpha,prior.o,theta,y,combos)$dose.rec
    combos <- c(combos, rep(idx,(n-length(y))))
    for (i in 1:(n-length(y))) {
      p<-runif(1)
      if (p<=r[idx]){
        y <- c(y,1)
      } else {
        y <- c(y,0)
      }
    }
  }
  #####
  
  MTD <- pocrm.imp(alpha,prior.o,theta,y,combos)$dose.rec
  ntox <- sum(y)
  if(r[MTD]==theta){
    correct <- 1
  }
  for (i in 1:length(r)) {
    if(r[i]==theta){
      npercent <- npercent + length(which(combos==i))
    }
  }
  npercent <- npercent/n
  
  for (i in 1:length(r)) {
    if(r[i]>theta){
      ptoxic <- ptoxic + length(which(combos==i))
    }
  }
  ptoxic <- ptoxic/n
  
  list(correct=correct, npercent=npercent, ntox=ntox, ptoxic=ptoxic)
}
