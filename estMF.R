library(MF)
library(ordinal)
library(bootstrap)
library(boot)
library(MASS)




nu_trt=2
b1=0.1725
size = 100
alpha=c(-1.5,-0.3,1,2,Inf)
alpha1=rep(alpha,nu_trt)
c=5
mlog_odds=rep(0,nu_trt*c)
treatment = rep( c(rep("trt", times = c), rep("con", times = c)))
dat = data.frame(treatment)
#View(dat)
for(j in 1:(nu_trt*c)){
  mlog_odds[j] = (alpha1[j]-b1*(dat$treatment[j] == "con"))
}
mcumproprobit=pnorm(mlog_odds)
#cumproplogit = plogis(log_odds)
mcell_pro=mcumproprobit
for(j in 2:(nu_trt*c)){72
  mcell_pro[j]=mcumproprobit[j]-mcumproprobit[j-1]
}
mcell_pro[seq(1,nu_trt*c,by=c)]=mcumproprobit[seq(1,nu_trt*c,by=c)]
dat$mcumpro=mcumproprobit
dat$mcelpro=mcell_pro
#View(dat)
#MF true value calculation
alpha=0
#p=c(0.02,0.28,0.1,0.2,0.4)
#v=c(0.3,0.25,0.15,0.25,0.05)
p=dat$mcelpro[6:10]
v=dat$mcelpro[1:5]
for(j in 1:length(p)){
  for(k in 1:length(v)){
    if(j==k){
      alpha=alpha+0.5*p[j]*v[k]
    }
    else if(j>k){
      alpha=alpha+p[j]*v[k]
    }
  }
}
true=2*alpha-1
#print(alpha)
#print(mf)
#true=0.081439394
count.probit=0
count.logit=0
count.loglog=0
count.sim=0
simN=1
simCI.probit=matrix(0,simN,2)
simCI.logit=matrix(0,simN,2)
simCI.loglog=matrix(0,simN,2)
simCI.cloglog=matrix(0,simN,2)
average_MF_over_5000.probit=rep(0,simN)
average_MF_over_5000.logit=rep(0,simN)
average_MF_over_5000.loglog=rep(0,simN)
average_MF_over_5000.cloglog=rep(0,simN)
log_odds=rep(0,nu_trt*c)
for(k in 1:simN){
  ###probit model confidence interval
  for(j in 1:(nu_trt*c)){
    log_odds[j] = alpha1[j]-b1*(dat$treatment[j] == "con")
  }
  cumproprobit=pnorm(log_odds)
  #cumproplogit = plogis(log_odds)
  cell_pro=cumproprobit
  for(j in 2:(nu_trt*c)){
    cell_pro[j]=cumproprobit[j]-cumproprobit[j-1]
  }
  cell_pro[seq(1,nu_trt*c,by=c)]=cumproprobit[seq(1,nu_trt*c,by=c)]
  dat$cumpro=cumproprobit
  dat$celpro=cell_pro
  #View(dat)
  fi=seq(1,(nu_trt*c-c+1), by=c)
  y=rep(0,nu_trt*c)
  for(i in fi ){
    #rmultinom(1, size = 110, prob = dat$celpro[i:(i+4)])
    #checkcell[i:(i+4)]=dat$celpro[i:(i+4)]
    y[i:(i+4)]= rmultinom(1, size = size, prob = dat$celpro[i:(i+4)])
  }
  dat$y=y
  #View(dat)
  #write.csv(dat, "ordinalsimulation2.csv", col.names = TRUE, row.names = FALSE)
  #size = 20#
  #ya=rep(0,size*nu_trt)
  trt_t=rep(0,size*nu_trt)
  #m=c(1,11)
  xt=1
  #i=11
  severity <- c(1,2,3,4,5)
  trt=dat$y[1:c]
  con=dat$y[(c+1):(c*nu_trt)]
  yb <- factor(c(rep(severity, con),
                 rep(severity, trt)), levels = severity)
  trt_t[xt:(xt+(size-1))]=1
  trt_t[((xt+(size-1))+1):(((xt+(size-1))+1)+(size-1))]=0
  newdata=data.frame(trt_t,yb)
  #View(newdata)
  y=factor(newdata$yb)
  trt=newdata$trt_t
  # blk=newdata$blk
  ## Cumulative link mixed model with two random terms:
  #probit.m <- clmm(y ~ trt + (1|blk), link = "probit",threshold = "equidistant")
  probit.mm <- clm(y ~ trt, link = "probit")
  summary(probit.mm)
  betap=tail(coef(probit.mm),n=1)
  alpha=pnorm(betap/sqrt(2),0,1)
  MF=2*alpha-1
  print(MF)
  ###probit model confidence interval
  #library(bootstrap)
  #library(boot)
  theta.hat=function(d,i)
  {
    yn<<-factor(d$yb[i])
    xn<<-d$trt_t[i]
    probit.m <- clm(yn ~ xn , link = "probit")
    betap=tail(coef(probit.m),n=1)
    MF=2*pnorm(betap/sqrt(2),0,1)-1
  }
  oboot=boot(data=newdata,statistic=theta.hat, R=5000)
  BCa.probit=boot.ci(oboot,conf=.95,type="bca")
  simCI.probit[k,]=BCa.probit$bca[4:5]
  if(true>=simCI.probit[k,1] && true<=simCI.probit[k,2]){
    count.probit=count.probit+1}
  average_MF_over_5000.probit[k]=oboot$t0
  count.sim=count.sim+1
}
print(count.probit)

library(tidyverse)
simdat <- newdata %>%
  group_by(trt_t) %>%
  mutate(id = row_number()) %>%
  mutate(ftrt = as.factor(trt_t)) %>%
  mutate(ftrt = recode_factor(ftrt, "0" = "con", "1" = "vac")) %>%
  mutate(nyb = as.numeric(as.character(yb))) %>%
  mutate(nyb2 = nyb*20)

simdat %>%
  ggplot(., aes(x=nyb2, y=as.factor(id), color=as.factor(trt_t))) +
    geom_point() + 
    facet_wrap(~paste("trt:", trt_t))



MFBoot(nyb2 ~ ftrt, data=simdat, compare=c("con", "vac"))






#Bias probit
library(ordinal)
library(bootstrap)
library(boot)
library(MASS)
nu_trt=2
b1=0.1725
size = 100
alpha=c(-1.5,-0.3,1,2,Inf)
alpha1=rep(alpha,nu_trt)
c=5
mlog_odds=rep(0,nu_trt*c)
treatment = rep( c(rep("trt", times = c), rep("con", times = c)))
dat = data.frame(treatment)
#View(dat)
for(j in 1:(nu_trt*c)){
  mlog_odds[j] = (alpha1[j]-b1*(dat$treatment[j] == "con"))
}
mcumproprobit=pnorm(mlog_odds)
#cumproplogit = plogis(log_odds)
mcell_pro=mcumproprobit
for(j in 2:(nu_trt*c)){
  mcell_pro[j]=mcumproprobit[j]-mcumproprobit[j-1]
}
mcell_pro[seq(1,nu_trt*c,by=c)]=mcumproprobit[seq(1,nu_trt*c,by=c)]
dat$mcumpro=mcumproprobit
dat$mcelpro=mcell_pro
#View(dat)
#MF true value calculation
alpha=0
#p=c(0.02,0.28,0.1,0.2,0.4)
#v=c(0.3,0.25,0.15,0.25,0.05)
p=dat$mcelpro[6:10]
v=dat$mcelpro[1:5]
for(j in 1:length(p)){
  for(k in 1:length(v)){
    if(j==k){
      alpha=alpha+0.5*p[j]*v[k]
    }
    else if(j>k){
      alpha=alpha+p[j]*v[k]
    }
  }
}
true=2*alpha-1
#print(alpha)
#print(mf)
#true=0.081439394
count.probit=0
count.logit=0
count.loglog=0
count.sim=0
simN=10000
simCI.probit=matrix(0,simN,2)
simCI.logit=matrix(0,simN,2)
simCI.loglog=matrix(0,simN,2)
simCI.cloglog=matrix(0,simN,2)
average_MF_over_5000.probit=rep(0,simN)
average_MF_over_5000.logit=rep(0,simN)
average_MF_over_5000.loglog=rep(0,simN)
average_MF_over_5000.cloglog=rep(0,simN)
log_odds=rep(0,nu_trt*c)
for(k in 1:simN){
  ###probit model confidence interval
  for(j in 1:(nu_trt*c)){
    log_odds[j] = alpha1[j]-b1*(dat$treatment[j] == "con")
  }
  cumproprobit=pnorm(log_odds)
  #cumproplogit = plogis(log_odds)
  cell_pro=cumproprobit
  for(j in 2:(nu_trt*c)){
    cell_pro[j]=cumproprobit[j]-cumproprobit[j-1]
  }
  cell_pro[seq(1,nu_trt*c,by=c)]=cumproprobit[seq(1,nu_trt*c,by=c)]
  dat$cumpro=cumproprobit
  dat$celpro=cell_pro
  #View(dat)
  fi=seq(1,(nu_trt*c-c+1), by=c)
  y=rep(0,nu_trt*c)
  for(i in fi ){
    #rmultinom(1, size = 110, prob = dat$celpro[i:(i+4)])
    #checkcell[i:(i+4)]=dat$celpro[i:(i+4)]
    y[i:(i+4)]= rmultinom(1, size = size, prob = dat$celpro[i:(i+4)])
  }
  dat$y=y
  #View(dat)
  #write.csv(dat, "ordinalsimulation2.csv", col.names = TRUE, row.names = FALSE)
  #size = 20#
  #ya=rep(0,size*nu_trt)
  trt_t=rep(0,size*nu_trt)
  #m=c(1,11)
  xt=1
  #i=11
  severity <- c(1,2,3,4,5)
  trt=dat$y[1:c]
  con=dat$y[(c+1):(c*nu_trt)]
  yb <- factor(c(rep(severity, con),
                 rep(severity, trt)), levels = severity)
  trt_t[xt:(xt+(size-1))]=1
  trt_t[((xt+(size-1))+1):(((xt+(size-1))+1)+(size-1))]=0
  newdata=data.frame(trt_t,yb)
  #View(newdata)
  y=factor(newdata$yb)
  trt=newdata$trt_t
  # blk=newdata$blk
  ## Cumulative link mixed model with two random terms:
  #probit.m <- clmm(y ~ trt + (1|blk), link = "probit",threshold = "equidistant")
  probit.mm <- clm(y ~ trt, link = "probit")
  summary(probit.mm)
  betap=tail(coef(probit.mm),n=1)
  alpha=pnorm(betap/sqrt(2),0,1)
  MF=2*alpha-1
  average_MF_over_5000.probit[k]=MF
}
#Bias Logit
library(ordinal)
library(bootstrap)
library(boot)
library(MASS)
nu_trt=2
b1=0.2842
size = 100
alpha=c(-1.5,-0.3,1,2,Inf)
alpha1=rep(alpha,nu_trt)
c=5
mlog_odds=rep(0,nu_trt*c)
treatment = rep( c(rep("trt", times = c), rep("con", times = c)))
dat = data.frame(treatment)
#View(dat)
for(j in 1:(nu_trt*c)){
  mlog_odds[j] = (alpha1[j]-b1*(dat$treatment[j] == "con"))
}
mcumproprobit=plogis(mlog_odds)
#cumproplogit = plogis(log_odds)
mcell_pro=mcumproprobit
for(j in 2:(nu_trt*c)){
  mcell_pro[j]=mcumproprobit[j]-mcumproprobit[j-1]
}
mcell_pro[seq(1,nu_trt*c,by=c)]=mcumproprobit[seq(1,nu_trt*c,by=c)]
dat$mcumpro=mcumproprobit
dat$mcelpro=mcell_pro
#View(dat)
#MF true value calculation
alpha=0
#p=c(0.02,0.28,0.1,0.2,0.4)
#v=c(0.3,0.25,0.15,0.25,0.05)
p=dat$mcelpro[6:10]
v=dat$mcelpro[1:5]
for(j in 1:length(p)){
  for(k in 1:length(v)){
    if(j==k){
      alpha=alpha+0.5*p[j]*v[k]
    }
    else if(j>k){
      alpha=alpha+p[j]*v[k]
    }
  }
}
true=2*alpha-1
print(true)
#print(mf)
#true=0.081439394
simN=10000
average_MF_over_5000.logit=rep(0,simN)
log_odds=rep(0,nu_trt*c)
for(k in 1:simN){
  ###probit model confidence interval
  for(j in 1:(nu_trt*c)){
    log_odds[j] = alpha1[j]-b1*(dat$treatment[j] == "con")
  }
  cumproprobit=plogis(log_odds)
  #cumproplogit = plogis(log_odds)
  cell_pro=cumproprobit
  for(j in 2:(nu_trt*c)){
    cell_pro[j]=cumproprobit[j]-cumproprobit[j-1]
  }
  cell_pro[seq(1,nu_trt*c,by=c)]=cumproprobit[seq(1,nu_trt*c,by=c)]
  dat$cumpro=cumproprobit
  dat$celpro=cell_pro
  #View(dat)
  fi=seq(1,(nu_trt*c-c+1), by=c)
  y=rep(0,nu_trt*c)
  for(i in fi ){
    #rmultinom(1, size = 110, prob = dat$celpro[i:(i+4)])
    #checkcell[i:(i+4)]=dat$celpro[i:(i+4)]
    y[i:(i+4)]= rmultinom(1, size = size, prob = dat$celpro[i:(i+4)])
  }
  dat$y=y
  #View(dat)
  #write.csv(dat, "ordinalsimulation2.csv", col.names = TRUE, row.names = FALSE)
  #size = 20#
  #ya=rep(0,size*nu_trt)
  trt_t=rep(0,size*nu_trt)
  #m=c(1,11)
  xt=1
  #i=11
  severity <- c(1,2,3,4,5)
  trt=dat$y[1:c]
  con=dat$y[(c+1):(c*nu_trt)]
  yb <- factor(c(rep(severity, con),
                 rep(severity, trt)), levels = severity)
  trt_t[xt:(xt+(size-1))]=1
  trt_t[((xt+(size-1))+1):(((xt+(size-1))+1)+(size-1))]=0
  newdata=data.frame(trt_t,yb)
  #View(newdata)
  y=factor(newdata$yb)
  trt=newdata$trt_t
  # blk=newdata$blk
  ## Cumulative link mixed model with two random terms:
  #probit.m <- clmm(y ~ trt + (1|blk), link = "probit",threshold = "equidistant")
  probit.m <- clm(y ~ trt, link = "logit")
  betap=tail(coef(probit.m),n=1)
  n1=10000
  y1=rlogis(n1, location = betap, scale = 1)
  y2=rlogis(n1, location = 0, scale = 1)
  MF=2*mean(y1>y2)-1
  average_MF_over_5000.logit[k]=MF}


