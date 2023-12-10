library(readxl)
library(dplyr)
library(VGAM)
#
data.use.contact=as.data.frame(read_xlsx("close contacts51786.xlsx"))
data_check_pair=as.data.frame(read_xlsx("total_transmission_pairs1.xlsx"))
data_case=as.data.frame(read.csv("total_case_data.csv"))
data_case_check=as.data.frame(read_xlsx("positive case 1139.xlsx"))
#classification of cases
#sporadic cases
all.sporadic.ID=unique(data_check_pair$infectee.case[data_check_pair$infectee.case==data_check_pair$infector.case])
sum(data_check_pair$cluster.risk=="sporadic cases")
#infector cases
all.infector.ID=unique(data_check_pair$infector.case)[!(unique(data_check_pair$infector.case)%in%all.sporadic.ID)]
#terminal cases
all.terminal.ID=unique(data_check_pair$infectee.case)[!(unique(data_check_pair$infectee.case)%in%c(all.sporadic.ID,all.infector.ID))]
#check
unique(data_all_contacts$infector.case_number)%in%all.sporadic.ID
##
data_case %>% 
  filter(caseid%in%all.infector.ID,
         age<18) %>%
  dplyr::pull(caseid) -> young.infector.ID
#
data_case %>% 
  filter(caseid%in%all.sporadic.ID,
         age<18) %>%
  dplyr::pull(caseid) -> young.sporadic.ID
#
data_case %>% 
  filter(caseid%in%all.terminal.ID,
         age<18) %>%
  dplyr::pull(caseid) -> young.terminal.ID
#
all_young_case.ID=c(young.infector.ID,young.sporadic.ID,young.terminal.ID)
#
pair=data_check_pair %>%
  mutate(infectee.case=ifelse(cluster.risk=="sporadic cases",NA,infectee.case),
         generation=as.factor(cluster.generation),
         infector.case=ifelse(cluster.risk=="sporadic cases",NA,infector.case))
#
#combine individual case information across datasets
offspring_data=data.frame(infector.ID=rep(NA,length(all_young_case.ID)),
                          offspring=rep(NA,length(all_young_case.ID)),
                          infector.age=rep(NA,length(all_young_case.ID)),
                          infector.sex=rep(NA,length(all_young_case.ID)),
                          infector.symptom.state=rep(NA,length(all_young_case.ID)),
                          infector.vac.dose=rep(NA,length(all_young_case.ID)),
                          case.class=rep(NA,length(all_young_case.ID)),
                          link=rep(NA,length(all_young_case.ID)),
                          household=rep(NA,length(all_young_case.ID)),
                          non_household=rep(NA,length(all_young_case.ID)),
                          #before=rep(NA,length(all_young_case.ID)),
                          #after=rep(NA,length(all_young_case.ID)),
                          terminal=rep(NA,length(all_young_case.ID)),
                          generation=rep(NA,length(all_young_case.ID)))
#
for (i in 1:nrow(offspring_data)){
  off.ID=pair$infectee.case[which(pair$infector.case==all_young_case.ID[i])]
  offspring_data$infector.ID[i]=all_young_case.ID[i]
  offspring_data$offspring[i]=nrow(filter(pair,infector.case==all_young_case.ID[i]))
  offspring_data$infector.age[i]=data_case$age[data_case$caseid==all_young_case.ID[i]]
  offspring_data$infector.sex[i]=data_other_check$`gender: Male Female`[data_other_check$ID==all_young_case.ID[i]]
  offspring_data$infector.onset.date[i]=data_case$onset.date[data_case$caseid==all_young_case.ID[i]]
  offspring_data$infector.symptom.state[i]=data_case$symptomatic[data_case$caseid==all_young_case.ID[i]]
  offspring_data$infector.vac.dose[i]=data_case$Vaccination[data_case$caseid==all_young_case.ID[i]]
  offspring_data$case.class[i]=data_case$imported[data_case$caseid==all_young_case.ID[i]]
  offspring_data$link[i]=ifelse(all_young_case.ID[i]%in%young.sporadic.ID,"Sporadic",
                                "Epi-linked")
  offspring_data$household[i]=sum(data_case$setting[data_case$caseid%in%off.ID]=="household")
  offspring_data$non_household[i]=sum(data_case$setting[data_case$caseid%in%off.ID]!="household")
  #offspring_data$before[i]=sum(data_other_check$npi[data_case$caseid%in%off.ID]=="before")
  #offspring_data$after[i]=sum(data_other_check$npi[data_case$caseid%in%off.ID]=="after")
  offspring_data$terminal[i]=ifelse(all_young_case.ID[i]%in%young.terminal.ID,"yes","no")
  offspring_data$generation[i]=max(as.numeric(pair%>%filter(infector.case==all_young_case.ID[i])%>%pull(generation)))
}
#
#
#offspring_data$generation=ifelse(offspring_data$generation==-Inf,NA,offspring_data$generation)
#offspring_data$generation=ifelse(offspring_data$link=="Sporadic",1,offspring_data$generation)
#
offspring_data$npi=ifelse(offspring_data$infector.ID%in%data_other_check$ID[which(data_other_check$npi=="before")],
                          "before","after")
#
offspring_data=offspring_data %>%
  mutate(age_group=case_when(
    infector.age<6~"0~5",
    infector.age<=12&infector.age>=6~"6~12",
    infector.age>12~"13~17"
  ))
#
####
#negative binomial model
nb.lik=\(parms,data){
  #R=parms[1]
  #k=parms[2]
  k=exp(parms[2])
  R=exp(parms[1])
  ll=-sum(log(dnbinom(data,mu=R,size=k)))
  return(ll)
}
#MCMC for R & k estimations
mcmc=\(burnin,M,dat,fn){
  like=fn
  R0=1
  k=1
  lb=0.000001
  ub=20
  R.lb=-5;R.ub=(2.5)
  k.lb=-5;k.ub=(3.5)
  #k.lb=lb;k.ub=ub
  #R.lb=lb;R.ub=ub
  sdR0=0.4
  sdk=0.2
  collect=c(0,R0,k,0,0)
  for (iter in 1:(burnin+M-1)){
    R0acpt=0
    kacpt=0
    R0new=rnorm(1,R0,sdR0)
    if(R0new>R.lb&R0new<R.ub){
      R0p=like(c(R0new,k),dat)
      if (is.finite(R0p)){
        testR0=exp(-1*(R0p-like(c(R0,k),dat)))
        if (runif(1)<min(1,testR0)) {
          R0=R0new
          R0acpt=1
        }
      }
    }
    knew=rnorm(1,k,sdk)
    if(knew>k.lb&knew<k.ub){
      kp=like(c(R0,knew),dat)
      if(is.finite(kp)){
        testk=exp(-1*(kp-like(c(R0,k),dat)))
        if(runif(1)<min(1,testk)){
          k=knew
          kacpt=1
        }
      }
    }
    collect=rbind(collect,c(iter,R0,k,R0acpt,kacpt))
  }
  return(collect)
}
burnin=40000
M=100000
#
set.seed(45782)
fit.all=mcmc(burnin,M,offspring_data$offspring,nb.lik)
hist(fit.all[(burnin+1):(burnin+M),2])
hist(fit.all[(burnin+1):(burnin+M),3])
quantile(exp(fit.all[(burnin+1):(burnin+M),2]),c(0.5,0.025,0.975))
quantile(exp(fit.all[(burnin+1):(burnin+M),3]),c(0.5,0.025,0.975))
#
#expected proportion of index case seeding 80% of all transmission, adapted from DOI: 10.12688/wellcomeopenres.15842.3
#
propresponsible=function(R0,k,prop){
  qm1=qnbinom(1-prop,k+1,mu=R0*(k+1)/k)
  remq=1-prop-pnbinom(qm1-1,k+1,mu=R0*(k+1)/k)
  remx=remq/dnbinom(qm1,k+1,mu=R0*(k+1)/k)
  q=qm1+1
  1-pnbinom(q-1,k,mu=R0)-dnbinom(q,k,mu=R0)*remx
}
#
prop80_all=quantile(propresponsible(exp(fit.all[(burnin+1):(burnin+M),2]),
                                    exp(fit.all[(burnin+1):(burnin+M),3]),
                                    0.8),c(0.5,0.025,0.975))
#
##construct data for SAR analyses
young_close_contact=offspring_data %>%
  filter(infector.ID%in%SAR.young.ID) %>%
  dplyr::select(-c(household,non_household,infector.onset.date,before,after))
#
#
data_SAR=data.use.contact %>%
  filter(infector.case_number%in%all_young_case.ID) %>%
  mutate(
    infector.lag=as.numeric(infector.lag),
    npi=case_when(
    infector.case_number%in%data_other_check$ID[which(data_other_check$npi=="after")]~"after",
    infector.case_number%in%data_other_check$ID[which(data_other_check$npi=="before")]~"before"
  ))
#
for (i in 1:nrow(young_close_contact)){
  young_close_contact$vac.lag[i]=unique(data_SAR$infector.lag[data_SAR$infector.case_number==young_close_contact$infector.ID[i]])
  young_close_contact$total.num.contact[i]=nrow(filter(data_SAR,infector.case_number==young_close_contact$infector.ID[i]))
  young_close_contact$positive.contact[i]=nrow(data_SAR %>% filter(outcome==1,infector.case_number==young_close_contact$infector.ID[i]))
  young_close_contact$household.total[i]=nrow(filter(data_SAR,infector.case_number==young_close_contact$infector.ID[i],
                                                     contact_setting=="household"))
  young_close_contact$household.positive[i]=nrow(filter(data_SAR,infector.case_number==young_close_contact$infector.ID[i],
                                                        contact_setting=="household",
                                                        outcome==1))
  young_close_contact$non_household.total[i]=nrow(filter(data_SAR,infector.case_number==young_close_contact$infector.ID[i],
                                                         contact_setting=="no-household"))
  young_close_contact$non_household.positive[i]=nrow(filter(data_SAR,infector.case_number==young_close_contact$infector.ID[i],
                                                            contact_setting=="no-household",
                                                            outcome==1))
  young_close_contact$before.total[i]=nrow(filter(data_SAR,infector.case_number==young_close_contact$infector.ID[i],
                                                  npi=="before"))
  young_close_contact$before.positive[i]=nrow(filter(data_SAR,infector.case_number==young_close_contact$infector.ID[i],
                                                     npi=="before",
                                                     outcome==1))
  young_close_contact$after.total[i]=nrow(filter(data_SAR,infector.case_number==young_close_contact$infector.ID[i],
                                                 npi=="after"))
  young_close_contact$after.positive[i]=nrow(filter(data_SAR,infector.case_number==young_close_contact$infector.ID[i],
                                                    npi=="after",
                                                    outcome==1))
  
}
#
##beta-binomial model
beta_bi_lik=\(parms,dat){
  total=dat$total.num.contact
  positive=dat$positive.contact
  a=parms[1]
  b=parms[2]
  ll=-sum(log(dbetabinom.ab(x = positive, 
                            size = total, 
                            shape1 = a, 
                            shape2 = b)))
  return(ll)
}
#mcmc for SAR estimation
mcmc=\(burnin,M,dat,fn){
  like=fn
  a=1
  b=1
  lb=0
  ub=100
  b.lb=lb;b.ub=ub
  a.lb=lb;a.ub=ub
  sda=0.4
  sdb=0.2
  collect=c(0,a,b,0,0)
  for (iter in 1:(burnin+M-1)){
    aacpt=0
    bacpt=0
    anew=rnorm(1,a,sda)
    if(anew>a.lb&anew<a.ub){
      ap=like(c(anew,b),dat)
      if (is.finite(ap)){
        testa=exp(-1*(ap-like(c(a,b),dat)))
        if (runif(1)<min(1,testa)) {
          a=anew
          aacpt=1
        }
      }
    }
    bnew=rnorm(1,b,sdb)
    if(bnew>b.lb&bnew<b.ub){
      bp=like(c(a,bnew),dat)
      if(is.finite(bp)){
        testb=exp(-1*(bp-like(c(a,b),dat)))
        if(runif(1)<min(1,testb)){
          b=bnew
          bacpt=1
        }
      }
    }
    collect=rbind(collect,c(iter,a,b,aacpt,bacpt))
  }
  return(collect)
}
#
set.seed(83102)
fit.SAR.all=mcmc(burnin,M,young_close_contact,beta_bi_lik)
hist(fit.SAR.all[(burnin+1):(burnin+M),2])
hist(fit.SAR.all[(burnin+1):(burnin+M),3])
a.all=fit.SAR.all[(burnin+1):(burnin+M),2]
b.all=fit.SAR.all[(burnin+1):(burnin+M),3]
#
mean.beta <- rep(0,length(a.all))
sd.beta <- rep(0,length(a.all))
median.beta <- rep(0,length(a.all))
ub.beta <- rep(0,length(a.all))
mode.beta <- rep(0,length(a.all))
CV.all <- rep(0,length(a.all))
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
for (i in 1:length(a.all)) {
  x <- rbeta(length(a.all),shape1 = a.all[i],shape2 = b.all[i])
  mean.beta[i] <- mean(x)
  sd.beta[i] <-sd(x)
  median.beta[i] <-median(x)
  ub.beta[i] <-quantile(x,0.95)
  mode.beta[i] <- Mode(x)
  CV.all[i] <- sd(x)/mean(x)
}
quantile(mean.beta,c(0.5,0.025,0.975))
quantile(median.beta,c(0.5,0.025,0.975))
quantile(sd.beta,c(0.5,0.025,0.975))
quantile(ub.beta,c(0.5,0.025,0.975))
quantile(mode.beta,c(0.5,0.025,0.975))
quantile(CV.all, c(0.5,0.025,0.975))













