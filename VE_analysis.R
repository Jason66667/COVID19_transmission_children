library(readxl)
library(dplyr)
library(arm)
#VE analyses
###dataset for VE against all infection
data.use.contact=read_xlsx("close contacts51786.xlsx")
outbreak.initial.date=as.Date("2022-07-31")
data_infect=data.use.contact %>%
  filter(infectee.age<18,
         !is.na(infectee.expose_date)) %>%
  mutate(
    infectee.expose_date=as.numeric(as.Date(infectee.expose_date)-outbreak.initial.date),
    infectee_lag=as.numeric(infectee_lag),
    age_group=case_when(
      infectee.age>=0&infectee.age<=12~"0~12",
      infectee.age>=13&infectee.age<=17~"13~17"
    )) %>%
  mutate(
    infector.vac_dose=case_when(
      infector.vac_dose==0~"0 dose",
      infector.vac_dose==1~"1 dose",
      infector.vac_dose==2~"2 dose",
      infector.vac_dose==3~"3 dose",
    ),
    infectee.vac_dose=case_when(
      infectee.vac_dose==0~"0 dose",
      infectee.vac_dose==1~"1 dose",
      infectee.vac_dose==2~"2 dose",
      infectee.vac_dose==3~"3 dose",
    ),
  ) %>%
  filter(infectee.vac_dose=="1 dose"|(infectee.vac_dose=="2 dose"&infectee_lag>=15))
#
data_infect_2_dose=data_infect %>%
  filter(infectee.vac_dose%in%c("1 dose","2 dose")) %>%
  mutate(infectee.lag.new=case_when(
    ((infectee_lag>365&infectee.vac_dose=="2 dose")|infectee.vac_dose=="1 dose")~"365~",
    (infectee_lag<=180&infectee.vac_dose=="2 dose")~"~180",
    (infectee_lag>180&infectee_lag<=365&infectee.vac_dose=="2 dose")~"180~365"
  ))
#
#
data_infect_2_dose=data_infect %>%
  filter(infectee.vac_dose%in%c("1 dose","2 dose")) %>%
  mutate(infectee.lag.new=case_when(
    ((infectee_lag>365&infectee.vac_dose=="2 dose")|infectee.vac_dose=="1 dose")~"365~",
    (infectee_lag<=365&infectee.vac_dose=="2 dose")~"~365"
  ))
#
##Bayesian generalized linear model
M_2dose_lag <- bayesglm(outcome ~ infectee.lag.new
                        +infectee.sex
                        +infectee.age
                        +infectee.expose_date
                        +infector.vac_dose
                        +contact_setting,
                        family=binomial(link="logit"),
                        data = data_infect_2_dose, 
                        prior.scale=Inf, prior.df=Inf)
##
summary(M_2dose_lag)
#
###dataset for VE against transmission
data_trans=data.use.contact %>%
  filter(infector.age<18,
         !is.na(infectee.expose_date)
         !is.na(infector.expose_date),
         ) %>%
  mutate(
    infectee.expose_date=as.numeric(as.Date(infectee.expose_date)-outbreak.initial.date),
    infector.expose_date=as.numeric(as.Date(infector.expose_date)-outbreak.initial.date),
    infector_lag=as.numeric(infectee_lag),
    age_group=case_when(
      infector.age>=0&infector.age<=12~"0~12",
      infector.age>=13&infector.age<=17~"13~17"
    )) %>%
  mutate(
    infector.vac_dose=case_when(
      infector.vac_dose==0~"0 dose",
      infector.vac_dose==1~"1 dose",
      infector.vac_dose==2~"2 dose",
      infector.vac_dose==3~"3 dose",
    ),
    infectee.vac_dose=case_when(
      infectee.vac_dose==0~"0 dose",
      infectee.vac_dose==1~"1 dose",
      infectee.vac_dose==2~"2 dose",
      infectee.vac_dose==3~"3 dose",
    ),
  ) %>%
  filter(infector.vac_dose=="1 dose"|(infector.vac_dose=="2 dose"&infector_lag>=15))
#
data_trans_2_dose=data_trans %>%
  filter(infector.vac_dose%in%c("1 dose","2 dose")) %>%
  mutate(infector.lag.new=case_when(
    ((infector_lag>365&infector.vac_dose=="2 dose")|infector.vac_dose=="1 dose")~"365~",
    (infector_lag<=180&infector.vac_dose=="2 dose")~"~180",
    (infector_lag>180&infector_lag<=365&infector.vac_dose=="2 dose")~"180~365"
  ))
#
data_trans_2_dose=data_trans %>%
  filter(infector.vac_dose%in%c("1 dose","2 dose")) %>%
  mutate(infector.lag.new=case_when(
    ((infector_lag>365&infector.vac_dose=="2 dose")|infector.vac_dose=="1 dose")~"365~",
    (infector_lag<=365&infector.vac_dose=="2 dose")~"~365"
  ))
#
##Bayesian generalized linear model
M_trans_2dose_lag <- bayesglm(outcome ~ infector.lag.new
                              +infector.age
                              +infector.sex
                              +infectee.expose_date
                              +infectee.vac_dose
                              +contact_setting, 
                              family=binomial(link="logit"),
                              data = data_trans_2_dose, 
                              prior.scale=Inf, prior.df=Inf))
#
summary(M_trans_2dose_lag)
#
###dataset for VE against symptomatic infection
data_case_check=as.data.frame(read_xlsx("positive case 1139.xlsx"))
young_case_ID=data_case_check$infectee.case_ID[data_case_check$age<18&(!is.na(data_case_check$`symptom onset date`))]
data_young_symptom=data.use.contact %>%
  filter(infectee.age<18,
         !is.na(infectee.expose_date)) %>%
  mutate(
    symptom.state=ifelse(infectee.case_ID %in% young_case_ID,1,0),
    infectee.expose_date=as.numeric(as.Date(infectee.expose_date)-outbreak.initial.date),
    infectee_lag=as.numeric(infectee_lag),
    age_group=case_when(
      infectee.age>=0&infectee.age<=12~"0~12",
      infectee.age>=13&infectee.age<=17~"13~17"
    )) %>%
  mutate(
    infector.vac_dose=case_when(
      infector.vac_dose==0~"0 dose",
      infector.vac_dose==1~"1 dose",
      infector.vac_dose==2~"2 dose",
      infector.vac_dose==3~"3 dose",
    ),
    infectee.vac_dose=case_when(
      infectee.vac_dose==0~"0 dose",
      infectee.vac_dose==1~"1 dose",
      infectee.vac_dose==2~"2 dose",
      infectee.vac_dose==3~"3 dose",
    ),
  ) %>%
  filter(infectee.vac_dose=="1 dose"|(infectee.vac_dose=="2 dose"&infectee_lag>=15))
#
data_young_symptom_lag=data_young_symptom %>%
  filter(infectee.vac_dose%in%c("1 dose","2 dose")) %>%
  mutate(infectee.lag.new=case_when(
    ((infectee_lag>365&infectee.vac_dose=="2 dose")|infectee.vac_dose=="1 dose")~"365~",
    (infectee_lag<=180&infectee.vac_dose=="2 dose")~"~180",
    (infectee_lag>180&infectee_lag<=365&infectee.vac_dose=="2 dose")~"180~365"
  ))
#
#
data_young_symptom_lag=data_young_symptom %>%
  filter(infectee.vac_dose%in%c("1 dose","2 dose")) %>%
  mutate(infectee.lag.new=case_when(
    ((infectee_lag>365&infectee.vac_dose=="2 dose")|infectee.vac_dose=="1 dose")~"~365",
    (infectee_lag<=365&infectee.vac_dose=="2 dose")~"~365"
  ))
#
##Bayesian generalized linear model
M3_lag <- bayesglm(symptom.state ~ infectee.lag.new
                   +infectee.age
                   +infectee.sex
                   +infectee.expose_date
                   +infector.vac_dose
                   +contact_setting, 
                   family=binomial(link="logit"),
                   data = data_young_symptom_lag, 
                   prior.scale=Inf, prior.df=Inf))
#
summary(M3_lag)

















