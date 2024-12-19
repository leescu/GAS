# ++++++++++++++++++++++++++++++++++++++++++++++++++++++ ####
# >>>>> Section 0. Packages and Functions used <<<<< ####
setwd("D:/PD_MDD/step 0. Data")
{#* section 0 Packages ####
  library(caret)
  library(car)
  library(cmprsk)
  library(dplyr)
  library(foreign)
  library(ggplot2)
  library(ggsci)
  library(ggrepel)
  library("ggthemes")
  library(lava)
  library(Matching)
  library(mediation)
  library(mice)
  library(pec)
  #install.packages("poLCA", dependencies = TRUE)
  library(poLCA)
  library(plyr)
  library(prodlim)
  library(reshape2)
  library(rms)
  library(riskRegression)
  library(survey)
  library(scales)
  library(survminer)
  library(survival)
  library(splines)
  library(timeROC)
  library(tableone)
  library(rms)
  library(withr)
  library(dplyr)
  library(doParallel)
}
# +================================================+ ####
# +====Section 0. Data arrangement=================+ ####
# +================================================+ #### 

load(file="2. Clinical/Interpolation_data_cross.Rdata")
#Interpolation_weighted<-subset(Interpolation_weighted,Race_ethnicity=="Non-Hispanic White")
Interpolation_data_cross<-Interpolation_data_cross[!is.na (Interpolation_data_cross$MDD),]
table(Interpolation_data_cross$Periodontal_disease)
table(Interpolation_data_cross$MDD)
table(Interpolation_data_cross$Periodontal_disease,Interpolation_data_cross$MDD)
Interpolation_data_cross$Periodontal_disease<-as.character(Interpolation_data_cross$Periodontal_disease)

Interpolation_data_cross$Periodontal_disease<-as.numeric(Interpolation_data_cross$Periodontal_disease)
table(Interpolation_data_cross$Periodontal_disease,Interpolation_data_cross$MDD)
colnames(Interpolation_data_cross)

colnames(Interpolation_data_cross)
#* Propensity Score Matching #####
#** moderately severe depression ##### 
#*** Model 0  ##### 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++ ####
# >>>>> section 1 PSM with PD  ####
#*** Model 1 Age_status,Sex,Race_ethnicity,SES,Marital_status ##### 
set.seed(0123)
psModel <- glm(formula =Periodontal_disease~Age+Sex+Ethnicity+Ethnicity+TDI_quantile
               +Smoke+Alcohol+PA+BMI_status+Diabetes+Hypertension,data =Interpolation_data_cross ,family = quasibinomial())
Interpolation_data_cross$pRhc <- predict(psModel, type = "response")
Interpolation_data_cross$pNoRhc <- 1 - Interpolation_data_cross$pRhc
Interpolation_data_cross$pAssign <- NA
Interpolation_data_cross$pAssign[Interpolation_data_cross$Periodontal_disease==1]<-Interpolation_data_cross$pRhc[Interpolation_data_cross$Periodontal_disease==1]
Interpolation_data_cross$pAssign[Interpolation_data_cross$Periodontal_disease==0]<-Interpolation_data_cross$pNoRhc[Interpolation_data_cross$Periodontal_disease==0]
Interpolation_data_cross$pMin <- pmin(Interpolation_data_cross$pRhc, Interpolation_data_cross$pNoRhc)
listMatch <- Match(Tr = (Interpolation_data_cross$Periodontal_disease == 1), X= log(Interpolation_data_cross$pRhc / Interpolation_data_cross$pNoRhc),M= 1,caliper  = 0.2,
                   replace  = FALSE,ties= TRUE,version  = "fast")
mb <- MatchBalance(psModel$formula, data=Interpolation_data_cross, match.out=listMatch, nboots=50)
Interpolation_MDD<- Interpolation_data_cross[unlist(listMatch[c("index.treated","index.control")]), ]
Interpolation_MDD_Cross<-Interpolation_MDD
table(Interpolation_MDD_Cross$MDD,Interpolation_MDD_Cross$Periodontal_disease)
table(Interpolation_MDD_Cross$MDD)
table(Interpolation_MDD_Cross$Periodontal_disease)
save(Interpolation_MDD_Cross, file="2. Clinical/Interpolation_MDD_Cross.Rdata")
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++ ####
# >>>>> section 1 PSM with MDD  ####
Interpolation_data_cross<-na.omit(Interpolation_data_cross)
set.seed(0123)
psModel <- glm(formula =MDD~Age+Sex+Ethnicity+Ethnicity+TDI_quantile
               +Smoke+Alcohol+PA+BMI_status+Diabetes+Hypertension,data =Interpolation_data_cross ,family = quasibinomial())
Interpolation_data_cross$pRhc <- predict(psModel, type = "response")
Interpolation_data_cross$pNoRhc <- 1 - Interpolation_data_cross$pRhc
Interpolation_data_cross$pAssign <- NA
Interpolation_data_cross$pAssign[Interpolation_data_cross$MDD==1]<-Interpolation_data_cross$pRhc[Interpolation_data_cross$MDD==1]
Interpolation_data_cross$pAssign[Interpolation_data_cross$MDD==0]<-Interpolation_data_cross$pNoRhc[Interpolation_data_cross$MDD==0]
Interpolation_data_cross$pMin <- pmin(Interpolation_data_cross$pRhc, Interpolation_data_cross$pNoRhc)
listMatch <- Match(Tr = (Interpolation_data_cross$MDD == 1), X= log(Interpolation_data_cross$pRhc / Interpolation_data_cross$pNoRhc),M= 1,caliper  = 0.2,
                   replace  = FALSE,ties= TRUE,version  = "fast")
mb <- MatchBalance(psModel$formula, data=Interpolation_data_cross, match.out=listMatch, nboots=50)
Interpolation_PD<- Interpolation_data_cross[unlist(listMatch[c("index.treated","index.control")]), ]
Interpolation_PD_Cross<-Interpolation_PD
table(Interpolation_PD_Cross$MDD,Interpolation_PD_Cross$Periodontal_disease)
table(Interpolation_PD_Cross$MDD)
table(Interpolation_PD_Cross$Periodontal_disease)

save(Interpolation_PD_Cross, file="2. Clinical/Interpolation_PD_Cross.Rdata")

# +================================================+ ####
# +====Section x. Corss sectional orginal==========+ ####
# +================================================+ #### 
load(file="2. Clinical/Interpolation_data_cross.Rdata")
Interpolation_data_cross<-na.omit(Interpolation_data_cross)
Interpolation_data_cross$Periodontal_disease<-as.character(Interpolation_data_cross$Periodontal_disease)

Interpolation_data_cross$Periodontal_disease<-as.numeric(Interpolation_data_cross$Periodontal_disease)
table(Interpolation_data_cross$Periodontal_disease)
table(Interpolation_data_cross$MDD)
colnames(Interpolation_data_cross)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++ ####
# >>>>> section 1 PD to MDD  ####
#*** Model 0  ##### 
MDD <- glm(MDD~Periodontal_disease,
           data=Interpolation_data_cross, family = binomial)
SUM<-summary(MDD)[["coefficients"]]
OR<-round(exp(SUM[2,1]),3)
CI5<-round(exp(SUM[2,1]-1.96*SUM[2,2]),3)
CI95<-round(exp(SUM[2,1]+1.96*SUM[2,2]),3)
P<-SUM[2,4]
result1 <- data.frame('OR'=OR,'lower .95'=CI5,'upper .95'=CI95,
                      'P value' =P,'model'="Model 0",'status'="PD to MDD cross")
#*** Model 1  ##### 
MDD <- glm(MDD~Periodontal_disease+
           Age+Sex+Ethnicity+Education+TDI_quantile,
           data=Interpolation_data_cross,
           family = binomial)
SUM<-summary(MDD)[["coefficients"]]
OR<-round(exp(SUM[2,1]),3)
CI5<-round(exp(SUM[2,1]-1.96*SUM[2,2]),3)
CI95<-round(exp(SUM[2,1]+1.96*SUM[2,2]),3)
P<-SUM[2,4]
result2 <- data.frame('OR'=OR,'lower .95'=CI5,'upper .95'=CI95,
                      'P value' =P,'model'="Model 1",'status'="PD to MDD cross")

#*** Model 2  ##### 
MDD <- glm(MDD~Periodontal_disease+
             Age+Sex+Ethnicity+Education+TDI_quantile+
             Smoke+Alcohol+PA
             ,data=Interpolation_data_cross,
           family = binomial)
SUM<-summary(MDD)[["coefficients"]]
OR<-round(exp(SUM[2,1]),3)
CI5<-round(exp(SUM[2,1]-1.96*SUM[2,2]),3)
CI95<-round(exp(SUM[2,1]+1.96*SUM[2,2]),3)
P<-SUM[2,4]
result3 <- data.frame('OR'=OR,'lower .95'=CI5,'upper .95'=CI95,
                      'P value' =P,'model'="Model 2",'status'="PD to MDD cross")

#*** Model 3  ##### 
MDD <- glm(MDD~Periodontal_disease+
             Age+Sex+Ethnicity+Education+TDI_quantile+
             Smoke+Alcohol+PA+
             BMI_status+Diabetes+Hypertension+Medication
           ,data=Interpolation_data_cross,
           family = binomial)
SUM<-summary(MDD)[["coefficients"]]
OR<-round(exp(SUM[2,1]),3)
CI5<-round(exp(SUM[2,1]-1.96*SUM[2,2]),3)
CI95<-round(exp(SUM[2,1]+1.96*SUM[2,2]),3)
P<-SUM[2,4]
result4 <- data.frame('OR'=OR,'lower .95'=CI5,'upper .95'=CI95,
                      'P value' =P,'model'="Model 3",'status'="PD to MDD cross")
result_cross_original_PD<-rbind(result1,result2,result3,result4)

result_cross_original_PD
save(result_cross_original_PD, file="2. Clinical/result_cross_original_PD.Rdata")
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++ ####
# >>>>> section2 MDD to PD  ####
PD <- glm(Periodontal_disease~MDD,
           data=Interpolation_data_cross, family = binomial)
SUM<-summary(MDD)[["coefficients"]]
OR<-round(exp(SUM[2,1]),3)
CI5<-round(exp(SUM[2,1]-1.96*SUM[2,2]),3)
CI95<-round(exp(SUM[2,1]+1.96*SUM[2,2]),3)
P<-SUM[2,4]
result0 <- data.frame('OR'=OR,'lower .95'=CI5,'upper .95'=CI95,
                      'P value' =P,'model'="Model 0",'status'="MDD to PD cross")
#*** Model 1  ##### 
MDD <- glm(Periodontal_disease~MDD+
             Age+Sex+Ethnicity+Education+TDI_quantile,
           data=Interpolation_data_cross,
           family = binomial)
SUM<-summary(MDD)[["coefficients"]]
OR<-round(exp(SUM[2,1]),3)
CI5<-round(exp(SUM[2,1]-1.96*SUM[2,2]),3)
CI95<-round(exp(SUM[2,1]+1.96*SUM[2,2]),3)
P<-SUM[2,4]
result1 <- data.frame('OR'=OR,'lower .95'=CI5,'upper .95'=CI95,
                      'P value' =P,'model'="Model 1",'status'="MDD to PD cross")

#*** Model 2  ##### 
MDD <- glm(Periodontal_disease~MDD+
             Age+Sex+Ethnicity+Education+TDI_quantile+
             Smoke+Alcohol+PA
           ,data=Interpolation_data_cross,
           family = binomial)
SUM<-summary(MDD)[["coefficients"]]
OR<-round(exp(SUM[2,1]),3)
CI5<-round(exp(SUM[2,1]-1.96*SUM[2,2]),3)
CI95<-round(exp(SUM[2,1]+1.96*SUM[2,2]),3)
P<-SUM[2,4]
result2 <- data.frame('OR'=OR,'lower .95'=CI5,'upper .95'=CI95,
                      'P value' =P,'model'="Model 2",'status'="MDD to PD cross")

#*** Model 3  ##### 
MDD <- glm(Periodontal_disease~MDD+
             Age+Sex+Ethnicity+Education+TDI_quantile+
             Smoke+Alcohol+PA+
             BMI_status+Diabetes+Hypertension+Medication
           ,data=Interpolation_data_cross,
           family = binomial)
SUM<-summary(MDD)[["coefficients"]]
OR<-round(exp(SUM[2,1]),3)
CI5<-round(exp(SUM[2,1]-1.96*SUM[2,2]),3)
CI95<-round(exp(SUM[2,1]+1.96*SUM[2,2]),3)
P<-SUM[2,4]
result3 <- data.frame('OR'=OR,'lower .95'=CI5,'upper .95'=CI95,
                      'P value' =P,'model'="Model 3",'status'="MDD to PD cross")
result_cross_original_MDD<-rbind(result0,result1,result2,result3)
result_cross_original_MDD
save(result_cross_original_MDD, file="2. Clinical/result_cross_original_MDD.Rdata")

# +================================================+ ####
# +====Section 0. Cohort===========================+ ####
# +================================================+ #### 
# >>>>> Total  ####
load(file="2. Clinical/Interpolation_data_cohort.Rdata")
#Interpolation_weighted<-subset(Interpolation_weighted,Race_ethnicity=="Non-Hispanic White")
table(Interpolation_data_cohort$MDD,Interpolation_data_cohort$Periodontal_disease)
table(Interpolation_data_cohort$Periodontal_disease)
table(Interpolation_data_cohort$MDD)
Interpolation_data_cohort$Periodontal_disease<-as.character(Interpolation_data_cohort$Periodontal_disease)

Interpolation_data_cohort$Periodontal_disease<-as.numeric(Interpolation_data_cohort$Periodontal_disease)
table(Interpolation_data_cohort$Periodontal_disease)
table(Interpolation_data_cohort$MDD)
colnames(Interpolation_data_cohort)
#* Propensity Score Matching #####
#** moderately severe depression ##### 
#*** Model 0  ##### 
MDD <- coxph(Surv(days,MDD)~Periodontal_disease
             ,data=Interpolation_data_cohort)
model_result<-summary(MDD)
P<-model_result[["coefficients"]][1,"Pr(>|z|)"]
HR<-as.numeric(model_result[["conf.int"]][1,c("exp(coef)","lower .95","upper .95")])
result0 <- data.frame('HR'=HR[1],'lower .95'=HR[2],'upper .95'=HR[3],
                      'P value' =P,'model'="Model 0",'status'="Overall MDD Risk")

MDD <- coxph(Surv(days,MDD)~Periodontal_disease+
              Age+Sex+Ethnicity+Education+TDI_quantile
             ,data=Interpolation_data_cohort)
model_result<-summary(MDD)
P<-model_result[["coefficients"]][1,"Pr(>|z|)"]
HR<-as.numeric(model_result[["conf.int"]][1,c("exp(coef)","lower .95","upper .95")])
result1 <- data.frame('HR'=HR[1],'lower .95'=HR[2],'upper .95'=HR[3],
                      'P value' =P,'model'="Model 1",'status'="Overall MDD Risk")

MDD <- coxph(Surv(days,MDD)~Periodontal_disease+
               Age+Sex+Ethnicity+Education+TDI_quantile+
               Smoke+Alcohol+PA
             ,data=Interpolation_data_cohort)
model_result<-summary(MDD)
P<-model_result[["coefficients"]][1,"Pr(>|z|)"]
HR<-as.numeric(model_result[["conf.int"]][1,c("exp(coef)","lower .95","upper .95")])
result2 <- data.frame('HR'=HR[1],'lower .95'=HR[2],'upper .95'=HR[3],
                      'P value' =P,'model'="Model 2",'status'="Overall MDD Risk")

MDD <- coxph(Surv(days,MDD)~Periodontal_disease+
                               Age+Sex+Ethnicity+Education+TDI_quantile+
                               Smoke+Alcohol+PA+
                               BMI_status+Diabetes+Hypertension+Medication
             ,data=Interpolation_data_cohort)
model_result<-summary(MDD)
P<-model_result[["coefficients"]][1,"Pr(>|z|)"]
HR<-as.numeric(model_result[["conf.int"]][1,c("exp(coef)","lower .95","upper .95")])
result3 <- data.frame('HR'=HR[1],'lower .95'=HR[2],'upper .95'=HR[3],
                                   'P value' =P,'model'="Model 3",'status'="Overall MDD Risk")

result_cohort<-rbind(result0,result1,result2,result3)
result_cohort
save(result_cohort, file="15. Cohort/result_cohort.Rdata")
# >>>>> 5 years  ####
load(file="15. Cohort/Interpolation_data_cohort.Rdata")

#Interpolation_weighted<-subset(Interpolation_weighted,Race_ethnicity=="Non-Hispanic White")
table(Interpolation_data_cohort$MDD,Interpolation_data_cohort$Periodontal_disease)
table(Interpolation_data_cohort$Periodontal_disease)
table(Interpolation_data_cohort$MDD)
Interpolation_data_cohort$Periodontal_disease<-as.character(Interpolation_data_cohort$Periodontal_disease)

Interpolation_data_cohort$Periodontal_disease<-as.numeric(Interpolation_data_cohort$Periodontal_disease)
table(Interpolation_data_cohort$Periodontal_disease)
table(Interpolation_data_cohort$MDD)
colnames(Interpolation_data_cohort)

Interpolation_data_cohort$days_5<- pmin(Interpolation_data_cohort$days, 5*365)  # 如果时间超过5年，设为5年
Interpolation_data_cohort$MDD_5 <- ifelse(Interpolation_data_cohort$days > 5*365, 0, Interpolation_data_cohort$MDD)  # 超过5年视为无事件

#*** Model 0  ##### 
MDD <- coxph(Surv(days_5,MDD_5)~Periodontal_disease
             ,data=Interpolation_data_cohort)
model_result<-summary(MDD)
P<-model_result[["coefficients"]][1,"Pr(>|z|)"]
HR<-as.numeric(model_result[["conf.int"]][1,c("exp(coef)","lower .95","upper .95")])
result0 <- data.frame('HR'=HR[1],'lower .95'=HR[2],'upper .95'=HR[3],
                      'P value' =P,'model'="Model 0",'status'="5-years MDD Risk")

MDD <- coxph(Surv(days_5,MDD_5)~Periodontal_disease+
               Age+Sex+Ethnicity+Education+TDI_quantile
             ,data=Interpolation_data_cohort)
model_result<-summary(MDD)
P<-model_result[["coefficients"]][1,"Pr(>|z|)"]
HR<-as.numeric(model_result[["conf.int"]][1,c("exp(coef)","lower .95","upper .95")])
result1 <- data.frame('HR'=HR[1],'lower .95'=HR[2],'upper .95'=HR[3],
                      'P value' =P,'model'="Model 1",'status'="5-years MDD Risk")

MDD <- coxph(Surv(days_5,MDD_5)~Periodontal_disease+
               Age+Sex+Ethnicity+Education+TDI_quantile+
               Smoke+Alcohol+PA
             ,data=Interpolation_data_cohort)
model_result<-summary(MDD)
P<-model_result[["coefficients"]][1,"Pr(>|z|)"]
HR<-as.numeric(model_result[["conf.int"]][1,c("exp(coef)","lower .95","upper .95")])
result2 <- data.frame('HR'=HR[1],'lower .95'=HR[2],'upper .95'=HR[3],
                      'P value' =P,'model'="Model 2",'status'="5-years MDD Risk")

MDD <- coxph(Surv(days_5,MDD_5)~Periodontal_disease+
               Age+Sex+Ethnicity+Education+TDI_quantile+
               Smoke+Alcohol+PA+
               BMI_status+Diabetes+Hypertension+Medication
             ,data=Interpolation_data_cohort)
model_result<-summary(MDD)
P<-model_result[["coefficients"]][1,"Pr(>|z|)"]
HR<-as.numeric(model_result[["conf.int"]][1,c("exp(coef)","lower .95","upper .95")])
result3 <- data.frame('HR'=HR[1],'lower .95'=HR[2],'upper .95'=HR[3],
                      'P value' =P,'model'="Model 3",'status'="5-years MDD Risk")

result_cohort_5<-rbind(result0,result1,result2,result3)
result_cohort_5
save(result_cohort_5, file="15. Cohort/result_cohort_5.Rdata")

# >>>>> 10 years  ####
load(file="15. Cohort/Interpolation_data_cohort.Rdata")

#Interpolation_weighted<-subset(Interpolation_weighted,Race_ethnicity=="Non-Hispanic White")
table(Interpolation_data_cohort$MDD,Interpolation_data_cohort$Periodontal_disease)
table(Interpolation_data_cohort$Periodontal_disease)
table(Interpolation_data_cohort$MDD)
Interpolation_data_cohort$Periodontal_disease<-as.character(Interpolation_data_cohort$Periodontal_disease)

Interpolation_data_cohort$Periodontal_disease<-as.numeric(Interpolation_data_cohort$Periodontal_disease)
table(Interpolation_data_cohort$Periodontal_disease)
table(Interpolation_data_cohort$MDD)
colnames(Interpolation_data_cohort)

Interpolation_data_cohort$days_10<- pmin(Interpolation_data_cohort$days, 10*365)  # 如果时间超过5年，设为5年
Interpolation_data_cohort$MDD_10 <- ifelse(Interpolation_data_cohort$days > 10*365, 0, Interpolation_data_cohort$MDD)  # 超过5年视为无事件
#* Propensity Score Matching #####
#** moderately severe depression ##### 
#*** Model 0  ##### 
MDD <- coxph(Surv(days_10,MDD_10)~Periodontal_disease
             ,data=Interpolation_data_cohort)
model_result<-summary(MDD)
P<-model_result[["coefficients"]][1,"Pr(>|z|)"]
HR<-as.numeric(model_result[["conf.int"]][1,c("exp(coef)","lower .95","upper .95")])
result0 <- data.frame('HR'=HR[1],'lower .95'=HR[2],'upper .95'=HR[3],
                      'P value' =P,'model'="Model 0",'status'="10-years MDD Risk")

MDD <- coxph(Surv(days_10,MDD_10)~Periodontal_disease+
               Age+Sex+Ethnicity+Education+TDI_quantile
             ,data=Interpolation_data_cohort)
model_result<-summary(MDD)
P<-model_result[["coefficients"]][1,"Pr(>|z|)"]
HR<-as.numeric(model_result[["conf.int"]][1,c("exp(coef)","lower .95","upper .95")])
result1 <- data.frame('HR'=HR[1],'lower .95'=HR[2],'upper .95'=HR[3],
                      'P value' =P,'model'="Model 1",'status'="10-years MDD Risk")

MDD <- coxph(Surv(days_10,MDD_10)~Periodontal_disease+
               Age+Sex+Ethnicity+Education+TDI_quantile+
               Smoke+Alcohol+PA
             ,data=Interpolation_data_cohort)
model_result<-summary(MDD)
P<-model_result[["coefficients"]][1,"Pr(>|z|)"]
HR<-as.numeric(model_result[["conf.int"]][1,c("exp(coef)","lower .95","upper .95")])
result2 <- data.frame('HR'=HR[1],'lower .95'=HR[2],'upper .95'=HR[3],
                      'P value' =P,'model'="Model 2",'status'="10-years MDD Risk")

MDD <- coxph(Surv(days_10,MDD_10)~Periodontal_disease+
               Age+Sex+Ethnicity+Education+TDI_quantile+
               Smoke+Alcohol+PA+
               BMI_status+Diabetes+Hypertension+Medication
             ,data=Interpolation_data_cohort)
model_result<-summary(MDD)
P<-model_result[["coefficients"]][1,"Pr(>|z|)"]
HR<-as.numeric(model_result[["conf.int"]][1,c("exp(coef)","lower .95","upper .95")])
result3 <- data.frame('HR'=HR[1],'lower .95'=HR[2],'upper .95'=HR[3],
                      'P value' =P,'model'="Model 3",'status'="10-years MDD Risk")

result_cohort_10<-rbind(result0,result1,result2,result3)
result_cohort_10
save(result_cohort_10, file="15. Cohort/result_cohort_10.Rdata")
# +================================================+ ####
# +====Section x. Corss sectional PD PSM===========+ ####
# +================================================+ #### 
load(file="2. Clinical/Interpolation_MDD_Cross.Rdata")
Interpolation_MDD_Cross<-na.omit(Interpolation_MDD_Cross)
Interpolation_MDD_Cross$Periodontal_disease<-as.character(Interpolation_MDD_Cross$Periodontal_disease)

Interpolation_MDD_Cross$Periodontal_disease<-as.numeric(Interpolation_MDD_Cross$Periodontal_disease)
table(Interpolation_MDD_Cross$Periodontal_disease)
table(Interpolation_MDD_Cross$MDD)
colnames(Interpolation_MDD_Cross)
#* Propensity Score Matching #####
#** moderately severe depression ##### 
#*** Model 0  ##### 
MDD <- glm(MDD~Periodontal_disease,
           data=Interpolation_MDD_Cross, family = binomial)
SUM<-summary(MDD)[["coefficients"]]
OR<-round(exp(SUM[2,1]),3)
CI5<-round(exp(SUM[2,1]-1.96*SUM[2,2]),3)
CI95<-round(exp(SUM[2,1]+1.96*SUM[2,2]),3)
P<-SUM[2,4]
result0 <- data.frame('OR'=OR,'lower .95'=CI5,'upper .95'=CI95,
                      'P value' =P,'model'="Model 0",'status'="PSM PD Cross")
#*** Model 1  ##### 
MDD <- glm(MDD~Periodontal_disease+
             Age+Sex+Ethnicity+Education+TDI_quantile,
           data=Interpolation_MDD_Cross,
           family = binomial)
SUM<-summary(MDD)[["coefficients"]]
OR<-round(exp(SUM[2,1]),3)
CI5<-round(exp(SUM[2,1]-1.96*SUM[2,2]),3)
CI95<-round(exp(SUM[2,1]+1.96*SUM[2,2]),3)
P<-SUM[2,4]
result1 <- data.frame('OR'=OR,'lower .95'=CI5,'upper .95'=CI95,
                      'P value' =P,'model'="Model 1",'status'="PSM PD Cross")

#*** Model 2  ##### 
MDD <- glm(MDD~Periodontal_disease+
             Age+Sex+Ethnicity+Education+TDI_quantile+
             Smoke+Alcohol+PA
           ,data=Interpolation_MDD_Cross,
           family = binomial)
SUM<-summary(MDD)[["coefficients"]]
OR<-round(exp(SUM[2,1]),3)
CI5<-round(exp(SUM[2,1]-1.96*SUM[2,2]),3)
CI95<-round(exp(SUM[2,1]+1.96*SUM[2,2]),3)
P<-SUM[2,4]
result2 <- data.frame('OR'=OR,'lower .95'=CI5,'upper .95'=CI95,
                      'P value' =P,'model'="Model 2",'status'="PSM PD Cross")

#*** Model 3  ##### 
MDD <- glm(MDD~Periodontal_disease+
             Age+Sex+Ethnicity+Education+TDI_quantile+
             Smoke+Alcohol+PA+
             BMI_status+Diabetes+Hypertension+Medication
           ,data=Interpolation_MDD_Cross,
           family = binomial)
SUM<-summary(MDD)[["coefficients"]]
OR<-round(exp(SUM[2,1]),3)
CI5<-round(exp(SUM[2,1]-1.96*SUM[2,2]),3)
CI95<-round(exp(SUM[2,1]+1.96*SUM[2,2]),3)
P<-SUM[2,4]
result3 <- data.frame('OR'=OR,'lower .95'=CI5,'upper .95'=CI95,
                      'P value' =P,'model'="Model 3",'status'="PSM PD Cross")
result_cross_PSM_PD<-rbind(result0,result1,result2,result3)
result_cross_PSM_PD
save(result_cross_PSM_PD, file="2. Clinical/result_cross_PSM_PD.Rdata")

# +================================================+ ####
# +====Section x. Corss sectional PD PSM===========+ ####
# +================================================+ #### 
load(file="2. Clinical/Interpolation_PD_Cross.Rdata")
Interpolation_PD_Cross<-na.omit(Interpolation_PD_Cross)
Interpolation_PD_Cross$Periodontal_disease<-as.character(Interpolation_PD_Cross$Periodontal_disease)

Interpolation_PD_Cross$Periodontal_disease<-as.numeric(Interpolation_PD_Cross$Periodontal_disease)
table(Interpolation_PD_Cross$Periodontal_disease)
table(Interpolation_PD_Cross$MDD)
colnames(Interpolation_PD_Cross)
#* Propensity Score Matching #####
#** moderately severe depression ##### 
#*** Model 0  ##### 
MDD <- glm(Periodontal_disease~MDD,
           data=Interpolation_PD_Cross, family = binomial)
SUM<-summary(MDD)[["coefficients"]]
OR<-round(exp(SUM[2,1]),3)
CI5<-round(exp(SUM[2,1]-1.96*SUM[2,2]),3)
CI95<-round(exp(SUM[2,1]+1.96*SUM[2,2]),3)
P<-SUM[2,4]
result0 <- data.frame('OR'=OR,'lower .95'=CI5,'upper .95'=CI95,
                      'P value' =P,'model'="Model 0",'status'="PSM MDD Cross")
#*** Model 1  ##### 
MDD <- glm(Periodontal_disease~MDD+
             Age+Sex+Ethnicity+Education+TDI_quantile,
           data=Interpolation_PD_Cross,
           family = binomial)
SUM<-summary(MDD)[["coefficients"]]
OR<-round(exp(SUM[2,1]),3)
CI5<-round(exp(SUM[2,1]-1.96*SUM[2,2]),3)
CI95<-round(exp(SUM[2,1]+1.96*SUM[2,2]),3)
P<-SUM[2,4]
result1 <- data.frame('OR'=OR,'lower .95'=CI5,'upper .95'=CI95,
                      'P value' =P,'model'="Model 1",'status'="PSM MDD Cross")

#*** Model 2  ##### 
MDD <- glm(Periodontal_disease~MDD+
             Age+Sex+Ethnicity+Education+TDI_quantile+
             Smoke+Alcohol+PA
           ,data=Interpolation_PD_Cross,
           family = binomial)
SUM<-summary(MDD)[["coefficients"]]
OR<-round(exp(SUM[2,1]),3)
CI5<-round(exp(SUM[2,1]-1.96*SUM[2,2]),3)
CI95<-round(exp(SUM[2,1]+1.96*SUM[2,2]),3)
P<-SUM[2,4]
result2 <- data.frame('OR'=OR,'lower .95'=CI5,'upper .95'=CI95,
                      'P value' =P,'model'="Model 2",'status'="PSM MDD Cross")

#*** Model 3  ##### 
MDD <- glm(Periodontal_disease~MDD+
             Age+Sex+Ethnicity+Education+TDI_quantile+
             Smoke+Alcohol+PA+
             BMI_status+Diabetes+Hypertension+Medication
           ,data=Interpolation_PD_Cross,
           family = binomial)
SUM<-summary(MDD)[["coefficients"]]
OR<-round(exp(SUM[2,1]),3)
CI5<-round(exp(SUM[2,1]-1.96*SUM[2,2]),3)
CI95<-round(exp(SUM[2,1]+1.96*SUM[2,2]),3)
P<-SUM[2,4]
result3 <- data.frame('OR'=OR,'lower .95'=CI5,'upper .95'=CI95,
                      'P value' =P,'model'="Model 3",'status'="PSM MDD Cross")
result_cross_PSM_MDD<-rbind(result0,result1,result2,result3)
result_cross_PSM_MDD
save(result_cross_PSM_MDD, file="2. Clinical/result_cross_PSM_MDD.Rdata")



# +================================================+ ####
# +====Section x. Ggplot1===========================+ ####
# +================================================+ #### 
load(file="2. Clinical/result_cross_original_PD.Rdata")
load(file="2. Clinical/result_cross_original_MDD.Rdata")

colnames(result_cross_original_PD)[1]<-"effect"
colnames(result_cross_original_MDD)[1]<-"effect"


Result<-rbind(result_cross_original_PD,result_cross_original_MDD)

# 定义一个函数将数值格式化为4位有效数字，并补齐不足位数
format_numeric <- function(x) {
  sprintf("%.3f", x)
}

format_p_value <- function(x) {
  format(x, scientific = TRUE, digits = 4)
}

Result$effect <- sapply(Result$effect, format_numeric)
Result$lower..95 <- sapply(Result$lower..95, format_numeric)
Result$upper..95 <- sapply(Result$upper..95, format_numeric)
Result$P.value <- sapply(Result$P.value, format_p_value)



# 新增列：显示效应值及95%置信区间范围
Result$CI <- paste0(
  Result$effect, " (", Result$lower..95, " - ", Result$upper..95, ")"
)



Result$effect <- as.numeric(Result$effect)
Result$lower..95 <- as.numeric(Result$lower..95)
Result$upper..95 <- as.numeric(Result$upper..95)
Result$model <- factor(Result$model, levels = c('Model 4', 'Model 3', 'Model 2', 'Model 1', 'Model 0'))
#创建自定义标签函数
 custom_labels <- c(
   "PD to MDD cross" = "PD to MDD Risk (Cross)",
   "MDD to PD cross" = "MDD to PD Risk (Cross)")
Result$status <- factor(
  Result$status, 
  levels = c("PD to MDD cross", "MDD to PD cross")
)
save(Result,file="2. Clinical/cross_original_Result.Rdata")
# 在每个分面图添加的左上角标签数据
label_data <- data.frame(
  status = c("PD to MDD cross", "MDD to PD cross"),
  label = c("Odds Ratio (95%CI)", "Odds Ratio (95%CI)"),
  effect = 0.91,  # 确保文本在x轴的左侧边缘附近
  model = "Effect size"  # 将文本放置在分面图的顶部
)
label_data$status <- factor(
  label_data$status, 
  levels = c("PD to MDD cross", "MDD to PD cross"))
# 绘制森林图
p<-ggplot(Result, aes(x = effect, y = model, color = status)) +
  geom_errorbarh(aes(xmin = lower..95, xmax = upper..95), height = 0.2, color = "black", linewidth = 0.7) +  # 黑色误差棒
  geom_point(size = 3) +  # 绘制点，颜色根据 status 区分
  theme_minimal() +  # 使用简洁主题
  labs(x = "Effect",y = "Model") +  # 设置标签
  theme(
    panel.grid.major.x = element_blank(),  # 去除 X 轴主网格线
    panel.grid.minor.x = element_blank(),  # 去除 X 轴次网格线
    panel.grid.minor.y = element_blank(),  # 去除 Y 轴次网格线
    legend.position = "right",  # 图例在右侧
    legend.text = element_text(size = 12),  # 增大图例字体大小
    legend.title = element_text(size = 14, face = "bold"),  # 图例标题加粗并增大字体
    axis.text.x = element_text(size = 14, color = "black"),  # X 轴刻度字体设为黑色
    axis.text.y = element_text(size = 14, color = "black"),  # Y 轴刻度字体设为黑色
    axis.title = element_text(size = 14, face = "bold", color = "black"),  # 坐标轴标题字体加粗且为黑色
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),  # X 轴刻度线
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),  # Y 轴刻度线
    strip.text = element_text(size = 12, face = "bold"),  # 增大分面标签字体并加粗
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # 添加图的边框
  ) +
  scale_color_manual(values = c("#4865a9", "#ef8a43")) +  # 分组颜色设置
  scale_x_continuous(limits = c(0.9, 1.7), breaks = seq(0.9, 1.7, by = 0.1)) +  # 设置 X 轴范围和刻度
  facet_grid(status ~ ., scales = "free_y", space = "free", labeller = as_labeller(custom_labels)) + # 自定义分面标签
  annotate("segment", x = 1, xend = 1, y = 0.5, yend = 4.5, color = "#e53b2c", linewidth = 1) +  # 部分红色垂直线
  geom_segment(data = Result, aes(x = 0.9, xend = 1.7, y = 5, yend = 5), 
               inherit.aes = FALSE, color = "white", linewidth = 0.5) +
  geom_text(data = label_data, aes(x = effect, y = model, label = label), 
            inherit.aes = FALSE, color = "black", size = 4,fontface = "bold",  hjust = 0)
#6*7.5

ggsave(
  filename = "2. Clinical/forest plot_cross.pdf",  # 文件名
  plot = p,  # 图形对象
  width = 7.5,  # 宽度 (英寸)
  height = 4.5,  # 高度 (英寸)
  dpi = 300,  # 分辨率设置为 300 DPI
  device = cairo_pdf  # 确保 PDF 的高质量输出
)




# +================================================+ ####
# +====Section x. Ggplot2===========================+ ####
# +================================================+ #### 
# +================================================+ ####
# +====Section x. Ggplot1===========================+ ####
# +================================================+ #### 

load(file="2. Clinical/result_cross_PSM_PD.Rdata")
load(file="2. Clinical/result_cross_PSM_MDD.Rdata")

colnames(result_cross_PSM_PD)[1]<-"effect"
colnames(result_cross_PSM_MDD)[1]<-"effect"


Result<-rbind(
              result_cross_PSM_PD,result_cross_PSM_MDD)

# 定义一个函数将数值格式化为4位有效数字，并补齐不足位数
format_numeric <- function(x) {
  sprintf("%.3f", x)
}

format_p_value <- function(x) {
  format(x, scientific = TRUE, digits = 4)
}

Result$effect <- sapply(Result$effect, format_numeric)
Result$lower..95 <- sapply(Result$lower..95, format_numeric)
Result$upper..95 <- sapply(Result$upper..95, format_numeric)
Result$P.value <- sapply(Result$P.value, format_p_value)



# 新增列：显示效应值及95%置信区间范围
Result$CI <- paste0(
  Result$effect, " (", Result$lower..95, " - ", Result$upper..95, ")"
)



Result$effect <- as.numeric(Result$effect)
Result$lower..95 <- as.numeric(Result$lower..95)
Result$upper..95 <- as.numeric(Result$upper..95)
Result$model <- factor(Result$model, levels = c('Model 4', 'Model 3', 'Model 2', 'Model 1', 'Model 0'))
#创建自定义标签函数
custom_labels <- c(
  "PSM PD Cross" = "PD to MDD Risk (PSM)",
  "PSM MDD Cross" = "MDD to PD Risk (PSM)")
Result$status <- factor(
  Result$status, 
  levels = c("PSM PD Cross", "PSM MDD Cross"
  )
)
save(Result,file="2. Clinical/cross_psm_Result.Rdata")
# 在每个分面图添加的左上角标签数据
label_data <- data.frame(
  status = c(
             "PSM PD Cross", "PSM MDD Cross"),
  label = c(
            "Odds Ratio (95%CI)", "Odds Ratio (95%CI)"),
  effect = 0.91,  # 确保文本在x轴的左侧边缘附近
  model = "Effect size"  # 将文本放置在分面图的顶部
)
label_data$status <- factor(
  label_data$status, 
  levels = c(
             "PSM PD Cross", "PSM MDD Cross"))
# 绘制森林图
p<-ggplot(Result, aes(x = effect, y = model, color = status)) +
  geom_errorbarh(aes(xmin = lower..95, xmax = upper..95), height = 0.2, color = "black", linewidth = 0.7) +  # 黑色误差棒
  geom_point(size = 3) +  # 绘制点，颜色根据 status 区分
  theme_minimal() +  # 使用简洁主题
  labs(x = "Effect",y = "Model") +  # 设置标签
  theme(
    panel.grid.major.x = element_blank(),  # 去除 X 轴主网格线
    panel.grid.minor.x = element_blank(),  # 去除 X 轴次网格线
    panel.grid.minor.y = element_blank(),  # 去除 Y 轴次网格线
    legend.position = "right",  # 图例在右侧
    legend.text = element_text(size = 12),  # 增大图例字体大小
    legend.title = element_text(size = 14, face = "bold"),  # 图例标题加粗并增大字体
    axis.text.x = element_text(size = 14, color = "black"),  # X 轴刻度字体设为黑色
    axis.text.y = element_text(size = 14, color = "black"),  # Y 轴刻度字体设为黑色
    axis.title = element_text(size = 14, face = "bold", color = "black"),  # 坐标轴标题字体加粗且为黑色
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),  # X 轴刻度线
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),  # Y 轴刻度线
    strip.text = element_text(size = 12, face = "bold"),  # 增大分面标签字体并加粗
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # 添加图的边框
  ) +
  scale_color_manual(values = c("#4865a9", "#ef8a43")) +  # 分组颜色设置
  scale_x_continuous(limits = c(0.9, 1.7), breaks = seq(0.9, 1.7, by = 0.1)) +  # 设置 X 轴范围和刻度
  facet_grid(status ~ ., scales = "free_y", space = "free", labeller = as_labeller(custom_labels)) + # 自定义分面标签
  annotate("segment", x = 1, xend = 1, y = 0.5, yend = 4.5, color = "#e53b2c", linewidth = 1) +  # 部分红色垂直线
  geom_segment(data = Result, aes(x = 0.9, xend = 1.7, y = 5, yend = 5), 
               inherit.aes = FALSE, color = "white", linewidth = 0.5) +
  geom_text(data = label_data, aes(x = effect, y = model, label = label), 
            inherit.aes = FALSE, color = "black", size = 4,fontface = "bold",  hjust = 0)
#6*7.5

ggsave(
  filename = "2. Clinical/forest plot psm.pdf",  # 文件名
  plot = p,  # 图形对象
  width = 7.5,  # 宽度 (英寸)
  height = 4.5,  # 高度 (英寸)
  dpi = 300,  # 分辨率设置为 300 DPI
  device = cairo_pdf  # 确保 PDF 的高质量输出
)




# +================================================+ ####
# +====Section x. Ggplot1===========================+ ####
# +================================================+ #### 

load(file="2. Clinical/result_cohort.Rdata")
load(file="2. Clinical/result_cohort_5.Rdata")
load(file="2. Clinical/result_cohort_10.Rdata")

colnames(result_cohort)[1]<-"effect"
colnames(result_cohort_5)[1]<-"effect"
colnames(result_cohort_10)[1]<-"effect"


Result<-rbind(
              result_cohort_5,result_cohort_10,result_cohort)

# 定义一个函数将数值格式化为4位有效数字，并补齐不足位数
format_numeric <- function(x) {
  sprintf("%.3f", x)
}

format_p_value <- function(x) {
  format(x, scientific = TRUE, digits = 4)
}

Result$effect <- sapply(Result$effect, format_numeric)
Result$lower..95 <- sapply(Result$lower..95, format_numeric)
Result$upper..95 <- sapply(Result$upper..95, format_numeric)
Result$P.value <- sapply(Result$P.value, format_p_value)



# 新增列：显示效应值及95%置信区间范围
Result$CI <- paste0(
  Result$effect, " (", Result$lower..95, " - ", Result$upper..95, ")"
)


Result$effect <- as.numeric(Result$effect)
Result$lower..95 <- as.numeric(Result$lower..95)
Result$upper..95 <- as.numeric(Result$upper..95)
Result$model <- factor(Result$model, levels = c('Model 4', 'Model 3', 'Model 2', 'Model 1', 'Model 0'))
#创建自定义标签函数
custom_labels <- c(
  "5-years MDD Risk" = "5-years MDD Risk",
  "10-years MDD Risk" = "10-years MDD Risk",
  "Overall MDD Risk" = "Overall MDD Risk")
Result$status <- factor(
  Result$status, 
  levels = c(
             "5-years MDD Risk","10-years MDD Risk" , "Overall MDD Risk"
  )
)
save(Result,file="2. Clinical/Cohort_Result.Rdata")
# 在每个分面图添加的左上角标签数据
label_data <- data.frame(
  status = c(
             "5-years MDD Risk","10-years MDD Risk" , "Overall MDD Risk"),
  label = c(
            "Hazard Ratio (95%CI)","Hazard Ratio (95%CI)","Hazard Ratio (95%CI)"),
  effect = 0.91,  # 确保文本在x轴的左侧边缘附近
  model = "Effect size"  # 将文本放置在分面图的顶部
)
label_data$status <- factor(
  label_data$status, 
  levels = c(
             "5-years MDD Risk","10-years MDD Risk" , "Overall MDD Risk"))
# 绘制森林图
p<-ggplot(Result, aes(x = effect, y = model, color = status)) +
  geom_errorbarh(aes(xmin = lower..95, xmax = upper..95), height = 0.2, color = "black", linewidth = 0.7) +  # 黑色误差棒
  geom_point(size = 3) +  # 绘制点，颜色根据 status 区分
  theme_minimal() +  # 使用简洁主题
  labs(x = "Effect",y = "Model") +  # 设置标签
  theme(
    panel.grid.major.x = element_blank(),  # 去除 X 轴主网格线
    panel.grid.minor.x = element_blank(),  # 去除 X 轴次网格线
    panel.grid.minor.y = element_blank(),  # 去除 Y 轴次网格线
    legend.position = "right",  # 图例在右侧
    legend.text = element_text(size = 12),  # 增大图例字体大小
    legend.title = element_text(size = 14, face = "bold"),  # 图例标题加粗并增大字体
    axis.text.x = element_text(size = 14, color = "black"),  # X 轴刻度字体设为黑色
    axis.text.y = element_text(size = 14, color = "black"),  # Y 轴刻度字体设为黑色
    axis.title = element_text(size = 14, face = "bold", color = "black"),  # 坐标轴标题字体加粗且为黑色
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),  # X 轴刻度线
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),  # Y 轴刻度线
    strip.text = element_text(size = 12, face = "bold"),  # 增大分面标签字体并加粗
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # 添加图的边框
  ) +
  scale_color_manual(values = c("#66448C","#e4e45f","#008F91")) +  # 分组颜色设置
  scale_x_continuous(limits = c(0.9, 1.7), breaks = seq(0.9, 1.7, by = 0.1)) +  # 设置 X 轴范围和刻度
  facet_grid(status ~ ., scales = "free_y", space = "free", labeller = as_labeller(custom_labels)) + # 自定义分面标签
  annotate("segment", x = 1, xend = 1, y = 0.5, yend = 4.5, color = "#e53b2c", linewidth = 1) +  # 部分红色垂直线
  geom_segment(data = Result, aes(x = 0.9, xend = 1.7, y = 5, yend = 5), 
               inherit.aes = FALSE, color = "white", linewidth = 0.5) +
  geom_text(data = label_data, aes(x = effect, y = model, label = label), 
            inherit.aes = FALSE, color = "black", size = 4,fontface = "bold",  hjust = 0)
#6*7.5

ggsave(
  filename = "2. Clinical/forest plot cohort.pdf",  # 文件名
  plot = p,  # 图形对象
  width = 7.67,  # 宽度 (英寸)
  height =6.5,  # 高度 (英寸)
  dpi = 300,  # 分辨率设置为 300 DPI
  device = cairo_pdf  # 确保 PDF 的高质量输出
)
#########################################################
# PD to MDD
data <- data.frame(
  Group = rep(c("Non-PD", "PD"), each = 2),
  MDD_Status = rep(c("Non-MDD", "MDD"), 2),
  Count = c(173936, 25450, 34326, 8105)
)



p<-ggplot(data, aes(x = Group, y = Count, fill = MDD_Status)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Non-MDD" = "#D8D9DA", "MDD" = "#4865a9")) +
  labs(x = "Group", y = "Count") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),  # 去除 X 轴主网格线
    panel.grid.minor.x = element_blank(),  # 去除 X 轴次网格线
    panel.grid.minor.y = element_blank(),  # 去除 Y 轴次网格线
    legend.text = element_text(size = 12),  # 增大图例字体大小
    axis.text.x = element_text(size = 14, color = "black"),  # X 轴刻度字体设为黑色
    axis.text.y = element_text(size = 14, color = "black"),  # Y 轴刻度字体设为黑色
    axis.title = element_text(size = 14, face = "bold", color = "black"),  # 坐标轴标题字体加粗且为黑色
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),  # X 轴刻度线
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),  # Y 轴刻度线
    strip.text = element_text(size = 12, face = "bold"),  # 增大分面标签字体并加粗
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # 添加图的边框
  )
ggsave(
  filename = "15. Cohort/Cross_PD.pdf",  # 文件名
  plot = p,  # 图形对象
  width = 3.95,  # 宽度 (英寸)
  height = 2.5,  # 高度 (英寸)
  dpi = 300,  # 分辨率设置为 300 DPI
  device = cairo_pdf  # 确保 PDF 的高质量输出
)

# MDD to PD
data
data$Group <- factor(data$Group, levels = c("PD","Non-PD"))
data$MDD_Status <- factor(data$MDD_Status, levels = c("Non-MDD","MDD"))
p<-ggplot(data, aes(x = MDD_Status, y = Count, fill =Group )) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Non-PD" = "#D8D9DA", "PD" = "#ef8a43")) +
  labs(x = "Group", y = "Count") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),  # 去除 X 轴主网格线
    panel.grid.minor.x = element_blank(),  # 去除 X 轴次网格线
    panel.grid.minor.y = element_blank(),  # 去除 Y 轴次网格线
    legend.text = element_text(size = 12),  # 增大图例字体大小
    axis.text.x = element_text(size = 14, color = "black"),  # X 轴刻度字体设为黑色
    axis.text.y = element_text(size = 14, color = "black"),  # Y 轴刻度字体设为黑色
    axis.title = element_text(size = 14, face = "bold", color = "black"),  # 坐标轴标题字体加粗且为黑色
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),  # X 轴刻度线
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),  # Y 轴刻度线
    strip.text = element_text(size = 12, face = "bold"),  # 增大分面标签字体并加粗
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # 添加图的边框
  )
ggsave(
  filename = "15. Cohort/Cross_MDD.pdf",  # 文件名
  plot = p,  # 图形对象
  width = 3.8,  # 宽度 (英寸)
  height = 2.5,  # 高度 (英寸)
  dpi = 300,  # 分辨率设置为 300 DPI
  device = cairo_pdf  # 确保 PDF 的高质量输出
)





# PD to MDD (PSM)
data <- data.frame(
  Group = rep(c("Non-PD", "PD"), each = 2),
  MDD_Status = rep(c("Non-MDD", "MDD"), 2),
  Count = c(36703, 5728, 34326, 8105)
)



p<-ggplot(data, aes(x = Group, y = Count, fill = MDD_Status)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Non-MDD" = "#D8D9DA", "MDD" = "#4865a9")) +
  labs(x = "Group", y = "Count") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),  # 去除 X 轴主网格线
    panel.grid.minor.x = element_blank(),  # 去除 X 轴次网格线
    panel.grid.minor.y = element_blank(),  # 去除 Y 轴次网格线
    legend.text = element_text(size = 12),  # 增大图例字体大小
    axis.text.x = element_text(size = 14, color = "black"),  # X 轴刻度字体设为黑色
    axis.text.y = element_text(size = 14, color = "black"),  # Y 轴刻度字体设为黑色
    axis.title = element_text(size = 14, face = "bold", color = "black"),  # 坐标轴标题字体加粗且为黑色
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),  # X 轴刻度线
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),  # Y 轴刻度线
    strip.text = element_text(size = 12, face = "bold"),  # 增大分面标签字体并加粗
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # 添加图的边框
  )
ggsave(
  filename = "15. Cohort/PSM_PD.pdf",  # 文件名
  plot = p,  # 图形对象
  width = 3.85,  # 宽度 (英寸)
  height = 2.5,  # 高度 (英寸)
  dpi = 300,  # 分辨率设置为 300 DPI
  device = cairo_pdf  # 确保 PDF 的高质量输出
)

# MDD to PD
data <- data.frame(
  Group = rep(c("Non-PD", "PD"), each = 2),
  MDD_Status = rep(c("Non-MDD", "MDD"), 2),
  Count = c(27754, 25450, 5801, 8105)
)
data$Group <- factor(data$Group, levels = c("PD","Non-PD"))
data$MDD_Status <- factor(data$MDD_Status, levels = c("Non-MDD","MDD"))
p<-ggplot(data, aes(x = MDD_Status, y = Count, fill =Group )) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Non-PD" = "#D8D9DA", "PD" = "#ef8a43")) +
  labs(x = "Group", y = "Count") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),  # 去除 X 轴主网格线
    panel.grid.minor.x = element_blank(),  # 去除 X 轴次网格线
    panel.grid.minor.y = element_blank(),  # 去除 Y 轴次网格线
    legend.text = element_text(size = 12),  # 增大图例字体大小
    axis.text.x = element_text(size = 14, color = "black"),  # X 轴刻度字体设为黑色
    axis.text.y = element_text(size = 14, color = "black"),  # Y 轴刻度字体设为黑色
    axis.title = element_text(size = 14, face = "bold", color = "black"),  # 坐标轴标题字体加粗且为黑色
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),  # X 轴刻度线
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),  # Y 轴刻度线
    strip.text = element_text(size = 12, face = "bold"),  # 增大分面标签字体并加粗
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # 添加图的边框
  )
ggsave(
  filename = "15. Cohort/PSM_MDD.pdf",  # 文件名
  plot = p,  # 图形对象
  width = 3.73,  # 宽度 (英寸)
  height = 2.5,  # 高度 (英寸)
  dpi = 300,  # 分辨率设置为 300 DPI
  device = cairo_pdf  # 确保 PDF 的高质量输出
)



# +================================================+ ####
# +====Section x. tableone===========================+ ####
# +================================================+ #### 

library(tableone)
load(file="2. Clinical/Interpolation_data_cross.Rdata")
Interpolation_data_cross<-Interpolation_data_cross[!is.na(Interpolation_data_cross$Periodontal_disease) 
                                                   & !is.na(Interpolation_data_cross$MDD), ]
# 定义变量列表
variables <- c("MDD", "Age", "Sex", "Ethnicity", "Education", "TDI", "Smoke", "Alcohol", "PA", 
               "BMI", "Systolic", "Diastolic", "Diabetes", "Medication", "TDI_quantile", "BMI_status", "Hypertension")

# 指定哪些变量是分类变量
factor_vars <- c("Sex", "Ethnicity", "Education", "Smoke", "Alcohol", "PA", "Diabetes", 
                 "Medication", "TDI_quantile", "BMI_status", "Hypertension", "MDD")

# 创建描述性统计表格
table <- CreateTableOne(vars = variables, 
                        strata = "Periodontal_disease",  # 根据Periodontal_disease分组
                        data = Interpolation_data_cross, 
                        factorVars = factor_vars,
                        test = TRUE)

# 打印表格，包含均值±标准误和显著性检验结果
table_df <- as.data.frame(print(table, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE))


write.csv(table_df, "2. Clinical/Interpolation_data_cross_summary.csv", row.names = TRUE)



# +================================================+ ####
# +====Section x. tableone cross===========================+ ####
# +================================================+ #### 




library(tableone)
load(file="2. Clinical/Interpolation_data_cross.Rdata")
Interpolation_data_cross<-Interpolation_data_cross[!is.na(Interpolation_data_cross$Periodontal_disease) 
                                                   & !is.na(Interpolation_data_cross$MDD), ]

mean_age <- mean(Interpolation_data_cross$Age, na.rm = TRUE)  # 计算age列的均值
n <- sum(!is.na(Interpolation_data_cross$Age))  # 计算非NA值的个数
se_age <- sd(Interpolation_data_cross$Age, na.rm = TRUE) / sqrt(n)  # 计算age列的标准误



# 定义变量列表
variables <- c("MDD", "Age", "Sex", "Ethnicity", "Education", "TDI", "Smoke", "Alcohol", "PA", 
               "BMI", "Systolic", "Diastolic", "Diabetes", "Medication", "TDI_quantile", "BMI_status", "Hypertension")

# 指定哪些变量是分类变量
factor_vars <- c("Sex", "Ethnicity", "Education", "Smoke", "Alcohol", "PA", "Diabetes", 
                 "Medication", "TDI_quantile", "BMI_status", "Hypertension", "MDD")

# 创建描述性统计表格
table <- CreateTableOne(vars = variables, 
                        strata = "Periodontal_disease",  # 根据Periodontal_disease分组
                        data = Interpolation_data_cross, 
                        factorVars = factor_vars,
                        test = TRUE)

# 打印表格，包含均值±标准误和显著性检验结果
table_df <- as.data.frame(print(table, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE))


write.csv(table_df, "2. Clinical/Interpolation_data_cross_summary.csv", row.names = TRUE)



# +================================================+ ####
# +====Section x. tableone orginal===========================+ ####
# +================================================+ #### 
# >>>>> section 4 Merge data  ####
load(file="1. Original/MDD_ICD.Rdata")
load(file="1. Original/MDD_CS.Rdata")
load(file="1. Original/Covariate_data_orginal.Rdata")
load(file="1. Original/Oral_data_final.Rdata")
load(file = "1. Original/Protein_final.Rdata")

MDD_gene_protein<-Protein_final[,1:3]
PD_gene_protein<-MDD_gene_protein
Protein_data<-as.data.frame(t(Protein_final[,-c(1:3)]))
colnames(Protein_data)<-Protein_final$Gene_symbol
rownames(Protein_data)
Protein_data$eid<-rownames(Protein_data)
Data_final_cohort<-merge(Oral_data_final,MDD_ICD,by = "eid",all= F)
Data_final_cross<-merge(Oral_data_final,MDD_CS,by = "eid",all= F)
Interpolation_data_cross<-merge(Data_final_cross,Covariate_data_orginal,by = "eid",all= F)
Interpolation_data_cross<-Interpolation_data_cross[!is.na(Interpolation_data_cross$Periodontal_disease) 
                                                   & !is.na(Interpolation_data_cross$MDD), ]

# 定义变量列表
variables <- c("MDD", "Age", "Sex", "Ethnicity", "Education", "TDI", "Smoke", "Alcohol", "PA", 
               "BMI", "Systolic", "Diastolic", "Diabetes", "Medication", "TDI_quantile", "BMI_status", "Hypertension")

# 指定哪些变量是分类变量
factor_vars <- c("Sex", "Ethnicity", "Education", "Smoke", "Alcohol", "PA", "Diabetes", 
                 "Medication", "TDI_quantile", "BMI_status", "Hypertension", "MDD")

# 创建描述性统计表格
table <- CreateTableOne(vars = variables, 
                        strata = "Periodontal_disease",  # 根据Periodontal_disease分组
                        data = Interpolation_data_cross, 
                        factorVars = factor_vars,
                        test = TRUE)

# 打印表格，包含均值±标准误和显著性检验结果
table_df <- as.data.frame(print(table, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE))


write.csv(table_df, "2. Clinical/Interpolation_data_crossorginal_summary.csv", row.names = TRUE)


# +================================================+ ####
# +====Section x. tableone psm PD===========================+ ####
# +================================================+ #### 
# >>>>> section 4 Merge data  ####

library(tableone)
load(file="2. Clinical/Interpolation_MDD_Cross.Rdata")
Interpolation_MDD_Cross<-Interpolation_MDD_Cross[!is.na(Interpolation_MDD_Cross$Periodontal_disease) 
                                                   & !is.na(Interpolation_MDD_Cross$MDD), ]

# 定义变量列表
variables <- c("MDD", "Age", "Sex", "Ethnicity", "Education", "TDI", "Smoke", "Alcohol", "PA", 
               "BMI", "Systolic", "Diastolic", "Diabetes", "Medication", "TDI_quantile", "BMI_status", "Hypertension")

# 指定哪些变量是分类变量
factor_vars <- c("Sex", "Ethnicity", "Education", "Smoke", "Alcohol", "PA", "Diabetes", 
                 "Medication", "TDI_quantile", "BMI_status", "Hypertension", "MDD")

# 创建描述性统计表格
table <- CreateTableOne(vars = variables, 
                        strata = "Periodontal_disease",  # 根据Periodontal_disease分组
                        data = Interpolation_MDD_Cross, 
                        factorVars = factor_vars,
                        test = TRUE)

# 打印表格，包含均值±标准误和显著性检验结果
table_df <- as.data.frame(print(table, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE))


write.csv(table_df, "2. Clinical/Interpolation_MDD_cross_summary.csv", row.names = TRUE)




# +================================================+ ####
# +====Section x. tableone psm MDD===========================+ ####
# +================================================+ #### 
# >>>>> section 4 Merge data  ####

library(tableone)
load(file="2. Clinical/Interpolation_PD_Cross.Rdata")
Interpolation_MDD_Cross<-Interpolation_PD_Cross[!is.na(Interpolation_PD_Cross$Periodontal_disease) 
                                                 & !is.na(Interpolation_PD_Cross$MDD), ]

# 定义变量列表
variables <- c("Periodontal_disease", "Age", "Sex", "Ethnicity", "Education", "TDI", "Smoke", "Alcohol", "PA", 
               "BMI", "Systolic", "Diastolic", "Diabetes", "Medication", "TDI_quantile", "BMI_status", "Hypertension")

# 指定哪些变量是分类变量
factor_vars <- c("Sex", "Ethnicity", "Education", "Smoke", "Alcohol", "PA", "Diabetes", 
                 "Medication", "TDI_quantile", "BMI_status", "Hypertension", "Periodontal_disease")

# 创建描述性统计表格
table <- CreateTableOne(vars = variables, 
                        strata = "MDD",  # 根据Periodontal_disease分组
                        data =Interpolation_PD_Cross, 
                        factorVars = factor_vars,
                        test = TRUE)

# 打印表格，包含均值±标准误和显著性检验结果
table_df <- as.data.frame(print(table, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE))


write.csv(table_df, "2. Clinical/Interpolation_PD_cross_summary.csv", row.names = TRUE)





# +================================================+ ####
# +====Section x. tableone cohort===========================+ ####
# +================================================+ #### 

library(tableone)
load(file="2. Clinical/Interpolation_data_cohort.Rdata")
Interpolation_data_cohort<-Interpolation_data_cohort[!is.na(Interpolation_data_cohort$Periodontal_disease) 
                                                   & !is.na(Interpolation_data_cohort$MDD)
                                                   & !is.na(Interpolation_data_cohort$days), ]
# 定义变量列表
variables <- c("MDD", "Age", "Sex", "Ethnicity", "Education", "TDI", "Smoke", "Alcohol", "PA", 
               "BMI", "Systolic", "Diastolic", "Diabetes", "Medication", "TDI_quantile", "BMI_status", "Hypertension")

# 指定哪些变量是分类变量
factor_vars <- c("Sex", "Ethnicity", "Education", "Smoke", "Alcohol", "PA", "Diabetes", 
                 "Medication", "TDI_quantile", "BMI_status", "Hypertension", "MDD")

# 创建描述性统计表格
table <- CreateTableOne(vars = variables, 
                        strata = "Periodontal_disease",  # 根据Periodontal_disease分组
                        data = Interpolation_data_cohort, 
                        factorVars = factor_vars,
                        test = TRUE)

# 打印表格，包含均值±标准误和显著性检验结果
table_df <- as.data.frame(print(table, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE))


write.csv(table_df, "2. Clinical/Interpolation_data_cohort_summary.csv", row.names = TRUE)



# +================================================+ ####
# +====Section x. MD/PD/total===========================+ ####
# +================================================+ #### 


load(file="2. Clinical/cross_original_Result.Rdata")

write.csv(Result, "2. Clinical/cross_original_Result.csv", row.names = TRUE)


load(file="2. Clinical/cross_psm_Result.Rdata")

write.csv(Result, "2. Clinical/cross_psm_Result.csv", row.names = TRUE)

load(file="2. Clinical/cohort_Result.Rdata")

write.csv(Result, "2. Clinical/cohort_Result.csv", row.names = TRUE)


load(file="2. Clinical/Interpolation_data_cross.Rdata")
Interpolation_data_cross<-Interpolation_data_cross[!is.na(Interpolation_data_cross$Periodontal_disease) 
                                                   & !is.na(Interpolation_data_cross$MDD), ]
table(Interpolation_data_cross$Periodontal_disease,Interpolation_data_cross$MDD)
table(Interpolation_data_cross$Periodontal_disease)
load(file="2. Clinical/Interpolation_MDD_Cross.Rdata")
table(Interpolation_MDD_Cross$Periodontal_disease,Interpolation_MDD_Cross$MDD)
table(Interpolation_MDD_Cross$Periodontal_disease)

load(file="2. Clinical/Interpolation_data_cross.Rdata")
table(Interpolation_data_cross$Periodontal_disease,Interpolation_data_cross$MDD)

load(file="2. Clinical/Interpolation_PD_Cross.Rdata")
table(Interpolation_PD_Cross$Periodontal_disease,Interpolation_PD_Cross$MDD)


load(file="2. Clinical/Interpolation_data_cohort.Rdata")
table(Interpolation_data_cohort$MDD,Interpolation_data_cohort$Periodontal_disease)


load(file="2. Clinical/Interpolation_data_cohort.Rdata")

Interpolation_data_cohort$days_5<- pmin(Interpolation_data_cohort$days, 5*365)  # 如果时间超过5年，设为5年
Interpolation_data_cohort$MDD_5 <- ifelse(Interpolation_data_cohort$days > 5*365, 0, Interpolation_data_cohort$MDD)  # 超过5年视为无事件

table(Interpolation_data_cohort$MDD_5,Interpolation_data_cohort$Periodontal_disease)

load(file="2. Clinical/Interpolation_data_cohort.Rdata")

Interpolation_data_cohort$days_5<- pmin(Interpolation_data_cohort$days, 10*365)  # 如果时间超过5年，设为5年
Interpolation_data_cohort$MDD_5 <- ifelse(Interpolation_data_cohort$days > 10*365, 0, Interpolation_data_cohort$MDD)  # 超过5年视为无事件

table(Interpolation_data_cohort$MDD_5,Interpolation_data_cohort$Periodontal_disease)


mean(Interpolation_data_cohort$days)/365

Table <-c("Type","Effect size", "95% CI", "P", NA, NA, NA)
Table <- rbind(Table, c("Cross-sectional study", "OR", "", "", NA, NA,NA, NA))
Table <- rbind(Table,c("PD to MDD", "", "", "", NA, NA, NA))
Table <- rbind(Table,c("Univariate Cox model",Results[1,]))
Table

Table <- rbind(Table,c("Cohort study", "HR", "", "", NA, NA, NA))
Table
Table <- rbind(Table,c("Univariate Cox model",Results[1,]))
Table <- rbind(Table,c("Multivariate Cox model",Results[2,]))
Table <- rbind(Table,c("PSM model",Results[3,]))
Table <- rbind(Table, c("Cross-sectional study", "OR", "", "", NA, NA,NA, NA))
Table <- rbind(Table,c("Univariate logistics model",Results[4,]))
Table <- rbind(Table,c("Multivariate logistics model",Results[5,]))
Table <- rbind(Table,c("PSM model",Results[6,]))
Table

colnames(Table)<-c('V1','V2','V3','V4','V5','V6','V7')
Table<-as.data.frame(Table)
Table
Table$V5 <- as.numeric(Table$V5)
Table$V6 <- as.numeric(Table$V6)
Table$V7 <- as.numeric(Table$V7)

forest <- Table
forest 
paste(c(rep("T",2),rep("F",3),"T",rep("F",3)), sep=",",collapse = ",")
a1<-paste(forest[,1])
a2<-paste(forest[,2])
a3<-paste(forest[,3])
a4<-paste(forest[,4])
labeltext=cbind(a1,a2,a3,a4)
labeltext[9,2]<-"1.470"
 library(forestplot)
#pdf("Fig. 3a .pdf",  height=3,width=11, onefile = FALSE)
forestplot(labeltext = labeltext,
           graphwidth = unit(45,'mm'),
           mean = forest$V5,
           lower = forest$V6,
           upper = forest$V7,
           col = fpColors(line = "#CC79A7", box="#D55E00"),
           is.summary = c(T,T,F,F,F,T,F,F,F),
           zero = 1,
           boxsize = 0.25,
           lineheight = unit(6,'mm'),
           colgap = unit(8,'mm'),
           xaxt = "n",
           xlab = "",
           cex.axis = 4,
           lwd.zero = 2,
           lwd.ci = 2,
           xticks = c(0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7),
           clip = c(-1,7),
           lwd.xaxis = 3,
           lty.ci = "solid",
           graph.pos = 4)
