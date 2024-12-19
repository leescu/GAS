# ++++++++++++++++++++++++++++++++++++++++++++++++++++++ ####
# >>>>> Section 0. Packages and Functions used <<<<< ####
setwd("D:/PD_MDD/step 0. Data")
library(data.table)
library(survival)
library(tableone)
library(Matching)
library(mediation)
library(survey)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(ggsci)
library(survminer)
library(tidyr)
library(forestplot)
library(mgcv)
library(plyr)
library(pheatmap)
library(car)

# 加载数据
load(file="16. Cross_sectional/Interpolation_data_cross.Rdata")
Interpolation_data_cross<-na.omit(Interpolation_data_cross)
load(file = "1. Original/Protein_data.Rdata")
Protein_data_Cross<-merge(Interpolation_data_cross,Protein_data,by = "eid",all = F)
data<-Protein_data_Cross
table(data$Periodontal_disease)
data$Periodontal_disease <- as.character(data$Periodontal_disease)
data$Periodontal_disease <- as.numeric(data$Periodontal_disease)



# +================================================+ ####
# +====Section 0. MDD to PD========================+ ####
# +================================================+ #### 

X <- "Periodontal_disease"
colnames(data)[20]
M <- colnames(data)[20:2942]  # 中介变量（蛋白质）
Y <- "MDD"

original_results <- data.frame(Type=character(), Estimate=numeric(), SE=numeric(), value=numeric(), P=numeric(),
                               Protein=character(), stringsAsFactors=FALSE)


mediation_results <- data.frame(Mediator=character(), ACME=numeric(), ACME_CI_Lower=numeric(), ACME_CI_Upper=numeric(), ACME_p_value=numeric(),
                                ADE=numeric(), ADE_CI_Lower=numeric(), ADE_CI_Upper=numeric(), ADE_p_value=numeric(),
                                Prop_Mediated=numeric(), Prop_Mediated_CI_Lower=numeric(), Prop_Mediated_CI_Upper=numeric(), Prop_Mediated_p_value=numeric(),
                                Total_Effect=numeric(), Total_Effect_CI_Lower=numeric(), Total_Effect_CI_Upper=numeric(), Total_Effect_p_value=numeric(), stringsAsFactors=FALSE)

remove_constant_vars <- function(data) {
  data[, sapply(data, function(col) length(unique(col)) > 1)]
}

total_start_time <- Sys.time()
for (i in seq_along(M)) {
  loop_start_time <- Sys.time()
  m <- M[i]
  print(paste("当前分析第",i ,"个。") )

  subset_data <- na.omit(data[, c(X, m, Y)])
  subset_data$Protein<-subset_data[,m]

  med.fit <- glm(Protein ~ Periodontal_disease, data=subset_data, family=gaussian())
  A<-summary(med.fit)[["coefficients"]]
  med<-A["Periodontal_disease",]
  out.fit <- glm(MDD ~ Protein+Periodontal_disease,data=subset_data, family = binomial(link = "probit"))
  B<-summary(out.fit)[["coefficients"]]
  out<-B["Protein",]
  
  result1<-as.data.frame(rbind(c("a",med),c("b",out)))
  original_results <- rbind(original_results, data.frame(Type=result1[1], Estimate=result1[2], SE=result1[3], value=result1[4], P=result1[5],
                                                         Protein=m))
  
  if (all(!is.na(c(A, B)))) {
    med.out <- mediate(med.fit, out.fit, treat=X, mediator="Protein", sims=100)

    result <- rbind(c("ACME", med.out$d.avg, med.out$d.avg.ci[1], med.out$d.avg.ci[2], med.out$d.avg.p),
                    c("ADE", med.out$z.avg, med.out$z.avg.ci[1], med.out$z.avg.ci[2], med.out$z.avg.p),
                    c("Prop. Mediated", med.out$n.avg, med.out$n.avg.ci[1], med.out$n.avg.ci[2], med.out$n.avg.p),
                    c("Total Effect", med.out$tau.coef, med.out$tau.ci[1], med.out$tau.ci[2], med.out$tau.p))

    mediation_results <- rbind(mediation_results, data.frame(Mediator=m,
                                                             ACME=result[1,2], ACME_CI_Lower=result[1,3], ACME_CI_Upper=result[1,4], ACME_p_value=result[1,5],
                                                             ADE=result[2,2], ADE_CI_Lower=result[2,3], ADE_CI_Upper=result[2,4], ADE_p_value=result[2,5],
                                                             Prop_Mediated=result[3,2], Prop_Mediated_CI_Lower=result[3,3], Prop_Mediated_CI_Upper=result[3,4], Prop_Mediated_p_value=result[3,5],
                                                             Total_Effect=result[4,2], Total_Effect_CI_Lower=result[4,3], Total_Effect_CI_Upper=result[4,4], Total_Effect_p_value=result[4,5]))
  } else {
    mediation_results <- rbind(mediation_results, data.frame(Mediator=m,
                                                             ACME=NA, ACME_CI_Lower=NA, ACME_CI_Upper=NA, ACME_p_value=NA,
                                                             ADE=NA, ADE_CI_Lower=NA, ADE_CI_Upper=NA, ADE_p_value=NA,
                                                             Prop_Mediated=NA, Prop_Mediated_CI_Lower=NA, Prop_Mediated_CI_Upper=NA, Prop_Mediated_p_value=NA,
                                                             Total_Effect=NA, Total_Effect_CI_Lower=NA, Total_Effect_CI_Upper=NA, Total_Effect_p_value=NA))
    print(paste("跳过第", i, "个中介变量，由于模型系数中存在 NA 值。"))
  }

  loop_end_time <- Sys.time()

  loop_time <- loop_end_time - loop_start_time
  print(paste("剩余进度",round((1-i/length(colnames(data)[20:2942]))*100,2),"%。  ","本次循环用时",round(loop_time, 2),"秒"))
}
write.csv(mediation_results,file="16. Cross_sectional/mediation_results_Cross_PD.csv")
save(original_results,file="16. Cross_sectional/original_results_Cross_PD.Rdata")

# +================================================+ ####
# +====Section 2. MDD to PD========================+ ####
# +================================================+ #### 

Y <- "Periodontal_disease"
colnames(data)[20]
M <- colnames(data)[20:2942]  # 中介变量（蛋白质）
X <- "MDD"

original_results <- data.frame(Type=character(), Estimate=numeric(), SE=numeric(), value=numeric(), P=numeric(),
                               Protein=character(), stringsAsFactors=FALSE)


mediation_results <- data.frame(Mediator=character(), ACME=numeric(), ACME_CI_Lower=numeric(), ACME_CI_Upper=numeric(), ACME_p_value=numeric(),
                                ADE=numeric(), ADE_CI_Lower=numeric(), ADE_CI_Upper=numeric(), ADE_p_value=numeric(),
                                Prop_Mediated=numeric(), Prop_Mediated_CI_Lower=numeric(), Prop_Mediated_CI_Upper=numeric(), Prop_Mediated_p_value=numeric(),
                                Total_Effect=numeric(), Total_Effect_CI_Lower=numeric(), Total_Effect_CI_Upper=numeric(), Total_Effect_p_value=numeric(), stringsAsFactors=FALSE)

remove_constant_vars <- function(data) {
  data[, sapply(data, function(col) length(unique(col)) > 1)]
}

total_start_time <- Sys.time()
for (i in seq_along(M)) {
  loop_start_time <- Sys.time()
  m <- M[i]
  print(paste("当前分析第",i ,"个。") )
  
  
  # 去除NA值
  subset_data <- na.omit(data[, c(X, m, Y)])
  subset_data$Protein<-subset_data[,m]
  
  # 模型1：自变量到中介变量
  med.fit <- glm(Protein ~MDD, data=subset_data, family=gaussian())
  A<-summary(med.fit)[["coefficients"]]
  med<-A["MDD",]
  out.fit <- glm(Periodontal_disease ~ Protein+MDD,data=subset_data, family = binomial(link = "probit"))
  B<-summary(out.fit)[["coefficients"]]
  out<-B["Protein",]
  
  result1<-as.data.frame(rbind(c("a",med),c("b",out)))
  original_results <- rbind(original_results, data.frame(Type=result1[1], Estimate=result1[2], SE=result1[3], value=result1[4], P=result1[5],
                                                         Protein=m))
  
  if (all(!is.na(c(A, B)))) {
    med.out <- mediate(med.fit, out.fit, treat=X, mediator="Protein", sims=100)
    
    # 提取中介分析结果
    result <- rbind(c("ACME", med.out$d.avg, med.out$d.avg.ci[1], med.out$d.avg.ci[2], med.out$d.avg.p),
                    c("ADE", med.out$z.avg, med.out$z.avg.ci[1], med.out$z.avg.ci[2], med.out$z.avg.p),
                    c("Prop. Mediated", med.out$n.avg, med.out$n.avg.ci[1], med.out$n.avg.ci[2], med.out$n.avg.p),
                    c("Total Effect", med.out$tau.coef, med.out$tau.ci[1], med.out$tau.ci[2], med.out$tau.p))
    
    # 将中介分析结果添加到数据框中
    mediation_results <- rbind(mediation_results, data.frame(Mediator=m,
                                                             ACME=result[1,2], ACME_CI_Lower=result[1,3], ACME_CI_Upper=result[1,4], ACME_p_value=result[1,5],
                                                             ADE=result[2,2], ADE_CI_Lower=result[2,3], ADE_CI_Upper=result[2,4], ADE_p_value=result[2,5],
                                                             Prop_Mediated=result[3,2], Prop_Mediated_CI_Lower=result[3,3], Prop_Mediated_CI_Upper=result[3,4], Prop_Mediated_p_value=result[3,5],
                                                             Total_Effect=result[4,2], Total_Effect_CI_Lower=result[4,3], Total_Effect_CI_Upper=result[4,4], Total_Effect_p_value=result[4,5]))
  } else {
    mediation_results <- rbind(mediation_results, data.frame(Mediator=m,
                                                             ACME=NA, ACME_CI_Lower=NA, ACME_CI_Upper=NA, ACME_p_value=NA,
                                                             ADE=NA, ADE_CI_Lower=NA, ADE_CI_Upper=NA, ADE_p_value=NA,
                                                             Prop_Mediated=NA, Prop_Mediated_CI_Lower=NA, Prop_Mediated_CI_Upper=NA, Prop_Mediated_p_value=NA,
                                                             Total_Effect=NA, Total_Effect_CI_Lower=NA, Total_Effect_CI_Upper=NA, Total_Effect_p_value=NA))
    print(paste("跳过第", i, "个中介变量，由于模型系数中存在 NA 值。"))
  }
  
  # 记录循环结束时间
  loop_end_time <- Sys.time()
  
  # 计算并打印循环用时
  loop_time <- loop_end_time - loop_start_time
  print(paste("剩余进度",round((1-i/length(colnames(data)[20:2942]))*100,2),"%。  ","本次循环用时",round(loop_time, 2),"秒"))
}
write.csv(mediation_results,file="16. Cross_sectional/mediation_results_Cross_MDD.csv")
save(original_results,file="16. Cross_sectional/original_results_Cross_MDD.Rdata")

