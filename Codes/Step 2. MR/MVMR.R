# +================================================+ ####
# +====Section 0. Packages and setwd===============+ ####
# +================================================+ #### 

setwd("D:/PD_MDD/Step 0. Data")
library(data.table)
library(dplyr)
library(TwoSampleMR)
library(devtools)
library(tidyverse)
library(dplyr)
library(ieugwasr)
library(MungeSumstats)
library(MendelianRandomization)
library(MVMR)
load( file ="2. MR/Periodontal_exposure.Rdata")
load(file ="2. MR/Periodontal_outcome.Rdata")
load( file ="2. MR/Depressive_exposure.Rdata")
load(file ="2. MR/Depressive_outcome.Rdata")
load( file ="2. MR/Drink_exposure.Rdata")
load( file ="2. MR/Cognitive_exposure.Rdata")
Periodontal_orginal_data<-fread("2. MR/Periodontal_orginal_data.csv")
Cognitive_orginal_data<-fread("2. MR/Cognitive_orginal_data.csv")
Drink_orginal_data<-fread("2. MR/Drink_orginal_data.csv")
Depressive_orginal_data<-fread("2. MR/Depressive_orginal_data.csv")
Depressive_outcome <- subset(Depressive_outcome,eaf.outcome>0.01)
Periodontal_exposure<- subset(Periodontal_exposure,eaf.exposure>0.01)
# +================================================+ ####
# +====Section 0. Packages and setwd===============+ ####
# +================================================+ #### 
# 定义协变量的列表，每个协变量包含exposure和original数据集
covariates <- list(Drink = list(exposure = Drink_exposure, original = Drink_orginal_data, name = "Drink"),
                   Cognitive = list(exposure = Cognitive_exposure, original = Cognitive_orginal_data, name = "Cognitive"))

# 初始化结果列表
results <- list()

# 循环处理1到5个协变量的组合
for (n in 1:2) {
  cov_comb <- combn(names(covariates), n, simplify = FALSE)
  
  for (comb in cov_comb) {
    # 排除同时包含 "Smoking" 和 "Smoke"，或 "Drinking" 和 "Drink" 的组合
    if (("Smoking" %in% comb && "Smoke" %in% comb) || ("Drinking" %in% comb && "Drink" %in% comb)) {
      next  # 跳过这个组合
    }
    
    print(paste("Processing combination with", n, "covariates:", paste(comb, collapse = "+")))
    
    # 1. 从每个Exposure数据集中获取所有显著的SNP
    significant_snps <- unique(unlist(lapply(comb, function(cov_name) covariates[[cov_name]]$exposure$SNP)))
    
    # 2. 获取所有SNP并去重，包括显著SNP和Periodontal_exposure中的SNP
    all_snps <- unique(c(Periodontal_exposure$SNP, significant_snps))
    
    # 3. 初始化一个新的dataframe，包含所有唯一的SNP
    rawdat_mvmr <- data.frame(SNP = all_snps)
    
    # 4. 从Periodontal_orginal_data提取beta和se，并合并到rawdat_mvmr
    periodontal_data_filtered <- Periodontal_orginal_data[Periodontal_orginal_data$SNP %in% all_snps, c("SNP", "beta", "se")]
    rawdat_mvmr <- merge(rawdat_mvmr, periodontal_data_filtered, by = "SNP", all.x = TRUE)
    colnames(rawdat_mvmr)[2:3] <- c("Periodontal_beta", "Periodontal_se")
    
    # 5. 从每个协变量的Original数据集提取beta和se，并合并到rawdat_mvmr
    for (i in seq_along(comb)) {
      cov_name <- comb[i]
      cov_data_filtered <- covariates[[cov_name]]$original[covariates[[cov_name]]$original$SNP %in% all_snps, c("SNP", "beta", "se")]
      rawdat_mvmr <- merge(rawdat_mvmr, cov_data_filtered, by = "SNP", all.x = TRUE)
      colnames(rawdat_mvmr)[(i*2 + 2):(i*2 + 3)] <- c(paste0(cov_name, "_beta"), paste0(cov_name, "_se"))
    }
    
    # 6. 从Depressive_outcome提取相同SNP的beta和se，并合并到rawdat_mvmr
    depressive_data_filtered <- Depressive_outcome[Depressive_outcome$SNP %in% all_snps, c("SNP", "beta.outcome", "se.outcome")]
    rawdat_mvmr <- merge(rawdat_mvmr, depressive_data_filtered, by = "SNP", all.x = TRUE)
    colnames(rawdat_mvmr)[(n*2 + 4):(n*2 + 5)] <- c("MDD_beta", "MDD_se")
    
    # 7. 删除含有NA的行
    rawdat_mvmr <- na.omit(rawdat_mvmr)
    
    # 8. 创建MRMVInputObject
    bx <- as.matrix(rawdat_mvmr[, paste0(c("Periodontal", comb), "_beta")])
    bxse <- as.matrix(rawdat_mvmr[, paste0(c("Periodontal", comb), "_se")])
    MDD_beta <- as.numeric(rawdat_mvmr$MDD_beta)
    MDD_se <- as.numeric(rawdat_mvmr$MDD_se)
    
    MRMVInputObject_1 <- mr_mvinput(bx = bx, bxse = bxse, by = MDD_beta, byse = MDD_se)
    
    # 9. 执行MVMR分析
    MRMVObject_ivw <- mr_mvivw(MRMVInputObject_1, model = "default", correl = FALSE, distribution = "normal", alpha = 0.05)
    MRMVObject_egger <- mr_mvegger(MRMVInputObject_1, orientate = 1, correl = FALSE, distribution = "normal", alpha = 0.05)
    MRMVObject_median <- mr_mvmedian(MRMVInputObject_1, distribution = "normal", alpha = 0.05, iterations = 100, seed = 0)
    MRMVObject_lasso <- mr_mvlasso(MRMVInputObject_1, orientate = 1, distribution = "normal", alpha = 0.05, lambda = numeric(0))
    
    # 10. 收集结果，并添加协变量信息
    result <- rbind(c(MRMVObject_ivw@Estimate[1],MRMVObject_ivw@StdError[1],MRMVObject_ivw@Pvalue[1],"IVW",paste(comb, collapse = "+")),
                    c(MRMVObject_egger@Estimate[1],MRMVObject_egger@StdError.Est[1],MRMVObject_egger@Pvalue.Est[1],'MR Egger',paste(comb, collapse = "+")),
                    c(MRMVObject_median@Estimate[1],MRMVObject_median@StdError[1],MRMVObject_median@Pvalue[1],'Weighted median',paste(comb, collapse = "+")),
                    c(MRMVObject_lasso@Estimate[1],MRMVObject_lasso@StdError[1],MRMVObject_lasso@Pvalue[1],"LASSO",paste(comb, collapse = "+")))  # 添加协变量组合信息
    
    # 保存结果
    results[[paste(comb, collapse = "+")]] <- result
  }
}

# 将结果转换为数据框并输出
results_df <- do.call(rbind, results)
colnames(results_df) <- c("Estimate", "StdError", "Pvalue", "Model","Covariates")

write.table(results_df, file ="2. MR/results_MVMR_eaf_PD.csv", sep ="," ,row.names =F, col.names =T)








# +================================================+ ####
# +====Section 0. Packages and setwd===============+ ####
# +================================================+ #### 

Periodontal_outcome <- subset(Periodontal_outcome,eaf.outcome>0.01)
Depressive_exposure<- subset(Depressive_exposure,eaf.exposure>0.01)
# +================================================+ ####
# +====Section 0. Packages and setwd===============+ ####
# +================================================+ #### 
# 定义协变量的列表，每个协变量包含exposure和original数据集
covariates <- list(Drink = list(exposure = Drink_exposure, original = Drink_orginal_data, name = "Drink"),
                   Cognitive = list(exposure = Cognitive_exposure, original = Cognitive_orginal_data, name = "Cognitive"))

# 初始化结果列表
results <- list()

# 循环处理1到5个协变量的组合
for (n in 1:2) {
  cov_comb <- combn(names(covariates), n, simplify = FALSE)
  
  for (comb in cov_comb) {
    # 排除同时包含 "Smoking" 和 "Smoke"，或 "Drinking" 和 "Drink" 的组合
    if (("Smoking" %in% comb && "Smoke" %in% comb) || ("Drinking" %in% comb && "Drink" %in% comb)) {
      next  # 跳过这个组合
    }
    
    print(paste("Processing combination with", n, "covariates:", paste(comb, collapse = "+")))
    
    # 1. 从每个Exposure数据集中获取所有显著的SNP
    significant_snps <- unique(unlist(lapply(comb, function(cov_name) covariates[[cov_name]]$exposure$SNP)))
    
    # 2. 获取所有SNP并去重，包括显著SNP和Periodontal_exposure中的SNP
    all_snps <- unique(c(Depressive_exposure$SNP, significant_snps))
    
    # 3. 初始化一个新的dataframe，包含所有唯一的SNP
    rawdat_mvmr <- data.frame(SNP = all_snps)
    
    # 4. 从Periodontal_orginal_data提取beta和se，并合并到rawdat_mvmr
    Depressive_data_filtered <- Depressive_orginal_data[Depressive_orginal_data$SNP %in% all_snps, c("SNP", "beta", "se")]
    rawdat_mvmr <- merge(rawdat_mvmr, Depressive_data_filtered, by = "SNP", all.x = TRUE)
    colnames(rawdat_mvmr)[2:3] <- c("Depressive_beta", "Depressive_se")
    
    # 5. 从每个协变量的Original数据集提取beta和se，并合并到rawdat_mvmr
    for (i in seq_along(comb)) {
      cov_name <- comb[i]
      cov_data_filtered <- covariates[[cov_name]]$original[covariates[[cov_name]]$original$SNP %in% all_snps, c("SNP", "beta", "se")]
      rawdat_mvmr <- merge(rawdat_mvmr, cov_data_filtered, by = "SNP", all.x = TRUE)
      colnames(rawdat_mvmr)[(i*2 + 2):(i*2 + 3)] <- c(paste0(cov_name, "_beta"), paste0(cov_name, "_se"))
    }
    
    # 6. 从Depressive_outcome提取相同SNP的beta和se，并合并到rawdat_mvmr
    Periodontal_data_filtered <- Periodontal_outcome[Periodontal_outcome$SNP %in% all_snps, c("SNP", "beta.outcome", "se.outcome")]
    rawdat_mvmr <- merge(rawdat_mvmr,Periodontal_data_filtered, by = "SNP", all.x = TRUE)
    colnames(rawdat_mvmr)[(n*2 + 4):(n*2 + 5)] <- c("MDD_beta", "MDD_se")
    
    # 7. 删除含有NA的行
    rawdat_mvmr <- na.omit(rawdat_mvmr)
    
    # 8. 创建MRMVInputObject
    bx <- as.matrix(rawdat_mvmr[, paste0(c("Depressive", comb), "_beta")])
    bxse <- as.matrix(rawdat_mvmr[, paste0(c("Depressive", comb), "_se")])
    MDD_beta <- as.numeric(rawdat_mvmr$MDD_beta)
    MDD_se <- as.numeric(rawdat_mvmr$MDD_se)
    
    MRMVInputObject_1 <- mr_mvinput(bx = bx, bxse = bxse, by = MDD_beta, byse = MDD_se)
    
    # 9. 执行MVMR分析
    MRMVObject_ivw <- mr_mvivw(MRMVInputObject_1, model = "default", correl = FALSE, distribution = "normal", alpha = 0.05)
    MRMVObject_egger <- mr_mvegger(MRMVInputObject_1, orientate = 1, correl = FALSE, distribution = "normal", alpha = 0.05)
    MRMVObject_median <- mr_mvmedian(MRMVInputObject_1, distribution = "normal", alpha = 0.05, iterations = 100, seed = 0)
    MRMVObject_lasso <- mr_mvlasso(MRMVInputObject_1, orientate = 1, distribution = "normal", alpha = 0.05, lambda = numeric(0))
    
    # 10. 收集结果，并添加协变量信息
    result <- rbind(c(MRMVObject_ivw@Estimate[1],MRMVObject_ivw@StdError[1],MRMVObject_ivw@Pvalue[1],"IVW",paste(comb, collapse = "+")),
                    c(MRMVObject_egger@Estimate[1],MRMVObject_egger@StdError.Est[1],MRMVObject_egger@Pvalue.Est[1],'MR Egger',paste(comb, collapse = "+")),
                    c(MRMVObject_median@Estimate[1],MRMVObject_median@StdError[1],MRMVObject_median@Pvalue[1],'Weighted median',paste(comb, collapse = "+")),
                    c(MRMVObject_lasso@Estimate[1],MRMVObject_lasso@StdError[1],MRMVObject_lasso@Pvalue[1],"LASSO",paste(comb, collapse = "+")))  # 添加协变量组合信息
    
    # 保存结果
    results[[paste(comb, collapse = "+")]] <- result
  }
}

# 将结果转换为数据框并输出
results_df <- do.call(rbind, results)
colnames(results_df) <- c("Estimate", "StdError", "Pvalue", "Model","Covariates")
results_df 
write.table(results_df, file ="2. MR/results_MVMR_eaf_MDD.csv", sep ="," ,row.names =F, col.names =T)

MVMR_PD<-fread( file ="2. MR/results_MVMR_eaf_PD.csv")
OR <- exp(MVMR_PD$Estimate)
CI_lower<- exp(MVMR_PD$Estimate - 1.96 *MVMR_PD$StdError)
CI_upper<- exp(MVMR_PD$Estimate + 1.96 *MVMR_PD$StdError)
Result<-as.data.frame(cbind(MVMR_PD$Model,MVMR_PD$Covariates,OR,CI_lower,CI_upper,MVMR_PD$Pvalue))
colnames(Result)<-c("Model","Covariates","OR","CI_lower","CI_upper","Pval")
Result[ , 3:6] <- lapply(Result[ , 3:6], as.numeric)
# 定义一个函数将数值格式化为4位有效数字，并补齐不足位数
format_numeric <- function(x) {
  sprintf("%.3f", x)
}

format_p_value <- function(x) {
  format(x, scientific = TRUE, digits = 4)
}

Result$OR <- sapply(Result$OR, format_numeric)
Result$CI_lower <- sapply(Result$CI_lower, format_numeric)
Result$CI_upper <- sapply(Result$CI_upper, format_numeric)
Result$Pval <- sapply(Result$Pval, format_p_value)
Result


# 新增列：显示效应值及95%置信区间范围
Result$CI <- paste0(
  Result$effect, " (", Result$lower..95, " - ", Result$upper..95, ")"
)
save(Result,file="2. MR/results_MVMR_PD.Rdata")


Result$OR <- as.numeric(Result$OR)
Result$CI_lower <- as.numeric(Result$CI_lower)
Result$CI_upper <- as.numeric(Result$CI_upper)
Result$Model <- factor(Result$Model, levels = c('LASSO', 'IVW', 'Weighted median', 'MR Egger','Effect size'))
#创建自定义标签函数
custom_labels <- c(
  "Drink" = "Drinking",
  "Cognitive" = "Cognitive",
  "Drink+Cognitive"="Drinking & Cognitive")
Result$Covariates <- factor(
  Result$Covariates, 
  levels = c("Drink","Cognitive","Drink+Cognitive")
)
# 在每个分面图添加的左上角标签数据
label_data <- data.frame(
  Covariates = c("Drink","Cognitive","Drink+Cognitive"),
  label = c("Odds Ratio (95%CI)","Odds Ratio (95%CI)","Odds Ratio (95%CI)"),
  OR = 0.91,  # 确保文本在x轴的左侧边缘附近
  Model = "Effect size"  # 将文本放置在分面图的顶部
)
label_data$Covariates <- factor(
  label_data$Covariates, 
  levels = c("Drink","Cognitive","Drink+Cognitive"))
label_data$Model <- factor(
  label_data$Model, 
  levels = c('LASSO', 'IVW', 'Weighted median', 'MR Egger','Effect size'))

# 绘制森林图
p<-ggplot(Result, aes(x = OR, y = Model, color = Covariates)) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2, color = "black", linewidth = 0.7) +  # 黑色误差棒
  geom_point(size = 3) +  # 绘制点，颜色根据 status 区分
  theme_minimal() +  # 使用简洁主题
  labs(x = "PD MVMR effect size",y = "Model") +  # 设置标签
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
  scale_x_continuous(limits = c(0.9, 2.0), breaks = seq(0.9, 2.0, by = 0.2)) +  # 设置 X 轴范围和刻度
  facet_grid(Covariates ~ ., scales = "free_y", space = "free", labeller = as_labeller(custom_labels)) + # 自定义分面标签
  annotate("segment", x = 1, xend = 1, y = 0.5, yend = 4.5, color = "#e53b2c", linewidth = 1) +  # 部分红色垂直线
  geom_segment(data = Result, aes(x = 0.9, xend = 1.7, y = 5, yend = 5), 
               inherit.aes = FALSE, color = "white", linewidth = 0.5) +
  geom_text(data = label_data, aes(x = OR , y = Model, label = label), 
            inherit.aes = FALSE, color = "black", size = 4,fontface = "bold",  hjust = 0)
#6*7.5
p
ggsave(
  filename = "2. MR/PD MVMR forest plot.pdf",  # 文件名
  plot = p,  # 图形对象
  width = 7.5,  # 宽度 (英寸)
  height = 6,  # 高度 (英寸)
  dpi = 300,  # 分辨率设置为 300 DPI
  device = cairo_pdf  # 确保 PDF 的高质量输出
)





MVMR_PD<-fread( file ="2. MR/results_MVMR_eaf_MDD.csv")
OR <- exp(MVMR_PD$Estimate)
CI_lower<- exp(MVMR_PD$Estimate - 1.96 *MVMR_PD$StdError)
CI_upper<- exp(MVMR_PD$Estimate + 1.96 *MVMR_PD$StdError)
Result<-as.data.frame(cbind(MVMR_PD$Model,MVMR_PD$Covariates,OR,CI_lower,CI_upper,MVMR_PD$Pvalue))
colnames(Result)<-c("Model","Covariates","OR","CI_lower","CI_upper","Pval")
Result[ , 3:6] <- lapply(Result[ , 3:6], as.numeric)
# 定义一个函数将数值格式化为4位有效数字，并补齐不足位数
format_numeric <- function(x) {
  sprintf("%.3f", x)
}

format_p_value <- function(x) {
  format(x, scientific = TRUE, digits = 4)
}

Result$OR <- sapply(Result$OR, format_numeric)
Result$CI_lower <- sapply(Result$CI_lower, format_numeric)
Result$CI_upper <- sapply(Result$CI_upper, format_numeric)
Result$Pval <- sapply(Result$Pval, format_p_value)
Result


# 新增列：显示效应值及95%置信区间范围
Result$CI <- paste0(
  Result$effect, " (", Result$lower..95, " - ", Result$upper..95, ")"
)
save(Result,file="2. MR/results_MVMR_MDD.Rdata")


Result$OR <- as.numeric(Result$OR)
Result$CI_lower <- as.numeric(Result$CI_lower)
Result$CI_upper <- as.numeric(Result$CI_upper)
Result$Model <- factor(Result$Model, levels = c('LASSO', 'IVW', 'Weighted median', 'MR Egger','Effect size'))
#创建自定义标签函数
custom_labels <- c(
  "Drink" = "Drinking",
  "Cognitive" = "Cognitive",
  "Drink+Cognitive"="Drinking & Cognitive")
Result$Covariates <- factor(
  Result$Covariates, 
  levels = c("Drink","Cognitive","Drink+Cognitive")
)
# 在每个分面图添加的左上角标签数据
label_data <- data.frame(
  Covariates = c("Drink","Cognitive","Drink+Cognitive"),
  label = c("Odds Ratio (95%CI)","Odds Ratio (95%CI)","Odds Ratio (95%CI)"),
  OR = 0.91,  # 确保文本在x轴的左侧边缘附近
  Model = "Effect size"  # 将文本放置在分面图的顶部
)
label_data$Covariates <- factor(
  label_data$Covariates, 
  levels = c("Drink","Cognitive","Drink+Cognitive"))
label_data$Model <- factor(
  label_data$Model, 
  levels = c('LASSO', 'IVW', 'Weighted median', 'MR Egger','Effect size'))

# 绘制森林图
p<-ggplot(Result, aes(x = OR, y = Model, color = Covariates)) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2, color = "black", linewidth = 0.7) +  # 黑色误差棒
  geom_point(size = 3) +  # 绘制点，颜色根据 status 区分
  theme_minimal() +  # 使用简洁主题
  labs(x = "MDD MVMR effect size",y = "Model") +  # 设置标签
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
  scale_x_continuous(limits = c(0.9, 2.0), breaks = seq(0.9, 2.0, by = 0.2)) +  # 设置 X 轴范围和刻度
  facet_grid(Covariates ~ ., scales = "free_y", space = "free", labeller = as_labeller(custom_labels)) + # 自定义分面标签
  annotate("segment", x = 1, xend = 1, y = 0.5, yend = 4.5, color = "#e53b2c", linewidth = 1) +  # 部分红色垂直线
  geom_segment(data = Result, aes(x = 0.9, xend = 1.7, y = 5, yend = 5), 
               inherit.aes = FALSE, color = "white", linewidth = 0.5) +
  geom_text(data = label_data, aes(x = OR , y = Model, label = label), 
            inherit.aes = FALSE, color = "black", size = 4,fontface = "bold",  hjust = 0)
#6*7.5
p
ggsave(
  filename = "2. MR/MDD MVMR forest plot.pdf",  # 文件名
  plot = p,  # 图形对象
  width = 7.5,  # 宽度 (英寸)
  height = 6,  # 高度 (英寸)
  dpi = 300,  # 分辨率设置为 300 DPI
  device = cairo_pdf  # 确保 PDF 的高质量输出
)