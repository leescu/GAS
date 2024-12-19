
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

library(lavaan)



# 加载数据
load(file="2. Clinical/Interpolation_data_cross.Rdata")
Interpolation_data_cross<-na.omit(Interpolation_data_cross)
load(file = "1. Original/Protein_data.Rdata")
Protein_data_eid<-Protein_data$eid

Protein_data$eid<-NULL
missing_percentage <- colMeans(is.na(Protein_data))  # 计算每列缺失值的比例

# 只保留缺失值比例小于 80% 的列
filtered <-as.data.frame( Protein_data[, missing_percentage < 0.2])
colnames(filtered)<-"Miss"
filter<-subset(filtered,Miss==T)
rownames(filter)
Protein_data <-as.data.frame( subset(Protein_data, select = rownames(filter)))
rownames(Protein_data)<-Protein_data_eid

missing_percentage <- rowMeans(is.na(Protein_data))
Protein_data<-as.data.frame( Protein_data[missing_percentage < 0.2, ])

colnames(Protein_data)<-toupper(colnames(Protein_data ))
Protein_ID<-colnames(Protein_data)
Protein_data_eid<-rownames(Protein_data)
Protein_matrix<-t(Protein_data)
# 确保 Protein_data_filtered 是矩阵
library(impute)
Protein_data_filled <- impute.knn(Protein_matrix)$data
# 检查插补后的数据是否还有缺失值
sum(is.na(Protein_data_filled))  
# 提取插补后的数据集（第一个插补结果）
Protein_data_filled<-as.data.frame(t(Protein_data_filled))
save(Protein_data_filled,file = "1. Original/Protein_data_filled.Rdata")
load(file = "1. Original/Protein_data_filled.Rdata")
Protein_data_filled$eid<-rownames(Protein_data_filled)
Total_data<-merge(Interpolation_data_cross,Protein_data_filled,by = "eid",all = F)
Total_Protein<-Total_data[,20:2930]
rownames(Total_Protein)<-Total_data$eid

library(psych)
# # 进行 KMO 测试和巴特利特球形检验
# kmo_result <- KMO(Protein_data_filled)
# print(kmo_result)
# 
# bartlett_test <- cortest.bartlett(protein_data)
# print(bartlett_test)
# scree(Protein_data_filled)  
cor_matrix <- cor(Total_Protein, use = "pairwise.complete.obs")
# 计算特征值
eigenvalues <- eigen(cor_matrix)$values
aparallel <- parallel(var = 2911, subject = 20979, rep = 100, cent = 0.95)$eigen$qevpea




eigen_df <- data.frame(
  Factor = 1:length(eigenvalues),
  Eigenvalue = eigenvalues
)

# 绘制碎石图
ggplot(eigen_df, aes(x = Factor, y = Eigenvalue)) +
  geom_point(size = 3, color = "blue") +
  geom_line() +
  labs(title = "Scree Plot", x = "Number of Factors", y = "Eigenvalue") +
  theme_minimal() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red")  # 添加水平线


# 计算相邻特征根的差值（即斜率）
slopes <- diff(eigenvalues)

# 计算斜率的变化幅度（第二阶导数）
second_slopes <- diff(slopes)

# 找到最大斜率变化的位置（折点）
elbow_point <- which.min(second_slopes) + 1  # +1 是因为 diff 函数减少了一个位置

# 绘制碎石图并标记折点
screeplot <- plot(1:length(eigenvalues), eigenvalues, geom="line") + 
  geom_vline(xintercept = elbow_point, linetype="dashed", color="red") + 
  labs(x = "因子数量", y = "特征根") + 
  theme_minimal()

print(screeplot)

# 输出折点位置
cat("建议的因子数量为：", elbow_point, "\n")














# 使用 nScree 函数
#results <- nScree(x = eigenvalues, aparallel = aparallel)

efa_result_4 <- fa(Total_Protein, nfactors = 4, rotate = "varimax", scores = "regression")
save(cor_matrix,results,eigenvalues,aparallel,efa_result_4,file="1. Original/eaf_result_4.Rdata")
load(file="1. Original/eaf_result.Rdata")
# 提取因子得分
factor_scores <- efa_result_4 $scores

head(factor_scores)
loadings_matrix <- as.matrix(efa_result_4 $loadings)

scree(Total_Protein)
# 查看前几行的分组结果
head(loadings_df)
loadings_matrix <- as.matrix(efa_result_4 [["loadings"]])
# 将矩阵转换为数据框
data <- as.data.frame(loadings_matrix)

# 检查转换后的矩阵维度
print(dim(loadings_matrix))

data <- as.data.frame(loadings_matrix)
write.table(loadings_matrix, file ="1. Original/EFA_orginal.csv", sep ="," ,row.names =T, col.names =T)
# 设置蛋白质名为行名，因子名为列名
data <- read.csv(file = "1. Original/EFA_orginal.csv", header = TRUE)

# 找到每个蛋白对应的主因子（即载荷值最大的因子）
data$Factor <- apply(data, 1, function(x) which.max(abs(x)))# 查看前几行结果

table(data$Factor)


# 将蛋白和对应的因子分组保存为数据框
grouped_proteins <- data.frame(
  Protein = rownames(loadings_matrix),
  Group = protein_groups
)
table(grouped_proteins$Group)



grouped_proteins
Protein_data_Cross<-merge(Interpolation_data_cross,Protein_data,by = "eid",all = F)
data<-Protein_data_Cross
table(data$Periodontal_disease)
data$Periodontal_disease<- as.character(data$Periodontal_disease)
colnames(data)[2]<-"PD"
data$PD <- as.numeric(data$PD)


# 选择 Protein 列（第24到2900列）
protein_cols <- colnames(data)[20:50]

# 创建潜变量模型，将这些蛋白列动态添加为模块
protein_module_formula <- paste("ProteinModule =~", paste(protein_cols, collapse = " + "))

# SEM模型定义，包括PD和MDD的路径关系
model <- paste0(
  protein_module_formula, "\n",   # 潜变量定义部分
  "PD ~ ProteinModule\n",          # ProteinModule 对 PD 的影响
  "MDD ~ ProteinModule\n",         # ProteinModule 对 MDD 的影响
  "PD ~~ MDD"                      # PD 和 MDD 之间的协方差
)

# 拟合模型
fit <- sem(model, data = data)



fit_summary <- summary(fit, standardized = TRUE, fit.measures = TRUE)
parameter_estimates <- parameterEstimates(fit, standardized = TRUE)

# 输出模型结果和拟合指标
summary(fit, fit.measures = TRUE)

# 提取路径系数、标准误和 P 值
fit_summary <- summary(fit, standardized = TRUE, fit.measures = TRUE)

# 从 summary 中提取感兴趣的部分
parameter_estimates <- parameterEstimates(fit, standardized = TRUE)

# 将结果转换为 dataframe
sem_results <- data.frame(
  Path = paste(parameter_estimates$lhs, 
               parameter_estimates$op, 
               parameter_estimates$rhs, sep = " "),
  Estimate = parameter_estimates$est,
  Std_Estimate = parameter_estimates$std.all,
  Std_Error = parameter_estimates$se,
  Z_value = parameter_estimates$z,
  P_value = parameter_estimates$pvalue
)

# 查看结果的前几行
fit_summary <- summary(fit, standardized = TRUE, fit.measures = TRUE)
parameter_estimates <- parameterEstimates(fit, standardized = TRUE)
protein_paths <- parameter_estimates[grepl("Protein", parameter_estimates$lhs) | grepl("Protein", parameter_estimates$rhs), ]
























# 创建空列表存储 AIC 和 BIC 值
aic_values <- c()
bic_values <- c()

# 遍历不同的因子数量进行 EFA
for (factors in 1:20) {  # 假设最多尝试 10 个因子
  efa_result <- fa(Total_Protein, nfactors = factors, rotate = "varimax")
  
  # 提取模型的对数似然值
  log_likelihood <- efa_result$STATISTIC
  
  # 计算 AIC 和 BIC
  n_params <- factors * (factors - 1) / 2 + factors * ncol(Total_Protein)  # 参数数量
  n_samples <- nrow(Total_Protein)  # 样本量
  
  aic <- -2 * log_likelihood + 2 * n_params
  bic <- -2 * log_likelihood + n_params * log(n_samples)
  
  # 保存 AIC 和 BIC 值
  aic_values[factors] <- aic
  bic_values[factors] <- bic
}

# 将 AIC 和 BIC 存储为数据框
criteria_df <- data.frame(Factors = 1:20, AIC = aic_values, BIC = bic_values)

# 找到 AIC 和 BIC 最小的因子数量
best_aic <- which.min(criteria_df$AIC)
best_bic <- which.min(criteria_df$BIC)

print(criteria_df)
cat("最佳因子数量（AIC）：", best_aic, "\n")
cat("最佳因子数量（BIC）：", best_bic, "\n")
