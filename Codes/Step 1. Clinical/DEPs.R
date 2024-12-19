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
library(ggrepel) 
# +================================================+ ####
# +====Section 1. PD to MDD========================+ ####
# +================================================+ #### 
load(file="16. Cross_sectional/Interpolation_data_cross.Rdata")
Interpolation_data_cross<-na.omit(Interpolation_data_cross)
load(file = "1. Original/Protein_data.Rdata")
Protein_data_Cross<-merge(Interpolation_data_cross,Protein_data,by = "eid",all = F)
data<-Protein_data_Cross
table(data$Periodontal_disease)
data$Periodontal_disease <- as.character(data$Periodontal_disease)
data$Periodontal_disease <- as.numeric(data$Periodontal_disease)

Uni_glm_model<- 
  function(x){
    FML<-as.formula(paste0("Periodontal_disease~",x,"+Age+Sex+Ethnicity+Education+TDI_quantile+
             Smoke+Alcohol+PA+BMI_status+Diabetes+Hypertension+Medication"))
    PD<-glm(FML,data=data,family = binomial)
    SUM<-summary(PD)[["coefficients"]]
    OR<-round(exp(SUM[2,1]),3)
    CI5<-round(exp(SUM[2,1]-1.96*SUM[2,2]),3)
    CI95<-round(exp(SUM[2,1]+1.96*SUM[2,2]),3)
    P<-SUM[2,4]
    Uni_glm_model <- data.frame('Characteristics'=x,
                                'OR' = OR,
                                'CIL' = CI5,
                                'CIH' = CI95,
                                'P' = P,
                                'status'="PD to MDD cross")           
    return(Uni_glm_model)
  }  

variable.names<- colnames(data)[20:2942]
Uni_glm<- lapply(variable.names, Uni_glm_model)
Uni_glm<- ldply(Uni_glm,data.frame)
rownames(Uni_glm)<-Uni_glm$Characteristics
Uni_glm <- Uni_glm[,-1]
PD_glm_protein<-Uni_glm
write.csv(Uni_glm,file="16. Cross_sectional/PD_glm_protein.csv")
save(PD_glm_protein,file="16. Cross_sectional/PD_glm_protein.Rdata")
# +================================================+ ####
# +====Section 1. MDD to PD========================+ ####
# +================================================+ #### 
load(file="16. Cross_sectional/Interpolation_data_cross.Rdata")
Interpolation_data_cross<-na.omit(Interpolation_data_cross)
load(file = "1. Original/Protein_data.Rdata")
Protein_data_Cross<-merge(Interpolation_data_cross,Protein_data,by = "eid",all = F)
data<-Protein_data_Cross
table(data$MDD)
data$MDD <- as.character(data$MDD)
data$MDD<- as.numeric(data$MDD)

Uni_glm_model<- 
  function(x){
    FML<-as.formula(paste0("MDD~",x,"+Age+Sex+Ethnicity+Education+TDI_quantile+
             Smoke+Alcohol+PA+BMI_status+Diabetes+Hypertension+Medication"))
    PD<-glm(FML,data=data,family = binomial)
    SUM<-summary(PD)[["coefficients"]]
    OR<-round(exp(SUM[2,1]),3)
    CI5<-round(exp(SUM[2,1]-1.96*SUM[2,2]),3)
    CI95<-round(exp(SUM[2,1]+1.96*SUM[2,2]),3)
    P<-SUM[2,4]
    Uni_glm_model <- data.frame('Characteristics'=x,
                                'OR' = OR,
                                'CIL' = CI5,
                                'CIH' = CI95,
                                'P' = P,
                                'status'="MDD to PD cross")           
    return(Uni_glm_model)
  }  

variable.names<- colnames(data)[20:2942]
Uni_glm<- lapply(variable.names, Uni_glm_model)
Uni_glm<- ldply(Uni_glm,data.frame)
rownames(Uni_glm)<-Uni_glm$Characteristics
Uni_glm <- Uni_glm[,-1]
MDD_glm_protein<-Uni_glm
write.csv(Uni_glm,file="16. Cross_sectional/MDD_glm_protein.csv")
save(MDD_glm_protein,file="16. Cross_sectional/MDD_glm_protein.Rdata")

load(file="16. Cross_sectional/MDD_glm_protein.Rdata")
load(file="16. Cross_sectional/PD_glm_protein.Rdata")
# +================================================+ ####
# +====Section 4. Volcano of DEGs==================+ ####
# +================================================+ ####

library(EnhancedVolcano)

colnames(PD_glm_protein)
PD_data<-PD_glm_protein
PD_data$Gene<-toupper(rownames(PD_data))
PD_data$Beta<- log(PD_data$OR)
PD_data$Padj<-p.adjust(PD_data$P,method ="fdr")
MDD_data<-MDD_glm_protein
MDD_data$Gene<-toupper(rownames(MDD_data))
MDD_data$Beta<- log(MDD_data$OR)
MDD_data$Padj<-p.adjust(MDD_data$P,method ="fdr")
Merge_data<-merge(PD_data,MDD_data,by="Gene",all=F, suffixes = c("_PD", "_MDD"))

Merge_sig_data<-subset(Merge_data,(Padj_PD<0.05&Padj_MDD<0.05)& #P<0.05
                        ( (OR_PD>1.5&OR_MDD>1.5)|
                         (OR_MDD<2/3&OR_PD<2/3) ))
Merge_sig_data<-subset(Merge_data,(Padj_PD<0.05&Padj_MDD<0.05))
vals<-Merge_sig_data$Gene

# 计算 -log10(Padj) 用于 Y 轴
PD_data$neg_log10_Padj <- -log10(PD_data$Padj)

# 根据 OR 和 Padj 创建颜色分组
PD_data$category <- with(PD_data, 
                         ifelse(Gene %in% vals, "Shared", 
                                ifelse(OR > 1 & Padj < 0.05, "Upregulated", 
                                       ifelse(OR < 1 & Padj < 0.05, "Downregulated", 
                                              "Normal"))))
PD_data$category<-factor(PD_data$category, levels = c("Normal", "Downregulated" , "Upregulated","Shared"))
category_counts <- as.data.frame(table(PD_data$category))
colnames(category_counts) <- c("Category", "Count")
category_counts$Percentage <- round(category_counts$Count / sum(category_counts$Count) * 100, 1)
category_counts$Label <- paste0(category_counts$Category, ": ", category_counts$Percentage, "%")

# 创建饼图
pie_chart <- ggplot(category_counts, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "black") +  # 添加黑色边框
  coord_polar(theta = "y") +  # 转换为饼图
  theme_void() +  # 去除背景和轴
  theme(legend.position = "none") +  # 去除图例
  scale_fill_manual(values = c("Upregulated" = "#ECA861", 
                               "Downregulated" = "#6091C9", 
                               "Shared" = "#A0332D", 
                               "Normal" = "#D8D9DA"))

# 将饼图转换为图层对象
pie_grob <- ggplotGrob(pie_chart)
PD_data<-subset(PD_data,OR<10)
# 计算 neg_log10_Padj 用于 Y 轴
PD_data$neg_log10_Padj <- -log10(PD_data$Padj)
top_genes <- PD_data %>% arrange(desc(OR)) %>% head(5)
# 创建散点图
p <- ggplot(PD_data, aes(x = Beta, y = neg_log10_Padj)) +
  geom_point(aes(fill = category), size = 2.5, alpha = 1, shape = 21, 
             color = "black", stroke = 0.3) +
  scale_fill_manual(values = c("Upregulated" = "#ECA861", 
                               "Downregulated" = "#6091C9", 
                               "Shared" = "#A0332D", 
                               "Normal" = "#D8D9DA")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "pink", linewidth = 0.88) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "pink", linewidth = 0.88) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  labs(x = "Beta", y = expression(-log[10](italic(P)[adj]))) +
  xlim(-0.52, 0.52) + ylim(0, 23) +  # 调整坐标轴范围
  
  # 标记 OR 值最高的 5 个基因，添加加粗字体和白色外框
  geom_text_repel(
    data = top_genes, aes(label = Gene), 
    fontface = "bold",  # 加粗字体
    box.padding = 0.5,  # 标签与点之间的间距
    max.overlaps = 10,  # 控制标签重叠
    segment.color = "black",  # 连线颜色
    segment.size = 0.5,  # 连线宽度
    color = "black",  # 字体颜色
    bg.color = "white",  # 标签背景色
    bg.r = 0.15  # 标签圆角半径
  )

# 将饼图添加到散点图的左上角
p<-p + annotation_custom(pie_grob, xmin = -0.5, xmax = -0.05, ymin = 12.1, ymax = 22)
p
ggsave(
  filename = "16. Cross_sectional/DEP PD.pdf",  # 文件名
  plot = p,  # 图形对象
  width = 6,  # 宽度 (英寸)
  height = 4.3,  # 高度 (英寸)
  dpi = 300,  # 分辨率设置为 300 DPI
  device = cairo_pdf  # 确保 PDF 的高质量输出
)

vals<-Merge_sig_data$Gene

# 计算 -log10(Padj) 用于 Y 轴
MDD_data$neg_log10_Padj <- -log10(MDD_data$Padj)

# 根据 OR 和 Padj 创建颜色分组
MDD_data$category <- with(MDD_data, 
                          ifelse(Gene %in% vals, "Shared", 
                                 ifelse(OR > 1 & Padj < 0.05, "Upregulated", 
                                        ifelse(OR < 1 & Padj < 0.05, "Downregulated", 
                                               "Normal"))))
MDD_data$category<-factor(MDD_data$category, levels = c("Normal", "Downregulated" , "Upregulated","Shared"))
category_counts <- as.data.frame(table(MDD_data$category))
colnames(category_counts) <- c("Category", "Count")
category_counts$Percentage <- round(category_counts$Count / sum(category_counts$Count) * 100, 1)
category_counts$Label <- paste0(category_counts$Category, ": ", category_counts$Percentage, "%")

# 创建饼图
pie_chart <- ggplot(category_counts, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "black") +  # 添加黑色边框
  coord_polar(theta = "y") +  # 转换为饼图
  theme_void() +  # 去除背景和轴
  theme(legend.position = "none") +  # 去除图例
  scale_fill_manual(values = c("Upregulated" = "#CC88B0", 
                               "Downregulated" = "#679DAC", 
                               "Shared" = "#A0332D", 
                               "Normal" = "#D8D9DA"))

# 将饼图转换为图层对象
pie_grob <- ggplotGrob(pie_chart)
MDD_data<-subset(MDD_data,OR<10)
# 计算 neg_log10_Padj 用于 Y 轴
MDD_data$neg_log10_Padj <- -log10(MDD_data$Padj)
top_genes <- MDD_data %>% arrange(desc(OR)) %>% head(5)
# 创建散点图
p <- ggplot(MDD_data, aes(x = Beta, y = neg_log10_Padj)) +
  geom_point(aes(fill = category), size = 2.5, alpha = 1, shape = 21, 
             color = "black", stroke = 0.3) +
  scale_fill_manual(values = c("Upregulated" = "#CC88B0", 
                               "Downregulated" = "#679DAC", 
                               "Shared" = "#A0332D", 
                               "Normal" = "#D8D9DA")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "pink", linewidth = 0.88) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "pink", linewidth = 0.88) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  labs(x = "Beta", y = expression(-log[10](italic(P)[adj]))) +
  xlim(-0.62, 0.62) + ylim(0, 19) +  # 调整坐标轴范围
  
  # 标记 OR 值最高的 5 个基因，添加加粗字体和白色外框
  geom_text_repel(
    data = top_genes, aes(label = Gene), 
    fontface = "bold",  # 加粗字体
    box.padding = 0.5,  # 标签与点之间的间距
    max.overlaps = 10,  # 控制标签重叠
    segment.color = "black",  # 连线颜色
    segment.size = 0.5,  # 连线宽度
    color = "black",  # 字体颜色
    bg.color = "white",  # 标签背景色
    bg.r = 0.15  # 标签圆角半径
  )

# 将饼图添加到散点图的左上角
p<-p + annotation_custom(pie_grob, xmin = -0.59615, xmax = -0.059615, ymin = 9.95276, ymax = 18.173913)
p
ggsave(
  filename = "16. Cross_sectional/DEP MDD.pdf",  # 文件名
  plot = p,  # 图形对象
  width = 6,  # 宽度 (英寸)
  height = 4.3,  # 高度 (英寸)
  dpi = 300,  # 分辨率设置为 300 DPI
  device = cairo_pdf  # 确保 PDF 的高质量输出
)
