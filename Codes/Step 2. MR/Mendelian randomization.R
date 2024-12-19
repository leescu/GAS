# +================================================+ ####
# +====Section 0. Packages and setwd===============+ ####
# +================================================+ #### 
library(TwoSampleMR)
library(data.table)
library(MRPRESSO)
setwd("D:/PD_MDD/Step 0. Data")
# +================================================+ ####
# +====Section 1. Periodontal & Depressive=========+ ####
# +================================================+ ####
load(file ="3. MR/Periodontal_exposure.Rdata")
load(file ="3. MR/Depressive_outcome.Rdata")

# >>>>> Section 1.1. Mendelian randomization <<<<< ####
Depressive_outcome <- subset(Depressive_outcome,eaf.outcome>0.01)
Periodontal_exposure<- subset(Periodontal_exposure,eaf.exposure>0.01)
Data <- harmonise_data(exposure_dat =Periodontal_exposure,
                       outcome_dat =Depressive_outcome)
set.seed(123)
res <- mr(Data)
OR <- exp(res$b)
CI_lower<- exp(res$b - 1.96 *res$se)
CI_upper<- exp(res$b + 1.96 *res$se)

PD_MR_result<-as.data.frame(cbind(res$method,OR,CI_lower,CI_upper,res$pval))
colnames(PD_MR_result)<-c("Name","OR","CI_lower","CI_upper","Pval")
write.csv(PD_MR_result,file="3. MR/PD_MR_result_total.csv")

PD_MR_qval<-mr_heterogeneity(Data)[2,8]


mr_heterogeneity(Data)
single <- mr_leaveoneout(Data)
OR <- exp(single$b)
CI_lower<- exp(single$b - 1.96 *single$se)
CI_upper<- exp(single$b + 1.96 *single$se)

PD_MR_single_result<-as.data.frame(cbind(single$SNP,OR,CI_lower,CI_upper,single$p))
colnames(PD_MR_single_result)<-c("Name","OR","CI_lower","CI_upper","Pval")
PD_MR_single_result<-PD_MR_single_result[1:6,]
write.csv(PD_MR_single_result,file="3. MR/PD_MR_result_single.csv")

res_single <- mr_singlesnp(Data)
res_single<-res_single[1:6,]
res_single$SNP<- paste0("Leave ", res_single$SNP)
OR <- exp(res_single$b)
CI_lower<- exp(res_single$b - 1.96 *res_single$se)
CI_upper<- exp(res_single$b + 1.96 *res_single$se)

PD_MR_leavesingle_result<-as.data.frame(cbind(res_single$SNP,OR,CI_lower,CI_upper,res_single$p))
colnames(PD_MR_leavesingle_result)<-c("Name","OR","CI_lower","CI_upper","Pval")
write.csv(PD_MR_leavesingle_result,file="3. MR/PD_MR_leavesingle_result.csv")
write.csv(res_single,file="3. MR/PD_MR_result_remove.csv")
# >>>>> Section 1.3. MRPRESSO <<<<< ####
set.seed(123)
mr_presso<-mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = Data, NbDistribution = 1000,  
          SignifThreshold = 0.05)
OR <- exp(mr_presso[["Main MR results"]]$`Causal Estimate`)
Data
CI_lower <- exp(mr_presso[["Main MR results"]]$`Causal Estimate`- 1.96 * mr_presso[["Main MR results"]]$Sd)
CI_upper <- exp(mr_presso[["Main MR results"]]$`Causal Estimate`+ 1.96 * mr_presso[["Main MR results"]]$Sd)
PD_mr_presso_result<-as.data.frame(cbind("MR PRESSO",OR,CI_lower,CI_upper,mr_presso[["Main MR results"]][,6,1]))
colnames(PD_mr_presso_result)<-c("Name","OR","CI_lower","CI_upper","Pval")
PD_MR<-rbind(PD_MR_result,PD_mr_presso_result[1,],PD_MR_single_result,PD_MR_leavesingle_result)
PD_MR[ , 2:5] <- lapply(PD_MR[ , 2:5], as.numeric)
# 定义一个函数将数值格式化为4位有效数字，并补齐不足位数
format_numeric <- function(x) {
  sprintf("%.3f", x)
}

format_p_value <- function(x) {
  format(x, scientific = TRUE, digits = 4)
}

PD_MR$OR <- sapply(PD_MR$OR, format_numeric)
PD_MR$CI_lower <- sapply(PD_MR$CI_lower, format_numeric)
PD_MR$CI_upper <- sapply(PD_MR$CI_upper, format_numeric)
PD_MR$Pval <- sapply(PD_MR$Pval, format_p_value)
PD_MR <- PD_MR[!PD_MR$Name %in% c("Simple mode"), ]

# 新增列：显示效应值及95%置信区间范围
PD_MR$CI <- paste0(
  PD_MR$OR, " (", PD_MR$CI_lower, " - ", PD_MR$CI_upper, ")"
)

PD_MR$Name[PD_MR$Name=="Inverse variance weighted"]<-"IVW"
write.csv(PD_MR,file="3. MR/PD_MR_Result.csv")
save(PD_MR,file="3. MR/PD_MR_Result.Rdata")

# >>>>> Section 1.2. pleio <<<<< ####
set.seed(123)
Data <- harmonise_data(exposure_dat =Periodontal_exposure,
                       outcome_dat =Depressive_outcome)
set.seed(123)

# mr_presso[["Main MR results"]][,'Causal Estimate'][1]
# mr_presso[["Main MR results"]][,'Sd'][1]
# mr_presso[["Main MR results"]][,'P-value'][1]
# SE_MR_PRESSO <- mr_presso[["Main MR results"]][,'Sd'][1]/ sqrt(6)

res <- mr(Data)
# res[6,]<-c(res[5,][1:4],"MR PRESSO",6,
#            mr_presso[["Main MR results"]][,'Causal Estimate'][1],
#            mr_presso[["Main MR results"]][,'Sd'][1],
#            mr_presso[["Main MR results"]][,'P-value'][1])
res <- res[res$method != "Simple mode", ]


mr_scatter<-mr_scatter_plot(res,Data)

mrres<-mr_scatter[["w2Cd43.gSJx5f"]][["plot_env"]][["mrres"]]
dat<-mr_scatter[["w2Cd43.gSJx5f"]][["plot_env"]][["d"]]
dat_plot<-dat[,c("beta.exposure","beta.outcome","se.exposure","se.outcome")]
dat_plot$xmin<-dat_plot$beta.exposure - dat_plot$se.exposure
dat_plot$xmax<-dat_plot$beta.exposure + dat_plot$se.exposure
dat_plot$ymin<-dat_plot$beta.outcome - dat_plot$se.outcome
dat_plot$ymax<-dat_plot$beta.outcome + dat_plot$se.outcome
mrres$method[mrres$method=="Inverse variance weighted"]<-"IVW"
# 提取 MR-Egger 和 IVW 的斜率（slope）

line_info <- data.frame(
  method =mrres$method,
  intercept = mrres$a,
  slope = mrres$b,
  color = c("#7da494", "#eab67a", "#884C45", "#d8a0c1"
            #,"#9f8db8"
            ),
  pval = mrres$pval
)
line_info$linetype <- ifelse(line_info$pval >= 0.05, "P >= 0.05", "P < 0.05")
p<-ggplot(dat_plot, aes(x = beta.exposure, y = beta.outcome)) +
  geom_abline(aes(intercept = intercept, slope = slope, color = method,  linetype =linetype), 
              data = line_info, size = 1) +
  # 绘制 SNP 的散点
 
  # 添加 Y 轴方向的误差条（se.outcome）
  geom_errorbar(aes(ymin = ymin, 
                    ymax = ymax), 
                width = 0.0003, color = "black", alpha = 1) +
  # 添加 X 轴方向的误差条（se.exposure）
  geom_errorbarh(aes(xmin =xmin, 
                     xmax =xmax), 
                 height = 0.0003, color = "black", alpha =1) +
  geom_point(size = 3, color = "#4865a9", alpha =1) +
  # 添加回归线：MR-Egger, Weighted Median, IVW, Weighted Mode
  
  # 设置 X 和 Y 轴标签
  labs(
    x = "SNP effect on PD", y = "SNP effect on MDD",
    color = "Method"  # 设置图例标题
  ) +
  # 设置 X 轴和 Y 轴范围
  xlim(0.01, 0.05) +
  ylim(-0.015, 0.02) +
  # 设置主题和图例位置
  theme_minimal() +  # 使用简洁主题
  theme(
   
    panel.grid.minor.x = element_blank(),  # 去除 X 轴次网格线
    panel.grid.minor.y = element_blank(),  # 去除 Y 轴次网格线
    legend.position = "left", 
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
  
  # 自定义线条颜色
  scale_color_manual(values = setNames(line_info$color, line_info$method))
ggsave(
  filename = "3. MR/PD plot.pdf",  # 文件名
  plot = p,  # 图形对象
  width = 5.5,  # 宽度 (英寸)
  height =3.2,  # 高度 (英寸)
  dpi = 300,  # 分辨率设置为 300 DPI
  device = cairo_pdf  # 确保 PDF 的高质量输出
)
# +================================================+ ####
# +====Section 2. Depressive & Periodontal=========+ ####
# +================================================+ ####
load(file ="3. MR/Periodontal_outcome.Rdata")
load(file ="3. MR/Depressive_exposure.Rdata")


# >>>>> Section 1.1. Mendelian randomization <<<<< ####
Periodontal_outcome <- subset(Periodontal_outcome,eaf.outcome>0.01)
Depressive_exposure<- subset(Depressive_exposure,eaf.exposure>0.01)
Data <- harmonise_data(exposure_dat =Depressive_exposure,
                       outcome_dat =Periodontal_outcome)

set.seed(123)
res <- mr(Data)
OR <- exp(res$b)
CI_lower<- exp(res$b - 1.96 *res$se)
CI_upper<- exp(res$b + 1.96 *res$se)

MDD_MR_result<-as.data.frame(cbind(res$method,OR,CI_lower,CI_upper,res$pval))
colnames(MDD_MR_result)<-c("Name","OR","CI_lower","CI_upper","Pval")
write.csv(MDD_MR_result,file="3. MR/MDD_MR_result_total.csv")
mr_heterogeneity(Data)
MDD_MR_qval<-mr_heterogeneity(Data)[2,8]



single <- mr_leaveoneout(Data)
OR <- exp(single$b)
CI_lower<- exp(single$b - 1.96 *single$se)
CI_upper<- exp(single$b + 1.96 *single$se)

MDD_MR_single_result<-as.data.frame(cbind(single$SNP,OR,CI_lower,CI_upper,single$p))
colnames(MDD_MR_single_result)<-c("Name","OR","CI_lower","CI_upper","Pval")
MDD_MR_single_result<-MDD_MR_single_result[1:39,]
write.csv(MDD_MR_single_result,file="3. MR/MDD_MR_result_single.csv")

res_single <- mr_singlesnp(Data)
res_single<-res_single[1:39,]
res_single$SNP<- paste0("Leave ", res_single$SNP)
OR <- exp(res_single$b)
CI_lower<- exp(res_single$b - 1.96 *res_single$se)
CI_upper<- exp(res_single$b + 1.96 *res_single$se)

MDD_MR_leavesingle_result<-as.data.frame(cbind(res_single$SNP,OR,CI_lower,CI_upper,res_single$p))
colnames(MDD_MR_leavesingle_result)<-c("Name","OR","CI_lower","CI_upper","Pval")
write.csv(MDD_MR_leavesingle_result,file="3. MR/MDD_MR_leavesingle_result.csv")
write.csv(res_single,file="3. MR/MDD_MR_result_remove.csv")
# >>>>> Section 1.3. MRPRESSO <<<<< ####
set.seed(123)
mr_presso<-mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                     OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = Data, NbDistribution = 1000,  
                     SignifThreshold = 0.05)
OR <- exp(mr_presso[["Main MR results"]]$`Causal Estimate`)
CI_lower <- exp(mr_presso[["Main MR results"]]$`Causal Estimate`- 1.96 * mr_presso[["Main MR results"]]$Sd)
CI_upper <- exp(mr_presso[["Main MR results"]]$`Causal Estimate`+ 1.96 * mr_presso[["Main MR results"]]$Sd)
MDD_mr_presso_result<-as.data.frame(cbind("MR PRESSO",OR,CI_lower,CI_upper,mr_presso[["Main MR results"]][,6,1]))
colnames(MDD_mr_presso_result)<-c("Name","OR","CI_lower","CI_upper","Pval")
MDD_MR<-rbind(MDD_MR_result,MDD_mr_presso_result[1,],MDD_MR_single_result,MDD_MR_leavesingle_result)
MDD_MR[ , 2:5] <- lapply(MDD_MR[ , 2:5], as.numeric)
# 定义一个函数将数值格式化为4位有效数字，并补齐不足位数
format_numeric <- function(x) {
  sprintf("%.3f", x)
}

format_p_value <- function(x) {
  format(x, scientific = TRUE, digits = 4)
}

MDD_MR$OR <- sapply(MDD_MR$OR, format_numeric)
MDD_MR$CI_lower <- sapply(MDD_MR$CI_lower, format_numeric)
MDD_MR$CI_upper <- sapply(MDD_MR$CI_upper, format_numeric)
MDD_MR$Pval <- sapply(MDD_MR$Pval, format_p_value)
MDD_MR <- MDD_MR[!MDD_MR$Name %in% c("Simple mode"), ]

# 新增列：显示效应值及95%置信区间范围
MDD_MR$CI <- paste0(
  MDD_MR$OR, " (", MDD_MR$CI_lower, " - ", MDD_MR$CI_upper, ")"
)

MDD_MR$Name[MDD_MR$Name=="Inverse variance weighted"]<-"IVW"
write.csv(MDD_MR,file="3. MR/MDD_MR_Result.csv")
save(MDD_MR,file="3. MR/MDD_MR_Result.Rdata")

# >>>>> Section 1.2. pleio <<<<< ####
set.seed(123)
Data <- harmonise_data(exposure_dat =Periodontal_exposure,
                       outcome_dat =Depressive_outcome)
set.seed(123)

# mr_presso[["Main MR results"]][,'Causal Estimate'][1]
# mr_presso[["Main MR results"]][,'Sd'][1]
# mr_presso[["Main MR results"]][,'P-value'][1]
# SE_MR_PRESSO <- mr_presso[["Main MR results"]][,'Sd'][1]/ sqrt(6)

res <- mr(Data)
# res[6,]<-c(res[5,][1:4],"MR PRESSO",6,
#            mr_presso[["Main MR results"]][,'Causal Estimate'][1],
#            mr_presso[["Main MR results"]][,'Sd'][1],
#            mr_presso[["Main MR results"]][,'P-value'][1])
res <- res[res$method != "Simple mode", ]


mr_scatter<-mr_scatter_plot(res,Data)

mrres<-mr_scatter[["iFHCQo.FCT8kd"]][["plot_env"]][["mrres"]]
dat<-mr_scatter[["iFHCQo.FCT8kd"]][["plot_env"]][["d"]]
dat_plot<-dat[,c("beta.exposure","beta.outcome","se.exposure","se.outcome")]
dat_plot$xmin<-dat_plot$beta.exposure - dat_plot$se.exposure
dat_plot$xmax<-dat_plot$beta.exposure + dat_plot$se.exposure
dat_plot$ymin<-dat_plot$beta.outcome - dat_plot$se.outcome
dat_plot$ymax<-dat_plot$beta.outcome + dat_plot$se.outcome
mrres$method[mrres$method=="Inverse variance weighted"]<-"IVW"
mrres
# 提取 MR-Egger 和 IVW 的斜率（slope）

line_info <- data.frame(
  method =mrres$method,
  intercept = mrres$a,
  slope = mrres$b,
  color = c("#7da494", "#eab67a", "#884C45", "#d8a0c1"
            #,"#9f8db8"
  ),
  pval = mrres$pval
)
line_info$linetype <- ifelse(line_info$pval >= 0.05, "P >= 0.05", "P < 0.05")
p<-ggplot(dat_plot, aes(x = beta.exposure, y = beta.outcome)) +
  geom_abline(aes(intercept = intercept, slope = slope, color = method,  linetype =linetype), 
              data = line_info, size = 1, alpha = 0.9) +
  # 绘制 SNP 的散点
  
  # 添加 Y 轴方向的误差条（se.outcome）
  geom_errorbar(aes(ymin = ymin, 
                    ymax = ymax), 
                width = 0.0003, color = "black", alpha = 1) +
  # 添加 X 轴方向的误差条（se.exposure）
  geom_errorbarh(aes(xmin =xmin, 
                     xmax =xmax), 
                 height = 0.0003, color = "black", alpha =1) +
  geom_point(size = 3, color = "#ef8a43", alpha =1) +
  # 添加回归线：MR-Egger, Weighted Median, IVW, Weighted Mode
  
  # 设置 X 和 Y 轴标签
  labs(
    x = "SNP effect on MDD", y = "SNP effect on PD",
    color = "Method"  # 设置图例标题
  ) +
  # 设置 X 轴和 Y 轴范围
  xlim(0.01, 0.09) +
  ylim(-0.03, 0.025) +
  # 设置主题和图例位置
  theme_minimal() +  # 使用简洁主题
  theme(
    
    panel.grid.minor.x = element_blank(),  # 去除 X 轴次网格线
    panel.grid.minor.y = element_blank(),  # 去除 Y 轴次网格线
    legend.position = "left", 
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
  
  # 自定义线条颜色
  scale_color_manual(values = setNames(line_info$color, line_info$method))
p
ggsave(
  filename = "3. MR/MMD plot.pdf",  # 文件名
  plot = p,  # 图形对象
  width = 5.5,  # 宽度 (英寸)
  height =3.2,  # 高度 (英寸)
  dpi = 300,  # 分辨率设置为 300 DPI
  device = cairo_pdf  # 确保 PDF 的高质量输出
)
# +================================================+ ####
# +====Section 3. forest plot PD===================+ ####
# +================================================+ ####
load(file="3. MR/PD_MR_Result.Rdata")
Result<-PD_MR[1:5,]
Result
Result$OR <- as.numeric(Result$OR)
Result$CI_lower <- as.numeric(Result$CI_lower)
Result$CI_upper <- as.numeric(Result$CI_upper)
Result$Name <- factor(Result$Name, levels = c('MR PRESSO', 'Weighted mode', 'IVW','Weighted median' , 'MR Egger'))
Y<-c(4,3,2,1)
# 在每个分面图添加的左上角标签数据
label_data <- data.frame(
  status = c(
    "PSM PD Cross"),
  label = c(
    "Odds Ratio (95%CI)"),
  effect = 0.91,  # 确保文本在x轴的左侧边缘附近
  model = "Effect size"  # 将文本放置在分面图的顶部
)
xlim_range <- c(0.9, 1.7)
# 绘制森林图
p<-ggplot(Result, aes(x = OR, y = Name)) +
  # 绘制截断的误差条，只在需要时添加箭头
  geom_segment(
    aes(
      x = pmax(CI_lower, xlim_range[1]),  # 左界限
      xend = pmin(CI_upper, xlim_range[2]),  # 右界限
      y = Name, yend = Name  # Y 轴位置不变
    ),
    color = "black", size = 0.7
  ) +
  # 添加箭头指示超出范围的误差条
  geom_segment(
    data = Result[Result$CI_lower < xlim_range[1], ],
    aes(
      x = xlim_range[1], xend = xlim_range[1] - 0.05, 
      y = Name, yend = Name
    ),
    color = "black", size = 0.7,
    arrow = arrow(angle = 18, type = "closed", length = unit(0.28, "cm"))
  ) +
  geom_segment(
    data = Result[Result$CI_upper > xlim_range[2], ],
    aes(
      x = xlim_range[2], xend = xlim_range[2] + 0.05, 
      y = Name, yend = Name
    ),
    color = "black", size = 0.7,
    arrow = arrow(angle = 18,type = "closed", length = unit(0.28, "cm"))
  ) +
  geom_segment(data = Result[Result$CI_lower >= xlim_range[1]&Result$CI_upper <= xlim_range[2], ],
    aes(
      x = pmax(CI_lower, xlim_range[1]),
      xend = pmax(CI_lower, xlim_range[1]),
      y =Y - 0.08, yend = Y + 0.08  # 设置竖线高度
    ),
    color = "black", size =1
  ) +
  
  # 添加右侧竖线
  geom_segment(data = Result[Result$CI_lower >= xlim_range[1]&Result$CI_upper <= xlim_range[2], ],
    aes(
      x = pmin(CI_upper, xlim_range[2]),
      xend = pmin(CI_upper, xlim_range[2]),
      y = Y - 0.08, yend = Y + 0.08  # 设置竖线高度
    ),
    color = "black", size = 1
  ) +
  # 绘制散点
  geom_point(size = 3,color="#4865a9") +  # 绘制点，颜色根据 status 区分
  
  theme_minimal() +  # 使用简洁主题
  labs(x = "PD MR effect size",y = "Model") +  # 设置标签
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
  scale_color_manual(values = c("#4865a9")) +  # 分组颜色设置
  scale_x_continuous(breaks = seq(0.9, 1.7, by = 0.1)) +  # 设置 X 轴范围和刻度
  annotate("segment", x = 1, xend = 1, y = 0.5, yend = 5.5, color = "#e53b2c", linewidth = 1) +  # 部分红色垂直线
  geom_segment(data = Result, aes(x = 0.85, xend = 1.7, y = 6, yend = 6), 
               inherit.aes = FALSE, color = "white", linewidth = 0.5) +
  geom_text(data = label_data, aes(x = effect, y = model, label = label), 
            inherit.aes = FALSE, color = "black", size = 4,fontface = "bold",  hjust = 0)
#6*7.5
p
ggsave(
  filename = "3. MR/forest plot PD.pdf",  # 文件名
  plot = p,  # 图形对象
  width = 5.5,  # 宽度 (英寸)
  height = 3.2,  # 高度 (英寸)
  dpi = 300,  # 分辨率设置为 300 DPI
  device = cairo_pdf  # 确保 PDF 的高质量输出
)
# +================================================+ ####
# +====Section 3. forest plot MDD==================+ ####
# +================================================+ ####
load(file="3. MR/MDD_MR_Result.Rdata")
Result<-MDD_MR[1:5,]
Result
Result$OR <- as.numeric(Result$OR)
Result$CI_lower <- as.numeric(Result$CI_lower)
Result$CI_upper <- as.numeric(Result$CI_upper)
Result$Name <- factor(Result$Name, levels = c('MR PRESSO', 'Weighted mode', 'IVW','Weighted median' , 'MR Egger'))
Y<-c(4,3,2,1)
# 在每个分面图添加的左上角标签数据
label_data <- data.frame(
  status = c(
    "PSM MDD Cross"),
  label = c(
    "Odds Ratio (95%CI)"),
  effect = 0.91,  # 确保文本在x轴的左侧边缘附近
  model = "Effect size"  # 将文本放置在分面图的顶部
)
m_range <- c(0.9, 1.7)
# 绘制森林图
p<-ggplot(Result, aes(x = OR, y = Name)) +
  # 绘制截断的误差条，只在需要时添加箭头
  geom_segment(
    aes(
      x = pmax(CI_lower, xlim_range[1]),  # 左界限
      xend = pmin(CI_upper, xlim_range[2]),  # 右界限
      y = Name, yend = Name  # Y 轴位置不变
    ),
    color = "black", size = 0.7
  ) +
  # 添加箭头指示超出范围的误差条
  geom_segment(
    data = Result[Result$CI_lower < xlim_range[1], ],
    aes(
      x = xlim_range[1], xend = xlim_range[1] - 0.05, 
      y = Name, yend = Name
    ),
    color = "black", size = 0.7,
    arrow = arrow(angle = 18, type = "closed", length = unit(0.28, "cm"))
  ) +
  geom_segment(
    data = Result[Result$CI_upper > xlim_range[2], ],
    aes(
      x = xlim_range[2], xend = xlim_range[2] + 0.05, 
      y = Name, yend = Name
    ),
    color = "black", size = 0.7,
    arrow = arrow(angle = 18,type = "closed", length = unit(0.28, "cm"))
  ) +
  geom_segment(data = Result[Result$CI_lower >= xlim_range[1]&Result$CI_upper <= xlim_range[2], ],
               aes(
                 x = pmax(CI_lower, xlim_range[1]),
                 xend = pmax(CI_lower, xlim_range[1]),
                 y =Y - 0.08, yend = Y + 0.08  # 设置竖线高度
               ),
               color = "black", size =1
  ) +
  
  # 添加右侧竖线
  geom_segment(data = Result[Result$CI_lower >= xlim_range[1]&Result$CI_upper <= xlim_range[2], ],
               aes(
                 x = pmin(CI_upper, xlim_range[2]),
                 xend = pmin(CI_upper, xlim_range[2]),
                 y = Y - 0.08, yend = Y + 0.08  # 设置竖线高度
               ),
               color = "black", size = 1
  ) +
  # 绘制散点
  geom_point(size = 3,color="#ef8a43") +  # 绘制点，颜色根据 status 区分
  
  theme_minimal() +  # 使用简洁主题
  labs(x = "MDD MR effect size",y = "Model") +  # 设置标签
  theme(
    panel.grid.major.x = element_blank(),  # 去除 X 轴主网格线
    panel.grid.minor.x = element_blank(),  # 去除 X 轴次网格线
    panel.grid.minor.y = element_blank(),  # 去除 Y 轴次网格线
    legend.position = "right",  # 图例在右侧
    legend.text = element_text(size = 12),  # 增大图例字体大小
    legend.title = element_text(size =14, face = "bold"),  # 图例标题加粗并增大字体
    axis.text.x = element_text(size = 14, color = "black"),  # X 轴刻度字体设为黑色
    axis.text.y = element_text(size = 14, color = "black"),  # Y 轴刻度字体设为黑色
    axis.title = element_text(size = 14, face = "bold", color = "black"),  # 坐标轴标题字体加粗且为黑色
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),  # X 轴刻度线
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),  # Y 轴刻度线
    strip.text = element_text(size = 12, face = "bold"),  # 增大分面标签字体并加粗
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # 添加图的边框
  ) +
  scale_color_manual(values = c("#ef8a43")) +  # 分组颜色设置
  scale_x_continuous(breaks = seq(0.9, 1.7, by = 0.1)) +  # 设置 X 轴范围和刻度
  annotate("segment", x = 1, xend = 1, y = 0.5, yend = 5.5, color = "#e53b2c", linewidth = 1) +  # 部分红色垂直线
  geom_segment(data = Result, aes(x = 0.85, xend = 1.7, y = 6, yend = 6), 
               inherit.aes = FALSE, color = "white", linewidth = 0.5) +
  geom_text(data = label_data, aes(x = effect, y = model, label = label), 
            inherit.aes = FALSE, color = "black", size = 4,fontface = "bold",  hjust = 0)
#6*7.5
p
ggsave(
  filename = "3. MR/forest plot MDD.pdf",  # 文件名
  plot = p,  # 图形对象
  width = 5.5,  # 宽度 (英寸)
  height = 3.2,  # 高度 (英寸)
  dpi = 300,  # 分辨率设置为 300 DPI
  device = cairo_pdf  # 确保 PDF 的高质量输出
)

