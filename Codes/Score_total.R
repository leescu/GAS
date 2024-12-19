# +================================================+ ####
# +====Section 0. Packages and setwd===============+ ####
# +================================================+ #### 
setwd("D:/PD_MDD/Step 0. Data")
library(data.table)
library(TwoSampleMR)
library(devtools)
library(tidyverse)
library(dplyr)
library(tidyr)
library(topsis)
library(biomaRt)
library(corrplot)
library(grid)
# +================================================+ ####
# +====Section 1. FUMA(PLACO) =====================+ ####
# +================================================+ #### 
FUMA_original<-fread("13. FUMA/FUMA_1e6/genes.txt")
FUMA_original<-FUMA_original[na.omit(FUMA_original$posMapMaxCADD)&FUMA_original$posMapMaxCADD!=0&FUMA_original$type=="protein_coding"]

# 连接到 Ensembl 数据库
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 获取 Entrez ID 到 Gene Symbol 的映射
entrez_result <- getBM(attributes = c("entrezgene_id", "hgnc_symbol"),
                       filters = "entrezgene_id",
                       values = FUMA_original$entrezID,
                       mart = ensembl)


# 获取 ENSG ID 到 Gene Symbol 的映射
ensg_result <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id",
                     values =FUMA_original$ensg,
                     mart = ensembl)

# 重命名列以便合并
colnames(entrez_result) <- c("id", "symbol_entrez")
colnames(ensg_result) <- c("id", "symbol_ensg")

entrez_data <- merge(FUMA_original, entrez_result, by.x = "entrezID",
                     by.y = "id", all.x = TRUE)
ensg_data <- merge(entrez_data, ensg_result, by.x = "ensg",
                   by.y = "id", all.x = TRUE)

final_data <- ensg_data %>%
  mutate(symbol = coalesce(symbol_entrez, symbol_ensg))


FUMA<-final_data[,c("symbol","posMapMaxCADD","minGwasP")]
FUMA<-na.omit(FUMA)
FUMA<- FUMA %>%
  group_by(symbol) %>%  # 按 symbol 分组
  slice_max(posMapMaxCADD, n = 1, with_ties = FALSE) %>%  # 取最大值
  ungroup()  # 取消分组
colnames(FUMA)
FUMA$PLACO_score <- 0.05 +(FUMA$posMapMaxCADD - min(FUMA$posMapMaxCADD)) / (max(FUMA$posMapMaxCADD) - min(FUMA$posMapMaxCADD))*(1 - 0.05)
FUMA$FUMA_P<-FUMA$minGwasP
PLACO_data<-FUMA[,c("symbol","PLACO_score")]


save(PLACO_data,file="23. Score_total/1. PLACO_data.Rdata")
# +================================================+ ####
# +====Section 2. CTMA ============================+ ####
# +================================================+ #### 
FUMA_original<-fread("12. CTMA/FUMA_1e6/genes.txt")
FUMA_original<-FUMA_original[na.omit(FUMA_original$posMapMaxCADD)&FUMA_original$posMapMaxCADD!=0&FUMA_original$type=="protein_coding"]

# 连接到 Ensembl 数据库
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 获取 Entrez ID 到 Gene Symbol 的映射
entrez_result <- getBM(attributes = c("entrezgene_id", "hgnc_symbol"),
                       filters = "entrezgene_id",
                       values = FUMA_original$entrezID,
                       mart = ensembl)

# 获取 ENSG ID 到 Gene Symbol 的映射
ensg_result <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id",
                     values =FUMA_original$ensg,
                     mart = ensembl)

# 重命名列以便合并
colnames(entrez_result) <- c("id", "symbol_entrez")
colnames(ensg_result) <- c("id", "symbol_ensg")

entrez_data <- merge(FUMA_original, entrez_result, by.x = "entrezID",
                     by.y = "id", all.x = TRUE)
ensg_data <- merge(entrez_data, ensg_result, by.x = "ensg",
                   by.y = "id", all.x = TRUE)
final_data <- ensg_data %>%
  mutate(symbol = coalesce(symbol_entrez, symbol_ensg)) 


FUMA<-final_data[,c("symbol","posMapMaxCADD","minGwasP")]
FUMA<-na.omit(FUMA)
FUMA<- FUMA %>%
  group_by(symbol) %>%  # 按 symbol 分组
  slice_max(posMapMaxCADD, n = 1, with_ties = FALSE) %>%  # 取最大值
  ungroup()  # 取消分组
colnames(FUMA)
FUMA$CTMA_score <- 0.05 +(FUMA$posMapMaxCADD - min(FUMA$posMapMaxCADD)) / (max(FUMA$posMapMaxCADD) - min(FUMA$posMapMaxCADD))*(1 - 0.05)
FUMA$CTMA_P<-FUMA$minGwasP
CTMA_data<-FUMA[,c("symbol","CTMA_score")]


save(CTMA_data,file="23. Score_total/2. CTMA_data.Rdata")
# +================================================+ ####
# +====Section 3. MR-JTI ==========================+ ####
# +================================================+ #### 
target_path <- "D:/PD_MDD/Step 0. Data/17. MR_JTI/JTI"  # 将此处替换为您的目标文件夹路径
files <- list.files(path = target_path, pattern = "MDD_JTI_.*\\.txt", full.names = TRUE)
 # 使用 lapply() 读取每个文件，并将它们存储在一个列表中
data_list <- lapply(files, read.table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
 # 使用 do.call() 和 rbind() 将所有数据框纵向合并为一个数据框
MDD_original<- do.call(rbind, data_list)
MDD<-subset(MDD_original,pvalue<0.05)
MDD <- MDD %>%
  group_by(genename) %>%  # 按 symbol 分组
  slice_max(abs(effect_size), n = 1, with_ties = FALSE) %>%  # 取绝对值最大的值
  ungroup()
MDD$MDD_effect_size<-MDD$effect_size

files <- list.files(path = target_path, pattern = "PD_JTI_.*\\.txt", full.names = TRUE)
 # 使用 lapply() 读取每个文件，并将它们存储在一个列表中
data_list <- lapply(files, read.table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
 # 使用 do.call() 和 rbind() 将所有数据框纵向合并为一个数据框
PD_original<- do.call(rbind, data_list)
PD<-subset(PD_original,pvalue<0.05)
PD <- PD %>%
   group_by(genename) %>%  # 按 symbol 分组
   slice_max(abs(effect_size), n = 1, with_ties = FALSE) %>%  # 取绝对值最大的值
   ungroup()
PD$PD_effect_size<-PD$effect_size
PD_MDD<-merge(MDD[,c("gene","MDD_effect_size")],PD[,c("gene","PD_effect_size")],by="gene",all=F)
PD_MDD$symbol<-PD_MDD$gene

PD_MDD<-PD_MDD  %>%
  mutate(effect_size = sqrt(MDD_effect_size * PD_effect_size))
PD_MDD <-na.omit(PD_MDD)
PD_MDD$JTI <- 0.05 +(PD_MDD$effect_size - min(PD_MDD$effect_size)) / (max(PD_MDD$effect_size) - min(PD_MDD$effect_size))*(1 - 0.05)

PD_MDD <- PD_MDD %>%
  mutate(JTI_score = if_else(MDD_effect_size < 0, -JTI, JTI))
JTI_data<-PD_MDD[,c("symbol","JTI_score")]
save(JTI_data,file="23. Score_total/3. JTI_data.Rdata")

# +================================================+ ####
# +====Section 4. DEGs ============================+ ####
# +================================================+ #### 
MDD_original<-read.table(file ="20. DEG/MDD_DEG.csv", header = T, sep =",")
MDD_original<-subset(MDD_original,P.Value<0.05)
MDD_original$MDD_effect_size<-MDD_original$logFC
MDD_original$symbol<-rownames(MDD_original)
PD_original<-read.table(file ="20. DEG/PD_DEG.csv", header = T, sep =",")
PD_original<-subset(PD_original,P.Value<0.05)
PD_original$PD_effect_size<-PD_original$logFC
PD_original$symbol<-rownames(PD_original)
PD_MDD<-merge(MDD_original[,c("symbol","MDD_effect_size")],PD_original[,c("symbol","PD_effect_size")],by="symbol",all=F)

PD_MDD <-PD_MDD %>%
  mutate(effect_size = sqrt(MDD_effect_size * PD_effect_size))
PD_MDD <-na.omit(PD_MDD) 
PD_MDD$DEG <- 0.05 +(PD_MDD$effect_size - min(PD_MDD$effect_size)) / (max(PD_MDD$effect_size) - min(PD_MDD$effect_size))*(1 - 0.05)

PD_MDD <- PD_MDD %>%
  mutate(DEG_score = if_else(MDD_effect_size < 0, -DEG, DEG))
DEG_data<-PD_MDD[,c("symbol","DEG_score")]
save(DEG_data,file="23. Score_total/4. DEG_data.Rdata")

# +================================================+ ####
# +====Section 5. WGCNA ===========================+ ####
# +================================================+ #### 
MDD_original_blue<-read.table(file ="21. WGCNA/MDD_blue_module.csv", header = T, sep =",")
MDD_original_turquoise<-read.table(file ="21. WGCNA/MDD_turquoise_module.csv", header = T, sep =",")

#MDD_original<-subset(MDD_original,p.MM.blue<0.05)
MDD_original_blue$MDD_effect_size<-MDD_original_blue$GS.Group
MDD_original_blue$symbol<-MDD_original_blue$Gene_symbol

MDD_original_turquoise$MDD_effect_size<-MDD_original_turquoise$GS.Group
MDD_original_turquoise$symbol<-MDD_original_turquoise$Gene_symbol
MDD_original<-rbind(MDD_original_turquoise,MDD_original_blue)
# 
# MDD_cleaned <- MDD_original %>%
#   group_by(symbol) %>% 
#   summarise(
#     MDD_effect_size = if (any(MDD_effect_size > 0) & any(MDD_effect_size < 0)) {
#       # 如果有一个 MDD_effect_size > 0，另一个 < 0，选择保留绝对值更大的
#       MDD_effect_size[which.max(abs(MDD_effect_size))]
#     } else if (all(MDD_effect_size > 0)) {
#       # 如果所有 MDD_effect_size 都大于 0，保留最大的
#       max(MDD_effect_size)
#     } else if (all(MDD_effect_size < 0)) {
#       # 如果所有 MDD_effect_size 都小于 0，保留最小的
#       min(MDD_effect_size)
#     } else {
#       # 如果只存在单一情况，保留唯一值
#       MDD_effect_size
#     }
#   ) %>%
#   ungroup()

PD_original<-read.table(file ="21. WGCNA/PD_yellow_module.csv", header = T, sep =",")
#PD_original<-subset(PD_original,p.MM.yellow<0.05)
PD_original$PD_effect_size<-PD_original$GS.Group
PD_original$symbol<-PD_original$Gene_symbol
PD_MDD<-merge(MDD_original[,c("symbol","MDD_effect_size")],PD_original[,c("symbol","PD_effect_size")],by="symbol",all=F)

PD_MDD <-PD_MDD %>%
  mutate(effect_size = sqrt(MDD_effect_size * PD_effect_size))
PD_MDD <-na.omit(PD_MDD) 
PD_MDD$WGCNA <- 0.05 +(PD_MDD$effect_size - min(PD_MDD$effect_size)) / (max(PD_MDD$effect_size) - min(PD_MDD$effect_size))*(1 - 0.05)

PD_MDD <- PD_MDD %>%
  mutate(WGCNA_score = if_else(MDD_effect_size < 0, -WGCNA, WGCNA))
WGCNA_data<-PD_MDD[,c("symbol","WGCNA_score")]
save(WGCNA_data,file="23. Score_total/5. WGCNA_data.Rdata")

# +================================================+ ####
# +====Section 6. eQTL ============================+ ####
# +================================================+ #### 
MDD_original_blood<-read.table(file ="16. SMR/MDD/trait_eSMR.merged_blood.tsv", header = T, sep ="\t")
MDD_original_brain<-read.table(file ="16. SMR/MDD/trait_eSMR.merged_brain.tsv", header = T, sep ="\t")
MDD_original<-rbind(MDD_original_blood[,c("Gene","probeID","b_eQTL","p_SMR","p_SMR","p_HEIDI")],
                    MDD_original_brain[,c("Gene","probeID","b_eQTL","p_SMR","p_SMR","p_HEIDI")])
colnames(MDD_original)
MDD_original<-subset(MDD_original,p_SMR<0.05&p_HEIDI>0.05)
MDD_original <- MDD_original%>%
  group_by(Gene) %>%  # 按 symbol 分组
  slice_max(abs(b_eQTL), n = 1, with_ties = FALSE) %>%  # 取绝对值最大的值
  ungroup()
MDD_original$MDD_effect_size<-MDD_original$b_eQTL
MDD_original$symbol<-MDD_original$Gene

PD_original_blood<-read.table(file ="16. SMR/PD/trait_eSMR.merged_blood.tsv", header = T, sep ="\t")
PD_original_brain<-read.table(file ="16. SMR/PD/trait_eSMR.merged_brain.tsv", header = T, sep ="\t")
PD_original<-rbind(PD_original_blood[,c("Gene","probeID","b_eQTL","p_SMR","p_SMR","p_HEIDI")],
                   PD_original_brain[,c("Gene","probeID","b_eQTL","p_SMR","p_SMR","p_HEIDI")])
colnames(PD_original)
PD_original<-subset(PD_original,p_SMR<0.05&p_HEIDI>0.05)
PD_original <- PD_original%>%
  group_by(Gene) %>%  # 按 symbol 分组
  slice_max(abs(b_eQTL), n = 1, with_ties = FALSE) %>%  # 取绝对值最大的值
  ungroup()
PD_original$PD_effect_size<-PD_original$b_eQTL
PD_original$symbol<-PD_original$Gene
PD_MDD<-merge(MDD_original[,c("symbol","MDD_effect_size")],PD_original[,c("symbol","PD_effect_size")],by="symbol",all=F)


PD_MDD <-PD_MDD %>%
  mutate(effect_size = sqrt(MDD_effect_size * PD_effect_size))
PD_MDD <-na.omit(PD_MDD) 
PD_MDD$eQTL <- 0.05 +(PD_MDD$effect_size - min(PD_MDD$effect_size)) / (max(PD_MDD$effect_size) - min(PD_MDD$effect_size))*(1 - 0.05)

PD_MDD <- PD_MDD %>%
  mutate(eQTL_score = if_else(MDD_effect_size < 0, -eQTL, eQTL))
eQTL_data<-PD_MDD[,c("symbol","eQTL_score")]
save(eQTL_data,file="23. Score_total/6. eQTL_data.Rdata")

# +================================================+ ####
# +====Section 7. pQTL ============================+ ####
# +================================================+ #### 
MDD_original<-read.table(file ="16. SMR/MDD/trait_pSMR.merged_blood.tsv", header = T, sep ="\t")
colnames(MDD_original)
MDD_original<-subset(MDD_original,p_SMR<0.05&p_HEIDI>0.05)

MDD_original$MDD_effect_size<-MDD_original$b_SMR
MDD_original$symbol<-MDD_original$index

MDD_original <- MDD_original %>%
  group_by(symbol) %>%  # 按 symbol 分组
  slice_max(abs(b_SMR), n = 1, with_ties = FALSE) %>%  # 取绝对值最大的值
  ungroup()

PD_original<-read.table(file ="16. SMR/PD/trait_pSMR.merged_blood.tsv", header = T, sep ="\t")
PD_original<-subset(PD_original,p_SMR<0.05&p_HEIDI>0.05)
PD_original$PD_effect_size<-PD_original$b_SMR
PD_original$symbol<-PD_original$Gene
PD_original <- PD_original %>%
  group_by(symbol) %>%  # 按 symbol 分组
  slice_max(abs(b_SMR), n = 1, with_ties = FALSE) %>%  # 取绝对值最大的值
  ungroup()

PD_MDD<-merge(MDD_original[,c("symbol","MDD_effect_size")],PD_original[,c("symbol","PD_effect_size")],by="symbol",all=F)


PD_MDD <-PD_MDD %>%
  mutate(effect_size = sqrt(MDD_effect_size * PD_effect_size))

PD_MDD <-na.omit(PD_MDD) 
PD_MDD$pQTL <- 0.05 +(PD_MDD$effect_size - min(PD_MDD$effect_size)) / (max(PD_MDD$effect_size) - min(PD_MDD$effect_size))*(1 - 0.05)

PD_MDD <- PD_MDD %>%
  mutate(pQTL_score = if_else(MDD_effect_size < 0, -pQTL, pQTL))
pQTL_data<-PD_MDD[,c("symbol","pQTL_score")]
save(pQTL_data,file="23. Score_total/7. pQTL_data.Rdata")

# +================================================+ ####
# +====Section 8. mQTL ============================+ ####
# +================================================+ #### 
MDD_original_blood<-read.table(file ="16. SMR/MDD/trait_mSMR.merged_blood.tsv", header = T, sep ="\t")
MDD_original_brain<-read.table(file ="16. SMR/MDD/trait_mSMR.merged_brain.tsv", header = T, sep ="\t")
MDD_original<-rbind(MDD_original_blood[,c("Gene","probeID","b_eQTL","p_SMR","p_SMR","p_HEIDI")],
                    MDD_original_brain[,c("Gene","probeID","b_eQTL","p_SMR","p_SMR","p_HEIDI")])
colnames(MDD_original)
MDD_original<-subset(MDD_original,p_SMR<0.05&p_HEIDI>0.05)
MDD_original <- MDD_original%>%
  group_by(Gene) %>%  # 按 symbol 分组
  slice_max(abs(b_eQTL), n = 1, with_ties = FALSE) %>%  # 取绝对值最大的值
  ungroup()
MDD_original$MDD_effect_size<-MDD_original$b_eQTL
MDD_original$symbol<-MDD_original$Gene

PD_original_blood<-read.table(file ="16. SMR/PD/trait_mSMR.merged_blood.tsv", header = T, sep ="\t")
PD_original_brain<-read.table(file ="16. SMR/PD/trait_mSMR.merged_brain.tsv", header = T, sep ="\t")
PD_original<-rbind(PD_original_blood[,c("Gene","probeID","b_eQTL","p_SMR","p_SMR","p_HEIDI")],
                   PD_original_brain[,c("Gene","probeID","b_eQTL","p_SMR","p_SMR","p_HEIDI")])
colnames(PD_original)
PD_original<-subset(PD_original,p_SMR<0.05&p_HEIDI>0.05)
PD_original <- PD_original%>%
  group_by(Gene) %>%  # 按 symbol 分组
  slice_max(abs(b_eQTL), n = 1, with_ties = FALSE) %>%  # 取绝对值最大的值
  ungroup()
PD_original$PD_effect_size<-PD_original$b_eQTL
PD_original$symbol<-PD_original$Gene
PD_MDD<-merge(MDD_original[,c("symbol","MDD_effect_size")],PD_original[,c("symbol","PD_effect_size")],by="symbol",all=F)


PD_MDD <-PD_MDD %>%
  mutate(effect_size = sqrt(MDD_effect_size * PD_effect_size))

PD_MDD <-na.omit(PD_MDD) 
PD_MDD$mQTL <- 0.05 +(PD_MDD$effect_size - min(PD_MDD$effect_size)) / (max(PD_MDD$effect_size) - min(PD_MDD$effect_size))*(1 - 0.05)

PD_MDD <- PD_MDD %>%
  mutate(mQTL_score = if_else(MDD_effect_size < 0, -mQTL, mQTL))
mQTL_data<-PD_MDD[,c("symbol","mQTL_score")]
save(mQTL_data,file="23. Score_total/8. mQTL_data.Rdata")

# +================================================+ ####
# +====Section 8. sQTL ============================+ ####
# +================================================+ #### 
MDD_original_blood<-read.table(file ="16. SMR/MDD/trait_sSMR.merged_blood.tsv", header = T, sep ="\t")
MDD_original_brain<-read.table(file ="16. SMR/MDD/trait_sSMR.merged_brain.tsv", header = T, sep ="\t")
MDD_original<-rbind(MDD_original_blood[,c("Gene","probeID","b_eQTL","p_SMR","p_SMR","p_HEIDI")],
                    MDD_original_brain[,c("Gene","probeID","b_eQTL","p_SMR","p_SMR","p_HEIDI")])
colnames(MDD_original)
MDD_original<-subset(MDD_original,p_SMR<0.05&p_HEIDI>0.05)
MDD_original <- MDD_original%>%
  group_by(Gene) %>%  # 按 symbol 分组
  slice_max(abs(b_eQTL), n = 1, with_ties = FALSE) %>%  # 取绝对值最大的值
  ungroup()
MDD_original$MDD_effect_size<-MDD_original$b_eQTL
MDD_original$symbol<-MDD_original$Gene

PD_original_blood<-read.table(file ="16. SMR/PD/trait_sSMR.merged_blood.tsv", header = T, sep ="\t")
PD_original_brain<-read.table(file ="16. SMR/PD/trait_sSMR.merged_brain.tsv", header = T, sep ="\t")
PD_original<-rbind(PD_original_blood[,c("Gene","probeID","b_eQTL","p_SMR","p_SMR","p_HEIDI")],
                   PD_original_brain[,c("Gene","probeID","b_eQTL","p_SMR","p_SMR","p_HEIDI")])
colnames(PD_original)
PD_original<-subset(PD_original,p_SMR<0.05&p_HEIDI>0.05)
PD_original <- PD_original%>%
  group_by(Gene) %>%  # 按 symbol 分组
  slice_max(abs(b_eQTL), n = 1, with_ties = FALSE) %>%  # 取绝对值最大的值
  ungroup()
PD_original$PD_effect_size<-PD_original$b_eQTL
PD_original$symbol<-PD_original$Gene
PD_MDD<-merge(MDD_original[,c("symbol","MDD_effect_size")],PD_original[,c("symbol","PD_effect_size")],by="symbol",all=F)


PD_MDD <-PD_MDD %>%
  mutate(effect_size = sqrt(MDD_effect_size * PD_effect_size))

PD_MDD <-na.omit(PD_MDD) 
PD_MDD$sQTL <- 0.05 +(PD_MDD$effect_size - min(PD_MDD$effect_size)) / (max(PD_MDD$effect_size) - min(PD_MDD$effect_size))*(1 - 0.05)

PD_MDD <- PD_MDD %>%
  mutate(sQTL_score = if_else(MDD_effect_size < 0, -sQTL, sQTL))
sQTL_data<-PD_MDD[,c("symbol","sQTL_score")]
save(sQTL_data,file="23. Score_total/9. sQTL_data.Rdata")

# +================================================+ ####
# +====Section 10. DEPs ============================+ ####
# +================================================+ #### 
load(file="2. Clinical/MDD_glm_protein.Rdata")
colnames(MDD_glm_protein)
MDD_original<-subset(MDD_glm_protein,P<0.05)
MDD_original$MDD_effect_size<-log(MDD_original$OR)
MDD_original$symbol<-toupper(rownames(MDD_original))
load(file="2. Clinical/PD_glm_protein.Rdata")
PD_original<-subset(PD_glm_protein,P<0.05)
PD_original$PD_effect_size<-log(PD_original$OR)
PD_original$symbol<-toupper(rownames(PD_original))
PD_MDD<-merge(MDD_original[,c("symbol","MDD_effect_size")],PD_original[,c("symbol","PD_effect_size")],by="symbol",all=F)
PD_MDD <-PD_MDD %>%
  mutate(effect_size = sqrt(MDD_effect_size * PD_effect_size))

PD_MDD <-na.omit(PD_MDD) 
PD_MDD$DEP <- 0.05 +(PD_MDD$effect_size - min(PD_MDD$effect_size)) / (max(PD_MDD$effect_size) - min(PD_MDD$effect_size))*(1 - 0.05)

PD_MDD <- PD_MDD %>%
  mutate(DEP_score = if_else(MDD_effect_size < 0, -DEP, DEP))
DEP_data<-PD_MDD[,c("symbol","DEP_score")]
save(DEP_data,file="23. Score_total/10. DEP_data.Rdata")

# +================================================+ ####
# +====Section 11. Mediation ======================+ ####
# +================================================+ #### 
MDD_original<-read.table(file ="2. Clinical/mediation_results_Cross_MDD.csv", header = T, sep =",")
colnames(MDD_original)
MDD_original<-subset(MDD_original,Prop_Mediated_p_value<0.05&ACME_p_value<0.05&Total_Effect_p_value<0.05)
MDD_original$MDD_effect_size<-MDD_original$Prop_Mediated
MDD_original$symbol<-toupper(MDD_original$Mediator)
PD_original<-read.table(file ="2. Clinical/mediation_results_Cross_PD.csv", header = T, sep =",")
colnames(PD_original)
PD_original<-subset(PD_original,Prop_Mediated_p_value<0.05&ACME_p_value<0.05&Total_Effect_p_value<0.05)
PD_original$PD_effect_size<-PD_original$Prop_Mediated
PD_original$symbol<-toupper(PD_original$Mediator)
PD_MDD<-merge(MDD_original[,c("symbol","MDD_effect_size")],PD_original[,c("symbol","PD_effect_size")],by="symbol",all=F)
PD_MDD <-PD_MDD %>%
  mutate(effect_size = sqrt(MDD_effect_size * PD_effect_size))

PD_MDD <-na.omit(PD_MDD) 
PD_MDD$Med<- 0.05 +(PD_MDD$effect_size - min(PD_MDD$effect_size)) / (max(PD_MDD$effect_size) - min(PD_MDD$effect_size))*(1 - 0.05)

PD_MDD <- PD_MDD %>%
  mutate(Med_score = if_else(MDD_effect_size < 0, -Med, Med))
Med_data<-PD_MDD[,c("symbol","Med_score")]
save(Med_data,file="23. Score_total/11. Med_data.Rdata")

# +================================================+ ####
# +====Section 11. Total ==========================+ ####
# +================================================+ #### 
load(file="23. Score_total/1. PLACO_data.Rdata")
#load(file="23. Score_total/2. CTMA_data.Rdata")
load(file="23. Score_total/3. JTI_data.Rdata")
load(file="23. Score_total/4. DEG_data.Rdata")
load(file="23. Score_total/5. WGCNA_data.Rdata")
load(file="23. Score_total/6. eQTL_data.Rdata")
load(file="23. Score_total/7. pQTL_data.Rdata")
load(file="23. Score_total/8. mQTL_data.Rdata")

load(file="23. Score_total/9. sQTL_data.Rdata")
load(file="23. Score_total/10. DEP_data.Rdata")
load(file="23. Score_total/11. Med_data.Rdata")
Merged_data <- PLACO_data %>%
  #full_join(CTMA_data, by = "symbol") %>%
  full_join(JTI_data, by = "symbol") %>%
  full_join(DEG_data, by = "symbol") %>%
  full_join(WGCNA_data, by = "symbol") %>%
  full_join(eQTL_data, by = "symbol") %>%
  full_join(pQTL_data, by = "symbol") %>%
  full_join(mQTL_data, by = "symbol") %>%
  full_join(sQTL_data, by = "symbol") %>%
  full_join(DEP_data, by = "symbol") %>%
  full_join(Med_data, by = "symbol")
non_na_counts <- colSums(!is.na(Merged_data[, c("PLACO_score","JTI_score", 
                                                "DEG_score",  "WGCNA_score", 
                                                "eQTL_score", "pQTL_score", "mQTL_score", "sQTL_score", 
                                                "DEP_score", "Med_score")]))

#Merged_data[is.na(Merged_data)]<-0
Merged_data <- Merged_data %>%
  mutate(Mean_score = rowMeans(across(c(PLACO_score, JTI_score, DEG_score, WGCNA_score,
                                        eQTL_score, pQTL_score, mQTL_score, sQTL_score, DEP_score, Med_score)), 
                               na.rm = TRUE))

Merged_data <- Merged_data %>%
  rowwise() %>%
  mutate(Up = sum(c_across(c(PLACO_score, JTI_score, DEG_score, WGCNA_score,
                             eQTL_score, pQTL_score, mQTL_score,  sQTL_score, DEP_score, Med_score)) > 0, na.rm = TRUE),
         Down = sum(c_across(c(JTI_score, DEG_score, WGCNA_score,
                               eQTL_score, pQTL_score, mQTL_score, sQTL_score, DEP_score, Med_score)) < 0, na.rm = TRUE),
         Sum  = sum(c_across(c(PLACO_score,JTI_score, DEG_score, WGCNA_score,
                               eQTL_score, pQTL_score, mQTL_score, sQTL_score, DEP_score, Med_score)) != 0, na.rm = TRUE) ) %>%
  ungroup()



# Merged_data$Up_score<-((Merged_data$PLACO_score+Merged_data$CTMA_score)+
#                         (abs(Merged_data$DEG_score)+abs(Merged_data$WGCNA_score)+
#                            abs(Merged_data$eQTL_score)+abs(Merged_data$pQTL_score)+abs(Merged_data$mQTL_score)+
#                            abs(Merged_data$DEP_score)+abs(Merged_data$Med_score)))


# Merged_data$Total_score<-((Merged_data$PLACO_score+Merged_data$PLACO_score)+
#                          (Merged_data$DEG_score+Merged_data$WGCNA_score+Merged_data$JTI_score+
#                           Merged_data$eQTL_score+Merged_data$pQTL_score+Merged_data$mQTL_score+
#                           Merged_data$DEP_score+Merged_data$Med_score)/(Merged_data$Up+0.001))

Merged_data[is.na(Merged_data)]<-0
data_matrix <- as.matrix(Merged_data[, c("PLACO_score",  "JTI_score", 
                                         "DEG_score", "WGCNA_score", 
                                         "eQTL_score", "pQTL_score", "mQTL_score","sQTL_score",
                                         "DEP_score", "Med_score")])


total_non_na <- sum(non_na_counts)


log_weights <- log(non_na_counts + 1)

# 归一化处理，确保权重和为 1
adjusted_weights <- log_weights / sum(log_weights)
# 对初步权重进行开根号处理
# +=============UP======================+ #### 
# 定义是否为收益属性（'+'表示收益属性，'-'表示成本属性）
impacts <- rep('+', ncol(data_matrix))  # 假设所有属性都是收益属性

# 进行 TOPSIS 分析，得到每个基因的评分
topsis_result <- topsis(data_matrix, adjusted_weights, impacts)

# 将评分结果加入到原始数据框中
Merged_data <- Merged_data %>%
  mutate(TOPSIS_score_up = topsis_result$rank)

# 查看结果
Merged_data$TOPSIS_weight_up<-Merged_data$TOPSIS_score_up/10^(Merged_data$Up-Merged_data$Down+1)
Merged_data$TOPSIS_score_weight_up<-rank(Merged_data$TOPSIS_weight_up)
Merged_data_up<-subset(Merged_data,Up>0&Mean_score>0)
# +=============DOWN======================+ #### 
# 定义是否为收益属性（'+'表示收益属性，'-'表示成本属性）
impacts <- c('+',rep('-', ncol(data_matrix)-1) ) # 假设所有属性都是收益属性

# 进行 TOPSIS 分析，得到每个基因的评分
topsis_result <- topsis(data_matrix, adjusted_weights, impacts)

# 将评分结果加入到原始数据框中
Merged_data <- Merged_data %>%
  mutate(TOPSIS_score_down = topsis_result$rank)

# 查看结果
print(Merged_data)
Merged_data$TOPSIS_weight_down<-Merged_data$TOPSIS_score_down/10^(Merged_data$Down-Merged_data$Up+1)
Merged_data$TOPSIS_score_weight_down<-rank(Merged_data$TOPSIS_weight_down)
Merged_data_down<-subset(Merged_data,Down>0&Mean_score<0)

Merged_data_down <- Merged_data_down[order(Merged_data_down$TOPSIS_score_weight_down), ]
Merged_data_up <- Merged_data_up[order(Merged_data_up$TOPSIS_score_weight_up), ]
write.csv(Merged_data_down, file = "23. Score_total/Merged_data_down.csv")
write.csv(Merged_data_up, file = "23. Score_total/Merged_data_up.csv")
write.csv(Merged_data, file = "23. Score_total/Merged_data.csv")


# +================================================+ ####
# +====Section 12. Plot ===========================+ ####
# +================================================+ #### 
Merged_data<-read.table("23. Score_total/Merged_data.csv",sep = ",",header = T,row.names = 1)

data_subset <- Merged_data[, c("DEP_score", "Med_score","PLACO_score", 
                               "JTI_score", "DEG_score", 
                               "WGCNA_score", "eQTL_score", "pQTL_score",  "sQTL_score", 
                               "mQTL_score" )]

colnames(data_subset)<-c("DEPs score", "Mediation score","PLACO score", 
                         "JTI score", "DEGs score", 
                         "WGCNA score", "eQTL score", "pQTL score",  "sQTL score", 
                         "mQTL score" )
Group<-as.data.frame(colnames(data_subset))
colnames(Group)<-"Score"
Group$Group<-c("Clinical data","Clinical data","GWAS data","GWAS data",
               "RNA-seq data","RNA-seq data",
               "QTL data","QTL data","QTL data","QTL data")
coltext <- c( "#A69900", "#8F1D00", "#056258", "#6A478F")
colbar <- colorRampPalette(c("#67001F", "#B2182B","#FFFFFF", "#2166AC", "#053061"))
Group$Color <- coltext[factor(Group$Group, 
                              levels = c("Clinical data","GWAS data", "RNA-seq data", "QTL data"))]
# 计算斯皮尔曼相关性系数矩阵
cor_matrix_spearman <- cor(data_subset, use = "pairwise.complete.obs", method = "spearman")
df_summary <- data_subset %>%
  mutate(`Clinical data`= rowSums(across(Group[Group$Group == "Clinical data", ]$Score))) %>%
  mutate(`GWAS data`= rowSums(across(Group[Group$Group == "GWAS data", ]$Score))) %>%
  mutate(`RNA-seq data`= rowSums(across(Group[Group$Group == "RNA-seq data", ]$Score))) %>%
  mutate(`QTL data`= rowSums(across(Group[Group$Group == "QTL data", ]$Score))) 
df_summary<-df_summary[,c("Clinical data","GWAS data", "RNA-seq data", "QTL data" )]
colnames(df_summary) <- c( "Clinical data","GWAS data", "RNA-seq data", "QTL data")
corr_summary <- cor(df_summary, use = "pairwise.complete.obs", method = "spearman")

dev.off()
par(oma = c(9, 9, 9, 9), mar = c(0, 0, 0, 0))

# 绘制第一个 corrplot (上半部分)
corrplot(cor_matrix_spearman, type = "upper", order = 'original',
         col = rev(colbar(200)), method = "color", 
         cl.pos = "n", outline = TRUE,
         tl.col = Group$Color, tl.cex = 1,
         diag = TRUE,
         tl.srt = 45, addCoef.col = "black")

# 在同一画布上绘制第二个 corrplot (下半部分)
par(new = TRUE)  # 在同一个图层上进行新的绘制
corrplot(corr_summary, type = "lower", order = 'original',
         col = rev(colbar(200)), method = "color", add = TRUE,
         cl.pos = "b", 
         tl.col = coltext, tl.cex = 1,
         tl.srt = 45, outline = TRUE, diag = TRUE,
         addCoef.col = "black")

# 绘制条形图
# 计算每列非零值的个数
non_zero_counts <- apply(data_subset, 2, function(x) sum(x != 0))
non_zero_counts <- rev(non_zero_counts)
color_palette <- colorRampPalette(c("#FFFFFF", "#c74546"))
color_values <- color_palette(100)  # 创建100个渐变颜色
value_bins <- cut(non_zero_counts, breaks = 100, labels = FALSE)  # 将数值分为100个区间
bar_colors <- color_values[value_bins]  # 为每个条形分配
par(new = TRUE)
par(fig = c(0.752, 0.95, 0.084, 0.605), mar = c(0,0, 0, 0))  # 调整 fig 参数以控制条形图的位置和大小
# 绘制条形图，使用颜色渐变并应用所需的调整
bar_positions <- barplot(non_zero_counts, horiz = TRUE, las = 1, col = bar_colors,
                         xlab = "Count", ylab = "",  names.arg = rep("", length(non_zero_counts)),
                         main = "", axes = FALSE, xlim = c(0, 1000), cex.lab = 1, font.lab = 2)

# 手动添加 X 轴和 Y 轴的刻度
axis(1, cex.axis =1, font.axis = 2, tck = -0.015, lwd = 2) # 缩短 X 轴刻度线长度，tck 控制刻度线长度

# 在条形图右侧显示每个条形的具体数值
text(x = non_zero_counts, y = bar_positions, labels = non_zero_counts, pos = 4, cex = 1, font = 2)
#20.5*16.7





