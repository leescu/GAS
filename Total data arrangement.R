# +================================================+ ####
# +====Section 0. Packages and setwd===============+ ####
# +================================================+ #### 
setwd("D:/PD_MDD/step 0. Data")
library(data.table)
library(TwoSampleMR)
library(devtools)
library(tidyverse)
library(dplyr)
library(tidyr)
#devtools::install_github("mrcieu/ieugwasr")
library(ieugwasr)
library(MungeSumstats)
library(mice)
library(biomaRt)
library(impute)
# +================================================+ ####
# +====Section 1. RAW data arrangement============+ ####
# +================================================+ ####
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++ ####
# >>>>> section 1 Cross-sectional MDD arrangement ####
# ** PHQ-9 ####
Depression<-fread("1. Original/Depression_UKB.csv")
colnames(Depression)
PHQ<-Depression[,c(1,13:21)]
table(PHQ$p120104)
PHQ[PHQ==''|PHQ=='Prefer not to answer']<-NA
PHQ[PHQ=='Not at all']<-0
PHQ[PHQ=='Several days']<-1
PHQ[PHQ=='More than half the days']<-2
PHQ[PHQ=='Nearly every day']<-3
PHQ<-na.omit(PHQ)
rownames(PHQ)<-PHQ$eid
PHQ$eid<-NULL

colnames(PHQ)
PHQ$p120104<-as.numeric(PHQ$p120104)
PHQ$p120105<-as.numeric(PHQ$p120105)
PHQ$p120106<-as.numeric(PHQ$p120106)
PHQ$p120107<-as.numeric(PHQ$p120107)
PHQ$p120108<-as.numeric(PHQ$p120108)
PHQ$p120109<-as.numeric(PHQ$p120109)
PHQ$p120110<-as.numeric(PHQ$p120110)
PHQ$p120111<-as.numeric(PHQ$p120111)
PHQ$p120112<-as.numeric(PHQ$p120112)

PHQ$PHQ<-rowSums(PHQ,na.rm = F)

PHQ$MDD1 <- ifelse(PHQ$PHQ >= 10, 1, 0)
PHQ$eid<-rownames(PHQ)
PHQ$eid<-as.character(PHQ$eid)
table(PHQ$MDD1)
# ** Smith et al. (2013) PLOS ONE. 8: e75362 ####
Smith<-Depression[,c("eid","p20126_i0")]

Smith[Smith==''|Smith=='Prefer not to answer']<-NA
table(Smith$p20126_i0)
Smith$MDD2[Smith$p20126_i0=="Probable Recurrent major depression (moderate)"|
             Smith$p20126_i0=="Probable Recurrent major depression (severe)"]<-1
Smith$MDD2[Smith$p20126_i0=="Bipolar I Disorder"|
             Smith$p20126_i0=="Bipolar II Disorder"|
             Smith$p20126_i0=="Single Probable major depression episode"|
             Smith$p20126_i0=="No Bipolar or Depression"]<-0
table(Smith$MDD2)
Smith$p20126_i0<-NULL
Smith$eid<-as.character(Smith$eid)
MDD_CS<-merge(PHQ,Smith,by = "eid",all = T)
MDD_CS$MDD[MDD_CS$MDD1==0|MDD_CS$MDD2==0]<-0
MDD_CS$MDD[MDD_CS$MDD1==1|MDD_CS$MDD2==1]<-1
table(MDD_CS$MDD)
colnames(MDD_CS)
MDD_CS<-MDD_CS[,c("eid","MDD")]
save(MDD_CS,file = "1. Original/MDD_CS.Rdata")
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++ ####
# >>>>> section 2 Cohort MDD arrangement ####
colnames(Depression)
Data_participant<-fread("1. Original/Oral_UKB.csv")
Data_death<-fread("1. Original/Death_UKB.csv")
Data_death$dnx_death_id<-NULL
colnames(Data_death)[2]<-"Death_date"
Data_death$Death_date<-as.Date(Data_death$Death_date)
Data_death<-Data_death[order(Data_death$Death_date),]

Data_death <- unique(Data_death)
Lost<-Data_participant[,c("eid","p191","p53_i0")]
ICD<-Depression[,c("eid","p130894","p130895","p130896","p130897")]
MDD_ICD<-merge(ICD,Lost,by = "eid",all =T)
MDD_ICD<-merge(MDD_ICD,Data_death,by = "eid",all =T)
MDD_ICD[MDD_ICD==""] <- NA
MDD_ICD$p130894[MDD_ICD$p130894=="Code has event date matching participant's date of birth"] <- NA
MDD_ICD$p130896[MDD_ICD$p130896=="Code has event date matching participant's date of birth"] <- NA

MDD_ICD$Depression_F32<-ifelse(is.na(MDD_ICD$p130894),0,1)
MDD_ICD$Depression_F33<-ifelse(is.na(MDD_ICD$p130896),0,1)
MDD_ICD$Depression_Lost<-ifelse(is.na(MDD_ICD$p191),0,1)
MDD_ICD$Depression_Death<-ifelse(is.na(MDD_ICD$date_of_death),0,1)
MDD_ICD$F32_date<-as.Date(MDD_ICD$p130894)
MDD_ICD$F33_date<-as.Date(MDD_ICD$p130896)
MDD_ICD$Lost_date<-as.Date(MDD_ICD$p191)
MDD_ICD$Last_date<-as.Date("2023-09-30")
MDD_ICD$Death_date<-as.Date(MDD_ICD$Death_date)
MDD_ICD$Final_date<-MDD_ICD%>%
  transmute(earliest_date = pmin(MDD_ICD$F32_date,
                                 MDD_ICD$F33_date,
                                 MDD_ICD$Lost_date,
                                 MDD_ICD$Death_date,
                                 MDD_ICD$Last_date,
                                 na.rm = T))
MDD_ICD$First_date<-as.Date(MDD_ICD$p53_i0)
MDD_ICD$days<-MDD_ICD$Final_date-MDD_ICD$First_date
MDD_ICD<-subset(MDD_ICD,days>0)
MDD_ICD$MDD[MDD_ICD$Depression_F32==0&MDD_ICD$Depression_F33==0]<-0
MDD_ICD$MDD[MDD_ICD$Depression_F32==1|MDD_ICD$Depression_F33==1]<-1
table(MDD_ICD$MDD)
colnames(MDD_ICD)
MDD_ICD<-MDD_ICD[,c("eid","MDD","days")]
save(MDD_ICD,file = "1. Original/MDD_ICD.Rdata")
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++ ####
# >>>>> section 3 Oral data arrangement ####
Data_participant<-fread("1. Original/Oral_UKB.csv")
Oral_data<-as.data.frame(Data_participant$eid)
colnames(Oral_data)<-"eid"
Oral_data$oral<-Data_participant$p6149_i0
Oral_data<-subset(Oral_data,oral!=""&oral!=" Prefer not to answer"&oral!="Prefer not to answer")
table(Oral_data$oral,useNA='ifany')
Oral_data$Mouth_ulcers<-grepl("*Mouth ulcers*",Oral_data$oral)
Oral_data$Painful_gums<-grepl("*Painful gums*",Oral_data$oral)
Oral_data$Bleeding_gums<-grepl("*Bleeding gums*",Oral_data$oral)
Oral_data$Loose_teeth<-grepl("*Loose teeth*",Oral_data$oral)
Oral_data$Toothache<-grepl("*Toothache*",Oral_data$oral)
Oral_data$Dentures<-grepl("*Dentures*",Oral_data$oral)
Oral_data[Oral_data==T]<-1
Oral_data[Oral_data==F]<-0
Oral_data$Periodontal_disease<-0
Oral_data$Periodontal_disease[Oral_data$Painful_gums==1|Oral_data$Bleeding_gums==1|Oral_data$Loose_teeth]<-1
Oral_data$Periodontal_disease<-as.factor(Oral_data$Periodontal_disease)
head(Oral_data)
Oral_data_final<-Oral_data[,c("eid","Periodontal_disease")]
table(Oral_data$Periodontal_disease)
save(Oral_data_final,file = "1. Original/Oral_data_final.Rdata")
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++ ####
# >>>>> section 4 Protein data  ####
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
Protein_data_filled <- impute.knn(Protein_matrix)$data
# 检查插补后的数据是否还有缺失值
sum(is.na(Protein_data_filled))  
# 提取插补后的数据集（第一个插补结果）
Protein_data_filled<-as.data.frame(Protein_data_filled)

ensembl <- useMart("ensembl")
gene_values<-rownames(Protein_data_filled)
ensembl_dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
gene_info <- getBM(attributes = c("external_gene_name","chromosome_name",  "entrezgene_id"),
                   filters = "external_gene_name",
                   values = rownames(Protein_data_filled), 
                   mart = ensembl_dataset)
gene_information<-gene_info[gene_info$chromosome_name%in% c(1:22), ]

gene_information<- gene_information[!duplicated(gene_information$external_gene_name), ]
colnames(gene_information)<-c("Gene_symbol","Chromosome","Gene_ID")
colnames(gene_information)
Protein_data_filled$Gene_symbol<-rownames(Protein_data_filled)
Protein_final<-merge(gene_information,Protein_data_filled,by="Gene_symbol",all =F)
save(Protein_final,file = "1. Original/Protein_final.Rdata")
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++ ####
# >>>>> section 4 Merge data  ####
load(file="1. Original/MDD_ICD.Rdata")
load(file="1. Original/MDD_CS.Rdata")
load(file="1. Original/Covariate_data_interpolation.Rdata")
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
Interpolation_data_cohort<-merge(Data_final_cohort,Covariate_data_interpolation,by = "eid",all= F)
Interpolation_data_cross<-merge(Data_final_cross,Covariate_data_interpolation,by = "eid",all= F)

save(Interpolation_data_cohort,file="15. Cohort/Interpolation_data_cohort.Rdata")
save(Interpolation_data_cross,file="16. Cross_sectional/Interpolation_data_cross.Rdata")
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++ ####
# >>>>> section 5 Merge Protein data  ####
Protein<-merge(Interpolation_data_cross,Protein_data,by="eid",all=F)
Protein$eid <- paste("eid", Protein$eid, sep = "")
Protein_PD<-Protein[!is.na(Protein$Periodontal_disease),]
Protein_MDD<-Protein[!is.na(Protein$MDD),]

MDD_pheno_protein<-Protein_MDD[,1:19]
MDD_exp_protein<-as.data.frame(t(Protein_MDD[,c(20:2797)]))
colnames(MDD_exp_protein)<-Protein_MDD$eid

PD_pheno_protein<-Protein_PD[,1:19]
PD_exp_protein<-as.data.frame(t(Protein_PD[,c(20:2797)]))
colnames(PD_exp_protein)<-Protein_PD$eid


write.table (PD_exp_protein, file ="1. Original/PD_exp_protein.txt", sep ="\t", row.names =T, col.names =T, quote =F)
write.table (PD_gene_protein, file ="1. Original/PD_gene_protein.txt", sep ="\t", row.names =F, col.names =T, quote =F)
write.table (PD_pheno_protein, file ="1. Original/PD_pheno_protein.txt", sep ="\t", row.names =F, col.names =T, quote =F)
write.table (MDD_exp_protein, file ="1. Original/MDD_exp_protein.txt", sep ="\t", row.names =T, col.names =T, quote =F)
write.table (MDD_gene_protein, file ="1. Original/MDD_gene_protein.txt", sep ="\t", row.names =F, col.names =T, quote =F)
write.table (MDD_pheno_protein, file ="1. Original/MDD_pheno_protein.txt", sep ="\t", row.names =F, col.names =T, quote =F)

# +================================================+ ####
# +====Section 6. FUMA data arrangement=============+ ####
# +================================================+ ####  
dir.create("6. FUMA")
# >>>>> section 6.1. P <5e-8  ####

PLACO_original<-fread("5. PLACO/PLACO_original.txt",data.table =T)
PLACO_subset<-subset(PLACO_original,p.placo<5e-8)
PLACO_subset$T.placo<-NULL
colnames(PLACO_subset)<-c("SNP","CHR","BP","P")
write.table(PLACO_subset, file = "5. PLACO/PLACO_sig_5e8.txt", sep = "\t", row.names = FALSE,quote=F)
write.table(PLACO_subset, file = "6. FUMA/PLACO_sig_5e8.txt", sep = "\t", row.names = FALSE,quote=F)


# >>>>> section 6.1. P <1e-6  ####
#PLACO_original<-fread("PLACO_original.txt")
PLACO_subset<-subset(PLACO_original,p.placo<1e-6)
PLACO_subset$T.placo<-NULL
colnames(PLACO_subset)<-c("SNP","CHR","BP","P")
write.table(PLACO_subset, file = "5. PLACO/PLACO_sig_1e6.txt", sep = "\t", row.names = FALSE,quote=F)
write.table(PLACO_subset, file = "6. FUMA/PLACO_sig_1e6.txt", sep = "\t", row.names = FALSE,quote=F)

# >>>>> section 6.1. P.adj <0.05  ####
#PLACO_original<-fread("PLACO_original.txt")
PLACO_original$p.adj<-p.adjust(PLACO_original$p.placo,  # P值列表
                               method ="BH"                       # FDR校正的方法
)
PLACO_subset<-subset(PLACO_original,p.adj<0.05)
PLACO_subset<- PLACO_subset[order(PLACO_subset$p.adj), ]
PLACO_subset$T.placo<-NULL
PLACO_subset$p.adj<-NULL
colnames(PLACO_subset)<-c("SNP","CHR","BP","P")
PLACO_subset<- PLACO_subset[PLACO_subset$CHR %in% c(1:22), ]
write.table(PLACO_subset, file = "5. PLACO/PLACO_sig_fdr.txt", sep = "\t", row.names = FALSE,quote=F)
write.table(PLACO_subset, file = "6. FUMA/PLACO_sig_fdr.txt", sep = "\t", row.names = FALSE,quote=F)

