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
#devtools::install_github("mrcieu/ieugwasr")
library(ieugwasr)
library(MungeSumstats)
library(GEOquery)
library(tibble)
library(biomaRt)
library(org.Hs.eg.db)
library(clusterProfiler)
getOption('timeout')
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(timeout=10000)
# +================================================+ ####
# +====Section 1. GWAS data arrangement============+ ####
# +================================================+ ####
# >>>>> Section 1.1. GWAS data arrangement <<<<< ####
#### Section 1.1.1 Merge SNP information ####
Depressive_orginal<-fread("0. Source/PGC_UKB_depression_genome-wide.txt")
SNP=fread("M:/MR_database/snp150_hg19.txt",data.table =T)
colnames(SNP)[2]<-'MarkerName'
Depressive_SNP<-merge(Depressive_orginal, SNP, by="MarkerName",all.x = T)

Depressive_SNP_last<-separate(Depressive_SNP,`chromosome:start`, c("chromosome", "base_pair_location"), ":", extra = "merge")
Depressive_orginal<-setorder(Depressive_SNP_last,P)
Depressive_orginal$N<-(170756 +329443)
Depressive_orginal%>%
  dplyr::select(SNP=MarkerName,CHR=chromosome,BP=base_pair_location,
                A1=A1,A2=A2,FRQ=Freq,
                BETA=LogOR,SE=StdErrLogOR,P=P,N=N) ->Depressive_orginal_MungeSumstats
head(Depressive_orginal)
Depressive_orginal_MungeSumstats$A1<- toupper(Depressive_orginal_MungeSumstats$A1)
Depressive_orginal_MungeSumstats$A2<- toupper(Depressive_orginal_MungeSumstats$A2)
Depressive_orginal_MungeSumstats$CHR<-as.numeric(Depressive_orginal_MungeSumstats$CHR)
Depressive_orginal_MungeSumstats$BP<-as.numeric(Depressive_orginal_MungeSumstats$BP)
Depressive_MungeSumstats<-format_sumstats(Depressive_orginal_MungeSumstats,
                                          ref_genome = "GRCh37",
                                          bi_allelic_filter =F,
                                          allele_flip_frq=F,
                                          return_data =T)
Depressive_sumstats <-Depressive_MungeSumstats
Depressive_sumstats%>%
  dplyr::select(SNP=SNP,chr=CHR,pos=BP,effect_allele=A1,other_allele=A2,
                eaf=FRQ,beta=BETA,se=SE,pval=P,samplesize=N
  ) ->Depressive_orginal_data
Depressive_orginal_data$Phenotype<-"Depressive Disorder"
Depressive_orginal_data$ncase<-170756
Depressive_orginal_data$ncontrol<-329443
Depressive_orginal_data$chr<-as.character(Depressive_orginal_data$chr)
Depressive_orginal_data$chr<-as.numeric(Depressive_orginal_data$chr)
Depressive_orginal_data$chr<-as.factor(Depressive_orginal_data$chr)

#### Section 1.1.2 Omit MHC SNP ####
Depressive_mhc_data<-subset(Depressive_orginal_data,chr!=6|(chr==6&(pos>33448354|pos<28477797)))
Depressive_mhc_data$chr<-as.character(Depressive_mhc_data$chr)
Depressive_mhc_data$chr<-as.factor(Depressive_mhc_data$chr)
Depressive_orginal_data<-Depressive_mhc_data[Depressive_mhc_data$chr %in% c(1:22), ]

#### Section 1.1.3 Omit duplicated SNP ####
Depressive_orginal_data_sorted <- Depressive_orginal_data[order(Depressive_orginal_data$pval), ]
Depressive_orginal_data_unique <- Depressive_orginal_data_sorted[!duplicated(Depressive_orginal_data_sorted$SNP), ]

#### Section 1.1.4 Save SNP ####
write.table(Depressive_orginal_data_unique, file ="1. Original/Depressive_orginal_data.csv", sep ="," ,row.names =F, col.names =T)
Depressive_orginal_data_exposure<-subset(Depressive_orginal_data,pval<5e-8)
write.table(Depressive_orginal_data_exposure, file ="3. MR/Depressive_orginal_data_exposure.csv", sep ="," ,row.names =F, col.names =T)
# +================================================+ ####
# +====Section 2. Gene data arrangement============+ ####
# +================================================+ ####
# >>>>> Section 1.2. Blood data arrangement <<<<< ####
#### Section 1.2.1 Download information ####
setwd("D:/PD_MDD/Step 0. Data/1. Original")
options( 'download.file.method.GEOquery' = 'libcurl' )
MDD_set <- getGEO(GEO = "GSE38206",AnnotGPL=TRUE,destdir = ".")
MDD_genset<- MDD_set[[1]]##ExpressionSet
MDD_gene_feature <- MDD_genset@featureData@data
MDD_gene_pheno <- MDD_genset@phenoData@data
head(MDD_gene_feature)
ENST_genes <- MDD_gene_feature[grepl("ENST", MDD_gene_feature$GeneName), ]
ensembl <- useMart("ensembl")
ensembl_dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
gene_info <- getBM(attributes = c("ensembl_transcript_id", "external_gene_name", "chromosome_name", "entrezgene_id"),
                   filters = "ensembl_transcript_id",
                   values = ENST_genes$GeneName, 
                   mart = ensembl_dataset)

gene_info <- gene_info %>%
  mutate(
    external_gene_name = ifelse(is.na(external_gene_name), "", external_gene_name),
    entrezgene_id = ifelse(is.na(entrezgene_id), "", entrezgene_id))
gene_info<-subset(gene_info,external_gene_name!=""&entrezgene_id!="")

ENST_genes <- ENST_genes %>%
  left_join(gene_info, by = c("GeneName" = "ensembl_transcript_id"),relationship = "many-to-many")

ACC_genes<-subset(MDD_gene_feature,GB_ACC!="")

gene_info <- getBM(attributes = c("external_gene_name","chromosome_name",  "entrezgene_id"),
                   filters = "external_gene_name",
                   values = ACC_genes$GeneName,  # 使用 ACC_genes 中的 GeneName 列
                   mart = ensembl_dataset)

ACC_gene <- merge(ACC_genes, gene_info, by.x = "GeneName", by.y = "external_gene_name", all =F)
ACC_genes_unique <- ACC_gene %>%
  distinct(ID, .keep_all = TRUE)
ENST_gene<-ENST_genes[,c("ID","external_gene_name","chromosome_name", "entrezgene_id")]
ACC_gene<-ACC_genes_unique[,c("ID","GeneName","chromosome_name", "entrezgene_id")]
colnames(ENST_gene)<-c("ID","Gene symbol","Chromosome","Gene ID")
colnames(ACC_gene)<-c("ID","Gene symbol","Chromosome", "Gene ID")
MDD_feature<-na.omit(rbind(ENST_gene,ACC_gene))

#### Section 1.2.2 Clearn data ####
colnames(MDD_gene_pheno)
table(MDD_gene_pheno$`sample collection:ch1`)

MDD_pheno<-subset(MDD_gene_pheno,`sample collection:ch1`=='1st week (during a severe episode of depression)')
rownames(MDD_pheno)
write.table (MDD_pheno, file ="MDD_blood_pheno.txt", sep ="\t", row.names =F , quote =F)
MDD_gene <-as.data.frame(exprs(MDD_genset)) 
MDD_gene$ID<-rownames(MDD_gene)
MDD_exp<-merge(MDD_gene,MDD_feature,by="ID",all = F)
MDD_exp<-MDD_exp[,c("Gene symbol","Gene ID","Chromosome", rownames(MDD_pheno))]
MDD_exp<-MDD_exp[MDD_exp$Chromosome%in% c(1:22), ]
#### Section 1.2.4 Expression result ####
MDD_exp_result <- MDD_exp %>%
  group_by(`Gene symbol`) %>%
  summarise(`Gene ID` = head(`Gene ID`, 1), Chromosome = head(`Chromosome`, 1), 
            across(starts_with("GSM"), \(x) mean(x, na.rm = TRUE)),.groups = 'drop')
MDD_exp_result<- as.data.frame(MDD_exp_result)
MDD_gene<-MDD_exp_result[,c(1:3)]
row.names(MDD_exp_result)<-MDD_exp_result$`Gene symbol`
MDD_exp<-MDD_exp_result[,-c(1:3)]

write.table (MDD_exp, file ="MDD_blood_exp_result.txt", sep ="\t", row.names =T, col.names =T, quote =F)
write.table (MDD_gene, file ="MDD_blood_gene.txt", sep ="\t", row.names =F, col.names =T, quote =F)
# >>>>> Section 1.2. tissues data arrangement <<<<< ####
MDD_set <- getGEO(GEO = "GSE53987",AnnotGPL=TRUE,destdir = ".")
MDD_genset<- MDD_set[[1]]##ExpressionSet
MDD_gene_feature <- MDD_genset@featureData@data
MDD_gene_pheno <- MDD_genset@phenoData@data
table(MDD_gene_pheno$characteristics_ch1.7)
MDD_gene_pheno<-subset(MDD_gene_pheno,characteristics_ch1.7=="disease state: control"|characteristics_ch1.7=="disease state: major depressive disorder")
write.table (MDD_gene_pheno, file ="MDD_tissues_pheno.txt", sep ="\t", row.names =F , quote =F)
MDD_gene <-as.data.frame(exprs(MDD_genset)) 
MDD_gene$ID<-rownames(MDD_gene)
MDD_gene_feature$`Gene symbol`<- gsub("///.*", "", MDD_gene_feature$`Gene symbol`)
MDD_gene_feature$`Gene ID`<- gsub("///.*", "", MDD_gene_feature$`Gene ID`)
MDD_gene_feature$Chromosome <- gsub(".*Chromosome ([0-9]+),.*", "\\1", MDD_gene_feature$`Chromosome annotation`)

MDD_feature<-na.omit(MDD_gene_feature[,c("ID","Gene symbol","Gene ID","Chromosome")])
MDD_feature<-subset(MDD_feature,`Gene symbol`!="")
MDD_exp<-merge(MDD_gene,MDD_feature,by="ID",all = F)
MDD_exp<-MDD_exp[,c("Gene symbol","Gene ID","Chromosome",rownames(MDD_gene_pheno))]
MDD_exp<-MDD_exp[MDD_exp$Chromosome%in% c(1:22), ]
#### Section 1.2.4 Expression result ####
MDD_exp_result <- MDD_exp %>%
  group_by(`Gene symbol`) %>%
  summarise(`Gene ID` =  head(`Gene ID`, 1),Chromosome =  head(Chromosome, 1), 
            across(starts_with("GSM"), \(x) mean(x, na.rm = TRUE)),.groups = 'drop')
MDD_exp_result<- as.data.frame(MDD_exp_result)
MDD_gene<-MDD_exp_result[,c(1:3)]
row.names(MDD_exp_result)<-MDD_exp_result$`Gene symbol`
MDD_exp<-MDD_exp_result[,-c(1:3)]

write.table (MDD_exp, file ="MDD_tissues_exp_result.txt", sep ="\t", row.names =T, col.names =T, quote =F)
write.table (MDD_gene, file ="MDD_tissues_gene.txt", sep ="\t", row.names =F, col.names =T, quote =F)

# +================================================+ ####
# +====Section 3. MR data arrangement==============+ ####
# +================================================+ ####  
# >>>>> Section 2.1. Outcome data <<<<< ####
Depressive_outcome<-read_outcome_data(
  filename="1. Original/Depressive_orginal_data.csv",
  snps = NULL,
  sep = ",",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  samplesize_col = "samplesize",
  gene_col = "gene",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE,
  chr_col = "chr",
  pos_col = "pos"
)
save(Depressive_outcome, file ="3. MR/Depressive_outcome.Rdata")
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++ ####
# >>>>> Section 2.2. Exposure data <<<<< ####
Depressive_exposure_orginal<-read_exposure_data(
  filename="3. MR/Depressive_orginal_data_exposure.csv",
  clump = FALSE,
  sep = ",",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  samplesize_col = "samplesize",
  gene_col = "gene",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE,
  chr_col = "chr",
  pos_col = "pos"
)
write.table(Depressive_exposure_orginal, file ="3. MR/Depressive_exposure_orginal.csv", sep ="," ,row.names =F, col.names =T)

Depressive_exposure_orginal<-fread("3. MR/Depressive_exposure_orginal.csv")
Depressive_exposure_orginal%>%
  dplyr::select(rsid=SNP,pval=pval.exposure,id=id.exposure) %>%
  filter(rsid!=".")->Depressive_exposure_5E8_clump
#Depressive<-clump_data(Depressive_exposure_orginal,pop = "EUR")
Rownames_local<- ld_clump_local(dat = Depressive_exposure_5E8_clump,
                                clump_kb = 10000,
                                clump_r2 = 0.001,
                                clump_p=5e-8,
                                bfile = "M:/MR_database/1kg.v3/EUR",
                                plink_bin = "M:/MR_database/plink/windows/plink.exe")
#select clumped data
Depressive_exposure<-filter(Depressive_exposure_orginal,SNP %in% Rownames_local$rsid)
write.table(Depressive_exposure, file ="3. MR/Depressive_exposure_covariate.csv", sep ="," ,row.names =F, col.names =T)
#Final exposure data
Depressive_exposure<-read.table( file ="3. MR/Depressive_exposure_covariate.csv", sep =",",header = T)
Depressive_exposure<-Depressive_exposure[-which(Depressive_exposure$SNP=="rs7551758"|
                                                  Depressive_exposure$SNP=="rs2568958"|
                                                  Depressive_exposure$SNP=="rs2232423"|
                                                  Depressive_exposure$SNP=="rs2214123"|
                                                  Depressive_exposure$SNP=="rs198457"|
                                                  Depressive_exposure$SNP=="rs1367635"|
                                                  Depressive_exposure$SNP=="rs12967143"),]

save(Depressive_exposure, file ="3. MR/Depressive_exposure.Rdata")
# +================================================+ ####
# +====Section 4. CAUSE data arrangement===========+ ####
# +================================================+ ####  
MDD_orginal_data<-fread("1. Original/Depressive_orginal_data.csv")
MDD_orginal_data%>%
  dplyr::select(snp=SNP,beta_hat=beta,se=se,A1=effect_allele,A2=other_allele,P_value=pval) ->MDD_summary
write.table (MDD_summary, file ="4. CAUSE/MDD_summary.txt", sep ="\t", row.names =F , quote =F)

# +================================================+ ####
# +====Section 5. GMSR data arrangement============+ ####
# +================================================+ ####  
#MDD_orginal_data<-fread("1. Original/Depressive_orginal_data.csv")
head(MDD_orginal_data)
MDD_orginal_data%>%
  dplyr::select(SNP=SNP,
                bzy=beta,bzy_se=se,bzy_pval=pval,bzy_n=ncontrol) ->MDD_summary
write.table (MDD_summary, file ="5. GSMR/MDD_summary.txt", sep ="\t", row.names =F , quote =F)

# +================================================+ ####
# +====Section 6. LDSC data arrangement============+ ####
# +================================================+ ####  
#MDD_orginal_data<-fread("1. Original/Depressive_orginal_data.csv")
MDD_orginal_data%>%
  dplyr::select(SNP=SNP,P=pval,N=samplesize,A1=effect_allele,A2=other_allele,N_CAS=ncase,N_CON=ncontrol,
                BETA=beta,FRQ=eaf) ->MDD_summary
write.table (MDD_summary, file ="6. LDSC/MDD_summary.txt", sep ="\t", row.names =F , quote =F)

# +================================================+ ####
# +====Section 7. HESS data arrangement============+ ####
# +================================================+ ####  
#MDD_orginal_data<-fread("1. Original/Depressive_orginal_data.csv")
MDD_orginal_data$Z_score<-MDD_orginal_data$beta/MDD_orginal_data$se

MDD_orginal_data%>%
  dplyr::select(SNP=SNP,CHR=chr,BP=pos,P=pval,N=samplesize,A1=effect_allele,A2=other_allele,Z=Z_score,N=samplesize) ->MDD_summary

write.table (MDD_summary, file ="7. HESS/MDD_summary.txt", sep ="\t", row.names =F , quote =F)
#GWAS hits
MDD_orginal_data%>%
  dplyr::select(CHR=chr,BP=pos) ->MDD_summary
write.table (MDD_summary, file ="7. HESS/MDD_GWAS_hits.txt", sep ="\t", row.names =F , quote =F)
# +================================================+ ####
# +====Section 8. GNOVA data arrangement===========+ ####
# +================================================+ ####  
#MDD_orginal_data<-fread("1. Original/Depressive_orginal_data.csv")
head(MDD_orginal_data)
MDD_orginal_data%>%
  dplyr::select(SNP=SNP,A1=effect_allele,A2=other_allele,N=samplesize,P=pval,BETA=beta) ->MDD_summary
write.table (MDD_summary, file ="8. GNOVA/MDD_summary.txt", sep ="\t", row.names =F , quote =F)

# +================================================+ ####
# +====Section 9. PLACO data arrangement===========+ ####
# +================================================+ ####  
#dir.create("5. PLACO")
#MDD_orginal_data<-fread("1. Original/Depressive_orginal_data.csv")
MDD_orginal_data$Z_score<-MDD_orginal_data$beta/MDD_orginal_data$se
head(MDD_orginal_data)
MDD_orginal_data%>%
  dplyr::select(SNP=SNP,P=pval,Z=Z_score) ->MDD_summary
write.table (MDD_summary, file ="9. PLACO/MDD_summary.txt", sep ="\t", row.names =F , quote =F)

# +================================================+ ####
# +====Section 10. MTAG data arrangement===========+ ####
# +================================================+ ####  
#MDD_orginal_data<-fread("1. Original/Depressive_orginal_data.csv")
MDD_orginal_data$Z_score<-MDD_orginal_data$beta/MDD_orginal_data$se
head(MDD_orginal_data)
MDD_orginal_data%>%
  dplyr::select(snpid=SNP,Chr=chr,bpos=pos,a1=effect_allele,
                a2=other_allele,freq=eaf,beta=beta,z=Z_score,pval=pval,n=samplesize) ->MDD_summary
write.table (MDD_summary, file ="10. MTAG/MDD_summary.txt", sep ="\t", row.names =F , quote =F)

# +================================================+ ####
# +====Section 11. CPASSOC data arrangement========+ ####
# +================================================+ ####  
#MDD_orginal_data<-fread("1. Original/Depressive_orginal_data.csv")
MDD_orginal_data$Z_score<-MDD_orginal_data$beta/MDD_orginal_data$se
head(MDD_orginal_data)
MDD_orginal_data%>%
  dplyr::select(SNP=SNP,Z=Z_score) ->MDD_summary
write.table (MDD_summary, file ="11. CPASSOC/MDD_summary.txt", sep ="\t", row.names =F , quote =F)

# +================================================+ ####
# +====Section 16. SMR data arrangement=============+ ####
# +================================================+ ####  
#dir.create("10. SMR")
#MDD_orginal_data<-fread("1. Original/Depressive_orginal_data.csv")
MDD_orginal_data$N<-500199
MDD_orginal_data%>%
  dplyr::select(SNP=SNP,A1=effect_allele,A2=other_allele,freq=eaf, b=beta, se=se, p=pval, n=N) ->MDD_summary
MDD_summary<-na.omit(MDD_summary)
write.table (MDD_summary, file ="16. SMR/MDD_summary.txt", sep ="\t", row.names =F , quote =F)

# +================================================+ ####
# +====Section 17. MR JTI data arrangement==============+ ####
# +================================================+ ####  
dir.create("17. MR_JTI")
MDD_orginal_data<-fread("1. Original/Depressive_orginal_data.csv")
head(MDD_orginal_data)
MDD_orginal_data%>%
  dplyr::select(rsid=SNP,beta=beta, se=se, eff_allele=effect_allele,ref_allele=other_allele) ->MDD_summary

head(MDD_summary)
write.table(MDD_summary, gzfile("17. MR_JTI/MDD_summary.gz"),sep = "\t", row.names = FALSE, quote = FALSE)



