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
library("GEOquery")
# +================================================+ ####
# +====Section 1. GWAS data arrangement============+ ####
# +================================================+ ####
# >>>>> Section 1.1. GWAS data arrangement <<<<< ####
#### Section 1.1.1 Merge SNP information ####
Periodontal_orginal<-fread("0. Source/finngen_R10_K11_GINGIVITIS_PERIODONTAL.gz")
colnames(Periodontal_orginal)
Periodontal_orginal<-setorder(Periodontal_orginal,pval)
Periodontal_orginal$rsids<-gsub("\\,.*", "",Periodontal_orginal$rsids)
Periodontal_orginal$N<-(97830+272252)
Periodontal_orginal%>%
  dplyr::select(SNP=rsids,CHR=`#chrom`,BP=pos,A1=alt,A2=ref,FRQ=af_alt,
                BETA=beta,SE=sebeta,P=pval,N=N,
                gene=nearest_genes) ->Periodontal_orginal_MungeSumstats
Periodontal_MungeSumstats<-format_sumstats(Periodontal_orginal_MungeSumstats,
                                           ref_genome = "GRCh38",
                                           bi_allelic_filter =F,
                                           allele_flip_frq=F,
                                           return_data =T)
Periodontal_sumstats <- liftover(sumstats_dt=Periodontal_MungeSumstats, 
                                 ref_genome = "hg38",
                                 convert_ref_genome="hg19")
Periodontal_sumstats%>%
  dplyr::select(SNP=SNP,chr=CHR,pos=BP,effect_allele=A1,other_allele=A2,
                eaf=FRQ,beta=BETA,se=SE,pval=P,samplesize=N,
                gene=GENE) ->Periodontal_orginal_data
Periodontal_orginal_data$Phenotype<-"Periodontal disease"
Periodontal_orginal_data$ncase<-97830
Periodontal_orginal_data$ncontrol<-272252

#### Section 1.1.2 Omit MHC SNP ####
Periodontal_orginal_data$chr<-as.character(Periodontal_orginal_data$chr)
Periodontal_orginal_data$chr<-as.numeric(Periodontal_orginal_data$chr)
Periodontal_mhc_data<-subset(Periodontal_orginal_data,chr!=6|(chr==6&(pos>33448354|pos<28477797)))
Periodontal_orginal_data<-Periodontal_mhc_data[Periodontal_mhc_data$chr %in% c(1:22), ]

#### Section 1.1.3 Omit duplicated SNP ####
Periodontal_orginal_data_sorted <- Periodontal_orginal_data[order(Periodontal_orginal_data$pval), ]
Periodontal_orginal_data_unique <- Periodontal_orginal_data_sorted[!duplicated(Periodontal_orginal_data_sorted$SNP), ]

#### Section 1.1.4 Save SNP ####
write.table(Periodontal_orginal_data_unique, file ="1. Original/Periodontal_orginal_data.csv", sep ="," ,row.names =F, col.names =T)
Periodontal_orginal_data_exposure<-subset(Periodontal_orginal_data,pval<1e-5)
write.table(Periodontal_orginal_data_exposure, file ="3. MR/Periodontal_orginal_data_exposure.csv", sep ="," ,row.names =F, col.names =T)
# +================================================+ ####
# +====Section 2. Gene data arrangement============+ ####
# +================================================+ #### 
# >>>>> Section 1.2. Blood Gene data arrangement <<<<< ####
#### Section 1.2.1 Download information ####
setwd("D:/PD_MDD/Step 0. Data/1. Original")
PD_set <- getGEO(GEO = "GSE156993",AnnotGPL=TRUE,destdir = ".")
PD_genset<- PD_set[[1]]##ExpressionSet
PD_gene_feature <- PD_genset@featureData@data
PD_gene_pheno <- PD_genset@phenoData@data

#### Section 1.2.2 Clearn data ####
table(PD_gene_pheno$source_name_ch1)
PD_pheno<-subset(PD_gene_pheno,source_name_ch1=="Healthy subjects, PBMC"|source_name_ch1=="Periodontitis, PBMC")
rownames(PD_pheno)
write.table (PD_pheno, file ="PD_blood_pheno.txt", sep ="\t", row.names =F , quote =F)
PD_gene <-as.data.frame(exprs(PD_genset)) 
PD_gene$ID<-rownames(PD_gene)
PD_gene_feature$`Gene symbol`<- gsub("///.*", "", PD_gene_feature$`Gene symbol`)
PD_gene_feature$`Gene ID`<- gsub("///.*", "", PD_gene_feature$`Gene ID`)
PD_gene_feature$Chromosome <- gsub(".*Chromosome ([0-9]+),.*", "\\1", PD_gene_feature$`Chromosome annotation`)

PD_feature<-na.omit(PD_gene_feature[,c("ID","Gene symbol","Gene ID","Chromosome")])
PD_feature<-subset(PD_feature,`Gene symbol`!="")
PD_exp<-merge(PD_gene,PD_feature,by="ID",all = F)
PD_exp<-PD_exp[,c("Gene symbol","Gene ID","Chromosome",rownames(PD_pheno))]
PD_exp<-PD_exp[PD_exp$Chromosome%in% c(1:22), ]
#### Section 1.2.4 Expression result ####
PD_exp_result <- PD_exp %>%
  group_by(`Gene symbol`) %>%
  summarise(`Gene ID` =  head(`Gene ID`, 1),Chromosome =  head(Chromosome, 1), 
    across(starts_with("GSM"), \(x) mean(x, na.rm = TRUE)),.groups = 'drop')
PD_exp_result<- as.data.frame(PD_exp_result)
PD_gene<-PD_exp_result[,c(1:3)]
row.names(PD_exp_result)<-PD_exp_result$`Gene symbol`
PD_exp<-PD_exp_result[,-c(1:3)]

write.table (PD_exp, file ="PD_blood_exp_result.txt", sep ="\t", row.names =T, col.names =T, quote =F)
write.table (PD_gene, file ="PD_blood_gene.txt", sep ="\t", row.names =F, col.names =T, quote =F)
# >>>>> Section 1.2. Sample Gene data arrangement <<<<< ####

PD_set <- getGEO(GEO = "GSE10334",AnnotGPL=TRUE,destdir = ".")
PD_genset<- PD_set[[1]]##ExpressionSet
PD_gene_feature <- PD_genset@featureData@data
PD_gene_pheno <- PD_genset@phenoData@data
table(PD_gene_pheno$characteristics_ch1)
write.table (PD_gene_pheno, file ="PD_tissues_pheno.txt", sep ="\t", row.names =F , quote =F)
PD_gene <-as.data.frame(exprs(PD_genset)) 
PD_gene$ID<-rownames(PD_gene)
PD_gene_feature$`Gene symbol`<- gsub("///.*", "", PD_gene_feature$`Gene symbol`)
PD_gene_feature$`Gene ID`<- gsub("///.*", "", PD_gene_feature$`Gene ID`)
PD_gene_feature$Chromosome <- gsub(".*Chromosome ([0-9]+),.*", "\\1", PD_gene_feature$`Chromosome annotation`)

PD_feature<-na.omit(PD_gene_feature[,c("ID","Gene symbol","Gene ID","Chromosome")])
PD_feature<-subset(PD_feature,`Gene symbol`!="")
PD_exp<-merge(PD_gene,PD_feature,by="ID",all = F)
PD_exp<-PD_exp[,c("Gene symbol","Gene ID","Chromosome",rownames(PD_gene_pheno))]
PD_exp<-PD_exp[PD_exp$Chromosome%in% c(1:22), ]
#### Section 1.2.4 Expression result ####
PD_exp_result <- PD_exp %>%
  group_by(`Gene symbol`) %>%
  summarise(`Gene ID` =  head(`Gene ID`, 1),Chromosome =  head(Chromosome, 1), 
            across(starts_with("GSM"), \(x) mean(x, na.rm = TRUE)),.groups = 'drop')
PD_exp_result<- as.data.frame(PD_exp_result)
PD_gene<-PD_exp_result[,c(1:3)]
row.names(PD_exp_result)<-PD_exp_result$`Gene symbol`
PD_exp<-PD_exp_result[,-c(1:3)]

write.table (PD_exp, file ="PD_tissues_exp_result.txt", sep ="\t", row.names =T, col.names =T, quote =F)
write.table (PD_gene, file ="PD_tissues_gene.txt", sep ="\t", row.names =F, col.names =T, quote =F)
# +================================================+ ####
# +====Section 3. MR data arrangement==============+ ####
# +================================================+ #### 
#dir.create("2. MR")
# >>>>> Section 2.1. Outcome data <<<<< ####
Periodontal_outcome<-read_outcome_data(
  filename="1. Original/Periodontal_orginal_data.csv",
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
save(Periodontal_outcome, file ="3. MR/Periodontal_outcome.Rdata")

# >>>>> Section 2.2. Exposure data <<<<< ####
Periodontal_exposure_orginal<-read_exposure_data(
  filename="3. MR/Periodontal_orginal_data_exposure.csv",
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
write.table(Periodontal_exposure_orginal, file ="3. MR/Periodontal_exposure_orginal.csv", sep ="," ,row.names =F, col.names =T)
Periodontal_exposure_orginal<-fread(file ="3. MR/Periodontal_exposure_orginal.csv")
Periodontal_exposure_orginal%>%
  dplyr::select(rsid=SNP,pval=pval.exposure,id=id.exposure) %>%
  filter(rsid!=".")->Periodontal_exposure_1E6_clump
#Periodontal<-clump_data(Periodontal_exposure_orginal,pop = "EUR")
Rownames_local<- ld_clump_local(dat = Periodontal_exposure_1E6_clump,
                                clump_kb = 10000,
                                clump_r2 = 0.001,
                                clump_p=1e-6,
                                bfile = "M:/MR_database/1kg.v3/EUR",
                                plink_bin = "M:/MR_database/plink/windows/plink.exe")
#select clumped data
Periodontal_exposure<-filter(Periodontal_exposure_orginal,SNP %in% Rownames_local$rsid)
write.table(Periodontal_exposure, file ="3. MR/Periodontal_exposure_covariate.csv", sep ="," ,row.names =F, col.names =T)
#Final exposure data
Periodontal_exposure<-read.table( file ="3. MR/Periodontal_exposure_covariate.csv", sep =",",header = T)
F_static<- (Periodontal_exposure$beta.exposure)^2/(Periodontal_exposure$se.exposure)^2
F_static

Periodontal_exposure<- Periodontal_exposure[-which(Periodontal_exposure$SNP=="rs7613444"|
                                                     Periodontal_exposure$SNP=="rs438811"),]
save(Periodontal_exposure, file ="3. MR/Periodontal_exposure.Rdata")

# +================================================+ ####
# +====Section 4. CAUSE data arrangement===========+ ####
# +================================================+ ####  
PD_orginal_data<-fread("1. Original/Periodontal_orginal_data.csv")
PD_orginal_data%>%
  dplyr::select(snp=SNP,beta_hat=beta,se=se,A1=effect_allele,A2=other_allele,P_value=pval) ->PD_summary
write.table (PD_summary, file ="4. CAUSE/PD_summary.txt", sep ="\t", row.names =F , quote =F)

# +================================================+ ####
# +====Section 5. GSMR data arrangement============+ ####
# +================================================+ ####  
#PD_orginal_data<-fread("1. Original/Periodontal_orginal_data.csv")
head(PD_orginal_data)
PD_orginal_data%>%
  dplyr::select(SNP=SNP,a1=effect_allele,a2=other_allele,a1_freq=eaf,
                bzx=beta,bzx_se=se,bzx_pval=pval,bzx_n=ncontrol) ->PD_summary
write.table (PD_summary, file ="5. GSMR/PD_summary.txt", sep ="\t", row.names =F , quote =F)

# +================================================+ ####
# +====Section 6. LDSC data arrangement============+ ####
# +================================================+ ####  
#PD_orginal_data<-fread("1. Original/Periodontal_orginal_data.csv")
PD_orginal_data%>%
  dplyr::select(SNP=SNP,P=pval,N=samplesize,A1=effect_allele,A2=other_allele,N_CAS=ncase,N_CON=ncontrol,
                BETA=beta,FRQ=eaf) ->PD_summary
write.table (PD_summary, file ="6. LDSC/PD_summary.txt", sep ="\t", row.names =F , quote =F)

# +================================================+ ####
# +====Section 7. HESS data arrangement============+ ####
# +================================================+ ####  
#PD_orginal_data<-fread("1. Original/Periodontal_orginal_data.csv")
PD_orginal_data$Z_score<-PD_orginal_data$beta/PD_orginal_data$se
PD_orginal_data%>%
  dplyr::select(SNP=SNP,CHR=chr,BP=pos,P=pval,N=samplesize,A1=effect_allele,A2=other_allele,Z=Z_score,N=samplesize) ->PD_summary

write.table (PD_summary, file ="7. HESS/PD_summary.txt", sep ="\t", row.names =F , quote =F)
#GWAS hits
PD_orginal_data%>%
  dplyr::select(CHR=chr,BP=pos) ->PD_summary
write.table (PD_summary, file ="7. HESS/PD_GWAS_hits.txt", sep ="\t", row.names =F , quote =F)

# +================================================+ ####
# +====Section 8. GNOVA data arrangement===========+ ####
# +================================================+ ####  
#PD_orginal_data<-fread("1. Original/Periodontal_orginal_data.csv")
head(PD_orginal_data)
PD_orginal_data%>%
  dplyr::select(SNP=SNP,A1=effect_allele,A2=other_allele,N=samplesize,P=pval,BETA=beta) ->PD_summary
write.table (PD_summary, file ="8. GNOVA/PD_summary.txt", sep ="\t", row.names =F , quote =F)

# +================================================+ ####
# +====Section 9. PLACO data arrangement==============+ ####
# +================================================+ ####  
#dir.create("5. PLACO")
#PD_orginal_data<-fread("1. Original/Periodontal_orginal_data.csv")
PD_orginal_data$Z_score<-PD_orginal_data$beta/PD_orginal_data$se
head(PD_orginal_data)
PD_orginal_data%>%
  dplyr::select(SNP=SNP,P=pval,Z=Z_score) ->PD_summary
write.table (PD_summary, file ="9. PLACO/PD_summary.txt", sep ="\t", row.names =F , quote =F)

# +================================================+ ####
# +====Section 10. MTAG data arrangement============+ ####
# +================================================+ #### 
#PD_orginal_data<-fread("1. Original/Periodontal_orginal_data.csv")
PD_orginal_data$Z_score<-PD_orginal_data$beta/PD_orginal_data$se
head(PD_orginal_data)
PD_orginal_data%>%
  dplyr::select(snpid=SNP,Chr=chr,bpos=pos,a1=effect_allele,
                a2=other_allele,freq=eaf,beta=beta,z=Z_score,pval=pval,n=samplesize) ->PD_summary
write.table (PD_summary, file ="10. MTAG/PD_summary.txt", sep ="\t", row.names =F , quote =F)

# +================================================+ ####
# +====Section 11. CPASSOC data arrangement========+ ####
# +================================================+ ####  
#PD_orginal_data<-fread("1. Original/Periodontal_orginal_data.csv")
PD_orginal_data$Z_score<-PD_orginal_data$beta/PD_orginal_data$se
head(PD_orginal_data)
PD_orginal_data%>%
  dplyr::select(SNP=SNP,Z=Z_score) ->PD_summary
write.table (PD_summary, file ="11. CPASSOC/PD_summary.txt", sep ="\t", row.names =F , quote =F)

# +================================================+ ####
# +====Section 16. SMR data arrangement=============+ ####
# +================================================+ ####  
#dir.create("10. SMR")
#PD_orginal_data<-fread("1. Original/Periodontal_orginal_data.csv")
head(PD_orginal_data)
PD_orginal_data%>%
  dplyr::select(SNP=SNP,A1=effect_allele,A2=other_allele,freq=eaf, b=beta, se=se, p=pval, n=samplesize) ->PD_summary
PD_summary<-na.omit(PD_summary)
write.table (PD_summary, file ="16. SMR/PD_summary.txt", sep ="\t", row.names =F , quote =F)

# +================================================+ ####
# +====Section 17. MR JTI data arrangement==============+ ####
# +================================================+ ####  
#dir.create("9. MR_JTI")
#PD_orginal_data<-fread("1. Original/Periodontal_orginal_data.csv")
head(PD_orginal_data)
PD_orginal_data%>%
  dplyr::select(rsid=SNP,beta=beta, se=se, eff_allele=effect_allele,ref_allele=other_allele) ->PD_summary
PD_summary<-na.omit(PD_summary)
head(PD_summary)
write.table(PD_summary, gzfile("17. MR_JTI/PD_summary.gz"), sep = "\t", row.names = FALSE, quote = FALSE)




