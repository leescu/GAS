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
#remotes::install_github("MRCIEU/TwoSampleMR")
#devtools::install_local("./TwoSampleMR-master.zip")
#BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")
#BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")

#BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh37")
#BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")


Depressive_orginal<-fread("2. MR/PGC_UKB_depression_genome-wide.txt")
head_original<-head(Depressive_orginal)
colnames(head_original)
#merge SNP information
SNP=fread("M:/MR_database/snp150_hg19.txt",data.table =T)
colnames(SNP)[2]<-'MarkerName'
Depressive_SNP<-merge(Depressive_orginal, SNP, by="MarkerName",all.x = T)

dd<-head(Depressive_SNP)
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
#Depressive_sumstats <- liftover(sumstats_dt=Depressive_MungeSumstats, 
#                               ref_genome = "hg19",
#                               convert_ref_genome="hg38")

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
Depressive_mhc_data<-subset(Depressive_orginal_data,chr!=6|(chr==6&(pos>33448354|pos<28477797)))
Depressive_mhc_data$chr<-as.character(Depressive_mhc_data$chr)
Depressive_mhc_data$chr<-as.factor(Depressive_mhc_data$chr)

write.table(Depressive_mhc_data, file ="data/Depressive_orginal_data.csv", sep ="," ,row.names =F, col.names =T)
Depressive_orginal_data<-fread("data/Depressive_orginal_data.csv")
#Depressive_orginal_data<-subset(Depressive_orginal_data,eaf>0.01)

Depressive_orginal_data_exposure<-subset(Depressive_orginal_data,pval<5e-8)
write.table(Depressive_orginal_data_exposure, file ="data/Depressive_orginal_data_exposure.csv", sep ="," ,row.names =F, col.names =T)

#outcome_data
Depressive_outcome<-read_outcome_data(
  filename="data/Depressive_orginal_data.csv",
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
save(Depressive_outcome, file ="data/Depressive_outcome.Rdata")
#exposure_data
Depressive_exposure_orginal<-read_exposure_data(
  filename="data/Depressive_orginal_data_exposure.csv",
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
write.table(Depressive_exposure_orginal, file ="data/Depressive_exposure_orginal.csv", sep ="," ,row.names =F, col.names =T)
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
write.table(Depressive_exposure, file ="data/Depressive_exposure_covariate.csv", sep ="," ,row.names =F, col.names =T)
#Final exposure data
Depressive_exposure<-read.table( file ="2. MR/Depressive_exposure_covariate.csv", sep =",",header = T)
 Depressive_exposure<-Depressive_exposure[-which(Depressive_exposure$SNP=="rs7551758"|
                                                Depressive_exposure$SNP=="rs2568958"|
                                                Depressive_exposure$SNP=="rs2232423"|
                                                Depressive_exposure$SNP=="rs2214123"|
                                                Depressive_exposure$SNP=="rs198457"|
                                                Depressive_exposure$SNP=="rs1367635"|
                                                Depressive_exposure$SNP=="rs12967143"),]

save(Depressive_exposure, file ="2. MR/Depressive_exposure.Rdata")






