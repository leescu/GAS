setwd("D:/PD_MDD/MR")
library(data.table)
library(dplyr)
library(TwoSampleMR)
library(devtools)
library(tidyverse)
library(dplyr)
#devtools::install_github("mrcieu/ieugwasr")
library(ieugwasr)
library(MungeSumstats)
#remotes::install_github("MRCIEU/TwoSampleMR")
#devtools::install_local("./TwoSampleMR-master.zip")
Periodontal_orginal<-fread("data/finngen_R10_K11_GINGIVITIS_PERIODONTAL.gz")
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

#	Optional 
Periodontal_orginal_data$Phenotype<-"Periodontal disease"
Periodontal_orginal_data$ncase<-97830
Periodontal_orginal_data$ncontrol<-272252
Periodontal_orginal_data$chr<-as.character(Periodontal_orginal_data$chr)
Periodontal_orginal_data$chr<-as.numeric(Periodontal_orginal_data$chr)
Periodontal_orginal_data$chr<-as.factor(Periodontal_orginal_data$chr)
Periodontal_mhc_data<-subset(Periodontal_orginal_data,chr!=6|(chr==6&(pos>33448354|pos<28477797)))
Periodontal_mhc_data$chr<-as.character(Periodontal_mhc_data$chr)
Periodontal_mhc_data$chr<-as.factor(Periodontal_mhc_data$chr)
write.table(Periodontal_mhc_data, file ="data/Periodontal_orginal_data.csv", sep ="," ,row.names =F, col.names =T)
#Periodontal_orginal_data<-fread("./Periodontal_orginal_data.csv")
#Periodontal_orginal_data<-subset(Periodontal_orginal_data,eaf>0.05)
Periodontal_orginal_data_exposure<-subset(Periodontal_orginal_data,pval<1e-5)
write.table(Periodontal_orginal_data_exposure, file ="data/Periodontal_orginal_data_exposure.csv", sep ="," ,row.names =F, col.names =T)
#outcome_data
Periodontal_outcome<-read_outcome_data(
  filename="data/Periodontal_orginal_data.csv",
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
save(Periodontal_outcome, file ="data/Periodontal_outcome.Rdata")
#exposure_data
Periodontal_exposure_orginal<-read_exposure_data(
  filename="data/Periodontal_orginal_data_exposure.csv",
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
write.table(Periodontal_exposure_orginal, file ="data/Periodontal_exposure_orginal.csv", sep ="," ,row.names =F, col.names =T)
Periodontal_exposure_orginal<-fread(file ="data/Periodontal_exposure_orginal.csv")
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
write.table(Periodontal_exposure, file ="data/Periodontal_exposure_covariate.csv", sep ="," ,row.names =F, col.names =T)
#Final exposure data
Periodontal_exposure<-read.table( file ="data/Periodontal_exposure_covariate.csv", sep =",",header = T)
F_static<- (Periodontal_exposure$beta.exposure)^2/(Periodontal_exposure$se.exposure)^2
F_static

Periodontal_exposure<- Periodontal_exposure[-which(Periodontal_exposure$SNP=="rs7613444"|
                                                     Periodontal_exposure$SNP=="rs438811"),]
save(Periodontal_exposure, file ="data/Periodontal_exposure.Rdata")



