setwd("I:/paper_9_PD&MDD")
library(data.table)
library(dplyr)
library(TwoSampleMR)
library(devtools)
library(tidyverse)
library(dplyr)

library(ieugwasr)
library(MungeSumstats)
Cognitive_orginal<-fread("data/GWAS_CP_all.txt")
colnames(Cognitive_orginal)
aa<-head(Cognitive_orginal)
Cognitive_orginal
Cognitive_orginal<-setorder(Cognitive_orginal,Pval)

Cognitive_orginal%>%
  dplyr::select(SNP=MarkerName,CHR=CHR,BP=POS,A1=A1,A2=A2,FRQ=EAF,
                BETA=Beta,SE=SE,P=Pval) ->Cognitive_orginal_MungeSumstats


Cognitive_orginal_MungeSumstats$CHR<-as.numeric(Cognitive_orginal_MungeSumstats$CHR)
Cognitive_orginal_MungeSumstats$BP<-as.numeric(Cognitive_orginal_MungeSumstats$BP)

Cognitive_MungeSumstats<-format_sumstats(Cognitive_orginal_MungeSumstats,
                                   ref_genome = "GRCh37",
                                   bi_allelic_filter =F,
                                   allele_flip_frq=F,
                                   return_data =T)
Cognitive_sumstats <-Cognitive_MungeSumstats
Cognitive_sumstats$N<-257841
Cognitive_sumstats%>%
  dplyr::select(SNP=SNP,chr=CHR,pos=BP,effect_allele=A1,other_allele=A2,
                eaf=FRQ,beta=BETA,se=SE,pval=P,samplesize=N) ->Cognitive_orginal_data

#	Optional 
Cognitive_orginal_data$Phenotype<-"Cognitive"

write.table(Cognitive_orginal_data, file ="data/Cognitive_orginal_data.csv", sep ="," ,row.names =F, col.names =T)
#Cognitive_orginal_data<-fread("./Cognitive_orginal_data.csv")
#Cognitive_orginal_data<-subset(Cognitive_orginal_data,eaf>0.05)
Cognitive_orginal_data_exposure<-subset(Cognitive_orginal_data,pval<1e-5)
write.table(Cognitive_orginal_data_exposure, file ="data/Cognitive_orginal_data_exposure.csv", sep ="," ,row.names =F, col.names =T)

#exposure_data
Cognitive_exposure_orginal<-read_exposure_data(
  filename="data/Cognitive_orginal_data_exposure.csv",
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
  samplesize_col = "samplesize",
  gene_col = "gene",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE,
  chr_col = "chr",
  pos_col = "pos"
)
write.table(Cognitive_exposure_orginal, file ="data/Cognitive_exposure_orginal.csv", sep ="," ,row.names =F, col.names =T)
Cognitive_exposure_orginal<-fread(file ="data/Cognitive_exposure_orginal.csv")
Cognitive_exposure_orginal%>%
  dplyr::select(rsid=SNP,pval=pval.exposure,id=id.exposure) %>%
  filter(rsid!=".")->Cognitive_exposure_1E6_clump
#Cognitive<-clump_data(Cognitive_exposure_orginal,pop = "EUR")
Rownames_local<- ld_clump_local(dat = Cognitive_exposure_1E6_clump,
                                clump_kb = 10000,
                                clump_r2 = 0.001,
                                clump_p=5e-8,
                                bfile = "M:/MR_database/1kg.v3/EUR",
                                plink_bin = "M:/MR_database/plink/windows/plink.exe")
#select clumped data
Cognitive_exposure<-filter(Cognitive_exposure_orginal,SNP %in% Rownames_local$rsid)
write.table(Cognitive_exposure, file ="data/Cognitive_exposure_covariate.csv", sep ="," ,row.names =F, col.names =T)
save(Cognitive_exposure, file ="data/Cognitive_exposure.Rdata")
Cognitive_exposure<-read.table( file ="data/Cognitive_exposure_covariate.csv", sep =",",header = T)
Cognitive_exposure<-Cognitive_exposure[-which(Cognitive_exposure$SNP=="rs10129426"|
                                                Cognitive_exposure$SNP=="rs1064608"|
                                                Cognitive_exposure$SNP=="rs11210871"|
                                                Cognitive_exposure$SNP=="rs1144593"|
                                                Cognitive_exposure$SNP=="rs11662271"|
                                                Cognitive_exposure$SNP=="rs11720523"|
                                                Cognitive_exposure$SNP=="rs11793831"|
                                                Cognitive_exposure$SNP=="rs12435486"|
                                                Cognitive_exposure$SNP=="rs12535854"|
                                                Cognitive_exposure$SNP=="rs12635303"|
                                                Cognitive_exposure$SNP=="rs13107325"|
                                                Cognitive_exposure$SNP=="rs13163336"|
                                                Cognitive_exposure$SNP=="rs13428598"|
                                                Cognitive_exposure$SNP=="rs148696809"|
                                                Cognitive_exposure$SNP=="rs1507010"|
                                                Cognitive_exposure$SNP=="rs1523048"|
                                                Cognitive_exposure$SNP=="rs159428"|
                                                Cognitive_exposure$SNP=="rs17049085"|
                                                Cognitive_exposure$SNP=="rs17106817"|
                                                Cognitive_exposure$SNP=="rs1906252"|
                                                Cognitive_exposure$SNP=="rs2239647"|
                                                Cognitive_exposure$SNP=="rs2352974"|
                                                Cognitive_exposure$SNP=="rs2426132"|
                                                Cognitive_exposure$SNP=="rs2652454"|
                                                Cognitive_exposure$SNP=="rs3128341"|
                                                Cognitive_exposure$SNP=="rs35390433"|
                                                Cognitive_exposure$SNP=="rs3843954"|
                                                Cognitive_exposure$SNP=="rs4342312"|
                                                Cognitive_exposure$SNP=="rs4463213"|
                                                Cognitive_exposure$SNP=="rs4470366"|
                                                Cognitive_exposure$SNP=="rs4744250"|
                                                Cognitive_exposure$SNP=="rs4976976"|
                                                Cognitive_exposure$SNP=="rs5751191"|
                                                Cognitive_exposure$SNP=="rs62065449"|
                                                Cognitive_exposure$SNP=="rs6798941"|
                                                Cognitive_exposure$SNP=="rs6975134"|
                                                Cognitive_exposure$SNP=="rs702222"|
                                                Cognitive_exposure$SNP=="rs7573001"|
                                                Cognitive_exposure$SNP=="rs78382112"|
                                                Cognitive_exposure$SNP=="rs9930063"),]

save(Cognitive_exposure, file ="data/Cognitive_exposure.Rdata")
