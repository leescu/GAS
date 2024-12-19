# +================================================+ ####
# +====Section 0. Packages and setwd===============+ ####
# +================================================+ #### 
setwd("D:/PD_MDD/step 0. Data")
library(base)
library(dplyr)
library(data.table)
library(tidyverse)
# +================================================+ ####
# +====Section 1. Data arrangement============+ ####
# +================================================+ ####  
MDD_orginal = fread('9. PLACO/MDD_summary.txt')
head(MDD_orginal)
PD_orginal = fread('9. PLACO/PD_summary.txt')
head(PD_orginal)
df_fianl=merge(MDD_orginal,PD_orginal,by='SNP',all=F)
df_fianl=na.omit(df_fianl)
Z.matrix=select(df_fianl,c('SNP','Z.x','Z.y'))
Z.matrix_SNP<-Z.matrix$SNP
row.names(Z.matrix)<-Z.matrix$SNP
Z.matrix=Z.matrix[,-1]
Z.matrix=as.matrix(Z.matrix)
row.names(Z.matrix)<-Z.matrix_SNP
P.matrix=select(df_fianl,c('SNP','P.x','P.y'))
P.matrix_SNP<-P.matrix$SNP
row.names(P.matrix)<-P.matrix$SNP
P.matrix=P.matrix[,-1]
P.matrix=as.matrix(P.matrix)
row.names(P.matrix)<-P.matrix_SNP

# +================================================+ ####
# +====Section 2. Function for calculation=========+ ####
# +================================================+ ####  
.pdfx<-function(x) besselK(x=abs(x),nu=0)/pi
.p.bessel<-function(z, varz, AbsTol=1e-13){
  p1<-2*as.double(integrate(Vectorize(.pdfx), abs(z[1]*z[2]/sqrt(varz[1])),Inf, abs.tol=AbsTol)$value)
  p2<-2*as.double(integrate(Vectorize(.pdfx), abs(z[1]*z[2]/sqrt(varz[2])),Inf, abs.tol=AbsTol)$value)
  p0<-2*as.double(integrate(Vectorize(.pdfx), abs(z[1]*z[2]),Inf, abs.tol=AbsTol)$value)
  pval.compnull<-p1+p2-p0
  return(pval.compnull)
}
# +================================================+ ####
# +====Section 3. Function for PLACO===============+ ####
# +================================================+ #### 

var.placo<-function(Z.matrix, P.matrix, p.threshold=1e-4){
  # Here Z.matrix is the pxk matrix of Z-scores where p is total no. of variants in the dataset, and k is the no. of traits 
  # Similarly, P.matrix is the corresponding pxk matrix of p-values where p is total no. of variants in the dataset, and k is the no. of traits
  # p.threshold determines which variants are deemed to have no association marginally (default: 1e-4)
  # checks
  k<-ncol(Z.matrix)
  if(k!=2) stop("This method is meant for 2 traits only. Columns correspond to traits.")
  ZP<-cbind(Z.matrix,P.matrix)
  ZP<-na.omit(ZP)
  
  rows.alt<-which(ZP[,3]<p.threshold & ZP[,4]<p.threshold)
  if(length(rows.alt)>0){
    ZP<-ZP[-rows.alt,]
    if(nrow(ZP)==0) stop(paste("No 'null' variant left at p-value threshold",p.threshold))
    if(nrow(ZP)<30) warning(paste("Too few 'null' variants at p-value threshold",p.threshold))
  }
  varz<-diag(var(ZP[,c(1,2)]))
  return(varz)
}
# +================================================+ ####
# +====Section 4. Function for matrix of the Z's===+ ####
# +================================================+ #### 

cor.pearson<-function(Z.matrix, P.matrix, p.threshold=1e-4){
  # Here Z.matrix is the pxk matrix of Z-scores where p is total no. of variants in the dataset, and k is the no. of traits 
  # Similarly, P.matrix is the corresponding pxk matrix of p-values where p is total no. of variants in the dataset, and k is the no. of traits
  # p.threshold determines which variants are deemed to have no association marginally (default: 1e-4)
  # checks
  k<-ncol(Z.matrix)
  if(k!=2) stop("This method is meant for 2 traits only.")
  # estimating correlation
  row.exclude<-which( apply(P.matrix, MARGIN = 1, function(x) any(x < p.threshold)) == TRUE )
  if(length(row.exclude)>0) Z.matrix<-Z.matrix[-row.exclude,]
  R<-cor(Z.matrix)
  return(R)
}

placo<-function(Z, VarZ, AbsTol=.Machine$double.eps^0.8){
  # Z: vector of Z-scores of size k=2 (i.e., collection of Z-scores of a particular SNP for k=2 traits)
  # VarZ: vector of variances of Z-scores (covariance assumed 0; so need to be independent traits)
  # AbsTol: absolute tolerance (accuracy paramater) for numerical integration.
  # checks				
  k<-length(Z)		
  if(k!=2) stop("This method is meant for 2 traits only.")
  if(length(VarZ)!=k) stop("Provide variance estimates for 2 traits as obtained using var.placo() function.")
  
  # test of pleiotropy: PLACO
  pvalue.b=.p.bessel(z=Z, varz=VarZ, AbsTol=AbsTol)
  return(list(T.placo=prod(Z), p.placo=pvalue.b))
}

# +================================================+ ####
# +====Section 5. Calculation for PLACO============+ ####
# +================================================+ #### 

set.seed(1)
VarZ=var.placo(Z.matrix = Z.matrix,P.matrix = P.matrix,p.threshold = 1e-04)
out=sapply(1:nrow(Z.matrix), function(i) placo(Z=Z.matrix[i,],VarZ=VarZ))
out.data<-as.data.frame(t(out))

out.data1<-out.data
out.data1$T.placo<-as.numeric(out.data1$T.placo)
out.data1$p.placo<-as.numeric(out.data1$p.placo)
write.table(out.data1,file='9. PLACO/out_data.txt',sep='\t')
str(out.data)  
PLACO_data<-as.data.frame(df_fianl$SNP)
colnames(PLACO_data)<-"SNP"

SNP=fread("M:/MR_database/snp150_hg19.txt",data.table =T)
head(SNP)
colnames(SNP)[2]<-"SNP"
PLACO_dat<-merge(PLACO_data,SNP,by="SNP",all.x=T)
PLACO_dat<- distinct(PLACO_dat, SNP, .keep_all = TRUE) 
PLACO_dat1<-separate(PLACO_dat,`chromosome:start`, c("chromosome", "BP"), ":", extra = "merge")
out.data$`T.placo`<-as.numeric(out.data$`T.placo`)
out.data$`p.placo`<-as.numeric(out.data$`p.placo`)
PLACO_original<-cbind(PLACO_dat1,out.data)
PLACO_original<- PLACO_original[order(PLACO_original$p.placo), ]
PLACO_original<- PLACO_original[!duplicated(PLACO_original$SNP), ]
write.table(PLACO_original, file = "9. PLACO/PLACO_original.txt", sep = "\t", row.names = FALSE)
write.table(Z.matrix,file='9. PLACO/Z.matrix.txt',sep='\t',quote = F)

# +================================================+ ####
# +====Section 6. Selection for PLACO SNP==========+ ####
# +================================================+ ####
# >>>>> section 6.1. P <5e-8  ####
PLACO_original<-fread("9. PLACO/PLACO_original.txt",data.table =T)
PLACO_subset<-subset(PLACO_original,p.placo<5e-8)
PLACO_subset$T.placo<-NULL
colnames(PLACO_subset)<-c("SNP","CHR","BP","P")
write.table(PLACO_subset, file = "9. PLACO/PLACO_sig_5e8.txt", sep = "\t", row.names = FALSE,quote=F)


# >>>>> section 6.1. P <1e-6  ####
#PLACO_original<-fread("PLACO_original.txt")
PLACO_subset<-subset(PLACO_original,p.placo<1e-6)
PLACO_subset$T.placo<-NULL
colnames(PLACO_subset)<-c("SNP","CHR","BP","P")
write.table(PLACO_subset, file = "9. PLACO/PLACO_sig_1e6.txt", sep = "\t", row.names = FALSE,quote=F)

# >>>>> section 6.1. P.adj <0.05  ####
#PLACO_original<-fread("PLACO_original.txt")
PLACO_original$p.adj<-p.adjust(PLACO_original$p.placo,  # P值列表
                               method ="fdr"                       # FDR校正的方法
)
PLACO_subset<-subset(PLACO_original,p.adj<0.05)
PLACO_subset<- PLACO_subset[order(PLACO_subset$p.adj), ]
PLACO_subset$T.placo<-NULL
PLACO_subset$p.adj<-NULL
colnames(PLACO_subset)<-c("SNP","CHR","BP","P")
PLACO_subset<- PLACO_subset[PLACO_subset$CHR %in% c(1:22), ]
write.table(PLACO_subset, file = "9. PLACO/PLACO_sig_fdr.txt", sep = "\t", row.names = FALSE,quote=F)




