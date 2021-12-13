#New PCA from filtered VCF of filteredVCF
library(vcfR)
library(ff)
library(ffbase)
library(pcaMethods)
library(ggbiplot)
library(tidyverse)
#library(GenomicRanges)

extract_gt_ff <- function(vcf){
  
  gt <- vcf@gt[,-1]
  names_gt <- colnames(gt)
  #maybe reverse?? to get number of ALT 
  gt[grepl("0/1|0\\|1", gt, perl = TRUE)] <- 1
  gt[grepl("1/1|1\\|1", gt, perl = TRUE)] <- 2
  gt[grepl("0/0|0\\|0", gt, perl = TRUE)] <- 0
  gt[grepl("./.|.\\|.", gt, perl = TRUE)] <- NA
  gt <- as.ffdf(as.ff(as.byte(gt), dim = dim(gt)))
  colnames(gt) <- names_gt
  return(gt)
}

##INPUT

vcf_file="/Users/fraca/data_Cgran/plastid/Cgl_chl.recode.vcf" ###Input vcfile
ind_pop_file='/Users/fraca/data_Cgran/plastid/li_plast_noCbp2' ##INput table first column is the sample_id, second column the population/group id for the PCA plot.
out_prefix='/Users/fraca/data_Cgran/plastid/Cgl_chl' #Output prefix

per_keep_ind=0.2 # individuals with more than X SNP missing will be removed
per_keep_snp=0.2 # SNPs with more than X individual missing will be removed
min_med_DP=10 # minimun median DP

######################################

vcf <- read.vcfR(vcf_file, limit = 1e+09, verbose = F, convertNA = FALSE)
vcf@fix[,"ID"] <- NA
chr_pos=paste0(vcf@fix[,'CHROM'],"_",vcf@fix[,'POS'])

##test the DP of each individual
dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)
gt_ff=extract_gt_ff(vcf)
gf_ff_ori=gt_ff
dp_ind=apply(dp,2,function(x) median(x,na.rm = T))

#remove the individual with less than min DP
# rem=which(dp_ind<min_med_DP)
# dp=dp[,-rem]
# gt_ff=gt_ff[,-rem]


#check and remvove ind with more than XX SNP missing
check_NA_ind=apply(gt_ff[,1:length(gt_ff)], 2, function(x) length(which(is.na(x))))
rem=which(check_NA_ind >= (dim(gt_ff)[1]*per_keep_ind))
if(length(snp_rem)!=0) {
  dp=dp[,-rem]
  gt_ff=gt_ff[,-rem]
}
# check and remove SNPs with more than XX indviduals missing
check_NA_snp=apply(gt_ff[,1:length(gt_ff)], 1, function(x) length(which(is.na(x))))
snp_rem=which(check_NA_snp >= (dim(gt_ff)[2]*per_keep_snp))

cat('Removed',length(snp_rem),' snp.\n')
#Removed 8848  snp.
if(length(snp_rem)!=0) {
  gt_ff=gt_ff[-snp_rem,]
  chr_pos=chr_pos[-snp_rem]
}

rownames(gt_ff)=chr_pos

#change ind names
ind_pop=read.table(ind_pop_file,header = F,stringsAsFactors = F)
colnames(ind_pop)=c("new_ID",	"Pop_id")
rownames(ind_pop)=ind_pop[,1]
ind_pop=ind_pop[colnames(gt_ff),]
pops=unique(ind_pop$Pop_id)
#ind_counts=as.tibble(ind_pop) %>% count(Pop_id,sort=T) 
colnames(gt_ff)=ind_pop$new_ID
head(ind_pop)

cat('Analized',dim(gt_ff)[1],'SNPs  in ',dim(gt_ff)[2],'populations\n')
# Analized 158145 SNPs  in  285 populations
cat('Percentage of missing data:',length(which(is.na(gt_ff[1:nrow(gt_ff),])))/(nrow(gt_ff)*ncol(gt_ff)),'\n')

#perform a probabilistic PCA
pc <- pca(gt_ff[1:nrow(gt_ff),], nPcs=3, method="ppca")
saveRDS(pc,file=paste0(out_prefix,"_pca.rds"))
pc=readRDS(file=paste0(out_prefix,"_pca.rds"))
summary(pc)

imputed <- completeObs(pc)
imputed2=round(imputed)
fixed=apply(imputed2,1,function(x) length(unique(x)))

bla=prcomp(t(imputed2), scale = TRUE)

pdf(file=paste0(out_prefix,"_PCA_PC12.pdf"),height = 20, width = 20)
g <- ggbiplot(bla, groups = ind_pop$Pop_id,choices = 1:2, ellipse = TRUE,var.axes = FALSE,labels=colnames(gt_ff))
print(g)

pdf(file=paste0(out_prefix,"_PCA_PC34.pdf"),height = 20, width = 20)
g <- ggbiplot(bla, groups = ind_pop$Pop_id,choices = 3:4, ellipse = TRUE,var.axes = FALSE,labels=colnames(gt_ff))
print(g)
dev.off()

pdf(file=paste0(out_prefix,"_PCA_PC56.pdf"),height = 20, width = 20)
g <- ggbiplot(bla, groups = ind_pop$Pop_id,choices = 5:6, ellipse = TRUE,var.axes = FALSE,labels=colnames(gt_ff))
print(g)
dev.off()