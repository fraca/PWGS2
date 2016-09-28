
###Select SNPs based on minmum MAF (total populations), and sequenced in at least NN populations
##create input X baypass
##crea file geno e vcf
#da fare dopo snp_merge3.R




############

###selezioni popolazioni

genoMAF2=function(name_tab,min_MAF,min_pop,pop_sel=NULL,name_out) {
  
  ff=file(paste(name_out,"_",min_MAF,"_",min_pop,"_stat3",sep=""),"w")
  tab_big=read.table(name_tab,header = TRUE,stringsAsFactors = FALSE)
  
  
  if(!is.null(pop_sel)) {
    sel2_alt=paste("A_",pop_sel,sep="")
    sel2_ref=paste("R_",pop_sel,sep="")
    cat(pop_sel,sep="\n",file=paste(name_out,"_",min_MAF,"_",min_pop,"_pops",sep=""))
    sel="chiave"
    for(i in 1:length(sel2_alt)) {
      sel=c(sel,sel2_ref[i],sel2_alt[i])
    }
    tab_big=tab_big[,sel]
  }
  
  
  cat("Before",dim(tab_big)[1],"\n",file=ff)  
  sel_alt=colnames(tab_big)[grep("A_",colnames(tab_big))]
  sel_ref=colnames(tab_big)[grep("R_",colnames(tab_big))]  
  

  
  n_NA=apply(tab_big[,sel_alt],1,function(x) length(which(!is.na(x))))

  sel=which(n_NA>=min_pop)

  tab_big=tab_big[sel,]
  cat("SNPs sequenced in at least in",min_pop,"pops:", length(sel),"\n",file=ff)
  
  tot_alt=apply(tab_big[,sel_alt],1,function(x) sum(x,na.rm = T))
  tot_ref=apply(tab_big[,sel_ref],1,function(x) sum(x,na.rm = T))

  cov_tot=tot_alt+tot_ref
  freq=tot_alt/cov_tot
  MAF=freq
  MAF[which(freq>0.5)]=1-freq[which(freq>0.5)]
  #hist(MAF,breaks = 100)
  rem_MAF=which(MAF<=min_MAF)
  if(length(rem_MAF)!=0) {
    tab_big=tab_big[-rem_MAF,]
    MAF=MAF[-rem_MAF]
    cov_tot=cov_tot[-rem_MAF]  
  }
  
  
  cat("SNPs with total MAF lower than",min_MAF,":", length(rem_MAF),"\n",file=ff)
  
  cat(tab_big[,1],sep="\n",file=paste(name_out,"_",min_MAF,"_",min_pop,".vcf",sep=""))
  tab_vcf=read.table(paste(name_out,"_",min_MAF,"_",min_pop,".vcf",sep=""),sep="-")
  tab_vcf=cbind(tab_vcf[,1:2],rep(".",dim(tab_vcf)[1]),tab_vcf[,3:4],rep("99",dim(tab_vcf)[1]),rep("PASS",dim(tab_vcf)[1]),paste("AF=",round(MAF,3),";DP=",cov_tot,sep=""))
  write.table(tab_vcf,file=paste(name_out,"_",min_MAF,"_",min_pop,".vcf",sep=""),quote=F,sep="\t",row.names = F,col.names = F)
  
  tab_big=tab_big[,-1]
  write.table(tab_big,file=paste(name_out,"_",min_MAF,"_",min_pop,".geno",sep=""),quote=F,sep=" ",row.names = F,col.names = F,na="0")
  
  cat("After:", dim(tab_big)[1],"\n",file=ff)
  
  close(ff)
}




#########################

###########################################################################################################################################

#####X treemix intergenic prese da total

# out_name="/home/marco/pool_first/merge_snp4/old"
# #pop_sel=c("p07C", "p07D", "p07E", "p07F", "p07G", "p07H", "p07J", "p07K", "p07L", "p07M")
# pop_sel=c("p07C","p07D","p07E","p07F","p07G","p07J","p07K","p07L","p07M","p07N","p07O","p07P","p07Q","p07R","p11C","p11AA","p11AB","p11AC","p11AE","p11AG","p11AH","p11AJ","p11B","p11G","p11H","p11J","p11K","p11L_F0","p11M","p11N","p11O","p11P","p11Q","p11S","p11T","p11U","p11W","p11X","p11Z","p11D","p11V","p11R","p11A","p11E","p11F","p11Y","p14A","p14B","p14C","p14D","p14E","p07H")
# tab_in="/home/marco/pool_first/merge_snp3/tot_f002_scaffold_"

tab_in="/home/marco/pool_first/merge_snp4/tot_p1_f15_scaffold_"



pop_sel=c("p07C","p07D","p07E","p07F","p07G","p07J","p07K","p07L","p07M","p07N","p07O","p07P","p07Q","p07R","p11C","p11AA","p11AB","p11AC","p11AE","p11AG","p11AH","p11AJ","p11B","p11G","p11H","p11J","p11K","p11L_F0","p11M","p11N","p11O","p11P","p11Q","p11S","p11T","p11U","p11W","p11X","p11Z","p11D","p11V","p11R","p11A","p11E","p11F","p11Y","p14A","p14B","p14C","p14D","p14E","p07H","pHa31")

out_name="/home/marco/pool_first/merge_snp4/int_0"
min_MAF=0

out_name="/home/marco/pool_first/merge_snp4/int_015"
min_MAF=0.015

# out_name="/home/marco/pool_first/merge_snp4/int_05"
# min_MAF=0.05

system(paste("mkdir ",out_name,sep=""))

min_pop=length(pop_sel)
i=8
for(i in 1:8) {
  name_tab=paste(tab_in,i,sep="")
  name_out=paste(out_name,"/scaffold_",i,sep="")
  genoMAF2(name_tab,min_MAF,min_pop,pop_sel,name_out)
  cat(i,"\t")
}

# ############
# ############
# ##fai scrpit /home/marco/pool_first/merge_snp4/snpEff_by_scaf.sh
# ############
# ############
is.odd <- function(x) x %% 2 != 0

prefix=paste(out_name,"/",sep="")
suffix=paste("_",min_MAF,"_",min_pop,sep="")
file_out=paste(out_name,"_x_treemix",sep="")
##add
# file_out=paste(out_name,"_x_treemix_norm",sep="")

i=1
pops=readLines(paste(prefix,"scaffold_",i,suffix,"_pops",sep=""))
aa=read.table(paste(prefix,"scaffold_",i,suffix,".geno",sep=""))
bb=read.table(paste(prefix,"scaffold_",i,suffix,".snpEff",sep=""))
sel=which(bb[,9]=="intergenic_region")
aa=aa[sel,]
ref=aa[,is.odd(1:dim(aa)[2])]
alt=aa[,!is.odd(1:dim(aa)[2])]
##add
# tot=ref+alt
# ref=round(ref/tot*100)
# alt=round(alt/tot*100)

for(j in 1:dim(ref)[2]) {
  ref[,j]=paste(ref[,j],alt[,j],sep=",")
}
colnames(ref)=pops
write.table(ref,file=file_out,quote=F,sep=" ",row.names=F)

for(i in 2:8) {
  aa=read.table(paste(prefix,"scaffold_",i,suffix,".geno",sep=""))
  bb=read.table(paste(prefix,"scaffold_",i,suffix,".snpEff",sep=""))
  sel=which(bb[,9]=="intergenic_region")
  aa=aa[sel,]
  ref=aa[,is.odd(1:dim(aa)[2])]
  alt=aa[,!is.odd(1:dim(aa)[2])]
  for(j in 1:dim(ref)[2]) {
    ref[,j]=paste(ref[,j],alt[,j],sep=",")
  }
  write.table(ref,file=file_out,quote=F,sep=" ",row.names=F,col.names=F,append = T)
  
}

#system(paste("rm -r ",out_name,sep=""))

###################################################################


#####X treemix total regions


tab_in="/home/marco/pool_first/merge_snp4/tot_p1_f15_scaffold_"

pop_sel=c("p07C","p07D","p07E","p07F","p07G","p07J","p07K","p07L","p07M","p07N","p07O","p07P","p07Q","p07R","p11C","p11AA","p11AB","p11AC","p11AE","p11AG","p11AH","p11AJ","p11B","p11G","p11H","p11J","p11K","p11L_F0","p11M","p11N","p11O","p11P","p11Q","p11S","p11T","p11U","p11W","p11X","p11Z","p11D","p11V","p11R","p11A","p11E","p11F","p11Y","p14A","p14B","p14C","p14D","p14E","p07H","pHa31")

out_name="/home/marco/pool_first/merge_snp4/tot_05"
min_MAF=0.05

out_name="/home/marco/pool_first/merge_snp4/tot_0"
min_MAF=0

system(paste("mkdir ",out_name,sep=""))

min_pop=length(pop_sel)
i=8
for(i in 1:8) {
  name_tab=paste(tab_in,i,sep="")
  name_out=paste(out_name,"/scaffold_",i,sep="")
  genoMAF2(name_tab,min_MAF,min_pop,pop_sel,name_out)
  cat(i,"\t")
}

is.odd <- function(x) x %% 2 != 0

prefix=paste(out_name,"/",sep="")
suffix=paste("_",min_MAF,"_",min_pop,sep="")
file_out=paste(out_name,"_x_treemix",sep="")

i=1
pops=readLines(paste(prefix,"scaffold_",i,suffix,"_pops",sep=""))
aa=read.table(paste(prefix,"scaffold_",i,suffix,".geno",sep=""))
ref=aa[,is.odd(1:dim(aa)[2])]
alt=aa[,!is.odd(1:dim(aa)[2])]

for(j in 1:dim(ref)[2]) {
  ref[,j]=paste(ref[,j],alt[,j],sep=",")
}
colnames(ref)=pops
write.table(ref,file=file_out,quote=F,sep=" ",row.names=F)

for(i in 2:8) {
  aa=read.table(paste(prefix,"scaffold_",i,suffix,".geno",sep=""))
  ref=aa[,is.odd(1:dim(aa)[2])]
  alt=aa[,!is.odd(1:dim(aa)[2])]
  for(j in 1:dim(ref)[2]) {
    ref[,j]=paste(ref[,j],alt[,j],sep=",")
  }
  write.table(ref,file=file_out,quote=F,sep=" ",row.names=F,col.names=F,append = T)
  
}


system(paste("rm -r ",out_name,sep=""))


#########################################
##su int 1000bp
tab_in="/home/marco/pool_first/merge_snp4/int_p53_f15_scaffold_"

pop_sel=c("p07C","p07D","p07E","p07F","p07G","p07J","p07K","p07L","p07M","p07N","p07O","p07P","p07Q","p07R","p11C","p11AA","p11AB","p11AC","p11AE","p11AG","p11AH","p11AJ","p11B","p11G","p11H","p11J","p11K","p11L_F0","p11M","p11N","p11O","p11P","p11Q","p11S","p11T","p11U","p11W","p11X","p11Z","p11D","p11V","p11R","p11A","p11E","p11F","p11Y","p14A","p14B","p14C","p14D","p14E","p07H","pHa31")

out_name="/home/marco/pool_first/merge_snp4/int_1000bp_05"
min_MAF=0.05

system(paste("mkdir ",out_name,sep=""))

min_pop=length(pop_sel)
i=8
for(i in 1:8) {
  name_tab=paste(tab_in,i,sep="")
  name_out=paste(out_name,"/scaffold_",i,sep="")
  genoMAF2(name_tab,min_MAF,min_pop,pop_sel,name_out)
  cat(i,"\t")
}


is.odd <- function(x) x %% 2 != 0

prefix=paste(out_name,"/",sep="")
suffix=paste("_",min_MAF,"_",min_pop,sep="")
file_out=paste(out_name,"_x_treemix",sep="")

i=1
pops=readLines(paste(prefix,"scaffold_",i,suffix,"_pops",sep=""))
aa=read.table(paste(prefix,"scaffold_",i,suffix,".geno",sep=""))

ref=aa[,is.odd(1:dim(aa)[2])]
alt=aa[,!is.odd(1:dim(aa)[2])]
##add
# tot=ref+alt
# ref=round(ref/tot*100)
# alt=round(alt/tot*100)

for(j in 1:dim(ref)[2]) {
  ref[,j]=paste(ref[,j],alt[,j],sep=",")
}
colnames(ref)=pops
write.table(ref,file=file_out,quote=F,sep=" ",row.names=F)

for(i in 2:8) {
  aa=read.table(paste(prefix,"scaffold_",i,suffix,".geno",sep=""))
  ref=aa[,is.odd(1:dim(aa)[2])]
  alt=aa[,!is.odd(1:dim(aa)[2])]
  for(j in 1:dim(ref)[2]) {
    ref[,j]=paste(ref[,j],alt[,j],sep=",")
  }
  write.table(ref,file=file_out,quote=F,sep=" ",row.names=F,col.names=F,append = T)
  
}


#############################
##scaffold separati

# pop_sel=c("p07C","p07D","p07E","p07F","p07G","p07J","p07K","p07L","p07M","p07N","p07O","p07P","p07Q","p07R","p11C","p11AA","p11AB","p11AC","p11AE","p11AG","p11AH","p11AJ","p11B","p11G","p11H","p11J","p11K","p11L_F0","p11M","p11N","p11O","p11P","p11Q","p11S","p11T","p11U","p11W","p11X","p11Z","p11D","p11V","p11R","p11A","p11E","p11F","p11Y","p14A","p14B","p14C","p14D","p14E","p07H","pHa31")
# 
# out_name="/home/marco/pool_first/merge_snp4/int_0"
# min_MAF=0
# min_pop=length(pop_sel)
# 
# is.odd <- function(x) x %% 2 != 0
# 
# prefix=paste(out_name,"/",sep="")
# suffix=paste("_",min_MAF,"_",min_pop,sep="")
# file_out=paste(out_name,"_x_treemix_scaf",sep="")
# 
# pops=readLines(paste(prefix,"scaffold_1",suffix,"_pops",sep=""))
# 
# for(i in 1:8) {
#   aa=read.table(paste(prefix,"scaffold_",i,suffix,".geno",sep=""))
#   bb=read.table(paste(prefix,"scaffold_",i,suffix,".snpEff",sep=""))
#   sel=which(bb[,9]=="intergenic_region")
#   aa=aa[sel,]
#   ref=aa[,is.odd(1:dim(aa)[2])]
#   alt=aa[,!is.odd(1:dim(aa)[2])]
#   for(j in 1:dim(ref)[2]) {
#     ref[,j]=paste(ref[,j],alt[,j],sep=",")
#   }
#   colnames(ref)=pops
#   write.table(ref,file=paste(file_out,i,sep=""),quote=F,sep=" ",row.names=F)
#   
# }



##############################
###solo scaf_1

out_name="/home/marco/pool_first/merge_snp4/scaf1_015"
min_MAF=0.015
out_name="/home/marco/pool_first/merge_snp4/scaf1_05"
min_MAF=0.05
#pop_sel=c("p07C", "p07D", "p07E", "p07F", "p07G", "p07H", "p07J", "p07K", "p07L", "p07M")
pop_sel=c("p07C","p07D","p07E","p07F","p07G","p07J","p07K","p07L","p07M","p07N","p07O","p07P","p07Q","p07R","p11C","p11AA","p11AB","p11AC","p11AE","p11AG","p11AH","p11AJ","p11B","p11G","p11H","p11J","p11K","p11L_F0","p11M","p11N","p11O","p11P","p11Q","p11S","p11T","p11U","p11W","p11X","p11Z","p11D","p11V","p11R","p11A","p11E","p11F","p11Y","p14A","p14B","p14C","p14D","p14E","p07H","pHa31")

out_name="/home/marco/pool_first/merge_snp4/scaf1_015_noadmx"
min_MAF=0.015
pop_sel=c("p07C","p07D","p07E","p07F","p07G","p07J","p07K","p07L","p07M","p07N","p07O","p07P","p07Q","p07R","p11C","p11AA","p11AB","p11AC","p11AE","p11AG","p11AH","p11AJ","p11B","p11G","p11H","p11J","p11K","p11L_F0","p11M","p11N","p11O","p11P","p11Q","p11S","p11T","p11U","p11W","p11X","p11Z","p11D","p11V","p11R","p11A","p11E","p11F","p11Y","p14A","p14C","pHa31")

tab_in="/home/marco/pool_first/merge_snp4/tot_p1_f15_scaffold_"
###

system(paste("mkdir ",out_name,sep=""))

min_pop=length(pop_sel)
i=1


name_tab=paste(tab_in,i,sep="")
name_out=paste(out_name,"/scaffold_",i,sep="")
genoMAF2(name_tab,min_MAF,min_pop,pop_sel,name_out)  

# ############
# ############
# ##fai scrpit /home/marco/pool_first/merge_snp4/snpEff_by_scaf.sh
# ############
# ############
is.odd <- function(x) x %% 2 != 0

prefix=paste(out_name,"/",sep="")
suffix=paste("_",min_MAF,"_",min_pop,sep="")
file_out=paste(out_name,"_x_treemix",sep="")


i=1
pops=readLines(paste(prefix,"scaffold_",i,suffix,"_pops",sep=""))
aa=read.table(paste(prefix,"scaffold_",i,suffix,".geno",sep=""))
bb=read.table(paste(prefix,"scaffold_",i,suffix,".snpEff",sep=""))
sel=which(bb[,9]=="intergenic_region")
aa=aa[sel,]
ref=aa[,is.odd(1:dim(aa)[2])]
alt=aa[,!is.odd(1:dim(aa)[2])]
for(j in 1:dim(ref)[2]) {
  ref[,j]=paste(ref[,j],alt[,j],sep=",")
}
colnames(ref)=pops
write.table(ref,file=file_out,quote=F,sep=" ",row.names=F)

aa=read.table(paste(prefix,"scaffold_",i,suffix,".geno",sep=""))
bb=read.table(paste(prefix,"scaffold_",i,suffix,".snpEff",sep=""))
ref=aa[,is.odd(1:dim(aa)[2])]
alt=aa[,!is.odd(1:dim(aa)[2])]
for(j in 1:dim(ref)[2]) {
  ref[,j]=paste(ref[,j],alt[,j],sep=",")
}
colnames(ref)=pops
write.table(ref,file=paste(file_out,"_all",sep=""),quote=F,sep=" ",row.names=F)

#system(paste("rm -r ",out_name,sep=""))

###########################################################################################################################################

##outcrossing


pop_sel=c("p07C", "p07D", "p07E", "p07F", "p07G", "p07H", "p07J", "p07K", "p07L", "p07M", "p07N", "p07O", "p07P", "p11A", "p11AA", "p11AB", "p11AE", "p11AG", "p11AH", "p11AJ", "p11D", "p11E", "p11F", "p11G", "p11H", "p11J", "p11L_F0", "p11M", "p11N", "p11O", "p11P", "p11Q", "p11R", "p11S", "p11T", "p11U", "p11V", "p11W", "p11X", "p11Y", "p11Z", "p14C")

out_name="/home/marco/pool_first/merge_snp4/outc_05"
min_pop=21
min_MAF=0.05 # total MAF

tab_in="/home/marco/pool_first/merge_snp4/tot_p1_f15_scaffold_"

system(paste("mkdir ",out_name,sep=""))

i=8
for(i in 1:8) {
  name_tab=paste(tab_in,i,sep="")
  name_out=paste(out_name,"/scaffold_",i,sep="")
  genoMAF2(name_tab,min_MAF,min_pop,pop_sel,name_out)  
}

# ############
# ############
# ##fai scrpit /home/marco/pool_first/merge_snp4/snpEff_by_scaf.sh
# ############
# ############

###selezione per omega solo intergenic

prefix="/home/marco/pool_first/merge_snp4/outc_05/"
min_pop=21
min_MAF=0.05 # total MAF
suffix=paste("_",min_MAF,"_",min_pop,sep="")
file_out="/home/marco/pool_first/merge_snp4/outc_int_05_omega.geno"


i=1
pops=readLines(paste(prefix,"scaffold_",i,suffix,"_pops",sep=""))
aa=read.table(paste(prefix,"scaffold_",i,suffix,".geno",sep=""))
bb=read.table(paste(prefix,"scaffold_",i,suffix,".snpEff",sep=""))
sel=which(bb[,9]=="intergenic_region")
aa=aa[sel,]
write.table(aa,file=file_out,quote=F,sep=" ",row.names=F,col.names=F)

for(i in 2:8) {
  aa=read.table(paste(prefix,"scaffold_",i,suffix,".geno",sep=""))
  write.table(aa,file=file_out,quote=F,sep=" ",row.names=F,col.names=F,append = T)
  
}


##############################
##X casa
###solo scaf_2 prova x baypass su tutto lo scaffold

pop_sel=c("p07C", "p07D", "p07E", "p07F", "p07G", "p07H", "p07J", "p07K", "p07L", "p07M", "p07N", "p07O", "p07P", "p11A", "p11AA", "p11AB", "p11AE", "p11AG", "p11AH", "p11AJ", "p11D", "p11E", "p11F", "p11G", "p11H", "p11J", "p11L_F0", "p11M", "p11N", "p11O", "p11P", "p11Q", "p11R", "p11S", "p11T", "p11U", "p11V", "p11W", "p11X", "p11Y", "p11Z", "p14C")

out_name="/home/marco/pool_first/merge_snp4/outc_05"
min_pop=21
min_MAF=0.05 # total MAF

tab_in="/home/marco/pool_first/merge_snp4/tot_p1_f15_scaffold_"

system(paste("mkdir ",out_name,sep=""))

i=8

name_tab=paste(tab_in,i,sep="")
name_out=paste(out_name,"/scaffold_",i,sep="")
genoMAF2(name_tab,min_MAF,min_pop,pop_sel,name_out)  

#################################################################################################################
#################################################################################################################

##create intergenic SNPs table away from genes


##intergenic far away 2000 bp

snpeff_file="/home/marco/pool_first/merge_snp4/outc_05_42_42.snpEff"

name_out="/home/marco/pool_first/merge_snp4/outc_05_42_42_int_2000bp"
int_2000="/home/marco/pool_first/create_BED/intergenic_v2_2000bp.BED"

name_out="/home/marco/pool_first/merge_snp4/outc_05_42_42_int_1000bp"
int_2000="/home/marco/pool_first/create_BED/intergenic_v2_1000bp.BED"

name_out="/home/marco/pool_first/merge_snp4/outc_05_42_42_int_500bp"
int_2000="/home/marco/pool_first/create_BED/intergenic_v2_500bp.BED"

snpeff_file="/home/marco/pool_first/merge_snp4/outc_05_21.snpEff"

name_out="/home/marco/pool_first/merge_snp4/outc_05_21_int_2000bp"
int_2000="/home/marco/pool_first/create_BED/intergenic_v2_2000bp.BED"

name_out="/home/marco/pool_first/merge_snp4/outc_05_21_int_1000bp"
int_2000="/home/marco/pool_first/create_BED/intergenic_v2_1000bp.BED"

name_out="/home/marco/pool_first/merge_snp4/outc_05_21_int_500bp"
int_2000="/home/marco/pool_first/create_BED/intergenic_v2_500bp.BED"

name_temp=paste(name_out,"_temp",sep="")
tab_snpeff=read.table(snpeff_file,stringsAsFactors = F)
tab_int=tab_snpeff[which(tab_snpeff[,9]=="intergenic_region"),c(1,2,2,3,4,5,6,7)]
write.table(tab_int,quote=F,row.names=F,col.names=F,sep="\t",file=name_temp)
system(paste("bedtools sort -i ",name_temp," | bedtools intersect -a stdin -b ",int_2000," -f 1 > ",name_out,sep=""))
system(paste("rm -r ",name_temp,sep=""))


