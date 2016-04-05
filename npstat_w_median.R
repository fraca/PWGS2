####################################################################################################
####################################################################################################
#function calculate maean, median quantiles and weighted from NPSTAT table
##INPUT

# file_in="/home/marco/pool_first/ge_int/p07C_5000_CDS_SNP_scaf/scaffold_1_snape.npstats"
# file_in="/home/marco/pool_first/ge_int/p14A_5000_CDS_SNP_scaf/scaffold_1_snape.npstats"
# sel_bp=1 # min bp in windows
# sel_S=0 #min SNP in windows
# tab=read.table(file_in,header = TRUE)
# # aa=ri_we2(tab,sel_bp,sel_S)
# 
# #OUTPUT
# #return
# #total bp
# #total number of SNPs, 
# #for Watterson Pi Tajima_D FayWu_H:
# # mean
# # weighted mean (number of SNPs in windows)
# # quantile c(0.05,0.25,0.5,0.75,0.95) 0.5 is the median
# 
# file_in="/home/marco/pool_first/snp_BED_end/p11AC_CDS_SNP_scaf/scaffold_1_varscan.npstats"
# # file_in="/home/marco/pool_first/snp_BED_end/p14C_CDS_SNP_scaf/scaffold_1_varscan.npstats"
# # file_in="/home/marco/pool_first/snp_BED_end/p07F_CDS_SNP_scaf/scaffold_1_varscan.npstats"
# # # sel_bp=500 # min bp in windows
# # # sel_S=10 #min SNP in windows
# sel_bp=1 # min bp in windows
# sel_S=0 #min SNP in windows
# # 
#tab=read.table(file_in,header = TRUE)
#bla=ri_we2(tab,sel_bp,sel_S)

# tab=tab[which(tab$length>=sel_bp),]
# tab=tab[which(tab$S>=sel_S),]
# tab[1:4,]
# hist(tab$Watterson,breaks=100)
# hist(tab$FayWu_H,breaks=100)

ri_we2=function(tab,sel_bp,sel_S) {
  require(Hmisc)
  
  tab=tab[which(tab$length>=sel_bp),]
  tab=tab[which(tab$S>=sel_S),]
  tot_bp=sum(tab$length)
  tot_S=sum(tab$S)
  tot_nonsyn_pol=sum(tab$nonsyn_pol)
  tot_syn_pol=sum(tab$syn_pol)
  tot_nonsyn_div=sum(tab$nonsyn_div)
  tot_syn_div=sum(tab$syn_div)
  ##alpha of NPStat is wrong revese
  
  NI=(tab$syn_div*tab$nonsyn_pol)/(tab$nonsyn_div*tab$syn_pol)
  alpha2=1-NI
  Dos=(tab$nonsyn_div/(tab$nonsyn_div+tab$syn_div))-(tab$nonsyn_pol/(tab$nonsyn_pol+tab$syn_pol))
  
  Wat_meanw=weighted.mean(tab$Watterson,tab$length,na.rm=TRUE)
  Wat_mean=mean(tab$Watterson,na.rm=TRUE)
  Wat_qua=quantile(tab$Watterson,c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=TRUE)
  Wat_qua2=wtd.quantile(tab$Watterson, weights=tab$length, probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
  
  Pi_meanw=weighted.mean(tab$Pi,tab$length,na.rm=TRUE)
  Pi_mean=mean(tab$Pi,na.rm=TRUE)
  Pi_qua=quantile(tab$Pi,c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=TRUE)
  Pi_qua2=wtd.quantile(tab$Pi, weights=tab$length, probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
  
  Tajima_D_meanw=weighted.mean(tab$Tajima_D,tab$length,na.rm=TRUE)
  Tajima_D_mean=mean(tab$Tajima_D,na.rm=TRUE)
  Tajima_D_qua=quantile(tab$Tajima_D,c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=TRUE)
  Tajima_D_qua2=wtd.quantile(tab$Tajima_D, weights=tab$length, probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
  
  FayWu_H_meanw=weighted.mean(tab$FayWu_H,tab$length,na.rm=TRUE)
  FayWu_H_mean=mean(tab$FayWu_H,na.rm=TRUE)
  FayWu_H_qua=quantile(tab$FayWu_H,c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=TRUE)
  FayWu_H_qua2=wtd.quantile(tab$FayWu_H, weights=tab$length, probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
  
  div_meanw=weighted.mean(tab$div,tab$length,na.rm=TRUE)
  div_mean=mean(tab$div,na.rm=TRUE)
  div_qua=quantile(tab$div,c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=TRUE)
  div_qua2=wtd.quantile(tab$div, weights=tab$length, probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
  
  alpha_meanw=weighted.mean(tab$alpha[which(tab$alpha!=Inf)],tab$length[which(tab$alpha!=Inf)],na.rm=TRUE)
  alpha_mean=mean(tab$alpha[which(tab$alpha!=Inf)],na.rm=TRUE)
  alpha_qua=quantile(tab$alpha[which(tab$alpha!=Inf)],c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=TRUE)
  alpha_qua2=wtd.quantile(tab$alpha[which(tab$alpha!=Inf)], weights=tab$length[which(tab$alpha!=Inf)], probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
  
  NI_meanw=weighted.mean(NI[which(NI!=Inf)],tab$length[which(NI!=Inf)],na.rm=TRUE)
  NI_mean=mean(NI[which(NI!=Inf)],na.rm=TRUE)
  NI_qua=quantile(NI[which(NI!=Inf)],c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=TRUE)
  NI_qua2=wtd.quantile(NI[which(NI!=Inf)], weights=tab$length[which(NI!=Inf)], probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
  
  alpha2_meanw=weighted.mean(alpha2[which(alpha2!=Inf)],tab$length[which(alpha2!=Inf)],na.rm=TRUE)
  alpha2_mean=mean(alpha2[which(alpha2!=Inf)],na.rm=TRUE)
  alpha2_qua=quantile(alpha2[which(alpha2!=Inf)],c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=TRUE)
  alpha2_qua2=wtd.quantile(alpha2[which(alpha2!=Inf)], weights=tab$length[which(alpha2!=Inf)], probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
  
  Dos_meanw=weighted.mean(Dos[which(Dos!=Inf)],tab$length[which(Dos!=Inf)],na.rm=TRUE)
  Dos_mean=mean(Dos[which(Dos!=Inf)],na.rm=TRUE)
  Dos_qua=quantile(Dos[which(Dos!=Inf)],c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=TRUE)
  Dos_qua2=wtd.quantile(Dos[which(Dos!=Inf)], weights=tab$length[which(Dos!=Inf)], probs=c(0,0.05,0.25,0.5,0.75,0.95,1))
  
  
  ris=c(tot_bp,tot_S,tot_nonsyn_pol,tot_syn_pol,tot_nonsyn_div,tot_syn_div,
        Wat_mean,Wat_meanw,Wat_qua,Wat_qua2,
        Pi_mean,Pi_meanw,Pi_qua,Pi_qua2,
        Tajima_D_mean,Tajima_D_meanw,Tajima_D_qua,Tajima_D_qua2,
        FayWu_H_mean,FayWu_H_meanw,FayWu_H_qua,FayWu_H_qua2,
        div_mean,div_meanw,div_qua,div_qua2,
        alpha_mean,alpha_meanw,alpha_qua,alpha_qua2,
        NI_mean,NI_meanw,NI_qua,NI_qua2,
        alpha2_mean,alpha2_meanw,alpha2_qua,alpha2_qua2,
        Dos_mean,Dos_meanw,Dos_qua,Dos_qua2)
  
  
  
  per_nomi=c("0","0.05","0.25","0.5","0.75","0.95","1")
  names(ris)=c("tot_bp","tot_S","tot_nonsyn_pol","tot_syn_pol","tot_nonsyn_div","tot_syn_div",
               "Wat_mean","Wat_mean_w",paste("Wat",per_nomi,sep="_"),paste("Wat_w",per_nomi,sep="_"),
               "Pi_mean","Pi_mean_w",paste("Pi",per_nomi,sep="_"),paste("Pi_w",per_nomi,sep="_"),
               "Tajima_D_mean","Tajima_D_mean_w",paste("Tajima_D",per_nomi,sep="_"),paste("Tajima_D_w",per_nomi,sep="_"),
               "FayWu_H_mean","FayWu_H_mean_w",paste("FayWu_H",per_nomi,sep="_"),paste("FayWu_H_w",per_nomi,sep="_"),
               "div_mean","div_mean_w",paste("div",per_nomi,sep="_"),paste("div_w",per_nomi,sep="_"),
               "alpha_mean","alpha_mean_w",paste("alpha",per_nomi,sep="_"),paste("alpha_w",per_nomi,sep="_"),
               "NI_mean","NI_mean_w",paste("NI",per_nomi,sep="_"),paste("NI_w",per_nomi,sep="_"),
               "alpha2_mean","alpha2_mean_w",paste("alpha2",per_nomi,sep="_"),paste("alpha2_w",per_nomi,sep="_"),
               "Dos_mean","Dos_mean_w",paste("Dos",per_nomi,sep="_"),paste("Dos_w",per_nomi,sep="_"))
  ris
  return(ris)
}



####################################################################################################
####################################################################################################
##call the function ri_we for population and scaffolds (work also on 1 scaffold)

##INPUT
# fold_pops=c( "/home/marco/pool_first/ge_int/p07C_5000_CDS_SNP_scaf/",
#              "/home/marco/pool_first/ge_int/p07D_5000_CDS_SNP_scaf/",
#              "/home/marco/pool_first/ge_int/p11A_5000_CDS_SNP_scaf/")
# 
# name_pops=c( "07C", "07D","11A")
# snp_call="snape"
# sel_bp=1
# sel_S=0
# name_out=paste("/home/marco/pool_first/gen_div/prova",snp_call,sel_bp,sel_S,sep="_")
# dist_np2(fold_pops,name_pops,snp_call,sel_bp,sel_S,name_out,scaffold=1:3)
##OUTPUT
# table with summary foreach pop.
###
dist_np2=function(fold_pops,name_pops,snp_call,sel_bp,sel_S,name_out,scaffold=1:8) {
  per_nomi=c("0","0.05","0.25","0.5","0.75","0.95","1")
  tab_tot=matrix(data=NA,ncol=150,nrow=length(fold_pops))
  colnames(tab_tot)=c("tot_bp","tot_S","tot_nonsyn_pol","tot_syn_pol","tot_nonsyn_div","tot_syn_div",
                      "Wat_mean","Wat_mean_w",paste("Wat",per_nomi,sep="_"),paste("Wat_w",per_nomi,sep="_"),
                      "Pi_mean","Pi_mean_w",paste("Pi",per_nomi,sep="_"),paste("Pi_w",per_nomi,sep="_"),
                      "Tajima_D_mean","Tajima_D_mean_w",paste("Tajima_D",per_nomi,sep="_"),paste("Tajima_D_w",per_nomi,sep="_"),
                      "FayWu_H_mean","FayWu_H_mean_w",paste("FayWu_H",per_nomi,sep="_"),paste("FayWu_H_w",per_nomi,sep="_"),
                      "div_mean","div_mean_w",paste("div",per_nomi,sep="_"),paste("div_w",per_nomi,sep="_"),
                      "alpha_mean","alpha_mean_w",paste("alpha",per_nomi,sep="_"),paste("alpha_w",per_nomi,sep="_"),
                      "NI_mean","NI_mean_w",paste("NI",per_nomi,sep="_"),paste("NI_w",per_nomi,sep="_"),
                      "alpha2_mean","alpha2_mean_w",paste("alpha2",per_nomi,sep="_"),paste("alpha2_w",per_nomi,sep="_"),
                      "Dos_mean","Dos_mean_w",paste("Dos",per_nomi,sep="_"),paste("Dos_w",per_nomi,sep="_"))
  
  rownames(tab_tot)=name_pops
  
  for(j in 1:length(fold_pops)) {
    tab_tm=matrix(data=NA,ncol=18,nrow=1)
    colnames(tab_tm)=c("window", "length", "length_outgroup", "read_depth", "S", "Watterson", "Pi",
                       "Tajima_D", "var_S", "var_Watterson", "unnorm_FayWu_H", "FayWu_H", "div", "nonsyn_pol",
                       "syn_pol", "nonsyn_div", "syn_div", "alpha")
    for(i in scaffold) {
      file_in=paste(fold_pops[j],"scaffold_",i,"_",snp_call,".npstats",sep="")
      tab=read.table(file_in,header = TRUE)
      tab_tm=rbind(tab_tm,tab)
    }
    tab_tm=tab_tm[-1,]
    riga=ri_we2(tab_tm,sel_bp,sel_S)
    tab_tot[name_pops[j],]=riga
  }
  write.table(tab_tot,sep="\t",quote=FALSE,file=name_out)
}


####################################################################################################
####################################################################################################
###New on all population, Varscan pv 0.15, min allele 3, min freq 0.015

#regions choosen: total, intergenic 1000 bp, intron, cds.
#min windows 1 or 2500 bp

name_pops=c("p07C","p07D","p07E","p07F","p07G","p07J","p07K","p07L","p07M","p07N","p07O","p07P","p07Q","p07R","p11C","p11AA","p11AB","p11AC","p11AE","p11AG","p11AH","p11AJ","p11B","p11G","p11H","p11J","p11K","p11L_F0","p11L_F1","p11M","p11N","p11O","p11P","p11Q","p11S","p11T","p11U","p11W","p11X","p11Z","p11D","p11V","p11R","p11A","p11E","p11F","p11Y","p14A","p14B","p14C","p14D","p14E","p07H","p13ARE","pHa21","pHa31")


folder="/home/marco/pool_first/snp_BED/"
nome="/home/marco/pool_first/npstat_output/sum"
snp_call="varscan"
sel_S=0
types=c("tot","intergenic_1000bp","Intron","CDS")


for(i in 1:length(types)) {
  type=types[i]
  fold_pops=paste(folder,name_pops,"_",type,"_SNP_scaf/",sep="")
  sel_bp=1
  name_out=paste(nome,type,snp_call,sel_bp,sel_S,sep="_")
  dist_np2(fold_pops,name_pops,snp_call,sel_bp,sel_S,name_out)
  sel_bp=2500
  name_out=paste(nome,type,snp_call,sel_bp,sel_S,sep="_")
  dist_np2(fold_pops,name_pops,snp_call,sel_bp,sel_S,name_out)
  
}


###number of snps

pops_sel=c("p07C", "p07D", "p07E", "p07F", "p07G", "p07H", "p07J", "p07K", "p07L", "p07M", "p07N", "p07O", "p07P", "p07Q", "p07R", "p11A", "p11AA", "p11AB", "p11AC", "p11AE", "p11AG", "p11AH", "p11AJ", "p11B", "p11C", "p11D", "p11E", "p11F", "p11G", "p11H", "p11J", "p11K", "p11L_F0", "p11M", "p11N", "p11O", "p11P", "p11Q", "p11R", "p11S", "p11T", "p11U", "p11V", "p11W", "p11X", "p11Y", "p11Z", "p14A", "p14B", "p14C", "p14D", "p14E")

scaf=paste("scaffold",1:8,sep="_")

folder="/home/marco/pool_first/snp_BED/"
region=c("tot","intergenic_1000bp","Intron","CDS")

tab_SNP=matrix(data=0,nrow=length(pops_sel),ncol=4)
tab_INDEL=matrix(data=0,nrow=length(pops_sel),ncol=4)
tab_bp=matrix(data=0,nrow=length(pops_sel),ncol=4)

for(i in 1:length(pops_sel)) {
  for(j in 1:length(region)) {
    for(z in 1:length(scaf)) {
      aa=readLines(paste(folder,pops_sel[i],"_",region[j],"_SNP_scaf/",scaf[z],"_stat",sep=""))
      tab_SNP[i,j]=tab_SNP[i,j]+as.numeric(unlist(strsplit(aa[12]," ")[[1]][1]))
      tab_INDEL[i,j]=tab_INDEL[i,j]+as.numeric(unlist(strsplit(aa[2]," ")[[1]][1]))
      bb=paste(folder,pops_sel[i],"_",region[j],"_SNP_scaf/",scaf[z],"_filt.BED",sep="")
      ri=system(paste("awk '{ris+=$3-$2} END {print ris}' ",bb,sep=""),intern=TRUE)
      tab_bp[i,j]=tab_bp[i,j]+as.numeric(ri)
    }
  }
}


colnames(tab_bp)=colnames(tab_SNP)=colnames(tab_INDEL)=c("tot","intergenic_1000bp","Intron","CDS")
rownames(tab_bp)=rownames(tab_SNP)=rownames(tab_INDEL)=pops_sel

write.table(tab_SNP,file="/home/marco/pool_first/npstat_output/tab_SNP_nf")
write.table(tab_INDEL,file="/home/marco/pool_first/npstat_output/tab_INDEL_nf")
write.table(tab_bp,file="/home/marco/pool_first/npstat_output/tab_bp_nf")


