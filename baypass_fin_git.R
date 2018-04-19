##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##These is to create The population variance matrix  (omega tables)
###create random samples of SNPs
#for doing the omega matrix
require(corrplot) ; require(ape)
source("/home/marco/program/baypass_2.1/utils/baypass_utils.R")

snp_rand=function(input,output,N_SNP,rep=3) {
  
  aa=read.table(input)

  for(i in 1:rep) {
    sel=sample(1:dim(aa)[1], n_SNP, replace=F)
    write.table(aa[sel,],quote=F, col.names = F,row.names = F,file=paste(output,"_rand",i,".geno",sep=""))
  }  
}

n_SNP=50000

input="/home/marco/pool_first/merge_snp4/outc_int_05_omega.geno"
output="/home/marco/pool_first/baypass_fin/omega_outc/ome_outc_int_05"

snp_rand(input,output,N_SNP)

cat(rep(50,42),sep=" ",file="/home/marco/pool_first/input_baypass/omega_outc/poolsize_42")


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
# Script used The population variance matrix  (omega tables) on SCIcore
# #! /bin/bash
# # Request "/bin/bash" as shell
# #$ -S /bin/bash
# #$ -o ome_1
# #$ -N ome_1
# #$ -cwd
# #$ -j y
# #$ -pe smp 6
# ##$ -q short.q
# #$ -q long.q
# ##$ -q very_long.q
# ##$ -q infinite.q
# 
# module load ifort
# bin_dir="/scicore/home/williy/fracasse/baypass_2.1/sources/"
# n_pop=42
# poolsize_file="poolsize_42"
# par_d0yij=4
# geno_file="ome_outc_int_05_rand1.geno"
# name_out="ome_outc1"
# $bin_dir"i_baypass" -npop $n_pop -gfile $geno_file -poolsizefile $poolsize_file -d0yij $par_d0yij -nthreads 6 -outprefix $name_out
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

##test similarity between matrix

mat_r1=as.matrix(read.table("/home/marco/pool_first/baypass_fin/omega_outc/ome_outc1_mat_omega.out"))
#mat_r2=as.matrix(read.table("/home/marco/pool_first/baypass_fin/omega_outc/ome_outc2_mat_omega.out"))
mat_r3=as.matrix(read.table("/home/marco/pool_first/baypass_fin/omega_outc/ome_outc3_mat_omega.out"))

#fmd.dist(mat_r1,mat_r2)
fmd.dist(mat_r1,mat_r3)
#fmd.dist(mat_r1,mat_r4)

##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
#Selection of environmental variables

cat(rep(50,42),sep=" ",file="/home/marco/pool_first/input_baypass/omega_outc/poolsize_42")


pop_sel=c("p07C", "p07D", "p07E", "p07F", "p07G", "p07H", "p07J", "p07K", "p07L", "p07M", "p07N", "p07O", "p07P", "p11A", "p11AA", "p11AB", "p11AE", "p11AG", "p11AH", "p11AJ", "p11D", "p11E", "p11F", "p11G", "p11H", "p11J", "p11L_F0", "p11M", "p11N", "p11O", "p11P", "p11Q", "p11R", "p11S", "p11T", "p11U", "p11V", "p11W", "p11X", "p11Y", "p11Z", "p14C")

aa=read.table("/home/marco/Dropbox/dottorato_Neuchatel/climate/tab_cli_160912.txt")

aa=aa[pop_sel,c("Sub","bio4","alpha","Tmin_ESp","bio16")]

##X outc
aa=t(aa)
write.table(aa,quote=FALSE,col.names = F,row.names = F,file="/home/marco/pool_first/baypass_fin/outc/cli")
cat(rownames(aa),sep="\n",file="/home/marco/pool_first/baypass_fin/outc/cli_row")

###X outc2

aa=aa[pop_sel,c("Sub","bio13","ai_yr","bio4","lgs","bio2","gcd_mix_anc")]
#cor(as.matrix(aa), method="spearman")
aa=t(aa)
write.table(aa,quote=FALSE,col.names = F,row.names = F,file="/home/marco/pool_first/baypass_fin/outc2/cli")
cat(rownames(aa),sep="\n",file="/home/marco/pool_first/baypass_fin/outc2/cli_row")


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
# This scriprt perform the analysis with the aux models on SCIcore
# #! /bin/bash
# #$ -S /bin/bash
# #$ -o work_aux.out
# #$ -N work_aux
# #$ -cwd
# #$ -j y
# #$ -pe smp 10
# ##$ -q short.q
# ##$ -q long.q
# ##$ -q very_long.q
# ##$ -q infinite.q
# module load ifort
# 
# bin_dir="/scicore/home/williy/fracasse/baypass_2.1/sources/"
# n_pop=42
# poolsize_file="poolsize_42"
# par_d0yij=4
# omega_file="ome_outc1_mat_omega.out"
# varcli_file="cli"
# geno_file=$dir"/b"${SGE_TASK_ID}
# name_out=$dir"_aux"$n"/b"${SGE_TASK_ID}
# $bin_dir"i_baypass" -npop $n_pop -gfile $geno_file -efile $varcli_file -poolsizefile $poolsize_file -d0yij $par_d0yij -scalecov -auxmodel -omegafile $omega_file -nthreads 10 -outprefix $name_out -seed $seed
# echo -e "Aux "${SGE_TASK_ID} >> $dir"_check_aux"$n


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
###merge output baypass
#INPUT
# ##folder with the baypass output
# folder="/home/marco/pool_first/input_baypass/mod_tot_f002_0.05_26/scaf_1_aux/"
##mode aux bas or std
# mode="aux"
##cli_rows rownmaes of climatic variables
# clirow="/home/marco/pool_first/input_baypass/mod_tot_f002_0.05_26/cli_fin_row"
# output directory
# dir_out="/home/marco/pool_first/input_baypass/mod_tot_f002_0.05_26/scaf_1_merge/"
# lim_wind
# to analize only some file not all

##OUTPUT
#directory with one file for each variable and XtX

################
merge_baypass=function(folder,mode,clirow,dir_out,lim_wind=NA) {
  
  clirow=readLines(clirow)
  
  #tab_vcf=read.table(vcf_file)
  
  # tab_ris=matrix(data=NA,nrow=dim(tab_vcf)[1],ncol=5+length(clirow))
  # tab_ris[,1:5]=tab_vcf[,c(1,2,4,5,8)]
  # tab_ris[1:4,]
  
  if(mode=="aux") {
    file_ris=(system(paste("ls ",folder),inter=T))
    file_ris=sort(file_ris[grep("_betai.out",file_ris)])
    if(!is.na(lim_wind)) {
      file_ris=file_ris[1:lim_wind]
    }
    ##primo
    
    i=1
    aa=read.table(paste(folder,file_ris[i],sep=""),header = TRUE)
    passo=dim(aa)[1]/length(clirow)
    
    uno=1
    due=passo
    for(j in 1:length(clirow)) {
      write.table(aa[uno:due,3:6],quote=F,row.names = F, col.names = F, file=paste(dir_out,clirow[j],"_",mode,sep=""))
      uno=due+1
      due=due+passo
      #cat(uno,due,"\n")
    }
    
    #in mezzo
    
    for(i in 2:length(file_ris)) {
      
      aa=read.table(paste(folder,file_ris[i],sep=""),header = TRUE)
      passo=dim(aa)[1]/length(clirow)
      uno=1
      due=passo
      
      for(j in 1:length(clirow)) {
        write.table(aa[uno:due,3:6],quote=F,row.names = F, col.names = F, file=paste(dir_out,clirow[j],"_",mode,sep=""),append = T)
        uno=due+1
        due=due+passo
        #cat(uno,due,"\n")
      }
    }
  }
  
  if(mode=="bas") {
    file_ris=(system(paste("ls ",folder),inter=T))
    file_ris=sort(file_ris[grep("_betai_reg.out",file_ris)])
    
    if(!is.na(lim_wind)) {
      file_ris=file_ris[1:lim_wind]
    }
    
    ##primo
    
    i=1
    aa=read.table(paste(folder,file_ris[i],sep=""),header = TRUE)
    passo=dim(aa)[1]/length(clirow)
    
    uno=1
    due=passo
    for(j in 1:length(clirow)) {
      write.table(aa[uno:due,3:8],quote=F,row.names = F, col.names = F, file=paste(dir_out,clirow[j],"_",mode,sep=""))
      uno=due+1
      due=due+passo
      #cat(uno,due,"\n")
    }
    
    #in mezzo
    
    for(i in 2:length(file_ris)) {
      
      aa=read.table(paste(folder,file_ris[i],sep=""),header = TRUE)
      passo=dim(aa)[1]/length(clirow)
      uno=1
      due=passo
      
      for(j in 1:length(clirow)) {
        write.table(aa[uno:due,3:8],quote=F,row.names = F, col.names = F, file=paste(dir_out,clirow[j],"_",mode,sep=""),append = T)
        uno=due+1
        due=due+passo
        #cat(uno,due,"\n")
      }
    }
  }
  
  if(mode=="std") {
    file_ris=(system(paste("ls ",folder),inter=T))
    file_ris=sort(file_ris[grep("_betai.out",file_ris)])
    
    if(!is.na(lim_wind)) {
      file_ris=file_ris[1:lim_wind]
    }
    ##primo
    
    i=1
    aa=read.table(paste(folder,file_ris[i],sep=""),header = TRUE)
    passo=dim(aa)[1]/length(clirow)
    
    uno=1
    due=passo
    for(j in 1:length(clirow)) {
      write.table(aa[uno:due,3:7],quote=F,row.names = F, col.names = F, file=paste(dir_out,clirow[j],"_",mode,sep=""))
      uno=due+1
      due=due+passo
      #cat(uno,due,"\n")
    }
    
    #in mezzo
    
    for(i in 2:length(file_ris)) {
      
      aa=read.table(paste(folder,file_ris[i],sep=""),header = TRUE)
      passo=dim(aa)[1]/length(clirow)
      uno=1
      due=passo
      
      for(j in 1:length(clirow)) {
        write.table(aa[uno:due,3:7],quote=F,row.names = F, col.names = F, file=paste(dir_out,clirow[j],"_",mode,sep=""),append = T)
        uno=due+1
        due=due+passo
        #cat(uno,due,"\n")
      }
    }
  }
  
  ##XtX
  
  file_ris=(system(paste("ls ",folder),inter=T))
  file_ris=sort(file_ris[grep("pi_xtx.out",file_ris)])
  
  if(!is.na(lim_wind)) {
    file_ris=file_ris[1:lim_wind]
  }
  ##primo
  
  i=1
  aa=read.table(paste(folder,file_ris[i],sep=""),header = TRUE)
  write.table(aa[,2:7],quote=F,row.names = F, col.names = F, file=paste(dir_out,"XtX_",mode,sep=""))
  #in mezzo
  
  for(i in 2:length(file_ris)) {
    aa=read.table(paste(folder,file_ris[i],sep=""),header = TRUE)
    write.table(aa[,2:7],quote=F,row.names = F, col.names = F, file=paste(dir_out,"XtX_",mode,sep=""),append = T)
  }
  
  
}

#################
merge_baypass2=function(clirow,scafs,mode,dir_out) {
  system(paste("mkdir ",dir_out))
  clirow=readLines(clirow)
  clirow=c("XtX",clirow)
  for(j in 1:length(clirow)) {
    
    sel_cli=clirow[j]  
    cat(sel_cli,"\t")
    i=1
    tab=read.table(paste(scafs[i],sel_cli,"_",mode,sep=""),stringsAsFactors = F)
    
    for(i in 2:length(scafs)) {
      tab=rbind(tab,read.table(paste(scafs[i],sel_cli,"_",mode,sep=""),stringsAsFactors = F))
    }
    write.table(tab,quote=F,row.names = F, col.names = F, file=paste(dir_out,sel_cli,"_",mode,sep=""))
    
  }
 
}

#################

# clirow="/home/marco/pool_first/baypass_fin/outc/cli_row"
# 
# for(i in 1:8) {
#   dir_out=paste("/home/marco/pool_first/baypass_fin/outc/var_scaf_",i,"_mod1/",sep="")
#   system(paste("mkdir ",dir_out))
#   mode="bas"
#   folder=paste("/home/marco/pool_first/baypass_fin/outc/scaf_",i,"_",mode,"1/",sep="")
#   merge_baypass(folder,mode,clirow,dir_out)  
#   mode="aux"
#   folder=paste("/home/marco/pool_first/baypass_fin/outc/scaf_",i,"_",mode,"1/",sep="")
#   merge_baypass(folder,mode,clirow,dir_out)  
#   
#   dir_out=paste("/home/marco/pool_first/baypass_fin/outc/var_scaf_",i,"_mod2/",sep="")
#   system(paste("mkdir ",dir_out))
#   mode="bas"
#   folder=paste("/home/marco/pool_first/baypass_fin/outc/scaf_",i,"_",mode,"2/",sep="")
#   merge_baypass(folder,mode,clirow,dir_out)  
#   mode="aux"
#   folder=paste("/home/marco/pool_first/baypass_fin/outc/scaf_",i,"_",mode,"2/",sep="")
#   merge_baypass(folder,mode,clirow,dir_out)  
#   
#   dir_out=paste("/home/marco/pool_first/baypass_fin/outc/var_scaf_",i,"_mod3/",sep="")
#   system(paste("mkdir ",dir_out))
#   mode="bas"
#   folder=paste("/home/marco/pool_first/baypass_fin/outc/scaf_",i,"_",mode,"3/",sep="")
#   merge_baypass(folder,mode,clirow,dir_out)  
#   mode="aux"
#   folder=paste("/home/marco/pool_first/baypass_fin/outc/scaf_",i,"_",mode,"3/",sep="")
#   merge_baypass(folder,mode,clirow,dir_out)  
#   
# }
# 
# 
# clirow="/home/marco/pool_first/baypass_fin/outc/cli_row"
# scafs=c("/home/marco/pool_first/baypass_fin/outc/var_scaf_1_mod1/",
#         "/home/marco/pool_first/baypass_fin/outc/var_scaf_2_mod1/",
#         "/home/marco/pool_first/baypass_fin/outc/var_scaf_3_mod1/",
#         "/home/marco/pool_first/baypass_fin/outc/var_scaf_4_mod1/",
#         "/home/marco/pool_first/baypass_fin/outc/var_scaf_5_mod1/",
#         "/home/marco/pool_first/baypass_fin/outc/var_scaf_6_mod1/",
#         "/home/marco/pool_first/baypass_fin/outc/var_scaf_7_mod1/",
#         "/home/marco/pool_first/baypass_fin/outc/var_scaf_8_mod1/")
# dir_out="/home/marco/pool_first/baypass_fin/outc/merged_mod1/"
# 
# system(paste("mkdir ",dir_out))
# mode="bas"
# merge_baypass2(clirow,scafs,mode,dir_out)
# mode="aux"
# merge_baypass2(clirow,scafs,mode,dir_out)
# ##delete scaffolds merged
# for(i in 1:length(scafs)) {
#   system(paste("rm -r ",scafs[i],sep=""))
# }
# 
# scafs=c("/home/marco/pool_first/baypass_fin/outc/var_scaf_1_mod2/",
#         "/home/marco/pool_first/baypass_fin/outc/var_scaf_2_mod2/",
#         "/home/marco/pool_first/baypass_fin/outc/var_scaf_3_mod2/",
#         "/home/marco/pool_first/baypass_fin/outc/var_scaf_4_mod2/",
#         "/home/marco/pool_first/baypass_fin/outc/var_scaf_5_mod2/",
#         "/home/marco/pool_first/baypass_fin/outc/var_scaf_6_mod2/",
#         "/home/marco/pool_first/baypass_fin/outc/var_scaf_7_mod2/",
#         "/home/marco/pool_first/baypass_fin/outc/var_scaf_8_mod2/")
# dir_out="/home/marco/pool_first/baypass_fin/outc/merged_mod2/"
# system(paste("mkdir ",dir_out))
# mode="bas"
# merge_baypass2(clirow,scafs,mode,dir_out)
# mode="aux"
# merge_baypass2(clirow,scafs,mode,dir_out)
# ##delete scaffolds merged
# for(i in 1:length(scafs)) {
#   system(paste("rm -r ",scafs[i],sep=""))
# }
# 
# scafs=c("/home/marco/pool_first/baypass_fin/outc/var_scaf_1_mod3/",
#         "/home/marco/pool_first/baypass_fin/outc/var_scaf_2_mod3/",
#         "/home/marco/pool_first/baypass_fin/outc/var_scaf_3_mod3/",
#         "/home/marco/pool_first/baypass_fin/outc/var_scaf_4_mod3/",
#         "/home/marco/pool_first/baypass_fin/outc/var_scaf_5_mod3/",
#         "/home/marco/pool_first/baypass_fin/outc/var_scaf_6_mod3/",
#         "/home/marco/pool_first/baypass_fin/outc/var_scaf_7_mod3/",
#         "/home/marco/pool_first/baypass_fin/outc/var_scaf_8_mod3/")
# dir_out="/home/marco/pool_first/baypass_fin/outc/merged_mod3/"
# system(paste("mkdir ",dir_out))
# mode="bas"
# merge_baypass2(clirow,scafs,mode,dir_out)
# mode="aux"
# merge_baypass2(clirow,scafs,mode,dir_out)
# ##delete scaffolds merged
# for(i in 1:length(scafs)) {
#   system(paste("rm -r ",scafs[i],sep=""))
# }
# 

######
##only bas more variables


clirow="/home/marco/pool_first/baypass_fin/outc2/cli_row"

for(i in 1:8) {
  dir_out=paste("/home/marco/pool_first/baypass_fin/outc2/var_scaf_",i,"_mod1/",sep="")
  system(paste("mkdir ",dir_out))
  mode="bas"
  folder=paste("/home/marco/pool_first/baypass_fin/outc2/scaf_",i,"_",mode,"1/",sep="")
  merge_baypass(folder,mode,clirow,dir_out)  
  
  dir_out=paste("/home/marco/pool_first/baypass_fin/outc2/var_scaf_",i,"_mod2/",sep="")
  system(paste("mkdir ",dir_out))
  mode="bas"
  folder=paste("/home/marco/pool_first/baypass_fin/outc2/scaf_",i,"_",mode,"2/",sep="")
  merge_baypass(folder,mode,clirow,dir_out)  
  
  dir_out=paste("/home/marco/pool_first/baypass_fin/outc2/var_scaf_",i,"_mod3/",sep="")
  system(paste("mkdir ",dir_out))
  mode="bas"
  folder=paste("/home/marco/pool_first/baypass_fin/outc2/scaf_",i,"_",mode,"3/",sep="")
  merge_baypass(folder,mode,clirow,dir_out)  
  
}


clirow="/home/marco/pool_first/baypass_fin/outc2/cli_row"
scafs=c("/home/marco/pool_first/baypass_fin/outc2/var_scaf_1_mod1/",
        "/home/marco/pool_first/baypass_fin/outc2/var_scaf_2_mod1/",
        "/home/marco/pool_first/baypass_fin/outc2/var_scaf_3_mod1/",
        "/home/marco/pool_first/baypass_fin/outc2/var_scaf_4_mod1/",
        "/home/marco/pool_first/baypass_fin/outc2/var_scaf_5_mod1/",
        "/home/marco/pool_first/baypass_fin/outc2/var_scaf_6_mod1/",
        "/home/marco/pool_first/baypass_fin/outc2/var_scaf_7_mod1/",
        "/home/marco/pool_first/baypass_fin/outc2/var_scaf_8_mod1/")
dir_out="/home/marco/pool_first/baypass_fin/outc2/merged_mod1/"

system(paste("mkdir ",dir_out))
mode="bas"
merge_baypass2(clirow,scafs,mode,dir_out)
##delete scaffolds merged
for(i in 1:length(scafs)) {
  system(paste("rm -r ",scafs[i],sep=""))
}

scafs=c("/home/marco/pool_first/baypass_fin/outc2/var_scaf_1_mod2/",
        "/home/marco/pool_first/baypass_fin/outc2/var_scaf_2_mod2/",
        "/home/marco/pool_first/baypass_fin/outc2/var_scaf_3_mod2/",
        "/home/marco/pool_first/baypass_fin/outc2/var_scaf_4_mod2/",
        "/home/marco/pool_first/baypass_fin/outc2/var_scaf_5_mod2/",
        "/home/marco/pool_first/baypass_fin/outc2/var_scaf_6_mod2/",
        "/home/marco/pool_first/baypass_fin/outc2/var_scaf_7_mod2/",
        "/home/marco/pool_first/baypass_fin/outc2/var_scaf_8_mod2/")
dir_out="/home/marco/pool_first/baypass_fin/outc2/merged_mod2/"
system(paste("mkdir ",dir_out))
mode="bas"
merge_baypass2(clirow,scafs,mode,dir_out)
##delete scaffolds merged
for(i in 1:length(scafs)) {
  system(paste("rm -r ",scafs[i],sep=""))
}

scafs=c("/home/marco/pool_first/baypass_fin/outc2/var_scaf_1_mod3/",
        "/home/marco/pool_first/baypass_fin/outc2/var_scaf_2_mod3/",
        "/home/marco/pool_first/baypass_fin/outc2/var_scaf_3_mod3/",
        "/home/marco/pool_first/baypass_fin/outc2/var_scaf_4_mod3/",
        "/home/marco/pool_first/baypass_fin/outc2/var_scaf_5_mod3/",
        "/home/marco/pool_first/baypass_fin/outc2/var_scaf_6_mod3/",
        "/home/marco/pool_first/baypass_fin/outc2/var_scaf_7_mod3/",
        "/home/marco/pool_first/baypass_fin/outc2/var_scaf_8_mod3/")
dir_out="/home/marco/pool_first/baypass_fin/outc2/merged_mod3/"
system(paste("mkdir ",dir_out))
mode="bas"
merge_baypass2(clirow,scafs,mode,dir_out)
##delete scaffolds merged
for(i in 1:length(scafs)) {
  system(paste("rm -r ",scafs[i],sep=""))
}


#############################


##only aux 5 variables


clirow="/home/marco/pool_first/baypass_fin/outc/cli_row"

for(i in 1:8) {
  dir_out=paste("/home/marco/pool_first/baypass_fin/outc/var_scaf_",i,"_mod1/",sep="")
  system(paste("mkdir ",dir_out))
  mode="aux"
  folder=paste("/home/marco/pool_first/baypass_fin/outc/scaf_",i,"_",mode,"1/",sep="")
  merge_baypass(folder,mode,clirow,dir_out)  
  
  dir_out=paste("/home/marco/pool_first/baypass_fin/outc/var_scaf_",i,"_mod2/",sep="")
  system(paste("mkdir ",dir_out))
  mode="aux"
  folder=paste("/home/marco/pool_first/baypass_fin/outc/scaf_",i,"_",mode,"2/",sep="")
  merge_baypass(folder,mode,clirow,dir_out)  
  
  dir_out=paste("/home/marco/pool_first/baypass_fin/outc/var_scaf_",i,"_mod3/",sep="")
  system(paste("mkdir ",dir_out))
  mode="aux"
  folder=paste("/home/marco/pool_first/baypass_fin/outc/scaf_",i,"_",mode,"3/",sep="")
  merge_baypass(folder,mode,clirow,dir_out)  
  
}


clirow="/home/marco/pool_first/baypass_fin/outc/cli_row"
scafs=c("/home/marco/pool_first/baypass_fin/outc/var_scaf_1_mod1/",
        "/home/marco/pool_first/baypass_fin/outc/var_scaf_2_mod1/",
        "/home/marco/pool_first/baypass_fin/outc/var_scaf_3_mod1/",
        "/home/marco/pool_first/baypass_fin/outc/var_scaf_4_mod1/",
        "/home/marco/pool_first/baypass_fin/outc/var_scaf_5_mod1/",
        "/home/marco/pool_first/baypass_fin/outc/var_scaf_6_mod1/",
        "/home/marco/pool_first/baypass_fin/outc/var_scaf_7_mod1/",
        "/home/marco/pool_first/baypass_fin/outc/var_scaf_8_mod1/")
dir_out="/home/marco/pool_first/baypass_fin/outc/merged_mod1/"

system(paste("mkdir ",dir_out))
mode="aux"
merge_baypass2(clirow,scafs,mode,dir_out)
##delete scaffolds merged
for(i in 1:length(scafs)) {
  system(paste("rm -r ",scafs[i],sep=""))
}

clirow="/home/marco/pool_first/baypass_fin/outc/cli_row"
scafs=c("/home/marco/pool_first/baypass_fin/outc/var_scaf_1_mod2/",
        "/home/marco/pool_first/baypass_fin/outc/var_scaf_2_mod2/",
        "/home/marco/pool_first/baypass_fin/outc/var_scaf_3_mod2/",
        "/home/marco/pool_first/baypass_fin/outc/var_scaf_4_mod2/",
        "/home/marco/pool_first/baypass_fin/outc/var_scaf_5_mod2/",
        "/home/marco/pool_first/baypass_fin/outc/var_scaf_6_mod2/",
        "/home/marco/pool_first/baypass_fin/outc/var_scaf_7_mod2/",
        "/home/marco/pool_first/baypass_fin/outc/var_scaf_8_mod2/")
dir_out="/home/marco/pool_first/baypass_fin/outc/merged_mod2/"

system(paste("mkdir ",dir_out))
mode="aux"
merge_baypass2(clirow,scafs,mode,dir_out)
##delete scaffolds merged
for(i in 1:length(scafs)) {
  system(paste("rm -r ",scafs[i],sep=""))
}

clirow="/home/marco/pool_first/baypass_fin/outc/cli_row"
scafs=c("/home/marco/pool_first/baypass_fin/outc/var_scaf_1_mod3/",
        "/home/marco/pool_first/baypass_fin/outc/var_scaf_2_mod3/",
        "/home/marco/pool_first/baypass_fin/outc/var_scaf_3_mod3/",
        "/home/marco/pool_first/baypass_fin/outc/var_scaf_4_mod3/",
        "/home/marco/pool_first/baypass_fin/outc/var_scaf_5_mod3/",
        "/home/marco/pool_first/baypass_fin/outc/var_scaf_6_mod3/",
        "/home/marco/pool_first/baypass_fin/outc/var_scaf_7_mod3/",
        "/home/marco/pool_first/baypass_fin/outc/var_scaf_8_mod3/")
dir_out="/home/marco/pool_first/baypass_fin/outc/merged_mod3/"

system(paste("mkdir ",dir_out))
mode="aux"
merge_baypass2(clirow,scafs,mode,dir_out)
##delete scaffolds merged
for(i in 1:length(scafs)) {
  system(paste("rm -r ",scafs[i],sep=""))
}




##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##new with BF>20
bay2cand_multi=function(snpeff_file,in_files,mode,out,sog=20,sel_scaf=NA,sel_eBPis=T) {
  tab_snpeff=read.table(snpeff_file,stringsAsFactors = F)
  if(!is.na(sel_scaf)) {
    sel=NULL
    for(i in 1:length(sel_scaf)) {
      sel=c(sel,which(tab_snpeff[,1]==sel_scaf[i]))
    }
    
    tab_snpeff=tab_snpeff[sel,]
  }
  
  
  tabs=NULL
  
  for(i in 1:length(in_files)) {
    tabs[[i]]=read.table(in_files[i])  
    
  }
  
  if(dim(tabs[[1]])[1]!=dim(tab_snpeff)[1]) {
    warning("Error different dimension tab and tab_snpEff")
  }
  
  
  if(mode=="aux") {
    
    bf_tot=NULL
    for(i in 1:length(in_files)) {
      bf_tot=c(bf_tot,tabs[[i]][,4])
    }
    
    sel=1:dim(tabs[[1]])[1]
    
    for(i in 1:length(in_files)) {
      sel=intersect(sel,which(tabs[[i]][,4]>=sog))
    }
    bf=sel
    for(i in 1:length(in_files)) {
      bf=cbind(bf,tabs[[i]][sel,4])
    }
  }
  
  if(mode=="bas") {
    
    bf_tot=NULL
    for(i in 1:length(in_files)) {
      bf_tot=c(bf_tot,tabs[[i]][,3])
    }


    
    sel=1:dim(tabs[[1]])[1]
    for(i in 1:length(in_files)) {
      if(sel_eBPis) {
        sel1=which(tabs[[i]][,3]>=sog)
        sel2=which(tabs[[i]][,6]>=3)
        sel=intersect(sel,intersect(sel1,sel2))
      } else {
        sel=intersect(sel,which(tabs[[i]][,3]>=sog))  
      } 
    }
    
    bf=sel
    for(i in 1:length(in_files)) {
      bf=cbind(bf,tabs[[i]][sel,3])
    }
    bf[1:3,]
  }
  
  ff=file(paste(out,"_",sog,"_stat",sep=""),"w")
  
  cat("There are",length(sel)," candidate SNPs. sog:",sog,"\n",file=ff)
  
  write.table(cbind(tab_snpeff[sel,],bf),quote=FALSE,sep="\t",col.names = F, row.names = F, file=paste(out,"_",sog,"_cand",sep=""))
  
  bla=tab_snpeff[sel,]
  
  regio=unique(tab_snpeff[sel,9])
  
  for(i in 1:length(regio)) {
    aa=length(which(bla[,9]==regio[i]))/dim(bla)[1]
    bb=length(which(tab_snpeff[,9]==regio[i]))/dim(tab_snpeff)[1]
    cat("SNPs in ",regio[i],": ",round(aa*100,2),"% (all ",round(bb*100,2),"%)\n",sep="",file=ff)
  }
  
  cat("\n\n",file=ff)
  
  aa=length(which(bla[,8]=="MODIFIER"))/dim(bla)[1]
  bb=length(which(tab_snpeff[,8]=="MODIFIER"))/dim(tab_snpeff)[1]
  cat("SNPs MODIFIER: ",round(aa*100,2),"% (all ",round(bb*100,2),"%)\n",sep="",file=ff)
  
  aa=length(which(bla[,8]=="LOW"))/dim(bla)[1]
  bb=length(which(tab_snpeff[,8]=="LOW"))/dim(tab_snpeff)[1]
  cat("SNPs LOW: ",round(aa*100,2),"% (",round(bb*100,2),"%)\n",sep="",file=ff)
  
  aa=length(which(bla[,8]=="MODERATE"))/dim(bla)[1]
  bb=length(which(tab_snpeff[,8]=="MODERATE"))/dim(tab_snpeff)[1]
  cat("SNPs MODERATE: ",round(aa*100,2),"% (",round(bb*100,2),"%)\n",sep="",file=ff)
  
  aa=length(which(bla[,8]=="HIGH"))/dim(bla)[1]
  bb=length(which(tab_snpeff[,8]=="HIGH"))/dim(tab_snpeff)[1]
  cat("SNPs HIGH: ",round(aa*100,2),"% (",round(bb*100,2),"%)\n",sep="",file=ff)
  
  close(ff)
}


bay2cand_multi2=function(int,snpeff_file,in_files,mode,out,sog=20,sel_scaf=NA,sel_eBPis=T) {
  tab_snpeff=read.table(snpeff_file,stringsAsFactors = F)
  tab_int=read.table(int,stringsAsFactors = F)
  if(!is.na(sel_scaf)) {
    sel=NULL
    for(i in 1:length(sel_scaf)) {
      sel=c(sel,which(tab_snpeff[,1]==sel_scaf[i]))
    }
    
    tab_snpeff=tab_snpeff[sel,]
  }
  
  
  tabs=NULL
  
  for(i in 1:length(in_files)) {
    tabs[[i]]=read.table(in_files[i])  
    
  }
  
  if(dim(tabs[[1]])[1]!=dim(tab_snpeff)[1]) {
    warning("Error different dimension tab and tab_snpEff")
  }
  
  
  if(mode=="aux") {
    
    bf_tot=NULL
    for(i in 1:length(in_files)) {
      bf_tot=c(bf_tot,tabs[[i]][,4])
    }
    
    sel=1:dim(tabs[[1]])[1]
    
    for(i in 1:length(in_files)) {
      sel=intersect(sel,which(tabs[[i]][,4]>=sog))
    }
    bf=sel
    for(i in 1:length(in_files)) {
      bf=cbind(bf,tabs[[i]][sel,4])
    }
  }
  
  if(mode=="bas") {
    
    bf_tot=NULL
    for(i in 1:length(in_files)) {
      bf_tot=c(bf_tot,tabs[[i]][,3])
    }
    
    
    
    sel=1:dim(tabs[[1]])[1]
    for(i in 1:length(in_files)) {
      if(sel_eBPis) {
        sel1=which(tabs[[i]][,3]>=sog)
        sel2=which(tabs[[i]][,6]>=3)
        sel=intersect(sel,intersect(sel1,sel2))
      } else {
        sel=intersect(sel,which(tabs[[i]][,3]>=sog))  
      } 
    }
    
    bf=sel
    for(i in 1:length(in_files)) {
      bf=cbind(bf,tabs[[i]][sel,3])
    }
    bf[1:3,]
  }
  
  ff=file(paste(out,"_",sog,"_stat",sep=""),"w")
  
  cat("There are",length(sel)," candidate SNPs. sog:",sog,"\n",file=ff)
  
  write.table(cbind(tab_snpeff[sel,],bf),quote=FALSE,sep="\t",col.names = F, row.names = F, file=paste(out,"_",sog,"_cand",sep=""))
  
  bla=tab_snpeff[sel,]
  
  regio=unique(tab_snpeff[sel,9])
  
  for(i in 1:length(regio)) {
    aa=length(which(bla[,9]==regio[i]))/dim(bla)[1]
    bb=length(which(tab_snpeff[,9]==regio[i]))/dim(tab_snpeff)[1]
    cat("SNPs in ",regio[i],": ",round(aa*100,2),"% (all ",round(bb*100,2),"%)\n",sep="",file=ff)
  }
  
  cat("\n\n",file=ff)
  
  aa=length(which(bla[,8]=="MODIFIER"))/dim(bla)[1]
  bb=length(which(tab_snpeff[,8]=="MODIFIER"))/dim(tab_snpeff)[1]
  cat("SNPs MODIFIER: ",round(aa*100,2),"% (all ",round(bb*100,2),"%)\n",sep="",file=ff)
  
  aa=length(which(bla[,8]=="LOW"))/dim(bla)[1]
  bb=length(which(tab_snpeff[,8]=="LOW"))/dim(tab_snpeff)[1]
  cat("SNPs LOW: ",round(aa*100,2),"% (",round(bb*100,2),"%)\n",sep="",file=ff)
  
  aa=length(which(bla[,8]=="MODERATE"))/dim(bla)[1]
  bb=length(which(tab_snpeff[,8]=="MODERATE"))/dim(tab_snpeff)[1]
  cat("SNPs MODERATE: ",round(aa*100,2),"% (",round(bb*100,2),"%)\n",sep="",file=ff)
  
  aa=length(which(bla[,8]=="HIGH"))/dim(bla)[1]
  bb=length(which(tab_snpeff[,8]=="HIGH"))/dim(tab_snpeff)[1]
  cat("SNPs HIGH: ",round(aa*100,2),"% (",round(bb*100,2),"%)\n",sep="",file=ff)
  
  int2=intersect(paste(bla[,1],bla[,2],sep=""),paste(tab_int[,1],tab_int[,2],sep=""))
  cat("SNPs intergenic (2000 bp): ",round((length(int2)/dim(bla)[1])*100,2),"% (",round((dim(tab_int)[1]/dim(tab_snpeff)[1])*100,2),"%)\n",sep="",file=ff)
  
  close(ff)
}
##########################################################################################################################################
##########################################################################################################################################


snpeff_file="/home/marco/pool_first/merge_snp4/outc_05_21.snpEff"

dir_out="/home/marco/pool_first/baypass_fin/outc/cand/"
system(paste("mkdir ",dir_out,sep=""))

clirow="/home/marco/pool_first/baypass_fin/outc/cli_row"
clirow=readLines(clirow)


i=1
sog=10
for(i in 1:length(clirow)) {
  
  mode="aux"
  
  in_files=c("/home/marco/pool_first/baypass_fin/outc/merged_mod1/","/home/marco/pool_first/baypass_fin/outc/merged_mod2/","/home/marco/pool_first/baypass_fin/outc/merged_mod3/")
  in_files=paste(in_files,clirow[i],"_",mode,sep="")
  out=paste(dir_out,clirow[i],"_",mode,sep="")
  #bay2cand_multi(snpeff_file,in_files,mode,out)
  bay2cand_multi(snpeff_file,in_files,mode,out,sog=sog)
  
  mode="bas"
  in_files=c("/home/marco/pool_first/baypass_fin/outc/merged_mod1/","/home/marco/pool_first/baypass_fin/outc/merged_mod2/","/home/marco/pool_first/baypass_fin/outc/merged_mod3/")
  in_files=paste(in_files,clirow[i],"_",mode,sep="")
  out=paste(dir_out,clirow[i],"_",mode,sep="")
  #bay2cand_multi(snpeff_file,in_files,mode,out)
  bay2cand_multi(snpeff_file,in_files,mode,out,sog=sog)
  
  
  cat(i,"\n")
}


##############
##only bas more variables

snpeff_file="/home/marco/pool_first/merge_snp4/outc_05_21.snpEff"

dir_out="/home/marco/pool_first/baypass_fin/outc2/cand/"
system(paste("mkdir ",dir_out,sep=""))
clirow="/home/marco/pool_first/baypass_fin/outc2/cli_row"
clirow=readLines(clirow)

i=1
sog=10
sog=20
sog=15

for(i in 1:length(clirow)) {
  
  mode="bas"
  in_files=c("/home/marco/pool_first/baypass_fin/outc2/merged_mod1/","/home/marco/pool_first/baypass_fin/outc2/merged_mod2/","/home/marco/pool_first/baypass_fin/outc2/merged_mod3/")
  in_files=paste(in_files,clirow[i],"_",mode,sep="")
  out=paste(dir_out,clirow[i],"_",mode,sep="")
  #bay2cand_multi(snpeff_file,in_files,mode,out)
  bay2cand_multi(snpeff_file,in_files,mode,out,sog=sog)
  
  
  cat(i,"\n")
}


#########################
##Aux 5 variable

snpeff_file="/home/marco/pool_first/merge_snp4/outc_05_21.snpEff"
int="/home/marco/pool_first/merge_snp4/outc_05_21_int_2000bp"
dir_out="/home/marco/pool_first/baypass_fin/outc/cand/"
system(paste("mkdir ",dir_out,sep=""))

clirow="/home/marco/pool_first/baypass_fin/outc/cli_row"
clirow=readLines(clirow)


i=1

for(i in 1:length(clirow)) {
  
  mode="aux"

  in_files=c("/home/marco/pool_first/baypass_fin/outc/merged_mod1/","/home/marco/pool_first/baypass_fin/outc/merged_mod2/","/home/marco/pool_first/baypass_fin/outc/merged_mod3/")
  in_files=paste(in_files,clirow[i],"_",mode,sep="")
  out=paste(dir_out,clirow[i],"_",mode,sep="")
  
  sog=10
  bay2cand_multi2(int,snpeff_file,in_files,mode,out,sog=sog)
  
  sog=20
  bay2cand_multi2(int,snpeff_file,in_files,mode,out,sog=sog)
  
  cat(i,"\n")
}

##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

##funziona :) 
run_snp2go=function(pos_cand,pos_nocand,out,gtf_file,go_file,gene_file,goall_file,orto_file,bp_ab=50,nruns=100000) {
  
  require(SNP2GO)
  pos_cand[,1]=gsub("scaffold_","",pos_cand[,1])
  pos_nocand[,1]=gsub("scaffold_","",pos_nocand[,1])
  
  cand=GRanges(seqnames=pos_cand[,1],ranges=IRanges(pos_cand[,2],pos_cand[,2]))
  nocand=GRanges(seqnames=pos_nocand[,1],ranges=IRanges(pos_nocand[,2],pos_nocand[,2]))
  
  
  
  cat("found",length(cand),"SNP candidates (",length(nocand),")\n")
  # Case 2: Using a GTF file + gene association file
  #nruns=1000
  y = snp2go(gtf=gtf_file,
             goFile=go_file,
             candidateSNPs=cand,
             noncandidateSNPs=nocand,
             FDR=0.05, 
             runs=nruns,
             extension=bp_ab)
  
  gff.significant.terms <- y$enriched$GO
  ris=cbind(as.character(gff.significant.terms),y$enriched$FDR)
  
  
  
  goall=read.table(goall_file,sep="\t",stringsAsFactors = FALSE)
  rownames(goall)=goall[,1]
  
  
  tab_ris=cbind(y$enriched,desc=goall[as.character(y$enriched$GO),2])
  
  
  if(dim(ris)[1]!=0) {
    write.table(ris,quote=FALSE,row.names = FALSE,col.names=FALSE,file=paste(out,"_snp2go_revigo",sep=""))    
    write.table(goall[ris[,1],],quote=TRUE,row.names = FALSE,col.names=FALSE,file=paste(out,"_snp2go_goterm",sep=""))
    tab_ris=tab_ris[order(tab_ris$FDR),]
    write.table(tab_ris,quote=TRUE,row.names = FALSE,col.names=TRUE,file=paste(out,"_snp2go_enrich",sep=""))
  }
  
  
  
  ##get all GO term associated with gene in SNP correlated
  go_t=read.table(go_file,stringsAsFactors = FALSE,sep="\t",header = TRUE)
  gene_t=read.table(gene_file,sep="\t",stringsAsFactors = FALSE)
  rownames(gene_t)=gene_t[,5]
  
  gene_t[,5]=0
  
  
  sel_tot=NULL
  go_all=NULL
  
  i=1
  for(i in 1:dim(pos_cand)[1]) {
    sel=which(gene_t[,1]==pos_cand[i,1])
    tab_p=gene_t[sel,]
    
    sel2=which(tab_p[,2]<=pos_cand[i,2])
    sel3=which(tab_p[,3]>=pos_cand[i,2])
    sel23=intersect(sel2,sel3)
    if(length(sel23)!=0) {
      gene_t[rownames(tab_p)[sel23],5]=gene_t[rownames(tab_p)[sel23],5]+1
      for(z in 1:length(sel23)) {
        sel_tot=c(sel_tot,rownames(tab_p)[sel23[z]])
        go_a=go_t[which(go_t[,1]==rownames(tab_p)[sel23[z]]),2]
        #go_a=cbind(go_a,rep(tot[sel_cand[i]],length(go_a)))
        go_all=c(go_all,go_a)   
      }
    }
  }
  
  go_all=go_all[which(go_all!="")]
  
  cat(go_all,sep="\n",file=paste(out,"_all_GO",sep=""))
  
  gene_t=gene_t[which(gene_t[,5]!=0),]
  
  gene_t=gene_t[order(gene_t[,5],decreasing = T),]
  
  
  if(dim(gene_t)[1]!=0) {
    orto=read.table(orto_file,fill = T,header = T,stringsAsFactors = F)
    atha=NULL
    i=1
    for(i in 1:dim(gene_t)[1]) {
      atha=c(atha,paste(unique(orto[which(orto[,1]==rownames(gene_t)[i]),2]),collapse=","))
    }
    gene_t=cbind(gene_t,atha)
    
    colnames(gene_t)=c("scaf","start","end","strand","nSNPs","thaliana")
    write.table(gene_t,quote=F,file=paste(out,"_genes",sep=""))
    
  }
  
  
  
  #   n_per=NULL
  #   uni_go_all=unique(go_all)
  #   
  #   #remove the three main GO CC MF BP
  #   uni_go_all=setdiff(uni_go_all,c("GO:0003674","GO:0005575","GO:0008150"))
  #   for(i in 1:length(uni_go_all)) {
  #     #     offspring = c(uni_go_all[i],GOMFOFFSPRING[[uni_go_all[i]]],GOBPOFFSPRING[[uni_go_all[i]]],GOCCOFFSPRING[[uni_go_all[i]]])
  #     #     nn=0
  #     #     for(z in 1:length(offspring)) {
  #     #       nn=nn+length(which(GO_c==offspring[z]))
  #     #     }
  #     nn=length(which(go_all==uni_go_all[i]))
  #     n_per=c(n_per,nn)
  #   }
  #   
  
  
  #   require(wordcloud)
  #   
  #   png(file=paste(out,"_cloud_all.png",sep=""), width=5000, height=5000) 
  #   wordcloud(goall[uni_go_all,2],n_per,
  #             scale=c(5,0.5), max.words=300, random.order=FALSE, min.freq = 5,
  #             rot.per=0, use.r.layout=FALSE, colors=brewer.pal(8, "Dark2"))
  #   dev.off()  
  #   
  #   if(dim(tab_ris)[1]!=0) {
  #     png(file=paste(out,"_cloud_snp2go_FDR.png",sep=""), width=5000, height=5000) 
  #     wordcloud( tab_ris$desc, -scale(tab_ris$FDR),
  #                scale=c(5,0.5), max.words=100, random.order=FALSE,
  #                rot.per=0, use.r.layout=FALSE, colors=brewer.pal(8, "Dark2"))
  #     dev.off()
  # 
  #     png(file=paste(out,"_cloud_snp2go_gG.png",sep=""), width=5000, height=5000) 
  #     wordcloud(tab_ris$desc, scale(tab_ris$g/tab_ris$G),
  #               scale=c(5,0.5), max.words=100, random.order=FALSE,
  #               rot.per=0, use.r.layout=FALSE, colors=brewer.pal(8, "Dark2"))
  #     dev.off()
  
  #     png(file=paste(out,"_cloud_snp2go_g.png",sep=""), width=5000, height=5000) 
  #     wordcloud(tab_ris$desc, tab_ris$g,
  #               scale=c(5,0.5), max.words=100, random.order=FALSE,
  #               rot.per=0, use.r.layout=FALSE, colors=brewer.pal(8, "Dark2"))
  #     dev.off()
  
  
  
  #     png(file=paste(out,"_cloud_snp2go_nc.png",sep=""), width=5000, height=5000) 
  #     wordcloud(tab_ris$desc, tab_ris$nc,
  #               scale=c(5,0.5), max.words=100, random.order=FALSE,
  #               rot.per=0, use.r.layout=FALSE, colors=brewer.pal(8, "Dark2"))
  #     dev.off()
  
  
}


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################


gtf_file="/home/marco/pool_first/Araly1_Genome/version2/Alyr_V2.gtf"
go_file="/home/marco/pool_first/Araly1_Genome/version2/mart_Alyrata_V2_BP"
goall_file="/home/marco/pool_first/annotation/goterms.txt"
orto_file="/home/marco/pool_first/Araly1_Genome/version2/Alyr_Atha_ortholog2"
snpeff_file="/home/marco/pool_first/merge_snp4/outc_05_21.snpEff"
dir_in="/home/marco/pool_first/baypass_fin/outc/cand/"
clirow="/home/marco/pool_first/baypass_fin/outc/cli_row"

bp_ap=50
gene_file="/home/marco/pool_first/Araly1_Genome/version2/Alyr_V2_genes_50bp"
dir_out="/home/marco/pool_first/baypass_fin/outc/snp2go/"

# bp_ap=500
# gene_file="/home/marco/pool_first/Araly1_Genome/version2/Alyr_V2_genes_500bp"
# dir_out="/home/marco/pool_first/baypass_fin/outc/snp2go_500/"
# bp_ap=1000
# gene_file="/home/marco/pool_first/Araly1_Genome/version2/Alyr_V2_genes_1000bp"
# dir_out="/home/marco/pool_first/baypass_fin/outc/snp2go_1000/"

# dir_in="/home/marco/pool_first/baypass_fin/outc2/cand/"
# clirow="/home/marco/pool_first/baypass_fin/outc2/cli_row"
# dir_out="/home/marco/pool_first/baypass_fin/outc2/snp2go/"

go_file="/home/marco/pool_first/Araly1_Genome/version2/mart_Alyrata_V2_BP_GOSLIM"
bp_ap=50
gene_file="/home/marco/pool_first/Araly1_Genome/version2/Alyr_V2_genes_50bp"
dir_out="/home/marco/pool_first/baypass_fin/outc/snp2go_goslim/"


go_file="/home/marco/pool_first/Araly1_Genome/version2/mart_Alyrata_V2_BP_161014"
bp_ap=50
gene_file="/home/marco/pool_first/Araly1_Genome/version2/Alyr_V2_genes_50bp"
dir_out="/home/marco/pool_first/baypass_fin/outc/snp2go_161014/"


sog=20
#sog=10

######
system(paste("mkdir ",dir_out,sep=""))
clirow=readLines(clirow)

tab_snpeff=read.table(snpeff_file,stringsAsFactors = F)
tab_snpeff=tab_snpeff[,1:2]
rownames(tab_snpeff)=paste(tab_snpeff[,1],tab_snpeff[,2],sep="-")

for(i in 1:length(clirow)) {
  mode="aux"
  pos_cand=read.table(paste(dir_in,clirow[i],"_",mode,"_",sog,"_cand",sep=""))
  sel=setdiff(rownames(tab_snpeff),paste(pos_cand[,1],pos_cand[,2],sep="-"))
  pos_nocand=tab_snpeff[sel,]
  out=paste(dir_out,clirow[i],"_",mode,"_",sog,sep="")
  run_snp2go(pos_cand,pos_nocand,out,gtf_file,go_file,gene_file,goall_file,orto_file,bp_ap)
  
#   mode="bas"
#   pos_cand=read.table(paste(dir_in,clirow[i],"_",mode,"_",sog,"_cand",sep=""))
#   sel=setdiff(rownames(tab_snpeff),paste(pos_cand[,1],pos_cand[,2],sep="-"))
#   pos_nocand=tab_snpeff[sel,]
#   out=paste(dir_out,clirow[i],"_",mode,"_",sog,sep="")
#   run_snp2go(pos_cand,pos_nocand,out,gtf_file,go_file,gene_file,goall_file,orto_file,bp_ap)
#   
}

###############################################################################
############################
#########merge gene
##create a table with the merged genes few between var and also with fischer and hancock.

##INPUT


clirow=c("Sub", "Tmin_ESp", "bio16", "bio4", "alpha")
clirow=c("Sub", "Tmin_ESp", "bio16", "alpha")
dir_in="/home/marco/Dropbox/dottorato_Neuchatel/baypass_aux/snp2go/"

suff="_aux_20_genes"
nome_out="/home/marco/Dropbox/dottorato_Neuchatel/hancock_fisher/aux_20"
###

#############


i=1
#aa=cbind(rownames(aa),aa[,1:4],aa[,6],aa[,5])
#colnames(aa)=c("chiave", "scaf","start","end","strand","thaliana",paste(clirow[i],"_nSNP",sep=""))

aa=read.table(paste(dir_in,clirow[i],suff,sep=""),header = T, fill=T,stringsAsFactors = F)
aa[,6]=gsub("^$","no",aa[,6])
chiave=paste(rownames(aa),aa[,1],aa[,2],aa[,3],aa[,4],aa[,6],sep="_")
aa=cbind(chiave,aa[,5])
colnames(aa)=c("chiave",clirow[i])
aa_tot=aa
aa[1:3,]

for(i in 2:length(clirow)) {
  aa=read.table(paste(dir_in,clirow[i],suff,sep=""),header = T, fill=T,stringsAsFactors = F)
  aa[,6]=gsub("^$","no",aa[,6])
  chiave=paste(rownames(aa),aa[,1],aa[,2],aa[,3],aa[,4],aa[,6],sep="_")
  aa=cbind(chiave,aa[,5])
  colnames(aa)=c("chiave",clirow[i])
  
  aa_tot=merge(x = aa_tot, y = aa, by = "chiave", all = TRUE)
}


aa_tot=as.matrix(aa_tot)
aa_tot[is.na(as.matrix(aa_tot))] <- 0

coli=colnames(aa_tot[,-1])
write.table(aa_tot,quote=F,sep="\t",col.names = F, row.names = F,file=paste(nome_out,"_temp",sep=""))
###########3
############
##a mano cambia _ in \t
###########
###########
aa_tot=read.table(paste(nome_out,"_temp",sep=""))

colnames(aa_tot)=c("ID_ly", "scaf","start","end","strand","ID_th",coli)
write.table(aa_tot,quote=F,sep="\t",file=nome_out)
# somma=apply(aa_tot[,7:dim(aa_tot)[2]],1,sum)
# bah=aa_tot[order(somma,decreasing = T),]





