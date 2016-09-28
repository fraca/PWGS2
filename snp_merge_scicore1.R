#########################################################################################################################
#########################################################################################################################
#####################
###nuovo
#crea input x baypass con missing data
#rimuove over max cov, N in ref, low that freq e triallelici nella pop e tra pops
###da fare
##parte lunga parallelizza.

###prova
# folder_in="/scicore/home/williy/fracasse/snp_BED_fin/"
# mod="tot"
# scafs=c("scaffold_2", "scaffold_8")
# name_pops=c("p07C", "p07D", "p11L_F0", "p07E","p07F","p14C","pHa21")
# min_freq=1.5
# min_pop=1
# nome_out="/scicore/home/williy/fracasse/snp_BED_fin/pro_p1_f2"

###ok
folder_in="/scicore/home/williy/fracasse/snp_BED_fin/"
mod="tot"
name_pops=c("p07C","p07D","p07E","p07F","p07G","p07J","p07K","p07L","p07M","p07N","p07O","p07P","p07Q","p07R","p11C","p11AA","p11AB","p11AC","p11AE","p11AG","p11AH","p11AJ","p11B","p11G","p11H","p11J","p11K","p11L_F0","p11M","p11N","p11O","p11P","p11Q","p11S","p11T","p11U","p11W","p11X","p11Z","p11D","p11V","p11R","p11A","p11E","p11F","p11Y","p14A","p14B","p14C","p14D","p14E","p07H","pHa31")
scafs=c("scaffold_1", "scaffold_2", "scaffold_3", "scaffold_4", "scaffold_5", "scaffold_6", "scaffold_7", "scaffold_8")
#scafs=c("scaffold_6", "scaffold_7", "scaffold_8")
min_freq=1.5
min_pop=1
nome_out="/scicore/home/williy/fracasse/snp_BED_fin/tot_p1_f15"

folder_in="/scicore/home/williy/fracasse/snp_BED_fin/"
#name_pops=c("p07C","p07D","pHa31")
#scafs=c("scaffold_6", "scaffold_7", "scaffold_8")
name_pops=c("p07C","p07D","p07E","p07F","p07G","p07J","p07K","p07L","p07M","p07N","p07O","p07P","p07Q","p07R","p11C","p11AA","p11AB","p11AC","p11AE","p11AG","p11AH","p11AJ","p11B","p11G","p11H","p11J","p11K","p11L_F0","p11M","p11N","p11O","p11P","p11Q","p11S","p11T","p11U","p11W","p11X","p11Z","p11D","p11V","p11R","p11A","p11E","p11F","p11Y","p14A","p14B","p14C","p14D","p14E","p07H","pHa31")
scafs=c("scaffold_1", "scaffold_2", "scaffold_3", "scaffold_4", "scaffold_5", "scaffold_6", "scaffold_7", "scaffold_8")
mod="intergenic_1000bp"
min_freq=1.5
min_pop=1
#nome_out="/scicore/home/williy/fracasse/snp_BED_fin/pro"
nome_out="/scicore/home/williy/fracasse/snp_BED_fin/int_p53_f15"



###

name_temp=paste(folder_in,"temp.BED",sep="")
temp_ri=paste(folder_in,"SNP_ri",sep="")
z=1

for(z in 1:length(scafs)) {

  
  scaf=scafs[z]  
  i=1
  
  input=paste(folder_in,name_pops[i],"_",mod,"_SNP_scaf/",scaf,sep="")
  tab=read.table(paste(input,".varscan",sep=""),stringsAsFactors=FALSE) ##vedi se cambiare varscan2 nome

  freq=as.numeric(gsub("%","",tab[,8]))
  rem_fr=which(freq<min_freq)
  if(length(rem_fr)!=0) {
    tab=tab[-rem_fr,]
  }
  
  max_cov=tab[,6]+tab[,7]
  mean_cov=round(mean(max_cov))

  blo=8-nchar(tab[,2])
  blo=unlist(lapply(blo,function(x) paste(rep("0",x),collapse="")))
  tab[,2]=paste(blo,tab[,2],sep="")
  
  col1=paste(tab[,1],tab[,2],tab[,4],tab[,20],sep="-")
  #a_alt=tab[,20]
  tab2=cbind(col1,tab[,6:7])
  colnames(tab2)=c("chiave",paste("R_",name_pops[i],sep=""),paste("A_",name_pops[i],sep=""))
  tab_tot=tab2
  mean_cov_tot=mean_cov
  cat("Merging",scaf,"...\n")
  i=2
  for(i in 2:length(name_pops)) {
    
    input=paste(folder_in,name_pops[i],"_",mod,"_SNP_scaf/",scaf,sep="")
    tab=read.table(paste(input,".varscan",sep=""),stringsAsFactors=FALSE) ##vedi se cambiare varscan2 nome
    
    freq=as.numeric(gsub("%","",tab[,8]))
    rem_fr=which(freq<min_freq)
    if(length(rem_fr)!=0) {
      tab=tab[-rem_fr,]
    }
    
    max_cov=tab[,6]+tab[,7]
    mean_cov=round(mean(max_cov))
    
    blo=8-nchar(tab[,2])
    blo=unlist(lapply(blo,function(x) paste(rep("0",x),collapse="")))
    tab[,2]=paste(blo,tab[,2],sep="")
    
    col1=paste(tab[,1],tab[,2],tab[,4],tab[,20],sep="-")
    a_alt=tab[,20]
    tab2=cbind(col1,tab[,6:7])
    colnames(tab2)=c("chiave",paste("R_",name_pops[i],sep=""),paste("A_",name_pops[i],sep=""))
    tab_tot=merge(x = tab_tot, y = tab2, by = "chiave", all = TRUE)
    mean_cov_tot=c(mean_cov_tot,mean_cov)
    cat(i,"\t")
  }
  names(mean_cov_tot)=name_pops
  
  
  ##############
  ##sel alt allele
  
  n_NA=apply(tab_tot[,grep("A_",colnames(tab_tot))],1,function(x) length(which(!is.na(x))))
  sel_alt=which(n_NA>=min_pop)
  
  ff=file(paste(nome_out,"_",scaf,"_stat2",sep=""),"w")
  cat("SNP tot before:",dim(tab_tot)[1],"\n",file=ff)
  cat("polimorphic in more than",min_pop,"pop:",length(sel_alt),"\n",file=ff)
  
  tab_tot=tab_tot[sel_alt,]
  
  
  ###ordina tabella
  bla=order(as.character(tab_tot[,1]))
  tab_tot=tab_tot[bla,]
  
  pos_tot=gsub(scaf,"",as.character(tab_tot[,1]))
  pos_tot=as.integer(gsub("[-ACTGN]","",pos_tot))
  ##in teoria con nuovo N non dovrebbe esserci
  #########
  ##rimuovo triallelici tra pops
  #length(unique(pos_tot))
  
  dup=which(duplicated(pos_tot) | duplicated(pos_tot[length(pos_tot):1])[length(pos_tot):1])
  if(length(dup)!=0) {
    tab_dup=tab_tot[dup,]
    #pos_dup=pos_tot[dup]
    
    tab_tot=tab_tot[-dup,]
    pos_tot=pos_tot[-dup]  
  } else {
    tab_dup=data.frame()
  }
  
  
  #after remove triallilic add rownames

  rownames(tab_tot)=paste("p",pos_tot,sep="")
  
  cat("total SNPs:",dim(tab_tot)[1],"\n",file=ff)
  cat("trialleic SNPs (between pops):",dim(tab_dup)[1],"\n",file=ff)
  close(ff)
  ####
  cat("wait...\n")
  

  for(i in 1:length(name_pops)) {
    
    input=paste(folder_in,name_pops[i],"_",mod,"_SNP_scaf/",scaf,"_filt2.BED",sep="")
    col_sel=grep(name_pops[i],colnames(tab_tot))
    
    #tab_tot
    a_NA=which(is.na(tab_tot[,col_sel[1]]))
    
    if(length(a_NA)!=0) {
      ##calling bedtools
      #tab_t=cbind(rep(scaf,length(a_NA)),as.integer(pos_tot[a_NA]),as.integer(pos_tot[a_NA]))
      tab_t=cbind(rep(scaf,length(a_NA)),pos_tot[a_NA],pos_tot[a_NA])
      
      write.table(tab_t,quote=F,row.names=F,col.names=F,sep="\t",file=name_temp)
      ##errore in file grande aggiunge NA in tab_t
      system(paste("bedtools sort -i ",name_temp," | bedtools intersect -a stdin -b ",input," -f 1 | awk \'{print $2}\' > ",temp_ri,sep=""))
      snp_ok=paste("p",readLines(temp_ri),sep="")
      
      ##problema snp_ok non presenti in tab_tot mettere as.integer in pos_tot
      
      if(length(snp_ok)!=0) {
        tab_tot[snp_ok,col_sel[1]]=mean_cov_tot[name_pops[i]]
        tab_tot[snp_ok,col_sel[2]]=0
      }
      
      if(dim(tab_tot)[1]!=length(pos_tot)) {
        save.image(file=paste(nome_out,".Rdata",sep=""))
        stop("Errore tab_tot cambia")
      }
    }

    
    

    
    cat(name_pops[i],"\t")
    
    
  }
  
  write.table(tab_tot,file=paste(nome_out,"_",scaf,sep=""),quote=F,sep=" ",row.names = F,col.names = T)
  write.table(tab_dup,file=paste(nome_out,"_",scaf,"_tripop",sep=""),quote=F,sep=" ",row.names = F,col.names = T)
  cat(mean_cov_tot,sep="\n",file=paste(nome_out,"_",scaf,"_meancov",sep=""))
  cat(scaf,"done\n")
  
  
}

#save.image(file=paste(nome_out,".Rdata",sep=""))
system(paste("rm ",name_temp,sep=""))
system(paste("rm ",temp_ri,sep=""))


#########################################################################################################################
