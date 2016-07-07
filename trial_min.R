#############################################################

#Selection of the SNP, remove triallelic.
#after merge_SNPs_pops.
#remove triallelic or more
#polymorphic at least in 2 (or more) populations
#better looking of the table

##INPUT

# #input tab from merge_SNPs_pops
# # tab_in="/home/marco/pool_first/SNP_database/prova3_tab_snape"
# tab_in="/home/marco/pool_first/SNP_database/p52_tab_snape"
# # # SNP polymorphic at least in n pops
# min_pop=2
# # 
# min_freq=0.05
# # # for now not done maybe later (to implement)#freq min to consider SNPs for triallelic:
# # # min_freq=0.05
# # righe="/home/marco/pool_first/SNP_database/prova3/SNP_tot_snape"
# righe="/home/marco/pool_first/SNP_database/p52/SNP_tot_snape"
# 
# 
# # #output file
# # nome_out="/home/marco/pool_first/SNP_database/prova3_tab_snape_clean"
# nome_out="/home/marco/pool_first/SNP_database/p52_tab_snape_clean"


######################################

trial_min=function(tab_in,min_pop,min_freq,righe,nome_out) {
#   tab_in="/tmp/tab_var"
#   min_freq=0.05
#   righe="/tmp/snp_var"
#   min_pop=2
#   nome_out="/tmp/out"
  
  ff_fixed=file(paste(nome_out,"_rfixed",sep=""),"w")
  ff_min_freq=file(paste(nome_out,"_rmin_freq",sep=""),"w")
  ff_min_pop=file(paste(nome_out,"_rmin_pop",sep=""),"w")
  #close(ff_fixed)
  
  system(paste("rm ", nome_out,"_rtri",sep=""))
  system(paste("rm ", nome_out,".BED",sep=""))
  
  
  righe=readLines(righe)
  colonne=readLines(paste(tab_in,"_col",sep=""))

  
  
  
  #da rifare
  #mem_limit=100000
  mem_limit=50000
  #mem_limit=5000
  j=1
  j_tot=ceiling(length(righe)/mem_limit)
  #cat(j_tot,"\n")  
  c1=1
  c2=mem_limit
  dim_dopo=0
  n_min_freq=0
  n_min_pop=0
  n_tri=0
  n_fixed=0

  alt=colonne[grep("ALT",colonne)]
  freq=setdiff(colonne,alt)

    
  for(j in 1:j_tot) {
    cat(j,"of",j_tot,"\n")
    
    #create 2 table
    aa=read.table(tab_in,skip=c1-1,nrows=mem_limit, stringsAsFactors = FALSE)
    colnames(aa)=colonne
    rownames(aa)=righe[c1:c2]
    bb=round(aa[,freq],4)
    aa=aa[,alt]   
    

    ###remove if freq is less than min_freq in all samples
    bb_b=apply(bb,1,function(x) x[!is.na(x)])
    
    bb_d=lapply(bb_b,function(x) length(which(x>=min_freq)))
    rem=which(bb_d==0)
    
    if(length(rem)!=0) {
      cat(rownames(aa)[rem],sep="\n",file=ff_min_freq)
      n_min_freq=n_min_freq+length(rem)
      aa=aa[-rem,]
      bb=bb[-rem,]
      
    }
    
    
    # put like 0 the freq in triallelic that are lower than 0.05 (does not matter for 0.95 beacuse I don't remove anyway.)
    # if number of alleles are the same put all to 0.05
    aa_b=apply(aa,1,function(x) unique(x[!is.na(x)]))
    aa_c=lapply(aa_b, length)
    sel_tri=which(aa_c>=2)

    
    for(i in 1:length(sel_tri)) {
      alle_ok=NULL
      for(h in 1:length(aa_b[[sel_tri[i]]])) {
        alle_ok=c(alle_ok,length(which(aa[sel_tri[i],]==aa_b[[sel_tri[i]]][h])))
      }
      alle_ok=which(alle_ok==max(alle_ok))
      
      sel_bb=which(bb[sel_tri[i],]<=min_freq)
      if(length(alle_ok)==1) {
        sel_bb=intersect(sel_bb,which(aa[sel_tri[i],]!=aa_b[[sel_tri[i]]][alle_ok]))
        #cat("entro",i,"\t")
      }
      aa[sel_tri[i],sel_bb]=NA
      bb[sel_tri[i],sel_bb]=NA
      #cat(sel_bb,"\n")
    }
    
    #remove triallelic
    aa_b=apply(aa,1,function(x) x[!is.na(x)])
    aa_b=lapply(aa_b, unique)
    aa_c=lapply(aa_b, length)
    sel_tri=which(aa_c>1)
    if(length(sel_tri)!=0) {
      tab_tri=rownames(aa)[sel_tri]
      for(i in 1:dim(aa)[2]) {
        tab_tri=cbind(tab_tri,aa[sel_tri,i],bb[sel_tri,i])
      }
      n_tri=n_tri+length(sel_tri)
      write.table(tab_tri,append=TRUE,col.names=FALSE,row.names=FALSE,file=paste(nome_out,"_rtri",sep=""),quote=FALSE)  
      aa=aa[-sel_tri,]
      bb=bb[-sel_tri,]
    }
    

    # toglie all fixed or NA
    sel_fix=which(apply(bb,1,function(x) length(unique(x)))==1)    
       
    if(length(sel_fix)!=0) {
      cat(rownames(aa)[sel_fix],sep="\n",file=ff_fixed)
      n_fixed=n_fixed+length(sel_fix)
      aa=aa[-sel_fix,]
      bb=bb[-sel_fix,]
      
    }

    
    # remove if SNP different from REF in less than min_pop
    #min_pop=2
    #min_pop=10
  
    aa_SNP=apply(aa,1,function(x) length(x[!is.na(x)]))

    sel_min_pop=which(aa_SNP<min_pop)
    
#     cat(min_pop,"\n")
#     cat(length(which(aa_SNP<2)),"2SNP\n")
#     cat(length(which(aa_SNP<10)),"10SNP\n")
#     dim_be=dim(aa)[1]
    
    if(length(sel_min_pop)!=0) {
      cat(rownames(aa)[sel_min_pop],sep="\n",file=ff_min_pop)
      n_min_pop=n_min_pop+length(sel_min_pop)
      aa=aa[-sel_min_pop,]
      bb=bb[-sel_min_pop,]
     
    }
    #cat(dim_be-dim(aa)[1],"==",length(sel_min_pop),"\n")
    
    #create like BED files
    aa_b=apply(aa,1,function(x) x[!is.na(x)])
    aa_b=unlist(lapply(aa_b, unique))

    perbed=strsplit(rownames(aa),"-")
    scaf=NULL
    pos=NULL
    ref=NULL
    
    for(i in 1:length(perbed)) {
      scaf=c(scaf,perbed[[i]][1])
      pos=c(pos,perbed[[i]][2])
      ref=c(ref,perbed[[i]][3])
    }

    for(i in 1:dim(bb)[2]) { #change NA to 0
      bb[which(is.na(bb[,i])),i]=0  
    }
    
    
    ris=cbind(SCAF=scaf,POS1=pos,POS2=pos,REF=ref,ALT=aa_b,bb)    
    write.table(ris,append=TRUE,col.names=FALSE,row.names=FALSE, file=paste(nome_out,".BED",sep=""),sep="\t",quote=FALSE)
    
    c1=c1+mem_limit
    c2=c2+mem_limit
    if(c2>length(righe))
      c2=length(righe) 
    dim_dopo=dim_dopo+dim(ris)[1]
  }
  
  cat(colnames(ris),sep="\n",file=paste(nome_out,"_col",sep=""))
  
  ff_stat=file(paste(nome_out,"stat",sep="_"),"w")
  
  cat("Before",length(righe),"SNPs\n",file=ff_stat)
  cat("All SNP freq less than",min_freq,"in all pop",n_min_freq,"\n",file=ff_stat)
  cat("Triallelic",n_tri,"\n",file=ff_stat)
  cat("Fixed or frequency less than 0.05",n_fixed,"\n",file=ff_stat)
  cat("SNP different from REF in less than",min_pop,"pop",n_min_pop,"\n",file=ff_stat)
  cat("After",dim_dopo,"SNPs\n",file=ff_stat)
  close(ff_stat)
  close(ff_fixed)
  close(ff_min_freq)
  close(ff_min_pop)
  
}

args <- commandArgs(trailingOnly = TRUE)

# cat(args[1],"\n")
# cat(args[2],"\n")
# cat(args[3],"\n")


cat("Removing triallelic, selecting SNPs in at least",args[2],"populations...\n")

trial_min(tab_in=args[1],min_pop=as.numeric(args[2]),min_freq=as.numeric(args[3]),righe=args[4],nome_out=args[5])
cat("done.\n")

################################################################