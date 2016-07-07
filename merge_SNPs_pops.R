
# #folder with file for each pop of SNPs
# folder="/home/marco/pool_first/SNP_database/prova3"
# folder="/home/marco/pool_first/SNP_database/p52"
# # #names of pops
# # fold_pops=c("07C_snape", "11L_F0_snape", "11M_snape")
# fold_pops=c("07C_snape", "07D_snape", "07E_snape", "07F_snape", "07G_snape", "07J_snape", "07K_snape", "07L_snape", "07M_snape", "07N_snape", "07O_snape", "07P_snape", "07Q_snape", "07R_snape", "11AA_snape", "11AB_snape", "11AC_snape", "11AE_snape", "11AG_snape", "11AH_snape", "11AJ_snape", "11B_snape", "11G_snape", "11H_snape", "11J_snape", "11K_snape", "11L_F0_snape", "11M_snape", "11N_snape", "11O_snape", "11P_snape", "11Q_snape", "11S_snape", "11T_snape", "11U_snape", "11W_snape", "11X_snape", "11Z_snape", "14A_snape", "14B_snape", "14C_snape", "14D_snape", "14E_snape", "11D_snape", "11V_snape", "11R_snape", "11A_snape", "11C_snape", "11E_snape", "11F_snape", "11Y_snape", "07H_snape")
# # #rownames of SNPs
# # righe="/home/marco/pool_first/SNP_database/prova3/SNP_tot_snape"
# righe="/home/marco/pool_first/SNP_database/p52/SNP_tot_snape"
# # #output name
# # nome_out="/home/marco/pool_first/SNP_database/prova3_tab_snape"
# nome_out="/home/marco/pool_first/SNP_database/p52_tab_snape"

#OUTPUT one table with all SNPs for all populations
###################################################################

merge_SNPs_pops=function(folder,fold_pops,righe,nome_out) {
  
  system(paste("rm ",nome_out,sep=""))

  righe=readLines(righe)
  #righe=righe[1:10000]
  colonne=NULL
  for(i in 1:length(fold_pops)) {
    colonne=c(colonne,paste("ALT_",fold_pops[i],sep=""),paste("freq_",fold_pops[i],sep=""))
  }
  cat(colonne,sep="\n",file=paste(nome_out,"_col",sep=""))

  #cat("fin qui2\n")
  mem_limit=100000
  #mem_limit=50000
  j_tot=ceiling(length(righe)/mem_limit)
  cat(j_tot,"\n")  
  c1=1
  c2=mem_limit
  for(j in 1:j_tot) {
    tab=matrix(nrow=length(righe[c1:c2]),ncol=length(fold_pops)*2)
    rownames(tab)=righe[c1:c2]
    colnames(tab)=colonne
    
    for(i in 1:length(fold_pops)) {
      cat(j,"of",j_tot,i,"of",length(fold_pops),"\n")
      bla=read.table(paste(folder,"/",fold_pops[i],sep=""),stringsAsFactors = FALSE)
      rownames(bla)=bla[,1]
      bla=bla[,2:3]    
      colnames(bla)=c(paste("ALT_",fold_pops[i],sep=""),paste("freq_",fold_pops[i],sep=""))
      
      rig_sel=intersect(rownames(tab),rownames(bla))
      bla=bla[rig_sel,]
      tab[rig_sel,paste("ALT_",fold_pops[i],sep="")]=bla[,paste("ALT_",fold_pops[i],sep="")]
      #cat("fin qui5 ",i,j,"\n")
      tab[rig_sel,paste("freq_",fold_pops[i],sep="")]=bla[,paste("freq_",fold_pops[i],sep="")]
      #cat("fin qui6 ",i,"\n")

    }    
    c1=c1+mem_limit
    c2=c2+mem_limit
    if(c2>length(righe))
      c2=length(righe) 
    write.table(tab,quote=FALSE,append=TRUE,row.names=FALSE,col.names=FALSE,file=nome_out)
    cat(dim(tab),"\n")
  }
  
}
#fine  


args <- commandArgs(trailingOnly = TRUE)

# cat(args[1],"\n")
# cat(args[2],"\n")
# cat(args[3],"\n")

folder=args[1]
righe=args[2]
nome_out=args[3]
fold_pops=readLines(paste(nome_out,"_pops",sep=""))

cat("Merging of all SNPs ...\n")
merge_SNPs_pops(folder,fold_pops,righe,nome_out)
cat("done.\n")

#############################################################