#######################################################
#Calculate frequencies from Varscan
###################################################

##support function for select_SNPs2.sh


###INPUT
#file_in contain each chromosome


###OUTPUT
#rewrite the same file
#
#chromosomes position position ref_allele alt_allele alt_allele_freq
# scaffold_1  141  141  C	T	1
# scaffold_1	351	351	C	A	0.424242424242424
# scaffold_1	845	845	T	G	0.0416666666666667
# scaffold_1	849	849	T	C	1
# scaffold_1	1303	1303	C	T	0.024
# scaffold_1	1398	1398	A	G	0.227513227513228
# scaffold_1	1593	1593	C	T	0.0202020202020202
# scaffold_1	1602	1602	C	T	0.370786516853933
# scaffold_1	1623	1623	C	T	0.0236686390532544
# scaffold_1	1825	1825	A	T	0.402173913043478



Varscan2Freq=function(file_in) {
  ##VARSCAN
  #file_in="/home/marco/pool_first/merge_snp/pro4/07C_varscan" 
  tab=read.table(file_in,stringsAsFactors=FALSE)
  tab[,8]=gsub("%","",tab[,8])
  mode(tab[,8])="numeric"
  tab[,8]=tab[,8]/100
  tab=tab[,c(1,2,3,4,20,8)]
  
  #remove triallelic SNP called by VarScan
  dup=!duplicated(tab[,2])
  tab=tab[dup,]
  tab[1:4,]
  rownames(tab)=apply(tab,1,function(x) gsub(" ","",paste(x[c(1,2,4)],collapse="-"))) 
  tab=tab[,5:6]
  write.table(tab,quote=FALSE,sep="\t",row.names=TRUE,col.names=FALSE,file=file_in)     
}


args <- commandArgs(trailingOnly = TRUE)

file_in=args[1]
Varscan2Freq(file_in)

