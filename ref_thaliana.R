###to get reference of thaliana

#first run prepmfa.sh

library(seqinr)
is.even <- function(x) x %% 2 == 0 #pari
is.odd <- function(x) x %% 2 != 0 #dispari

########################################3
###input

setwd("/home/marco/Araly1_Genome/A_lyrata_Vs_TAIR10/")
#chr=c("ChrC", "ChrM")
#name_output="All_Plastid"
chr=c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
name_output="All_regions"


#############################################################################
##select region of A. thaliana that align on the A. lyrata scaffold with A. lyrata + strand and coordinates
##first run /home/marco/Araly1_Genome/A_lyrata_Vs_TAIR10/prepmfa.sh
sca_tot=NULL
start_tot=NULL
end_tot=NULL
gap_th=NULL # gap position in thaliana alignement
gap_ly=NULL
seq_new=NULL
nome_tot=NULL
score_tot=NULL
mono_tot=NULL

for(j in 1:length(chr)) {
  cat(i,"\n")
  seq=read.fasta(paste(chr[j],"_mod.mfa",sep=""))
  l_seq=length(seq)
  seq_ly=seq[c(1:l_seq)[is.even(1:l_seq)]]
  seq_th=seq[c(1:l_seq)[is.odd(1:l_seq)]]
  score=read.table(paste(chr[j],".score",sep=""),stringsAsFactors = FALSE)
  for(i in 1:length(seq_ly)) {
    nome=unlist(strsplit(names(seq_ly[i]),":"))
    nome2=paste(strsplit(names(seq_th[i]),":")[[1]][4:5],collapse=":")
    sca_tot=c(sca_tot,nome[3])
    names(seq_th)[i]=paste(paste(nome[3],nome[4],sep=":"),nome2,sep="&")
    start_tot=c(start_tot,strsplit(nome[4],"-")[[1]][1])
    end_tot=c(end_tot,strsplit(nome[4],"-")[[1]][2])
 
    #gap in lyr
    gap_sel=which(seq_ly[i][[1]]=="-")
    gap_ly=c(gap_ly,length(gap_sel))
    if(length(gap_sel)!=0)
      seq_th[i][[1]]=seq_th[i][[1]][-gap_sel]
    
    if(nome[5]=="(-)") { #if - do revcomp
      seq_th[i][[1]]=rev(comp(seq_th[i][[1]]))
      gap=which(is.na(seq_th[i][[1]])) ##gap in th
      if(length(gap)!=0)
        seq_th[i][[1]][gap]="n"
      gap_th=c(gap_th,length(gap))
    } else {
      gap=which(seq_th[i][[1]]=="-") ##gap in th
      if(length(gap)!=0)
        seq_th[i][[1]][gap]="n"
      gap_th=c(gap_th,length(gap))
    }
  }
  score_tot=c(score_tot,score[,1])
  mono_tot=c(mono_tot,score[,2])
  nome_tot=c(nome_tot,names(seq_th))
  seq_new=c(seq_new,seq_th)

}
tab_fin=cbind(sca_tot,start_tot,end_tot,score_tot,gap_th,gap_ly,mono_tot)
rownames(tab_fin)=nome_tot
write.table(tab_fin, quote=FALSE, file= paste(name_output,"_tab",sep=""))
write.fasta(sequences = seq_new, names = names(seq_new), nbchar = 80, file.out = paste(name_output,".fasta",sep=""))

####################################################################################

###create different fasta file one for each A. lyrata scaffold
##sort the alignement in base of score-(gap_th+gap_ly) ascending order

scafs=c("scaffold_1", "scaffold_2", "scaffold_3", "scaffold_4", "scaffold_5", "scaffold_6", "scaffold_7", "scaffold_8")
tab_fin=read.table("All_regions_tab",stringsAsFactors = FALSE)
align=read.fasta("All_regions.fasta")



for (i in 1:length(scafs)) {
  nomi=rownames(tab_fin)[which(tab_fin$sca_tot==scafs[i])]
  tab=tab_fin[nomi,]
  score=tab$score_tot-(tab$gap_th+tab$gap_ly)
  score=order(score)
  nomi=rownames(tab[score,])
  fasta=align[nomi]
  write.fasta(sequences =fasta, names = names(fasta), nbchar = 80, file.out = paste(scafs[i],"_th.fasta",sep=""))
  cat(i,"\n")
}


##########################################################################
###create hybrid chromosomes inser each alignement on A. lrata chromosomes

library(seqinr)
setwd("/home/marco/Araly1_Genome/A_lyrata_Vs_TAIR10/")



scafs=c("scaffold_1", "scaffold_2", "scaffold_3", "scaffold_4", "scaffold_5", "scaffold_6", "scaffold_7", "scaffold_8")
tab_fin=read.table("All_regions_tab",stringsAsFactors = FALSE)

cat("inizio\n")
for(i in 1:length(scafs)) {
  ins=read.fasta(paste(scafs[i],"_th.fasta",sep=""))
  chr=read.fasta(paste("lyr_th_scaf/",scafs[i],".fa",sep=""))
  for(j in 1:length(ins)) {
    start=tab_fin[names(ins[j]),"start_tot"]
    end=tab_fin[names(ins[j]),"end_tot"]
    chr[[1]][start:end]=ins[[j]]
    cat(i,j,"\n")
  }
  write.fasta(sequences =chr, names = paste(scafs[i],"_th",sep=""), nbchar = 80, file.out = paste("lyr_th_scaf/",scafs[i],"_th.fa",sep=""))
  
}



######################################################################
#select group of alignment
#get all region done
setwd("/home/marco/Araly1_Genome/A_lyrata_Vs_TAIR10/")
tab_fin=read.table("All_regions_tab",stringsAsFactors = FALSE)
tab_fin[1:4,]

tab_sca=tab_fin[which(tab_fin$sca_tot=="scaffold_1"),]


int_start=467605
int_end=467838

int_start=467488
int_end=468492

int_start=456557
int_end=467501
#>A. lyrata scaffold_1:456557-467501 (-)
sel=intersect(which(tab_sca$start_tot <= int_end),which(tab_sca$end_tot >= int_start))


tab_sca[sel,]

###separate dirty.R
#separe fasta file dirty

ff=readLines("/home/marco/program/bwa_genomes/Araly1_assembly_scaffolds.fasta")

nn=grep(">scaffold_9$",ff)

cat(ff[1:(nn-1)],sep="\n",file="/home/marco/program/bwa_genomes/Ara18.fasta")

cat(ff[nn:length(ff)],sep="\n",file="/home/marco/program/bwa_genomes/Ara_small.fasta")

