PWGS2
=====

Pipeline for analyzing output of VarScan (after PWGS pipeline).


- Create reference genome of *Arabidopsis thaliana*.
- Annotation with SnpEff.
- Get weighted median for linear model.
- Phylogeny with TreeMix.


Software used:

- SnpEff (Cingolani et al. 2012)
- SeqinR (Charif & Lobry et al. 2007)



##Create reference genome of *Arabidopsis thaliana*.   
It uses the multiple genome alignment (http://pipeline.lbl.gov/downloads.shtml) between *A. thaliana* (TIAR 10) and *A. lyrata*  (version 1) done with GenomeVISTA (Dubchak et al. 2009). First it sorts the aligned regions in an ascending order according to alignment score minus the number of gap in *A. thaliana* and *A. lyrata*. After it replace the 8 scaffold of *A. lyrata* with the aligned regions of *A. thaliana*. It start from the aligned region with lower score. The regions without alignment are encoded with N.  

It use the functions in the R script **ref_thaliana.R** and the script **prepmfa.sh**.

##Annotation with SnpEff.  
It annotates the INDEL and SNPs called by VarScan with the new annotation of *A. lyrata* (Rawat et al. 2015). The script **PWGS2_snpEff_genome_build.sh** creates the database for SnpEff. The script **PWGS2_snpEff.sh** annotates the INDEL and SNPs for each population.  

INPUT  
bin_dir= directory with the executable files  
scaf_num=names of the analyzed scaffolds  
folder=folder with the VarScan SNPs and INDEL  
fold_pops=names of the populations  
anno_SnpEff=SnpEff database  
out=output directory  

OUTPUT  
out_SNP directory with the annotated SNPs  
out_INDEL directory with the annotated INDEL  




##Get weighted median for linear model.  
It calculates the median and the mean of different statistics  of NPStat (Ferretti et al 2013) and NI (Rand and Kann 1996) and DoS (Stoletzki and Eyre-Walker 2011). It calculates the quantiles on the weighted median based on the number of the bp sequenced in the window analyzed. Run the function dist_np2 in the R script npstat_w_median.R  

INPUT  
fold_pops=c( "/home/marco/pool_first/ge_int/p07C_5000_CDS_SNP_scaf/", "/home/marco/pool_first/ge_int/p07D_5000_CDS_SNP_scaf/", "/home/marco/pool_first/ge_int/p11A_5000_CDS_SNP_scaf/")
 
name_pops=c( "07C", "07D","11A")
snp_call="snape"
sel_bp=1
sel_S=0
name_out=paste("/home/marco/pool_first/gen_div/prova",snp_call,sel_bp,sel_S,sep="_")
dist_np2(fold_pops,name_pops,snp_call,sel_bp,sel_S,name_out,scaffold=1:3)

OUTPUT  
table with summary foreach pop.  


##Phylogeny with TreeMix.  



