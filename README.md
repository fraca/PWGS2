PWGS2
=====

It analyzes the output of the PWGS pipeline.


- Create reference genome of *Arabidopsis thaliana*.
- Annotation with SnpEff.
- Weighted median of NPStat output.
- Phylogeny with TreeMix.


Software used:

- SnpEff (Cingolani et al. 2012)
- SeqinR (Charif & Lobry et al. 2007)
- bedtools (BEDTools, Quinlan et al. 2010)
- R Core Team (2016)
- TreeMix (Pickrell & Pritchard 2012)

##Create reference genome of *Arabidopsis thaliana*.   
It uses the multiple genome alignment (http://pipeline.lbl.gov/downloads.shtml) between *A. thaliana* (TIAR 10) and *A. lyrata*  (version 1) done with GenomeVISTA (Dubchak et al. 2009). First it sorts the aligned regions in an ascending order according to alignment score minus the number of gap in *A. thaliana* and *A. lyrata*. After it replace the 8 scaffold of *A. lyrata* with the aligned regions of *A. thaliana*. It start from the aligned region with lower score. The regions without alignment are encoded with N.  

It uses the functions in the R script **ref_thaliana.R** and the script **prepmfa.sh**.

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




##Weighted median of NPStat output.  
It calculates the median and the mean of different statistics  of NPStat (Ferretti et al 2013). It calculates the quantiles on the weighted median based on the number of the bp sequenced in the window analyzed. Run the function dist_np2 in the R script **npstat_w_median.R**  

INPUT  
fold_pops=folders with the VarScan SNPs  
 
name_pops= names of the populations  
snp_call= SNPcalling software used (snape or varscan)  
sel_bp= minimum number of base pairs to take into account a window  
sel_S= minimum number of SNPs to take into account a window  
name_out=output name  

OUTPUT  
table with summary foreach pop.  


##Phylogeny with TreeMix.  
The script **merge_Treemix.sh** create the input files for TreeMix (Pickrell & Pritchard 2012). The script **treemix_ini.sh** run TreeMix program. This pipeline is designed to run on sciCORE HPC with qsub command.  

**merge_Treemix.sh**  

It selects the biallelic SNPs present in regions sequenced in all populations and create the input file for TreeMix. It uses the functions in the R script **merge_SNPs_pops.R** and **trial_min.R**  

INPUT  

fold_pops= names of the populations  
folder=folder with the VarScan SNPs  
type= type of region to analyze  
nome_out=output name  
min_pop= minimum number of populations where the SNP is present  
min_freq= minimum frequency for a SNPs  

OUTPUT  

nome_out_varscan.TreeMix.gz input file for TreeMix  
nome_out_varscan_clean.BED table with biallelic SNPs present in regions sequenced in all populations  
nome_out_varscan_clean_col colnames of the table above  
nome_out_varscan table with all the SNPs  
nome_out_varscan_col colnames of the table above  
nome_out_varscan_clean_rfixed table with SNPs fixed in all populations  
nome_out_varscan_clean_rmin_freq table with SNPs that have frequency less than min_freq in all populations  
nome_out_varscan_clean_rmin_pop table with SNPs that are present in less than min_pop  
nome_out_varscan_clean_rtri table with triallelic SNPs  
nome_out_varscan_clean_stat print the number of SNPs removed at each step  




**treemix_ini.sh**  
It run maximum likelihood trees done with TreeMix. It use the **graph_Treemix.R**, **plotting_funcs.R** and **treemix_work.sh** scripts.

INPUT  

bin_dir= directory with the executable files  
in_file= table with biallelic SNPs present in regions sequenced in all populations  
in_pops= colnames of the table above  
ord_pop= list of the names of the populations in the order you would like
them to be plotted  
link= number of SNPs for each block (500)  
n_sim= number of tree to generate  
migs= number of migration events (1 2 3 4 5 6)  
name_out= output name  

OUTPUT  
directory with the TreeMix trees and the graph trees  
name_out_llik It print the inital and final likelihood for each tree  


