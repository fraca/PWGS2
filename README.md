PWGS2
=====

Severals scripts that perform analysis on pool-seq data obtained with the [PWGS pipeline](https://github.com/fraca/PWGS). They have been used in:

- [ M. Fracassetti Genetic diversity and climate adaption in *Arabidopsis lyrata*](https://edoc.unibas.ch/44796/)
- [Willi Y. et al. Accumulation of Mutational Load at the Edges of a Species Range](https://academic.oup.com/mbe/article/35/4/781/4810769)

Software used:

- SnpEff ([Cingolani et al. 2012](http://snpeff.sourceforge.net/))
- SeqinR ([Charif & Lobry et al. 2007](https://cran.r-project.org/web/packages/seqinr/index.html))
- BEDTools ([BEDTools, Quinlan et al. 2010](http://bedtools.readthedocs.io/en/latest/))
- TreeMix ([Pickrell & Pritchard 2012](https://bitbucket.org/nygcresearch/treemix/wiki/Home))
- Baypass ([Gautier 2015](http://www1.montpellier.inra.fr/CBGP/software/baypass/))
- SIFT4G ([Vaser et al. 2016](http://sift.bii.a-star.edu.sg/sift4g/))

## Create reference genome of *Arabidopsis thaliana*.   
It uses the multiple genome alignment (http://pipeline.lbl.gov/downloads.shtml) between *A. thaliana* (TIAR 10) and *A. lyrata*  (version 1) done with GenomeVISTA ([Dubchak et al. 2009](http://genome.lbl.gov/vista/index.shtml)). First it sorts the aligned regions in an ascending order according to alignment score minus the number of gap in *A. thaliana* and *A. lyrata*. After it replace the 8 scaffold of *A. lyrata* with the aligned regions of *A. thaliana*. It start from the aligned region with lower score. The regions without alignment are encoded with N.  

It uses the functions in the R script **ref_thaliana.R** and the script **prepmfa.sh**.

## Annotation with SnpEff.  
It annotates the INDEL and SNPs called by VarScan with the new annotation of *A. lyrata* ([Rawat et al. 2015](https://doi.org/10.1371/journal.pone.0137391.s002)). The script **PWGS2_snpEff_genome_build.sh** creates the database for SnpEff. The script **PWGS2_snpEff.sh** annotates the INDEL and SNPs for each population.  

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




## Weighted median of NPStat output.  
It calculates the median and the mean of different statistics  of NPStat ([Ferretti et al 2013](https://github.com/lucaferretti/npstat)). It calculates the quantiles on the weighted median based on the number of the bp sequenced in the window analyzed. Run the function dist_np2 in the R script **npstat_w_median.R**  

INPUT  
fold_pops=folders with the VarScan SNPs  
 
name_pops= names of the populations  
snp_call= SNPcalling software used (snape or varscan)  
sel_bp= minimum number of base pairs to take into account a window  
sel_S= minimum number of SNPs to take into account a window  
name_out=output name  

OUTPUT  
table with summary foreach pop.  


## Explanatory variables for linear model  
The script **linear_model_input_table_fin.R** calculates the explanatory for the linear model.


## Merging of the SNP.  
It merges the SNPs called by Varscan. It retains only the SNP biallelic (between populations), present in min_pop populations, that has min_freq allele frequency and that has min_MAF MAF (minor allele frequency). It prepares the input file for Treemix and Baypass. It uses BEDTools,**snp_merge_scicore1.R** and **snp_merge3_MAF.R**  

**snp_merge_scicore1.R**  
INPUT  
folder_in= folder with the output of PWGS  
mod= genomic region (tot, intergenic, cds ...)  
name_pops= selected populations  
scafs= scaffold  
min_freq= minimum frequency to retain a SNP  
min_pop= minimum number of populations where a SNP have to be present.
nome_out= Output name.


OUTPUT  
nome_out_scaf_stat2 info file  
nome_out_scaf table with all the SNPs  
nome_out_scaf_tripop table with triallelic SNPs  
nome_out_scaf_meancov Mean read depth  

**snp_merge3_MAF.R**  
INPUT  
name_tab= Input table from **snp_merge_scicore1.R**  
min_MAF= minimum total MAF (across all populations) to retain a SNPs  
min_pop= minimum number of populations where a SNP have to be present  
pop_sel= selected populations  
name_out= Output name  

OUTPUT  
name_out_min_MAF_min_pop_stat3=info file  
name_out_min_MAF_min_pop.vcf= vcf file with popsitions  
name_out_min_MAF_min_pop.geno= table with number of reference and alternative alleles  
name_out_min_MAF_min_pop_pops= names of populations analyzed  


## Relatedness tree with TreeMix.  
It creates relatedness trees with migration events. The script **treemix_ini.sh** run TreeMix program and call **treemix_work.sh** and **treemix_work_no_mig.sh**. This pipeline is designed to run on sciCORE HPC with qsub command.  

**treemix_ini.sh**  
It run maximum likelihood trees done with TreeMix. It uses the **graph_Treemix.R**, **plotting_funcs.R**, **treemix_work.sh**, **treemix_work_no_mig.sh**.


INPUT  

ord_pop= order of populations  
migs= number of migration events  
in_file= input file for treemix  
link= block of SNPs (for linkage disequilibrium)  
name_out= output name  

n_sim= number of simulations  

OUTPUT  

directory with the TreeMix trees and the graph trees  
name_out_llik It print the inital and final likelihood for each tree  

## Baypass analysis  
Baypass ([Gautier 2015](http://www1.montpellier.inra.fr/CBGP/software/baypass/)) is a program for the Environmental Association Analysis EAA. The input files are done with **snp_merge3_MAF.R**. First the population variant matrix is calculated and after the association with between the SNPs and the environmental variable is performed. Finally a Gene Ontology analysis is done with SNP2GO (Szkiba et al. 2014). In the script **bypass_fin_git.R** there are some function in R that I used. 


## SIFT4G analysis  
SIFT4G predict whether an amino acid substitustion is deleterious. The script **SIFT_ini.sh** run SIFT program and call **SIFT_work.sh**, after run **SIFT_fin.sh**. This pipeline is designed to run on sciCORE HPC with qsub command.   



INPUT  

scaf_num= scaffold names   
fold_pops= popolation names  

n_threads= number cores  
folder= input folder  
type= type of region (tot, int, CDS)  

out_dir= output directory

sift_pro= location SIFT4G_Annotator_v2.4.jar
database= database SIFT4G (v.1.0.23)

OUTPUT  

The SIFT4G_Annotator_v2 output files in the output directory  



