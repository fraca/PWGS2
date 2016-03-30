#! /bin/bash

###############################################################################
#work in local

#first modify the file /home/marco/program/snpEff/snpEff.config
#adding

## Arabidopsis lyrata version2
#alyrata_v2.genome : Arabidopsis_lyrata

## Arabidopsis lyrata version2
#athalyr_v2.genome : Arabidopsis_thaliana_on_lyrata


cat /home/marco/pool_first/Araly1_Genome/A_lyrata_Vs_TAIR10/lyr_th_scaf/scaffold_1_th.fa /home/marco/pool_first/Araly1_Genome/A_lyrata_Vs_TAIR10/lyr_th_scaf/scaffold_2_th.fa /home/marco/pool_first/Araly1_Genome/A_lyrata_Vs_TAIR10/lyr_th_scaf/scaffold_3_th.fa /home/marco/pool_first/Araly1_Genome/A_lyrata_Vs_TAIR10/lyr_th_scaf/scaffold_4_th.fa /home/marco/pool_first/Araly1_Genome/A_lyrata_Vs_TAIR10/lyr_th_scaf/scaffold_5_th.fa /home/marco/pool_first/Araly1_Genome/A_lyrata_Vs_TAIR10/lyr_th_scaf/scaffold_6_th.fa /home/marco/pool_first/Araly1_Genome/A_lyrata_Vs_TAIR10/lyr_th_scaf/scaffold_7_th.fa /home/marco/pool_first/Araly1_Genome/A_lyrata_Vs_TAIR10/lyr_th_scaf/scaffold_8_th.fa > /home/marco/program/snpEff/data/genomes/athalyr_v2.fa

#select gff file same for both
awk '$3=="CDS" || $3=="gene" || $3=="exon" ||$3=="intron" || $3=="transcript" || $3=="start_codon" || $3=="stop_codon" || $3=="miRNA_gene" || $3=="tRNA_gene"' /home/marco/pool_first/Araly1_Genome/version2/Additional_file-2_10Feb2015.gff > /home/marco/program/snpEff/data/alyrata_v2/genes.gff

# per lyrata
java -Xmx20g -Xss1g -jar snpEff.jar build -gff3 -v alyrata_v2 2>bla


java -Xmx10g -Xss1g -jar snpEff.jar build -gff3 -v athalyr_v2 2>bla





