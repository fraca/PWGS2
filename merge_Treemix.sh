#! /bin/bash
################################################
##ALL the three scripts togheter


fold_pops=(07C 07D 07E 07F 07G 07J 07K 07L 07M 07N 07O 07P 07Q 07R 11C 11AA 11AB 11AC 11AE 11AG 11AH 11AJ 11B 11G 11H 11J 11K 11L_F0 11M 11N 11O 11P 11Q 11S 11T 11U 11W 11X 11Z 11D 11V 11R 11A 11E 11F 11Y 14A 14B 14C 14D 14E 07H Ha31)

folder="/home/marco/pool_first/snp_BED"
type="intergenic_1000bp"
nome_out="intergenic_1000bp_f002_p1"

min_pop=1
min_freq=0.02



################
##merge different positions
#used multiIntersectBed (from bedtools 2.22 copy /home/marco/program/bedtools2/bin to /usr/bin)
# merge_position2.sh

scaf_tot=(scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 scaffold_6 scaffold_7 scaffold_8) 


unicum_m=""
i=0
len=${#fold_pops[*]}

len2=${#scaf_tot[*]}

while [ $i -lt $len ]; do
  #echo "$i: ${fold_pops[$i]}"
  j=0
  while [ $j -lt $len2 ]; do
    #echo "$i $j"
    bedtools sort -i $folder"/p"${fold_pops[$i]}"_"$type"_SNP_scaf/"${scaf_tot[$j]}"_filt.BED" > $nome_out"tmp.BED"
    mv $nome_out"tmp.BED" $folder"/p"${fold_pops[$i]}"_"$type"_SNP_scaf/"${scaf_tot[$j]}"_filt.BED"
    
    unicum_m=$unicum_m$folder"/p"${fold_pops[$i]}"_"$type"_SNP_scaf/"${scaf_tot[$j]}"_filt.BED "
    let j++
  done
  let i++
done
#unicum=${unicum%?} X togliere ultimo char

#echo $unicum_m

multiIntersectBed -i $unicum_m > $nome_out"_multi.BED"
awk '$4=='$len' {print $1"\t"$2"\t"$3}' $nome_out"_multi.BED" > $nome_out".BED"

rm $nome_out"_multi.BED" ##big file to remove.

#####################################################################################
# select_SNP2.sh
################

dir_out="/home/marco/pool_first/merge_snp/"$nome_out

int="/home/marco/pool_first/merge_snp/"$nome_out".BED"

###

rm -r $dir_out
mkdir $dir_out

scaf_num=(scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 scaffold_6 scaffold_7 scaffold_8)



i=0
len=${#fold_pops[*]}
while [ $i -lt $len ]; do
  echo "$i: ${fold_pops[$i]}"
  
  touch $dir_out"/"${fold_pops[$i]}"_varscan"

  j=0
  len2=${#scaf_num[*]}
  while [ $j -lt $len2 ]; do
    /home/marco/program/bedtools2/bin/bedtools intersect -u -a $folder"/p"${fold_pops[$i]}"_"$type"_SNP_scaf/"${scaf_num[$j]}".varscan" -b $int -f 1  >> $dir_out"/"${fold_pops[$i]}"_varscan"
  
  let j++
  done
  
  ##R scripts
  Rscript Varscan2Freq.R $dir_out"/"${fold_pops[$i]}"_varscan"
  
  cut -f 1 $dir_out"/"${fold_pops[$i]}"_varscan" >> $dir_out"/SNP_tot_varscan_t"
  
  
let i++
done

#selction shared SNPs
sort -u $dir_out"/SNP_tot_varscan_t" > $dir_out"/SNP_tot_varscan"
rm $dir_out"/SNP_tot_varscan_t"


#####################################################################################
# merge_SNP2.sh
################


#lascia cosi' di meno da problemi x trial_min.R
ulimit -v 3500000 -m 3500000 

rm $nome_out"_varscan_pops"

i=0
len=${#fold_pops[*]}
while [ $i -lt $len ]; do

  echo ${fold_pops[$i]}"_varscan" >> $nome_out"_varscan_pops"
  let i++
  
done


righe=$dir_out"/SNP_tot_varscan"
Rscript merge_SNPs_pops.R $nome_out"/" $righe $nome_out"_varscan"
Rscript trial_min.R $nome_out"_varscan" $min_pop $min_freq $righe $nome_out"_varscan_clean"

#Selection of the SNP, remove triallelic.

##input TreeMix
awk '{OFS=""; ORS=""; OFMT = "%.0f"; for(i=6;i<=NF;i++) print $i*50,",",50-($i*50)," ";print "\n"}' $nome_out"_varscan_clean.BED" > $nome_out"_varscan.TreeMix"
sed 's/_varscan//g' $nome_out"_varscan_pops" | sed 's/^/p/g' | paste -s -d ' ' > $nome_out"_varscan_pops2"
echo -e  "$(cat $nome_out"_varscan_pops2")\n$(cat $nome_out"_varscan.TreeMix")" | gzip -c > $nome_out"_varscan.TreeMix.gz"
rm $nome_out"_varscan.TreeMix"
rm $nome_out"_varscan_pops"
rm $nome_out"_varscan_pops2"

