#! /bin/bash

###############################################################################
#work in local
##INPUT



scaf_num=(scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 scaffold_6 scaffold_7 scaffold_8)

folder="/home/marco/pool_first/snp_BED_tot/"

fold_pops=(07C 07D 07E 07F 07G 07J 07K 07L 07M 07N 07O 07P 07Q 07R 11C 11AA 11AB 11AC 11AE 11AG 11AH 11AJ 11B 11G 11H 11J 11K 11L_F0 11L_F1 11M 11N 11O 11P 11Q 11S 11T 11U 11W 11X 11Z 11D 11V 11R 11A 11E 11F 11Y 14A 14B 14C 14D 14E 07H 13ARE Ha21 Ha31)



#fold_pops=(07C 07D Ha31)

#anno_SnpEff="alyrata_v2"
#out="SnpEff_all"

anno_SnpEff="athalyr_v2"
out="SnpEff_all_th"

#############################################################################

mkdir $out"_SNP"
mkdir $out"_INDEL"
temp="temp_"$out

j=0
len_pop=${#fold_pops[*]}
while [ $j -lt $len_pop ]; do
  
  echo ${fold_pops[$j]}
  
  echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"  > $temp"_SNP.vcf"
  
  echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"  > $temp"_INDEL.vcf"
  
  i=0
  len=${#scaf_num[*]}
  while [ $i -lt $len ]; do
    
    awk '{freq=$8; gsub("%","",freq); freq=freq/100; rede=$6+$7; if($5=="Y") {if($4=="C") {print $1"\t"$2"\t.\t"$4"\tT\t99\tPASS\tAF="freq";DP="rede} else {print $1"\t"$2"\t.\t"$4"\tC\t99\tPASS\tAF="freq";DP="rede}} else if($5=="R") {if($4=="G") {print $1"\t"$2"\t.\t"$4"\tA\t99\tPASS\tAF="freq";DP="rede} else {print $1"\t"$2"\t.\t"$4"\tG\t99\tPASS\tAF="freq";DP="rede}} else if($5=="W") {if($4=="A") {print $1"\t"$2"\t.\t"$4"\tT\t99\tPASS\tAF="freq";DP="rede} else {print $1"\t"$2"\t.\t"$4"\tA\t99\tPASS\tAF="freq";DP="rede}} else if($5=="S") {if($4=="C") {print $1"\t"$2"\t.\t"$4"\tG\t99\tPASS\tAF="freq";DP="rede} else {print $1"\t"$2"\t.\t"$4"\tC\t99\tPASS\tAF="freq";DP="rede}} else if($5=="K") {if($4=="G") {print $1"\t"$2"\t.\t"$4"\tT\t99\tPASS\tAF="freq";DP="rede} else {print $1"\t"$2"\t.\t"$4"\tG\t99\tPASS\tAF="freq";DP="rede}} else if($5=="M") {if($4=="C") {print $1"\t"$2"\t.\t"$4"\tA\t99\tPASS\tAF="freq";DP="rede} else {print $1"\t"$2"\t.\t"$4"\tC\t99\tPASS\tAF="freq";DP="rede}} else {print $1"\t"$2"\t.\t"$4"\t"$5"\t99\tPASS\tAF="freq";DP="rede}}' $folder"p"${fold_pops[$j]}"_tot_SNP_scaf/"${scaf_num[$i]}".varscan" >> $temp"_SNP.vcf"
    
    awk '{freq=$8; gsub("%","",freq); freq=freq/100; rede=$6+$7; if($5=="Y") {if($4=="C") {print $1"\t"$2"\t.\t"$4"\tT\t99\tPASS\tAF="freq";DP="rede} else {print $1"\t"$2"\t.\t"$4"\tC\t99\tPASS\tAF="freq";DP="rede}} else if($5=="R") {if($4=="G") {print $1"\t"$2"\t.\t"$4"\tA\t99\tPASS\tAF="freq";DP="rede} else {print $1"\t"$2"\t.\t"$4"\tG\t99\tPASS\tAF="freq";DP="rede}} else if($5=="W") {if($4=="A") {print $1"\t"$2"\t.\t"$4"\tT\t99\tPASS\tAF="freq";DP="rede} else {print $1"\t"$2"\t.\t"$4"\tA\t99\tPASS\tAF="freq";DP="rede}} else if($5=="S") {if($4=="C") {print $1"\t"$2"\t.\t"$4"\tG\t99\tPASS\tAF="freq";DP="rede} else {print $1"\t"$2"\t.\t"$4"\tC\t99\tPASS\tAF="freq";DP="rede}} else if($5=="K") {if($4=="G") {print $1"\t"$2"\t.\t"$4"\tT\t99\tPASS\tAF="freq";DP="rede} else {print $1"\t"$2"\t.\t"$4"\tG\t99\tPASS\tAF="freq";DP="rede}} else if($5=="M") {if($4=="C") {print $1"\t"$2"\t.\t"$4"\tA\t99\tPASS\tAF="freq";DP="rede} else {print $1"\t"$2"\t.\t"$4"\tC\t99\tPASS\tAF="freq";DP="rede}} else {print $1"\t"$2"\t.\t"$4"\t"$5"\t99\tPASS\tAF="freq";DP="rede}}' $folder"p"${fold_pops[$j]}"_tot_SNP_scaf/"${scaf_num[$i]}"_indel.varscan" >> $temp"_INDEL.vcf"

    
  let i++
  done
  
  java -Xmx8g -jar /home/marco/program/snpEff/snpEff.jar -c /home/marco/program/snpEff/snpEff.config -stats $out"_SNP/p"${fold_pops[$j]}"_genes" $anno_SnpEff $temp"_SNP.vcf" > $out"_SNP/p"${fold_pops[$j]}

  java -Xmx8g -jar /home/marco/program/snpEff/snpEff.jar -c /home/marco/program/snpEff/snpEff.config -stats $out"_INDEL/p"${fold_pops[$j]}"_genes" $anno_SnpEff $temp"_INDEL.vcf" > $out"_INDEL/p"${fold_pops[$j]}

  rm $out"_INDEL/p"${fold_pops[$j]}"_genes"
  rm $out"_SNP/p"${fold_pops[$j]}"_genes"
  rm $temp"_SNP.vcf"
  rm $temp"_INDEL.vcf"
  
let j++
done



