#! /bin/bash

###############################################################################
##prepare 
##INPUT


scaf_num=(scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 scaffold_6 scaffold_7 scaffold_8)


fold_pops=(07C 07D 07E 07F 07G 07J 07K 07L 07M 07N 07O 07P 07Q 07R 11C 11AA 11AB 11AC 11AE 11AG 11AH 11AJ 11B 11G 11H 11J 11K 11L_F0 11L_F1 11M 11N 11O 11P 11Q 11S 11T 11U 11W 11X 11Z 11D 11V 11R 11A 11E 11F 11Y 14A 14B 14C 14D 14E 07H 13ARE Ha21 Ha31)

#fold_pops=(07C 07D 07E 07F 07G 07J 07K 07L 07M 07N 07O 07P 07Q 07R 11C 11AA 11AB 11AC 11AE 11AG 11AH 11AJ 11B 11G 11H 11J 11K 11L_F0 11L_F1 11M 11N 11O 11P 11Q 11S 11T 11U 11W 11X 11Z 11D 11V 11R 11A 11E 11F 11Y 14A 14B 14C 14D 14E 07H)

n_threads=1
folder="/scicore/home/williy/GROUP/snp_BED_50X/"
type="tot"

out_dir="/scicore/home/williy/fracasse/sift_50x/"

sift_pro="/scicore/home/williy/fracasse/sift_50x/SIFT4G_Annotator_v2.4.jar"
database="/scicore/home/williy/fracasse/v.1.0.23/"

#scaf_num=(scaffold_1 scaffold_2)
#fold_pops=(07C 07D)

mkdir $out_dir"ris_sift"


j=0
len_pop=${#fold_pops[*]}
while [ $j -lt $len_pop ]; do
  
  #############################################################################

 
  i=0
  len=${#scaf_num[*]}
  while [ $i -lt $len ]; do
    echo ${scaf_num[$i]}"--"${fold_pops[$j]}
    qsub -v scaf=${scaf_num[$i]},pop=${fold_pops[$j]},folder=$folder,type=$type,out_dir=$out_dir,$sift_pro=sift_pro,database=$database -pe smp $n_threads -o "p"${fold_pops[$j]}"_"${scaf_num[$i]}".out" -N "p"${fold_pops[$j]}"_"${scaf_num[$i]} SIFT_work.sh
    
  let i++
  done


  
let j++
done






