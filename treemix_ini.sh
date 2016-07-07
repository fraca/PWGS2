#! /bin/bash

###############################################################################




##INPUT

bin_dir="/scicore/home/williy/fracasse/pipe_bin/"
#bin_dir="/scicore/home/williy/fracasse/treemix-1.12/src/"

in_file="intergenic_1000bp_f002_p1_varscan_clean.BED"
in_pops="intergenic_1000bp_f002_p1_varscan_clean_col"
ord_pop="end_ord"
link=500
n_sim=100
migs=(1 2 3 4 5 6)

name_out="inter_500"


##########################


awk '{OFS=""; ORS=""; OFMT = "%.0f"; for(i=6;i<=NF;i++) print $i*50,",",50-($i*50)," ";print "\n"}' $in_file > $name_out"_pro"
grep varscan $in_pops | sed 's/freq_//g' | sed 's/_varscan//g' | sed 's/^/p/g' | paste -s -d ' ' > $name_out"_pro2"

cat $name_out"_pro2" $name_out"_pro" > $name_out"_pro3"

gzip -c $name_out"_pro3" > $name_out"_pro3.gz"



echo "Fatto..."

rm -r $name_out"_nomig"
rm $name_out"_nomig.out"
mkdir $name_out"_nomig"
echo "Likelihoods:" > $name_out"_nomig_llik"
  
qsub -v bin_dir=$bin_dir,link=$link,name_out=$name_out,ord_pop=$ord_pop -t 1:100 -tc 50 -o $name_out"_nomig.out" -N $name_out"_nomig" treemix_work_no_mig.sh

##for each scaffold



j=0
len_migs=${#migs[*]}
while [ $j -lt $len_migs ]; do
  
  rm -r $name_out"_m"${migs[$j]}  
  rm $name_out"_m"${migs[$j]}".out"  
  mkdir $name_out"_m"${migs[$j]}  
  echo "Likelihoods:" > $name_out"_m"${migs[$j]}"_llik"
  
  qsub -v bin_dir=$bin_dir,link=$link,mig=${migs[$j]},name_out=$name_out,ord_pop=$ord_pop -t 1:$n_sim -tc 50 -q long.q -o $name_out"_m"${migs[$j]}".out" -N $name_out"_m"${migs[$j]} treemix_work.sh

 
let j++
done



rm $name_out"_pro"
rm $name_out"_pro2"
rm $name_out"_pro3"


