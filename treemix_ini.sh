#! /bin/bash

###############################################################################




##INPUT

bin_dir="/scicore/home/williy/fracasse/pipe_bin/"
#bin_dir="/scicore/home/williy/fracasse/treemix-1.12/src/"
ord_pop="end_ord"
migs=(3 4 5 6 7)

in_file="int_0_x_treemix.gz"

link=500
name_out="int_0"

link=1000
name_out="int_0_k1000"



in_file="int_0_x_treemix_norm.gz"

link=500
name_out="int_0_norm"

in_file="int_015_x_treemix.gz"

link=500
name_out="int_015"



in_file="int_1000bp_05_x_treemix.gz"
link=500
name_out="int_1000bp_05"

ord_pop="end_ord_no11L"
migs=(3 5 7)
in_file="int_1000bp_05_x_treemix_no11L.gz"
link=500
name_out="int_1000bp_05_x_treemix_no11L"

ord_pop="end_ord"

migs=(3 5 7)
in_file="tot_05_x_treemix.gz"
link=500
name_out="tot_05"

ord_pop="end_ord"

migs=(3 5 6 7)

in_file="tot_0_x_treemix.gz"
link=500
name_out="tot_0"

n_sim=100

##########################


mkdir $name_out"_nomig"
echo "Likelihoods:" > $name_out"_nomig_llik"
  
qsub -v in_file=$in_file,bin_dir=$bin_dir,link=$link,name_out=$name_out,ord_pop=$ord_pop -t 1:$n_sim -tc 50 -o $name_out"_nomig.out" -N $name_out"_nomig" treemix_work_no_mig.sh

##for each scaffold


j=0
len_migs=${#migs[*]}
while [ $j -lt $len_migs ]; do
  
  rm -r $name_out"_m"${migs[$j]}  
  rm $name_out"_m"${migs[$j]}".out"  
  mkdir $name_out"_m"${migs[$j]}  
  echo "Likelihoods:" > $name_out"_m"${migs[$j]}"_llik"
  
  qsub -v in_file=$in_file,bin_dir=$bin_dir,link=$link,mig=${migs[$j]},name_out=$name_out,ord_pop=$ord_pop -t 1:$n_sim -tc 50 -q long.q -o $name_out"_m"${migs[$j]}".out" -N $name_out"_m"${migs[$j]} treemix_work.sh
  
 
let j++
done



#rm $name_out"_pro"
#rm $name_out"_pro2"
#rm $name_out"_pro3"


