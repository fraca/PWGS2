#! /bin/bash

# Request "/bin/bash" as shell
#$ -S /bin/bash

#$ -pe smp 1

# Start the job from the current working directory
#$ -cwd

# Merge standard output and standard error
#$ -j y



###############################################################################


# to load gsl library fro NPStat

module load GSL/1.16-intel-2015.02
module load Boost/1.57.0-goolf-1.4.10
module load R/3.2.4-goolf-1.7.20
module load GCC/4.9.2

 
$bin_dir"treemix" -i $name_out"_pro3.gz" -root pHa31 -k $link -m $mig -se -global -o $name_out"_m"$mig"/t_"${SGE_TASK_ID}
  
Rscript graph_Treemix.R $name_out"_m"$mig"/t_"${SGE_TASK_ID} $ord_pop
  
echo -e "t_"${SGE_TASK_ID} >> $name_out"_m"$mig"_llik"
cat $name_out"_m"$mig"/t_"${SGE_TASK_ID}".llik" >> $name_out"_m"$mig"_llik"



echo "Fine."

