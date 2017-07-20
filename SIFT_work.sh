#! /bin/bash

# Request "/bin/bash" as shell
#$ -S /bin/bash

# Start the job from the current working directory
#$ -cwd

# Merge standard output and standard error
#$ -j y

#  how much memory? DEFAULT = 2G 
########################
#$ -l membycore=16G

###############################################################################


# to load gsl library fro NPStat

awk '{freq=$8; gsub("%","",freq); uno=$1; gsub("scaffold_","",uno); freq=freq/100; rede=$6+$7; if($5=="Y") {if($4=="C") {print uno"\t"$2"\t.\t"$4"\tT\t99\tPASS\tAF="freq";DP="rede} else {print uno"\t"$2"\t.\t"$4"\tC\t99\tPASS\tAF="freq";DP="rede}} else if($5=="R") {if($4=="G") {print uno"\t"$2"\t.\t"$4"\tA\t99\tPASS\tAF="freq";DP="rede} else {print uno"\t"$2"\t.\t"$4"\tG\t99\tPASS\tAF="freq";DP="rede}} else if($5=="W") {if($4=="A") {print uno"\t"$2"\t.\t"$4"\tT\t99\tPASS\tAF="freq";DP="rede} else {print uno"\t"$2"\t.\t"$4"\tA\t99\tPASS\tAF="freq";DP="rede}} else if($5=="S") {if($4=="C") {print uno"\t"$2"\t.\t"$4"\tG\t99\tPASS\tAF="freq";DP="rede} else {print uno"\t"$2"\t.\t"$4"\tC\t99\tPASS\tAF="freq";DP="rede}} else if($5=="K") {if($4=="G") {print uno"\t"$2"\t.\t"$4"\tT\t99\tPASS\tAF="freq";DP="rede} else {print uno"\t"$2"\t.\t"$4"\tG\t99\tPASS\tAF="freq";DP="rede}} else if($5=="M") {if($4=="C") {print uno"\t"$2"\t.\t"$4"\tA\t99\tPASS\tAF="freq";DP="rede} else {print uno"\t"$2"\t.\t"$4"\tC\t99\tPASS\tAF="freq";DP="rede}} else {print uno"\t"$2"\t.\t"$4"\t"$5"\t99\tPASS\tAF="freq";DP="rede}}' $folder"p"$pop"_"$type"_SNP_scaf/"$scaf".varscan" > $out_dir"p"$pop"_"$scaf"_tm.vcf"

echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" | cat - $out_dir"p"$pop"_"$scaf"_tm.vcf" > $out_dir"p"$pop"_"$scaf".vcf"

rm $out_dir"p"$pop"_"$scaf"_tm.vcf"
 
echo "vcf create..."


java -jar SIFT4G_Annotator_v2.4.jar -c -i $out_dir"p"$pop"_"$scaf".vcf" -d $database -r $out_dir"/ris_sift/"


echo "Done"