#! /bin/bash


chr=(Chr1 Chr2 Chr3 Chr4 Chr5 ChrC ChrM )


for ((j=0; j<${#chr[*]}; j++))
do
  echo $j ${chr[$j]}
  grep "^=" ${chr[$j]}.mfa | sed 's/= score = //g' | sed 's/  type = /\t/g' > ${chr[$j]}.score
  sed -e 's/ /\:/g' ${chr[$j]}.mfa | grep "^=" -v > ${chr[$j]}_mod.mfa
done
