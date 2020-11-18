
cd /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/scripts/05_emase

for F in $(cat ../prefix.txt)
do 
  a=`grep "CRS:" emase_dump_${F}.txt | grep -o -E '[0-9]+'`

  for((i=0;i<=$a;i=i+100));
  do
    b=$(($i + 100))
    if [ "$b" -ge "$a" ]
    then
        j=$a
    else
        j=$(($i+100))
    fi	
    
  echo ${a}_${i}_${j}

  done
done
