cd /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/scripts/05_emase


for F in $(cat ../prefix.txt)
do 
	a=`grep "CRS:" emase_dump_${F}.txt | grep -o -E '[0-9]+'`
	
	mkdir /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/05_emase/${F}

for((i=0;i<=$a;i=i+100));
do
	if [ $(($i + 100)) -gt $a ]
	then
		j=$a
	else
		j=$(($i+100))
    fi	
    
	echo "
#PBS -l nodes=1:ppn=8,walltime=50:00:00
#PBS -q batch

module load Anaconda/4.2.0
source activate scTools

cd /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/scripts/05_emase

emase-zero -m 2 -o /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/05_emase/${F}/${i}_${j} -s ${i}:${j} -g /projects/howell-lab/yangs/resources/emase.gene2transcripts.tsv -t 0.000001 -i 9999 -v /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/04_alntools/${F}.bin " >> emase_${F}_${i}_${j}.sh

qsub emase_${F}_${i}_${j}.sh

done

#emase-zero -m 2 -o emase/emase_zero -s 0:6115 -g /projects/howell-lab/yangs/resources/emase.gene2transcripts.tsv -t 0.000001 -i 9999 -v aln/alntools_outbase.bin

done
