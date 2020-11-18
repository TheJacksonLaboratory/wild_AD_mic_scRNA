# ${PREFIX_STRAIN} and ${STRAIN} are input from qsub -v
# STRAIN = C57BL6J CAST_EiJ PWK_PhJ WSB_EiJ
# PREFIX_STRAIN = B6J.txt CAST.txt PWK.txt WSB.txt 
#STRAIN


# example of qsub
# separate by the variable list by comma
#qsub 03_bowtie.sh -v STRAIN=C57BL6J,PREFIX_STRAIN=B6J.txt
#qsub 03_bowtie.sh -v STRAIN=CAST_EiJ,PREFIX_STRAIN=CAST.txt
#qsub 03_bowtie.sh -v STRAIN=PWK_PhJ,PREFIX_STRAIN=PWK.txt
#qsub 03_bowtie.sh -v STRAIN=WSB_EiJ,PREFIX_STRAIN=WSB.txt

cd /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/scripts


for F in $(cat ${PREFIX_STRAIN})
 do 
 	a=`find ../cfastq/ -name "${F}*" |wc -l`
 	for ((i=0; i<$a; i++))
 	do
 	echo "
#PBS -l nodes=1:ppn=8,walltime=50:00:00
#PBS -m e
#PBS -q batch
	
cd /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/scripts
module load bowtie/1.0.0
module load samtools/1.8

bowtie -q -p 8 -a --best --strata --sam -v 3 /projects/howell-lab/yangs/resources/genomes/${STRAIN}/bowtie1/transcripts ../cfastq/${F}_${i}.fastq | samtools view -bS - > ../03_bowtie/${F}_${i}.bam" >> 03_bowtie/bowtie_${F}_${i}.sh

qsub 03_bowtie/bowtie_${F}_${i}.sh
	done
done



