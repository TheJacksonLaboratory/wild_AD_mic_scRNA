
cd /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/scripts/04_alntools

for F in $(cat ../prefix.txt)
	do echo "
	
#PBS -l nodes=1:ppn=8,walltime=50:00:00
#PBS -m e
#PBS -q batch

cd $PBS_O_WORKDIR

module load Anaconda/4.2.0
source activate scTools

alntools bam2both -c 16 --multisample /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/03_bowtie/${F} /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/04_alntools/${F}.bin /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/04_alntools/${F}.h5 -v" >>alntools_${F}.sh

qsub alntools_${F}.sh

done



