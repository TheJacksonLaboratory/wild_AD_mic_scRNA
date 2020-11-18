#copied this code directly into helix console to batch qsub cfastq files

cd /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/scripts
for F in $(cat prefix.txt)
 do 
echo "

#PBS -l nodes=1:ppn=1,walltime=20:00:00
#PBS -m e
#PBS -q batch
#PBS -M Stanley.Yang@jax.org

cd /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/scripts

module load Anaconda
source activate scTools

dir_F=/projects/howell-lab/00_fastq/yangs/2019_04_scRNA_CD11b

cfastq3 -d -c 5000000 -o ../cfastq/$F -e 1 $dir_F/${F}_I1_001.fastq.gz $dir_F/${F}_R1_001.fastq.gz $dir_F/${F}_R2_001.fastq.gz" >> cfastq_${F}.sh

qsub cfastq_${F}.sh

done
