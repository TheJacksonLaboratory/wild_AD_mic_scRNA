cd /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/scripts/06_scbase

for F in $(cat ../prefix.txt)
	do echo "
	
#PBS -l nodes=1:ppn=8,walltime=50:00:00
#PBS -m e
#PBS -q batch

cd $PBS_O_WORKDIR

module load Anaconda/4.2.0
source activate scbase

mkdir /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/05_emase/${F}

scbase collate --counts --name *.gene.counts -t /projects/howell-lab/yangs/resources/emase.gene2transcripts.tsv -v /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/05_emase/${F} /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/06_scbase/${F}.loom" >>scbase_${F}.sh

qsub scbase_${F}.sh

done