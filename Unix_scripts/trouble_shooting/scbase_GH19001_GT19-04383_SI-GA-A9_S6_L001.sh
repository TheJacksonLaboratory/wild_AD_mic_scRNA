
	
#PBS -l nodes=1:ppn=8,walltime=50:00:00
#PBS -m e
#PBS -q batch

cd /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/scripts/06_scbase

module load Anaconda/4.2.0
source activate scbase

scbase collate --counts --name *.gene.counts -t /projects/howell-lab/yangs/resources/emase.transcripts.info -v /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/05_emase/GH19001_GT19-04383_SI-GA-A9_S6_L001 /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/06_scbase/GH19001_GT19-04383_SI-GA-A9_S6_L001.loom
