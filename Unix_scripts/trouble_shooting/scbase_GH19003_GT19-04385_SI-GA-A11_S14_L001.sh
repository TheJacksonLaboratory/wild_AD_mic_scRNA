
	
#PBS -l nodes=1:ppn=8,walltime=50:00:00
#PBS -m e
#PBS -q batch

cd /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/scripts/06_scbase

module load Anaconda/4.2.0
source activate scbase

mkdir /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/05_emase/GH19003_GT19-04385_SI-GA-A11_S14_L001

scbase collate --counts --name *.gene.counts -t /projects/howell-lab/yangs/resources/emase.gene2transcripts.tsv -v /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/05_emase/GH19003_GT19-04385_SI-GA-A11_S14_L001 /projects/howell-lab/yangs/projects/2019_04_scRNA_CD11b/GH19001_GH19008/06_scbase/GH19003_GT19-04385_SI-GA-A11_S14_L001.loom
