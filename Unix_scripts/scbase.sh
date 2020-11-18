#PBS -l nodes=1:ppn=8,walltime=50:00:00
#PBS -m e
#PBS -q batch
#PBS -M Stanley.Yang@jax.org

cd /projects/howell-lab/yangs/projects/2018_10_scRNA_test1/test1

module load Anaconda/4.2.0
source activate scbase

scbase collate --counts --name emase_*.gene.counts -t /projects/howell-lab/yangs/resources/GeneIDs.tsv -v emase_output scbase/myeloid.loom
