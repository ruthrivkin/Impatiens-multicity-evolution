#!/bin/bash
#SBATCH --account=grdi_genarcc
#SBATCH --partition=standard
#SBATCH --job-name=demux34
#SBATCH --output=logs/demux34
#SBATCH --error=logs/demux34
#SBATCH --time=10000
#SBATCH --nodes=1
#SBATCH --mem=200GB
#SBATCH --cpus-per-task=1
#SBATCH --comment="image=registry.maze.science.gc.ca/ssc-hpcs/generic-job:ubuntu22.04,ssh=true,tmpfs_size=20G"
#sbatch scripts/bcftools_chr_intervals_abalone.sh 

#Set base directory
BASEDIR=/gpfs/fs7/grdi/genarcc/wp3/PolarBears/ICurb
#Source environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate
conda activate $BASEDIR/env/preprocess

#Demultiplex Plate 1.2 samples

axe-demux -m0 -z 3 -c -v \
-t $BASEDIR/Demultiplex/3.4.demult.txt \
-b $BASEDIR/seqs/Raw_Sequence_3_4/3.4_ApeKI-188PlexKey.txt \
-f $BASEDIR/seqs/Raw_Sequence_3_4/40580-3.4_R1.fastq.gz \
-r $BASEDIR/seqs/Raw_Sequence_3_4/40580-3.4_R2.fastq.gz \
-F $BASEDIR/Demultiplex/Plate3.4_ \
-R $BASEDIR/Demultiplex/Plate3.4_
