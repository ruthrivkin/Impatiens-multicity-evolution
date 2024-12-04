#!/bin/bash
#SBATCH --account=grdi_genarcc
#SBATCH --partition=standard
#SBATCH --job-name=demuxtor
#SBATCH --output=logs/demuxtor
#SBATCH --error=logs/demuxtor
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
-t $BASEDIR/Demultiplex/Toronto.demult.txt \
-b $BASEDIR/seqs/Raw_TorontoSeqs/TO_ApeK1-Plate2Key.txt \
-f $BASEDIR/seqs/Raw_TorontoSeqs/HF2NMCCX2_1_200813_FD07777933_Other__R_200605_ROBELS1_LIBX10_M006_R1.fastq.gz \
-r $BASEDIR/seqs/Raw_TorontoSeqs/HF2NMCCX2_1_200813_FD07777933_Other__R_200605_ROBELS1_LIBX10_M006_R2.fastq.gz \
-F $BASEDIR/Demultiplex/Toronto_ \
-R $BASEDIR/Demultiplex/Toronto_

