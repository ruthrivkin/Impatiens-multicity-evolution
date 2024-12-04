#!/bin/bash
#SBATCH --account=grdi_genarcc
#SBATCH --partition=standard
#SBATCH --job-name=trimg
#SBATCH --output=logs/trimg
#SBATCH --error=logs/trimg
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
conda activate $BASEDIR/env/preprocess/cutadaptenv

#Trim files 
cd $BASEDIR/Demultiplex
perl batch_trim.pl keyfile.txt

