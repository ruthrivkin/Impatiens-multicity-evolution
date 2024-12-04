#!/bin/bash
#SBATCH --account=grdi_genarcc
#SBATCH --partition=standard
#SBATCH --job-name=pixy
#SBATCH --output=logs/pixy
#SBATCH --error=logs/pixy
#SBATCH --time=10000
#SBATCH --nodes=1
#SBATCH --mem=200GB
#SBATCH --cpus-per-task=1
#SBATCH --comment="image=registry.maze.science.gc.ca/ssc-hpcs/generic-job:ubuntu22.04,ssh=true,tmpfs_size=20G"

#Set base directory
BASEDIR=/gpfs/fs7/grdi/genarcc/wp3/PolarBears/ICurb
#Source environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate
conda activate /gpfs/fs7/grdi/genarcc/wp3/PolarBears/env/pixy


pixy --stats pi fst dxy \
--vcf $BASEDIR/Analysis/Pixy/pixy.vcf.gz \
--populations $BASEDIR/Analysis/Pixy/PixyPopFile.txt \
--window_size 10000 \
--output_prefix noldprune \
--n_cores 8 
