#!/bin/bash
#SBATCH --account=grdi_genarcc
#SBATCH --partition=standard
#SBATCH --job-name=readgroups
#SBATCH --output=logs/readgroups
#SBATCH --error=logs/readgroups
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
conda activate $BASEDIR/env/assembly

for i in `cat $BASEDIR/bamsamples.txt` 
do gatk AddOrReplaceReadGroups \
    I=$BASEDIR/Bams/$i.bam \
    O=$BASEDIR/Bams/$i.RG.bam \
    SORT_ORDER=coordinate \
    RGLB=IC \
    RGPU=x \
    RGPL=illumina \
    RGSM=$i
done
