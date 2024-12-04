#!/bin/bash
#SBATCH --account=grdi_genarcc
#SBATCH --open-mode=append
#SBATCH --partition=standard
#SBATCH --job-name=assemble
#SBATCH --output=logs/assemble
#SBATCH --error=logs/assemble
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
conda activate $BASEDIR/env/assembly

# Set reference genome 
REF=$BASEDIR/Reference/Impatiens_capensis_hirise.fasta

# declare variables
IND=$1
FORWARD=$BASEDIR/Trimmed/${IND}_R1.fq.gz
REVERSE=$BASEDIR/Trimmed/${IND}_R2.fq.gz
OUTPUT=$BASEDIR/Bams/${IND}.bam

# then align and sort
echo "Aligning $IND with bwa"
bwa mem -M -t 4 $REF $FORWARD $REVERSE | \
samtools view -b | \
samtools sort -T ${IND} > $OUTPUT
#end script
