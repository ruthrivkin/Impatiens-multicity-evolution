#!/bin/bash
#SBATCH --account=grdi_genarcc
#SBATCH --partition=standard
#SBATCH --job-name=callsnps
#SBATCH --output=logs/callsnps.allsites
#SBATCH --error=logs/callsnps.allsites
#SBATCH --time=48:00:00
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
BAM=$BASEDIR/SeedSamples.txt

#Variant only
#bcftools mpileup -Ou -b $BAM -f $REF -I -a AD,DP,SP,ADF,ADR -d 50 | bcftools call -a GQ,GP -mv -Ov > SNPs/variants.vcf

#All sites
bcftools mpileup -Ou -b $BAM -f $REF -I -a AD,DP,SP,ADF,ADR -d 50 | bcftools call -a GQ,GP -m -Ov > SNPs/allsites.vcf
