#!/bin/bash
#SBATCH -J fastq-dump
#SBATCH -p dna
#SBATCH -N 1
#SBATCH --mem=10G
#SBATCH --cpus-per-task=2
#SBATCH -o slurm.%N.%j.%x_%A_%a.out        # STDOUT
#SBATCH -e slurm.%N.%j.%x_%A_%a.err        # STDERR
#SBATCH --mail-type=END
#SBATCH --mail-user=jxsu22@m.fudan.edu.cn

fastq-dump --split-3 --gzip $1 --outdir $2
