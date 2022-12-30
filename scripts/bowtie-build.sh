#!/bin/bash
#SBATCH -J bowtie-build
#SBATCH -p dna
#SBATCH -N 4
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH -o slurm.%j.%x.out        # STDOUT
#SBATCH -e slurm.%j.%x.err        # STDERR
#SBATCH --mail-type=END  # 发送哪一种email通知：BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jxsu22@m.fudan.edu.cn

bowtie2-build m8.fa m8