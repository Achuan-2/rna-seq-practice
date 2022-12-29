#!/bin/bash
#SBATCH -J trim_galore
#SBATCH -p dna
#SBATCH -N 1
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH -o slurm.%j.%x.out        # STDOUT
#SBATCH -e slurm.%j.%x.err        # STDERR
#SBATCH --mail-type=END  # 发送哪一种email通知：BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jxsu22@m.fudan.edu.cn

mode=$1
SRR=$2
if [ "$mode" == "single" ];then
    # 如果是单端
    trim_galore --illumina --fastqc  $SRR.fastq.gz
elif [ "$mode" == "paired" ];then
    # 如果是多端
    trim_galore --illumina --fastqc --paired ${SRR}_1.fastq.gz ${SRR}_2.fastq.gz
fi
