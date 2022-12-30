#!/bin/bash
#SBATCH -J cuffdiff
#SBATCH -p dna
#SBATCH -N 1
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH -o slurm.%j.%x.out   # STDOUT
#SBATCH -e slurm.%j.%x.err   # STDERR
#SBATCH --mail-type=END  # 发送哪一种email通知：BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jxsu22@m.fudan.edu.cn

outdir=$1
gtf=$2
bam_dir=$3
label1=$4
label2=$5

sample1=$(ls $bam_dir/*$label1* | xargs | tr ' ' ',')
sample2=$(ls $bam_dir/*$label2* | xargs | tr ' ' ',')
cuffdiff -p 2 -o $outdir \
    -L $label1,$label2 \
    $gtf \
    $sample1 \
    $sample2
