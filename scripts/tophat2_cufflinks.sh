#!/bin/bash
#SBATCH -J tophat2_cufflinks
#SBATCH -p dna
#SBATCH -N 4
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH -o slurm.%j.%x.out        # STDOUT
#SBATCH -e slurm.%j.%x.err        # STDERR
#SBATCH --mail-type=END  # 发送哪一种email通知：BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jxsu22@m.fudan.edu.cn

echo "Usage:"
echo "   tophat2_cufflinks.sh {mode} {threads} {transcriptome-index} {bowtie2-index} {SRR} {fq1} [{fq2}] "
echo ""


INDEX=$PWD/00_index
DATA=$PWD/01_rawdata
RESULT=$PWD/02_result

mode=$1
if test -z $mode # 检测字符是否为空
then
   echo "please input the mode(single or paired)"
   exit
fi
threads=$2
if test -z $threads
then
   echo "please input the number of threads"
   exit
fi

Annotation=$3
if test -z $Annotation
then
   echo "please input transcriptome-index(/share/Genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/hg19_genes/genes.gff)" 
   exit
fi

bowtie2Index=$4
if test -z $bowtie2Index
then
   echo "please input bowtie2-index(/share/Genomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome)" 
   exit
fi

SRR=$5
if test -z $SRR
then
   echo "please input SRR id"
   exit
fi

fq1=$6
if test -z $fq1
then
   echo "please input fasta1"
   exit
fi


if [ "$mode" == "paired" ];then

   fq2=$7
   if test -z $fq2
   then
      echo "please input fasta1"
      exit
   fi
fi

#====================================================
echo "Running info"
echo "Project:    "$PWD
echo "Read:       "$SRR
echo "Annotation: "$Annotation
echo "Genome:     "$bowtie2Index
echo " "

########################RUN##############################
mkdir -p ${RESULT}/tophat2/${SRR}
mkdir -p ${RESULT}/cufflinks/${SRR}

if [ "$mode" == "single" ];then
    # 如果是单端
    tophat2 -p ${threads} -o ${RESULT}/tophat2/${SRR} ${INDEX}/${bowtie2Index} ${DATA}/${fq1}
elif [ "$mode" == "paired" ];then
    # 如果是多端
    tophat2 -p ${threads} -o ${RESULT}/tophat2/${SRR} ${INDEX}/${bowtie2Index} ${DATA}/${fq1} ${DATA}/${fq2}
fi

cufflinks -p ${threads} -o ${RESULT}/cufflinks/${SRR} -G ${INDEX}/$Annotation ${RESULT}/tophat2/${SRR}/accepted_hits.bam

echo " "
echo " Running ${SRR} is compeleted."
echo " "
