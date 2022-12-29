
# 1 Upstream processing

## 1.1 Data download

slurm: [./scripts/trim_galore.sh](./scripts/trim_galore.sh)


submit: 

```bash
PROJECT=/home/u22211520038/workplace/mini
cd $PROJECT
cat $PROJECT/SRR_Acc_List.txt| while read id;do
    sbatch ~/scripts/fastq-dump.sh ${id} $PROJECT/01_rawdata/ 
done
```

## 1.2 trim_galore: adapter and quality trimming

slurm: [./scripts/trim_galore.sh](./scripts/trim_galore.sh)


submit: 

```bash
PROJECT=/home/u22211520038/workplace/mini
mode=paired

cd $PROJECT/01_rawdata/
cat $PROJECT/SRR_Acc_List.txt| while read SRR;do 
    sbatch ~/scripts/trim_galore.sh $mode $SRR
done
```

## 1.3 Prepare reference genome and annotation files

genome: 

```bash
PROJECT=/home/u22211520038/workplace/mini
mkdir -p $PROJECT/00_index
cd $PROJECT/00_index
# download
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz
# uncompress and rename
gunzip GRCh37.p13.genome.fa.gz
mkdir hg19 
mv GRCh37.p13.genome.fa hg19/hg19.fa
cd hg19/
# bowtie2-build
bowtie2-build hg19.fa hg19
```

gtf

```bash
PROJECT=/home/u22211520038/workplace/mini
cd $PROJECT/00_index
# download
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
# uncompress and rename
gunzip gencode.v19.annotation.gtf.gz
mv gencode.v19.annotation.gtf_withproteinids hg19_genes.gtf
```

## 1.4 Tophat2 & cufflinks

slurm: [./scripts/tophat2_cufflinks.sh](./scripts/tophat2_cufflinks.sh)

submit:

```bash
PROJECT=/home/u22211520038/workplace/mini
mode=paired
cd $PROJECT
cat $PROJECT/SRR_Acc_List.txt | while read SRR;do
    sbatch ~/scripts/tophat2_cufflinks.sh \
        $mode \
        4 \
        hg19_genes.gtf \
        hg19_Bowtie2Index/hg19 \
        ${SRR} \
        ${SRR}_1_val_1.fq.gz \
        ${SRR}_2_val_2.fq.gz 
done
```

‍

‍