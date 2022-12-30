# 1 Upstream processing

## 1.1 Data download

Slurm: [./scripts/trim_galore.sh](./scripts/trim_galore.sh)

Run:

```bash
PROJECT=/home/u22211520038/workplace/mini

cd $PROJECT
cat $PROJECT/00_metadata/SRR_Acc_List.txt| while read id;do
    sbatch ~/scripts/fastq-dump.sh ${id} $PROJECT/01_rawdata/ 
done
```

## 1.2 trim_galore: adapter and quality trimming

Slurm: [./scripts/trim_galore.sh](./scripts/trim_galore.sh)

Run:

```bash
PROJECT=/home/u22211520038/workplace/mini
mode=single

cd $PROJECT/01_rawdata/
cat $PROJECT/00_metadata/SRR_Acc_List.txt| while read SRR;do 
    sbatch ~/scripts/trim_galore.sh $mode $SRR
done
```

## 1.3 Prepare reference genome and annotation files

Gencode website page：[Mouse Release M8(GRCm38.p4)](https://www.gencodegenes.org/mouse/release_M8.html)

Download Genome:

```bash
PROJECT=/home/u22211520038/workplace/mini
mkdir -p $PROJECT/00_index
cd $PROJECT/00_index
# download
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M8/GRCm38.p4.genome.fa.gz
# uncompress and rename
gunzip GRCm38.p4.genome.fa.gz
mkdir -p m8_Bowtie2Index
mv GRCm38.p4.genome.fa m8_Bowtie2Index/m8.fa
cd m8_Bowtie2Index/
# bowtie2-build
bowtie2-build m8.fa m8 --threads 5
```

Download GTF:

```bash
PROJECT=/home/u22211520038/workplace/mini
cd $PROJECT/00_index

# download
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M8/gencode.vM8.annotation.gtf.gz

# uncompress and rename
gunzip gencode.vM8.annotation.gtf.gz
mv gencode.vM8.annotation.gtf m8_genes.gtf
```

## 1.4 Tophat2 & cufflinks

Slurm: [./scripts/tophat2_cufflinks.sh](./scripts/tophat2_cufflinks.sh)

Run:

```bash
PROJECT=/home/u22211520038/workplace/mini
mode=single
cd $PROJECT
cat $PROJECT/00_metadata/SRR_Acc_List.txt | while read SRR;do
    sbatch ~/scripts/tophat2_cufflinks.sh \
        $mode \
        4 \
        m8_genes.gtf \
         m8_Bowtie2Index/m8 \
        ${SRR} \
        ${SRR}_trimmed.fq.gz
done
```

## ‍1.5 FPKM to TPM

Rscript: [./scripts/FPKM2TPM.R](./scripts/FPKM2TPM.R)

Run：

```bash
PROJECT=/home/u22211520038/workplace/mini
mkdir -p $PROJECT/02_result/fpkm
mkdir -p $PROJECT/02_result/FPKM2TPM
cd $PROJECT/02_result

# Create symbolic links, and name it as a group category.
ln -s $PROJECT/02_result/cufflinks/SRR10123774/genes.fpkm_tracking $PROJECT/02_result/fpkm/Ctrl_18M_A.fpkm
ln -s $PROJECT/02_result/cufflinks/SRR10123775/genes.fpkm_tracking $PROJECT/02_result/fpkm/Ctrl_18M_B.fpkm
ln -s $PROJECT/02_result/cufflinks/SRR10123776/genes.fpkm_tracking $PROJECT/02_result/fpkm/Ctrl_18M_C.fpkm
ln -s $PROJECT/02_result/cufflinks/SRR10123786/genes.fpkm_tracking $PROJECT/02_result/fpkm/MCAO_18M_A.fpkm
ln -s $PROJECT/02_result/cufflinks/SRR10123787/genes.fpkm_tracking $PROJECT/02_result/fpkm/MCAO_18M_B.fpkm
ln -s $PROJECT/02_result/cufflinks/SRR10123788/genes.fpkm_tracking $PROJECT/02_result/fpkm/MCAO_18M_C.fpkm

# Transform FPKM into TPM through Rscript
~/scripts/FPKM2TPM.R -f \
	Ctrl_18M_A.fpkm,Ctrl_18M_B.fpkm,Ctrl_18M_C.fpkm,MCAO_18M_A.fpkm,MCAO_18M_B.fpkm,MCAO_18M_C.fpkm \
	-o $PROJECT/02_result/FPKM2TPM/Expression
```

## 1.6 Cuffdiff

Slurm：[./scripts/cuffdiff.sh](./scripts/cuffdiff.sh)

Run:

```bash
PROJECT=/home/u22211520038/workplace/mini
mkdir -p $PROJECT/02_result/bam
mkdir -p $PROJECT/02_result/cuffdiff


# Create symbolic links, and name it as a group category.
ln -s $PROJECT/02_result/tophat2/SRR10123774/accepted_hits.bam $PROJECT/02_result/bam/Ctrl_18M_A.bam
ln -s $PROJECT/02_result/tophat2/SRR10123775/accepted_hits.bam $PROJECT/02_result/bam/Ctrl_18M_B.bam
ln -s $PROJECT/02_result/tophat2/SRR10123776/accepted_hits.bam $PROJECT/02_result/bam/Ctrl_18M_C.bam
ln -s $PROJECT/02_result/tophat2/SRR10123786/accepted_hits.bam $PROJECT/02_result/bam/MCAO_18M_A.bam
ln -s $PROJECT/02_result/tophat2/SRR10123787/accepted_hits.bam $PROJECT/02_result/bam/MCAO_18M_B.bam
ln -s $PROJECT/02_result/tophat2/SRR10123788/accepted_hits.bam $PROJECT/02_result/bam/MCAO_18M_C.bam


cd $PROJECT/02_result/cuffdiff

# Run cuffdiff for differential gene analysis.
sbatch ~/scripts/cuffdiff.sh \
    $PROJECT/02_result/cuffdiff \
    $PROJECT/00_index/m8_genes.gtf \
    $PROJECT/02_result/bam \
    Ctrl \
    MCAO


```
