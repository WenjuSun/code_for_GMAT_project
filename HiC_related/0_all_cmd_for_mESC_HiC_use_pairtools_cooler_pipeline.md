

[TOC]



```bash
mamba create --name hicenv bwa-mem2 python=3.10.14 pairtools cooler=0.9.3 cooltools coolpuppy hicexplorer coolbox pygenometracks jupyter ipykernel 
mamba activate hicenv
# 
python --version
# Python 3.10.14
mamba list |grep pair
# pairix                    0.3.8           py310h1af8fb7_2    bioconda
# pairtools                 1.1.0           py310h8dfefeb_1    bioconda
mamba list |grep cool
# coolbox                   0.3.8              pyhdfd78af_0    bioconda
# cooler                    0.9.3              pyhdfd78af_0    bioconda
# coolpuppy                 1.1.0              pyh086e186_0    bioconda
# cooltools                 0.7.0           py310h581d4b6_2    bioconda
# hic2cool                  1.0.1              pyh7cba7a3_0    bioconda
mamba list |grep hic
# # packages in environment at /home/sunwenju/miniforge3/envs/hicenv:
# hic2cool                  1.0.1              pyh7cba7a3_0    bioconda
# hicexplorer               3.7.4              pyhdfd78af_0    bioconda
# hicmatrix                 17.2               pyhdfd78af_0    bioconda

```



# Data Preparation

```bash
# GSE98671 Nora mESC WT and CTCF IAA HiC data
# /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE98671_mESC_HiC_data
# GSE94452 Ren mESC WT and CTCF IAA HiC data
# /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE94452_mESC_HiC_data

# 
# part1: /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC
# part2: /data1/sunwenju/940_raw_data_backup/2024041201_WN_BPK25_HiC_FLG_ChromSeq

# 
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/
# # WN BPK25 HiC rep1/2 part1 data links
# ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/0_raw_data_download/20240305/Hic_BPK25/bpk25_hic_mesc_S0_L000_R1_000.fastq.gz WN_mESC_BPK25_HiC_rep1_R1.part1.fastq.gz
# ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/0_raw_data_download/20240305/Hic_BPK25/bpk25_hic_mesc_S0_L000_R2_000.fastq.gz WN_mESC_BPK25_HiC_rep1_R2.part1.fastq.gz
# ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/0_raw_data_download/20240305/Hic_BPK25/bpk25_hic_mesc_rep_S0_L000_R1_000.fastq.gz WN_mESC_BPK25_HiC_rep2_R1.part1.fastq.gz
# ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/0_raw_data_download/20240305/Hic_BPK25/bpk25_hic_mesc_rep_S0_L000_R2_000.fastq.gz WN_mESC_BPK25_HiC_rep2_R2.part1.fastq.gz
# # WN BPK25 HiC rep1/2 part2 data links
# ln -s /data1/sunwenju/940_raw_data_backup/2024041201_WN_BPK25_HiC_FLG_ChromSeq/0_raw_data_download/20240412/Hic_BPK25/bpk25_hic_mesc_S0_L000_R1_000.fastq.gz WN_mESC_BPK25_HiC_rep1_R1.part2.fastq.gz
# ln -s /data1/sunwenju/940_raw_data_backup/2024041201_WN_BPK25_HiC_FLG_ChromSeq/0_raw_data_download/20240412/Hic_BPK25/bpk25_hic_mesc_S0_L000_R2_000.fastq.gz WN_mESC_BPK25_HiC_rep1_R2.part2.fastq.gz
# ln -s /data1/sunwenju/940_raw_data_backup/2024041201_WN_BPK25_HiC_FLG_ChromSeq/0_raw_data_download/20240412/Hic_BPK25/bpk25_hic_mesc_rep_S0_L000_R1_000.fastq.gz WN_mESC_BPK25_HiC_rep2_R1.part2.fastq.gz
# ln -s /data1/sunwenju/940_raw_data_backup/2024041201_WN_BPK25_HiC_FLG_ChromSeq/0_raw_data_download/20240412/Hic_BPK25/bpk25_hic_mesc_rep_S0_L000_R2_000.fastq.gz WN_mESC_BPK25_HiC_rep2_R2.part2.fastq.gz
# ls -1 
# # WN BPK25 HiC rep1/2  part12 merge
# cat WN_mESC_BPK25_HiC_rep1_R1.part1.fastq.gz WN_mESC_BPK25_HiC_rep1_R1.part2.fastq.gz >WN_mESC_BPK25_HiC_rep1_R1.fastq.gz
# cat WN_mESC_BPK25_HiC_rep1_R2.part1.fastq.gz WN_mESC_BPK25_HiC_rep1_R2.part2.fastq.gz >WN_mESC_BPK25_HiC_rep1_R2.fastq.gz
# cat WN_mESC_BPK25_HiC_rep2_R1.part1.fastq.gz WN_mESC_BPK25_HiC_rep2_R1.part2.fastq.gz >WN_mESC_BPK25_HiC_rep2_R1.fastq.gz
# cat WN_mESC_BPK25_HiC_rep2_R2.part1.fastq.gz WN_mESC_BPK25_HiC_rep2_R2.part2.fastq.gz >WN_mESC_BPK25_HiC_rep2_R2.fastq.gz
# rm WN_mESC_BPK25_HiC*part[12].fastq.gz

# WN BPK25 HiC rep1/2 part12 merge fastq 
ln -s /data1/sunwenju/projects_on_940xa/2024041701_WN_mESC_BPK25_HiC_add/WN_mESC_BPK25_HiC_rep1_R1.fastq.gz WN_mESC_BPK25_HiC_rep1_R1.fastq.gz
ln -s /data1/sunwenju/projects_on_940xa/2024041701_WN_mESC_BPK25_HiC_add/WN_mESC_BPK25_HiC_rep1_R2.fastq.gz WN_mESC_BPK25_HiC_rep1_R2.fastq.gz
ln -s /data1/sunwenju/projects_on_940xa/2024041701_WN_mESC_BPK25_HiC_add/WN_mESC_BPK25_HiC_rep2_R1.fastq.gz WN_mESC_BPK25_HiC_rep2_R1.fastq.gz
ln -s /data1/sunwenju/projects_on_940xa/2024041701_WN_mESC_BPK25_HiC_add/WN_mESC_BPK25_HiC_rep2_R2.fastq.gz WN_mESC_BPK25_HiC_rep2_R2.fastq.gz

# Nora WT & IAA HiC rep1/2 data links
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE98671_mESC_HiC_data/SRR5633683_1.fastq.gz Nora_mESC_WT_HiC_rep1_R1.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE98671_mESC_HiC_data/SRR5633683_2.fastq.gz Nora_mESC_WT_HiC_rep1_R2.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE98671_mESC_HiC_data/SRR5633685_1.fastq.gz Nora_mESC_WT_HiC_rep2_R1.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE98671_mESC_HiC_data/SRR5633685_2.fastq.gz Nora_mESC_WT_HiC_rep2_R2.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE98671_mESC_HiC_data/SRR5633688_1.fastq.gz Nora_mESC_IAA_HiC_rep1_R1.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE98671_mESC_HiC_data/SRR5633688_2.fastq.gz Nora_mESC_IAA_HiC_rep1_R2.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE98671_mESC_HiC_data/SRR5633689_1.fastq.gz Nora_mESC_IAA_HiC_rep2_R1.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE98671_mESC_HiC_data/SRR5633689_2.fastq.gz Nora_mESC_IAA_HiC_rep2_R2.fastq.gz

# Ren WT & IAA HiC rep1/2 data links
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE94452_mESC_HiC_data/SRR11245481_1.fastq.gz Ren_mESC_WT_HiC_rep1_R1.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE94452_mESC_HiC_data/SRR11245481_2.fastq.gz Ren_mESC_WT_HiC_rep1_R2.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE94452_mESC_HiC_data/SRR11245482_1.fastq.gz Ren_mESC_WT_HiC_rep2_R1.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE94452_mESC_HiC_data/SRR11245482_2.fastq.gz Ren_mESC_WT_HiC_rep2_R2.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE94452_mESC_HiC_data/SRR11245484_1.fastq.gz Ren_mESC_IAA_HiC_rep1_R1.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE94452_mESC_HiC_data/SRR11245484_2.fastq.gz Ren_mESC_IAA_HiC_rep1_R2.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE94452_mESC_HiC_data/SRR11245485_1.fastq.gz Ren_mESC_IAA_HiC_rep2_R1.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE94452_mESC_HiC_data/SRR11245485_2.fastq.gz Ren_mESC_IAA_HiC_rep2_R2.fastq.gz

ls -1 *.fastq.gz
ls -1 *_R1.fastq.gz

```

## raw data fastqc

```bash
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline
mkdir 1_mESC_HiC_data_fastqc_result
nohup fastqc -outdir 1_mESC_HiC_data_fastqc_result --threads 72 *mESC*.fastq.gz >1_mESC_HiC_data_fastqc.runlog 2>&1 &
tail -f 1_mESC_HiC_data_fastqc.runlog

#
cd 1_mESC_HiC_data_fastqc_result
multiqc ./*fastqc.zip --ignore *.html --no-megaqc-upload
cd ../

```

```bash
# 5 组 并行运行，为各组创建文件夹并创建fastq数据文件连接
# 注意，创建的fastq文件连接，必须符合 *_R*.fastq* ,否则提示找不到fastq文件
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline
mkdir -p Nora_WT/fastq
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_WT/fastq
# Nora WT HiC rep1/2 data links
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE98671_mESC_HiC_data/SRR5633683_1.fastq.gz Nora_mESC_WT_HiC_rep1_R1.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE98671_mESC_HiC_data/SRR5633683_2.fastq.gz Nora_mESC_WT_HiC_rep1_R2.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE98671_mESC_HiC_data/SRR5633685_1.fastq.gz Nora_mESC_WT_HiC_rep2_R1.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE98671_mESC_HiC_data/SRR5633685_2.fastq.gz Nora_mESC_WT_HiC_rep2_R2.fastq.gz
ls -1

cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline
mkdir -p Nora_IAA/fastq
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_IAA/fastq
# Nora IAA HiC rep1/2 data links
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE98671_mESC_HiC_data/SRR5633688_1.fastq.gz Nora_mESC_IAA_HiC_rep1_R1.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE98671_mESC_HiC_data/SRR5633688_2.fastq.gz Nora_mESC_IAA_HiC_rep1_R2.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE98671_mESC_HiC_data/SRR5633689_1.fastq.gz Nora_mESC_IAA_HiC_rep2_R1.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE98671_mESC_HiC_data/SRR5633689_2.fastq.gz Nora_mESC_IAA_HiC_rep2_R2.fastq.gz
ls -1

cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline
mkdir -p Ren_WT/fastq
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_WT/fastq
# Ren WT HiC rep1/2 data links
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE94452_mESC_HiC_data/SRR11245481_1.fastq.gz Ren_mESC_WT_HiC_rep1_R1.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE94452_mESC_HiC_data/SRR11245481_2.fastq.gz Ren_mESC_WT_HiC_rep1_R2.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE94452_mESC_HiC_data/SRR11245482_1.fastq.gz Ren_mESC_WT_HiC_rep2_R1.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE94452_mESC_HiC_data/SRR11245482_2.fastq.gz Ren_mESC_WT_HiC_rep2_R2.fastq.gz
ls -1

cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline
mkdir -p Ren_IAA/fastq
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_IAA/fastq
# Ren IAA HiC rep1/2 data links
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE94452_mESC_HiC_data/SRR11245484_1.fastq.gz Ren_mESC_IAA_HiC_rep1_R1.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE94452_mESC_HiC_data/SRR11245484_2.fastq.gz Ren_mESC_IAA_HiC_rep1_R2.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE94452_mESC_HiC_data/SRR11245485_1.fastq.gz Ren_mESC_IAA_HiC_rep2_R1.fastq.gz
ln -s /data1/sunwenju/940_raw_data_backup/2024030701_WN_mESC_BPK25_HiC/GSE94452_mESC_HiC_data/SRR11245485_2.fastq.gz Ren_mESC_IAA_HiC_rep2_R2.fastq.gz
ls -1


cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline
mkdir -p WN_BPK25/fastq
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/WN_BPK25/fastq
# WN BPK25 HiC rep1/2 data links
ln -s /data1/sunwenju/projects_on_940xa/2024041701_WN_mESC_BPK25_HiC_add/WN_mESC_BPK25_HiC_rep1_R1.fastq.gz WN_mESC_BPK25_HiC_rep1_R1.fastq.gz
ln -s /data1/sunwenju/projects_on_940xa/2024041701_WN_mESC_BPK25_HiC_add/WN_mESC_BPK25_HiC_rep1_R2.fastq.gz WN_mESC_BPK25_HiC_rep1_R2.fastq.gz
ln -s /data1/sunwenju/projects_on_940xa/2024041701_WN_mESC_BPK25_HiC_add/WN_mESC_BPK25_HiC_rep2_R1.fastq.gz WN_mESC_BPK25_HiC_rep2_R1.fastq.gz
ln -s /data1/sunwenju/projects_on_940xa/2024041701_WN_mESC_BPK25_HiC_add/WN_mESC_BPK25_HiC_rep2_R2.fastq.gz WN_mESC_BPK25_HiC_rep2_R2.fastq.gz
ls -1

```



# bwa-mem2 reads mapping to mm10

> ref：https://pairtools.readthedocs.io/en/latest/examples/pairtools_walkthrough.html
>
> and https://github.com/open2c/pairtools/discussions/118
>
> 

```bash
mamba activate hicenv
# 
cd /home/sunwenju/3_genomeData/mouse/mm10
mkdir bwamem2index
cd /home/sunwenju/3_genomeData/mouse/mm10/bwamem2index
ln -s ../mm10.fa mm10.fa
ln -s ../mm10.fa.fai mm10.fa.fai
# 
nohup time bwa-mem2 index -p bwamem2idx mm10.fa >bwamem2_index_building.runlog 2>&1 &
tail -f bwamem2_index_building.runlog

```



## Nora_WT

```bash
# Nora_WT
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_WT/
nohup time bwa-mem2 mem -t 16 -SP -L 0 /home/sunwenju/3_genomeData/mouse/mm10/bwamem2index/bwamem2idx ./fastq/Nora_mESC_WT_HiC_rep1_R1.fastq.gz ./fastq/Nora_mESC_WT_HiC_rep1_R2.fastq.gz -o Nora_mESC_WT_HiC_rep1.sam >Nora_mESC_WT_HiC_rep1.bwamem2.runlog 2>&1 &
nohup time bwa-mem2 mem -t 16 -SP -L 0 /home/sunwenju/3_genomeData/mouse/mm10/bwamem2index/bwamem2idx ./fastq/Nora_mESC_WT_HiC_rep2_R1.fastq.gz ./fastq/Nora_mESC_WT_HiC_rep2_R2.fastq.gz -o Nora_mESC_WT_HiC_rep2.sam >Nora_mESC_WT_HiC_rep2.bwamem2.runlog 2>&1 &
tail -f Nora_mESC_WT_HiC_rep1.bwamem2.runlog
tail -f Nora_mESC_WT_HiC_rep2.bwamem2.runlog

```


## Nora_IAA

```bash
# Nora_IAA
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_IAA/
nohup time bwa-mem2 mem -t 16 -SP -L 0 /home/sunwenju/3_genomeData/mouse/mm10/bwamem2index/bwamem2idx ./fastq/Nora_mESC_IAA_HiC_rep1_R1.fastq.gz ./fastq/Nora_mESC_IAA_HiC_rep1_R2.fastq.gz -o Nora_mESC_IAA_HiC_rep1.sam >Nora_mESC_IAA_HiC_rep1.bwamem2.runlog 2>&1 &
nohup time bwa-mem2 mem -t 16 -SP -L 0 /home/sunwenju/3_genomeData/mouse/mm10/bwamem2index/bwamem2idx ./fastq/Nora_mESC_IAA_HiC_rep2_R1.fastq.gz ./fastq/Nora_mESC_IAA_HiC_rep2_R2.fastq.gz -o Nora_mESC_IAA_HiC_rep2.sam >Nora_mESC_IAA_HiC_rep2.bwamem2.runlog 2>&1 &
tail -f Nora_mESC_IAA_HiC_rep1.bwamem2.runlog
tail -f Nora_mESC_IAA_HiC_rep2.bwamem2.runlog

```

## Ren_WT

```bash
# Ren_WT
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_WT/
nohup time bwa-mem2 mem -t 16 -SP -L 0 /home/sunwenju/3_genomeData/mouse/mm10/bwamem2index/bwamem2idx ./fastq/Ren_mESC_WT_HiC_rep1_R1.fastq.gz ./fastq/Ren_mESC_WT_HiC_rep1_R2.fastq.gz -o Ren_mESC_WT_HiC_rep1.sam >Ren_mESC_WT_HiC_rep1.bwamem2.runlog 2>&1 &
nohup time bwa-mem2 mem -t 16 -SP -L 0 /home/sunwenju/3_genomeData/mouse/mm10/bwamem2index/bwamem2idx ./fastq/Ren_mESC_WT_HiC_rep2_R1.fastq.gz ./fastq/Ren_mESC_WT_HiC_rep2_R2.fastq.gz -o Ren_mESC_WT_HiC_rep2.sam >Ren_mESC_WT_HiC_rep2.bwamem2.runlog 2>&1 &
tail -f Ren_mESC_WT_HiC_rep1.bwamem2.runlog
tail -f Ren_mESC_WT_HiC_rep2.bwamem2.runlog

```

## Ren_IAA

```bash
# Ren_IAA
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_IAA/
nohup time bwa-mem2 mem -t 16 -SP -L 0 /home/sunwenju/3_genomeData/mouse/mm10/bwamem2index/bwamem2idx ./fastq/Ren_mESC_IAA_HiC_rep1_R1.fastq.gz ./fastq/Ren_mESC_IAA_HiC_rep1_R2.fastq.gz -o Ren_mESC_IAA_HiC_rep1.sam >Ren_mESC_IAA_HiC_rep1.bwamem2.runlog 2>&1 &
nohup time bwa-mem2 mem -t 16 -SP -L 0 /home/sunwenju/3_genomeData/mouse/mm10/bwamem2index/bwamem2idx ./fastq/Ren_mESC_IAA_HiC_rep2_R1.fastq.gz ./fastq/Ren_mESC_IAA_HiC_rep2_R2.fastq.gz -o Ren_mESC_IAA_HiC_rep2.sam >Ren_mESC_IAA_HiC_rep2.bwamem2.runlog 2>&1 &
tail -f Ren_mESC_IAA_HiC_rep1.bwamem2.runlog
tail -f Ren_mESC_IAA_HiC_rep2.bwamem2.runlog

```



## WN_BPK25

```bash
# WN_BPK25
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/WN_BPK25/
nohup time bwa-mem2 mem -t 16 -SP -L 0 /home/sunwenju/3_genomeData/mouse/mm10/bwamem2index/bwamem2idx ./fastq/WN_mESC_BPK25_HiC_rep1_R1.fastq.gz ./fastq/WN_mESC_BPK25_HiC_rep1_R2.fastq.gz -o WN_mESC_BPK25_HiC_rep1.sam >WN_mESC_BPK25_HiC_rep1.bwamem2.runlog 2>&1 &
nohup time bwa-mem2 mem -t 16 -SP -L 0 /home/sunwenju/3_genomeData/mouse/mm10/bwamem2index/bwamem2idx ./fastq/WN_mESC_BPK25_HiC_rep2_R1.fastq.gz ./fastq/WN_mESC_BPK25_HiC_rep2_R2.fastq.gz -o WN_mESC_BPK25_HiC_rep2.sam >WN_mESC_BPK25_HiC_rep2.bwamem2.runlog 2>&1 &
tail -f WN_mESC_BPK25_HiC_rep1.bwamem2.runlog
tail -f WN_mESC_BPK25_HiC_rep2.bwamem2.runlog

```



```bash
# test pairtools parse use raw sam
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_WT
mkdir pairtools_parse_test
cd pairtools_parse_test
head -n 500000 ../Nora_mESC_WT_HiC_rep1.sam >Nora_mESC_WT_HiC_rep1.head50w.sam
pairtools parse -c /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes -o test_pairs_use_raw_sam.pairs.gz Nora_mESC_WT_HiC_rep1.head50w.sam

# test sorted bam (no index)
samtools view --threads 16 -bS Nora_mESC_WT_HiC_rep1.head50w.sam | samtools sort --threads 16 -o Nora_mESC_WT_HiC_rep1.head50w.sorted.bam
pairtools parse -c /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes -o test_pairs_use_sorted_noidx_bam.pairs.gz Nora_mESC_WT_HiC_rep1.head50w.sorted.bam

# test sorted bam (with index)
samtools index Nora_mESC_WT_HiC_rep1.head50w.sorted.bam
pairtools parse -c /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes -o test_pairs_use_sorted_idx_bam.pairs.gz Nora_mESC_WT_HiC_rep1.head50w.sorted.bam


zcat test_pairs_use_raw_sam.pairs.gz |wc -l
# 208058
zcat test_pairs_use_sorted_noidx_bam.pairs.gz |wc -l
# 420998
zcat test_pairs_use_sorted_idx_bam.pairs.gz |wc -l
# 420998

#
zcat test_pairs_use_raw_sam.pairs.gz |grep -v "#" | cut -f 8 |sort |uniq -c 
  #  7927 MM
  #  9240 MR
  # 13458 MU
  #  1627 NM
  #   287 NN
  # 19202 NR
  #  1197 NU
  # 30147 RU
  # 30228 UR
  # 90500 UU
zcat test_pairs_use_sorted_noidx_bam.pairs.gz |grep -v "#" | cut -f 8 |sort |uniq -c 
#    248 MM
#     22 MR
#    593 MU
#    563 NM
#    296 NN
#  49891 NR
#   1158 NU
#    745 RU
#    779 UR
#  21305 UU
# 340352 XX

zcat test_pairs_use_sorted_idx_bam.pairs.gz |grep -v "#" | cut -f 8 |sort |uniq -c 
#    248 MM
#     22 MR
#    593 MU
#    563 NM
#    296 NN
#  49891 NR
#   1158 NU
#    745 RU
#    779 UR
#  21305 UU
# 340352 XX

```



# pairtools parse/sort/dedup/select and to cool

> restriction enzyme: 
>
> Nora ：129/Ola mESC： WT&IAA HiC （HindIII）； PE100
>
> Ren ：F123 mESC ：WT&IAA HiC （MboI）； PE51
>
> WN: F123 mESC : BPK25 HiC （DpnII）； PE150
>
> ref：https://pairtools.readthedocs.io/en/latest/parsing.html
>
> However, in complex walks (two fragments on both reads or more than two fragments on any side) you need specialized functionality that will report all the deduplicated pairs in the complex walks. This is especially relevant if you have the reads length > 100 bp, since more than 20% or all restriction fragments in the genome are then shorter than the read length. 



> https://pairtools.readthedocs.io/en/latest/cli_tools.html#pairtools-parse
>
> If the path ends with .bam, the input is decompressed from bam with samtools. By default, the input is read from stdin.
>

## Nora_WT

```bash
# 
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_WT
nohup time bash ../3_sam_to_sorted_bam_and_index.sh >3_sam_to_sorted_bam_and_index.sh.runlog 2>&1 &
tail -f 3_sam_to_sorted_bam_and_index.sh.runlog

```

### pairtools parse and sort

```bash
# pairtools parse hic pairs
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_WT
nohup time pairtools --verbose parse -c /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes -o Nora_mESC_WT_HiC_rep1.pairs.gz --assembly mm10 --min-mapq 20 --drop-seq --drop-sam --add-pair-index --add-columns mapq --output-stats Nora_mESC_WT_HiC_rep1.parse.stats.txt --walks-policy 5unique --no-flip Nora_mESC_WT_HiC_rep1.sam >Nora_mESC_WT_HiC_rep1.pairtools.parse_pairs.runlog 2>&1 &
nohup time pairtools --verbose parse -c /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes -o Nora_mESC_WT_HiC_rep2.pairs.gz --assembly mm10 --min-mapq 20 --drop-seq --drop-sam --add-pair-index --add-columns mapq --output-stats Nora_mESC_WT_HiC_rep2.parse.stats.txt --walks-policy 5unique --no-flip Nora_mESC_WT_HiC_rep2.sam >Nora_mESC_WT_HiC_rep2.pairtools.parse_pairs.runlog 2>&1 &
tail -f Nora_mESC_WT_HiC_rep1.pairtools.parse_pairs.runlog
tail -f Nora_mESC_WT_HiC_rep2.pairtools.parse_pairs.runlog

# Sort pairs using pairtools sort:
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_WT
nohup time pairtools --verbose sort --nproc 8 --tmpdir ./ -o Nora_mESC_WT_HiC_rep1.pairs.sorted.gz Nora_mESC_WT_HiC_rep1.pairs.gz >Nora_mESC_WT_HiC_rep1.pairtools.sort.runlog 2>&1 &
nohup time pairtools --verbose sort --nproc 8 --tmpdir ./ -o Nora_mESC_WT_HiC_rep2.pairs.sorted.gz Nora_mESC_WT_HiC_rep2.pairs.gz >Nora_mESC_WT_HiC_rep2.pairtools.sort.runlog 2>&1 &
tail -f Nora_mESC_WT_HiC_rep1.pairtools.sort.runlog
tail -f Nora_mESC_WT_HiC_rep2.pairtools.sort.runlog

# raw pairs stats
nohup time pairtools stats Nora_mESC_WT_HiC_rep1.pairs.gz -o Nora_mESC_WT_HiC_rep1.pairs.stats.tsv >Nora_mESC_WT_HiC_rep1.pairs.stats.runlog 2>&1 &
nohup time pairtools stats Nora_mESC_WT_HiC_rep2.pairs.gz -o Nora_mESC_WT_HiC_rep2.pairs.stats.tsv >Nora_mESC_WT_HiC_rep2.pairs.stats.runlog 2>&1 &
tail -f Nora_mESC_WT_HiC_rep1.pairs.stats.runlog
tail -f Nora_mESC_WT_HiC_rep2.pairs.stats.runlog

pairtools stats -o Nora_mESC_WT_HiC.rep12.pairs.stats.merge.tsv --merge Nora_mESC_WT_HiC_rep1.pairs.stats.tsv Nora_mESC_WT_HiC_rep2.pairs.stats.tsv
wc -l *.tsv
# 2429 Nora_mESC_WT_HiC.rep12.pairs.stats.merge.tsv
# 2323 Nora_mESC_WT_HiC_rep1.pairs.stats.tsv
# 2362 Nora_mESC_WT_HiC_rep2.pairs.stats.tsv

head *.tsv
# ==> Nora_mESC_WT_HiC.rep12.pairs.stats.merge.tsv <==
# total	706832052
# total_unmapped	50428003
# total_single_sided_mapped	159955118
# total_mapped	496448931
# total_dups	0
# total_nodups	496448931
# cis	387879367
# trans	108569564
# pair_types/NU	982109
# pair_types/NM	3570807

# ==> Nora_mESC_WT_HiC_rep1.pairs.stats.tsv <==
# total	359216310
# total_unmapped	22914648
# total_single_sided_mapped	82303444
# total_mapped	253998218
# total_dups	0
# total_nodups	253998218
# cis	195649217
# trans	58349001
# pair_types/UU	154048101
# pair_types/RU	49941044

# ==> Nora_mESC_WT_HiC_rep2.pairs.stats.tsv <==
# total	347615742
# total_unmapped	27513355
# total_single_sided_mapped	77651674
# total_mapped	242450713
# total_dups	0
# total_nodups	242450713
# cis	192230150
# trans	50220563
# pair_types/UU	154079122
# pair_types/RU	44887394

```

### pairtools  dedup

```bash
## Detect and remove duplicates using pairtools dedup and generate statistics
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_WT
nohup time pairtools --verbose dedup -p 16 \
 --output Nora_mESC_WT_HiC_rep1.nodups.pairs.gz \
 --output-dups Nora_mESC_WT_HiC_rep1.dups.pairs.gz \
 --output-unmapped Nora_mESC_WT_HiC_rep1.unmapped.pairs.gz \
 --output-stats Nora_mESC_WT_HiC_rep1.dedup.stats.txt \
 Nora_mESC_WT_HiC_rep1.pairs.sorted.gz \
 >Nora_mESC_WT_HiC_rep1.pairtools.dedup.runlog 2>&1 &
nohup time pairtools --verbose dedup -p 16 \
 --output Nora_mESC_WT_HiC_rep2.nodups.pairs.gz \
 --output-dups Nora_mESC_WT_HiC_rep2.dups.pairs.gz \
 --output-unmapped Nora_mESC_WT_HiC_rep2.unmapped.pairs.gz \
 --output-stats Nora_mESC_WT_HiC_rep2.dedup.stats.txt \
 Nora_mESC_WT_HiC_rep2.pairs.sorted.gz \
 >Nora_mESC_WT_HiC_rep2.pairtools.dedup.runlog 2>&1 &
tail -f Nora_mESC_WT_HiC_rep1.pairtools.dedup.runlog
tail -f Nora_mESC_WT_HiC_rep2.pairtools.dedup.runlog

```

```bash
# pairs type stat
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_WT
# zcat Nora_mESC_WT_HiC_rep1.pairs.gz |grep -v "#" |cut -f 8 |head -n 100 |sort |uniq -c |sort -k1,1nr
zcat Nora_mESC_WT_HiC_rep1.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 154048101 UU
#  50009073 UR
#  49941044 RU
#  18728296 MM
#  16832653 RN
#  15764303 NR
#  14426043 UM
#  13993452 MU
#   9498959 RM
#   9448511 MR
#   2020751 MN
#   1869231 UN
#   1722288 NM
#    470292 NU
#    443313 NN
zcat Nora_mESC_WT_HiC_rep1.pairs.sorted.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 154048101 UU
#  50009073 UR
#  49941044 RU
#  18728296 MM
#  16832653 RN
#  15764303 NR
#  14426043 UM
#  13993452 MU
#   9498959 RM
#   9448511 MR
#   2020751 MN
#   1869231 UN
#   1722288 NM
#    470292 NU
#    443313 NN
zcat Nora_mESC_WT_HiC_rep1.dups.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 34360766 DD
zcat Nora_mESC_WT_HiC_rep1.nodups.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 129296779 UU
#  45196597 UR
#  45144076 RU

zcat Nora_mESC_WT_HiC_rep2.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 154079122 UU
#  44887394 RU
#  43484197 UR
#  22080624 MM
#  15108590 RN
#  14783028 NR
#  13985857 UM
#  13602597 MU
#   8900262 RM
#   8552215 MR
#   2207308 UN
#   2104864 MN
#   1848519 NM
#   1479348 NN
#    511817 NU
zcat Nora_mESC_WT_HiC_rep2.pairs.sorted.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 154079122 UU
#  44887394 RU
#  43484197 UR
#  22080624 MM
#  15108590 RN
#  14783028 NR
#  13985857 UM
#  13602597 MU
#   8900262 RM
#   8552215 MR
#   2207308 UN
#   2104864 MN
#   1848519 NM
#   1479348 NN
#    511817 NU
zcat Nora_mESC_WT_HiC_rep2.dups.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 37633584 DD
zcat Nora_mESC_WT_HiC_rep2.nodups.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 127493588 UU
#  39273913 RU
#  38049628 UR

```

### cooler cload pairs (pairs to cool)

```bash
# ref https://pairtools.readthedocs.io/en/latest/protocols_pipelines.html
# Recommended pairtools parameters for standard Hi-C protocols
# https://pairtools.readthedocs.io/en/latest/examples/pairtools_walkthrough.html
# select MAPQ1/2 >= 30 
nohup time pairtools --verbose select "(mapq1>=30) and (mapq2>=30)" Nora_mESC_WT_HiC_rep1.nodups.pairs.gz -o Nora_mESC_WT_HiC_rep1.nodups.pairs.Qge30.gz >Nora_mESC_WT_HiC_rep1.nodups.pairs.select_Qge30.runlog 2>&1 &
nohup time pairtools --verbose select "(mapq1>=30) and (mapq2>=30)" Nora_mESC_WT_HiC_rep2.nodups.pairs.gz -o Nora_mESC_WT_HiC_rep2.nodups.pairs.Qge30.gz >Nora_mESC_WT_HiC_rep2.nodups.pairs.select_Qge30.runlog 2>&1 &
tail -f Nora_mESC_WT_HiC_rep1.nodups.pairs.select_Qge30.runlog
tail -f Nora_mESC_WT_HiC_rep2.nodups.pairs.select_Qge30.runlog

# pairs to 1k cool
nohup time cooler cload pairs --assembly mm10 -c1 2 -p1 3 -c2 4 -p2 5 --zero-based --temp-dir . ../mm10.chrom.sizes:1000 Nora_mESC_WT_HiC_rep1.nodups.pairs.Qge30.gz Nora_mESC_WT_HiC_rep1.nodups.pairs.1k.cool >Nora_mESC_WT_HiC_rep1.nodups.pairs.cooler_cload.runlog 2>&1 &
nohup time cooler cload pairs --assembly mm10 -c1 2 -p1 3 -c2 4 -p2 5 --zero-based --temp-dir . ../mm10.chrom.sizes:1000 Nora_mESC_WT_HiC_rep2.nodups.pairs.Qge30.gz Nora_mESC_WT_HiC_rep2.nodups.pairs.1k.cool >Nora_mESC_WT_HiC_rep2.nodups.pairs.cooler_cload.runlog 2>&1 &
tail -f Nora_mESC_WT_HiC_rep1.nodups.pairs.cooler_cload.runlog
tail -f Nora_mESC_WT_HiC_rep2.nodups.pairs.cooler_cload.runlog

# rep1/2 raw 1k cool merge
nohup time cooler merge Nora_mESC_WT_HiC.rep12.merged.1k.cool Nora_mESC_WT_HiC_rep1.nodups.pairs.1k.cool Nora_mESC_WT_HiC_rep2.nodups.pairs.1k.cool >Nora_mESC_WT_HiC.rep12.1k.cool.merge.runlog 2>&1 &
tail -f Nora_mESC_WT_HiC.rep12.1k.cool.merge.runlog

cooler info Nora_mESC_WT_HiC_rep1.nodups.pairs.1k.cool
# {
#     "bin-size": 1000,
#     "bin-type": "fixed",
#     "creation-date": "2024-06-12T19:55:51.460724",
#     "format": "HDF5::Cooler",
#     "format-url": "https://github.com/open2c/cooler",
#     "format-version": 3,
#     "generated-by": "cooler-0.10.0",
#     "genome-assembly": "mm10",
#     "nbins": 2725532,
#     "nchroms": 21,
#     "nnz": 156056115,
#     "storage-mode": "symmetric-upper",
#     "sum": 213472977
# }


cooler info Nora_mESC_WT_HiC_rep2.nodups.pairs.1k.cool
# {
#     "bin-size": 1000,
#     "bin-type": "fixed",
#     "creation-date": "2024-06-12T19:55:11.649771",
#     "format": "HDF5::Cooler",
#     "format-url": "https://github.com/open2c/cooler",
#     "format-version": 3,
#     "generated-by": "cooler-0.10.0",
#     "genome-assembly": "mm10",
#     "nbins": 2725532,
#     "nchroms": 21,
#     "nnz": 129716272,
#     "storage-mode": "symmetric-upper",
#     "sum": 198724615
# }

cooler info Nora_mESC_WT_HiC.rep12.merged.1k.cool
# {
#     "bin-size": 1000,
#     "bin-type": "fixed",
#     "creation-date": "2024-06-12T20:00:13.248513",
#     "format": "HDF5::Cooler",
#     "format-url": "https://github.com/open2c/cooler",
#     "format-version": 3,
#     "generated-by": "cooler-0.10.0",
#     "genome-assembly": "mm10",
#     "nbins": 2725532,
#     "nchroms": 21,
#     "nnz": 267316479,
#     "storage-mode": "symmetric-upper",
#     "sum": 412197592
# }

```

### cooler zoomify cool to mcool

```bash
# zoomify 1k cool to mcool
# -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000
nohup time cooler zoomify -p 16 -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000 --balance -o Nora_mESC_WT_HiC_rep1.nodups.pairs.1k.zoomify.mcool Nora_mESC_WT_HiC_rep1.nodups.pairs.1k.cool >Nora_mESC_WT_HiC_rep1.nodups.pairs.1k.cool.zoomify.runlog 2>&1 &
nohup time cooler zoomify -p 16 -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000 --balance -o Nora_mESC_WT_HiC_rep2.nodups.pairs.1k.zoomify.mcool Nora_mESC_WT_HiC_rep2.nodups.pairs.1k.cool >Nora_mESC_WT_HiC_rep2.nodups.pairs.1k.cool.zoomify.runlog 2>&1 &
tail -f Nora_mESC_WT_HiC_rep1.nodups.pairs.1k.cool.zoomify.runlog
tail -f Nora_mESC_WT_HiC_rep2.nodups.pairs.1k.cool.zoomify.runlog

# zoomify 1k cool to mcool for rep12 merged 1k cool
# -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000
nohup time cooler zoomify -p 16 -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000 --balance -o Nora_mESC_WT_HiC.rep12.merged.1k.zoomify.mcool Nora_mESC_WT_HiC.rep12.merged.1k.cool >Nora_mESC_WT_HiC.rep12.merged.1k.cool.zoomify.runlog 2>&1 &
tail -f Nora_mESC_WT_HiC.rep12.merged.1k.cool.zoomify.runlog

# for i in `ls -1 *.mcool`;do echo $i; cooler info $i"::/resolutions/1000";done

```



## Nora_IAA

```bash
# 
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_IAA
nohup time bash ../3_sam_to_sorted_bam_and_index.sh >3_sam_to_sorted_bam_and_index.sh.runlog 2>&1 &
tail -f 3_sam_to_sorted_bam_and_index.sh.runlog

```

### pairtools parse and sort
```bash
# pairtools parse hic pairs
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_IAA
nohup time pairtools --verbose parse -c /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes -o Nora_mESC_IAA_HiC_rep1.pairs.gz --assembly mm10 --min-mapq 20 --drop-seq --drop-sam --add-pair-index --add-columns mapq --output-stats Nora_mESC_IAA_HiC_rep1.parse.stats.txt --walks-policy 5unique --no-flip Nora_mESC_IAA_HiC_rep1.sam >Nora_mESC_IAA_HiC_rep1.pairtools.parse_pairs.runlog 2>&1 &
nohup time pairtools --verbose parse -c /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes -o Nora_mESC_IAA_HiC_rep2.pairs.gz --assembly mm10 --min-mapq 20 --drop-seq --drop-sam --add-pair-index --add-columns mapq --output-stats Nora_mESC_IAA_HiC_rep2.parse.stats.txt --walks-policy 5unique --no-flip Nora_mESC_IAA_HiC_rep2.sam >Nora_mESC_IAA_HiC_rep2.pairtools.parse_pairs.runlog 2>&1 &
tail -f Nora_mESC_IAA_HiC_rep1.pairtools.parse_pairs.runlog
tail -f Nora_mESC_IAA_HiC_rep2.pairtools.parse_pairs.runlog

# Sort pairs using pairtools sort:
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_IAA
nohup time pairtools --verbose sort --nproc 8 --tmpdir ./ -o Nora_mESC_IAA_HiC_rep1.pairs.sorted.gz Nora_mESC_IAA_HiC_rep1.pairs.gz >Nora_mESC_IAA_HiC_rep1.pairtools.sort.runlog 2>&1 &
nohup time pairtools --verbose sort --nproc 8 --tmpdir ./ -o Nora_mESC_IAA_HiC_rep2.pairs.sorted.gz Nora_mESC_IAA_HiC_rep2.pairs.gz >Nora_mESC_IAA_HiC_rep2.pairtools.sort.runlog 2>&1 &
tail -f Nora_mESC_IAA_HiC_rep1.pairtools.sort.runlog
tail -f Nora_mESC_IAA_HiC_rep2.pairtools.sort.runlog

# raw pairs stats
nohup time pairtools stats Nora_mESC_IAA_HiC_rep1.pairs.gz -o Nora_mESC_IAA_HiC_rep1.pairs.stats.tsv >Nora_mESC_IAA_HiC_rep1.pairs.stats.runlog 2>&1 &
nohup time pairtools stats Nora_mESC_IAA_HiC_rep2.pairs.gz -o Nora_mESC_IAA_HiC_rep2.pairs.stats.tsv >Nora_mESC_IAA_HiC_rep2.pairs.stats.runlog 2>&1 &
tail -f Nora_mESC_IAA_HiC_rep1.pairs.stats.runlog
tail -f Nora_mESC_IAA_HiC_rep2.pairs.stats.runlog

pairtools stats -o Nora_mESC_IAA_HiC.rep12.pairs.stats.merge.tsv --merge Nora_mESC_IAA_HiC_rep1.pairs.stats.tsv Nora_mESC_IAA_HiC_rep2.pairs.stats.tsv
wc -l *.tsv
# 2390 Nora_mESC_IAA_HiC.rep12.pairs.stats.merge.tsv
# 2333 Nora_mESC_IAA_HiC_rep1.pairs.stats.tsv
# 2265 Nora_mESC_IAA_HiC_rep2.pairs.stats.tsv

head *.tsv
# ==> Nora_mESC_IAA_HiC.rep12.pairs.stats.merge.tsv <==
# total	548506141
# total_unmapped	35763547
# total_single_sided_mapped	126753598
# total_mapped	385988996
# total_dups	0
# total_nodups	385988996
# cis	283929114
# trans	102059882
# pair_types/RU	72268788
# pair_types/UU	241111854

# ==> Nora_mESC_IAA_HiC_rep1.pairs.stats.tsv <==
# total	364742417
# total_unmapped	23751362
# total_single_sided_mapped	83236588
# total_mapped	257754467
# total_dups	0
# total_nodups	257754467
# cis	189612210
# trans	68142257
# pair_types/UU	160465019
# pair_types/UR	48865996

# ==> Nora_mESC_IAA_HiC_rep2.pairs.stats.tsv <==
# total	183763724
# total_unmapped	12012185
# total_single_sided_mapped	43517010
# total_mapped	128234529
# total_dups	0
# total_nodups	128234529
# cis	94316904
# trans	33917625
# pair_types/UU	80646835
# pair_types/RU	23845336


```

### pairtools dedup
```bash
## Detect and remove duplicates using pairtools dedup and generate statistics
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_IAA
nohup time pairtools --verbose dedup -p 16 \
 --output Nora_mESC_IAA_HiC_rep1.nodups.pairs.gz \
 --output-dups Nora_mESC_IAA_HiC_rep1.dups.pairs.gz \
 --output-unmapped Nora_mESC_IAA_HiC_rep1.unmapped.pairs.gz \
 --output-stats Nora_mESC_IAA_HiC_rep1.dedup.stats.txt \
 Nora_mESC_IAA_HiC_rep1.pairs.sorted.gz \
 >Nora_mESC_IAA_HiC_rep1.pairtools.dedup.runlog 2>&1 &
nohup time pairtools --verbose dedup -p 16 \
 --output Nora_mESC_IAA_HiC_rep2.nodups.pairs.gz \
 --output-dups Nora_mESC_IAA_HiC_rep2.dups.pairs.gz \
 --output-unmapped Nora_mESC_IAA_HiC_rep2.unmapped.pairs.gz \
 --output-stats Nora_mESC_IAA_HiC_rep2.dedup.stats.txt \
 Nora_mESC_IAA_HiC_rep2.pairs.sorted.gz \
 >Nora_mESC_IAA_HiC_rep2.pairtools.dedup.runlog 2>&1 &
tail -f Nora_mESC_IAA_HiC_rep1.pairtools.dedup.runlog
tail -f Nora_mESC_IAA_HiC_rep2.pairtools.dedup.runlog

```

```bash
# pairs type stat
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_IAA
# zcat Nora_mESC_IAA_HiC_rep1.pairs.gz |grep -v "#" |cut -f 8 |head -n 100 |sort |uniq -c |sort -k1,1nr
zcat Nora_mESC_IAA_HiC_rep1.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 160465019 UU
#  48865996 UR
#  48423452 RU
#  19528635 MM
#  16398447 RN
#  15578058 UM
#  15120687 MU
#  14923622 NR
#   9424572 RM
#   9422909 MR
#   2054574 MN
#   1878143 UN
#   1702233 NM
#    490150 NU
#    465920 NN
zcat Nora_mESC_IAA_HiC_rep1.pairs.sorted.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 160465019 UU
#  48865996 UR
#  48423452 RU
#  19528635 MM
#  16398447 RN
#  15578058 UM
#  15120687 MU
#  14923622 NR
#   9424572 RM
#   9422909 MR
#   2054574 MN
#   1878143 UN
#   1702233 NM
#    490150 NU
#    465920 NN
zcat Nora_mESC_IAA_HiC_rep1.dups.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 41326218 DD
zcat Nora_mESC_IAA_HiC_rep1.nodups.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 132193134 UU
#  42316450 UR
#  41918665 RU

zcat Nora_mESC_IAA_HiC_rep2.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 80646835 UU
# 23845336 RU
# 23742358 UR
#  9768643 MM
#  8185821 NR
#  8010221 RN
#  7864582 MU
#  7818054 UM
#  4602388 RM
#  4593704 MR
#  1349888 NU
#  1092352 UN
#  1077104 NM
#  1026293 MN
#   140145 NN
zcat Nora_mESC_IAA_HiC_rep2.pairs.sorted.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 80646835 UU
# 23845336 RU
# 23742358 UR
#  9768643 MM
#  8185821 NR
#  8010221 RN
#  7864582 MU
#  7818054 UM
#  4602388 RM
#  4593704 MR
#  1349888 NU
#  1092352 UN
#  1077104 NM
#  1026293 MN
#   140145 NN
zcat Nora_mESC_IAA_HiC_rep2.dups.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 3653965 DD
zcat Nora_mESC_IAA_HiC_rep2.nodups.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 77192742 UU
# 23744578 RU
# 23643244 UR

```

### cooler cload pairs (pairs to cool)

```bash
# ref https://pairtools.readthedocs.io/en/latest/protocols_pipelines.html
# Recommended pairtools parameters for standard Hi-C protocols
# https://pairtools.readthedocs.io/en/latest/examples/pairtools_walkthrough.html
# select MAPQ1/2 >= 30 
mamba activate hicenv
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_IAA
nohup time pairtools --verbose select "(mapq1>=30) and (mapq2>=30)" Nora_mESC_IAA_HiC_rep1.nodups.pairs.gz -o Nora_mESC_IAA_HiC_rep1.nodups.pairs.Qge30.gz >Nora_mESC_IAA_HiC_rep1.nodups.pairs.select_Qge30.runlog 2>&1 &
nohup time pairtools --verbose select "(mapq1>=30) and (mapq2>=30)" Nora_mESC_IAA_HiC_rep2.nodups.pairs.gz -o Nora_mESC_IAA_HiC_rep2.nodups.pairs.Qge30.gz >Nora_mESC_IAA_HiC_rep2.nodups.pairs.select_Qge30.runlog 2>&1 &
tail -f Nora_mESC_IAA_HiC_rep1.nodups.pairs.select_Qge30.runlog
tail -f Nora_mESC_IAA_HiC_rep2.nodups.pairs.select_Qge30.runlog

# pairs to 1k cool
nohup time cooler cload pairs --assembly mm10 -c1 2 -p1 3 -c2 4 -p2 5 --zero-based --temp-dir . ../mm10.chrom.sizes:1000 Nora_mESC_IAA_HiC_rep1.nodups.pairs.Qge30.gz Nora_mESC_IAA_HiC_rep1.nodups.pairs.1k.cool >Nora_mESC_IAA_HiC_rep1.nodups.pairs.cooler_cload.runlog 2>&1 &
nohup time cooler cload pairs --assembly mm10 -c1 2 -p1 3 -c2 4 -p2 5 --zero-based --temp-dir . ../mm10.chrom.sizes:1000 Nora_mESC_IAA_HiC_rep2.nodups.pairs.Qge30.gz Nora_mESC_IAA_HiC_rep2.nodups.pairs.1k.cool >Nora_mESC_IAA_HiC_rep2.nodups.pairs.cooler_cload.runlog 2>&1 &
tail -f Nora_mESC_IAA_HiC_rep1.nodups.pairs.cooler_cload.runlog
tail -f Nora_mESC_IAA_HiC_rep2.nodups.pairs.cooler_cload.runlog

# rep1/2 raw 1k cool merge
nohup time cooler merge Nora_mESC_IAA_HiC.rep12.merged.1k.cool Nora_mESC_IAA_HiC_rep1.nodups.pairs.1k.cool Nora_mESC_IAA_HiC_rep2.nodups.pairs.1k.cool >Nora_mESC_IAA_HiC.rep12.1k.cool.merge.runlog 2>&1 &
tail -f Nora_mESC_IAA_HiC.rep12.1k.cool.merge.runlog

cooler info Nora_mESC_IAA_HiC_rep1.nodups.pairs.1k.cool
# {
#     "bin-size": 1000,
#     "bin-type": "fixed",
#     "creation-date": "2024-06-12T20:17:04.131914",
#     "format": "HDF5::Cooler",
#     "format-url": "https://github.com/open2c/cooler",
#     "format-version": 3,
#     "generated-by": "cooler-0.10.0",
#     "genome-assembly": "mm10",
#     "nbins": 2725532,
#     "nchroms": 21,
#     "nnz": 155395818,
#     "storage-mode": "symmetric-upper",
#     "sum": 210371683
# }

cooler info Nora_mESC_IAA_HiC_rep2.nodups.pairs.1k.cool
# {
#     "bin-size": 1000,
#     "bin-type": "fixed",
#     "creation-date": "2024-06-12T20:13:20.521425",
#     "format": "HDF5::Cooler",
#     "format-url": "https://github.com/open2c/cooler",
#     "format-version": 3,
#     "generated-by": "cooler-0.10.0",
#     "genome-assembly": "mm10",
#     "nbins": 2725532,
#     "nchroms": 21,
#     "nnz": 91796053,
#     "storage-mode": "symmetric-upper",
#     "sum": 121061264
# }

cooler info Nora_mESC_IAA_HiC.rep12.merged.1k.cool
# {
#     "bin-size": 1000,
#     "bin-type": "fixed",
#     "creation-date": "2024-06-12T20:34:11.082776",
#     "format": "HDF5::Cooler",
#     "format-url": "https://github.com/open2c/cooler",
#     "format-version": 3,
#     "generated-by": "cooler-0.10.0",
#     "genome-assembly": "mm10",
#     "nbins": 2725532,
#     "nchroms": 21,
#     "nnz": 233818855,
#     "storage-mode": "symmetric-upper",
#     "sum": 331432947
# }

```

### cooler zoomify cool to mcool

```bash
# zoomify 1k cool to mcool
# -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000
nohup time cooler zoomify -p 16 -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000 --balance -o Nora_mESC_IAA_HiC_rep1.nodups.pairs.1k.zoomify.mcool Nora_mESC_IAA_HiC_rep1.nodups.pairs.1k.cool >Nora_mESC_IAA_HiC_rep1.nodups.pairs.1k.cool.zoomify.runlog 2>&1 &
nohup time cooler zoomify -p 16 -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000 --balance -o Nora_mESC_IAA_HiC_rep2.nodups.pairs.1k.zoomify.mcool Nora_mESC_IAA_HiC_rep2.nodups.pairs.1k.cool >Nora_mESC_IAA_HiC_rep2.nodups.pairs.1k.cool.zoomify.runlog 2>&1 &
tail -f Nora_mESC_IAA_HiC_rep1.nodups.pairs.1k.cool.zoomify.runlog
tail -f Nora_mESC_IAA_HiC_rep2.nodups.pairs.1k.cool.zoomify.runlog

# zoomify 1k cool to mcool for rep12 merged 1k cool
# -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000
nohup time cooler zoomify -p 16 -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000 --balance -o Nora_mESC_IAA_HiC.rep12.merged.1k.zoomify.mcool Nora_mESC_IAA_HiC.rep12.merged.1k.cool >Nora_mESC_IAA_HiC.rep12.merged.1k.cool.zoomify.runlog 2>&1 &
tail -f Nora_mESC_IAA_HiC.rep12.merged.1k.cool.zoomify.runlog

# for i in `ls -1 *.mcool`;do echo $i; cooler info $i"::/resolutions/1000";done

```



## Ren_WT

```bash
#
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_WT
nohup time bash ../3_sam_to_sorted_bam_and_index.sh >3_sam_to_sorted_bam_and_index.sh.runlog 2>&1 &
tail -f 3_sam_to_sorted_bam_and_index.sh.runlog

```

### pairtools parse and sort
```bash
# pairtools parse hic pairs
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_WT
nohup time pairtools --verbose parse -c /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes -o Ren_mESC_WT_HiC_rep1.pairs.gz --assembly mm10 --min-mapq 20 --drop-seq --drop-sam --add-pair-index --add-columns mapq --output-stats Ren_mESC_WT_HiC_rep1.parse.stats.txt --walks-policy 5unique --no-flip Ren_mESC_WT_HiC_rep1.sam >Ren_mESC_WT_HiC_rep1.pairtools.parse_pairs.runlog 2>&1 &
nohup time pairtools --verbose parse -c /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes -o Ren_mESC_WT_HiC_rep2.pairs.gz --assembly mm10 --min-mapq 20 --drop-seq --drop-sam --add-pair-index --add-columns mapq --output-stats Ren_mESC_WT_HiC_rep2.parse.stats.txt --walks-policy 5unique --no-flip Ren_mESC_WT_HiC_rep2.sam >Ren_mESC_WT_HiC_rep2.pairtools.parse_pairs.runlog 2>&1 &
tail -f Ren_mESC_WT_HiC_rep1.pairtools.parse_pairs.runlog
tail -f Ren_mESC_WT_HiC_rep2.pairtools.parse_pairs.runlog

# Sort pairs using pairtools sort:
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_WT
nohup time pairtools --verbose sort --nproc 8 --tmpdir ./ -o Ren_mESC_WT_HiC_rep1.pairs.sorted.gz Ren_mESC_WT_HiC_rep1.pairs.gz >Ren_mESC_WT_HiC_rep1.pairtools.sort.runlog 2>&1 &
nohup time pairtools --verbose sort --nproc 8 --tmpdir ./ -o Ren_mESC_WT_HiC_rep2.pairs.sorted.gz Ren_mESC_WT_HiC_rep2.pairs.gz >Ren_mESC_WT_HiC_rep2.pairtools.sort.runlog 2>&1 &
tail -f Ren_mESC_WT_HiC_rep1.pairtools.sort.runlog
tail -f Ren_mESC_WT_HiC_rep2.pairtools.sort.runlog

# raw pairs stats
nohup time pairtools stats Ren_mESC_WT_HiC_rep1.pairs.gz -o Ren_mESC_WT_HiC_rep1.pairs.stats.tsv >Ren_mESC_WT_HiC_rep1.pairs.stats.runlog 2>&1 &
nohup time pairtools stats Ren_mESC_WT_HiC_rep2.pairs.gz -o Ren_mESC_WT_HiC_rep2.pairs.stats.tsv >Ren_mESC_WT_HiC_rep2.pairs.stats.runlog 2>&1 &
tail -f Ren_mESC_WT_HiC_rep1.pairs.stats.runlog
tail -f Ren_mESC_WT_HiC_rep2.pairs.stats.runlog

pairtools stats -o Ren_mESC_WT_HiC.rep12.pairs.stats.merge.tsv --merge Ren_mESC_WT_HiC_rep1.pairs.stats.tsv Ren_mESC_WT_HiC_rep2.pairs.stats.tsv
wc -l *.tsv
# 2559 Ren_mESC_WT_HiC.rep12.pairs.stats.merge.tsv
# 2449 Ren_mESC_WT_HiC_rep1.pairs.stats.tsv
# 2489 Ren_mESC_WT_HiC_rep2.pairs.stats.tsv

head *.tsv
# ==> Ren_mESC_WT_HiC.rep12.pairs.stats.merge.tsv <==
# total	2444199876
# total_unmapped	329359071
# total_single_sided_mapped	708081624
# total_mapped	1406759181
# total_dups	0
# total_nodups	1406759181
# cis	1245904141
# trans	160855040
# pair_types/MN	37042382
# pair_types/RM	140075

# ==> Ren_mESC_WT_HiC_rep1.pairs.stats.tsv <==
# total	1297636490
# total_unmapped	219715215
# total_single_sided_mapped	358929941
# total_mapped	718991334
# total_dups	0
# total_nodups	718991334
# cis	644768590
# trans	74222744
# pair_types/UU	718919702
# pair_types/MM	179514913

# ==> Ren_mESC_WT_HiC_rep2.pairs.stats.tsv <==
# total	1146563386
# total_unmapped	109643856
# total_single_sided_mapped	349151683
# total_mapped	687767847
# total_dups	0
# total_nodups	687767847
# cis	601135551
# trans	86632296
# pair_types/UU	687758792
# pair_types/MU	122895799

```

### pairtools dedup
```bash
## Detect and remove duplicates using pairtools dedup and generate statistics
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_WT
nohup time pairtools --verbose dedup -p 16 \
 --output Ren_mESC_WT_HiC_rep1.nodups.pairs.gz \
 --output-dups Ren_mESC_WT_HiC_rep1.dups.pairs.gz \
 --output-unmapped Ren_mESC_WT_HiC_rep1.unmapped.pairs.gz \
 --output-stats Ren_mESC_WT_HiC_rep1.dedup.stats.txt \
 Ren_mESC_WT_HiC_rep1.pairs.sorted.gz \
 >Ren_mESC_WT_HiC_rep1.pairtools.dedup.runlog 2>&1 &
nohup time pairtools --verbose dedup -p 16 \
 --output Ren_mESC_WT_HiC_rep2.nodups.pairs.gz \
 --output-dups Ren_mESC_WT_HiC_rep2.dups.pairs.gz \
 --output-unmapped Ren_mESC_WT_HiC_rep2.unmapped.pairs.gz \
 --output-stats Ren_mESC_WT_HiC_rep2.dedup.stats.txt \
 Ren_mESC_WT_HiC_rep2.pairs.sorted.gz \
 >Ren_mESC_WT_HiC_rep2.pairtools.dedup.runlog 2>&1 &
tail -f Ren_mESC_WT_HiC_rep1.pairtools.dedup.runlog
tail -f Ren_mESC_WT_HiC_rep2.pairtools.dedup.runlog

```

```bash
# pairs type stat
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_WT
# zcat Ren_mESC_WT_HiC_rep1.pairs.gz |grep -v "#" |cut -f 8 |head -n 100 |sort |uniq -c |sort -k1,1nr
zcat Ren_mESC_WT_HiC_rep1.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 718919702 UU
# 179514913 MM
# 127866001 UM
# 124115057 MU
#  70441819 UN
#  33239150 NU
#  23025061 MN
#  11167107 NM
#   6008134 NN
#   1842884 RN
#   1298979 NR
#     85757 RM
#     67101 RU
#     40294 MR
#      4531 UR
zcat Ren_mESC_WT_HiC_rep1.pairs.sorted.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 718919702 UU
# 179514913 MM
# 127866001 UM
# 124115057 MU
#  70441819 UN
#  33239150 NU
#  23025061 MN
#  11167107 NM
#   6008134 NN
#   1842884 RN
#   1298979 NR
#     85757 RM
#     67101 RU
#     40294 MR
#      4531 UR
zcat Ren_mESC_WT_HiC_rep1.dups.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 293099412 DD
zcat Ren_mESC_WT_HiC_rep1.nodups.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 425853040 UU
#     36079 RU
#      2803 UR


zcat Ren_mESC_WT_HiC_rep2.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 687758792 UU
# 125499324 UM
# 122895799 MU
#  82185914 MM
#  57497333 UN
#  39894464 NU
#  14017321 MN
#   9728649 NM
#   3711972 NN
#   1677703 RN
#   1590649 NR
#     54318 RM
#     42093 MR
#      7863 RU
#      1192 UR
zcat Ren_mESC_WT_HiC_rep2.pairs.sorted.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 687758792 UU
# 125499324 UM
# 122895799 MU
#  82185914 MM
#  57497333 UN
#  39894464 NU
#  14017321 MN
#   9728649 NM
#   3711972 NN
#   1677703 RN
#   1590649 NR
#     54318 RM
#     42093 MR
#      7863 RU
#      1192 UR
zcat Ren_mESC_WT_HiC_rep2.dups.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 62546820 DD
zcat Ren_mESC_WT_HiC_rep2.nodups.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 625212729 UU
#      7210 RU
#      1088 UR

```

### cooler cload pairs (pairs to cool)

```bash
# ref https://pairtools.readthedocs.io/en/latest/protocols_pipelines.html
# Recommended pairtools parameters for standard Hi-C protocols
# https://pairtools.readthedocs.io/en/latest/examples/pairtools_walkthrough.html
# select MAPQ1/2 >= 30 
mamba activate hicenv
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_WT
nohup time pairtools --verbose select "(mapq1>=30) and (mapq2>=30)" Ren_mESC_WT_HiC_rep1.nodups.pairs.gz -o Ren_mESC_WT_HiC_rep1.nodups.pairs.Qge30.gz >Ren_mESC_WT_HiC_rep1.nodups.pairs.select_Qge30.runlog 2>&1 &
nohup time pairtools --verbose select "(mapq1>=30) and (mapq2>=30)" Ren_mESC_WT_HiC_rep2.nodups.pairs.gz -o Ren_mESC_WT_HiC_rep2.nodups.pairs.Qge30.gz >Ren_mESC_WT_HiC_rep2.nodups.pairs.select_Qge30.runlog 2>&1 &
tail -f Ren_mESC_WT_HiC_rep1.nodups.pairs.select_Qge30.runlog
tail -f Ren_mESC_WT_HiC_rep2.nodups.pairs.select_Qge30.runlog

# pairs to 1k cool
nohup time cooler cload pairs --assembly mm10 -c1 2 -p1 3 -c2 4 -p2 5 --zero-based --temp-dir . ../mm10.chrom.sizes:1000 Ren_mESC_WT_HiC_rep1.nodups.pairs.Qge30.gz Ren_mESC_WT_HiC_rep1.nodups.pairs.1k.cool >Ren_mESC_WT_HiC_rep1.nodups.pairs.cooler_cload.runlog 2>&1 &
nohup time cooler cload pairs --assembly mm10 -c1 2 -p1 3 -c2 4 -p2 5 --zero-based --temp-dir . ../mm10.chrom.sizes:1000 Ren_mESC_WT_HiC_rep2.nodups.pairs.Qge30.gz Ren_mESC_WT_HiC_rep2.nodups.pairs.1k.cool >Ren_mESC_WT_HiC_rep2.nodups.pairs.cooler_cload.runlog 2>&1 &
tail -f Ren_mESC_WT_HiC_rep1.nodups.pairs.cooler_cload.runlog
tail -f Ren_mESC_WT_HiC_rep2.nodups.pairs.cooler_cload.runlog

# rep1/2 raw 1k cool merge
nohup time cooler merge Ren_mESC_WT_HiC.rep12.merged.1k.cool Ren_mESC_WT_HiC_rep1.nodups.pairs.1k.cool Ren_mESC_WT_HiC_rep2.nodups.pairs.1k.cool >Ren_mESC_WT_HiC.rep12.1k.cool.merge.runlog 2>&1 &
tail -f Ren_mESC_WT_HiC.rep12.1k.cool.merge.runlog

cooler info Ren_mESC_WT_HiC_rep1.nodups.pairs.1k.cool
# {
#     "bin-size": 1000,
#     "bin-type": "fixed",
#     "creation-date": "2024-06-12T20:25:44.059581",
#     "format": "HDF5::Cooler",
#     "format-url": "https://github.com/open2c/cooler",
#     "format-version": 3,
#     "generated-by": "cooler-0.10.0",
#     "genome-assembly": "mm10",
#     "nbins": 2725532,
#     "nchroms": 21,
#     "nnz": 235690949,
#     "storage-mode": "symmetric-upper",
#     "sum": 411507394
# }

cooler info Ren_mESC_WT_HiC_rep2.nodups.pairs.1k.cool
# {
#     "bin-size": 1000,
#     "bin-type": "fixed",
#     "creation-date": "2024-06-12T20:34:19.502419",
#     "format": "HDF5::Cooler",
#     "format-url": "https://github.com/open2c/cooler",
#     "format-version": 3,
#     "generated-by": "cooler-0.10.0",
#     "genome-assembly": "mm10",
#     "nbins": 2725532,
#     "nchroms": 21,
#     "nnz": 432716116,
#     "storage-mode": "symmetric-upper",
#     "sum": 605856320
# }

cooler info Ren_mESC_WT_HiC.rep12.merged.1k.cool
# {
#     "bin-size": 1000,
#     "bin-type": "fixed",
#     "creation-date": "2024-06-12T20:47:05.913921",
#     "format": "HDF5::Cooler",
#     "format-url": "https://github.com/open2c/cooler",
#     "format-version": 3,
#     "generated-by": "cooler-0.10.0",
#     "genome-assembly": "mm10",
#     "nbins": 2725532,
#     "nchroms": 21,
#     "nnz": 622246225,
#     "storage-mode": "symmetric-upper",
#     "sum": 1017363714
# }

```

### cooler zoomify cool to mcool

```bash
# zoomify 1k cool to mcool
# -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000
nohup time cooler zoomify -p 16 -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000 --balance -o Ren_mESC_WT_HiC_rep1.nodups.pairs.1k.zoomify.mcool Ren_mESC_WT_HiC_rep1.nodups.pairs.1k.cool >Ren_mESC_WT_HiC_rep1.nodups.pairs.1k.cool.zoomify.runlog 2>&1 &
nohup time cooler zoomify -p 16 -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000 --balance -o Ren_mESC_WT_HiC_rep2.nodups.pairs.1k.zoomify.mcool Ren_mESC_WT_HiC_rep2.nodups.pairs.1k.cool >Ren_mESC_WT_HiC_rep2.nodups.pairs.1k.cool.zoomify.runlog 2>&1 &
tail -f Ren_mESC_WT_HiC_rep1.nodups.pairs.1k.cool.zoomify.runlog
tail -f Ren_mESC_WT_HiC_rep2.nodups.pairs.1k.cool.zoomify.runlog

# zoomify 1k cool to mcool for rep12 merged 1k cool
# -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000
nohup time cooler zoomify -p 16 -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000 --balance -o Ren_mESC_WT_HiC.rep12.merged.1k.zoomify.mcool Ren_mESC_WT_HiC.rep12.merged.1k.cool >Ren_mESC_WT_HiC.rep12.merged.1k.cool.zoomify.runlog 2>&1 &
tail -f Ren_mESC_WT_HiC.rep12.merged.1k.cool.zoomify.runlog

# for i in `ls -1 *.mcool`;do echo $i; cooler info $i"::/resolutions/1000";done

```



## Ren_IAA

```bash
#
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_IAA
nohup time bash ../3_sam_to_sorted_bam_and_index.sh >3_sam_to_sorted_bam_and_index.sh.runlog 2>&1 &
tail -f 3_sam_to_sorted_bam_and_index.sh.runlog

```

### pairtools parse and sort
```bash
# pairtools parse hic pairs
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_IAA
nohup time pairtools --verbose parse -c /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes -o Ren_mESC_IAA_HiC_rep1.pairs.gz --assembly mm10 --min-mapq 20 --drop-seq --drop-sam --add-pair-index --add-columns mapq --output-stats Ren_mESC_IAA_HiC_rep1.parse.stats.txt --walks-policy 5unique --no-flip Ren_mESC_IAA_HiC_rep1.sam >Ren_mESC_IAA_HiC_rep1.pairtools.parse_pairs.runlog 2>&1 &
nohup time pairtools --verbose parse -c /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes -o Ren_mESC_IAA_HiC_rep2.pairs.gz --assembly mm10 --min-mapq 20 --drop-seq --drop-sam --add-pair-index --add-columns mapq --output-stats Ren_mESC_IAA_HiC_rep2.parse.stats.txt --walks-policy 5unique --no-flip Ren_mESC_IAA_HiC_rep2.sam >Ren_mESC_IAA_HiC_rep2.pairtools.parse_pairs.runlog 2>&1 &
tail -f Ren_mESC_IAA_HiC_rep1.pairtools.parse_pairs.runlog
tail -f Ren_mESC_IAA_HiC_rep2.pairtools.parse_pairs.runlog

# Sort pairs using pairtools sort:
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_IAA
nohup time pairtools --verbose sort --nproc 8 --tmpdir ./ -o Ren_mESC_IAA_HiC_rep1.pairs.sorted.gz Ren_mESC_IAA_HiC_rep1.pairs.gz >Ren_mESC_IAA_HiC_rep1.pairtools.sort.runlog 2>&1 &
nohup time pairtools --verbose sort --nproc 8 --tmpdir ./ -o Ren_mESC_IAA_HiC_rep2.pairs.sorted.gz Ren_mESC_IAA_HiC_rep2.pairs.gz >Ren_mESC_IAA_HiC_rep2.pairtools.sort.runlog 2>&1 &
tail -f Ren_mESC_IAA_HiC_rep1.pairtools.sort.runlog
tail -f Ren_mESC_IAA_HiC_rep2.pairtools.sort.runlog

# raw pairs stats
nohup time pairtools stats Ren_mESC_IAA_HiC_rep1.pairs.gz -o Ren_mESC_IAA_HiC_rep1.pairs.stats.tsv >Ren_mESC_IAA_HiC_rep1.pairs.stats.runlog 2>&1 &
nohup time pairtools stats Ren_mESC_IAA_HiC_rep2.pairs.gz -o Ren_mESC_IAA_HiC_rep2.pairs.stats.tsv >Ren_mESC_IAA_HiC_rep2.pairs.stats.runlog 2>&1 &
tail -f Ren_mESC_IAA_HiC_rep1.pairs.stats.runlog
tail -f Ren_mESC_IAA_HiC_rep2.pairs.stats.runlog

pairtools stats -o Ren_mESC_IAA_HiC.rep12.pairs.stats.merge.tsv --merge Ren_mESC_IAA_HiC_rep1.pairs.stats.tsv Ren_mESC_IAA_HiC_rep2.pairs.stats.tsv
wc -l *.tsv
# 2497 Ren_mESC_IAA_HiC.rep12.pairs.stats.merge.tsv
# 2397 Ren_mESC_IAA_HiC_rep1.pairs.stats.tsv
# 2427 Ren_mESC_IAA_HiC_rep2.pairs.stats.tsv

head *.tsv
# ==> Ren_mESC_IAA_HiC.rep12.pairs.stats.merge.tsv <==
# total	1347065014
# total_unmapped	160695884
# total_single_sided_mapped	399285076
# total_mapped	787084054
# total_dups	0
# total_nodups	787084054
# cis	678240735
# trans	108843319
# pair_types/MM	123530697
# pair_types/RM	80378

# ==> Ren_mESC_IAA_HiC_rep1.pairs.stats.tsv <==
# total	630346337
# total_unmapped	87552374
# total_single_sided_mapped	176417541
# total_mapped	366376422
# total_dups	0
# total_nodups	366376422
# cis	314045439
# trans	52330983
# pair_types/UU	366350187
# pair_types/MM	70408531

# ==> Ren_mESC_IAA_HiC_rep2.pairs.stats.tsv <==
# total	716718677
# total_unmapped	73143510
# total_single_sided_mapped	222867535
# total_mapped	420707632
# total_dups	0
# total_nodups	420707632
# cis	364195296
# trans	56512336
# pair_types/UU	420702490
# pair_types/UM	76464046

```

### pairtools dedup
```bash
## Detect and remove duplicates using pairtools dedup and generate statistics
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_IAA
nohup time pairtools --verbose dedup -p 16 \
 --output Ren_mESC_IAA_HiC_rep1.nodups.pairs.gz \
 --output-dups Ren_mESC_IAA_HiC_rep1.dups.pairs.gz \
 --output-unmapped Ren_mESC_IAA_HiC_rep1.unmapped.pairs.gz \
 --output-stats Ren_mESC_IAA_HiC_rep1.dedup.stats.txt \
 Ren_mESC_IAA_HiC_rep1.pairs.sorted.gz \
 >Ren_mESC_IAA_HiC_rep1.pairtools.dedup.runlog 2>&1 &
nohup time pairtools --verbose dedup -p 16 \
 --output Ren_mESC_IAA_HiC_rep2.nodups.pairs.gz \
 --output-dups Ren_mESC_IAA_HiC_rep2.dups.pairs.gz \
 --output-unmapped Ren_mESC_IAA_HiC_rep2.unmapped.pairs.gz \
 --output-stats Ren_mESC_IAA_HiC_rep2.dedup.stats.txt \
 Ren_mESC_IAA_HiC_rep2.pairs.sorted.gz \
 >Ren_mESC_IAA_HiC_rep2.pairtools.dedup.runlog 2>&1 &
tail -f Ren_mESC_IAA_HiC_rep1.pairtools.dedup.runlog
tail -f Ren_mESC_IAA_HiC_rep2.pairtools.dedup.runlog

```

```bash
# pairs type stat
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_IAA
# zcat Ren_mESC_IAA_HiC_rep1.pairs.gz |grep -v "#" |cut -f 8 |head -n 100 |sort |uniq -c |sort -k1,1nr
zcat Ren_mESC_IAA_HiC_rep1.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 366350187 UU
#  70408531 MM
#  62998641 UM
#  61274376 MU
#  32714867 UN
#  17758987 NU
#   9505766 MN
#   5096996 NM
#   2541081 NN
#    908168 RN
#    696296 NR
#     45418 RM
#     24506 RU
#     20788 MR
#      1729 UR
zcat Ren_mESC_IAA_HiC_rep1.pairs.sorted.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 366350187 UU
#  70408531 MM
#  62998641 UM
#  61274376 MU
#  32714867 UN
#  17758987 NU
#   9505766 MN
#   5096996 NM
#   2541081 NN
#    908168 RN
#    696296 NR
#     45418 RM
#     24506 RU
#     20788 MR
#      1729 UR
zcat Ren_mESC_IAA_HiC_rep1.dups.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 86869062 DD
zcat Ren_mESC_IAA_HiC_rep1.nodups.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 279486339 UU
#     19629 RU
#      1392 UR


zcat Ren_mESC_IAA_HiC_rep2.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 420702490 UU
#  76464046 UM
#  74257727 MU
#  53122166 MM
#  44440787 UN
#  25489712 NU
#  10810615 MN
#   6254492 NM
#   2956237 NN
#   1095886 RN
#   1059596 NR
#     34960 RM
#     24821 MR
#      4426 RU
#       716 UR
zcat Ren_mESC_IAA_HiC_rep2.pairs.sorted.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 420702490 UU
#  76464046 UM
#  74257727 MU
#  53122166 MM
#  44440787 UN
#  25489712 NU
#  10810615 MN
#   6254492 NM
#   2956237 NN
#   1095886 RN
#   1059596 NR
#     34960 RM
#     24821 MR
#      4426 RU
#       716 UR
zcat Ren_mESC_IAA_HiC_rep2.dups.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 30568476 DD
zcat Ren_mESC_IAA_HiC_rep2.nodups.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 390134282 UU
#      4194 RU
#       680 UR

```

### cooler cload pairs (pairs to cool)

```bash
# ref https://pairtools.readthedocs.io/en/latest/protocols_pipelines.html
# Recommended pairtools parameters for standard Hi-C protocols
# https://pairtools.readthedocs.io/en/latest/examples/pairtools_walkthrough.html
# select MAPQ1/2 >= 30 
mamba activate hicenv
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_IAA
nohup time pairtools --verbose select "(mapq1>=30) and (mapq2>=30)" Ren_mESC_IAA_HiC_rep1.nodups.pairs.gz -o Ren_mESC_IAA_HiC_rep1.nodups.pairs.Qge30.gz >Ren_mESC_IAA_HiC_rep1.nodups.pairs.select_Qge30.runlog 2>&1 &
nohup time pairtools --verbose select "(mapq1>=30) and (mapq2>=30)" Ren_mESC_IAA_HiC_rep2.nodups.pairs.gz -o Ren_mESC_IAA_HiC_rep2.nodups.pairs.Qge30.gz >Ren_mESC_IAA_HiC_rep2.nodups.pairs.select_Qge30.runlog 2>&1 &
tail -f Ren_mESC_IAA_HiC_rep1.nodups.pairs.select_Qge30.runlog
tail -f Ren_mESC_IAA_HiC_rep2.nodups.pairs.select_Qge30.runlog

# pairs to 1k cool
nohup time cooler cload pairs --assembly mm10 -c1 2 -p1 3 -c2 4 -p2 5 --zero-based --temp-dir . ../mm10.chrom.sizes:1000 Ren_mESC_IAA_HiC_rep1.nodups.pairs.Qge30.gz Ren_mESC_IAA_HiC_rep1.nodups.pairs.1k.cool >Ren_mESC_IAA_HiC_rep1.nodups.pairs.cooler_cload.runlog 2>&1 &
nohup time cooler cload pairs --assembly mm10 -c1 2 -p1 3 -c2 4 -p2 5 --zero-based --temp-dir . ../mm10.chrom.sizes:1000 Ren_mESC_IAA_HiC_rep2.nodups.pairs.Qge30.gz Ren_mESC_IAA_HiC_rep2.nodups.pairs.1k.cool >Ren_mESC_IAA_HiC_rep2.nodups.pairs.cooler_cload.runlog 2>&1 &
tail -f Ren_mESC_IAA_HiC_rep1.nodups.pairs.cooler_cload.runlog
tail -f Ren_mESC_IAA_HiC_rep2.nodups.pairs.cooler_cload.runlog

# rep1/2 raw 1k cool merge
nohup time cooler merge Ren_mESC_IAA_HiC.rep12.merged.1k.cool Ren_mESC_IAA_HiC_rep1.nodups.pairs.1k.cool Ren_mESC_IAA_HiC_rep2.nodups.pairs.1k.cool >Ren_mESC_IAA_HiC.rep12.1k.cool.merge.runlog 2>&1 &
tail -f Ren_mESC_IAA_HiC.rep12.1k.cool.merge.runlog

cooler info Ren_mESC_IAA_HiC_rep1.nodups.pairs.1k.cool
# {
#     "bin-size": 1000,
#     "bin-type": "fixed",
#     "creation-date": "2024-06-12T20:21:42.609368",
#     "format": "HDF5::Cooler",
#     "format-url": "https://github.com/open2c/cooler",
#     "format-version": 3,
#     "generated-by": "cooler-0.10.0",
#     "genome-assembly": "mm10",
#     "nbins": 2725532,
#     "nchroms": 21,
#     "nnz": 163226995,
#     "storage-mode": "symmetric-upper",
#     "sum": 270698316
# }

cooler info Ren_mESC_IAA_HiC_rep2.nodups.pairs.1k.cool
# {
#     "bin-size": 1000,
#     "bin-type": "fixed",
#     "creation-date": "2024-06-12T20:26:35.902143",
#     "format": "HDF5::Cooler",
#     "format-url": "https://github.com/open2c/cooler",
#     "format-version": 3,
#     "generated-by": "cooler-0.10.0",
#     "genome-assembly": "mm10",
#     "nbins": 2725532,
#     "nchroms": 21,
#     "nnz": 281952978,
#     "storage-mode": "symmetric-upper",
#     "sum": 377961858
# }

cooler info Ren_mESC_IAA_HiC.rep12.merged.1k.cool
# {
#     "bin-size": 1000,
#     "bin-type": "fixed",
#     "creation-date": "2024-06-12T20:39:21.355340",
#     "format": "HDF5::Cooler",
#     "format-url": "https://github.com/open2c/cooler",
#     "format-version": 3,
#     "generated-by": "cooler-0.10.0",
#     "genome-assembly": "mm10",
#     "nbins": 2725532,
#     "nchroms": 21,
#     "nnz": 418264161,
#     "storage-mode": "symmetric-upper",
#     "sum": 648660174
# }

```

### cooler zoomify cool to mcool

```bash
# zoomify 1k cool to mcool
# -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000
nohup time cooler zoomify -p 16 -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000 --balance -o Ren_mESC_IAA_HiC_rep1.nodups.pairs.1k.zoomify.mcool Ren_mESC_IAA_HiC_rep1.nodups.pairs.1k.cool >Ren_mESC_IAA_HiC_rep1.nodups.pairs.1k.cool.zoomify.runlog 2>&1 &
nohup time cooler zoomify -p 16 -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000 --balance -o Ren_mESC_IAA_HiC_rep2.nodups.pairs.1k.zoomify.mcool Ren_mESC_IAA_HiC_rep2.nodups.pairs.1k.cool >Ren_mESC_IAA_HiC_rep2.nodups.pairs.1k.cool.zoomify.runlog 2>&1 &
tail -f Ren_mESC_IAA_HiC_rep1.nodups.pairs.1k.cool.zoomify.runlog
tail -f Ren_mESC_IAA_HiC_rep2.nodups.pairs.1k.cool.zoomify.runlog

# zoomify 1k cool to mcool for rep12 merged 1k cool
# -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000
nohup time cooler zoomify -p 16 -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000 --balance -o Ren_mESC_IAA_HiC.rep12.merged.1k.zoomify.mcool Ren_mESC_IAA_HiC.rep12.merged.1k.cool >Ren_mESC_IAA_HiC.rep12.merged.1k.cool.zoomify.runlog 2>&1 &
tail -f Ren_mESC_IAA_HiC.rep12.merged.1k.cool.zoomify.runlog

# for i in `ls -1 *.mcool`;do echo $i; cooler info $i"::/resolutions/1000";done

```



## WN_BPK25

```bash
# 
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/WN_BPK25
nohup time bash ../3_sam_to_sorted_bam_and_index.sh >3_sam_to_sorted_bam_and_index.sh.runlog 2>&1 &
tail -f 3_sam_to_sorted_bam_and_index.sh.runlog

```

### pairtools parse and sort
```bash
# pairtools parse hic pairs
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/WN_BPK25
nohup time pairtools --verbose parse -c /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes -o WN_mESC_BPK25_HiC_rep1.pairs.gz --assembly mm10 --min-mapq 20 --drop-seq --drop-sam --add-pair-index --add-columns mapq --output-stats WN_mESC_BPK25_HiC_rep1.parse.stats.txt --walks-policy 5unique --no-flip WN_mESC_BPK25_HiC_rep1.sam >WN_mESC_BPK25_HiC_rep1.pairtools.parse_pairs.runlog 2>&1 &
nohup time pairtools --verbose parse -c /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes -o WN_mESC_BPK25_HiC_rep2.pairs.gz --assembly mm10 --min-mapq 20 --drop-seq --drop-sam --add-pair-index --add-columns mapq --output-stats WN_mESC_BPK25_HiC_rep2.parse.stats.txt --walks-policy 5unique --no-flip WN_mESC_BPK25_HiC_rep2.sam >WN_mESC_BPK25_HiC_rep2.pairtools.parse_pairs.runlog 2>&1 &
tail -f WN_mESC_BPK25_HiC_rep1.pairtools.parse_pairs.runlog
tail -f WN_mESC_BPK25_HiC_rep2.pairtools.parse_pairs.runlog

# Sort pairs using pairtools sort:
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/WN_BPK25
nohup time pairtools --verbose sort --nproc 8 --tmpdir ./ -o WN_mESC_BPK25_HiC_rep1.pairs.sorted.gz WN_mESC_BPK25_HiC_rep1.pairs.gz >WN_mESC_BPK25_HiC_rep1.pairtools.sort.runlog 2>&1 &
nohup time pairtools --verbose sort --nproc 8 --tmpdir ./ -o WN_mESC_BPK25_HiC_rep2.pairs.sorted.gz WN_mESC_BPK25_HiC_rep2.pairs.gz >WN_mESC_BPK25_HiC_rep2.pairtools.sort.runlog 2>&1 &
tail -f WN_mESC_BPK25_HiC_rep1.pairtools.sort.runlog
tail -f WN_mESC_BPK25_HiC_rep2.pairtools.sort.runlog

# raw pairs stats
nohup time pairtools stats WN_mESC_BPK25_HiC_rep1.pairs.gz -o WN_mESC_BPK25_HiC_rep1.pairs.stats.tsv >WN_mESC_BPK25_HiC_rep1.pairs.stats.runlog 2>&1 &
nohup time pairtools stats WN_mESC_BPK25_HiC_rep2.pairs.gz -o WN_mESC_BPK25_HiC_rep2.pairs.stats.tsv >WN_mESC_BPK25_HiC_rep2.pairs.stats.runlog 2>&1 &
tail -f WN_mESC_BPK25_HiC_rep1.pairs.stats.runlog
tail -f WN_mESC_BPK25_HiC_rep2.pairs.stats.runlog

pairtools stats -o WN_mESC_BPK25_HiC.rep12.pairs.stats.merge.tsv --merge WN_mESC_BPK25_HiC_rep1.pairs.stats.tsv WN_mESC_BPK25_HiC_rep2.pairs.stats.tsv
wc -l *.tsv
# 2559 WN_mESC_BPK25_HiC.rep12.pairs.stats.merge.tsv
# 2514 WN_mESC_BPK25_HiC_rep1.pairs.stats.tsv
# 2405 WN_mESC_BPK25_HiC_rep2.pairs.stats.tsv

head *.tsv
# ==> WN_mESC_BPK25_HiC.rep12.pairs.stats.merge.tsv <==
# total	975875572
# total_unmapped	101068539
# total_single_sided_mapped	218881698
# total_mapped	655925335
# total_dups	0
# total_nodups	655925335
# cis	476036353
# trans	179888982
# pair_types/MN	6041437
# pair_types/MM	76634146

# ==> WN_mESC_BPK25_HiC_rep1.pairs.stats.tsv <==
# total	597286793
# total_unmapped	56364641
# total_single_sided_mapped	139050690
# total_mapped	401871462
# total_dups	0
# total_nodups	401871462
# cis	282628720
# trans	119242742
# pair_types/UU	229725157
# pair_types/RU	85960368

# ==> WN_mESC_BPK25_HiC_rep2.pairs.stats.tsv <==
# total	378588779
# total_unmapped	44703898
# total_single_sided_mapped	79831008
# total_mapped	254053873
# total_dups	0
# total_nodups	254053873
# cis	193407633
# trans	60646240
# pair_types/UU	158055216
# pair_types/RU	48112073

```

### pairtools dedup
```bash
## Detect and remove duplicates using pairtools dedup and generate statistics
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/WN_BPK25
nohup time pairtools --verbose dedup -p 16 \
 --output WN_mESC_BPK25_HiC_rep1.nodups.pairs.gz \
 --output-dups WN_mESC_BPK25_HiC_rep1.dups.pairs.gz \
 --output-unmapped WN_mESC_BPK25_HiC_rep1.unmapped.pairs.gz \
 --output-stats WN_mESC_BPK25_HiC_rep1.dedup.stats.txt \
 WN_mESC_BPK25_HiC_rep1.pairs.sorted.gz \
 >WN_mESC_BPK25_HiC_rep1.pairtools.dedup.runlog 2>&1 &
nohup time pairtools --verbose dedup -p 16 \
 --output WN_mESC_BPK25_HiC_rep2.nodups.pairs.gz \
 --output-dups WN_mESC_BPK25_HiC_rep2.dups.pairs.gz \
 --output-unmapped WN_mESC_BPK25_HiC_rep2.unmapped.pairs.gz \
 --output-stats WN_mESC_BPK25_HiC_rep2.dedup.stats.txt \
 WN_mESC_BPK25_HiC_rep2.pairs.sorted.gz \
 >WN_mESC_BPK25_HiC_rep2.pairtools.dedup.runlog 2>&1 &
tail -f WN_mESC_BPK25_HiC_rep1.pairtools.dedup.runlog
tail -f WN_mESC_BPK25_HiC_rep2.pairtools.dedup.runlog

```

```bash
# pairs type stat
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/WN_BPK25
# zcat WN_mESC_BPK25_HiC_rep1.pairs.gz |grep -v "#" |cut -f 8 |head -n 100 |sort |uniq -c |sort -k1,1nr
zcat WN_mESC_BPK25_HiC_rep1.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 229725157 UU
#  86185937 UR
#  85960368 RU
#  44841588 MM
#  25641411 UM
#  21973985 RM
#  21886152 MU
#  21551974 MR
#  20419354 RN
#  18549271 NR
#   7647264 UN
#   4975557 NM
#   3802035 MN
#   2745461 NN
#   1381279 NU
zcat WN_mESC_BPK25_HiC_rep1.pairs.sorted.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 229725157 UU
#  86185937 UR
#  85960368 RU
#  44841588 MM
#  25641411 UM
#  21973985 RM
#  21886152 MU
#  21551974 MR
#  20419354 RN
#  18549271 NR
#   7647264 UN
#   4975557 NM
#   3802035 MN
#   2745461 NN
#   1381279 NU
zcat WN_mESC_BPK25_HiC_rep1.dups.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 135712450 DD
zcat WN_mESC_BPK25_HiC_rep1.nodups.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 151386313 UU
#  57393476 RU
#  57379223 UR


zcat WN_mESC_BPK25_HiC_rep2.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 158055216 UU
#  48112073 RU
#  47886584 UR
#  31792558 MM
#  16550353 UM
#  13193753 MU
#  12276892 RM
#  11753005 MR
#  10827593 RN
#   9479638 NR
#   7642349 NM
#   4870556 UN
#   3029589 NN
#   2239402 MN
#    879218 NU
zcat WN_mESC_BPK25_HiC_rep2.pairs.sorted.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 158055216 UU
#  48112073 RU
#  47886584 UR
#  31792558 MM
#  16550353 UM
#  13193753 MU
#  12276892 RM
#  11753005 MR
#  10827593 RN
#   9479638 NR
#   7642349 NM
#   4870556 UN
#   3029589 NN
#   2239402 MN
#    879218 NU
zcat WN_mESC_BPK25_HiC_rep2.dups.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 99460894 DD
zcat WN_mESC_BPK25_HiC_rep2.nodups.pairs.gz |grep -v "#" |cut -f 8 |sort |uniq -c |sort -k1,1nr
# 95533714 UU
# 29640911 RU
# 29418354 UR

```

### cooler cload pairs (pairs to cool)

```bash
# ref https://pairtools.readthedocs.io/en/latest/protocols_pipelines.html
# Recommended pairtools parameters for standard Hi-C protocols
# https://pairtools.readthedocs.io/en/latest/examples/pairtools_walkthrough.html
# select MAPQ1/2 >= 30 
mamba activate hicenv
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/WN_BPK25
nohup time pairtools --verbose select "(mapq1>=30) and (mapq2>=30)" WN_mESC_BPK25_HiC_rep1.nodups.pairs.gz -o WN_mESC_BPK25_HiC_rep1.nodups.pairs.Qge30.gz >WN_mESC_BPK25_HiC_rep1.nodups.pairs.select_Qge30.runlog 2>&1 &
nohup time pairtools --verbose select "(mapq1>=30) and (mapq2>=30)" WN_mESC_BPK25_HiC_rep2.nodups.pairs.gz -o WN_mESC_BPK25_HiC_rep2.nodups.pairs.Qge30.gz >WN_mESC_BPK25_HiC_rep2.nodups.pairs.select_Qge30.runlog 2>&1 &
tail -f WN_mESC_BPK25_HiC_rep1.nodups.pairs.select_Qge30.runlog
tail -f WN_mESC_BPK25_HiC_rep2.nodups.pairs.select_Qge30.runlog

# pairs to 1k cool
nohup time cooler cload pairs --assembly mm10 -c1 2 -p1 3 -c2 4 -p2 5 --zero-based --temp-dir . ../mm10.chrom.sizes:1000 WN_mESC_BPK25_HiC_rep1.nodups.pairs.Qge30.gz WN_mESC_BPK25_HiC_rep1.nodups.pairs.1k.cool >WN_mESC_BPK25_HiC_rep1.nodups.pairs.cooler_cload.runlog 2>&1 &
nohup time cooler cload pairs --assembly mm10 -c1 2 -p1 3 -c2 4 -p2 5 --zero-based --temp-dir . ../mm10.chrom.sizes:1000 WN_mESC_BPK25_HiC_rep2.nodups.pairs.Qge30.gz WN_mESC_BPK25_HiC_rep2.nodups.pairs.1k.cool >WN_mESC_BPK25_HiC_rep2.nodups.pairs.cooler_cload.runlog 2>&1 &
tail -f WN_mESC_BPK25_HiC_rep1.nodups.pairs.cooler_cload.runlog
tail -f WN_mESC_BPK25_HiC_rep2.nodups.pairs.cooler_cload.runlog

# rep1/2 raw 1k cool merge
nohup time cooler merge WN_mESC_BPK25_HiC.rep12.merged.1k.cool WN_mESC_BPK25_HiC_rep1.nodups.pairs.1k.cool WN_mESC_BPK25_HiC_rep2.nodups.pairs.1k.cool >WN_mESC_BPK25_HiC.rep12.1k.cool.merge.runlog 2>&1 &
tail -f WN_mESC_BPK25_HiC.rep12.1k.cool.merge.runlog

cooler info WN_mESC_BPK25_HiC_rep1.nodups.pairs.1k.cool
# {
#     "bin-size": 1000,
#     "bin-type": "fixed",
#     "creation-date": "2024-06-12T20:23:03.554739",
#     "format": "HDF5::Cooler",
#     "format-url": "https://github.com/open2c/cooler",
#     "format-version": 3,
#     "generated-by": "cooler-0.10.0",
#     "genome-assembly": "mm10",
#     "nbins": 2725532,
#     "nchroms": 21,
#     "nnz": 186536134,
#     "storage-mode": "symmetric-upper",
#     "sum": 257280453
# }

cooler info WN_mESC_BPK25_HiC_rep2.nodups.pairs.1k.cool
# {
#     "bin-size": 1000,
#     "bin-type": "fixed",
#     "creation-date": "2024-06-12T20:17:49.417259",
#     "format": "HDF5::Cooler",
#     "format-url": "https://github.com/open2c/cooler",
#     "format-version": 3,
#     "generated-by": "cooler-0.10.0",
#     "genome-assembly": "mm10",
#     "nbins": 2725532,
#     "nchroms": 21,
#     "nnz": 94892694,
#     "storage-mode": "symmetric-upper",
#     "sum": 149814135
# }

cooler info WN_mESC_BPK25_HiC.rep12.merged.1k.cool
# {
#     "bin-size": 1000,
#     "bin-type": "fixed",
#     "creation-date": "2024-06-12T20:37:53.485310",
#     "format": "HDF5::Cooler",
#     "format-url": "https://github.com/open2c/cooler",
#     "format-version": 3,
#     "generated-by": "cooler-0.10.0",
#     "genome-assembly": "mm10",
#     "nbins": 2725532,
#     "nchroms": 21,
#     "nnz": 271853397,
#     "storage-mode": "symmetric-upper",
#     "sum": 407094588
# }

```

### cooler zoomify cool to mcool

```bash
# zoomify 1k cool to mcool
# -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000
nohup time cooler zoomify -p 16 -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000 --balance -o WN_mESC_BPK25_HiC_rep1.nodups.pairs.1k.zoomify.mcool WN_mESC_BPK25_HiC_rep1.nodups.pairs.1k.cool >WN_mESC_BPK25_HiC_rep1.nodups.pairs.1k.cool.zoomify.runlog 2>&1 &
nohup time cooler zoomify -p 16 -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000 --balance -o WN_mESC_BPK25_HiC_rep2.nodups.pairs.1k.zoomify.mcool WN_mESC_BPK25_HiC_rep2.nodups.pairs.1k.cool >WN_mESC_BPK25_HiC_rep2.nodups.pairs.1k.cool.zoomify.runlog 2>&1 &
tail -f WN_mESC_BPK25_HiC_rep1.nodups.pairs.1k.cool.zoomify.runlog
tail -f WN_mESC_BPK25_HiC_rep2.nodups.pairs.1k.cool.zoomify.runlog

# zoomify 1k cool to mcool for rep12 merged 1k cool
# -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000
nohup time cooler zoomify -p 16 -r 1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000,75000,100000,150000,200000,250000,500000 --balance -o WN_mESC_BPK25_HiC.rep12.merged.1k.zoomify.mcool WN_mESC_BPK25_HiC.rep12.merged.1k.cool >WN_mESC_BPK25_HiC.rep12.merged.1k.cool.zoomify.runlog 2>&1 &
tail -f WN_mESC_BPK25_HiC.rep12.merged.1k.cool.zoomify.runlog

# for i in `ls -1 *.mcool`;do echo $i; cooler info $i"::/resolutions/1000";done

```

# replot Nora paper heatmap Figs

> ref  https://cooltools.readthedocs.io/en/latest/notebooks/viz.html
>

## Nora_WT

```bash
# 
# 4_cooltools_plot_cmd_use_Nora_WT_mcool.r20k.ipynb
# 4_cooltools_plot_cmd_use_Nora_WT_mcool.r10k.ipynb
# 4_cooltools_plot_cmd_use_Nora_WT_mcool.r5k.ipynb
# 4_cooltools_plot_cmd_use_Nora_WT_mcool.r1k.ipynb

```

## Nora_IAA

```bash
# 
# 4_cooltools_plot_cmd_use_Nora_IAA_mcool.r20k.ipynb
# 4_cooltools_plot_cmd_use_Nora_IAA_mcool.r10k.ipynb
# 4_cooltools_plot_cmd_use_Nora_IAA_mcool.r5k.ipynb

```

## Ren_WT

```bash
# 
# 4_cooltools_plot_cmd_use_Ren_WT_mcool.r20k.ipynb
# 4_cooltools_plot_cmd_use_Ren_WT_mcool.r10k.ipynb
# 4_cooltools_plot_cmd_use_Ren_WT_mcool.r5k.ipynb

```

## Ren_IAA

```bash
# 
# 4_cooltools_plot_cmd_use_Ren_IAA_mcool.r20k.ipynb
# 4_cooltools_plot_cmd_use_Ren_IAA_mcool.r10k.ipynb
# 4_cooltools_plot_cmd_use_Ren_IAA_mcool.r5k.ipynb

```

## WN_BPK25

```bash
# 
# 4_cooltools_plot_cmd_use_WN_BPK25_mcool.r20k.ipynb
# 4_cooltools_plot_cmd_use_WN_BPK25_mcool.r10k.ipynb
# 4_cooltools_plot_cmd_use_WN_BPK25_mcool.r5k.ipynb

```

## Nora_WT_IAA_WN_BPK25

```bash
# 
# 4_cooltools_plot_cmd_use_Nora_WN_mcool.r20k.ipynb
# 4_cooltools_plot_cmd_use_Nora_WN_mcool.r10k.ipynb

```

## Ren_WT_IAA_WN_BPK25

```bash
# 
# 4_cooltools_plot_cmd_use_Ren_WN_mcool.r20k.ipynb
# 4_cooltools_plot_cmd_use_Ren_WN_mcool.r10k.ipynb

```

## Nora_Ren_WN_all

```bash
# 
# 4_Nora_Ren_WN.all.rep12merged.fall.heatmap.summary.ipynb

# 
# 5_cooltools_plot_cmd_for_chromsome_view.ipynb

```



#  AB compartment anlysis and saddleplots

## Nora_WT

```bash
# test
# 5_cooltools_cmd_for_expect_AB_TAD_loops.under_diff_resolution.Nora_WT_mcool.test.ipynb

# 
# /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_WT
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r200k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r100k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r50k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r20k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r10k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r5k.ipynb

```

## Nora_IAA

```bash
# 
# /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_IAA
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r200k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r100k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r50k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r20k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r10k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r5k.ipynb

```

## Ren_WT

```bash
# 
# /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_WT
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r200k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r100k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r50k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r20k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r10k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r5k.ipynb

```

## Ren_IAA

```bash
# 
# /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_IAA
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r200k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r100k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r50k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r20k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r10k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r5k.ipynb

```

## WN_BPK25

```bash
# 
# /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/WN_BPK25
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r200k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r100k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r50k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r20k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r10k.ipynb
# 5_cooltools_cmd_for_saddleplot_of_expect_and_AB_compartment.under_r5k.ipynb

```



# insulation boundaries and TADs analysis and plot

## Nora_WT

```bash
# 
# /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_WT
# 6_cooltools_cmd_for_insulation_boundaries_and_TADs.under_r50k.ipynb
# 6_cooltools_cmd_for_insulation_boundaries_and_TADs.under_r20k.ipynb
# 6_cooltools_cmd_for_insulation_boundaries_and_TADs.under_r10k.ipynb
# 6_cooltools_cmd_for_insulation_boundaries_and_TADs.under_r5k.ipynb

```

## Nora_IAA

```bash
# 
# /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_IAA
# 6_cooltools_cmd_for_insulation_boundaries_and_TADs.under_r50k.ipynb
# 6_cooltools_cmd_for_insulation_boundaries_and_TADs.under_r20k.ipynb
# 6_cooltools_cmd_for_insulation_boundaries_and_TADs.under_r10k.ipynb
# 6_cooltools_cmd_for_insulation_boundaries_and_TADs.under_r5k.ipynb

```

## Ren_WT

```bash
# 
# /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_WT
# 6_cooltools_cmd_for_insulation_boundaries_and_TADs.under_r50k.ipynb
# 6_cooltools_cmd_for_insulation_boundaries_and_TADs.under_r20k.ipynb
# 6_cooltools_cmd_for_insulation_boundaries_and_TADs.under_r10k.ipynb
# 6_cooltools_cmd_for_insulation_boundaries_and_TADs.under_r5k.ipynb

```

## Ren_IAA

```bash
# 
# /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_IAA
# 6_cooltools_cmd_for_insulation_boundaries_and_TADs.under_r50k.ipynb
# 6_cooltools_cmd_for_insulation_boundaries_and_TADs.under_r20k.ipynb
# 6_cooltools_cmd_for_insulation_boundaries_and_TADs.under_r10k.ipynb
# 6_cooltools_cmd_for_insulation_boundaries_and_TADs.under_r5k.ipynb

```

## WN_BPK25

```bash
# 
# /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/WN_BPK25
# 6_cooltools_cmd_for_insulation_boundaries_and_TADs.under_r50k.ipynb
# 6_cooltools_cmd_for_insulation_boundaries_and_TADs.under_r20k.ipynb
# 6_cooltools_cmd_for_insulation_boundaries_and_TADs.under_r10k.ipynb
# 6_cooltools_cmd_for_insulation_boundaries_and_TADs.under_r5k.ipynb

```

## Nora_WT_IAA_WN_BPK25

```bash
# 
# 6_cooltools_plot_cmd_use_Nora_WN_mcool.r20k.with_TAD_and_loops.ipynb
# 6_cooltools_plot_cmd_use_Nora_WN_mcool.r10k.with_TAD_and_loops.ipynb

```

## Ren_WT_IAA_WN_BPK25

```bash
# 
# 6_cooltools_plot_cmd_use_Ren_WN_mcool.r20k.with_TAD_and_loops.ipynb
# 6_cooltools_plot_cmd_use_Ren_WN_mcool.r10k.with_TAD_and_loops.ipynb

```

## Nora_Ren_WN_all

```bash
# 
# 6_Nora_Ren_WN.all.rep12merged.fall.heatmap.with_tads_and_loops.summary.ipynb

```



# dots / loops analysis and plot


## Nora_WT

```bash
# 
# /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_WT
# 7_cooltools_cmd_for_dots_and_loops_calling.under_r20k.ipynb
# 7_cooltools_cmd_for_dots_and_loops_calling.under_r10k.ipynb
# 7_cooltools_cmd_for_dots_and_loops_calling.under_r5k.ipynb

```

## Nora_IAA

```bash
# 
# /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_IAA
# 7_cooltools_cmd_for_dots_and_loops_calling.under_r20k.ipynb
# 7_cooltools_cmd_for_dots_and_loops_calling.under_r10k.ipynb
# 7_cooltools_cmd_for_dots_and_loops_calling.under_r5k.ipynb

```

## Ren_WT

```bash
# 
# /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_WT
# 7_cooltools_cmd_for_dots_and_loops_calling.under_r20k.ipynb
# 7_cooltools_cmd_for_dots_and_loops_calling.under_r10k.ipynb
# 7_cooltools_cmd_for_dots_and_loops_calling.under_r5k.ipynb

```

## Ren_IAA

```bash
# 
# /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_IAA
# 7_cooltools_cmd_for_dots_and_loops_calling.under_r20k.ipynb
# 7_cooltools_cmd_for_dots_and_loops_calling.under_r10k.ipynb
# 7_cooltools_cmd_for_dots_and_loops_calling.under_r5k.ipynb

```

## WN_BPK25

```bash
# 
# /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/WN_BPK25
# 7_cooltools_cmd_for_dots_and_loops_calling.under_r20k.ipynb
# 7_cooltools_cmd_for_dots_and_loops_calling.under_r10k.ipynb
# 7_cooltools_cmd_for_dots_and_loops_calling.under_r5k.ipynb

```

## Nora_WT_IAA_WN_BPK25

```bash
# 
# 7_cooltools_plot_cmd_use_Nora_WN_mcool.r20k.with_loops.ipynb
# 7_cooltools_plot_cmd_use_Nora_WN_mcool.r10k.with_loops.ipynb

```

## Ren_WT_IAA_WN_BPK25

```bash
# 
# 7_cooltools_plot_cmd_use_Ren_WN_mcool.r20k.with_loops.ipynb
# 7_cooltools_plot_cmd_use_Ren_WN_mcool.r10k.with_loops.ipynb

```

## Nora_Ren_WN_all

```bash
# 
# 7_Nora_Ren_WN.all.rep12merged.fall.heatmap.with_loops.summary.ipynb

```



# coolpuppy pileup analysis

## WT d0246 mESC CTCF summit flank pileup

### CTCF summit bed3+FC

```bash
# mESC WT CTCF narrowPeak 
cd /home/sunwenju/projects_on_940xa/2022121401_WN_mESC_d0246_CTCF_MBD3_ChIP_Seq/macs2_peak_calling_2_for_mESC_d0246_WT_CTCF_rep123_together
find . -name *narrowPeak |sort 
# ./mESC_d0_WT_IP_CTCF_rep123_TvsC/mESC_d0_WT_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# ./mESC_d2_WT_IP_CTCF_rep123_TvsC/mESC_d2_WT_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# ./mESC_d4_WT_IP_CTCF_rep123_TvsC/mESC_d4_WT_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# ./mESC_d6_WT_IP_CTCF_rep123_TvsC/mESC_d6_WT_IP_CTCF_rep123_TvsC_peaks.narrowPeak
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline
ln -s /home/sunwenju/projects_on_940xa/2022121401_WN_mESC_d0246_CTCF_MBD3_ChIP_Seq/macs2_peak_calling_2_for_mESC_d0246_WT_CTCF_rep123_together/mESC_d0_WT_IP_CTCF_rep123_TvsC/mESC_d0_WT_IP_CTCF_rep123_TvsC_peaks.narrowPeak mESC_d0_WT_CTCF.narrowPeak
ln -s /home/sunwenju/projects_on_940xa/2022121401_WN_mESC_d0246_CTCF_MBD3_ChIP_Seq/macs2_peak_calling_2_for_mESC_d0246_WT_CTCF_rep123_together/mESC_d2_WT_IP_CTCF_rep123_TvsC/mESC_d2_WT_IP_CTCF_rep123_TvsC_peaks.narrowPeak mESC_d2_WT_CTCF.narrowPeak
ln -s /home/sunwenju/projects_on_940xa/2022121401_WN_mESC_d0246_CTCF_MBD3_ChIP_Seq/macs2_peak_calling_2_for_mESC_d0246_WT_CTCF_rep123_together/mESC_d4_WT_IP_CTCF_rep123_TvsC/mESC_d4_WT_IP_CTCF_rep123_TvsC_peaks.narrowPeak mESC_d4_WT_CTCF.narrowPeak
ln -s /home/sunwenju/projects_on_940xa/2022121401_WN_mESC_d0246_CTCF_MBD3_ChIP_Seq/macs2_peak_calling_2_for_mESC_d0246_WT_CTCF_rep123_together/mESC_d6_WT_IP_CTCF_rep123_TvsC/mESC_d6_WT_IP_CTCF_rep123_TvsC_peaks.narrowPeak mESC_d6_WT_CTCF.narrowPeak
wc -l mESC*CTCF.narrowPeak
# 50544 mESC_d0_WT_CTCF.narrowPeak
# 28732 mESC_d2_WT_CTCF.narrowPeak
# 37518 mESC_d4_WT_CTCF.narrowPeak
# 38907 mESC_d6_WT_CTCF.narrowPeak

# 
for i in `ls -1 mESC*CTCF.narrowPeak|sort`;do echo $i; awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2+$10,$2+$10+1}' $i |wc -l;done
# mESC_d0_WT_CTCF.narrowPeak
# 50544
# mESC_d2_WT_CTCF.narrowPeak
# 28732
# mESC_d4_WT_CTCF.narrowPeak
# 37518
# mESC_d6_WT_CTCF.narrowPeak
# 38907
for i in `ls -1 mESC*CTCF.narrowPeak|sort`;do echo $i; awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2+$10,$2+$10+1,$7}' $i |wc -l;done
# mESC_d0_WT_CTCF.narrowPeak
# 50544
# mESC_d2_WT_CTCF.narrowPeak
# 28732
# mESC_d4_WT_CTCF.narrowPeak
# 37518
# mESC_d6_WT_CTCF.narrowPeak
# 38907
for i in `ls -1 mESC*CTCF.narrowPeak|sort`;do echo $i; awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2+$10,$2+$10+1,$7}' $i |sort -k1,1 -k2,2n |uniq |wc -l;done
# mESC_d0_WT_CTCF.narrowPeak
# 50544
# mESC_d2_WT_CTCF.narrowPeak
# 28732
# mESC_d4_WT_CTCF.narrowPeak
# 37518
# mESC_d6_WT_CTCF.narrowPeak
# 38907

# 使用如下命令构造 summit bed3 + fold_enrichment 
for i in `ls -1 mESC*CTCF.narrowPeak|sort`;do echo $i; awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2+$10,$2+$10+1,$7}' $i |sort -k1,1 -k2,2n >${i%.narrowPeak}.summit.sorted.bed;done
wc -l mESC_*_CTCF.summit.sorted.bed
# 50544 mESC_d0_WT_CTCF.summit.sorted.bed
# 28732 mESC_d2_WT_CTCF.summit.sorted.bed
# 37518 mESC_d4_WT_CTCF.summit.sorted.bed
# 38907 mESC_d6_WT_CTCF.summit.sorted.bed

awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2+$10-25,$2+$10+1+25,$7,"d0"}' mESC_d0_WT_CTCF.narrowPeak |head
awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2+$10-25,$2+$10+1+25,$7,"d0"}' mESC_d0_WT_CTCF.narrowPeak |sort -k1,1 -k2,2n >mESC_d0_WT_CTCF.summit.eB25.sorted.bed
awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2+$10-25,$2+$10+1+25,$7,"d2"}' mESC_d2_WT_CTCF.narrowPeak |sort -k1,1 -k2,2n >mESC_d2_WT_CTCF.summit.eB25.sorted.bed
awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2+$10-25,$2+$10+1+25,$7,"d4"}' mESC_d4_WT_CTCF.narrowPeak |sort -k1,1 -k2,2n >mESC_d4_WT_CTCF.summit.eB25.sorted.bed
awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2+$10-25,$2+$10+1+25,$7,"d6"}' mESC_d6_WT_CTCF.narrowPeak |sort -k1,1 -k2,2n >mESC_d6_WT_CTCF.summit.eB25.sorted.bed

ls -1 mESC_*_CTCF.summit.eB25.sorted.bed |sort
cat `ls -1 mESC_*_CTCF.summit.eB25.sorted.bed |sort` |wc -l
# 155701
cat `ls -1 mESC_*_CTCF.summit.eB25.sorted.bed |sort` |sort -k1,1 -k2,2n |bedtools merge -c 4 -o mean -i stdin |head
cat `ls -1 mESC_*_CTCF.summit.eB25.sorted.bed |sort` |sort -k1,1 -k2,2n |bedtools merge -c 4 -o mean -i stdin |awk 'BEGIN{FS="\t";OFS="\t"}{ct=int(($2+$3)/2);print $1,ct,ct+1}' |head
cat `ls -1 mESC_*_CTCF.summit.eB25.sorted.bed |sort` |sort -k1,1 -k2,2n |bedtools merge -c 4 -o mean -i stdin |awk 'BEGIN{FS="\t";OFS="\t"}{ct=int(($2+$3)/2);print $1,ct,ct+1}' |wc -l
# 58252
cat `ls -1 mESC_*_CTCF.summit.eB25.sorted.bed |sort` |sort -k1,1 -k2,2n |bedtools merge -c 4 -o mean -i stdin |awk 'BEGIN{FS="\t";OFS="\t"}{ct=int(($2+$3)/2);print $1,ct,ct+1}' |sort -k1,1 -k2,2n |uniq |wc -l
# 58252
cat `ls -1 mESC_*_CTCF.summit.eB25.sorted.bed |sort` |sort -k1,1 -k2,2n |bedtools merge -c 4 -o mean -i stdin |awk 'BEGIN{FS="\t";OFS="\t"}{ct=int(($2+$3)/2);print $1,ct,ct+1,$4}' |head
cat `ls -1 mESC_*_CTCF.summit.eB25.sorted.bed |sort` |sort -k1,1 -k2,2n |bedtools merge -c 4 -o mean -i stdin |awk 'BEGIN{FS="\t";OFS="\t"}{ct=int(($2+$3)/2);print $1,ct,ct+1,$4}' |sort -k1,1 -k2,2n >mESC_WT_d0246_union_CTCF.summit.sorted.bed
wc -l mESC_WT_d0246_union_CTCF.summit.sorted.bed
# 58252 mESC_WT_d0246_union_CTCF.summit.sorted.bed

```

### pileup use 5hic r20k data

```bash
# 
# 8.1_coolpuppy_pileup_plot_for_mESC_WT_CTCF_BS.use_Nora_WN_Ren_mcool.r20k.ipynb

```

### pileup use 5hic r10k data

```bash
# 
# 8.1_coolpuppy_pileup_plot_for_mESC_WT_CTCF_BS.use_Nora_WN_Ren_mcool.r10k.ipynb

```

### pileup use 5hic r5k data

```bash
# 
# 8.1_coolpuppy_pileup_plot_for_mESC_WT_CTCF_BS.use_Nora_WN_Ren_mcool.r5k.ipynb

```



## WT d0246 mESC MBD3 summit flank pileup

### 构造 MBD3 summit bed3+FC

```bash
# mESC WT MBD3 narrowPeak 
# /home/sunwenju/projects_on_940xa/2022121401_WN_mESC_d0246_CTCF_MBD3_ChIP_Seq/macs2_peak_calling_2_for_mESC_d0246_WT_MBD3_rep123_together/
cd /home/sunwenju/projects_on_940xa/2022121401_WN_mESC_d0246_CTCF_MBD3_ChIP_Seq/macs2_peak_calling_2_for_mESC_d0246_WT_MBD3_rep123_together/
find . -name *narrowPeak |sort 
# ./mESC_d0_WT_IP_MBD3_rep123_TvsC/mESC_d0_WT_IP_MBD3_rep123_TvsC_peaks.narrowPeak
# ./mESC_d2_WT_IP_MBD3_rep123_TvsC/mESC_d2_WT_IP_MBD3_rep123_TvsC_peaks.narrowPeak
# ./mESC_d4_WT_IP_MBD3_rep123_TvsC/mESC_d4_WT_IP_MBD3_rep123_TvsC_peaks.narrowPeak
# ./mESC_d6_WT_IP_MBD3_rep123_TvsC/mESC_d6_WT_IP_MBD3_rep123_TvsC_peaks.narrowPeak
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline
ln -s /home/sunwenju/projects_on_940xa/2022121401_WN_mESC_d0246_CTCF_MBD3_ChIP_Seq/macs2_peak_calling_2_for_mESC_d0246_WT_MBD3_rep123_together/mESC_d0_WT_IP_MBD3_rep123_TvsC/mESC_d0_WT_IP_MBD3_rep123_TvsC_peaks.narrowPeak mESC_d0_WT_MBD3.narrowPeak
ln -s /home/sunwenju/projects_on_940xa/2022121401_WN_mESC_d0246_CTCF_MBD3_ChIP_Seq/macs2_peak_calling_2_for_mESC_d0246_WT_MBD3_rep123_together/mESC_d2_WT_IP_MBD3_rep123_TvsC/mESC_d2_WT_IP_MBD3_rep123_TvsC_peaks.narrowPeak mESC_d2_WT_MBD3.narrowPeak
ln -s /home/sunwenju/projects_on_940xa/2022121401_WN_mESC_d0246_CTCF_MBD3_ChIP_Seq/macs2_peak_calling_2_for_mESC_d0246_WT_MBD3_rep123_together/mESC_d4_WT_IP_MBD3_rep123_TvsC/mESC_d4_WT_IP_MBD3_rep123_TvsC_peaks.narrowPeak mESC_d4_WT_MBD3.narrowPeak
ln -s /home/sunwenju/projects_on_940xa/2022121401_WN_mESC_d0246_CTCF_MBD3_ChIP_Seq/macs2_peak_calling_2_for_mESC_d0246_WT_MBD3_rep123_together/mESC_d6_WT_IP_MBD3_rep123_TvsC/mESC_d6_WT_IP_MBD3_rep123_TvsC_peaks.narrowPeak mESC_d6_WT_MBD3.narrowPeak
wc -l mESC*MBD3.narrowPeak
# 10915 mESC_d0_WT_MBD3.narrowPeak
#  2252 mESC_d2_WT_MBD3.narrowPeak
#  9437 mESC_d4_WT_MBD3.narrowPeak
# 11230 mESC_d6_WT_MBD3.narrowPeak

# 
for i in `ls -1 mESC*MBD3.narrowPeak|sort`;do echo $i; awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2+$10,$2+$10+1}' $i |wc -l;done
# mESC_d0_WT_MBD3.narrowPeak
# 10915
# mESC_d2_WT_MBD3.narrowPeak
# 2252
# mESC_d4_WT_MBD3.narrowPeak
# 9437
# mESC_d6_WT_MBD3.narrowPeak
# 11230
for i in `ls -1 mESC*MBD3.narrowPeak|sort`;do echo $i; awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2+$10,$2+$10+1,$7}' $i |wc -l;done
# mESC_d0_WT_MBD3.narrowPeak
# 10915
# mESC_d2_WT_MBD3.narrowPeak
# 2252
# mESC_d4_WT_MBD3.narrowPeak
# 9437
# mESC_d6_WT_MBD3.narrowPeak
# 11230
for i in `ls -1 mESC*MBD3.narrowPeak|sort`;do echo $i; awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2+$10,$2+$10+1,$7}' $i |sort -k1,1 -k2,2n |uniq |wc -l;done
# mESC_d0_WT_MBD3.narrowPeak
# 10915
# mESC_d2_WT_MBD3.narrowPeak
# 2252
# mESC_d4_WT_MBD3.narrowPeak
# 9437
# mESC_d6_WT_MBD3.narrowPeak
# 11230

# summit bed3 + fold_enrichment 
for i in `ls -1 mESC*MBD3.narrowPeak|sort`;do echo $i; awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2+$10,$2+$10+1,$7}' $i |sort -k1,1 -k2,2n >${i%.narrowPeak}.summit.sorted.bed;done
wc -l mESC_*_MBD3.summit.sorted.bed
# 10915 mESC_d0_WT_MBD3.summit.sorted.bed
#  2252 mESC_d2_WT_MBD3.summit.sorted.bed
#  9437 mESC_d4_WT_MBD3.summit.sorted.bed
# 11230 mESC_d6_WT_MBD3.summit.sorted.bed

awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2+$10-25,$2+$10+1+25,$7,"d0"}' mESC_d0_WT_MBD3.narrowPeak |head
awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2+$10-25,$2+$10+1+25,$7,"d0"}' mESC_d0_WT_MBD3.narrowPeak |sort -k1,1 -k2,2n >mESC_d0_WT_MBD3.summit.eB25.sorted.bed
awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2+$10-25,$2+$10+1+25,$7,"d2"}' mESC_d2_WT_MBD3.narrowPeak |sort -k1,1 -k2,2n >mESC_d2_WT_MBD3.summit.eB25.sorted.bed
awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2+$10-25,$2+$10+1+25,$7,"d4"}' mESC_d4_WT_MBD3.narrowPeak |sort -k1,1 -k2,2n >mESC_d4_WT_MBD3.summit.eB25.sorted.bed
awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2+$10-25,$2+$10+1+25,$7,"d6"}' mESC_d6_WT_MBD3.narrowPeak |sort -k1,1 -k2,2n >mESC_d6_WT_MBD3.summit.eB25.sorted.bed

ls -1 mESC_*_MBD3.summit.eB25.sorted.bed |sort
cat `ls -1 mESC_*_MBD3.summit.eB25.sorted.bed |sort` |wc -l
# 33834
cat `ls -1 mESC_*_MBD3.summit.eB25.sorted.bed |sort` |sort -k1,1 -k2,2n |bedtools merge -c 4 -o mean -i stdin |head
cat `ls -1 mESC_*_MBD3.summit.eB25.sorted.bed |sort` |sort -k1,1 -k2,2n |bedtools merge -c 4 -o mean -i stdin |awk 'BEGIN{FS="\t";OFS="\t"}{ct=int(($2+$3)/2);print $1,ct,ct+1}' |head
cat `ls -1 mESC_*_MBD3.summit.eB25.sorted.bed |sort` |sort -k1,1 -k2,2n |bedtools merge -c 4 -o mean -i stdin |awk 'BEGIN{FS="\t";OFS="\t"}{ct=int(($2+$3)/2);print $1,ct,ct+1}' |wc -l
# 20043
cat `ls -1 mESC_*_MBD3.summit.eB25.sorted.bed |sort` |sort -k1,1 -k2,2n |bedtools merge -c 4 -o mean -i stdin |awk 'BEGIN{FS="\t";OFS="\t"}{ct=int(($2+$3)/2);print $1,ct,ct+1}' |sort -k1,1 -k2,2n |uniq |wc -l
# 20043
cat `ls -1 mESC_*_MBD3.summit.eB25.sorted.bed |sort` |sort -k1,1 -k2,2n |bedtools merge -c 4 -o mean -i stdin |awk 'BEGIN{FS="\t";OFS="\t"}{ct=int(($2+$3)/2);print $1,ct,ct+1,$4}' |head
cat `ls -1 mESC_*_MBD3.summit.eB25.sorted.bed |sort` |sort -k1,1 -k2,2n |bedtools merge -c 4 -o mean -i stdin |awk 'BEGIN{FS="\t";OFS="\t"}{ct=int(($2+$3)/2);print $1,ct,ct+1,$4}' |sort -k1,1 -k2,2n >mESC_WT_d0246_union_MBD3.summit.sorted.bed
wc -l mESC_WT_d0246_union_MBD3.summit.sorted.bed
# 20043 mESC_WT_d0246_union_MBD3.summit.sorted.bed

```

### pileup use 5hic r20k data

```bash
# 
# 8.1_coolpuppy_pileup_plot_for_mESC_WT_MBD3_BS.use_Nora_WN_Ren_mcool.r20k.ipynb

```

### pileup use 5hic r10k data

```bash
# 
# 8.1_coolpuppy_pileup_plot_for_mESC_WT_MBD3_BS.use_Nora_WN_Ren_mcool.r10k.ipynb

```

### pileup use 5hic r5k data

```bash
# 
# 8.1_coolpuppy_pileup_plot_for_mESC_WT_MBD3_BS.use_Nora_WN_Ren_mcool.r5k.ipynb

```



## Nora WT loop anchors pileup

```bash
#
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline
bedtools slop -b 150 -i mESC_WT_d0246_union_CTCF.summit.sorted.bed -g /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes |sort -k1,1 -k2,2n >mESC_WT_d0246_union_CTCF.summit.eB150.sorted.bed

sed '1d' nora_wt_r20k_loops_df.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"loop_"NR,"left"; print $4,$5,$6,"loop_"NR,"right"}' | bedtools intersect -a stdin -b mESC_WT_d0246_union_CTCF.summit.eB150.sorted.bed -wao |head
# chr1	4500000	4520000	loop_1	left	chr1	4516597	4516898	9.385175	301
# chr1	4760000	4780000	loop_1	right	chr1	4768432	4768733	12.067015	301
# chr1	4760000	4780000	loop_1	right	chr1	4769914	4770215	6.0398025	301
# chr1	5220000	5240000	loop_2	left	.	-1	-1	.	0
# chr1	5380000	5400000	loop_2	right	.	-1	-1	.	0

sed '1d' nora_wt_r20k_loops_df.txt |wc -l
# 2628

# CBS_type:O
sed '1d' nora_wt_r20k_loops_df.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"loop_"NR,"left"; print $4,$5,$6,"loop_"NR,"right"}' | bedtools intersect -a stdin -b mESC_WT_d0246_union_CTCF.summit.eB150.sorted.bed -wao |awk '{if($10==0){print $4}}' |head -n 5
# loop_2
# loop_2
# loop_3
# loop_3
# loop_4
sed '1d' nora_wt_r20k_loops_df.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"loop_"NR,"left"; print $4,$5,$6,"loop_"NR,"right"}' | bedtools intersect -a stdin -b mESC_WT_d0246_union_CTCF.summit.eB150.sorted.bed -wao |awk '{if($10==0){print $4}}' |sort |uniq |wc -l
# 1365

# CBS_type:LRB
sed '1d' nora_wt_r20k_loops_df.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"loop_"NR,"left"; print $4,$5,$6,"loop_"NR,"right"}' | bedtools intersect -a stdin -b mESC_WT_d0246_union_CTCF.summit.eB150.sorted.bed -wao |awk '{if($10!=0){print $4}}' |sort |uniq |wc -l
# 1781

# CBS_type:L|R
# 1365 + 1781 - 2628 = 518
# CBS_type: B   1781 - 518 = 1263

sed '1d' nora_wt_r20k_loops_df.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"loop_"NR,"left"}' | bedtools intersect -a stdin -b mESC_WT_d0246_union_CTCF.summit.eB150.sorted.bed -wao |awk '{if($10!=0){print $4}}' |sort |uniq |wc -l
# 1541
# CBS_type:L 1541 - 1263 = 278
sed '1d' nora_wt_r20k_loops_df.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $4,$5,$6,"loop_"NR,"right"}' | bedtools intersect -a stdin -b mESC_WT_d0246_union_CTCF.summit.eB150.sorted.bed -wao |awk '{if($10!=0){print $4}}' |sort |uniq |wc -l
# 1503
# CBS_type:R 1503 - 1263 = 240

# CBS_type:O  2628 - 1263 - 278 - 240 = 847

```

### r20k / r10k / r5k loop anchors pileup analysis

```bash
# 
# 8.3_coolpuppy_pileup_plot_for_Nora_WT_r20k_loop_anchors.use_Nora_WN_Ren_mcool.r20k.ipynb
# 8.3_coolpuppy_pileup_plot_for_Nora_WT_r10k_loop_anchors.use_Nora_WN_Ren_mcool.r10k.ipynb
# 8.3_coolpuppy_pileup_plot_for_Nora_WT_r5k_loop_anchors.use_Nora_WN_Ren_mcool.r5k.ipynb

```





## Ren WT loop anchors pileup

```bash
#
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline
bedtools slop -b 150 -i mESC_WT_d0246_union_CTCF.summit.sorted.bed -g /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes |sort -k1,1 -k2,2n >mESC_WT_d0246_union_CTCF.summit.eB150.sorted.bed

sed '1d' ren_wt_r20k_loops_df.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"loop_"NR,"left"; print $4,$5,$6,"loop_"NR,"right"}' | bedtools intersect -a stdin -b mESC_WT_d0246_union_CTCF.summit.eB150.sorted.bed -wao |head
# chr1	4500000	4520000	loop_1	left	chr1	4516597	4516898	9.385175	301
# chr1	4760000	4780000	loop_1	right	chr1	4768432	4768733	12.067015	301
# chr1	4760000	4780000	loop_1	right	chr1	4769914	4770215	6.0398025	301
# chr1	5220000	5240000	loop_2	left	.	-1	-1	.	0
# chr1	5380000	5400000	loop_2	right	.	-1	-1	.	0

sed '1d' ren_wt_r20k_loops_df.txt |wc -l
# 2427

# CBS_type:O
sed '1d' ren_wt_r20k_loops_df.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"loop_"NR,"left"; print $4,$5,$6,"loop_"NR,"right"}' | bedtools intersect -a stdin -b mESC_WT_d0246_union_CTCF.summit.eB150.sorted.bed -wao |awk '{if($10==0){print $4}}' |head -n 5
# loop_6
# loop_18
# loop_29
# loop_32
# loop_33
sed '1d' ren_wt_r20k_loops_df.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"loop_"NR,"left"; print $4,$5,$6,"loop_"NR,"right"}' | bedtools intersect -a stdin -b mESC_WT_d0246_union_CTCF.summit.eB150.sorted.bed -wao |awk '{if($10==0){print $4}}' |sort |uniq |wc -l
# 375

# CBS_type:LRB
sed '1d' ren_wt_r20k_loops_df.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"loop_"NR,"left"; print $4,$5,$6,"loop_"NR,"right"}' | bedtools intersect -a stdin -b mESC_WT_d0246_union_CTCF.summit.eB150.sorted.bed -wao |awk '{if($10!=0){print $4}}' |sort |uniq |wc -l
# 2334

# CBS_type:L|R
# 375 + 2334 - 2427 = 282
# CBS_type: B   2334 - 282 = 2052

sed '1d' ren_wt_r20k_loops_df.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"loop_"NR,"left"}' | bedtools intersect -a stdin -b mESC_WT_d0246_union_CTCF.summit.eB150.sorted.bed -wao |awk '{if($10!=0){print $4}}' |sort |uniq |wc -l
# 2203
# CBS_type:L 2203 - 2052 = 151
sed '1d' ren_wt_r20k_loops_df.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $4,$5,$6,"loop_"NR,"right"}' | bedtools intersect -a stdin -b mESC_WT_d0246_union_CTCF.summit.eB150.sorted.bed -wao |awk '{if($10!=0){print $4}}' |sort |uniq |wc -l
# 2183
# CBS_type:R 2183 - 2052 = 131

# CBS_type:O  2427 - 2052 - 151 - 131 = 93

```
### r20k / r10k / r5k loop anchors pileup analysis

```bash
# 
# 8.4_coolpuppy_pileup_plot_for_Ren_WT_r20k_loop_anchors.use_Nora_WN_Ren_mcool.r20k.ipynb
# 8.4_coolpuppy_pileup_plot_for_Ren_WT_r10k_loop_anchors.use_Nora_WN_Ren_mcool.r10k.ipynb
# 8.4_coolpuppy_pileup_plot_for_Ren_WT_r5k_loop_anchors.use_Nora_WN_Ren_mcool.r5k.ipynb

```



## Nora WT insulation boundaries pileup

```bash
# 8.5_coolpuppy_pileup_plot_for_Nora_WT_r20k_insulation_boundaries.use_Nora_WN_Ren_mcool.r20k.ipynb
# 8.5_coolpuppy_pileup_plot_for_Nora_WT_r10k_insulation_boundaries.use_Nora_WN_Ren_mcool.r10k.ipynb
# 8.5_coolpuppy_pileup_plot_for_Nora_WT_r5k_insulation_boundaries.use_Nora_WN_Ren_mcool.r5k.ipynb

cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_WT
pwd
# /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Nora_WT
# r10k 
wc -l Nora_mESC_WT_HiC.rep12.merged.mcool_r10k.w15r_insulation_boundary.bed
# 6482 Nora_mESC_WT_HiC.rep12.merged.mcool_r10k.w15r_insulation_boundary.bed
bedtools intersect -a Nora_mESC_WT_HiC.rep12.merged.mcool_r10k.w15r_insulation_boundary.bed -b ../mESC_WT_d0246_union_CTCF.summit.sorted.bed -wao |awk '{if($8!=0){print $1,$2,$3}}' |sort -k1,2 -k2,2n |uniq |wc -l
# 2155  (with CTCF BS)
bedtools intersect -a Nora_mESC_WT_HiC.rep12.merged.mcool_r10k.w15r_insulation_boundary.bed -b ../mESC_WT_d0246_union_MBD3.summit.sorted.bed -wao |awk '{if($8!=0){print $1,$2,$3}}' |sort -k1,2 -k2,2n |uniq |wc -l
# 1044  (with MBD3 BS)

# r20k
wc -l Nora_mESC_WT_HiC.rep12.merged.mcool_r20k.w5r_insulation_boundary.bed
# 6830 Nora_mESC_WT_HiC.rep12.merged.mcool_r20k.w5r_insulation_boundary.bed
bedtools intersect -a Nora_mESC_WT_HiC.rep12.merged.mcool_r20k.w5r_insulation_boundary.bed -b ../mESC_WT_d0246_union_CTCF.summit.sorted.bed -wao |awk '{if($8!=0){print $1,$2,$3}}' |sort -k1,2 -k2,2n |uniq |wc -l
# 3680  (with CTCF BS)
bedtools intersect -a Nora_mESC_WT_HiC.rep12.merged.mcool_r20k.w5r_insulation_boundary.bed -b ../mESC_WT_d0246_union_MBD3.summit.sorted.bed -wao |awk '{if($8!=0){print $1,$2,$3}}' |sort -k1,2 -k2,2n |uniq |wc -l
# 1956  (with MBD3 BS)


# Nroa WT 
# r10k: 6482 boundaries;  2155 with CTCF BS; 1044 with MBD3 BS;
# r20k: 6830 boundaries;  3680 with CTCF BS; 1956 with MBD3 BS;

```

## Ren WT insulation boundaries pileup

```bash
# 8.6_coolpuppy_pileup_plot_for_Ren_WT_r20k_insulation_boundaries.use_Nora_WN_Ren_mcool.r20k.ipynb
# 8.6_coolpuppy_pileup_plot_for_Ren_WT_r10k_insulation_boundaries.use_Nora_WN_Ren_mcool.r10k.ipynb
# 8.6_coolpuppy_pileup_plot_for_Ren_WT_r5k_insulation_boundaries.use_Nora_WN_Ren_mcool.r5k.ipynb

cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_WT
pwd
# /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline/Ren_WT
# r10k
wc -l Ren_mESC_WT_HiC.rep12.merged.mcool_r10k.w15r_insulation_boundary.bed
# 6184 Ren_mESC_WT_HiC.rep12.merged.mcool_r10k.w15r_insulation_boundary.bed
bedtools intersect -a Ren_mESC_WT_HiC.rep12.merged.mcool_r10k.w15r_insulation_boundary.bed -b ../mESC_WT_d0246_union_CTCF.summit.sorted.bed -wao |awk '{if($8!=0){print $1,$2,$3}}' |sort -k1,2 -k2,2n |uniq |wc -l
# 2438  (with CTCF BS)
bedtools intersect -a Ren_mESC_WT_HiC.rep12.merged.mcool_r10k.w15r_insulation_boundary.bed -b ../mESC_WT_d0246_union_MBD3.summit.sorted.bed -wao |awk '{if($8!=0){print $1,$2,$3}}' |sort -k1,2 -k2,2n |uniq |wc -l
# 1145  (with MBD3 BS)

# r20k
wc -l Ren_mESC_WT_HiC.rep12.merged.mcool_r20k.w5r_insulation_boundary.bed
# 6178 Ren_mESC_WT_HiC.rep12.merged.mcool_r20k.w5r_insulation_boundary.bed
bedtools intersect -a Ren_mESC_WT_HiC.rep12.merged.mcool_r20k.w5r_insulation_boundary.bed -b ../mESC_WT_d0246_union_CTCF.summit.sorted.bed -wao |awk '{if($8!=0){print $1,$2,$3}}' |sort -k1,2 -k2,2n |uniq |wc -l
# 3789  (with CTCF BS)
bedtools intersect -a Ren_mESC_WT_HiC.rep12.merged.mcool_r20k.w5r_insulation_boundary.bed -b ../mESC_WT_d0246_union_MBD3.summit.sorted.bed -wao |awk '{if($8!=0){print $1,$2,$3}}' |sort -k1,2 -k2,2n |uniq |wc -l
# 2048  (with MBD3 BS)

# Ren WT 
# r10k: 6184 boundaries;  2438 with CTCF BS; 1145 with MBD3 BS;
# r20k: 6178 boundaries;  3789 with CTCF BS; 2048 with MBD3 BS;

```

## Nora IAA insulation boundaries pileup

```bash
# 8.7_coolpuppy_pileup_plot_for_Nora_IAA_r20k_insulation_boundaries.use_Nora_WN_Ren_mcool.r20k.ipynb
# 8.7_coolpuppy_pileup_plot_for_Nora_IAA_r10k_insulation_boundaries.use_Nora_WN_Ren_mcool.r10k.ipynb
# 8.7_coolpuppy_pileup_plot_for_Nora_IAA_r5k_insulation_boundaries.use_Nora_WN_Ren_mcool.r5k.ipynb

```
## Ren IAA insulation boundaries pileup

```bash
# 8.8_coolpuppy_pileup_plot_for_Ren_IAA_r20k_insulation_boundaries.use_Nora_WN_Ren_mcool.r20k.ipynb
# 8.8_coolpuppy_pileup_plot_for_Ren_IAA_r10k_insulation_boundaries.use_Nora_WN_Ren_mcool.r10k.ipynb
# 8.8_coolpuppy_pileup_plot_for_Ren_IAA_r5k_insulation_boundaries.use_Nora_WN_Ren_mcool.r5k.ipynb

```
## WN BPK25 insulation boundaries pileup

```bash
# 8.9_coolpuppy_pileup_plot_for_WN_BPK25_r20k_insulation_boundaries.use_Nora_WN_Ren_mcool.r20k.ipynb
# 8.9_coolpuppy_pileup_plot_for_WN_BPK25_r10k_insulation_boundaries.use_Nora_WN_Ren_mcool.r10k.ipynb
# 8.9_coolpuppy_pileup_plot_for_WN_BPK25_r5k_insulation_boundaries.use_Nora_WN_Ren_mcool.r5k.ipynb

```




## Nora WT TADs pileup

```bash
# 8.10_coolpuppy_pileup_plot_for_Nora_WT_r20k_TADs.use_Nora_WN_Ren_mcool.r20k.ipynb
# 8.10_coolpuppy_pileup_plot_for_Nora_WT_r10k_TADs.use_Nora_WN_Ren_mcool.r10k.ipynb
# 8.10_coolpuppy_pileup_plot_for_Nora_WT_r5k_TADs.use_Nora_WN_Ren_mcool.r5k.ipynb

```

## Ren WT TADs pileup

```bash
# 8.11_coolpuppy_pileup_plot_for_Ren_WT_r20k_TADs.use_Nora_WN_Ren_mcool.r20k.ipynb
# 8.11_coolpuppy_pileup_plot_for_Ren_WT_r10k_TADs.use_Nora_WN_Ren_mcool.r10k.ipynb
# 8.11_coolpuppy_pileup_plot_for_Ren_WT_r5k_TADs.use_Nora_WN_Ren_mcool.r5k.ipynb

```
## Nora IAA TADs pileup

```bash
# 8.12_coolpuppy_pileup_plot_for_Nora_IAA_r20k_TADs.use_Nora_WN_Ren_mcool.r20k.ipynb
# 8.12_coolpuppy_pileup_plot_for_Nora_IAA_r10k_TADs.use_Nora_WN_Ren_mcool.r10k.ipynb
# 8.12_coolpuppy_pileup_plot_for_Nora_IAA_r5k_TADs.use_Nora_WN_Ren_mcool.r5k.ipynb

```

## Ren IAA TADs pileup

```bash
# 8.13_coolpuppy_pileup_plot_for_Ren_IAA_r20k_TADs.use_Nora_WN_Ren_mcool.r20k.ipynb
# 8.13_coolpuppy_pileup_plot_for_Ren_IAA_r10k_TADs.use_Nora_WN_Ren_mcool.r10k.ipynb
# 8.13_coolpuppy_pileup_plot_for_Ren_IAA_r5k_TADs.use_Nora_WN_Ren_mcool.r5k.ipynb

```
## WN BPK25 TADs pileup

```bash
# 8.14_coolpuppy_pileup_plot_for_WN_BPK25_r20k_TADs.use_Nora_WN_Ren_mcool.r20k.ipynb
# 8.14_coolpuppy_pileup_plot_for_WN_BPK25_r10k_TADs.use_Nora_WN_Ren_mcool.r10k.ipynb
# 8.14_coolpuppy_pileup_plot_for_WN_BPK25_r5k_TADs.use_Nora_WN_Ren_mcool.r5k.ipynb

```



## WT d0246 mESC CTCF & MBD3 intervene 3types summits flank pileup (20240920)

### CTCF & MBD3  eB150 intervene

```bash
cd /data1/sunwenju/projects_on_940xa/2024041703_mESC_HiC_use_pairtools_cooler_pipeline
ls -1 mESC_WT_d0246_*.bed
# mESC_WT_d0246_union_CTCF.summit.eB150.sorted.bed
# mESC_WT_d0246_union_CTCF.summit.sorted.bed
# mESC_WT_d0246_union_MBD3.summit.sorted.bed
wc -l  mESC_WT_d0246_*.bed
# 58252 mESC_WT_d0246_union_CTCF.summit.eB150.sorted.bed
# 58252 mESC_WT_d0246_union_CTCF.summit.sorted.bed
# 20043 mESC_WT_d0246_union_MBD3.summit.sorted.bed

# 
bedtools slop -b 150 -i mESC_WT_d0246_union_MBD3.summit.sorted.bed -g /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes |sort -k1,1 -k2,2n |uniq >mESC_WT_d0246_union_MBD3.summit.eB150.sorted.bed
wc -l  mESC_WT_d0246_*.sorted.bed
# 58252 mESC_WT_d0246_union_CTCF.summit.eB150.sorted.bed
# 58252 mESC_WT_d0246_union_CTCF.summit.sorted.bed
# 20043 mESC_WT_d0246_union_MBD3.summit.eB150.sorted.bed
# 20043 mESC_WT_d0246_union_MBD3.summit.sorted.bed

# intervene 
intervene venn --type genomic -i \
 mESC_WT_d0246_union_CTCF.summit.eB150.sorted.bed \
 mESC_WT_d0246_union_MBD3.summit.eB150.sorted.bed \
 --names=mESC_WT_CTCF_PseB150,mESC_WT_MBD3_PseB150 \
 --output mESC_WT_CTCF_MBD3_PSeB150_intervene_venn \
 --project mESC_WT_CTCF_MBD3_PSeB150_intervene_venn \
 --title mESC_WT_CTCF_MBD3_PSeB150_intervene_venn \
 --save-overlaps \
 --figtype pdf \
 --figsize 6 6

# 3c summits
awk 'BEGIN{FS="\t";OFS="\t"}{center=int(($2+$3)/2);print $1,center,center+1,$4}' 01_mESC_WT_MBD3_PseB150.bed |head
# chr1	4435638	4435639	4.7787
# chr1	5296835	5296836	4.15304
awk 'BEGIN{FS="\t";OFS="\t"}{center=int(($2+$3)/2);print $1,center,center+1,$4}' ./mESC_WT_CTCF_MBD3_PSeB150_intervene_venn/sets/01_mESC_WT_MBD3_PseB150.bed |sort -k1,1 -k2,2n |uniq >mESC_WT_CvM3.MBD3_only_PS.sorted.bed
awk 'BEGIN{FS="\t";OFS="\t"}{center=int(($2+$3)/2);print $1,center,center+1,$4}' ./mESC_WT_CTCF_MBD3_PSeB150_intervene_venn/sets/10_mESC_WT_CTCF_PseB150.bed |sort -k1,1 -k2,2n |uniq >mESC_WT_CvM3.CTCF_only_PS.sorted.bed
awk 'BEGIN{FS="\t";OFS="\t"}{center=int(($2+$3)/2);print $1,center,center+1,$4}' ./mESC_WT_CTCF_MBD3_PSeB150_intervene_venn/sets/11_mESC_WT_CTCF_PseB150_mESC_WT_MBD3_PseB150.bed |sort -k1,1 -k2,2n |uniq >mESC_WT_CvM3.CTCF_MBD3_common_PS.sorted.bed

```

### pileup use 5hic r20k data

```bash
# 
# 8.15_coolpuppy_pileup_plot_for_mESC_WT_CvM3_CTCF_only_BS.use_Nora_WN_Ren_mcool.r10k.ipynb
# 8.15_coolpuppy_pileup_plot_for_mESC_WT_CvM3_CTCF_only_BS.use_Nora_WN_Ren_mcool.r10k.summary.pptx

```

### pileup use 5hic r10k data

```bash
# 
# 8.15_coolpuppy_pileup_plot_for_mESC_WT_CvM3_MBD3_only_BS.use_Nora_WN_Ren_mcool.r10k.ipynb
# 8.15_coolpuppy_pileup_plot_for_mESC_WT_CvM3_MBD3_only_BS.use_Nora_WN_Ren_mcool.r10k.summary.pptx

```

### pileup use 5hic r5k data

```bash
# 
# 8.15_coolpuppy_pileup_plot_for_mESC_WT_CvM3_CTCF_MBD3_common_BS.use_Nora_WN_Ren_mcool.r10k.ipynb
# 8.15_coolpuppy_pileup_plot_for_mESC_WT_CvM3_CTCF_MBD3_common_BS.use_Nora_WN_Ren_mcool.r10k.summary.pptx

# 
# 8.15_coolpuppy_pileup_plot_for_mESC_WT_CvM3_BS.use_Nora_WN_Ren_mcool.r10k.summary.pptx
```









