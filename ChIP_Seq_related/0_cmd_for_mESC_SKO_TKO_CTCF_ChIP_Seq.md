





> SKO: DNMT1a/b/2 and TET1/2/3 Six genes KO，（low methylation level, ~<20%）
>
> TKO: TET1/2/3 Three genes KO,  （high methylation level， >75%）
>
> WT: Normal mESC （Normal methylation level， 50% ± 20%）



# QC

```bash
cd /data1/sunwenju/projects_on_940xa/
mkdir 2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
# 根据样本信息表构造数据连接创建脚本，上传并运行
bash 0_make_and_rename_soft_links_for_WN_mESC_SKO_TKO_CTCF_ChIP_Seq_raw_data.sh
ls -1

# 原始数据进行质量评估
mkdir 1_all_ChIP_Seq_raw_data_fastqc_result
nohup fastqc -outdir 1_all_ChIP_Seq_raw_data_fastqc_result --threads 72 *.fq.gz >1_all_ChIP_Seq_raw_data_fastqc.runlog 2>&1 &
tail -f 1_all_ChIP_Seq_raw_data_fastqc.runlog

# 并用 multiqc 进行结果汇总
# 提取原始数据reads数量
cd 1_all_ChIP_Seq_raw_data_fastqc_result
multiqc ./*fastqc.zip --ignore *.html --no-megaqc-upload
cd ../
# PE150，有接头， illumina universal adapter

```

# cutadapt and QC

## cutadapt

```bash
# 根据 fastqc 的结果提示，在 illumina 的接头序列文件中找到如下序列，并对其中一个文件进行了grep统计，确认是接头
# Read 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
# Read 2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

# zcat grep 可以看到有接头且reads前两nt含N的reads较多所以同时增加 -u 2 -U 2 删除 reads 1/2 开始的2nt
# 使用如下命令生成去接头脚本并运行
for i in `ls -1 *rep[123].R1.fq.gz`;do echo 'cutadapt -j 128 -q 20,20 -u 2 -U 2 --trim-n -m 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o '${i%.R1.fq.gz}.R1.cutadapt.fq.gz' -p '${i%.R1.fq.gz}.R2.cutadapt.fq.gz' '$i' '${i%.R1.fq.gz}.R2.fq.gz;done >2_all_mESC_SKO_TKO_WT_CTCF_ChIP_Seq_raw_data_run_cutadapt.sh
cat 2_all_mESC_SKO_TKO_WT_CTCF_ChIP_Seq_raw_data_run_cutadapt.sh
nohup bash 2_all_mESC_SKO_TKO_WT_CTCF_ChIP_Seq_raw_data_run_cutadapt.sh >2_all_mESC_SKO_TKO_WT_CTCF_ChIP_Seq_raw_data_run_cutadapt.sh.runlog 2>&1 &
tail -f 2_all_mESC_SKO_TKO_WT_CTCF_ChIP_Seq_raw_data_run_cutadapt.sh.runlog

```

## QC after cutadapt

```bash
# 去接头后再次进行质量评估
mkdir 2_all_ChIP_Seq_cutadapt_fastqc_result
nohup fastqc -outdir 2_all_ChIP_Seq_cutadapt_fastqc_result --threads 128 *.cutadapt.fq.gz >2_all_ChIP_Seq_cutadapt_fastqc.runlog 2>&1 &
tail -f 2_all_ChIP_Seq_cutadapt_fastqc.runlog

# 并用 multiqc 进行结果汇总
# 提取原始数据reads数量
cd 2_all_ChIP_Seq_cutadapt_fastqc_result
multiqc ./*fastqc.zip --ignore *.html --no-megaqc-upload
cd ../

# 删除无用的原始数据链接
rm *rep[123].R[12].fq.gz

```

# Genome mapping

```bash
# 可使用如下命令构造比对脚本， 然后运行
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
for i in `ls -1 *rep[123].R1.cutadapt.fq.gz`;do echo 'bowtie2 -t -p 128 --no-unal -x /home/sunwenju/3_genomeData/mouse/mm10/bowtie2Index/mm10bt2idx -1 '$i' -2 '${i%.R1.cutadapt.fq.gz}.R2.cutadapt.fq.gz' -S '${i%.R1.cutadapt.fq.gz}.bowtie2.sam' --un-gz '${i%.R1.cutadapt.fq.gz}.notAlign.unpairedReads.fq.gz' --un-conc-gz '${i%.R1.cutadapt.fq.gz}.notAlign.pairedReads.fq.gz' >'${i%.R1.cutadapt.fq.gz}.bowtie2.runlog' 2>&1';done >3_all_mESC_SKO_TKO_WT_CTCF_ChIP_Seq_bowtie2_map.sh
cat 3_all_mESC_SKO_TKO_WT_CTCF_ChIP_Seq_bowtie2_map.sh
# 运行比对脚本
nohup bash 3_all_mESC_SKO_TKO_WT_CTCF_ChIP_Seq_bowtie2_map.sh >3_all_mESC_SKO_TKO_WT_CTCF_ChIP_Seq_bowtie2_map.sh.runlog 2>&1 &
# tail -f 3_all_mESC_SKO_TKO_WT_CTCF_ChIP_Seq_bowtie2_map.sh.runlog
tail -f `ls -1 *bowtie2.runlog|tail -n 1`

# 使用 multiqc 对 bowtie2 比对情况进行汇总
multiqc --no-megaqc-upload -m bowtie2 -o 3_all_mESC_SKO_TKO_WT_CTCF_ChIP_Seq_bowtie2_map_results_summary ./*.bowtie2.runlog

```

# Convert sam to sorted bam

## sam to sorted bam

```bash
# 将比对sam结果转换为 按座位排序的 bam 并建立索引
# for i in `ls -1 *.bowtie2.sam|sort`; do echo $i; samtools view --threads 36 -bS $i | samtools sort --threads 36 -o ${i%.sam}.sorted.bam; samtools index ${i%.sam}.sorted.bam; done
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
nohup bash 4_mESC_TKO_SKO_CTCF_ChIP_Seq_data_bowtie2_sam_to_sorted_bam.sh >4_mESC_TKO_SKO_CTCF_ChIP_Seq_data_bowtie2_sam_to_sorted_bam.sh.runlog 2>&1 &
tail -f 4_mESC_TKO_SKO_CTCF_ChIP_Seq_data_bowtie2_sam_to_sorted_bam.sh.runlog
# rm *.bowtie2.sam

```

## mapped reads stat

```bash
# 使用 samtools flagstat 统计比对情况
# 每个bam文件不到一分钟
for i in `ls -1 *bowtie2.sorted.bam`;do echo $i; samtools flagstat --threads 16 $i;done >3_mESC_TKO_SKO_CTCF_ChIP_Seq_bowtie2_map_samtools_flagstat.txt

# 统计比对到各个染色体上的 总read数量 
for i in `ls -1 *bowtie2.sorted.bam|sort`;do echo $i; samtools view $i |grep -v -e '^@' |cut -f 3 |sort |uniq -c |awk 'BEGIN{OFS="\t"}{if($0 !~ /_/){print $2,$1}}' |sort -k1,1 >${i%.sam}.chr_mapped_reads_stat.txt;done
# 统计比对到各个染色体上的 uniq read数量 
for i in `ls -1 *bowtie2.sorted.bam|sort`;do echo $i; samtools view $i |grep -v -e '^@' |cut -f 3,4 |sort |uniq |cut -f 1 |sort |uniq -c |awk 'BEGIN{OFS="\t"}{if($0 !~ /_/){print $2,$1}}' |sort -k1,1 >${i%.sam}.chr_mapped_uniq_reads_stat.txt;done

# 比对到各个染色体上的 read数量 统计结果合并
paste `ls -1 *chr_mapped_reads_stat.txt|sort` |head
# 检测各个统计结果的染色体数量和顺序是否一致
paste `ls -1 *chr_mapped_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"}{if($1!=$3 || $1!=$5 || $1!=$7 || $1!=$9 || $1!=$11 || $1!=$13 || $1!=$15 || $1!=$17 || $1!=$19 || $1!=$21 || $1!=$23 || $1!=$25 || $1!=$27 || $1!=$29 || $1!=$31 || $1!=$33 || $1!=$35){print $1,$2}}' |head
# 染色体数量和顺序一致，则合并数据列
paste `ls -1 *chr_mapped_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"; print "chrom","mESC_SKO_INPUT_rep1","mESC_SKO_INPUT_rep2","mESC_SKO_INPUT_rep3","mESC_SKO_IP_CTCF_rep1","mESC_SKO_IP_CTCF_rep2","mESC_SKO_IP_CTCF_rep3","mESC_TKO_INPUT_rep1","mESC_TKO_INPUT_rep2","mESC_TKO_INPUT_rep3","mESC_TKO_IP_CTCF_rep1","mESC_TKO_IP_CTCF_rep2","mESC_TKO_IP_CTCF_rep3","mESC_WT_INPUT_rep1","mESC_WT_INPUT_rep2","mESC_WT_INPUT_rep3","mESC_WT_IP_CTCF_rep1","mESC_WT_IP_CTCF_rep2","mESC_WT_IP_CTCF_rep3"}{if($1==$3 && $1==$5 && $1==$7 && $1==$9 && $1==$11 && $1==$13 && $1==$15 && $1==$17 && $1==$19 && $1==$21 && $1==$23 && $1==$25 && $1==$27 && $1==$29 && $1==$31 && $1==$33 && $1==$35){print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36}}' |head
paste `ls -1 *chr_mapped_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"; print "chrom","mESC_SKO_INPUT_rep1","mESC_SKO_INPUT_rep2","mESC_SKO_INPUT_rep3","mESC_SKO_IP_CTCF_rep1","mESC_SKO_IP_CTCF_rep2","mESC_SKO_IP_CTCF_rep3","mESC_TKO_INPUT_rep1","mESC_TKO_INPUT_rep2","mESC_TKO_INPUT_rep3","mESC_TKO_IP_CTCF_rep1","mESC_TKO_IP_CTCF_rep2","mESC_TKO_IP_CTCF_rep3","mESC_WT_INPUT_rep1","mESC_WT_INPUT_rep2","mESC_WT_INPUT_rep3","mESC_WT_IP_CTCF_rep1","mESC_WT_IP_CTCF_rep2","mESC_WT_IP_CTCF_rep3"}{if($1==$3 && $1==$5 && $1==$7 && $1==$9 && $1==$11 && $1==$13 && $1==$15 && $1==$17 && $1==$19 && $1==$21 && $1==$23 && $1==$25 && $1==$27 && $1==$29 && $1==$31 && $1==$33 && $1==$35){print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36}}' >all_samples_chr_mapped_reads_stat.summary.txt
cat all_samples_chr_mapped_reads_stat.summary.txt

# 比对到各个染色体上的 uniq read数量 统计结果合并
paste `ls -1 *chr_mapped_uniq_reads_stat.txt|sort` |head
# 检测各个统计结果的染色体数量和顺序是否一致
paste `ls -1 *chr_mapped_uniq_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"}{if($1!=$3 || $1!=$5 || $1!=$7 || $1!=$9 || $1!=$11 || $1!=$13 || $1!=$15 || $1!=$17 || $1!=$19 || $1!=$21 || $1!=$23 || $1!=$25 || $1!=$27 || $1!=$29 || $1!=$31 || $1!=$33 || $1!=$35){print $1,$2}}' |head
# 染色体数量和顺序一致，则合并数据列
paste `ls -1 *chr_mapped_uniq_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"; print "chrom","mESC_SKO_INPUT_rep1","mESC_SKO_INPUT_rep2","mESC_SKO_INPUT_rep3","mESC_SKO_IP_CTCF_rep1","mESC_SKO_IP_CTCF_rep2","mESC_SKO_IP_CTCF_rep3","mESC_TKO_INPUT_rep1","mESC_TKO_INPUT_rep2","mESC_TKO_INPUT_rep3","mESC_TKO_IP_CTCF_rep1","mESC_TKO_IP_CTCF_rep2","mESC_TKO_IP_CTCF_rep3","mESC_WT_INPUT_rep1","mESC_WT_INPUT_rep2","mESC_WT_INPUT_rep3","mESC_WT_IP_CTCF_rep1","mESC_WT_IP_CTCF_rep2","mESC_WT_IP_CTCF_rep3"}{if($1==$3 && $1==$5 && $1==$7 && $1==$9 && $1==$11 && $1==$13 && $1==$15 && $1==$17 && $1==$19 && $1==$21 && $1==$23 && $1==$25 && $1==$27 && $1==$29 && $1==$31 && $1==$33 && $1==$35){print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36}}' |head
paste `ls -1 *chr_mapped_uniq_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"; print "chrom","mESC_SKO_INPUT_rep1","mESC_SKO_INPUT_rep2","mESC_SKO_INPUT_rep3","mESC_SKO_IP_CTCF_rep1","mESC_SKO_IP_CTCF_rep2","mESC_SKO_IP_CTCF_rep3","mESC_TKO_INPUT_rep1","mESC_TKO_INPUT_rep2","mESC_TKO_INPUT_rep3","mESC_TKO_IP_CTCF_rep1","mESC_TKO_IP_CTCF_rep2","mESC_TKO_IP_CTCF_rep3","mESC_WT_INPUT_rep1","mESC_WT_INPUT_rep2","mESC_WT_INPUT_rep3","mESC_WT_IP_CTCF_rep1","mESC_WT_IP_CTCF_rep2","mESC_WT_IP_CTCF_rep3"}{if($1==$3 && $1==$5 && $1==$7 && $1==$9 && $1==$11 && $1==$13 && $1==$15 && $1==$17 && $1==$19 && $1==$21 && $1==$23 && $1==$25 && $1==$27 && $1==$29 && $1==$31 && $1==$33 && $1==$35){print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36}}' >all_samples_chr_mapped_uniq_reads_stat.summary.txt
cat all_samples_chr_mapped_uniq_reads_stat.summary.txt

```

## rep1/2/3 bam merge

```bash
# 绘制reads 分布热图时要用到 rep123合并的bam文件
# 创建 bam merge 脚本 4.3_make_all_mESC_TKO_SKO_CTCF_ChIP_Seq_rep123_merged_bam_and_index.sh 内容如下：
# # 合并各重复的 bam
# samtools merge --threads 72 mESC_SKO_INPUT.rep123.merged.sorted.bam mESC_SKO_INPUT_rep1.bowtie2.sorted.bam mESC_SKO_INPUT_rep2.bowtie2.sorted.bam mESC_SKO_INPUT_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 mESC_SKO_IP_CTCF.rep123.merged.sorted.bam mESC_SKO_IP_CTCF_rep1.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep2.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 mESC_TKO_INPUT.rep123.merged.sorted.bam mESC_TKO_INPUT_rep1.bowtie2.sorted.bam mESC_TKO_INPUT_rep2.bowtie2.sorted.bam mESC_TKO_INPUT_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 mESC_TKO_IP_CTCF.rep123.merged.sorted.bam mESC_TKO_IP_CTCF_rep1.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep2.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 mESC_WT_INPUT.rep123.merged.sorted.bam mESC_WT_INPUT_rep1.bowtie2.sorted.bam mESC_WT_INPUT_rep2.bowtie2.sorted.bam mESC_WT_INPUT_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 mESC_WT_IP_CTCF.rep123.merged.sorted.bam mESC_WT_IP_CTCF_rep1.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep2.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep3.bowtie2.sorted.bam
# # 为所有 merged bam 建立索引
# for i in `ls -1 *.rep123.merged.sorted.bam`; do echo $i; samtools index -@ 72 $i;done

cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
nohup bash 4.3_make_all_mESC_TKO_SKO_CTCF_ChIP_Seq_rep123_merged_bam_and_index.sh >4.3_make_all_mESC_TKO_SKO_CTCF_ChIP_Seq_rep123_merged_bam_and_index.sh.runlog 2>&1 &
tail -f 4.3_make_all_mESC_TKO_SKO_CTCF_ChIP_Seq_rep123_merged_bam_and_index.sh.runlog

```

# Detection of enrichment and repeatability

## plotFingerprint

```bash
# 参考： https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html
# 可以判断实验的成败、富集程度的好坏等
# mESC_SKO
plotFingerprint -b mESC_SKO_INPUT_rep1.bowtie2.sorted.bam mESC_SKO_INPUT_rep2.bowtie2.sorted.bam mESC_SKO_INPUT_rep3.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep1.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep2.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep3.bowtie2.sorted.bam --labels mESC_SKO_INPUT_rep1 mESC_SKO_INPUT_rep2 mESC_SKO_INPUT_rep3 mESC_SKO_IP_CTCF_rep1 mESC_SKO_IP_CTCF_rep2 mESC_SKO_IP_CTCF_rep3 -o plotFingerprint_for_mESC_SKO_CTCF_ChIP_Seq_samples.pdf --extendReads --ignoreDuplicates --centerReads --minMappingQuality 30 --numberOfSamples 1000000 --plotTitle "Fingerprints of different samples for mESC SKO" --skipZeros --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed -p "max"
# mESC_TKO
plotFingerprint -b mESC_TKO_INPUT_rep1.bowtie2.sorted.bam mESC_TKO_INPUT_rep2.bowtie2.sorted.bam mESC_TKO_INPUT_rep3.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep1.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep2.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep3.bowtie2.sorted.bam --labels mESC_TKO_INPUT_rep1 mESC_TKO_INPUT_rep2 mESC_TKO_INPUT_rep3 mESC_TKO_IP_CTCF_rep1 mESC_TKO_IP_CTCF_rep2 mESC_TKO_IP_CTCF_rep3 -o plotFingerprint_for_mESC_TKO_CTCF_ChIP_Seq_samples.pdf --extendReads --ignoreDuplicates --centerReads --minMappingQuality 30 --numberOfSamples 1000000 --plotTitle "Fingerprints of different samples for mESC TKO" --skipZeros --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed -p "max"
# mESC_WT
plotFingerprint -b mESC_WT_INPUT_rep1.bowtie2.sorted.bam mESC_WT_INPUT_rep2.bowtie2.sorted.bam mESC_WT_INPUT_rep3.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep1.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep2.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep3.bowtie2.sorted.bam --labels mESC_WT_INPUT_rep1 mESC_WT_INPUT_rep2 mESC_WT_INPUT_rep3 mESC_WT_IP_CTCF_rep1 mESC_WT_IP_CTCF_rep2 mESC_WT_IP_CTCF_rep3 -o plotFingerprint_for_mESC_WT_CTCF_ChIP_Seq_samples.pdf --extendReads --ignoreDuplicates --centerReads --minMappingQuality 30 --numberOfSamples 1000000 --plotTitle "Fingerprints of different samples for mESC WT" --skipZeros --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed -p "max"

# 图太多，使用linux的 pdftoppm 将 pdf 转 png 方便 vs code image gallery 预览缩略图
for i in `ls -1 plotFingerprint*.pdf`;do echo $i; pdftoppm -png $i ${i%.pdf};done

```

## plotCorrelation & PCA

```bash
# 使用 deeptools 分析 样本相似度
# 先对全基因组分bin统计 uniq reads数量
multiBamSummary bins --bamfiles mESC_SKO_INPUT_rep1.bowtie2.sorted.bam mESC_SKO_INPUT_rep2.bowtie2.sorted.bam mESC_SKO_INPUT_rep3.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep1.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep2.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep3.bowtie2.sorted.bam mESC_TKO_INPUT_rep1.bowtie2.sorted.bam mESC_TKO_INPUT_rep2.bowtie2.sorted.bam mESC_TKO_INPUT_rep3.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep1.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep2.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep3.bowtie2.sorted.bam mESC_WT_INPUT_rep1.bowtie2.sorted.bam mESC_WT_INPUT_rep2.bowtie2.sorted.bam mESC_WT_INPUT_rep3.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep1.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep2.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep3.bowtie2.sorted.bam -out multiBamSummary_results_for_all_bam_files.npz --labels mESC_SKO_INPUT_rep1 mESC_SKO_INPUT_rep2 mESC_SKO_INPUT_rep3 mESC_SKO_IP_CTCF_rep1 mESC_SKO_IP_CTCF_rep2 mESC_SKO_IP_CTCF_rep3 mESC_TKO_INPUT_rep1 mESC_TKO_INPUT_rep2 mESC_TKO_INPUT_rep3 mESC_TKO_IP_CTCF_rep1 mESC_TKO_IP_CTCF_rep2 mESC_TKO_IP_CTCF_rep3 mESC_WT_INPUT_rep1 mESC_WT_INPUT_rep2 mESC_WT_INPUT_rep3 mESC_WT_IP_CTCF_rep1 mESC_WT_IP_CTCF_rep2 mESC_WT_IP_CTCF_rep3 --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# 用 plotCorrelation 和 plotPCA 分析样本的相似度
plotCorrelation --corData multiBamSummary_results_for_all_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_all_samples.pdf --plotTitle "Pearson Correlation heatmap for all samples" --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_all_samples.txt --colorMap coolwarm --plotNumbers
plotCorrelation --corData multiBamSummary_results_for_all_bam_files.npz --corMethod pearson --whatToPlot scatterplot --plotFile plotCorrelation_scatterplot_plot_for_all_samples.pdf --plotTitle "Pearson Correlation scatterplot for all samples" --skipZeros --removeOutliers
plotPCA --corData multiBamSummary_results_for_all_bam_files.npz --plotFile plotPCA_plot_for_all_samples.pdf --plotTitle "plotPCA_plot_for_all_samples" --ntop 0 --log2 --outFileNameData plotPCA_stat_data_for_all_samples.txt
plotPCA --rowCenter --corData multiBamSummary_results_for_all_bam_files.npz --plotFile plotPCA_plot_for_all_samples.rowCenter.pdf --plotTitle "plotPCA_plot_for_all_samples" --ntop 0 --log2 --outFileNameData plotPCA_stat_data_for_all_samples.rowCenter.txt

# mESC_SKO 样本的相似度
multiBamSummary bins --bamfiles mESC_SKO_INPUT_rep1.bowtie2.sorted.bam mESC_SKO_INPUT_rep2.bowtie2.sorted.bam mESC_SKO_INPUT_rep3.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep1.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep2.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep3.bowtie2.sorted.bam -out multiBamSummary_results_for_mESC_SKO_bam_files.npz --labels mESC_SKO_INPUT_rep1 mESC_SKO_INPUT_rep2 mESC_SKO_INPUT_rep3 mESC_SKO_IP_CTCF_rep1 mESC_SKO_IP_CTCF_rep2 mESC_SKO_IP_CTCF_rep3 --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
plotCorrelation --corData multiBamSummary_results_for_mESC_SKO_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_mESC_SKO_samples.pdf --plotTitle "Pearson Correlation heatmap for mESC SKO samples" --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_mESC_SKO_samples.txt --colorMap coolwarm --plotNumbers
plotCorrelation --corData multiBamSummary_results_for_mESC_SKO_bam_files.npz --corMethod pearson --whatToPlot scatterplot --plotFile plotCorrelation_scatterplot_plot_for_mESC_SKO_samples.pdf --plotTitle "Pearson Correlation scatterplot for mESC SKO samples" --skipZeros --removeOutliers
plotPCA --corData multiBamSummary_results_for_mESC_SKO_bam_files.npz --plotFile plotPCA_plot_for_mESC_SKO_samples.pdf --plotTitle "plotPCA_plot_for_mESC_SKO_samples" --ntop 0 --log2 --outFileNameData plotPCA_stat_data_for_mESC_SKO_samples.txt
plotPCA --rowCenter --corData multiBamSummary_results_for_mESC_SKO_bam_files.npz --plotFile plotPCA_plot_for_mESC_SKO_samples.rowCenter.pdf --plotTitle "plotPCA_plot_for_mESC_SKO_samples" --ntop 0 --log2 --outFileNameData plotPCA_stat_data_for_mESC_SKO_samples.rowCenter.txt

# mESC_TKO 样本的相似度
multiBamSummary bins --bamfiles mESC_TKO_INPUT_rep1.bowtie2.sorted.bam mESC_TKO_INPUT_rep2.bowtie2.sorted.bam mESC_TKO_INPUT_rep3.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep1.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep2.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep3.bowtie2.sorted.bam -out multiBamSummary_results_for_mESC_TKO_bam_files.npz --labels mESC_TKO_INPUT_rep1 mESC_TKO_INPUT_rep2 mESC_TKO_INPUT_rep3 mESC_TKO_IP_CTCF_rep1 mESC_TKO_IP_CTCF_rep2 mESC_TKO_IP_CTCF_rep3 --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
plotCorrelation --corData multiBamSummary_results_for_mESC_TKO_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_mESC_TKO_samples.pdf --plotTitle "Pearson Correlation heatmap for mESC TKO samples" --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_mESC_TKO_samples.txt --colorMap coolwarm --plotNumbers
plotCorrelation --corData multiBamSummary_results_for_mESC_TKO_bam_files.npz --corMethod pearson --whatToPlot scatterplot --plotFile plotCorrelation_scatterplot_plot_for_mESC_TKO_samples.pdf --plotTitle "Pearson Correlation scatterplot for mESC TKO samples" --skipZeros --removeOutliers
plotPCA --corData multiBamSummary_results_for_mESC_TKO_bam_files.npz --plotFile plotPCA_plot_for_mESC_TKO_samples.pdf --plotTitle "plotPCA_plot_for_mESC_TKO_samples" --ntop 0 --log2 --outFileNameData plotPCA_stat_data_for_mESC_TKO_samples.txt
plotPCA --rowCenter --corData multiBamSummary_results_for_mESC_TKO_bam_files.npz --plotFile plotPCA_plot_for_mESC_TKO_samples.rowCenter.pdf --plotTitle "plotPCA_plot_for_mESC_TKO_samples" --ntop 0 --log2 --outFileNameData plotPCA_stat_data_for_mESC_TKO_samples.rowCenter.txt

# mESC_WT 样本的相似度
multiBamSummary bins --bamfiles mESC_WT_INPUT_rep1.bowtie2.sorted.bam mESC_WT_INPUT_rep2.bowtie2.sorted.bam mESC_WT_INPUT_rep3.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep1.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep2.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep3.bowtie2.sorted.bam -out multiBamSummary_results_for_mESC_WT_bam_files.npz --labels mESC_WT_INPUT_rep1 mESC_WT_INPUT_rep2 mESC_WT_INPUT_rep3 mESC_WT_IP_CTCF_rep1 mESC_WT_IP_CTCF_rep2 mESC_WT_IP_CTCF_rep3 --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
plotCorrelation --corData multiBamSummary_results_for_mESC_WT_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_mESC_WT_samples.pdf --plotTitle "Pearson Correlation heatmap for mESC WT samples" --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_mESC_WT_samples.txt --colorMap coolwarm --plotNumbers
plotCorrelation --corData multiBamSummary_results_for_mESC_WT_bam_files.npz --corMethod pearson --whatToPlot scatterplot --plotFile plotCorrelation_scatterplot_plot_for_mESC_WT_samples.pdf --plotTitle "Pearson Correlation scatterplot for mESC WT samples" --skipZeros --removeOutliers
plotPCA --corData multiBamSummary_results_for_mESC_WT_bam_files.npz --plotFile plotPCA_plot_for_mESC_WT_samples.pdf --plotTitle "plotPCA_plot_for_mESC_WT_samples" --ntop 0 --log2 --outFileNameData plotPCA_stat_data_for_mESC_WT_samples.txt
plotPCA --rowCenter --corData multiBamSummary_results_for_mESC_WT_bam_files.npz --plotFile plotPCA_plot_for_mESC_WT_samples.rowCenter.pdf --plotTitle "plotPCA_plot_for_mESC_WT_samples" --ntop 0 --log2 --outFileNameData plotPCA_stat_data_for_mESC_WT_samples.rowCenter.txt

# 图太多，使用linux的 pdftoppm 将 pdf 转 png 方便 vs code image gallery 预览缩略图
for i in `ls -1 plotCorrelation*.pdf`;do echo $i; pdftoppm -png $i ${i%.pdf};done
for i in `ls -1 plotPCA*.pdf`;do echo $i; pdftoppm -png $i ${i%.pdf};done

```

# macs2 peak calling

## peak calling for each repetition separately 

###  mESC WT CTCF ChIP-Seq peak calling

```bash
# data preperation
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
mkdir 6_macs2_peak_calling_1_for_mESC_WT_CTCF_each_repeat
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_1_for_mESC_WT_CTCF_each_repeat
# 创建所需的bam文件连接
# bam
ln -s ../mESC_WT_INPUT_rep1.bowtie2.sorted.bam mESC_WT_INPUT_rep1.bowtie2.sorted.bam
ln -s ../mESC_WT_INPUT_rep2.bowtie2.sorted.bam mESC_WT_INPUT_rep2.bowtie2.sorted.bam
ln -s ../mESC_WT_INPUT_rep3.bowtie2.sorted.bam mESC_WT_INPUT_rep3.bowtie2.sorted.bam
ln -s ../mESC_WT_IP_CTCF_rep1.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../mESC_WT_IP_CTCF_rep2.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../mESC_WT_IP_CTCF_rep3.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep3.bowtie2.sorted.bam
# bam.bai
ln -s ../mESC_WT_INPUT_rep1.bowtie2.sorted.bam.bai mESC_WT_INPUT_rep1.bowtie2.sorted.bam.bai
ln -s ../mESC_WT_INPUT_rep2.bowtie2.sorted.bam.bai mESC_WT_INPUT_rep2.bowtie2.sorted.bam.bai
ln -s ../mESC_WT_INPUT_rep3.bowtie2.sorted.bam.bai mESC_WT_INPUT_rep3.bowtie2.sorted.bam.bai
ln -s ../mESC_WT_IP_CTCF_rep1.bowtie2.sorted.bam.bai mESC_WT_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../mESC_WT_IP_CTCF_rep2.bowtie2.sorted.bam.bai mESC_WT_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../mESC_WT_IP_CTCF_rep3.bowtie2.sorted.bam.bai mESC_WT_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *.bam

# mESC_WT
nohup macs2 callpeak -t mESC_WT_IP_CTCF_rep1.bowtie2.sorted.bam -c mESC_WT_INPUT_rep1.bowtie2.sorted.bam --format BAMPE -g mm --keep-dup 1 --outdir mESC_WT_IP_CTCF_rep1_TvsC -n mESC_WT_IP_CTCF_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_mESC_WT_IP_CTCF_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t mESC_WT_IP_CTCF_rep2.bowtie2.sorted.bam -c mESC_WT_INPUT_rep2.bowtie2.sorted.bam --format BAMPE -g mm --keep-dup 1 --outdir mESC_WT_IP_CTCF_rep2_TvsC -n mESC_WT_IP_CTCF_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_mESC_WT_IP_CTCF_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t mESC_WT_IP_CTCF_rep3.bowtie2.sorted.bam -c mESC_WT_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g mm --keep-dup 1 --outdir mESC_WT_IP_CTCF_rep3_TvsC -n mESC_WT_IP_CTCF_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_mESC_WT_IP_CTCF_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_mESC_WT_IP_CTCF_rep1_TvsC.runlog
tail -f macs2_callpeak_for_mESC_WT_IP_CTCF_rep2_TvsC.runlog
tail -f macs2_callpeak_for_mESC_WT_IP_CTCF_rep3_TvsC.runlog

# 统计peak数量
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_1_for_mESC_WT_CTCF_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./mESC_WT_IP_CTCF_rep1_TvsC/mESC_WT_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# 48787
# ./mESC_WT_IP_CTCF_rep2_TvsC/mESC_WT_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# 48479
# ./mESC_WT_IP_CTCF_rep3_TvsC/mESC_WT_IP_CTCF_rep3_TvsC_peaks.narrowPeak
# 48181

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./mESC_WT_IP_CTCF_rep1_TvsC/mESC_WT_IP_CTCF_rep1_TvsC_summits.bed
# 49545
# ./mESC_WT_IP_CTCF_rep2_TvsC/mESC_WT_IP_CTCF_rep2_TvsC_summits.bed
# 49302
# ./mESC_WT_IP_CTCF_rep3_TvsC/mESC_WT_IP_CTCF_rep3_TvsC_summits.bed
# 49012

cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_1_for_mESC_WT_CTCF_each_repeat
# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_mESC_WT_IP_CTCF_ChIP_Seq_samples_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bigwig|sort`;do echo $i;done

```


#### rep1/2/3 peaks/summits venn plot

```bash
# 由于三个重复的peak数量相当，用韦恩图看一下peak的重复性如何
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_1_for_mESC_WT_CTCF_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i;done
# ./mESC_WT_IP_CTCF_rep1_TvsC/mESC_WT_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# ./mESC_WT_IP_CTCF_rep2_TvsC/mESC_WT_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# ./mESC_WT_IP_CTCF_rep3_TvsC/mESC_WT_IP_CTCF_rep3_TvsC_peaks.narrowPeak
for i in `find . -name **summits.bed|sort`;do echo $i;done
# ./mESC_WT_IP_CTCF_rep1_TvsC/mESC_WT_IP_CTCF_rep1_TvsC_summits.bed
# ./mESC_WT_IP_CTCF_rep2_TvsC/mESC_WT_IP_CTCF_rep2_TvsC_summits.bed
# ./mESC_WT_IP_CTCF_rep3_TvsC/mESC_WT_IP_CTCF_rep3_TvsC_summits.bed

for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.sorted.bed;done
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i%.bed}.sorted.bed ;done
mv `find . -name *sorted.bed` ./
# 构造 summits eB150 即 summits 两侧各延伸150nt
for i in `ls -1 mESC*TvsC_summits.sorted.bed|sort`;do echo $i; bedtools slop -b 150 -i $i -g /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes |sort -k1,1 -k2,2n |uniq >${i%.sorted.bed}.eB150.sorted.bed;done

# intervene plot
## mESC_WT CTCF rep1/2/3 peaks venn
ls -1 *narrowPeak.sorted.bed
# mESC_WT_IP_CTCF_rep1_TvsC_peaks.narrowPeak.sorted.bed
# mESC_WT_IP_CTCF_rep2_TvsC_peaks.narrowPeak.sorted.bed
# mESC_WT_IP_CTCF_rep3_TvsC_peaks.narrowPeak.sorted.bed
intervene venn -i mESC_WT_IP_CTCF_rep1_TvsC_peaks.narrowPeak.sorted.bed mESC_WT_IP_CTCF_rep2_TvsC_peaks.narrowPeak.sorted.bed mESC_WT_IP_CTCF_rep3_TvsC_peaks.narrowPeak.sorted.bed --names=mESC_WT_IP_CTCF_rep1_peaks,mESC_WT_IP_CTCF_rep2_peaks,mESC_WT_IP_CTCF_rep3_peaks -o 0_mESC_WT_IP_CTCF_rep123_peaks_venn --save-overlaps

## mESC_WT CTCF rep1/2/3 summits eB150 venn
ls -1 *summits.eB150.sorted.bed
# mESC_WT_IP_CTCF_rep1_TvsC_summits.eB150.sorted.bed
# mESC_WT_IP_CTCF_rep2_TvsC_summits.eB150.sorted.bed
# mESC_WT_IP_CTCF_rep3_TvsC_summits.eB150.sorted.bed
intervene venn -i mESC_WT_IP_CTCF_rep1_TvsC_summits.eB150.sorted.bed mESC_WT_IP_CTCF_rep2_TvsC_summits.eB150.sorted.bed mESC_WT_IP_CTCF_rep3_TvsC_summits.eB150.sorted.bed --names=mESC_WT_IP_CTCF_rep1_seB150,mESC_WT_IP_CTCF_rep2_seB150,mESC_WT_IP_CTCF_rep3_seB150 -o 0_mESC_WT_IP_CTCF_rep123_summits_eB150_venn --save-overlaps

# 面积权重韦恩图统一在项目根目录 0_make_venn_plot.ipynb 中绘制

```

###  mESC SKO CTCF ChIP-Seq peak calling

```bash
# data preperation
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
mkdir 6_macs2_peak_calling_1_for_mESC_SKO_CTCF_each_repeat
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_1_for_mESC_SKO_CTCF_each_repeat
# 创建所需的bam文件连接
# bam
ln -s ../mESC_SKO_INPUT_rep1.bowtie2.sorted.bam mESC_SKO_INPUT_rep1.bowtie2.sorted.bam
ln -s ../mESC_SKO_INPUT_rep2.bowtie2.sorted.bam mESC_SKO_INPUT_rep2.bowtie2.sorted.bam
ln -s ../mESC_SKO_INPUT_rep3.bowtie2.sorted.bam mESC_SKO_INPUT_rep3.bowtie2.sorted.bam
ln -s ../mESC_SKO_IP_CTCF_rep1.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../mESC_SKO_IP_CTCF_rep2.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../mESC_SKO_IP_CTCF_rep3.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep3.bowtie2.sorted.bam
# bam.bai
ln -s ../mESC_SKO_INPUT_rep1.bowtie2.sorted.bam.bai mESC_SKO_INPUT_rep1.bowtie2.sorted.bam.bai
ln -s ../mESC_SKO_INPUT_rep2.bowtie2.sorted.bam.bai mESC_SKO_INPUT_rep2.bowtie2.sorted.bam.bai
ln -s ../mESC_SKO_INPUT_rep3.bowtie2.sorted.bam.bai mESC_SKO_INPUT_rep3.bowtie2.sorted.bam.bai
ln -s ../mESC_SKO_IP_CTCF_rep1.bowtie2.sorted.bam.bai mESC_SKO_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../mESC_SKO_IP_CTCF_rep2.bowtie2.sorted.bam.bai mESC_SKO_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../mESC_SKO_IP_CTCF_rep3.bowtie2.sorted.bam.bai mESC_SKO_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *.bam

# mESC_SKO
nohup macs2 callpeak -t mESC_SKO_IP_CTCF_rep1.bowtie2.sorted.bam -c mESC_SKO_INPUT_rep1.bowtie2.sorted.bam --format BAMPE -g mm --keep-dup 1 --outdir mESC_SKO_IP_CTCF_rep1_TvsC -n mESC_SKO_IP_CTCF_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_mESC_SKO_IP_CTCF_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t mESC_SKO_IP_CTCF_rep2.bowtie2.sorted.bam -c mESC_SKO_INPUT_rep2.bowtie2.sorted.bam --format BAMPE -g mm --keep-dup 1 --outdir mESC_SKO_IP_CTCF_rep2_TvsC -n mESC_SKO_IP_CTCF_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_mESC_SKO_IP_CTCF_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t mESC_SKO_IP_CTCF_rep3.bowtie2.sorted.bam -c mESC_SKO_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g mm --keep-dup 1 --outdir mESC_SKO_IP_CTCF_rep3_TvsC -n mESC_SKO_IP_CTCF_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_mESC_SKO_IP_CTCF_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_mESC_SKO_IP_CTCF_rep1_TvsC.runlog
tail -f macs2_callpeak_for_mESC_SKO_IP_CTCF_rep2_TvsC.runlog
tail -f macs2_callpeak_for_mESC_SKO_IP_CTCF_rep3_TvsC.runlog

# 统计peak数量
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_1_for_mESC_SKO_CTCF_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./mESC_SKO_IP_CTCF_rep1_TvsC/mESC_SKO_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# 37909
# ./mESC_SKO_IP_CTCF_rep2_TvsC/mESC_SKO_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# 44296
# ./mESC_SKO_IP_CTCF_rep3_TvsC/mESC_SKO_IP_CTCF_rep3_TvsC_peaks.narrowPeak
# 43228

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./mESC_SKO_IP_CTCF_rep1_TvsC/mESC_SKO_IP_CTCF_rep1_TvsC_summits.bed
# 38448
# ./mESC_SKO_IP_CTCF_rep2_TvsC/mESC_SKO_IP_CTCF_rep2_TvsC_summits.bed
# 44935
# ./mESC_SKO_IP_CTCF_rep3_TvsC/mESC_SKO_IP_CTCF_rep3_TvsC_summits.bed
# 43875

cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_1_for_mESC_SKO_CTCF_each_repeat
# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_mESC_SKO_IP_CTCF_ChIP_Seq_samples_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bigwig|sort`;do echo $i;done

```


#### rep1/2/3 peaks/summits  venn plot

```bash
# 由于三个重复的peak数量相当，用韦恩图看一下peak的重复性如何
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_1_for_mESC_SKO_CTCF_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i;done
# ./mESC_SKO_IP_CTCF_rep1_TvsC/mESC_SKO_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# ./mESC_SKO_IP_CTCF_rep2_TvsC/mESC_SKO_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# ./mESC_SKO_IP_CTCF_rep3_TvsC/mESC_SKO_IP_CTCF_rep3_TvsC_peaks.narrowPeak
for i in `find . -name **summits.bed|sort`;do echo $i;done
# ./mESC_SKO_IP_CTCF_rep1_TvsC/mESC_SKO_IP_CTCF_rep1_TvsC_summits.bed
# ./mESC_SKO_IP_CTCF_rep2_TvsC/mESC_SKO_IP_CTCF_rep2_TvsC_summits.bed
# ./mESC_SKO_IP_CTCF_rep3_TvsC/mESC_SKO_IP_CTCF_rep3_TvsC_summits.bed

for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.sorted.bed;done
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i%.bed}.sorted.bed ;done
mv `find . -name *sorted.bed` ./
# 构造 summits eB150 即 summits 两侧各延伸150nt
for i in `ls -1 mESC*TvsC_summits.sorted.bed|sort`;do echo $i; bedtools slop -b 150 -i $i -g /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes |sort -k1,1 -k2,2n |uniq >${i%.sorted.bed}.eB150.sorted.bed;done

# intervene plot
## mESC_SKO CTCF rep1/2/3 peaks venn
ls -1 *narrowPeak.sorted.bed
# mESC_SKO_IP_CTCF_rep1_TvsC_peaks.narrowPeak.sorted.bed
# mESC_SKO_IP_CTCF_rep2_TvsC_peaks.narrowPeak.sorted.bed
# mESC_SKO_IP_CTCF_rep3_TvsC_peaks.narrowPeak.sorted.bed
intervene venn -i mESC_SKO_IP_CTCF_rep1_TvsC_peaks.narrowPeak.sorted.bed mESC_SKO_IP_CTCF_rep2_TvsC_peaks.narrowPeak.sorted.bed mESC_SKO_IP_CTCF_rep3_TvsC_peaks.narrowPeak.sorted.bed --names=mESC_SKO_IP_CTCF_rep1_peaks,mESC_SKO_IP_CTCF_rep2_peaks,mESC_SKO_IP_CTCF_rep3_peaks -o 0_mESC_SKO_IP_CTCF_rep123_peaks_venn --save-overlaps

## mESC_SKO CTCF rep1/2/3 summits eB150 venn
ls -1 *summits.eB150.sorted.bed
# mESC_SKO_IP_CTCF_rep1_TvsC_summits.eB150.sorted.bed
# mESC_SKO_IP_CTCF_rep2_TvsC_summits.eB150.sorted.bed
# mESC_SKO_IP_CTCF_rep3_TvsC_summits.eB150.sorted.bed
intervene venn -i mESC_SKO_IP_CTCF_rep1_TvsC_summits.eB150.sorted.bed mESC_SKO_IP_CTCF_rep2_TvsC_summits.eB150.sorted.bed mESC_SKO_IP_CTCF_rep3_TvsC_summits.eB150.sorted.bed --names=mESC_SKO_IP_CTCF_rep1_seB150,mESC_SKO_IP_CTCF_rep2_seB150,mESC_SKO_IP_CTCF_rep3_seB150 -o 0_mESC_SKO_IP_CTCF_rep123_summits_eB150_venn --save-overlaps

# 面积权重韦恩图统一在项目根目录 0_make_venn_plot.ipynb 中绘制

```



###  mESC TKO CTCF ChIP-Seq peak calling

```bash
# data preperation
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
mkdir 6_macs2_peak_calling_1_for_mESC_TKO_CTCF_each_repeat
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_1_for_mESC_TKO_CTCF_each_repeat
# 创建所需的bam文件连接
# bam
ln -s ../mESC_TKO_INPUT_rep1.bowtie2.sorted.bam mESC_TKO_INPUT_rep1.bowtie2.sorted.bam
ln -s ../mESC_TKO_INPUT_rep2.bowtie2.sorted.bam mESC_TKO_INPUT_rep2.bowtie2.sorted.bam
ln -s ../mESC_TKO_INPUT_rep3.bowtie2.sorted.bam mESC_TKO_INPUT_rep3.bowtie2.sorted.bam
ln -s ../mESC_TKO_IP_CTCF_rep1.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../mESC_TKO_IP_CTCF_rep2.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../mESC_TKO_IP_CTCF_rep3.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep3.bowtie2.sorted.bam
# bam.bai
ln -s ../mESC_TKO_INPUT_rep1.bowtie2.sorted.bam.bai mESC_TKO_INPUT_rep1.bowtie2.sorted.bam.bai
ln -s ../mESC_TKO_INPUT_rep2.bowtie2.sorted.bam.bai mESC_TKO_INPUT_rep2.bowtie2.sorted.bam.bai
ln -s ../mESC_TKO_INPUT_rep3.bowtie2.sorted.bam.bai mESC_TKO_INPUT_rep3.bowtie2.sorted.bam.bai
ln -s ../mESC_TKO_IP_CTCF_rep1.bowtie2.sorted.bam.bai mESC_TKO_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../mESC_TKO_IP_CTCF_rep2.bowtie2.sorted.bam.bai mESC_TKO_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../mESC_TKO_IP_CTCF_rep3.bowtie2.sorted.bam.bai mESC_TKO_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *.bam

# mESC_TKO
nohup macs2 callpeak -t mESC_TKO_IP_CTCF_rep1.bowtie2.sorted.bam -c mESC_TKO_INPUT_rep1.bowtie2.sorted.bam --format BAMPE -g mm --keep-dup 1 --outdir mESC_TKO_IP_CTCF_rep1_TvsC -n mESC_TKO_IP_CTCF_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_mESC_TKO_IP_CTCF_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t mESC_TKO_IP_CTCF_rep2.bowtie2.sorted.bam -c mESC_TKO_INPUT_rep2.bowtie2.sorted.bam --format BAMPE -g mm --keep-dup 1 --outdir mESC_TKO_IP_CTCF_rep2_TvsC -n mESC_TKO_IP_CTCF_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_mESC_TKO_IP_CTCF_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t mESC_TKO_IP_CTCF_rep3.bowtie2.sorted.bam -c mESC_TKO_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g mm --keep-dup 1 --outdir mESC_TKO_IP_CTCF_rep3_TvsC -n mESC_TKO_IP_CTCF_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_mESC_TKO_IP_CTCF_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_mESC_TKO_IP_CTCF_rep1_TvsC.runlog
tail -f macs2_callpeak_for_mESC_TKO_IP_CTCF_rep2_TvsC.runlog
tail -f macs2_callpeak_for_mESC_TKO_IP_CTCF_rep3_TvsC.runlog

# 统计peak数量
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_1_for_mESC_TKO_CTCF_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./mESC_TKO_IP_CTCF_rep1_TvsC/mESC_TKO_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# 41852
# ./mESC_TKO_IP_CTCF_rep2_TvsC/mESC_TKO_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# 42347
# ./mESC_TKO_IP_CTCF_rep3_TvsC/mESC_TKO_IP_CTCF_rep3_TvsC_peaks.narrowPeak
# 41534

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./mESC_TKO_IP_CTCF_rep1_TvsC/mESC_TKO_IP_CTCF_rep1_TvsC_summits.bed
# 42588
# ./mESC_TKO_IP_CTCF_rep2_TvsC/mESC_TKO_IP_CTCF_rep2_TvsC_summits.bed
# 43074
# ./mESC_TKO_IP_CTCF_rep3_TvsC/mESC_TKO_IP_CTCF_rep3_TvsC_summits.bed
# 42246

cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_1_for_mESC_TKO_CTCF_each_repeat
# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_mESC_TKO_IP_CTCF_ChIP_Seq_samples_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bigwig|sort`;do echo $i;done

```


#### rep1/2/3 peaks/summits  venn plot

```bash
# 由于三个重复的peak数量相当，用韦恩图看一下peak的重复性如何
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_1_for_mESC_TKO_CTCF_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i;done
# ./mESC_TKO_IP_CTCF_rep1_TvsC/mESC_TKO_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# ./mESC_TKO_IP_CTCF_rep2_TvsC/mESC_TKO_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# ./mESC_TKO_IP_CTCF_rep3_TvsC/mESC_TKO_IP_CTCF_rep3_TvsC_peaks.narrowPeak
for i in `find . -name **summits.bed|sort`;do echo $i;done
# ./mESC_TKO_IP_CTCF_rep1_TvsC/mESC_TKO_IP_CTCF_rep1_TvsC_summits.bed
# ./mESC_TKO_IP_CTCF_rep2_TvsC/mESC_TKO_IP_CTCF_rep2_TvsC_summits.bed
# ./mESC_TKO_IP_CTCF_rep3_TvsC/mESC_TKO_IP_CTCF_rep3_TvsC_summits.bed

for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.sorted.bed;done
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i%.bed}.sorted.bed ;done
mv `find . -name *sorted.bed` ./
# 构造 summits eB150 即 summits 两侧各延伸150nt
for i in `ls -1 mESC*TvsC_summits.sorted.bed|sort`;do echo $i; bedtools slop -b 150 -i $i -g /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes |sort -k1,1 -k2,2n |uniq >${i%.sorted.bed}.eB150.sorted.bed;done

# intervene plot
## mESC_TKO CTCF rep1/2/3 peaks venn
ls -1 *narrowPeak.sorted.bed
# mESC_TKO_IP_CTCF_rep1_TvsC_peaks.narrowPeak.sorted.bed
# mESC_TKO_IP_CTCF_rep2_TvsC_peaks.narrowPeak.sorted.bed
# mESC_TKO_IP_CTCF_rep3_TvsC_peaks.narrowPeak.sorted.bed
intervene venn -i mESC_TKO_IP_CTCF_rep1_TvsC_peaks.narrowPeak.sorted.bed mESC_TKO_IP_CTCF_rep2_TvsC_peaks.narrowPeak.sorted.bed mESC_TKO_IP_CTCF_rep3_TvsC_peaks.narrowPeak.sorted.bed --names=mESC_TKO_IP_CTCF_rep1_peaks,mESC_TKO_IP_CTCF_rep2_peaks,mESC_TKO_IP_CTCF_rep3_peaks -o 0_mESC_TKO_IP_CTCF_rep123_peaks_venn --save-overlaps

## mESC_TKO CTCF rep1/2/3 summits eB150 venn
ls -1 *summits.eB150.sorted.bed
# mESC_TKO_IP_CTCF_rep1_TvsC_summits.eB150.sorted.bed
# mESC_TKO_IP_CTCF_rep2_TvsC_summits.eB150.sorted.bed
# mESC_TKO_IP_CTCF_rep3_TvsC_summits.eB150.sorted.bed
intervene venn -i mESC_TKO_IP_CTCF_rep1_TvsC_summits.eB150.sorted.bed mESC_TKO_IP_CTCF_rep2_TvsC_summits.eB150.sorted.bed mESC_TKO_IP_CTCF_rep3_TvsC_summits.eB150.sorted.bed --names=mESC_TKO_IP_CTCF_rep1_seB150,mESC_TKO_IP_CTCF_rep2_seB150,mESC_TKO_IP_CTCF_rep3_seB150 -o 0_mESC_TKO_IP_CTCF_rep123_summits_eB150_venn --save-overlaps

# 面积权重韦恩图统一在项目根目录 0_make_venn_plot.ipynb 中绘制

```



## peak calling for rep123 merged

### mESC WT CTCF  rep123 merged peak calling
```bash
# data preperation
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
mkdir 6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together
# 创建所需的bam文件连接
# bam
ln -s ../mESC_WT_INPUT_rep1.bowtie2.sorted.bam mESC_WT_INPUT_rep1.bowtie2.sorted.bam
ln -s ../mESC_WT_INPUT_rep2.bowtie2.sorted.bam mESC_WT_INPUT_rep2.bowtie2.sorted.bam
ln -s ../mESC_WT_INPUT_rep3.bowtie2.sorted.bam mESC_WT_INPUT_rep3.bowtie2.sorted.bam
ln -s ../mESC_WT_IP_CTCF_rep1.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../mESC_WT_IP_CTCF_rep2.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../mESC_WT_IP_CTCF_rep3.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep3.bowtie2.sorted.bam
# bam.bai
ln -s ../mESC_WT_INPUT_rep1.bowtie2.sorted.bam.bai mESC_WT_INPUT_rep1.bowtie2.sorted.bam.bai
ln -s ../mESC_WT_INPUT_rep2.bowtie2.sorted.bam.bai mESC_WT_INPUT_rep2.bowtie2.sorted.bam.bai
ln -s ../mESC_WT_INPUT_rep3.bowtie2.sorted.bam.bai mESC_WT_INPUT_rep3.bowtie2.sorted.bam.bai
ln -s ../mESC_WT_IP_CTCF_rep1.bowtie2.sorted.bam.bai mESC_WT_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../mESC_WT_IP_CTCF_rep2.bowtie2.sorted.bam.bai mESC_WT_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../mESC_WT_IP_CTCF_rep3.bowtie2.sorted.bam.bai mESC_WT_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *.bam

# mESC_WT_IP_CTCF_rep123_TvsC
nohup macs2 callpeak -t mESC_WT_IP_CTCF_rep1.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep2.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep3.bowtie2.sorted.bam -c mESC_WT_INPUT_rep1.bowtie2.sorted.bam mESC_WT_INPUT_rep2.bowtie2.sorted.bam mESC_WT_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g mm --keep-dup 1 --outdir mESC_WT_IP_CTCF_rep123_TvsC -n mESC_WT_IP_CTCF_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_mESC_WT_IP_CTCF_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_mESC_WT_IP_CTCF_rep123_TvsC.runlog

# 统计peak数量
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./mESC_WT_IP_CTCF_rep123_TvsC/mESC_WT_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# 71075

# 统计 peak summit 数量
for i in `find . -name *summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./mESC_WT_IP_CTCF_rep123_TvsC/mESC_WT_IP_CTCF_rep123_TvsC_summits.bed
# 73125

cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together
# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_mESC_WT_IP_CTCF_ChIP_Seq_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes ${i%.bdg}.bigwig;done

```
```bash
# 进行 CvsT 比较，根据结果数量，排查是否存在标签标反的情况
# data preperation
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
mkdir 6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together_CvsT
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together_CvsT
# 创建所需的bam文件连接
# bam
ln -s ../mESC_WT_INPUT_rep1.bowtie2.sorted.bam mESC_WT_INPUT_rep1.bowtie2.sorted.bam
ln -s ../mESC_WT_INPUT_rep2.bowtie2.sorted.bam mESC_WT_INPUT_rep2.bowtie2.sorted.bam
ln -s ../mESC_WT_INPUT_rep3.bowtie2.sorted.bam mESC_WT_INPUT_rep3.bowtie2.sorted.bam
ln -s ../mESC_WT_IP_CTCF_rep1.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../mESC_WT_IP_CTCF_rep2.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../mESC_WT_IP_CTCF_rep3.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep3.bowtie2.sorted.bam
# bam.bai
ln -s ../mESC_WT_INPUT_rep1.bowtie2.sorted.bam.bai mESC_WT_INPUT_rep1.bowtie2.sorted.bam.bai
ln -s ../mESC_WT_INPUT_rep2.bowtie2.sorted.bam.bai mESC_WT_INPUT_rep2.bowtie2.sorted.bam.bai
ln -s ../mESC_WT_INPUT_rep3.bowtie2.sorted.bam.bai mESC_WT_INPUT_rep3.bowtie2.sorted.bam.bai
ln -s ../mESC_WT_IP_CTCF_rep1.bowtie2.sorted.bam.bai mESC_WT_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../mESC_WT_IP_CTCF_rep2.bowtie2.sorted.bam.bai mESC_WT_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../mESC_WT_IP_CTCF_rep3.bowtie2.sorted.bam.bai mESC_WT_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *.bam

# mESC_WT_IP_CTCF_rep123_CvsT
nohup macs2 callpeak -c mESC_WT_IP_CTCF_rep1.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep2.bowtie2.sorted.bam mESC_WT_IP_CTCF_rep3.bowtie2.sorted.bam -t mESC_WT_INPUT_rep1.bowtie2.sorted.bam mESC_WT_INPUT_rep2.bowtie2.sorted.bam mESC_WT_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g mm --keep-dup 1 --outdir mESC_WT_IP_CTCF_rep123_CvsT -n mESC_WT_IP_CTCF_rep123_CvsT -B --call-summits --cutoff-analysis >macs2_callpeak_for_mESC_WT_IP_CTCF_rep123_CvsT.runlog 2>&1 &
tail -f macs2_callpeak_for_mESC_WT_IP_CTCF_rep123_CvsT.runlog

# 统计peak数量
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together_CvsT
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./mESC_WT_IP_CTCF_rep123_CvsT/mESC_WT_IP_CTCF_rep123_CvsT_peaks.narrowPeak
# 1777

# 统计 peak summit 数量
for i in `find . -name *summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./mESC_WT_IP_CTCF_rep123_CvsT/mESC_WT_IP_CTCF_rep123_CvsT_summits.bed
# 2033

# 删除 bdg 文件
du -hs
for i in `find . -name *CvsT*.bdg|sort`;do echo $i;done
for i in `find . -name *CvsT*.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *CvsT*.bdg|sort`;do echo $i;done
du -hs

```
### mESC SKO CTCF  rep123 merged peak calling
```bash
# data preperation
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
mkdir 6_macs2_peak_calling_2_for_mESC_SKO_CTCF_rep123_together
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_SKO_CTCF_rep123_together
# 创建所需的bam文件连接
# bam
ln -s ../mESC_SKO_INPUT_rep1.bowtie2.sorted.bam mESC_SKO_INPUT_rep1.bowtie2.sorted.bam
ln -s ../mESC_SKO_INPUT_rep2.bowtie2.sorted.bam mESC_SKO_INPUT_rep2.bowtie2.sorted.bam
ln -s ../mESC_SKO_INPUT_rep3.bowtie2.sorted.bam mESC_SKO_INPUT_rep3.bowtie2.sorted.bam
ln -s ../mESC_SKO_IP_CTCF_rep1.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../mESC_SKO_IP_CTCF_rep2.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../mESC_SKO_IP_CTCF_rep3.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep3.bowtie2.sorted.bam
# bam.bai
ln -s ../mESC_SKO_INPUT_rep1.bowtie2.sorted.bam.bai mESC_SKO_INPUT_rep1.bowtie2.sorted.bam.bai
ln -s ../mESC_SKO_INPUT_rep2.bowtie2.sorted.bam.bai mESC_SKO_INPUT_rep2.bowtie2.sorted.bam.bai
ln -s ../mESC_SKO_INPUT_rep3.bowtie2.sorted.bam.bai mESC_SKO_INPUT_rep3.bowtie2.sorted.bam.bai
ln -s ../mESC_SKO_IP_CTCF_rep1.bowtie2.sorted.bam.bai mESC_SKO_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../mESC_SKO_IP_CTCF_rep2.bowtie2.sorted.bam.bai mESC_SKO_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../mESC_SKO_IP_CTCF_rep3.bowtie2.sorted.bam.bai mESC_SKO_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *.bam

# mESC_SKO_IP_CTCF_rep123_TvsC
nohup macs2 callpeak -t mESC_SKO_IP_CTCF_rep1.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep2.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep3.bowtie2.sorted.bam -c mESC_SKO_INPUT_rep1.bowtie2.sorted.bam mESC_SKO_INPUT_rep2.bowtie2.sorted.bam mESC_SKO_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g mm --keep-dup 1 --outdir mESC_SKO_IP_CTCF_rep123_TvsC -n mESC_SKO_IP_CTCF_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_mESC_SKO_IP_CTCF_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_mESC_SKO_IP_CTCF_rep123_TvsC.runlog

# 统计peak数量
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_SKO_CTCF_rep123_together
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./mESC_SKO_IP_CTCF_rep123_TvsC/mESC_SKO_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# 59372

# 统计 peak summit 数量
for i in `find . -name *summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./mESC_SKO_IP_CTCF_rep123_TvsC/mESC_SKO_IP_CTCF_rep123_TvsC_summits.bed
# 61027

cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_SKO_CTCF_rep123_together
# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_mESC_SKO_IP_CTCF_ChIP_Seq_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_SKO_CTCF_rep123_together
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes ${i%.bdg}.bigwig;done

```
```bash
# 进行 CvsT 比较，根据结果数量，排查是否存在标签标反的情况
# data preperation
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
mkdir 6_macs2_peak_calling_2_for_mESC_SKO_CTCF_rep123_together_CvsT
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_SKO_CTCF_rep123_together_CvsT
# 创建所需的bam文件连接
# bam
ln -s ../mESC_SKO_INPUT_rep1.bowtie2.sorted.bam mESC_SKO_INPUT_rep1.bowtie2.sorted.bam
ln -s ../mESC_SKO_INPUT_rep2.bowtie2.sorted.bam mESC_SKO_INPUT_rep2.bowtie2.sorted.bam
ln -s ../mESC_SKO_INPUT_rep3.bowtie2.sorted.bam mESC_SKO_INPUT_rep3.bowtie2.sorted.bam
ln -s ../mESC_SKO_IP_CTCF_rep1.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../mESC_SKO_IP_CTCF_rep2.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../mESC_SKO_IP_CTCF_rep3.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep3.bowtie2.sorted.bam
# bam.bai
ln -s ../mESC_SKO_INPUT_rep1.bowtie2.sorted.bam.bai mESC_SKO_INPUT_rep1.bowtie2.sorted.bam.bai
ln -s ../mESC_SKO_INPUT_rep2.bowtie2.sorted.bam.bai mESC_SKO_INPUT_rep2.bowtie2.sorted.bam.bai
ln -s ../mESC_SKO_INPUT_rep3.bowtie2.sorted.bam.bai mESC_SKO_INPUT_rep3.bowtie2.sorted.bam.bai
ln -s ../mESC_SKO_IP_CTCF_rep1.bowtie2.sorted.bam.bai mESC_SKO_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../mESC_SKO_IP_CTCF_rep2.bowtie2.sorted.bam.bai mESC_SKO_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../mESC_SKO_IP_CTCF_rep3.bowtie2.sorted.bam.bai mESC_SKO_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *.bam

# mESC_SKO_IP_CTCF_rep123_CvsT
nohup macs2 callpeak -c mESC_SKO_IP_CTCF_rep1.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep2.bowtie2.sorted.bam mESC_SKO_IP_CTCF_rep3.bowtie2.sorted.bam -t mESC_SKO_INPUT_rep1.bowtie2.sorted.bam mESC_SKO_INPUT_rep2.bowtie2.sorted.bam mESC_SKO_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g mm --keep-dup 1 --outdir mESC_SKO_IP_CTCF_rep123_CvsT -n mESC_SKO_IP_CTCF_rep123_CvsT -B --call-summits --cutoff-analysis >macs2_callpeak_for_mESC_SKO_IP_CTCF_rep123_CvsT.runlog 2>&1 &
tail -f macs2_callpeak_for_mESC_SKO_IP_CTCF_rep123_CvsT.runlog

# 统计peak数量
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_SKO_CTCF_rep123_together_CvsT
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./mESC_SKO_IP_CTCF_rep123_CvsT/mESC_SKO_IP_CTCF_rep123_CvsT_peaks.narrowPeak
# 955

# 统计 peak summit 数量
for i in `find . -name *summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./mESC_SKO_IP_CTCF_rep123_CvsT/mESC_SKO_IP_CTCF_rep123_CvsT_summits.bed
# 1149

# 删除 bdg 文件
du -hs
for i in `find . -name *CvsT*.bdg|sort`;do echo $i;done
for i in `find . -name *CvsT*.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *CvsT*.bdg|sort`;do echo $i;done
du -hs

```
### mESC TKO CTCF  rep123 merged peak calling
```bash
# data preperation
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
mkdir 6_macs2_peak_calling_2_for_mESC_TKO_CTCF_rep123_together
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_TKO_CTCF_rep123_together
# 创建所需的bam文件连接
# bam
ln -s ../mESC_TKO_INPUT_rep1.bowtie2.sorted.bam mESC_TKO_INPUT_rep1.bowtie2.sorted.bam
ln -s ../mESC_TKO_INPUT_rep2.bowtie2.sorted.bam mESC_TKO_INPUT_rep2.bowtie2.sorted.bam
ln -s ../mESC_TKO_INPUT_rep3.bowtie2.sorted.bam mESC_TKO_INPUT_rep3.bowtie2.sorted.bam
ln -s ../mESC_TKO_IP_CTCF_rep1.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../mESC_TKO_IP_CTCF_rep2.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../mESC_TKO_IP_CTCF_rep3.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep3.bowtie2.sorted.bam
# bam.bai
ln -s ../mESC_TKO_INPUT_rep1.bowtie2.sorted.bam.bai mESC_TKO_INPUT_rep1.bowtie2.sorted.bam.bai
ln -s ../mESC_TKO_INPUT_rep2.bowtie2.sorted.bam.bai mESC_TKO_INPUT_rep2.bowtie2.sorted.bam.bai
ln -s ../mESC_TKO_INPUT_rep3.bowtie2.sorted.bam.bai mESC_TKO_INPUT_rep3.bowtie2.sorted.bam.bai
ln -s ../mESC_TKO_IP_CTCF_rep1.bowtie2.sorted.bam.bai mESC_TKO_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../mESC_TKO_IP_CTCF_rep2.bowtie2.sorted.bam.bai mESC_TKO_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../mESC_TKO_IP_CTCF_rep3.bowtie2.sorted.bam.bai mESC_TKO_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *.bam

# mESC_TKO_IP_CTCF_rep123_TvsC
nohup macs2 callpeak -t mESC_TKO_IP_CTCF_rep1.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep2.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep3.bowtie2.sorted.bam -c mESC_TKO_INPUT_rep1.bowtie2.sorted.bam mESC_TKO_INPUT_rep2.bowtie2.sorted.bam mESC_TKO_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g mm --keep-dup 1 --outdir mESC_TKO_IP_CTCF_rep123_TvsC -n mESC_TKO_IP_CTCF_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_mESC_TKO_IP_CTCF_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_mESC_TKO_IP_CTCF_rep123_TvsC.runlog

# 统计peak数量
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_TKO_CTCF_rep123_together
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./mESC_TKO_IP_CTCF_rep123_TvsC/mESC_TKO_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# 62238

# 统计 peak summit 数量
for i in `find . -name *summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./mESC_TKO_IP_CTCF_rep123_TvsC/mESC_TKO_IP_CTCF_rep123_TvsC_summits.bed
# 64121

cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_TKO_CTCF_rep123_together
# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_mESC_TKO_IP_CTCF_ChIP_Seq_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_TKO_CTCF_rep123_together
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes ${i%.bdg}.bigwig;done

```
```bash
# 进行 CvsT 比较，根据结果数量，排查是否存在标签标反的情况
# data preperation
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
mkdir 6_macs2_peak_calling_2_for_mESC_TKO_CTCF_rep123_together_CvsT
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_TKO_CTCF_rep123_together_CvsT
# 创建所需的bam文件连接
# bam
ln -s ../mESC_TKO_INPUT_rep1.bowtie2.sorted.bam mESC_TKO_INPUT_rep1.bowtie2.sorted.bam
ln -s ../mESC_TKO_INPUT_rep2.bowtie2.sorted.bam mESC_TKO_INPUT_rep2.bowtie2.sorted.bam
ln -s ../mESC_TKO_INPUT_rep3.bowtie2.sorted.bam mESC_TKO_INPUT_rep3.bowtie2.sorted.bam
ln -s ../mESC_TKO_IP_CTCF_rep1.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../mESC_TKO_IP_CTCF_rep2.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../mESC_TKO_IP_CTCF_rep3.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep3.bowtie2.sorted.bam
# bam.bai
ln -s ../mESC_TKO_INPUT_rep1.bowtie2.sorted.bam.bai mESC_TKO_INPUT_rep1.bowtie2.sorted.bam.bai
ln -s ../mESC_TKO_INPUT_rep2.bowtie2.sorted.bam.bai mESC_TKO_INPUT_rep2.bowtie2.sorted.bam.bai
ln -s ../mESC_TKO_INPUT_rep3.bowtie2.sorted.bam.bai mESC_TKO_INPUT_rep3.bowtie2.sorted.bam.bai
ln -s ../mESC_TKO_IP_CTCF_rep1.bowtie2.sorted.bam.bai mESC_TKO_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../mESC_TKO_IP_CTCF_rep2.bowtie2.sorted.bam.bai mESC_TKO_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../mESC_TKO_IP_CTCF_rep3.bowtie2.sorted.bam.bai mESC_TKO_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *.bam

# mESC_TKO_IP_CTCF_rep123_CvsT
nohup macs2 callpeak -c mESC_TKO_IP_CTCF_rep1.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep2.bowtie2.sorted.bam mESC_TKO_IP_CTCF_rep3.bowtie2.sorted.bam -t mESC_TKO_INPUT_rep1.bowtie2.sorted.bam mESC_TKO_INPUT_rep2.bowtie2.sorted.bam mESC_TKO_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g mm --keep-dup 1 --outdir mESC_TKO_IP_CTCF_rep123_CvsT -n mESC_TKO_IP_CTCF_rep123_CvsT -B --call-summits --cutoff-analysis >macs2_callpeak_for_mESC_TKO_IP_CTCF_rep123_CvsT.runlog 2>&1 &
tail -f macs2_callpeak_for_mESC_TKO_IP_CTCF_rep123_CvsT.runlog

# 统计peak数量
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_TKO_CTCF_rep123_together_CvsT
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./mESC_TKO_IP_CTCF_rep123_CvsT/mESC_TKO_IP_CTCF_rep123_CvsT_peaks.narrowPeak
# 1070

# 统计 peak summit 数量
for i in `find . -name *summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./mESC_TKO_IP_CTCF_rep123_CvsT/mESC_TKO_IP_CTCF_rep123_CvsT_summits.bed
# 1252

# 删除 bdg 文件
du -hs
for i in `find . -name *CvsT*.bdg|sort`;do echo $i;done
for i in `find . -name *CvsT*.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *CvsT*.bdg|sort`;do echo $i;done
du -hs

```

### WT & SKO & TKO peaks / summits eB150 intervene

```bash
# SKO raw narrowPeak 和 summit 改造为 sorted bed3
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_SKO_CTCF_rep123_together/mESC_SKO_IP_CTCF_rep123_TvsC/
## raw narrowPeaks to sorted bed3
cut -f 1,2,3 mESC_SKO_IP_CTCF_rep123_TvsC_peaks.narrowPeak |head
cut -f 1,2,3 mESC_SKO_IP_CTCF_rep123_TvsC_peaks.narrowPeak |sort -k1,1 -k2,2n |uniq >mESC_SKO_IP_CTCF.narrowPeaks.sorted.bed
## peak summits to summits eB150 sorted beds
cut -f 1,2,3 mESC_SKO_IP_CTCF_rep123_TvsC_summits.bed |head
cut -f 1,2,3 mESC_SKO_IP_CTCF_rep123_TvsC_summits.bed |sort -k1,1 -k2,2n |uniq |bedtools slop -b 150 -i stdin -g /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes |sort -k1,1 -k2,2n |uniq >mESC_SKO_IP_CTCF.peakSummits.eB150.sorted.bed

# TKO raw narrowPeak 和 summit 改造为 sorted bed3
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_TKO_CTCF_rep123_together/mESC_TKO_IP_CTCF_rep123_TvsC/
## raw narrowPeaks to sorted bed3
cut -f 1,2,3 mESC_TKO_IP_CTCF_rep123_TvsC_peaks.narrowPeak |head
cut -f 1,2,3 mESC_TKO_IP_CTCF_rep123_TvsC_peaks.narrowPeak |sort -k1,1 -k2,2n |uniq >mESC_TKO_IP_CTCF.narrowPeaks.sorted.bed
## peak summits to summits eB150 sorted beds
cut -f 1,2,3 mESC_TKO_IP_CTCF_rep123_TvsC_summits.bed |head
cut -f 1,2,3 mESC_TKO_IP_CTCF_rep123_TvsC_summits.bed |sort -k1,1 -k2,2n |uniq |bedtools slop -b 150 -i stdin -g /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes |sort -k1,1 -k2,2n |uniq >mESC_TKO_IP_CTCF.peakSummits.eB150.sorted.bed

# WT raw narrowPeak 和 summit 改造为 sorted bed3
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/mESC_WT_IP_CTCF_rep123_TvsC/
## raw narrowPeaks to sorted bed3
cut -f 1,2,3 mESC_WT_IP_CTCF_rep123_TvsC_peaks.narrowPeak |head
cut -f 1,2,3 mESC_WT_IP_CTCF_rep123_TvsC_peaks.narrowPeak |sort -k1,1 -k2,2n |uniq >mESC_WT_IP_CTCF.narrowPeaks.sorted.bed
## peak summits to summits eB150 sorted beds
cut -f 1,2,3 mESC_WT_IP_CTCF_rep123_TvsC_summits.bed |head
cut -f 1,2,3 mESC_WT_IP_CTCF_rep123_TvsC_summits.bed |sort -k1,1 -k2,2n |uniq |bedtools slop -b 150 -i stdin -g /home/sunwenju/3_genomeData/mouse/mm10/mm10.chrom.sizes |sort -k1,1 -k2,2n |uniq >mESC_WT_IP_CTCF.peakSummits.eB150.sorted.bed

cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together
# 创建 WT peak and summit eB150 bed file links
ln -s ./mESC_WT_IP_CTCF_rep123_TvsC/mESC_WT_IP_CTCF.narrowPeaks.sorted.bed ./mESC_WT_IP_CTCF.narrowPeaks.sorted.bed
ln -s ./mESC_WT_IP_CTCF_rep123_TvsC/mESC_WT_IP_CTCF.peakSummits.eB150.sorted.bed ./mESC_WT_IP_CTCF.peakSummits.eB150.sorted.bed
# 创建 SKO peak and summit eB150 bed file links
ln -s ../6_macs2_peak_calling_2_for_mESC_SKO_CTCF_rep123_together/mESC_SKO_IP_CTCF_rep123_TvsC/mESC_SKO_IP_CTCF.narrowPeaks.sorted.bed ./mESC_SKO_IP_CTCF.narrowPeaks.sorted.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_SKO_CTCF_rep123_together/mESC_SKO_IP_CTCF_rep123_TvsC/mESC_SKO_IP_CTCF.peakSummits.eB150.sorted.bed ./mESC_SKO_IP_CTCF.peakSummits.eB150.sorted.bed
# 创建 TKO peak and summit eB150 bed file links
ln -s ../6_macs2_peak_calling_2_for_mESC_TKO_CTCF_rep123_together/mESC_TKO_IP_CTCF_rep123_TvsC/mESC_TKO_IP_CTCF.narrowPeaks.sorted.bed ./mESC_TKO_IP_CTCF.narrowPeaks.sorted.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_TKO_CTCF_rep123_together/mESC_TKO_IP_CTCF_rep123_TvsC/mESC_TKO_IP_CTCF.peakSummits.eB150.sorted.bed ./mESC_TKO_IP_CTCF.peakSummits.eB150.sorted.bed
ls -1 *.bed
wc -l *.bed
# 59372 mESC_SKO_IP_CTCF.narrowPeaks.sorted.bed
# 61027 mESC_SKO_IP_CTCF.peakSummits.eB150.sorted.bed
# 62238 mESC_TKO_IP_CTCF.narrowPeaks.sorted.bed
# 64121 mESC_TKO_IP_CTCF.peakSummits.eB150.sorted.bed
# 71075 mESC_WT_IP_CTCF.narrowPeaks.sorted.bed
# 73125 mESC_WT_IP_CTCF.peakSummits.eB150.sorted.bed

# intervene plot
## mESC WT/STO/TKO CTCF peaks venn
ls -1 *narrowPeaks.sorted.bed
# mESC_WT_IP_CTCF.narrowPeaks.sorted.bed
# mESC_SKO_IP_CTCF.narrowPeaks.sorted.bed
# mESC_TKO_IP_CTCF.narrowPeaks.sorted.bed
intervene venn -i mESC_WT_IP_CTCF.narrowPeaks.sorted.bed mESC_SKO_IP_CTCF.narrowPeaks.sorted.bed mESC_TKO_IP_CTCF.narrowPeaks.sorted.bed --names=mESC_WT_IP_CTCF_peaks,mESC_SKO_IP_CTCF_peaks,mESC_TKO_IP_CTCF_peaks -o 0_mESC_WT_SKO_TKO_IP_CTCF_peaks_venn --save-overlaps

## mESC WT/STO/TKO CTCF summits eB150 venn
ls -1 *peakSummits.eB150.sorted.bed
# mESC_WT_IP_CTCF.peakSummits.eB150.sorted.bed
# mESC_SKO_IP_CTCF.peakSummits.eB150.sorted.bed
# mESC_TKO_IP_CTCF.peakSummits.eB150.sorted.bed
intervene venn -i mESC_WT_IP_CTCF.peakSummits.eB150.sorted.bed mESC_SKO_IP_CTCF.peakSummits.eB150.sorted.bed mESC_TKO_IP_CTCF.peakSummits.eB150.sorted.bed --names=mESC_WT_IP_CTCF_seB150,mESC_SKO_IP_CTCF_seB150,mESC_TKO_IP_CTCF_seB150 -o 0_mESC_WT_SKO_TKO_IP_CTCF_summits_eB150_venn --save-overlaps

# 面积权重韦恩图统一在项目根目录 0_make_venn_plot.ipynb 中绘制

```

### bedtools intersection for peaks/summits eB150 from the three groups of experiments

#### SKO_vs_WT

##### peaks intersect

```bash
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together
# common from WT
bedtools intersect -a mESC_WT_IP_CTCF.narrowPeaks.sorted.bed -b mESC_SKO_IP_CTCF.narrowPeaks.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7!=0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 55504
# common from SKO
bedtools intersect -a mESC_SKO_IP_CTCF.narrowPeaks.sorted.bed -b mESC_WT_IP_CTCF.narrowPeaks.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7!=0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 55461

# only in WT
bedtools intersect -a mESC_WT_IP_CTCF.narrowPeaks.sorted.bed -b mESC_SKO_IP_CTCF.narrowPeaks.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 15571
# 
# common
bedtools intersect -a mESC_WT_IP_CTCF.narrowPeaks.sorted.bed -b mESC_SKO_IP_CTCF.narrowPeaks.sorted.bed | sort -k1,1 -k2,2n |uniq |wc -l
# 55629
# 
# only in SKO
bedtools intersect -a mESC_SKO_IP_CTCF.narrowPeaks.sorted.bed -b mESC_WT_IP_CTCF.narrowPeaks.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 3911

# 输出三类 region 及 centers bed
## only in WT
bedtools intersect -a mESC_WT_IP_CTCF.narrowPeaks.sorted.bed -b mESC_SKO_IP_CTCF.narrowPeaks.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq >mESC_WT_SKO_intersect.WT_only.regions.bed
awk 'BEGIN{FS="\t";OFS="\t"}{center=int(($2+$3)/2); print $1, center, center+1}' mESC_WT_SKO_intersect.WT_only.regions.bed >mESC_WT_SKO_intersect.WT_only.region.centers.bed
# common
bedtools intersect -a mESC_WT_IP_CTCF.narrowPeaks.sorted.bed -b mESC_SKO_IP_CTCF.narrowPeaks.sorted.bed | sort -k1,1 -k2,2n |uniq >mESC_WT_SKO_intersect.common.regions.bed
awk 'BEGIN{FS="\t";OFS="\t"}{center=int(($2+$3)/2); print $1, center, center+1}' mESC_WT_SKO_intersect.common.regions.bed >mESC_WT_SKO_intersect.common.region.centers.bed
# only in SKO
bedtools intersect -a mESC_SKO_IP_CTCF.narrowPeaks.sorted.bed -b mESC_WT_IP_CTCF.narrowPeaks.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq >mESC_WT_SKO_intersect.SKO_only.regions.bed
awk 'BEGIN{FS="\t";OFS="\t"}{center=int(($2+$3)/2); print $1, center, center+1}' mESC_WT_SKO_intersect.SKO_only.regions.bed >mESC_WT_SKO_intersect.SKO_only.region.centers.bed

```

##### summits eB150 intersect

```bash
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together
# common from WT
bedtools intersect -a mESC_WT_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_SKO_IP_CTCF.peakSummits.eB150.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7!=0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 56932
# common from SKO
bedtools intersect -a mESC_SKO_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_WT_IP_CTCF.peakSummits.eB150.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7!=0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 56826

# only in WT
bedtools intersect -a mESC_WT_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_SKO_IP_CTCF.peakSummits.eB150.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 16193
# 
# common
bedtools intersect -a mESC_WT_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_SKO_IP_CTCF.peakSummits.eB150.sorted.bed | sort -k1,1 -k2,2n |uniq |wc -l
# 57376
# 
# only in SKO
bedtools intersect -a mESC_SKO_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_WT_IP_CTCF.peakSummits.eB150.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 4201

# 输出三类 summits bed
## only in WT
bedtools intersect -a mESC_WT_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_SKO_IP_CTCF.peakSummits.eB150.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |awk 'BEGIN{FS="\t";OFS="\t"}{center=int(($2+$3)/2); print $1, center, center+1}' >mESC_WT_SKO_PSeB150_intersect.WT_only.summits.bed
# common
bedtools intersect -a mESC_WT_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_SKO_IP_CTCF.peakSummits.eB150.sorted.bed | sort -k1,1 -k2,2n |uniq |awk 'BEGIN{FS="\t";OFS="\t"}{center=int(($2+$3)/2); print $1, center, center+1}'>mESC_WT_SKO_PSeB150_intersect.common.summits.bed
# only in SKO
bedtools intersect -a mESC_SKO_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_WT_IP_CTCF.peakSummits.eB150.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |awk 'BEGIN{FS="\t";OFS="\t"}{center=int(($2+$3)/2); print $1, center, center+1}' >mESC_WT_SKO_PSeB150_intersect.SKO_only.summits.bed

```



#### TKO_vs_WT

##### peaks intersect

```bash
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together
# common from WT
bedtools intersect -a mESC_WT_IP_CTCF.narrowPeaks.sorted.bed -b mESC_TKO_IP_CTCF.narrowPeaks.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7!=0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 57583
# common from TKO
bedtools intersect -a mESC_TKO_IP_CTCF.narrowPeaks.sorted.bed -b mESC_WT_IP_CTCF.narrowPeaks.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7!=0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 57575

# only in WT
bedtools intersect -a mESC_WT_IP_CTCF.narrowPeaks.sorted.bed -b mESC_TKO_IP_CTCF.narrowPeaks.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 13492
# 
# common
bedtools intersect -a mESC_WT_IP_CTCF.narrowPeaks.sorted.bed -b mESC_TKO_IP_CTCF.narrowPeaks.sorted.bed | sort -k1,1 -k2,2n |uniq |wc -l
# 57713
# 
# only in TKO
bedtools intersect -a mESC_TKO_IP_CTCF.narrowPeaks.sorted.bed -b mESC_WT_IP_CTCF.narrowPeaks.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 4663

# 输出三类 region 及 centers bed
## only in WT
bedtools intersect -a mESC_WT_IP_CTCF.narrowPeaks.sorted.bed -b mESC_TKO_IP_CTCF.narrowPeaks.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq >mESC_WT_TKO_intersect.WT_only.regions.bed
awk 'BEGIN{FS="\t";OFS="\t"}{center=int(($2+$3)/2); print $1, center, center+1}' mESC_WT_TKO_intersect.WT_only.regions.bed >mESC_WT_TKO_intersect.WT_only.region.centers.bed
# common
bedtools intersect -a mESC_WT_IP_CTCF.narrowPeaks.sorted.bed -b mESC_TKO_IP_CTCF.narrowPeaks.sorted.bed | sort -k1,1 -k2,2n |uniq >mESC_WT_TKO_intersect.common.regions.bed
awk 'BEGIN{FS="\t";OFS="\t"}{center=int(($2+$3)/2); print $1, center, center+1}' mESC_WT_TKO_intersect.common.regions.bed >mESC_WT_TKO_intersect.common.region.centers.bed
# only in TKO
bedtools intersect -a mESC_TKO_IP_CTCF.narrowPeaks.sorted.bed -b mESC_WT_IP_CTCF.narrowPeaks.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq >mESC_WT_TKO_intersect.TKO_only.regions.bed
awk 'BEGIN{FS="\t";OFS="\t"}{center=int(($2+$3)/2); print $1, center, center+1}' mESC_WT_TKO_intersect.TKO_only.regions.bed >mESC_WT_TKO_intersect.TKO_only.region.centers.bed

```

##### summits eB150 intersect

```bash
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together
# common from WT
bedtools intersect -a mESC_WT_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_TKO_IP_CTCF.peakSummits.eB150.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7!=0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 59124
# common from TKO
bedtools intersect -a mESC_TKO_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_WT_IP_CTCF.peakSummits.eB150.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7!=0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 59096

# only in WT
bedtools intersect -a mESC_WT_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_TKO_IP_CTCF.peakSummits.eB150.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 14001
# 
# common
bedtools intersect -a mESC_WT_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_TKO_IP_CTCF.peakSummits.eB150.sorted.bed | sort -k1,1 -k2,2n |uniq |wc -l
# 59708
# 
# only in TKO
bedtools intersect -a mESC_TKO_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_WT_IP_CTCF.peakSummits.eB150.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 5025

# 输出三类 summits bed
## only in WT
bedtools intersect -a mESC_WT_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_TKO_IP_CTCF.peakSummits.eB150.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |awk 'BEGIN{FS="\t";OFS="\t"}{center=int(($2+$3)/2); print $1, center, center+1}' >mESC_WT_TKO_PSeB150_intersect.WT_only.summits.bed
# common
bedtools intersect -a mESC_WT_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_TKO_IP_CTCF.peakSummits.eB150.sorted.bed | sort -k1,1 -k2,2n |uniq |awk 'BEGIN{FS="\t";OFS="\t"}{center=int(($2+$3)/2); print $1, center, center+1}'>mESC_WT_TKO_PSeB150_intersect.common.summits.bed
# only in TKO
bedtools intersect -a mESC_TKO_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_WT_IP_CTCF.peakSummits.eB150.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |awk 'BEGIN{FS="\t";OFS="\t"}{center=int(($2+$3)/2); print $1, center, center+1}' >mESC_WT_TKO_PSeB150_intersect.TKO_only.summits.bed

```



#### SKO_vs_TKO

##### peaks intersect

```bash
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together
# common from SKO
bedtools intersect -a mESC_SKO_IP_CTCF.narrowPeaks.sorted.bed -b mESC_TKO_IP_CTCF.narrowPeaks.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7!=0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 53207
# common from TKO
bedtools intersect -a mESC_TKO_IP_CTCF.narrowPeaks.sorted.bed -b mESC_SKO_IP_CTCF.narrowPeaks.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7!=0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 53242

# only in SKO
bedtools intersect -a mESC_SKO_IP_CTCF.narrowPeaks.sorted.bed -b mESC_TKO_IP_CTCF.narrowPeaks.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 6165
# 
# common
bedtools intersect -a mESC_SKO_IP_CTCF.narrowPeaks.sorted.bed -b mESC_TKO_IP_CTCF.narrowPeaks.sorted.bed | sort -k1,1 -k2,2n |uniq |wc -l
# 53357
# 
# only in TKO
bedtools intersect -a mESC_TKO_IP_CTCF.narrowPeaks.sorted.bed -b mESC_SKO_IP_CTCF.narrowPeaks.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 8996

# 输出三类 region 及 centers bed
## only in SKO
bedtools intersect -a mESC_SKO_IP_CTCF.narrowPeaks.sorted.bed -b mESC_TKO_IP_CTCF.narrowPeaks.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq >mESC_SKO_TKO_intersect.SKO_only.regions.bed
awk 'BEGIN{FS="\t";OFS="\t"}{center=int(($2+$3)/2); print $1, center, center+1}' mESC_SKO_TKO_intersect.SKO_only.regions.bed >mESC_SKO_TKO_intersect.SKO_only.region.centers.bed
# common
bedtools intersect -a mESC_SKO_IP_CTCF.narrowPeaks.sorted.bed -b mESC_TKO_IP_CTCF.narrowPeaks.sorted.bed | sort -k1,1 -k2,2n |uniq >mESC_SKO_TKO_intersect.common.regions.bed
awk 'BEGIN{FS="\t";OFS="\t"}{center=int(($2+$3)/2); print $1, center, center+1}' mESC_SKO_TKO_intersect.common.regions.bed >mESC_SKO_TKO_intersect.common.region.centers.bed
# only in TKO
bedtools intersect -a mESC_TKO_IP_CTCF.narrowPeaks.sorted.bed -b mESC_SKO_IP_CTCF.narrowPeaks.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq >mESC_SKO_TKO_intersect.TKO_only.regions.bed
awk 'BEGIN{FS="\t";OFS="\t"}{center=int(($2+$3)/2); print $1, center, center+1}' mESC_SKO_TKO_intersect.TKO_only.regions.bed >mESC_SKO_TKO_intersect.TKO_only.region.centers.bed

```

##### summits eB150 intersect

```bash
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together
# common from SKO
bedtools intersect -a mESC_SKO_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_TKO_IP_CTCF.peakSummits.eB150.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7!=0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 54509
# common from TKO
bedtools intersect -a mESC_TKO_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_SKO_IP_CTCF.peakSummits.eB150.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7!=0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 54583

# only in SKO
bedtools intersect -a mESC_SKO_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_TKO_IP_CTCF.peakSummits.eB150.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 6518
# 
# common
bedtools intersect -a mESC_SKO_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_TKO_IP_CTCF.peakSummits.eB150.sorted.bed | sort -k1,1 -k2,2n |uniq |wc -l
# 54994
# 
# only in TKO
bedtools intersect -a mESC_TKO_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_SKO_IP_CTCF.peakSummits.eB150.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |wc -l
# 9538

# 输出三类 summits bed
## only in SKO
bedtools intersect -a mESC_SKO_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_TKO_IP_CTCF.peakSummits.eB150.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |awk 'BEGIN{FS="\t";OFS="\t"}{center=int(($2+$3)/2); print $1, center, center+1}' >mESC_SKO_TKO_PSeB150_intersect.SKO_only.summits.bed
# common
bedtools intersect -a mESC_SKO_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_TKO_IP_CTCF.peakSummits.eB150.sorted.bed | sort -k1,1 -k2,2n |uniq |awk 'BEGIN{FS="\t";OFS="\t"}{center=int(($2+$3)/2); print $1, center, center+1}'>mESC_SKO_TKO_PSeB150_intersect.common.summits.bed
# only in TKO
bedtools intersect -a mESC_TKO_IP_CTCF.peakSummits.eB150.sorted.bed -b mESC_SKO_IP_CTCF.peakSummits.eB150.sorted.bed -wao |awk 'BEGIN{FS="\t";OFS="\t"}{if($7==0){print $1,$2,$3}}' |sort -k1,1 -k2,2n |uniq |awk 'BEGIN{FS="\t";OFS="\t"}{center=int(($2+$3)/2); print $1, center, center+1}' >mESC_SKO_TKO_PSeB150_intersect.TKO_only.summits.bed

```



# use bamCompare to reate log2 transformed bigwig files

```bash
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
mkdir 8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap
ls -1 ../*rep[123].bowtie2.sorted.bam
ls -1 ../*rep123.merged.sorted.bam
# 根据要用的 bam 文件列表构造 bamCompare 运行脚本，上传并运行 
# 8.1_run_bamCompare_for_WT_SKO_TKO_each_repeat.sh
# 8.2_run_bamCompare_for_WT_SKO_TKO_rep123_merged_bam.sh

# 由于使用原始的 encode_mm10_black_list.bed 软件运行报错
# Your blacklist file(s) has (have) regions that overlap. 
# Proceeding with such a file would result in deepTools incorrectly calculating scaling factors. 
# As such, you MUST fix this issue before being able to proceed.
# 使用如下命令先对 black_list 的重叠区间进行合并然排序去重，最后在 8.1 和 8.2 脚本中使用新构造的 encode_mm10_black_list.sorted.bed
cd /home/sunwenju/3_genomeData/mouse/mm10
bedtools merge -i encode_mm10_black_list.bed |sort -k1,1 -k2,2n |uniq >encode_mm10_black_list.sorted.bed

# 分别运行 8.1 和 8.2 构造 TvsC bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap
nohup time bash 8.1_run_bamCompare_for_WT_SKO_TKO_each_repeat.sh >8.1_run_bamCompare_for_WT_SKO_TKO_each_repeat.sh.runlog 2>&1 &
tail -f 8.1_run_bamCompare_for_WT_SKO_TKO_each_repeat.sh.runlog

cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap
nohup time bash 8.2_run_bamCompare_for_WT_SKO_TKO_rep123_merged_bam.sh >8.2_run_bamCompare_for_WT_SKO_TKO_rep123_merged_bam.sh.runlog 2>&1 &
tail -f 8.2_run_bamCompare_for_WT_SKO_TKO_rep123_merged_bam.sh.runlog

```

# peak summits reads enrichment heatmap

## plotHeatMap_for_all_use_mESC_WT_CTCF_summits

### data preperation

```bash
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
mkdir 9_plotHeatMap_for_all_use_mESC_WT_CTCF_summits
cd 9_plotHeatMap_for_all_use_mESC_WT_CTCF_summits
# 创建bigwig文件连接
## rep1/2/3 individual
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF_rep1.TvsC.bigwig mESC_WT_IP_CTCF_rep1.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF_rep2.TvsC.bigwig mESC_WT_IP_CTCF_rep2.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF_rep3.TvsC.bigwig mESC_WT_IP_CTCF_rep3.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF_rep1.TvsC.bigwig mESC_SKO_IP_CTCF_rep1.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF_rep2.TvsC.bigwig mESC_SKO_IP_CTCF_rep2.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF_rep3.TvsC.bigwig mESC_SKO_IP_CTCF_rep3.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF_rep1.TvsC.bigwig mESC_TKO_IP_CTCF_rep1.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF_rep2.TvsC.bigwig mESC_TKO_IP_CTCF_rep2.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF_rep3.TvsC.bigwig mESC_TKO_IP_CTCF_rep3.TvsC.bigwig
## rep123 merged
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig
ls -1 *rep[123].TvsC.bigwig
ls -1 *rep123_merged.TvsC.bigwig

# 创建 mESC_WT_CTCF_summits bed 文件连接
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/mESC_WT_IP_CTCF_rep123_TvsC/mESC_WT_IP_CTCF_rep123_TvsC_summits.bed mESC_WT_IP_CTCF_summits.bed
head mESC_WT_IP_CTCF_summits.bed
wc -l mESC_WT_IP_CTCF_summits.bed
# 73125 mESC_WT_IP_CTCF_summits.bed
cut -f 1,2,3 mESC_WT_IP_CTCF_summits.bed |head
cut -f 1,2,3 mESC_WT_IP_CTCF_summits.bed |sort -k1,1 -k2,2n |uniq >mESC_WT_IP_CTCF_summits.sorted.bed
head mESC_WT_IP_CTCF_summits.sorted.bed
wc -l mESC_WT_IP_CTCF_summits.sorted.bed
# 73125 mESC_WT_IP_CTCF_summits.sorted.bed

```

### computeMatrix and plotHeatmap

```bash
## rep1/2/3 individual bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/9_plotHeatMap_for_all_use_mESC_WT_CTCF_summits
computeMatrix reference-point -p "max" -R mESC_WT_IP_CTCF_summits.sorted.bed -S mESC_WT_IP_CTCF_rep1.TvsC.bigwig mESC_WT_IP_CTCF_rep2.TvsC.bigwig mESC_WT_IP_CTCF_rep3.TvsC.bigwig mESC_SKO_IP_CTCF_rep1.TvsC.bigwig mESC_SKO_IP_CTCF_rep2.TvsC.bigwig mESC_SKO_IP_CTCF_rep3.TvsC.bigwig mESC_TKO_IP_CTCF_rep1.TvsC.bigwig mESC_TKO_IP_CTCF_rep2.TvsC.bigwig mESC_TKO_IP_CTCF_rep3.TvsC.bigwig --referencePoint center -b 1000 -a 1000 --sortRegions descend --sortUsingSamples 1 2 3 --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3" -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.mESC_WT_CTCF_summits.gz --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.mESC_WT_CTCF_summits.matrix --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.mESC_WT_CTCF_summits.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.mESC_WT_CTCF_summits.gz -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.mESC_WT_CTCF_summits.plot.pdf --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" --regionsLabel "mESC_WT_CTCF_summits(73125)" --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3"

## rep123 merged bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/9_plotHeatMap_for_all_use_mESC_WT_CTCF_summits
computeMatrix reference-point -p "max" -R mESC_WT_IP_CTCF_summits.sorted.bed -S mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig --referencePoint center -b 1000 -a 1000 --sortRegions descend  --sortUsingSamples 1 --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF" -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.mESC_WT_CTCF_summits.gz --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.mESC_WT_CTCF_summits.matrix --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.mESC_WT_CTCF_summits.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.mESC_WT_CTCF_summits.gz -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.mESC_WT_CTCF_summits.plot.pdf --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" --regionsLabel "mESC_WT_CTCF_summits(73125)" --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF"

```



## plotHeatMap_for_all_use_mESC_SKO_CTCF_summits

### data preperation

```bash
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
mkdir 9_plotHeatMap_for_all_use_mESC_SKO_CTCF_summits
cd 9_plotHeatMap_for_all_use_mESC_SKO_CTCF_summits
# 创建bigwig文件连接
## rep1/2/3 individual
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF_rep1.TvsC.bigwig mESC_WT_IP_CTCF_rep1.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF_rep2.TvsC.bigwig mESC_WT_IP_CTCF_rep2.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF_rep3.TvsC.bigwig mESC_WT_IP_CTCF_rep3.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF_rep1.TvsC.bigwig mESC_SKO_IP_CTCF_rep1.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF_rep2.TvsC.bigwig mESC_SKO_IP_CTCF_rep2.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF_rep3.TvsC.bigwig mESC_SKO_IP_CTCF_rep3.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF_rep1.TvsC.bigwig mESC_TKO_IP_CTCF_rep1.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF_rep2.TvsC.bigwig mESC_TKO_IP_CTCF_rep2.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF_rep3.TvsC.bigwig mESC_TKO_IP_CTCF_rep3.TvsC.bigwig
## rep123 merged
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig
ls -1 *rep[123].TvsC.bigwig
ls -1 *rep123_merged.TvsC.bigwig

# 创建 mESC_SKO_CTCF_summits bed 文件连接
ln -s ../6_macs2_peak_calling_2_for_mESC_SKO_CTCF_rep123_together/mESC_SKO_IP_CTCF_rep123_TvsC/mESC_SKO_IP_CTCF_rep123_TvsC_summits.bed mESC_SKO_IP_CTCF_summits.bed
head mESC_SKO_IP_CTCF_summits.bed
wc -l mESC_SKO_IP_CTCF_summits.bed
# 61027 mESC_SKO_IP_CTCF_summits.bed
cut -f 1,2,3 mESC_SKO_IP_CTCF_summits.bed |head
cut -f 1,2,3 mESC_SKO_IP_CTCF_summits.bed |sort -k1,1 -k2,2n |uniq >mESC_SKO_IP_CTCF_summits.sorted.bed
head mESC_SKO_IP_CTCF_summits.sorted.bed
wc -l mESC_SKO_IP_CTCF_summits.sorted.bed
# 61027 mESC_SKO_IP_CTCF_summits.sorted.bed

```

### computeMatrix and plotHeatmap

```bash
## rep1/2/3 individual bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/9_plotHeatMap_for_all_use_mESC_SKO_CTCF_summits
computeMatrix reference-point -p "max" -R mESC_SKO_IP_CTCF_summits.sorted.bed -S mESC_WT_IP_CTCF_rep1.TvsC.bigwig mESC_WT_IP_CTCF_rep2.TvsC.bigwig mESC_WT_IP_CTCF_rep3.TvsC.bigwig mESC_SKO_IP_CTCF_rep1.TvsC.bigwig mESC_SKO_IP_CTCF_rep2.TvsC.bigwig mESC_SKO_IP_CTCF_rep3.TvsC.bigwig mESC_TKO_IP_CTCF_rep1.TvsC.bigwig mESC_TKO_IP_CTCF_rep2.TvsC.bigwig mESC_TKO_IP_CTCF_rep3.TvsC.bigwig --referencePoint center -b 1000 -a 1000 --sortRegions descend --sortUsingSamples 1 2 3 --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3" -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.mESC_SKO_CTCF_summits.gz --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.mESC_SKO_CTCF_summits.matrix --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.mESC_SKO_CTCF_summits.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.mESC_SKO_CTCF_summits.gz -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.mESC_SKO_CTCF_summits.plot.pdf --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" --regionsLabel "mESC_SKO_CTCF_summits(61027)" --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3"

## rep123 merged bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/9_plotHeatMap_for_all_use_mESC_SKO_CTCF_summits
computeMatrix reference-point -p "max" -R mESC_SKO_IP_CTCF_summits.sorted.bed -S mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig --referencePoint center -b 1000 -a 1000 --sortRegions descend  --sortUsingSamples 1 --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF" -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.mESC_SKO_CTCF_summits.gz --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.mESC_SKO_CTCF_summits.matrix --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.mESC_SKO_CTCF_summits.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.mESC_SKO_CTCF_summits.gz -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.mESC_SKO_CTCF_summits.plot.pdf --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" --regionsLabel "mESC_SKO_CTCF_summits(61027)" --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF"

```



## plotHeatMap_for_all_use_mESC_TKO_CTCF_summits

### data preperation

```bash
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
mkdir 9_plotHeatMap_for_all_use_mESC_TKO_CTCF_summits
cd 9_plotHeatMap_for_all_use_mESC_TKO_CTCF_summits
# 创建bigwig文件连接
## rep1/2/3 individual
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF_rep1.TvsC.bigwig mESC_WT_IP_CTCF_rep1.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF_rep2.TvsC.bigwig mESC_WT_IP_CTCF_rep2.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF_rep3.TvsC.bigwig mESC_WT_IP_CTCF_rep3.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF_rep1.TvsC.bigwig mESC_SKO_IP_CTCF_rep1.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF_rep2.TvsC.bigwig mESC_SKO_IP_CTCF_rep2.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF_rep3.TvsC.bigwig mESC_SKO_IP_CTCF_rep3.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF_rep1.TvsC.bigwig mESC_TKO_IP_CTCF_rep1.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF_rep2.TvsC.bigwig mESC_TKO_IP_CTCF_rep2.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF_rep3.TvsC.bigwig mESC_TKO_IP_CTCF_rep3.TvsC.bigwig
## rep123 merged
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig
ls -1 *rep[123].TvsC.bigwig
ls -1 *rep123_merged.TvsC.bigwig

# 创建 mESC_TKO_CTCF_summits bed 文件连接
ln -s ../6_macs2_peak_calling_2_for_mESC_TKO_CTCF_rep123_together/mESC_TKO_IP_CTCF_rep123_TvsC/mESC_TKO_IP_CTCF_rep123_TvsC_summits.bed mESC_TKO_IP_CTCF_summits.bed
head mESC_TKO_IP_CTCF_summits.bed
wc -l mESC_TKO_IP_CTCF_summits.bed
# 64121 mESC_TKO_IP_CTCF_summits.bed
cut -f 1,2,3 mESC_TKO_IP_CTCF_summits.bed |head
cut -f 1,2,3 mESC_TKO_IP_CTCF_summits.bed |sort -k1,1 -k2,2n |uniq >mESC_TKO_IP_CTCF_summits.sorted.bed
head mESC_TKO_IP_CTCF_summits.sorted.bed
wc -l mESC_TKO_IP_CTCF_summits.sorted.bed
# 64121 mESC_TKO_IP_CTCF_summits.sorted.bed

```

### computeMatrix and plotHeatmap

```bash
## rep1/2/3 individual bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/9_plotHeatMap_for_all_use_mESC_TKO_CTCF_summits
computeMatrix reference-point -p "max" -R mESC_TKO_IP_CTCF_summits.sorted.bed -S mESC_WT_IP_CTCF_rep1.TvsC.bigwig mESC_WT_IP_CTCF_rep2.TvsC.bigwig mESC_WT_IP_CTCF_rep3.TvsC.bigwig mESC_SKO_IP_CTCF_rep1.TvsC.bigwig mESC_SKO_IP_CTCF_rep2.TvsC.bigwig mESC_SKO_IP_CTCF_rep3.TvsC.bigwig mESC_TKO_IP_CTCF_rep1.TvsC.bigwig mESC_TKO_IP_CTCF_rep2.TvsC.bigwig mESC_TKO_IP_CTCF_rep3.TvsC.bigwig --referencePoint center -b 1000 -a 1000 --sortRegions descend --sortUsingSamples 1 2 3 --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3" -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.mESC_TKO_CTCF_summits.gz --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.mESC_TKO_CTCF_summits.matrix --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.mESC_TKO_CTCF_summits.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.mESC_TKO_CTCF_summits.gz -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.mESC_TKO_CTCF_summits.plot.pdf --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" --regionsLabel "mESC_TKO_CTCF_summits(64121)" --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3"

## rep123 merged bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/9_plotHeatMap_for_all_use_mESC_TKO_CTCF_summits
computeMatrix reference-point -p "max" -R mESC_TKO_IP_CTCF_summits.sorted.bed -S mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig --referencePoint center -b 1000 -a 1000 --sortRegions descend  --sortUsingSamples 1 --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF" -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.mESC_TKO_CTCF_summits.gz --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.mESC_TKO_CTCF_summits.matrix --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.mESC_TKO_CTCF_summits.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.mESC_TKO_CTCF_summits.gz -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.mESC_TKO_CTCF_summits.plot.pdf --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" --regionsLabel "mESC_TKO_CTCF_summits(64121)" --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF"

```

# 三组实验两种方法各自三类 peak center 周围信号富集热图

## WT_vs_SKO

### 创建 TvsC_bigwig 与 两种方法各自三类 peak center 文件连接

```bash
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
mkdir 10_plotHeatMap_for_all_use_mESC_WTvsSKO_CTCF_summits
cd 10_plotHeatMap_for_all_use_mESC_WTvsSKO_CTCF_summits
# 创建bigwig文件连接
## rep1/2/3 individual
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF_rep1.TvsC.bigwig mESC_WT_IP_CTCF_rep1.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF_rep2.TvsC.bigwig mESC_WT_IP_CTCF_rep2.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF_rep3.TvsC.bigwig mESC_WT_IP_CTCF_rep3.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF_rep1.TvsC.bigwig mESC_SKO_IP_CTCF_rep1.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF_rep2.TvsC.bigwig mESC_SKO_IP_CTCF_rep2.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF_rep3.TvsC.bigwig mESC_SKO_IP_CTCF_rep3.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF_rep1.TvsC.bigwig mESC_TKO_IP_CTCF_rep1.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF_rep2.TvsC.bigwig mESC_TKO_IP_CTCF_rep2.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF_rep3.TvsC.bigwig mESC_TKO_IP_CTCF_rep3.TvsC.bigwig
## rep123 merged
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig
ls -1 *rep[123].TvsC.bigwig
ls -1 *rep123_merged.TvsC.bigwig

#  创建三类 peak center 文件连接 (bedtools intersect)
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/mESC_WT_SKO_intersect.WT_only.region.centers.bed mESC_WT_SKO_intersect.WT_only.region.centers.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/mESC_WT_SKO_intersect.common.region.centers.bed mESC_WT_SKO_intersect.common.region.centers.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/mESC_WT_SKO_intersect.SKO_only.region.centers.bed mESC_WT_SKO_intersect.SKO_only.region.centers.bed
wc -l *.region.centers.bed
# 15571 mESC_WT_SKO_intersect.WT_only.region.centers.bed
# 55629 mESC_WT_SKO_intersect.common.region.centers.bed
#  3911 mESC_WT_SKO_intersect.SKO_only.region.centers.bed

#  创建三类 peak center 文件连接 (bdgdiff)
ln -s ../7_macs2_peak_calling_3_for_SKO_vs_WT_bdgdiff/RIB_WT_SKO.CTCF.summit.sorted.bed RIB_WT_SKO.CTCF.summit.sorted.bed
ln -s ../7_macs2_peak_calling_3_for_SKO_vs_WT_bdgdiff/ROI_SKO.CTCF.summit.sorted.bed ROI_SKO.CTCF.summit.sorted.bed
ln -s ../7_macs2_peak_calling_3_for_SKO_vs_WT_bdgdiff/ROI_WT.CTCF.summit.sorted.bed ROI_WT.CTCF.summit.sorted.bed
wc -l *.summit.sorted.bed
#  3379 ROI_WT.CTCF.summit.sorted.bed
# 48821 RIB_WT_SKO.CTCF.summit.sorted.bed
#  1711 ROI_SKO.CTCF.summit.sorted.bed

```

### computeMatrix and plotHeatmap

#### use bedtools intersect 三类 peak center

```bash
# use bedtools intersect 三类 peak center
## rep123 merged bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/10_plotHeatMap_for_all_use_mESC_WTvsSKO_CTCF_summits
computeMatrix reference-point -p "max" -R mESC_WT_SKO_intersect.WT_only.region.centers.bed mESC_WT_SKO_intersect.common.region.centers.bed mESC_WT_SKO_intersect.SKO_only.region.centers.bed \
    -S mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig \
    --referencePoint center -b 1000 -a 1000 --sortRegions descend  --sortUsingSamples 1 --missingDataAsZero --skipZeros \
    --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed \
    --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF" \
    -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKO.intersect_3cPeakCenter.gz \
    --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKO.intersect_3cPeakCenter.matrix \
    --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKO.intersect_3cPeakCenter.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKO.intersect_3cPeakCenter.gz \
    -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKO.intersect_3cPeakCenter.plot.pdf \
    --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 \
    --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" \
    --regionsLabel "WT_only (15571)" "common (55629)" "SKO_only (3911)" --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF"

## rep1/2/3 individual bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/10_plotHeatMap_for_all_use_mESC_WTvsSKO_CTCF_summits
computeMatrix reference-point -p "max" -R mESC_WT_SKO_intersect.WT_only.region.centers.bed mESC_WT_SKO_intersect.common.region.centers.bed mESC_WT_SKO_intersect.SKO_only.region.centers.bed \
    -S mESC_WT_IP_CTCF_rep1.TvsC.bigwig mESC_WT_IP_CTCF_rep2.TvsC.bigwig mESC_WT_IP_CTCF_rep3.TvsC.bigwig \
    mESC_SKO_IP_CTCF_rep1.TvsC.bigwig mESC_SKO_IP_CTCF_rep2.TvsC.bigwig mESC_SKO_IP_CTCF_rep3.TvsC.bigwig \
    mESC_TKO_IP_CTCF_rep1.TvsC.bigwig mESC_TKO_IP_CTCF_rep2.TvsC.bigwig mESC_TKO_IP_CTCF_rep3.TvsC.bigwig \
    --referencePoint center -b 1000 -a 1000 --sortRegions descend --sortUsingSamples 1 2 3 --missingDataAsZero --skipZeros \
    --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed \
    --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" \
    "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" \
    "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3" \
    -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKO.intersect_3cPeakCenter.gz \
    --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKO.intersect_3cPeakCenter.matrix \
    --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKO.intersect_3cPeakCenter.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKO.intersect_3cPeakCenter.gz \
    -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKO.intersect_3cPeakCenter.plot.pdf \
    --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 \
    --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" \
    --regionsLabel "WT_only (15571)" "common (55629)" "SKO_only (3911)" \
    --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" \
    "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" \
    "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3"

```

#### use bdgdiff 三类 peak center

```bash
# use bdgdiff 三类 peak center
## rep123 merged bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/10_plotHeatMap_for_all_use_mESC_WTvsSKO_CTCF_summits
computeMatrix reference-point -p "max" -R ROI_WT.CTCF.summit.sorted.bed RIB_WT_SKO.CTCF.summit.sorted.bed ROI_SKO.CTCF.summit.sorted.bed \
    -S mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig \
    --referencePoint center -b 1000 -a 1000 --sortRegions descend  --sortUsingSamples 1 --missingDataAsZero --skipZeros \
    --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed \
    --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF" \
    -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKO.bdgdiff_3cPeakCenter.gz \
    --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKO.bdgdiff_3cPeakCenter.matrix \
    --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKO.bdgdiff_3cPeakCenter.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKO.bdgdiff_3cPeakCenter.gz \
    -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKO.bdgdiff_3cPeakCenter.plot.pdf \
    --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 \
    --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" \
    --regionsLabel "ROI_WT (3379)" "common (48821)" "ROI_SKO (1711)" --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF"

## rep1/2/3 individual bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/10_plotHeatMap_for_all_use_mESC_WTvsSKO_CTCF_summits
computeMatrix reference-point -p "max" -R ROI_WT.CTCF.summit.sorted.bed RIB_WT_SKO.CTCF.summit.sorted.bed ROI_SKO.CTCF.summit.sorted.bed \
    -S mESC_WT_IP_CTCF_rep1.TvsC.bigwig mESC_WT_IP_CTCF_rep2.TvsC.bigwig mESC_WT_IP_CTCF_rep3.TvsC.bigwig \
    mESC_SKO_IP_CTCF_rep1.TvsC.bigwig mESC_SKO_IP_CTCF_rep2.TvsC.bigwig mESC_SKO_IP_CTCF_rep3.TvsC.bigwig \
    mESC_TKO_IP_CTCF_rep1.TvsC.bigwig mESC_TKO_IP_CTCF_rep2.TvsC.bigwig mESC_TKO_IP_CTCF_rep3.TvsC.bigwig \
    --referencePoint center -b 1000 -a 1000 --sortRegions descend --sortUsingSamples 1 2 3 --missingDataAsZero --skipZeros \
    --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed \
    --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" \
    "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" \
    "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3" \
    -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKO.bdgdiff_3cPeakCenter.gz \
    --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKO.bdgdiff_3cPeakCenter.matrix \
    --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKO.bdgdiff_3cPeakCenter.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKO.bdgdiff_3cPeakCenter.gz \
    -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKO.bdgdiff_3cPeakCenter.plot.pdf \
    --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 \
    --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" \
    --regionsLabel "ROI_WT (3379)" "common (48821)" "ROI_SKO (1711)" \
    --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" \
    "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" \
    "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3"

# 最后为pdf图生成对应的png图，方便预览
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/10_plotHeatMap_for_all_use_mESC_WTvsSKO_CTCF_summits
for i in `ls -1 bamCompare_matrix_of*.pdf`;do echo $i; pdftoppm -png $i ${i%.pdf};done

```



## WT_vs_TKO

### 创建 TvsC_bigwig 与 两种方法各自三类 peak center 文件连接
```bash
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
mkdir 10_plotHeatMap_for_all_use_mESC_WTvsTKO_CTCF_summits
cd 10_plotHeatMap_for_all_use_mESC_WTvsTKO_CTCF_summits
# 创建bigwig文件连接
## rep1/2/3 individual
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF_rep1.TvsC.bigwig mESC_WT_IP_CTCF_rep1.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF_rep2.TvsC.bigwig mESC_WT_IP_CTCF_rep2.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF_rep3.TvsC.bigwig mESC_WT_IP_CTCF_rep3.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF_rep1.TvsC.bigwig mESC_SKO_IP_CTCF_rep1.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF_rep2.TvsC.bigwig mESC_SKO_IP_CTCF_rep2.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF_rep3.TvsC.bigwig mESC_SKO_IP_CTCF_rep3.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF_rep1.TvsC.bigwig mESC_TKO_IP_CTCF_rep1.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF_rep2.TvsC.bigwig mESC_TKO_IP_CTCF_rep2.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF_rep3.TvsC.bigwig mESC_TKO_IP_CTCF_rep3.TvsC.bigwig
## rep123 merged
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig
ls -1 *rep[123].TvsC.bigwig
ls -1 *rep123_merged.TvsC.bigwig

#  创建三类 peak center 文件连接 (bedtools intersect)
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/mESC_WT_TKO_intersect.WT_only.region.centers.bed mESC_WT_TKO_intersect.WT_only.region.centers.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/mESC_WT_TKO_intersect.common.region.centers.bed mESC_WT_TKO_intersect.common.region.centers.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/mESC_WT_TKO_intersect.TKO_only.region.centers.bed mESC_WT_TKO_intersect.TKO_only.region.centers.bed
wc -l *.region.centers.bed
# 13492 mESC_WT_TKO_intersect.WT_only.region.centers.bed
# 57713 mESC_WT_TKO_intersect.common.region.centers.bed
#  4663 mESC_WT_TKO_intersect.TKO_only.region.centers.bed

#  创建三类 peak center 文件连接 (bdgdiff)
ln -s ../7_macs2_peak_calling_3_for_TKO_vs_WT_bdgdiff/RIB_WT_TKO.CTCF.summit.sorted.bed RIB_WT_TKO.CTCF.summit.sorted.bed
ln -s ../7_macs2_peak_calling_3_for_TKO_vs_WT_bdgdiff/ROI_TKO.CTCF.summit.sorted.bed ROI_TKO.CTCF.summit.sorted.bed
ln -s ../7_macs2_peak_calling_3_for_TKO_vs_WT_bdgdiff/ROI_WT.CTCF.summit.sorted.bed ROI_WT.CTCF.summit.sorted.bed
wc -l *.summit.sorted.bed
#  6436 ROI_WT.CTCF.summit.sorted.bed
# 50966 RIB_WT_TKO.CTCF.summit.sorted.bed
#   872 ROI_TKO.CTCF.summit.sorted.bed

```

### computeMatrix and plotHeatmap

#### use bedtools intersect 三类 peak center

```bash
# use bedtools intersect 三类 peak center
## rep123 merged bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/10_plotHeatMap_for_all_use_mESC_WTvsTKO_CTCF_summits
computeMatrix reference-point -p "max" -R mESC_WT_TKO_intersect.WT_only.region.centers.bed mESC_WT_TKO_intersect.common.region.centers.bed mESC_WT_TKO_intersect.TKO_only.region.centers.bed \
    -S mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig \
    --referencePoint center -b 1000 -a 1000 --sortRegions descend  --sortUsingSamples 1 --missingDataAsZero --skipZeros \
    --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed \
    --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF" \
    -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsTKO.intersect_3cPeakCenter.gz \
    --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsTKO.intersect_3cPeakCenter.matrix \
    --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsTKO.intersect_3cPeakCenter.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsTKO.intersect_3cPeakCenter.gz \
    -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsTKO.intersect_3cPeakCenter.plot.pdf \
    --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 \
    --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" \
    --regionsLabel "WT_only (13492)" "common (57713)" "SKO_only (4663)" --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF"

## rep1/2/3 individual bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/10_plotHeatMap_for_all_use_mESC_WTvsTKO_CTCF_summits
computeMatrix reference-point -p "max" -R mESC_WT_TKO_intersect.WT_only.region.centers.bed mESC_WT_TKO_intersect.common.region.centers.bed mESC_WT_TKO_intersect.TKO_only.region.centers.bed \
    -S mESC_WT_IP_CTCF_rep1.TvsC.bigwig mESC_WT_IP_CTCF_rep2.TvsC.bigwig mESC_WT_IP_CTCF_rep3.TvsC.bigwig \
    mESC_SKO_IP_CTCF_rep1.TvsC.bigwig mESC_SKO_IP_CTCF_rep2.TvsC.bigwig mESC_SKO_IP_CTCF_rep3.TvsC.bigwig \
    mESC_TKO_IP_CTCF_rep1.TvsC.bigwig mESC_TKO_IP_CTCF_rep2.TvsC.bigwig mESC_TKO_IP_CTCF_rep3.TvsC.bigwig \
    --referencePoint center -b 1000 -a 1000 --sortRegions descend --sortUsingSamples 1 2 3 --missingDataAsZero --skipZeros \
    --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed \
    --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" \
    "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" \
    "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3" \
    -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsTKO.intersect_3cPeakCenter.gz \
    --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsTKO.intersect_3cPeakCenter.matrix \
    --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsTKO.intersect_3cPeakCenter.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsTKO.intersect_3cPeakCenter.gz \
    -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsTKO.intersect_3cPeakCenter.plot.pdf \
    --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 \
    --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" \
    --regionsLabel "WT_only (13492)" "common (57713)" "SKO_only (4663)" \
    --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" \
    "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" \
    "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3"

```

#### use bdgdiff 三类 peak center

```bash
# use bdgdiff 三类 peak center
## rep123 merged bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/10_plotHeatMap_for_all_use_mESC_WTvsTKO_CTCF_summits
computeMatrix reference-point -p "max" -R ROI_WT.CTCF.summit.sorted.bed RIB_WT_TKO.CTCF.summit.sorted.bed ROI_TKO.CTCF.summit.sorted.bed \
    -S mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig \
    --referencePoint center -b 1000 -a 1000 --sortRegions descend  --sortUsingSamples 1 --missingDataAsZero --skipZeros \
    --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed \
    --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF" \
    -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsTKO.bdgdiff_3cPeakCenter.gz \
    --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsTKO.bdgdiff_3cPeakCenter.matrix \
    --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsTKO.bdgdiff_3cPeakCenter.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsTKO.bdgdiff_3cPeakCenter.gz \
    -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsTKO.bdgdiff_3cPeakCenter.plot.pdf \
    --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 \
    --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" \
    --regionsLabel "ROI_WT (6436)" "common (50966)" "ROI_SKO (872)" --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF"

## rep1/2/3 individual bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/10_plotHeatMap_for_all_use_mESC_WTvsTKO_CTCF_summits
computeMatrix reference-point -p "max" -R ROI_WT.CTCF.summit.sorted.bed RIB_WT_TKO.CTCF.summit.sorted.bed ROI_TKO.CTCF.summit.sorted.bed \
    -S mESC_WT_IP_CTCF_rep1.TvsC.bigwig mESC_WT_IP_CTCF_rep2.TvsC.bigwig mESC_WT_IP_CTCF_rep3.TvsC.bigwig \
    mESC_SKO_IP_CTCF_rep1.TvsC.bigwig mESC_SKO_IP_CTCF_rep2.TvsC.bigwig mESC_SKO_IP_CTCF_rep3.TvsC.bigwig \
    mESC_TKO_IP_CTCF_rep1.TvsC.bigwig mESC_TKO_IP_CTCF_rep2.TvsC.bigwig mESC_TKO_IP_CTCF_rep3.TvsC.bigwig \
    --referencePoint center -b 1000 -a 1000 --sortRegions descend --sortUsingSamples 1 2 3 --missingDataAsZero --skipZeros \
    --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed \
    --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" \
    "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" \
    "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3" \
    -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsTKO.bdgdiff_3cPeakCenter.gz \
    --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsTKO.bdgdiff_3cPeakCenter.matrix \
    --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsTKO.bdgdiff_3cPeakCenter.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsTKO.bdgdiff_3cPeakCenter.gz \
    -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsTKO.bdgdiff_3cPeakCenter.plot.pdf \
    --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 \
    --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" \
    --regionsLabel "ROI_WT (6436)" "common (50966)" "ROI_SKO (872)" \
    --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" \
    "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" \
    "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3"

# 最后为pdf图生成对应的png图，方便预览
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/10_plotHeatMap_for_all_use_mESC_WTvsTKO_CTCF_summits
for i in `ls -1 bamCompare_matrix_of*.pdf`;do echo $i; pdftoppm -png $i ${i%.pdf};done

```



## SKO_vs_TKO

### 创建 TvsC_bigwig 与 两种方法各自三类 peak center 文件连接
```bash
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
mkdir 10_plotHeatMap_for_all_use_mESC_SKOvsTKO_CTCF_summits
cd 10_plotHeatMap_for_all_use_mESC_SKOvsTKO_CTCF_summits
# 创建bigwig文件连接
## rep1/2/3 individual
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF_rep1.TvsC.bigwig mESC_WT_IP_CTCF_rep1.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF_rep2.TvsC.bigwig mESC_WT_IP_CTCF_rep2.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF_rep3.TvsC.bigwig mESC_WT_IP_CTCF_rep3.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF_rep1.TvsC.bigwig mESC_SKO_IP_CTCF_rep1.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF_rep2.TvsC.bigwig mESC_SKO_IP_CTCF_rep2.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF_rep3.TvsC.bigwig mESC_SKO_IP_CTCF_rep3.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF_rep1.TvsC.bigwig mESC_TKO_IP_CTCF_rep1.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF_rep2.TvsC.bigwig mESC_TKO_IP_CTCF_rep2.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF_rep3.TvsC.bigwig mESC_TKO_IP_CTCF_rep3.TvsC.bigwig
## rep123 merged
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig
ls -1 *rep[123].TvsC.bigwig
ls -1 *rep123_merged.TvsC.bigwig

#  创建三类 peak center 文件连接 (bedtools intersect)
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/mESC_SKO_TKO_intersect.SKO_only.region.centers.bed mESC_SKO_TKO_intersect.SKO_only.region.centers.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/mESC_SKO_TKO_intersect.common.region.centers.bed mESC_SKO_TKO_intersect.common.region.centers.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/mESC_SKO_TKO_intersect.TKO_only.region.centers.bed mESC_SKO_TKO_intersect.TKO_only.region.centers.bed
wc -l *.region.centers.bed
#  6165 mESC_SKO_TKO_intersect.SKO_only.region.centers.bed
# 53357 mESC_SKO_TKO_intersect.common.region.centers.bed
#  8996 mESC_SKO_TKO_intersect.TKO_only.region.centers.bed

#  创建三类 peak center 文件连接 (bdgdiff)
ln -s ../7_macs2_peak_calling_3_for_SKO_vs_TKO_bdgdiff/ROI_SKO.CTCF.summit.sorted.bed ROI_SKO.CTCF.summit.sorted.bed
ln -s ../7_macs2_peak_calling_3_for_SKO_vs_TKO_bdgdiff/RIB_SKO_TKO.CTCF.summit.sorted.bed RIB_SKO_TKO.CTCF.summit.sorted.bed
ln -s ../7_macs2_peak_calling_3_for_SKO_vs_TKO_bdgdiff/ROI_TKO.CTCF.summit.sorted.bed ROI_TKO.CTCF.summit.sorted.bed
wc -l *.summit.sorted.bed
#  3705 ROI_SKO.CTCF.summit.sorted.bed
# 47428 RIB_SKO_TKO.CTCF.summit.sorted.bed
#   922 ROI_TKO.CTCF.summit.sorted.bed

```

### computeMatrix and plotHeatmap

#### use bedtools intersect 三类 peak center

```bash
# use bedtools intersect 三类 peak center
## rep123 merged bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/10_plotHeatMap_for_all_use_mESC_SKOvsTKO_CTCF_summits
computeMatrix reference-point -p "max" -R mESC_SKO_TKO_intersect.SKO_only.region.centers.bed mESC_SKO_TKO_intersect.common.region.centers.bed mESC_SKO_TKO_intersect.TKO_only.region.centers.bed \
    -S mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig \
    --referencePoint center -b 1000 -a 1000 --sortRegions descend  --sortUsingSamples 1 --missingDataAsZero --skipZeros \
    --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed \
    --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF" \
    -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.SKOvsTKO.intersect_3cPeakCenter.gz \
    --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.SKOvsTKO.intersect_3cPeakCenter.matrix \
    --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.SKOvsTKO.intersect_3cPeakCenter.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.SKOvsTKO.intersect_3cPeakCenter.gz \
    -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.SKOvsTKO.intersect_3cPeakCenter.plot.pdf \
    --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 \
    --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" \
    --regionsLabel "SKO_only (6165)" "common (53357)" "TKO_only (8996)" --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF"

## rep1/2/3 individual bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/10_plotHeatMap_for_all_use_mESC_SKOvsTKO_CTCF_summits
computeMatrix reference-point -p "max" -R mESC_SKO_TKO_intersect.SKO_only.region.centers.bed mESC_SKO_TKO_intersect.common.region.centers.bed mESC_SKO_TKO_intersect.TKO_only.region.centers.bed \
    -S mESC_WT_IP_CTCF_rep1.TvsC.bigwig mESC_WT_IP_CTCF_rep2.TvsC.bigwig mESC_WT_IP_CTCF_rep3.TvsC.bigwig \
    mESC_SKO_IP_CTCF_rep1.TvsC.bigwig mESC_SKO_IP_CTCF_rep2.TvsC.bigwig mESC_SKO_IP_CTCF_rep3.TvsC.bigwig \
    mESC_TKO_IP_CTCF_rep1.TvsC.bigwig mESC_TKO_IP_CTCF_rep2.TvsC.bigwig mESC_TKO_IP_CTCF_rep3.TvsC.bigwig \
    --referencePoint center -b 1000 -a 1000 --sortRegions descend --sortUsingSamples 1 2 3 --missingDataAsZero --skipZeros \
    --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed \
    --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" \
    "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" \
    "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3" \
    -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.SKOvsTKO.intersect_3cPeakCenter.gz \
    --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.SKOvsTKO.intersect_3cPeakCenter.matrix \
    --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.SKOvsTKO.intersect_3cPeakCenter.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.SKOvsTKO.intersect_3cPeakCenter.gz \
    -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.SKOvsTKO.intersect_3cPeakCenter.plot.pdf \
    --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 \
    --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" \
    --regionsLabel "SKO_only (6165)" "common (53357)" "TKO_only (8996)" \
    --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" \
    "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" \
    "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3"

```

#### use bdgdiff 三类 peak center

```bash
# use bdgdiff 三类 peak center
## rep123 merged bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/10_plotHeatMap_for_all_use_mESC_SKOvsTKO_CTCF_summits
computeMatrix reference-point -p "max" -R ROI_SKO.CTCF.summit.sorted.bed RIB_SKO_TKO.CTCF.summit.sorted.bed ROI_TKO.CTCF.summit.sorted.bed \
    -S mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig \
    --referencePoint center -b 1000 -a 1000 --sortRegions descend  --sortUsingSamples 1 --missingDataAsZero --skipZeros \
    --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed \
    --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF" \
    -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.SKOvsTKO.bdgdiff_3cPeakCenter.gz \
    --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.SKOvsTKO.bdgdiff_3cPeakCenter.matrix \
    --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.SKOvsTKO.bdgdiff_3cPeakCenter.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.SKOvsTKO.bdgdiff_3cPeakCenter.gz \
    -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.SKOvsTKO.bdgdiff_3cPeakCenter.plot.pdf \
    --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 \
    --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" \
    --regionsLabel "ROI_SKO (3705)" "common (47428)" "ROI_TKO (922)" --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF"

## rep1/2/3 individual bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/10_plotHeatMap_for_all_use_mESC_SKOvsTKO_CTCF_summits
computeMatrix reference-point -p "max" -R ROI_SKO.CTCF.summit.sorted.bed RIB_SKO_TKO.CTCF.summit.sorted.bed ROI_TKO.CTCF.summit.sorted.bed \
    -S mESC_WT_IP_CTCF_rep1.TvsC.bigwig mESC_WT_IP_CTCF_rep2.TvsC.bigwig mESC_WT_IP_CTCF_rep3.TvsC.bigwig \
    mESC_SKO_IP_CTCF_rep1.TvsC.bigwig mESC_SKO_IP_CTCF_rep2.TvsC.bigwig mESC_SKO_IP_CTCF_rep3.TvsC.bigwig \
    mESC_TKO_IP_CTCF_rep1.TvsC.bigwig mESC_TKO_IP_CTCF_rep2.TvsC.bigwig mESC_TKO_IP_CTCF_rep3.TvsC.bigwig \
    --referencePoint center -b 1000 -a 1000 --sortRegions descend --sortUsingSamples 1 2 3 --missingDataAsZero --skipZeros \
    --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed \
    --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" \
    "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" \
    "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3" \
    -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.SKOvsTKO.bdgdiff_3cPeakCenter.gz \
    --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.SKOvsTKO.bdgdiff_3cPeakCenter.matrix \
    --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.SKOvsTKO.bdgdiff_3cPeakCenter.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.SKOvsTKO.bdgdiff_3cPeakCenter.gz \
    -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.SKOvsTKO.bdgdiff_3cPeakCenter.plot.pdf \
    --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 \
    --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" \
    --regionsLabel "ROI_SKO (3705)" "common (47428)" "ROI_TKO (922)" \
    --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" \
    "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" \
    "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3"

# 最后为pdf图生成对应的png图，方便预览
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/10_plotHeatMap_for_all_use_mESC_SKOvsTKO_CTCF_summits
for i in `ls -1 bamCompare_matrix_of*.pdf`;do echo $i; pdftoppm -png $i ${i%.pdf};done

```

# 三组实验 intervene 结果的7类 region centers  周围信号富集热图

## 创建 bigwig 文件连接并构造7类 region centers 

```bash
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
mkdir 11_plotHeatMap_for_all_use_mESC_WTvsSKOvsTKO_CTCF_summits
cd 11_plotHeatMap_for_all_use_mESC_WTvsSKOvsTKO_CTCF_summits
# 创建bigwig文件连接
## rep1/2/3 individual
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF_rep1.TvsC.bigwig mESC_WT_IP_CTCF_rep1.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF_rep2.TvsC.bigwig mESC_WT_IP_CTCF_rep2.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF_rep3.TvsC.bigwig mESC_WT_IP_CTCF_rep3.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF_rep1.TvsC.bigwig mESC_SKO_IP_CTCF_rep1.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF_rep2.TvsC.bigwig mESC_SKO_IP_CTCF_rep2.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF_rep3.TvsC.bigwig mESC_SKO_IP_CTCF_rep3.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF_rep1.TvsC.bigwig mESC_TKO_IP_CTCF_rep1.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF_rep2.TvsC.bigwig mESC_TKO_IP_CTCF_rep2.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF_rep3.TvsC.bigwig mESC_TKO_IP_CTCF_rep3.TvsC.bigwig
## rep123 merged
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig
ln -s ../8_bamCompare_TvsC_bigwig_for_computeMatrix_and_plotHeatmap/mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig
ls -1 *rep[123].TvsC.bigwig
ls -1 *rep123_merged.TvsC.bigwig

#  创建7类 region centers 文件 (peaks intervene results)
#  peaks intervene 7 类区间结果在如下目录
#  ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/0_mESC_WT_SKO_TKO_IP_CTCF_peaks_venn/sets/
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/0_mESC_WT_SKO_TKO_IP_CTCF_peaks_venn/sets/001_mESC_TKO_IP_CTCF_peaks.bed peaks.001.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/0_mESC_WT_SKO_TKO_IP_CTCF_peaks_venn/sets/010_mESC_SKO_IP_CTCF_peaks.bed peaks.010.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/0_mESC_WT_SKO_TKO_IP_CTCF_peaks_venn/sets/011_mESC_SKO_IP_CTCF_peaks_mESC_TKO_IP_CTCF_peaks.bed peaks.011.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/0_mESC_WT_SKO_TKO_IP_CTCF_peaks_venn/sets/100_mESC_WT_IP_CTCF_peaks.bed peaks.100.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/0_mESC_WT_SKO_TKO_IP_CTCF_peaks_venn/sets/101_mESC_WT_IP_CTCF_peaks_mESC_TKO_IP_CTCF_peaks.bed peaks.101.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/0_mESC_WT_SKO_TKO_IP_CTCF_peaks_venn/sets/110_mESC_WT_IP_CTCF_peaks_mESC_SKO_IP_CTCF_peaks.bed peaks.110.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/0_mESC_WT_SKO_TKO_IP_CTCF_peaks_venn/sets/111_mESC_WT_IP_CTCF_peaks_mESC_SKO_IP_CTCF_peaks_mESC_TKO_IP_CTCF_peaks.bed peaks.111.bed
wc -l peaks*.bed
# 构造7类区间中心点 bed
awk 'BEGIN{"\t";OFS="\t"}{center=int(($2+$3)/2);print $1,center,center+1}' peaks.001.bed |sort -k1,1 -k2,2n |uniq >PKIc.001.sorted.bed
awk 'BEGIN{"\t";OFS="\t"}{center=int(($2+$3)/2);print $1,center,center+1}' peaks.010.bed |sort -k1,1 -k2,2n |uniq >PKIc.010.sorted.bed
awk 'BEGIN{"\t";OFS="\t"}{center=int(($2+$3)/2);print $1,center,center+1}' peaks.011.bed |sort -k1,1 -k2,2n |uniq >PKIc.011.sorted.bed
awk 'BEGIN{"\t";OFS="\t"}{center=int(($2+$3)/2);print $1,center,center+1}' peaks.100.bed |sort -k1,1 -k2,2n |uniq >PKIc.100.sorted.bed
awk 'BEGIN{"\t";OFS="\t"}{center=int(($2+$3)/2);print $1,center,center+1}' peaks.101.bed |sort -k1,1 -k2,2n |uniq >PKIc.101.sorted.bed
awk 'BEGIN{"\t";OFS="\t"}{center=int(($2+$3)/2);print $1,center,center+1}' peaks.110.bed |sort -k1,1 -k2,2n |uniq >PKIc.110.sorted.bed
awk 'BEGIN{"\t";OFS="\t"}{center=int(($2+$3)/2);print $1,center,center+1}' peaks.111.bed |sort -k1,1 -k2,2n |uniq >PKIc.111.sorted.bed
wc -l PKIc*.sorted.bed
#  9934 PKIc.100.sorted.bed
#  2610 PKIc.010.sorted.bed
#  3359 PKIc.001.sorted.bed
# 51946 PKIc.111.sorted.bed
#  3558 PKIc.110.sorted.bed
#  1301 PKIc.011.sorted.bed
#  5637 PKIc.101.sorted.bed

#  创建7类 region centers 文件 (peak summits eB150 intervene results)
#  peak summits eB150 intervene 7 类区间结果在如下目录
#  ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/0_mESC_WT_SKO_TKO_IP_CTCF_summits_eB150_venn/sets
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/0_mESC_WT_SKO_TKO_IP_CTCF_summits_eB150_venn/sets/001_mESC_TKO_IP_CTCF_seB150.bed eB150.001.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/0_mESC_WT_SKO_TKO_IP_CTCF_summits_eB150_venn/sets/010_mESC_SKO_IP_CTCF_seB150.bed eB150.010.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/0_mESC_WT_SKO_TKO_IP_CTCF_summits_eB150_venn/sets/011_mESC_SKO_IP_CTCF_seB150_mESC_TKO_IP_CTCF_seB150.bed eB150.011.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/0_mESC_WT_SKO_TKO_IP_CTCF_summits_eB150_venn/sets/100_mESC_WT_IP_CTCF_seB150.bed eB150.100.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/0_mESC_WT_SKO_TKO_IP_CTCF_summits_eB150_venn/sets/101_mESC_WT_IP_CTCF_seB150_mESC_TKO_IP_CTCF_seB150.bed eB150.101.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/0_mESC_WT_SKO_TKO_IP_CTCF_summits_eB150_venn/sets/110_mESC_WT_IP_CTCF_seB150_mESC_SKO_IP_CTCF_seB150.bed eB150.110.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/0_mESC_WT_SKO_TKO_IP_CTCF_summits_eB150_venn/sets/111_mESC_WT_IP_CTCF_seB150_mESC_SKO_IP_CTCF_seB150_mESC_TKO_IP_CTCF_seB150.bed eB150.111.bed
wc -l eB150*.bed
# 构造7类区间中心点 bed
awk 'BEGIN{"\t";OFS="\t"}{center=int(($2+$3)/2);print $1,center,center+1}' eB150.001.bed |sort -k1,1 -k2,2n |uniq >PSEIc.001.sorted.bed
awk 'BEGIN{"\t";OFS="\t"}{center=int(($2+$3)/2);print $1,center,center+1}' eB150.010.bed |sort -k1,1 -k2,2n |uniq >PSEIc.010.sorted.bed
awk 'BEGIN{"\t";OFS="\t"}{center=int(($2+$3)/2);print $1,center,center+1}' eB150.011.bed |sort -k1,1 -k2,2n |uniq >PSEIc.011.sorted.bed
awk 'BEGIN{"\t";OFS="\t"}{center=int(($2+$3)/2);print $1,center,center+1}' eB150.100.bed |sort -k1,1 -k2,2n |uniq >PSEIc.100.sorted.bed
awk 'BEGIN{"\t";OFS="\t"}{center=int(($2+$3)/2);print $1,center,center+1}' eB150.101.bed |sort -k1,1 -k2,2n |uniq >PSEIc.101.sorted.bed
awk 'BEGIN{"\t";OFS="\t"}{center=int(($2+$3)/2);print $1,center,center+1}' eB150.110.bed |sort -k1,1 -k2,2n |uniq >PSEIc.110.sorted.bed
awk 'BEGIN{"\t";OFS="\t"}{center=int(($2+$3)/2);print $1,center,center+1}' eB150.111.bed |sort -k1,1 -k2,2n |uniq >PSEIc.111.sorted.bed
wc -l PSEIc*.sorted.bed
# 10304 PSEIc.100.sorted.bed
#  2828 PSEIc.010.sorted.bed
#  3645 PSEIc.001.sorted.bed
# 53235 PSEIc.111.sorted.bed
#  3697 PSEIc.110.sorted.bed
#  1373 PSEIc.011.sorted.bed
#  5889 PSEIc.101.sorted.bed

```

## 计算并绘制热图

### use 7cRC from peaks intervene results

```bash
# use peaks intervene results 7类 region centers
## rep123 merged bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/11_plotHeatMap_for_all_use_mESC_WTvsSKOvsTKO_CTCF_summits
computeMatrix reference-point -p "max" -R PKIc.100.sorted.bed PKIc.010.sorted.bed PKIc.001.sorted.bed PKIc.111.sorted.bed PKIc.110.sorted.bed PKIc.011.sorted.bed PKIc.101.sorted.bed \
    -S mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig \
    --referencePoint center -b 1000 -a 1000 --sortRegions descend  --sortUsingSamples 1 --missingDataAsZero --skipZeros \
    --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed \
    --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF" \
    -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKOvsTKO.peaksIntervene_7cRegionCenters.gz \
    --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKOvsTKO.peaksIntervene_7cRegionCenters.matrix \
    --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKOvsTKO.peaksIntervene_7cRegionCenters.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKOvsTKO.peaksIntervene_7cRegionCenters.gz \
    -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKOvsTKO.peaksIntervene_7cRegionCenters.plot.pdf \
    --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 \
    --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" \
    --regionsLabel "WOO (9934)" "OSO (2610)" "OOT (3359)" "WST (51946)" "WSO (3558)" "OST (1301)" "WOT (5637)" \
    --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF"

## rep1/2/3 individual bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/11_plotHeatMap_for_all_use_mESC_WTvsSKOvsTKO_CTCF_summits
computeMatrix reference-point -p "max" -R PKIc.100.sorted.bed PKIc.010.sorted.bed PKIc.001.sorted.bed PKIc.111.sorted.bed PKIc.110.sorted.bed PKIc.011.sorted.bed PKIc.101.sorted.bed \
    -S mESC_WT_IP_CTCF_rep1.TvsC.bigwig mESC_WT_IP_CTCF_rep2.TvsC.bigwig mESC_WT_IP_CTCF_rep3.TvsC.bigwig \
    mESC_SKO_IP_CTCF_rep1.TvsC.bigwig mESC_SKO_IP_CTCF_rep2.TvsC.bigwig mESC_SKO_IP_CTCF_rep3.TvsC.bigwig \
    mESC_TKO_IP_CTCF_rep1.TvsC.bigwig mESC_TKO_IP_CTCF_rep2.TvsC.bigwig mESC_TKO_IP_CTCF_rep3.TvsC.bigwig \
    --referencePoint center -b 1000 -a 1000 --sortRegions descend --sortUsingSamples 1 2 3 --missingDataAsZero --skipZeros \
    --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed \
    --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" \
    "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" \
    "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3" \
    -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKOvsTKO.peaksIntervene_7cRegionCenters.gz \
    --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKOvsTKO.peaksIntervene_7cRegionCenters.matrix \
    --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKOvsTKO.peaksIntervene_7cRegionCenters.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKOvsTKO.peaksIntervene_7cRegionCenters.gz \
    -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKOvsTKO.peaksIntervene_7cRegionCenters.plot.pdf \
    --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 \
    --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" \
    --regionsLabel "WOO (9934)" "OSO (2610)" "OOT (3359)" "WST (51946)" "WSO (3558)" "OST (1301)" "WOT (5637)" \
    --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" \
    "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" \
    "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3"

```

#### use 6cRC from peaks intervene results (remove WST common)

```bash
# use peaks intervene results 6类 region centers  (remove WST common)
## rep123 merged bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/11_plotHeatMap_for_all_use_mESC_WTvsSKOvsTKO_CTCF_summits
computeMatrix reference-point -p "max" -R PKIc.100.sorted.bed PKIc.010.sorted.bed PKIc.001.sorted.bed PKIc.110.sorted.bed PKIc.011.sorted.bed PKIc.101.sorted.bed \
    -S mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig \
    --referencePoint center -b 1000 -a 1000 --sortRegions descend  --sortUsingSamples 1 --missingDataAsZero --skipZeros \
    --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed \
    --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF" \
    -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKOvsTKO.peaksIntervene_6cRegionCenters.gz \
    --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKOvsTKO.peaksIntervene_6cRegionCenters.matrix \
    --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKOvsTKO.peaksIntervene_6cRegionCenters.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKOvsTKO.peaksIntervene_6cRegionCenters.gz \
    -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKOvsTKO.peaksIntervene_6cRegionCenters.plot.pdf \
    --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 \
    --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" \
    --regionsLabel "WOO (9934)" "OSO (2610)" "OOT (3359)" "WSO (3558)" "OST (1301)" "WOT (5637)" \
    --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF"

## rep1/2/3 individual bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/11_plotHeatMap_for_all_use_mESC_WTvsSKOvsTKO_CTCF_summits
computeMatrix reference-point -p "max" -R PKIc.100.sorted.bed PKIc.010.sorted.bed PKIc.001.sorted.bed PKIc.110.sorted.bed PKIc.011.sorted.bed PKIc.101.sorted.bed \
    -S mESC_WT_IP_CTCF_rep1.TvsC.bigwig mESC_WT_IP_CTCF_rep2.TvsC.bigwig mESC_WT_IP_CTCF_rep3.TvsC.bigwig \
    mESC_SKO_IP_CTCF_rep1.TvsC.bigwig mESC_SKO_IP_CTCF_rep2.TvsC.bigwig mESC_SKO_IP_CTCF_rep3.TvsC.bigwig \
    mESC_TKO_IP_CTCF_rep1.TvsC.bigwig mESC_TKO_IP_CTCF_rep2.TvsC.bigwig mESC_TKO_IP_CTCF_rep3.TvsC.bigwig \
    --referencePoint center -b 1000 -a 1000 --sortRegions descend --sortUsingSamples 1 2 3 --missingDataAsZero --skipZeros \
    --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed \
    --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" \
    "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" \
    "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3" \
    -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKOvsTKO.peaksIntervene_6cRegionCenters.gz \
    --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKOvsTKO.peaksIntervene_6cRegionCenters.matrix \
    --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKOvsTKO.peaksIntervene_6cRegionCenters.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKOvsTKO.peaksIntervene_6cRegionCenters.gz \
    -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKOvsTKO.peaksIntervene_6cRegionCenters.plot.pdf \
    --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 \
    --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" \
    --regionsLabel "WOO (9934)" "OSO (2610)" "OOT (3359)" "WSO (3558)" "OST (1301)" "WOT (5637)" \
    --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" \
    "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" \
    "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3"

```



### use 7cRC from peak summits eB150 intervene results

```bash
# use peak summits eB150 intervene results 7类 region centers
## rep123 merged bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/11_plotHeatMap_for_all_use_mESC_WTvsSKOvsTKO_CTCF_summits
computeMatrix reference-point -p "max" -R PSEIc.100.sorted.bed PSEIc.010.sorted.bed PSEIc.001.sorted.bed PSEIc.111.sorted.bed PSEIc.110.sorted.bed PSEIc.011.sorted.bed PSEIc.101.sorted.bed \
    -S mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig \
    --referencePoint center -b 1000 -a 1000 --sortRegions descend  --sortUsingSamples 1 --missingDataAsZero --skipZeros \
    --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed \
    --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF" \
    -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKOvsTKO.PSeB150Intervene_7cRegionCenters.gz \
    --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKOvsTKO.PSeB150Intervene_7cRegionCenters.matrix \
    --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKOvsTKO.PSeB150Intervene_7cRegionCenters.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKOvsTKO.PSeB150Intervene_7cRegionCenters.gz \
    -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKOvsTKO.PSeB150Intervene_7cRegionCenters.plot.pdf \
    --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 \
    --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" \
    --regionsLabel "WOO (10304)" "OSO (2828)" "OOT (3645)" "WST (53235)" "WSO (3697)" "OST (1373)" "WOT (5889)" \
    --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF"

## rep1/2/3 individual bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/11_plotHeatMap_for_all_use_mESC_WTvsSKOvsTKO_CTCF_summits
computeMatrix reference-point -p "max" -R PSEIc.100.sorted.bed PSEIc.010.sorted.bed PSEIc.001.sorted.bed PSEIc.111.sorted.bed PSEIc.110.sorted.bed PSEIc.011.sorted.bed PSEIc.101.sorted.bed \
    -S mESC_WT_IP_CTCF_rep1.TvsC.bigwig mESC_WT_IP_CTCF_rep2.TvsC.bigwig mESC_WT_IP_CTCF_rep3.TvsC.bigwig \
    mESC_SKO_IP_CTCF_rep1.TvsC.bigwig mESC_SKO_IP_CTCF_rep2.TvsC.bigwig mESC_SKO_IP_CTCF_rep3.TvsC.bigwig \
    mESC_TKO_IP_CTCF_rep1.TvsC.bigwig mESC_TKO_IP_CTCF_rep2.TvsC.bigwig mESC_TKO_IP_CTCF_rep3.TvsC.bigwig \
    --referencePoint center -b 1000 -a 1000 --sortRegions descend --sortUsingSamples 1 2 3 --missingDataAsZero --skipZeros \
    --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed \
    --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" \
    "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" \
    "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3" \
    -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKOvsTKO.PSeB150Intervene_7cRegionCenters.gz \
    --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKOvsTKO.PSeB150Intervene_7cRegionCenters.matrix \
    --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKOvsTKO.PSeB150Intervene_7cRegionCenters.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKOvsTKO.PSeB150Intervene_7cRegionCenters.gz \
    -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKOvsTKO.PSeB150Intervene_7cRegionCenters.plot.pdf \
    --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 \
    --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" \
    --regionsLabel "WOO (10304)" "OSO (2828)" "OOT (3645)" "WST (53235)" "WSO (3697)" "OST (1373)" "WOT (5889)" \
    --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" \
    "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" \
    "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3"

```

#### use 6cRC from peak summits eB150 intervene results  (remove WST common)

```bash
# use peak summits eB150 intervene results 6类 region centers  (remove WST common)
## rep123 merged bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/11_plotHeatMap_for_all_use_mESC_WTvsSKOvsTKO_CTCF_summits
computeMatrix reference-point -p "max" -R PSEIc.100.sorted.bed PSEIc.010.sorted.bed PSEIc.001.sorted.bed PSEIc.110.sorted.bed PSEIc.011.sorted.bed PSEIc.101.sorted.bed \
    -S mESC_WT_IP_CTCF.rep123_merged.TvsC.bigwig mESC_SKO_IP_CTCF.rep123_merged.TvsC.bigwig mESC_TKO_IP_CTCF.rep123_merged.TvsC.bigwig \
    --referencePoint center -b 1000 -a 1000 --sortRegions descend  --sortUsingSamples 1 --missingDataAsZero --skipZeros \
    --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed \
    --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF" \
    -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKOvsTKO.PSeB150Intervene_6cRegionCenters.gz \
    --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKOvsTKO.PSeB150Intervene_6cRegionCenters.matrix \
    --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKOvsTKO.PSeB150Intervene_6cRegionCenters.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKOvsTKO.PSeB150Intervene_6cRegionCenters.gz \
    -o bamCompare_matrix_of_WT_SKO_TKO_rep123_merged_for_plotHeatMap.WTvsSKOvsTKO.PSeB150Intervene_6cRegionCenters.plot.pdf \
    --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 \
    --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" \
    --regionsLabel "WOO (10304)" "OSO (2828)" "OOT (3645)" "WSO (3697)" "OST (1373)" "WOT (5889)" \
    --samplesLabel "mESC_WT_IP_CTCF" "mESC_SKO_IP_CTCF" "mESC_TKO_IP_CTCF"

## rep1/2/3 individual bigwig
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/11_plotHeatMap_for_all_use_mESC_WTvsSKOvsTKO_CTCF_summits
computeMatrix reference-point -p "max" -R PSEIc.100.sorted.bed PSEIc.010.sorted.bed PSEIc.001.sorted.bed PSEIc.110.sorted.bed PSEIc.011.sorted.bed PSEIc.101.sorted.bed \
    -S mESC_WT_IP_CTCF_rep1.TvsC.bigwig mESC_WT_IP_CTCF_rep2.TvsC.bigwig mESC_WT_IP_CTCF_rep3.TvsC.bigwig \
    mESC_SKO_IP_CTCF_rep1.TvsC.bigwig mESC_SKO_IP_CTCF_rep2.TvsC.bigwig mESC_SKO_IP_CTCF_rep3.TvsC.bigwig \
    mESC_TKO_IP_CTCF_rep1.TvsC.bigwig mESC_TKO_IP_CTCF_rep2.TvsC.bigwig mESC_TKO_IP_CTCF_rep3.TvsC.bigwig \
    --referencePoint center -b 1000 -a 1000 --sortRegions descend --sortUsingSamples 1 2 3 --missingDataAsZero --skipZeros \
    --blackListFileName /home/sunwenju/3_genomeData/mouse/mm10/encode_mm10_black_list.sorted.bed \
    --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" \
    "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" \
    "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3" \
    -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKOvsTKO.PSeB150Intervene_6cRegionCenters.gz \
    --outFileNameMatrix bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKOvsTKO.PSeB150Intervene_6cRegionCenters.matrix \
    --outFileSortedRegions bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKOvsTKO.PSeB150Intervene_6cRegionCenters.sorted_regions.bed
plotHeatmap -m bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKOvsTKO.PSeB150Intervene_6cRegionCenters.gz \
    -o bamCompare_matrix_of_WT_SKO_TKO_all_repeats_for_plotHeatMap.WTvsSKOvsTKO.PSeB150Intervene_6cRegionCenters.plot.pdf \
    --colorMap 'coolwarm' --zMin "auto" --zMax "auto" --heatmapHeight 18 --heatmapWidth 5 \
    --whatToShow "plot, heatmap and colorbar" --refPointLabel "0" --xAxisLabel "distance to CTCF_summits" \
    --regionsLabel "WOO (10304)" "OSO (2828)" "OOT (3645)" "WSO (3697)" "OST (1373)" "WOT (5889)" \
    --samplesLabel "mESC_WT_IP_CTCF_rep1" "mESC_WT_IP_CTCF_rep2" "mESC_WT_IP_CTCF_rep3" \
    "mESC_SKO_IP_CTCF_rep1" "mESC_SKO_IP_CTCF_rep2" "mESC_SKO_IP_CTCF_rep3" \
    "mESC_TKO_IP_CTCF_rep1" "mESC_TKO_IP_CTCF_rep2" "mESC_TKO_IP_CTCF_rep3"

```

```bash
# 最后为pdf图生成对应的png图，方便预览
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/11_plotHeatMap_for_all_use_mESC_WTvsSKOvsTKO_CTCF_summits
for i in `ls -1 bamCompare_matrix_of*.pdf`;do echo $i; pdftoppm -png $i ${i%.pdf};done

```

# 绘制WST三组peak两两intersect生成的三类peak summit f2k的甲基化分布图

> mESC WT / SKO / TKO 三种细胞的 WGBS数据已经分析生成三种细胞全基因组CG位点的甲基化状态文件 bed6+4 （与之前绘制甲基化分布的输入格式一致，可直接作为 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py 的甲基化位点输入文件），文件在如下目录
>
> /data1/sunwenju/projects_on_940xa/2024031102_WN_mESC_BS_Seq_add/2_all_samples_genome_wide_CpG_sites_info

```bash
cd /data1/sunwenju/projects_on_940xa/2024031102_WN_mESC_BS_Seq_add/2_all_samples_genome_wide_CpG_sites_info/
ls -1
wc -l mESC_*all_CpG_result.txt
# 40219388 mESC_SKO_BS_Seq_rep123.all_CpG_result.txt
# 40309292 mESC_TKO_BS_Seq_rep123.all_CpG_result.txt
# 40129991 mESC_WT_BS_Seq_rep123.all_CpG_result.txt

```

## WT_vs_SKO

### 创建三类 peak center 和甲基化位点文件链接

```bash
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
mkdir 12_mESC_WTvsSKO_3cCTCF_summits_f2k_mCpG_stat_plot
cd 12_mESC_WTvsSKO_3cCTCF_summits_f2k_mCpG_stat_plot
#  创建三类 peak center 文件连接 (bedtools intersect)
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/mESC_WT_SKO_intersect.WT_only.region.centers.bed mESC_WT_SKO_intersect.WT_only.region.centers.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/mESC_WT_SKO_intersect.common.region.centers.bed mESC_WT_SKO_intersect.common.region.centers.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/mESC_WT_SKO_intersect.SKO_only.region.centers.bed mESC_WT_SKO_intersect.SKO_only.region.centers.bed
wc -l *.region.centers.bed
# 15571 mESC_WT_SKO_intersect.WT_only.region.centers.bed
# 55629 mESC_WT_SKO_intersect.common.region.centers.bed
#  3911 mESC_WT_SKO_intersect.SKO_only.region.centers.bed

# 创建CpG甲基化位点结果文件连接
ln -s /data1/sunwenju/projects_on_940xa/2024031102_WN_mESC_BS_Seq_add/2_all_samples_genome_wide_CpG_sites_info/mESC_SKO_BS_Seq_rep123.all_CpG_result.txt mESC_SKO_BS_Seq_rep123.all_CpG_result.txt
ln -s /data1/sunwenju/projects_on_940xa/2024031102_WN_mESC_BS_Seq_add/2_all_samples_genome_wide_CpG_sites_info/mESC_TKO_BS_Seq_rep123.all_CpG_result.txt mESC_TKO_BS_Seq_rep123.all_CpG_result.txt
ln -s /data1/sunwenju/projects_on_940xa/2024031102_WN_mESC_BS_Seq_add/2_all_samples_genome_wide_CpG_sites_info/mESC_WT_BS_Seq_rep123.all_CpG_result.txt mESC_WT_BS_Seq_rep123.all_CpG_result.txt
ls -1 mESC_*all_CpG_result.txt

```

### 统计三类 CTCF peak center f2k mCpG 分布并绘图

```bash
# 统计三类 CTCF peak centers 上下游甲基化位点的分布
# 总结之前的经验，可以全部单独运行
## mESC_WTvSKO.WT_only_PKc.f2k  mCpG stat
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_WT_SKO_intersect.WT_only.region.centers.bed -b mESC_WT_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WTvSKO.WT_only_PKc.f2k.WT_mCpG.stat.txt >mESC_WTvSKO.WT_only_PKc.f2k.WT_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_WT_SKO_intersect.WT_only.region.centers.bed -b mESC_SKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WTvSKO.WT_only_PKc.f2k.SKO_mCpG.stat.txt >mESC_WTvSKO.WT_only_PKc.f2k.SKO_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_WT_SKO_intersect.WT_only.region.centers.bed -b mESC_TKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WTvSKO.WT_only_PKc.f2k.TKO_mCpG.stat.txt >mESC_WTvSKO.WT_only_PKc.f2k.TKO_mCpG.stat.runlog 2>&1 &
tail -f mESC_WTvSKO.WT_only_PKc.f2k.WT_mCpG.stat.runlog
tail -f mESC_WTvSKO.WT_only_PKc.f2k.SKO_mCpG.stat.runlog
tail -f mESC_WTvSKO.WT_only_PKc.f2k.TKO_mCpG.stat.runlog
## mESC_WTvSKO.common_PKc.f2k  mCpG stat
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_WT_SKO_intersect.common.region.centers.bed -b mESC_WT_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WTvSKO.common_PKc.f2k.WT_mCpG.stat.txt >mESC_WTvSKO.common_PKc.f2k.WT_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_WT_SKO_intersect.common.region.centers.bed -b mESC_SKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WTvSKO.common_PKc.f2k.SKO_mCpG.stat.txt >mESC_WTvSKO.common_PKc.f2k.SKO_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_WT_SKO_intersect.common.region.centers.bed -b mESC_TKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WTvSKO.common_PKc.f2k.TKO_mCpG.stat.txt >mESC_WTvSKO.common_PKc.f2k.TKO_mCpG.stat.runlog 2>&1 &
tail -f mESC_WTvSKO.common_PKc.f2k.WT_mCpG.stat.runlog
tail -f mESC_WTvSKO.common_PKc.f2k.SKO_mCpG.stat.runlog
tail -f mESC_WTvSKO.common_PKc.f2k.TKO_mCpG.stat.runlog
## mESC_WTvSKO.SKO_only_PKc.f2k  mCpG stat
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_WT_SKO_intersect.SKO_only.region.centers.bed -b mESC_WT_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WTvSKO.SKO_only_PKc.f2k.WT_mCpG.stat.txt >mESC_WTvSKO.SKO_only_PKc.f2k.WT_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_WT_SKO_intersect.SKO_only.region.centers.bed -b mESC_SKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WTvSKO.SKO_only_PKc.f2k.SKO_mCpG.stat.txt >mESC_WTvSKO.SKO_only_PKc.f2k.SKO_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_WT_SKO_intersect.SKO_only.region.centers.bed -b mESC_TKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WTvSKO.SKO_only_PKc.f2k.TKO_mCpG.stat.txt >mESC_WTvSKO.SKO_only_PKc.f2k.TKO_mCpG.stat.runlog 2>&1 &
tail -f mESC_WTvSKO.SKO_only_PKc.f2k.WT_mCpG.stat.runlog
tail -f mESC_WTvSKO.SKO_only_PKc.f2k.SKO_mCpG.stat.runlog
tail -f mESC_WTvSKO.SKO_only_PKc.f2k.TKO_mCpG.stat.runlog

# 由于绘图仅需要前几列数据，所以仅输出前几列用于后续作图分析
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/12_mESC_WTvsSKO_3cCTCF_summits_f2k_mCpG_stat_plot
cut -f 1,2,3,4,5,6 mESC_WTvSKO.WT_only_PKc.f2k.SKO_mCpG.stat.txt >mESC_WTvSKO.WT_only_PKc.f2k.SKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WTvSKO.WT_only_PKc.f2k.TKO_mCpG.stat.txt >mESC_WTvSKO.WT_only_PKc.f2k.TKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WTvSKO.WT_only_PKc.f2k.WT_mCpG.stat.txt >mESC_WTvSKO.WT_only_PKc.f2k.WT_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WTvSKO.common_PKc.f2k.SKO_mCpG.stat.txt >mESC_WTvSKO.common_PKc.f2k.SKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WTvSKO.common_PKc.f2k.TKO_mCpG.stat.txt >mESC_WTvSKO.common_PKc.f2k.TKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WTvSKO.common_PKc.f2k.WT_mCpG.stat.txt >mESC_WTvSKO.common_PKc.f2k.WT_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WTvSKO.SKO_only_PKc.f2k.SKO_mCpG.stat.txt >mESC_WTvSKO.SKO_only_PKc.f2k.SKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WTvSKO.SKO_only_PKc.f2k.TKO_mCpG.stat.txt >mESC_WTvSKO.SKO_only_PKc.f2k.TKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WTvSKO.SKO_only_PKc.f2k.WT_mCpG.stat.txt >mESC_WTvSKO.SKO_only_PKc.f2k.WT_mCpG.stat.cut1_6.txt

# 统计结果绘图
# 测试发现 平滑曲线参数不固定，因此需要每个结果单独编写脚本
# 因此可以使用 jupyter notebook，方便调试参数

```



## WT_vs_TKO

### 创建三类 peak center 和甲基化位点文件链接

```bash
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
mkdir 12_mESC_WTvsTKO_3cCTCF_summits_f2k_mCpG_stat_plot
cd 12_mESC_WTvsTKO_3cCTCF_summits_f2k_mCpG_stat_plot
#  创建三类 peak center 文件连接 (bedtools intersect)
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/mESC_WT_TKO_intersect.WT_only.region.centers.bed mESC_WT_TKO_intersect.WT_only.region.centers.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/mESC_WT_TKO_intersect.common.region.centers.bed mESC_WT_TKO_intersect.common.region.centers.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/mESC_WT_TKO_intersect.TKO_only.region.centers.bed mESC_WT_TKO_intersect.TKO_only.region.centers.bed
wc -l *.region.centers.bed
# 13492 mESC_WT_TKO_intersect.WT_only.region.centers.bed
# 57713 mESC_WT_TKO_intersect.common.region.centers.bed
#  4663 mESC_WT_TKO_intersect.TKO_only.region.centers.bed

# 创建CpG甲基化位点结果文件连接
ln -s /data1/sunwenju/projects_on_940xa/2024031102_WN_mESC_BS_Seq_add/2_all_samples_genome_wide_CpG_sites_info/mESC_SKO_BS_Seq_rep123.all_CpG_result.txt mESC_SKO_BS_Seq_rep123.all_CpG_result.txt
ln -s /data1/sunwenju/projects_on_940xa/2024031102_WN_mESC_BS_Seq_add/2_all_samples_genome_wide_CpG_sites_info/mESC_TKO_BS_Seq_rep123.all_CpG_result.txt mESC_TKO_BS_Seq_rep123.all_CpG_result.txt
ln -s /data1/sunwenju/projects_on_940xa/2024031102_WN_mESC_BS_Seq_add/2_all_samples_genome_wide_CpG_sites_info/mESC_WT_BS_Seq_rep123.all_CpG_result.txt mESC_WT_BS_Seq_rep123.all_CpG_result.txt
ls -1 mESC_*all_CpG_result.txt

```

### 统计三类 CTCF peak center f2k mCpG 分布并绘图

```bash
# 统计三类 CTCF peak centers 上下游甲基化位点的分布
# 总结之前的经验，可以全部单独运行
## mESC_WTvTKO.WT_only_PKc.f2k  mCpG stat
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/12_mESC_WTvsTKO_3cCTCF_summits_f2k_mCpG_stat_plot
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_WT_TKO_intersect.WT_only.region.centers.bed -b mESC_WT_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WTvTKO.WT_only_PKc.f2k.WT_mCpG.stat.txt >mESC_WTvTKO.WT_only_PKc.f2k.WT_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_WT_TKO_intersect.WT_only.region.centers.bed -b mESC_SKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WTvTKO.WT_only_PKc.f2k.SKO_mCpG.stat.txt >mESC_WTvTKO.WT_only_PKc.f2k.SKO_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_WT_TKO_intersect.WT_only.region.centers.bed -b mESC_TKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WTvTKO.WT_only_PKc.f2k.TKO_mCpG.stat.txt >mESC_WTvTKO.WT_only_PKc.f2k.TKO_mCpG.stat.runlog 2>&1 &
tail -f mESC_WTvTKO.WT_only_PKc.f2k.WT_mCpG.stat.runlog
tail -f mESC_WTvTKO.WT_only_PKc.f2k.SKO_mCpG.stat.runlog
tail -f mESC_WTvTKO.WT_only_PKc.f2k.TKO_mCpG.stat.runlog
## mESC_WTvTKO.common_PKc.f2k  mCpG stat
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/12_mESC_WTvsTKO_3cCTCF_summits_f2k_mCpG_stat_plot
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_WT_TKO_intersect.common.region.centers.bed -b mESC_WT_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WTvTKO.common_PKc.f2k.WT_mCpG.stat.txt >mESC_WTvTKO.common_PKc.f2k.WT_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_WT_TKO_intersect.common.region.centers.bed -b mESC_SKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WTvTKO.common_PKc.f2k.SKO_mCpG.stat.txt >mESC_WTvTKO.common_PKc.f2k.SKO_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_WT_TKO_intersect.common.region.centers.bed -b mESC_TKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WTvTKO.common_PKc.f2k.TKO_mCpG.stat.txt >mESC_WTvTKO.common_PKc.f2k.TKO_mCpG.stat.runlog 2>&1 &
tail -f mESC_WTvTKO.common_PKc.f2k.WT_mCpG.stat.runlog
tail -f mESC_WTvTKO.common_PKc.f2k.SKO_mCpG.stat.runlog
tail -f mESC_WTvTKO.common_PKc.f2k.TKO_mCpG.stat.runlog
## mESC_WTvTKO.TKO_only_PKc.f2k  mCpG stat
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/12_mESC_WTvsTKO_3cCTCF_summits_f2k_mCpG_stat_plot
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_WT_TKO_intersect.TKO_only.region.centers.bed -b mESC_WT_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WTvTKO.TKO_only_PKc.f2k.WT_mCpG.stat.txt >mESC_WTvTKO.TKO_only_PKc.f2k.WT_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_WT_TKO_intersect.TKO_only.region.centers.bed -b mESC_SKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WTvTKO.TKO_only_PKc.f2k.SKO_mCpG.stat.txt >mESC_WTvTKO.TKO_only_PKc.f2k.SKO_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_WT_TKO_intersect.TKO_only.region.centers.bed -b mESC_TKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WTvTKO.TKO_only_PKc.f2k.TKO_mCpG.stat.txt >mESC_WTvTKO.TKO_only_PKc.f2k.TKO_mCpG.stat.runlog 2>&1 &
tail -f mESC_WTvTKO.TKO_only_PKc.f2k.WT_mCpG.stat.runlog
tail -f mESC_WTvTKO.TKO_only_PKc.f2k.SKO_mCpG.stat.runlog
tail -f mESC_WTvTKO.TKO_only_PKc.f2k.TKO_mCpG.stat.runlog

# 由于绘图仅需要前几列数据，所以仅输出前几列用于后续作图分析
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/12_mESC_WTvsTKO_3cCTCF_summits_f2k_mCpG_stat_plot
cut -f 1,2,3,4,5,6 mESC_WTvTKO.WT_only_PKc.f2k.SKO_mCpG.stat.txt >mESC_WTvTKO.WT_only_PKc.f2k.SKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WTvTKO.WT_only_PKc.f2k.TKO_mCpG.stat.txt >mESC_WTvTKO.WT_only_PKc.f2k.TKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WTvTKO.WT_only_PKc.f2k.WT_mCpG.stat.txt >mESC_WTvTKO.WT_only_PKc.f2k.WT_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WTvTKO.common_PKc.f2k.SKO_mCpG.stat.txt >mESC_WTvTKO.common_PKc.f2k.SKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WTvTKO.common_PKc.f2k.TKO_mCpG.stat.txt >mESC_WTvTKO.common_PKc.f2k.TKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WTvTKO.common_PKc.f2k.WT_mCpG.stat.txt >mESC_WTvTKO.common_PKc.f2k.WT_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WTvTKO.TKO_only_PKc.f2k.SKO_mCpG.stat.txt >mESC_WTvTKO.TKO_only_PKc.f2k.SKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WTvTKO.TKO_only_PKc.f2k.TKO_mCpG.stat.txt >mESC_WTvTKO.TKO_only_PKc.f2k.TKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WTvTKO.TKO_only_PKc.f2k.WT_mCpG.stat.txt >mESC_WTvTKO.TKO_only_PKc.f2k.WT_mCpG.stat.cut1_6.txt

# 统计结果绘图
# 测试发现 平滑曲线参数不固定，因此需要每个结果单独编写脚本
# 因此可以使用 jupyter notebook，方便调试参数

```



## SKO_vs_TKO

### 创建三类 peak center 和甲基化位点文件链接

```bash
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
mkdir 12_mESC_SKOvsTKO_3cCTCF_summits_f2k_mCpG_stat_plot
cd 12_mESC_SKOvsTKO_3cCTCF_summits_f2k_mCpG_stat_plot
#  创建三类 peak center 文件连接 (bedtools intersect)
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/mESC_SKO_TKO_intersect.SKO_only.region.centers.bed mESC_SKO_TKO_intersect.SKO_only.region.centers.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/mESC_SKO_TKO_intersect.common.region.centers.bed mESC_SKO_TKO_intersect.common.region.centers.bed
ln -s ../6_macs2_peak_calling_2_for_mESC_WT_CTCF_rep123_together/mESC_SKO_TKO_intersect.TKO_only.region.centers.bed mESC_SKO_TKO_intersect.TKO_only.region.centers.bed
wc -l *.region.centers.bed
#  6165 mESC_SKO_TKO_intersect.SKO_only.region.centers.bed
# 53357 mESC_SKO_TKO_intersect.common.region.centers.bed
#  8996 mESC_SKO_TKO_intersect.TKO_only.region.centers.bed

# 创建CpG甲基化位点结果文件连接
ln -s /data1/sunwenju/projects_on_940xa/2024031102_WN_mESC_BS_Seq_add/2_all_samples_genome_wide_CpG_sites_info/mESC_SKO_BS_Seq_rep123.all_CpG_result.txt mESC_SKO_BS_Seq_rep123.all_CpG_result.txt
ln -s /data1/sunwenju/projects_on_940xa/2024031102_WN_mESC_BS_Seq_add/2_all_samples_genome_wide_CpG_sites_info/mESC_TKO_BS_Seq_rep123.all_CpG_result.txt mESC_TKO_BS_Seq_rep123.all_CpG_result.txt
ln -s /data1/sunwenju/projects_on_940xa/2024031102_WN_mESC_BS_Seq_add/2_all_samples_genome_wide_CpG_sites_info/mESC_WT_BS_Seq_rep123.all_CpG_result.txt mESC_WT_BS_Seq_rep123.all_CpG_result.txt
ls -1 mESC_*all_CpG_result.txt

```

### 统计三类 CTCF peak center f2k mCpG 分布并绘图

```bash
# 统计三类 CTCF peak centers 上下游甲基化位点的分布
# 总结之前的经验，可以全部单独运行
## mESC_SKOvTKO.SKO_only_PKc.f2k  mCpG stat
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/12_mESC_SKOvsTKO_3cCTCF_summits_f2k_mCpG_stat_plot
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_SKO_TKO_intersect.SKO_only.region.centers.bed -b mESC_WT_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_SKOvTKO.SKO_only_PKc.f2k.WT_mCpG.stat.txt >mESC_SKOvTKO.SKO_only_PKc.f2k.WT_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_SKO_TKO_intersect.SKO_only.region.centers.bed -b mESC_SKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_SKOvTKO.SKO_only_PKc.f2k.SKO_mCpG.stat.txt >mESC_SKOvTKO.SKO_only_PKc.f2k.SKO_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_SKO_TKO_intersect.SKO_only.region.centers.bed -b mESC_TKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_SKOvTKO.SKO_only_PKc.f2k.TKO_mCpG.stat.txt >mESC_SKOvTKO.SKO_only_PKc.f2k.TKO_mCpG.stat.runlog 2>&1 &
tail -f mESC_SKOvTKO.SKO_only_PKc.f2k.WT_mCpG.stat.runlog
tail -f mESC_SKOvTKO.SKO_only_PKc.f2k.SKO_mCpG.stat.runlog
tail -f mESC_SKOvTKO.SKO_only_PKc.f2k.TKO_mCpG.stat.runlog
## mESC_SKOvTKO.common_PKc.f2k  mCpG stat
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/12_mESC_SKOvsTKO_3cCTCF_summits_f2k_mCpG_stat_plot
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_SKO_TKO_intersect.common.region.centers.bed -b mESC_WT_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_SKOvTKO.common_PKc.f2k.WT_mCpG.stat.txt >mESC_SKOvTKO.common_PKc.f2k.WT_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_SKO_TKO_intersect.common.region.centers.bed -b mESC_SKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_SKOvTKO.common_PKc.f2k.SKO_mCpG.stat.txt >mESC_SKOvTKO.common_PKc.f2k.SKO_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_SKO_TKO_intersect.common.region.centers.bed -b mESC_TKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_SKOvTKO.common_PKc.f2k.TKO_mCpG.stat.txt >mESC_SKOvTKO.common_PKc.f2k.TKO_mCpG.stat.runlog 2>&1 &
tail -f mESC_SKOvTKO.common_PKc.f2k.WT_mCpG.stat.runlog
tail -f mESC_SKOvTKO.common_PKc.f2k.SKO_mCpG.stat.runlog
tail -f mESC_SKOvTKO.common_PKc.f2k.TKO_mCpG.stat.runlog
## mESC_SKOvTKO.TKO_only_PKc.f2k  mCpG stat
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/12_mESC_SKOvsTKO_3cCTCF_summits_f2k_mCpG_stat_plot
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_SKO_TKO_intersect.TKO_only.region.centers.bed -b mESC_WT_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_SKOvTKO.TKO_only_PKc.f2k.WT_mCpG.stat.txt >mESC_SKOvTKO.TKO_only_PKc.f2k.WT_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_SKO_TKO_intersect.TKO_only.region.centers.bed -b mESC_SKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_SKOvTKO.TKO_only_PKc.f2k.SKO_mCpG.stat.txt >mESC_SKOvTKO.TKO_only_PKc.f2k.SKO_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a mESC_SKO_TKO_intersect.TKO_only.region.centers.bed -b mESC_TKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_SKOvTKO.TKO_only_PKc.f2k.TKO_mCpG.stat.txt >mESC_SKOvTKO.TKO_only_PKc.f2k.TKO_mCpG.stat.runlog 2>&1 &
tail -f mESC_SKOvTKO.TKO_only_PKc.f2k.WT_mCpG.stat.runlog
tail -f mESC_SKOvTKO.TKO_only_PKc.f2k.SKO_mCpG.stat.runlog
tail -f mESC_SKOvTKO.TKO_only_PKc.f2k.TKO_mCpG.stat.runlog

# 由于绘图仅需要前几列数据，所以仅输出前几列用于后续作图分析
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/12_mESC_SKOvsTKO_3cCTCF_summits_f2k_mCpG_stat_plot
cut -f 1,2,3,4,5,6 mESC_SKOvTKO.SKO_only_PKc.f2k.SKO_mCpG.stat.txt >mESC_SKOvTKO.SKO_only_PKc.f2k.SKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_SKOvTKO.SKO_only_PKc.f2k.TKO_mCpG.stat.txt >mESC_SKOvTKO.SKO_only_PKc.f2k.TKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_SKOvTKO.SKO_only_PKc.f2k.WT_mCpG.stat.txt >mESC_SKOvTKO.SKO_only_PKc.f2k.WT_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_SKOvTKO.common_PKc.f2k.SKO_mCpG.stat.txt >mESC_SKOvTKO.common_PKc.f2k.SKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_SKOvTKO.common_PKc.f2k.TKO_mCpG.stat.txt >mESC_SKOvTKO.common_PKc.f2k.TKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_SKOvTKO.common_PKc.f2k.WT_mCpG.stat.txt >mESC_SKOvTKO.common_PKc.f2k.WT_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_SKOvTKO.TKO_only_PKc.f2k.SKO_mCpG.stat.txt >mESC_SKOvTKO.TKO_only_PKc.f2k.SKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_SKOvTKO.TKO_only_PKc.f2k.TKO_mCpG.stat.txt >mESC_SKOvTKO.TKO_only_PKc.f2k.TKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_SKOvTKO.TKO_only_PKc.f2k.WT_mCpG.stat.txt >mESC_SKOvTKO.TKO_only_PKc.f2k.WT_mCpG.stat.cut1_6.txt

# 统计结果绘图
# 测试发现 平滑曲线参数不固定，因此需要每个结果单独编写脚本
# 因此可以使用 jupyter notebook，方便调试参数

```

# 三组实验 intervene 结果7类 peak centers  f2k 甲基化分布统计绘图

## 创建7类 peak center 和 甲基化位点文件链接

> 三组实验 intervene 结果7类 peak centers 在 11_plotHeatMap_for_all_use_mESC_WTvsSKOvsTKO_CTCF_summits 中已经构造过

```bash
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq
mkdir 13_mESC_WST_intervene_7cCTCF_summits_f2k_mCpG_stat_plot
cd 13_mESC_WST_intervene_7cCTCF_summits_f2k_mCpG_stat_plot
#  创建7类 peak center 文件连接 (WST peaks intervene)
ln -s ../11_plotHeatMap_for_all_use_mESC_WTvsSKOvsTKO_CTCF_summits/PKIc.100.sorted.bed WST.PKIc.100.sorted.bed
ln -s ../11_plotHeatMap_for_all_use_mESC_WTvsSKOvsTKO_CTCF_summits/PKIc.010.sorted.bed WST.PKIc.010.sorted.bed
ln -s ../11_plotHeatMap_for_all_use_mESC_WTvsSKOvsTKO_CTCF_summits/PKIc.001.sorted.bed WST.PKIc.001.sorted.bed
ln -s ../11_plotHeatMap_for_all_use_mESC_WTvsSKOvsTKO_CTCF_summits/PKIc.110.sorted.bed WST.PKIc.110.sorted.bed
ln -s ../11_plotHeatMap_for_all_use_mESC_WTvsSKOvsTKO_CTCF_summits/PKIc.011.sorted.bed WST.PKIc.011.sorted.bed
ln -s ../11_plotHeatMap_for_all_use_mESC_WTvsSKOvsTKO_CTCF_summits/PKIc.101.sorted.bed WST.PKIc.101.sorted.bed
ln -s ../11_plotHeatMap_for_all_use_mESC_WTvsSKOvsTKO_CTCF_summits/PKIc.111.sorted.bed WST.PKIc.111.sorted.bed
ls -1
wc -l WST.PKIc.*.bed
#  9934 WST.PKIc.100.sorted.bed
#  2610 WST.PKIc.010.sorted.bed
#  3359 WST.PKIc.001.sorted.bed
#  3558 WST.PKIc.110.sorted.bed
#  1301 WST.PKIc.011.sorted.bed
#  5637 WST.PKIc.101.sorted.bed
# 51946 WST.PKIc.111.sorted.bed

# 创建CpG甲基化位点结果文件连接
ln -s /data1/sunwenju/projects_on_940xa/2024031102_WN_mESC_BS_Seq_add/2_all_samples_genome_wide_CpG_sites_info/mESC_SKO_BS_Seq_rep123.all_CpG_result.txt mESC_SKO_BS_Seq_rep123.all_CpG_result.txt
ln -s /data1/sunwenju/projects_on_940xa/2024031102_WN_mESC_BS_Seq_add/2_all_samples_genome_wide_CpG_sites_info/mESC_TKO_BS_Seq_rep123.all_CpG_result.txt mESC_TKO_BS_Seq_rep123.all_CpG_result.txt
ln -s /data1/sunwenju/projects_on_940xa/2024031102_WN_mESC_BS_Seq_add/2_all_samples_genome_wide_CpG_sites_info/mESC_WT_BS_Seq_rep123.all_CpG_result.txt mESC_WT_BS_Seq_rep123.all_CpG_result.txt
ls -1 mESC_*all_CpG_result.txt

```

## 统计7类 peak center f2k 甲基化分布

```bash
# 统计三类 CTCF peak centers 上下游甲基化位点的分布
# 总结之前的经验，可以全部单独运行
## mESC_WST_PKIc.100.PKc.f2k  mCpG stat
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a WST.PKIc.100.sorted.bed -b mESC_WT_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WST_PKIc.100.PKc.f2k.WT_mCpG.stat.txt >mESC_WST_PKIc.100.PKc.f2k.WT_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a WST.PKIc.100.sorted.bed -b mESC_SKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WST_PKIc.100.PKc.f2k.SKO_mCpG.stat.txt >mESC_WST_PKIc.100.PKc.f2k.SKO_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a WST.PKIc.100.sorted.bed -b mESC_TKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WST_PKIc.100.PKc.f2k.TKO_mCpG.stat.txt >mESC_WST_PKIc.100.PKc.f2k.TKO_mCpG.stat.runlog 2>&1 &
tail -f mESC_WST_PKIc.100.PKc.f2k.WT_mCpG.stat.runlog
tail -f mESC_WST_PKIc.100.PKc.f2k.SKO_mCpG.stat.runlog
tail -f mESC_WST_PKIc.100.PKc.f2k.TKO_mCpG.stat.runlog
## mESC_WST_PKIc.010.PKc.f2k  mCpG stat
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a WST.PKIc.010.sorted.bed -b mESC_WT_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WST_PKIc.010.PKc.f2k.WT_mCpG.stat.txt >mESC_WST_PKIc.010.PKc.f2k.WT_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a WST.PKIc.010.sorted.bed -b mESC_SKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WST_PKIc.010.PKc.f2k.SKO_mCpG.stat.txt >mESC_WST_PKIc.010.PKc.f2k.SKO_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a WST.PKIc.010.sorted.bed -b mESC_TKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WST_PKIc.010.PKc.f2k.TKO_mCpG.stat.txt >mESC_WST_PKIc.010.PKc.f2k.TKO_mCpG.stat.runlog 2>&1 &
tail -f mESC_WST_PKIc.010.PKc.f2k.WT_mCpG.stat.runlog
tail -f mESC_WST_PKIc.010.PKc.f2k.SKO_mCpG.stat.runlog
tail -f mESC_WST_PKIc.010.PKc.f2k.TKO_mCpG.stat.runlog
## mESC_WST_PKIc.001.PKc.f2k  mCpG stat
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a WST.PKIc.001.sorted.bed -b mESC_WT_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WST_PKIc.001.PKc.f2k.WT_mCpG.stat.txt >mESC_WST_PKIc.001.PKc.f2k.WT_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a WST.PKIc.001.sorted.bed -b mESC_SKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WST_PKIc.001.PKc.f2k.SKO_mCpG.stat.txt >mESC_WST_PKIc.001.PKc.f2k.SKO_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a WST.PKIc.001.sorted.bed -b mESC_TKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WST_PKIc.001.PKc.f2k.TKO_mCpG.stat.txt >mESC_WST_PKIc.001.PKc.f2k.TKO_mCpG.stat.runlog 2>&1 &
tail -f mESC_WST_PKIc.001.PKc.f2k.WT_mCpG.stat.runlog
tail -f mESC_WST_PKIc.001.PKc.f2k.SKO_mCpG.stat.runlog
tail -f mESC_WST_PKIc.001.PKc.f2k.TKO_mCpG.stat.runlog
## mESC_WST_PKIc.110.PKc.f2k  mCpG stat
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a WST.PKIc.110.sorted.bed -b mESC_WT_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WST_PKIc.110.PKc.f2k.WT_mCpG.stat.txt >mESC_WST_PKIc.110.PKc.f2k.WT_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a WST.PKIc.110.sorted.bed -b mESC_SKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WST_PKIc.110.PKc.f2k.SKO_mCpG.stat.txt >mESC_WST_PKIc.110.PKc.f2k.SKO_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a WST.PKIc.110.sorted.bed -b mESC_TKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WST_PKIc.110.PKc.f2k.TKO_mCpG.stat.txt >mESC_WST_PKIc.110.PKc.f2k.TKO_mCpG.stat.runlog 2>&1 &
tail -f mESC_WST_PKIc.110.PKc.f2k.WT_mCpG.stat.runlog
tail -f mESC_WST_PKIc.110.PKc.f2k.SKO_mCpG.stat.runlog
tail -f mESC_WST_PKIc.110.PKc.f2k.TKO_mCpG.stat.runlog
## mESC_WST_PKIc.011.PKc.f2k  mCpG stat
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a WST.PKIc.011.sorted.bed -b mESC_WT_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WST_PKIc.011.PKc.f2k.WT_mCpG.stat.txt >mESC_WST_PKIc.011.PKc.f2k.WT_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a WST.PKIc.011.sorted.bed -b mESC_SKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WST_PKIc.011.PKc.f2k.SKO_mCpG.stat.txt >mESC_WST_PKIc.011.PKc.f2k.SKO_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a WST.PKIc.011.sorted.bed -b mESC_TKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WST_PKIc.011.PKc.f2k.TKO_mCpG.stat.txt >mESC_WST_PKIc.011.PKc.f2k.TKO_mCpG.stat.runlog 2>&1 &
tail -f mESC_WST_PKIc.011.PKc.f2k.WT_mCpG.stat.runlog
tail -f mESC_WST_PKIc.011.PKc.f2k.SKO_mCpG.stat.runlog
tail -f mESC_WST_PKIc.011.PKc.f2k.TKO_mCpG.stat.runlog
## mESC_WST_PKIc.101.PKc.f2k  mCpG stat
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a WST.PKIc.101.sorted.bed -b mESC_WT_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WST_PKIc.101.PKc.f2k.WT_mCpG.stat.txt >mESC_WST_PKIc.101.PKc.f2k.WT_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a WST.PKIc.101.sorted.bed -b mESC_SKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WST_PKIc.101.PKc.f2k.SKO_mCpG.stat.txt >mESC_WST_PKIc.101.PKc.f2k.SKO_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a WST.PKIc.101.sorted.bed -b mESC_TKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WST_PKIc.101.PKc.f2k.TKO_mCpG.stat.txt >mESC_WST_PKIc.101.PKc.f2k.TKO_mCpG.stat.runlog 2>&1 &
tail -f mESC_WST_PKIc.101.PKc.f2k.WT_mCpG.stat.runlog
tail -f mESC_WST_PKIc.101.PKc.f2k.SKO_mCpG.stat.runlog
tail -f mESC_WST_PKIc.101.PKc.f2k.TKO_mCpG.stat.runlog
## mESC_WST_PKIc.111.PKc.f2k  mCpG stat
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a WST.PKIc.111.sorted.bed -b mESC_WT_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WST_PKIc.111.PKc.f2k.WT_mCpG.stat.txt >mESC_WST_PKIc.111.PKc.f2k.WT_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a WST.PKIc.111.sorted.bed -b mESC_SKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WST_PKIc.111.PKc.f2k.SKO_mCpG.stat.txt >mESC_WST_PKIc.111.PKc.f2k.SKO_mCpG.stat.runlog 2>&1 &
nohup time python3 /home/sunwenju/1_softwares/0_my_useful_scripts/make_meta_profile_for_point_frank.py -a WST.PKIc.111.sorted.bed -b mESC_TKO_BS_Seq_rep123.all_CpG_result.txt -c 1 -d 1000 -o mESC_WST_PKIc.111.PKc.f2k.TKO_mCpG.stat.txt >mESC_WST_PKIc.111.PKc.f2k.TKO_mCpG.stat.runlog 2>&1 &
tail -f mESC_WST_PKIc.111.PKc.f2k.WT_mCpG.stat.runlog
tail -f mESC_WST_PKIc.111.PKc.f2k.SKO_mCpG.stat.runlog
tail -f mESC_WST_PKIc.111.PKc.f2k.TKO_mCpG.stat.runlog

```

## 提取绘图数据并在jupyter里绘图

```bash
# 由于绘图仅需要前几列数据，所以仅输出前几列用于后续作图分析
cd /data1/sunwenju/projects_on_940xa/2024082801_WN_mESC_TKO_SKO_CTCF_ChIP_Seq/13_mESC_WST_intervene_7cCTCF_summits_f2k_mCpG_stat_plot
#
cut -f 1,2,3,4,5,6 mESC_WST_PKIc.100.PKc.f2k.SKO_mCpG.stat.txt >mESC_WST_PKIc.100.PKc.f2k.SKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WST_PKIc.100.PKc.f2k.TKO_mCpG.stat.txt >mESC_WST_PKIc.100.PKc.f2k.TKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WST_PKIc.100.PKc.f2k.WT_mCpG.stat.txt >mESC_WST_PKIc.100.PKc.f2k.WT_mCpG.stat.cut1_6.txt
#
cut -f 1,2,3,4,5,6 mESC_WST_PKIc.010.PKc.f2k.SKO_mCpG.stat.txt >mESC_WST_PKIc.010.PKc.f2k.SKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WST_PKIc.010.PKc.f2k.TKO_mCpG.stat.txt >mESC_WST_PKIc.010.PKc.f2k.TKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WST_PKIc.010.PKc.f2k.WT_mCpG.stat.txt >mESC_WST_PKIc.010.PKc.f2k.WT_mCpG.stat.cut1_6.txt
#
cut -f 1,2,3,4,5,6 mESC_WST_PKIc.001.PKc.f2k.SKO_mCpG.stat.txt >mESC_WST_PKIc.001.PKc.f2k.SKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WST_PKIc.001.PKc.f2k.TKO_mCpG.stat.txt >mESC_WST_PKIc.001.PKc.f2k.TKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WST_PKIc.001.PKc.f2k.WT_mCpG.stat.txt >mESC_WST_PKIc.001.PKc.f2k.WT_mCpG.stat.cut1_6.txt
#
cut -f 1,2,3,4,5,6 mESC_WST_PKIc.110.PKc.f2k.SKO_mCpG.stat.txt >mESC_WST_PKIc.110.PKc.f2k.SKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WST_PKIc.110.PKc.f2k.TKO_mCpG.stat.txt >mESC_WST_PKIc.110.PKc.f2k.TKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WST_PKIc.110.PKc.f2k.WT_mCpG.stat.txt >mESC_WST_PKIc.110.PKc.f2k.WT_mCpG.stat.cut1_6.txt
#
cut -f 1,2,3,4,5,6 mESC_WST_PKIc.011.PKc.f2k.SKO_mCpG.stat.txt >mESC_WST_PKIc.011.PKc.f2k.SKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WST_PKIc.011.PKc.f2k.TKO_mCpG.stat.txt >mESC_WST_PKIc.011.PKc.f2k.TKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WST_PKIc.011.PKc.f2k.WT_mCpG.stat.txt >mESC_WST_PKIc.011.PKc.f2k.WT_mCpG.stat.cut1_6.txt
#
cut -f 1,2,3,4,5,6 mESC_WST_PKIc.101.PKc.f2k.SKO_mCpG.stat.txt >mESC_WST_PKIc.101.PKc.f2k.SKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WST_PKIc.101.PKc.f2k.TKO_mCpG.stat.txt >mESC_WST_PKIc.101.PKc.f2k.TKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WST_PKIc.101.PKc.f2k.WT_mCpG.stat.txt >mESC_WST_PKIc.101.PKc.f2k.WT_mCpG.stat.cut1_6.txt
#
cut -f 1,2,3,4,5,6 mESC_WST_PKIc.111.PKc.f2k.SKO_mCpG.stat.txt >mESC_WST_PKIc.111.PKc.f2k.SKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WST_PKIc.111.PKc.f2k.TKO_mCpG.stat.txt >mESC_WST_PKIc.111.PKc.f2k.TKO_mCpG.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 mESC_WST_PKIc.111.PKc.f2k.WT_mCpG.stat.txt >mESC_WST_PKIc.111.PKc.f2k.WT_mCpG.stat.cut1_6.txt

# 统计结果绘图
# 测试发现 平滑曲线参数不固定，因此需要每个结果单独编写脚本
# 因此可以使用 jupyter notebook，方便调试参数

```



