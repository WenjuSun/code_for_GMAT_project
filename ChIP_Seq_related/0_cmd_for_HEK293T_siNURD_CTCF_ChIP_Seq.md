# QC

```bash
cd /home/sunwenju/projects_on_940xa/
mkdir 2023081601_WN_siNURD_CTCF_ChIP_Seq
# 根据吴楠提供的样本信息表 
# BGI-20230716-IP-CTCF-1-96-C2样本信息表.xlsx
# BGI-20230716-IP-CTCF-2-96-C3样本信息表 .xlsx
# 构造数据连接创建脚本
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq
bash 0_make_and_rename_soft_links_for_WN_siNURD_CTCF_ChIP_Seq_raw_data.sh
ls -1

# 原始数据进行质量评估
mkdir 1_all_293T_siNURD_CTCF_ChIP_Seq_data_fastqc_result
nohup fastqc -outdir 1_all_293T_siNURD_CTCF_ChIP_Seq_data_fastqc_result --threads 72 *.fq.gz >1_all_293T_siNURD_CTCF_ChIP_Seq_data_fastqc.runlog 2>&1 &
tail -f 1_all_293T_siNURD_CTCF_ChIP_Seq_data_fastqc.runlog

# 并用 multiqc 进行结果汇总
# 提取原始数据reads数量
cd 1_all_293T_siNURD_CTCF_ChIP_Seq_data_fastqc_result
multiqc ./*fastqc.zip --ignore *.html --no-megaqc-upload
cd ../
# PE100，无接头
```

# Genome mapping

```bash
# 可使用如下命令构造比对脚本， 然后运行
for i in `ls -1 *rep[123].R1.fq.gz`;do echo 'bowtie2 -t -p 72 --no-unal -x /home/sunwenju/3_genomeData/human/hg19/bowtie2Index/hg19bt2idx -1 '$i' -2 '${i%.R1.fq.gz}.R2.fq.gz' -S '${i%.R1.fq.gz}.bowtie2.sam' --un-gz '${i%.R1.fq.gz}.notAlign.unpairedReads.fq.gz' --un-conc-gz '${i%.R1.fq.gz}.notAlign.pairedReads.fq.gz' >'${i%.R1.fq.gz}.bowtie2.runlog' 2>&1';done >3_all_HEK293T_siNURD_CTCF_ChIP_Seq_bowtie2_map.sh
nohup bash 3_all_HEK293T_siNURD_CTCF_ChIP_Seq_bowtie2_map.sh >3_all_HEK293T_siNURD_CTCF_ChIP_Seq_bowtie2_map.sh.runlog 2>&1 &
# tail -f 3_all_HEK293T_siNURD_CTCF_ChIP_Seq_bowtie2_map.sh.runlog
tail -f `ls -1 *bowtie2.runlog|tail -n 1`

# 使用 multiqc 对 bowtie2 比对情况进行汇总
multiqc --no-megaqc-upload -m bowtie2 -o 3_all_HEK293T_siNURD_CTCF_ChIP_Seq_bowtie2_map_results_summary ./*.bowtie2.runlog

```



# Convert sam to sorted bam

```bash
# 将比对sam结果转换为 按座位排序的 bam 并建立索引
# for i in `ls -1 *.bowtie2.sam|sort`; do echo $i; samtools view --threads 36 -bS $i | samtools sort --threads 36 -o ${i%.sam}.sorted.bam; samtools index ${i%.sam}.sorted.bam; done
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq
nohup bash 3_WN_siNURD_CTCF_ChIP_Seq_data_bowtie2_sam_to_sorted_bam.sh >3_WN_siNURD_CTCF_ChIP_Seq_data_bowtie2_sam_to_sorted_bam.sh.runlog 2>&1 &
tail -f 3_WN_siNURD_CTCF_ChIP_Seq_data_bowtie2_sam_to_sorted_bam.sh.runlog
# rm *.bowtie2.sam

# 使用samtools flagstat 统计比对情况
# 每个bam文件不到一分钟
for i in `ls -1 *bowtie2.sorted.bam`;do echo $i; samtools flagstat --threads 16 $i;done >3_293T_IAA72h_ChIP_Seq_bowtie2_map_samtools_flagstat.txt

# 统计比对到各个染色体上的 总read数量 
for i in `ls -1 *bowtie2.sorted.bam|sort`;do echo $i; samtools view $i |grep -v -e '^@' |cut -f 3 |sort |uniq -c |awk 'BEGIN{OFS="\t"}{if($0 !~ /_/){print $2,$1}}' |sort -k1,1 >${i%.sam}.chr_mapped_reads_stat.txt;done
# 统计比对到各个染色体上的 uniq read数量 
for i in `ls -1 *bowtie2.sorted.bam|sort`;do echo $i; samtools view $i |grep -v -e '^@' |cut -f 3,4 |sort |uniq |cut -f 1 |sort |uniq -c |awk 'BEGIN{OFS="\t"}{if($0 !~ /_/){print $2,$1}}' |sort -k1,1 >${i%.sam}.chr_mapped_uniq_reads_stat.txt;done

# 比对到各个染色体上的 read数量 统计结果合并
paste `ls -1 *chr_mapped_reads_stat.txt|sort` |head
# 检测各个统计结果的染色体数量和顺序是否一致
paste `ls -1 *chr_mapped_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"}{if($1!=$3 || $1!=$5 || $1!=$7 || $1!=$9 || $1!=$11 || $1!=$13 || $1!=$15 || $1!=$17 || $1!=$19 || $1!=$21 || $1!=$23 || $1!=$25 || $1!=$27 || $1!=$29 || $1!=$31 || $1!=$33 || $1!=$35 || $1!=$37 || $1!=$39 || $1!=$41 || $1!=$43 || $1!=$45 || $1!=$47 || $1!=$49 || $1!=$51 || $1!=$53 || $1!=$55 || $1!=$57 || $1!=$59 || $1!=$61 || $1!=$63 || $1!=$65 || $1!=$67 || $1!=$69 || $1!=$71){print $1,$2}}' |head
# 染色体数量和顺序一致，则合并数据列
paste `ls -1 *chr_mapped_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"; print "chrom","siCHD3_input_rep1","siCHD3_input_rep2","siCHD3_input_rep3","siCHD3_IP_CTCF_rep1","siCHD3_IP_CTCF_rep2","siCHD3_IP_CTCF_rep3","siCHD4_input_rep1","siCHD4_input_rep2","siCHD4_input_rep3","siCHD4_IP_CTCF_rep1","siCHD4_IP_CTCF_rep2","siCHD4_IP_CTCF_rep3","siCHD5_input_rep1","siCHD5_input_rep2","siCHD5_input_rep3","siCHD5_IP_CTCF_rep1","siCHD5_IP_CTCF_rep2","siCHD5_IP_CTCF_rep3","siMBD2_input_rep1","siMBD2_input_rep2","siMBD2_input_rep3","siMBD2_IP_CTCF_rep1","siMBD2_IP_CTCF_rep2","siMBD2_IP_CTCF_rep3","siMBD3_input_rep1","siMBD3_input_rep2","siMBD3_input_rep3","siMBD3_IP_CTCF_rep1","siMBD3_IP_CTCF_rep2","siMBD3_IP_CTCF_rep3","siNC_input_rep1","siNC_input_rep2","siNC_input_rep3","siNC_IP_CTCF_rep1","siNC_IP_CTCF_rep2","siNC_IP_CTCF_rep3"}{if($1==$3 && $1==$5 && $1==$7 && $1==$9 && $1==$11 && $1==$13 && $1==$15 && $1==$17 && $1==$19 && $1==$21 && $1==$23 && $1==$25 && $1==$27 && $1==$29 && $1==$31 && $1==$33 && $1==$35 && $1==$37 && $1==$39 && $1==$41 && $1==$43 && $1==$45 && $1==$47 && $1==$49 && $1==$51 && $1==$53 && $1==$55 && $1==$57 && $1==$59 && $1==$61 && $1==$63 && $1==$65 && $1==$67 && $1==$69 && $1==$71){print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36,$38,$40,$42,$44,$46,$48,$50,$52,$54,$56,$58,$60,$62,$64,$66,$68,$70,$72}}' |head
paste `ls -1 *chr_mapped_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"; print "chrom","siCHD3_input_rep1","siCHD3_input_rep2","siCHD3_input_rep3","siCHD3_IP_CTCF_rep1","siCHD3_IP_CTCF_rep2","siCHD3_IP_CTCF_rep3","siCHD4_input_rep1","siCHD4_input_rep2","siCHD4_input_rep3","siCHD4_IP_CTCF_rep1","siCHD4_IP_CTCF_rep2","siCHD4_IP_CTCF_rep3","siCHD5_input_rep1","siCHD5_input_rep2","siCHD5_input_rep3","siCHD5_IP_CTCF_rep1","siCHD5_IP_CTCF_rep2","siCHD5_IP_CTCF_rep3","siMBD2_input_rep1","siMBD2_input_rep2","siMBD2_input_rep3","siMBD2_IP_CTCF_rep1","siMBD2_IP_CTCF_rep2","siMBD2_IP_CTCF_rep3","siMBD3_input_rep1","siMBD3_input_rep2","siMBD3_input_rep3","siMBD3_IP_CTCF_rep1","siMBD3_IP_CTCF_rep2","siMBD3_IP_CTCF_rep3","siNC_input_rep1","siNC_input_rep2","siNC_input_rep3","siNC_IP_CTCF_rep1","siNC_IP_CTCF_rep2","siNC_IP_CTCF_rep3"}{if($1==$3 && $1==$5 && $1==$7 && $1==$9 && $1==$11 && $1==$13 && $1==$15 && $1==$17 && $1==$19 && $1==$21 && $1==$23 && $1==$25 && $1==$27 && $1==$29 && $1==$31 && $1==$33 && $1==$35 && $1==$37 && $1==$39 && $1==$41 && $1==$43 && $1==$45 && $1==$47 && $1==$49 && $1==$51 && $1==$53 && $1==$55 && $1==$57 && $1==$59 && $1==$61 && $1==$63 && $1==$65 && $1==$67 && $1==$69 && $1==$71){print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36,$38,$40,$42,$44,$46,$48,$50,$52,$54,$56,$58,$60,$62,$64,$66,$68,$70,$72}}' >all_samples_chr_mapped_reads_stat.summary.txt

# 比对到各个染色体上的 uniq read数量 统计结果合并
paste `ls -1 *chr_mapped_uniq_reads_stat.txt|sort` |head
# 检测各个统计结果的染色体数量和顺序是否一致
paste `ls -1 *chr_mapped_uniq_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"}{if($1!=$3 || $1!=$5 || $1!=$7 || $1!=$9 || $1!=$11 || $1!=$13 || $1!=$15 || $1!=$17 || $1!=$19 || $1!=$21 || $1!=$23 || $1!=$25 || $1!=$27 || $1!=$29 || $1!=$31 || $1!=$33 || $1!=$35 || $1!=$37 || $1!=$39 || $1!=$41 || $1!=$43 || $1!=$45 || $1!=$47 || $1!=$49 || $1!=$51 || $1!=$53 || $1!=$55 || $1!=$57 || $1!=$59 || $1!=$61 || $1!=$63 || $1!=$65 || $1!=$67 || $1!=$69 || $1!=$71){print $1,$2}}' |head
# 染色体数量和顺序一致，则合并数据列
paste `ls -1 *chr_mapped_uniq_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"; print "chrom","siCHD3_input_rep1","siCHD3_input_rep2","siCHD3_input_rep3","siCHD3_IP_CTCF_rep1","siCHD3_IP_CTCF_rep2","siCHD3_IP_CTCF_rep3","siCHD4_input_rep1","siCHD4_input_rep2","siCHD4_input_rep3","siCHD4_IP_CTCF_rep1","siCHD4_IP_CTCF_rep2","siCHD4_IP_CTCF_rep3","siCHD5_input_rep1","siCHD5_input_rep2","siCHD5_input_rep3","siCHD5_IP_CTCF_rep1","siCHD5_IP_CTCF_rep2","siCHD5_IP_CTCF_rep3","siMBD2_input_rep1","siMBD2_input_rep2","siMBD2_input_rep3","siMBD2_IP_CTCF_rep1","siMBD2_IP_CTCF_rep2","siMBD2_IP_CTCF_rep3","siMBD3_input_rep1","siMBD3_input_rep2","siMBD3_input_rep3","siMBD3_IP_CTCF_rep1","siMBD3_IP_CTCF_rep2","siMBD3_IP_CTCF_rep3","siNC_input_rep1","siNC_input_rep2","siNC_input_rep3","siNC_IP_CTCF_rep1","siNC_IP_CTCF_rep2","siNC_IP_CTCF_rep3"}{if($1==$3 && $1==$5 && $1==$7 && $1==$9 && $1==$11 && $1==$13 && $1==$15 && $1==$17 && $1==$19 && $1==$21 && $1==$23 && $1==$25 && $1==$27 && $1==$29 && $1==$31 && $1==$33 && $1==$35 && $1==$37 && $1==$39 && $1==$41 && $1==$43 && $1==$45 && $1==$47 && $1==$49 && $1==$51 && $1==$53 && $1==$55 && $1==$57 && $1==$59 && $1==$61 && $1==$63 && $1==$65 && $1==$67 && $1==$69 && $1==$71){print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36,$38,$40,$42,$44,$46,$48,$50,$52,$54,$56,$58,$60,$62,$64,$66,$68,$70,$72}}' |head
paste `ls -1 *chr_mapped_uniq_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"; print "chrom","siCHD3_input_rep1","siCHD3_input_rep2","siCHD3_input_rep3","siCHD3_IP_CTCF_rep1","siCHD3_IP_CTCF_rep2","siCHD3_IP_CTCF_rep3","siCHD4_input_rep1","siCHD4_input_rep2","siCHD4_input_rep3","siCHD4_IP_CTCF_rep1","siCHD4_IP_CTCF_rep2","siCHD4_IP_CTCF_rep3","siCHD5_input_rep1","siCHD5_input_rep2","siCHD5_input_rep3","siCHD5_IP_CTCF_rep1","siCHD5_IP_CTCF_rep2","siCHD5_IP_CTCF_rep3","siMBD2_input_rep1","siMBD2_input_rep2","siMBD2_input_rep3","siMBD2_IP_CTCF_rep1","siMBD2_IP_CTCF_rep2","siMBD2_IP_CTCF_rep3","siMBD3_input_rep1","siMBD3_input_rep2","siMBD3_input_rep3","siMBD3_IP_CTCF_rep1","siMBD3_IP_CTCF_rep2","siMBD3_IP_CTCF_rep3","siNC_input_rep1","siNC_input_rep2","siNC_input_rep3","siNC_IP_CTCF_rep1","siNC_IP_CTCF_rep2","siNC_IP_CTCF_rep3"}{if($1==$3 && $1==$5 && $1==$7 && $1==$9 && $1==$11 && $1==$13 && $1==$15 && $1==$17 && $1==$19 && $1==$21 && $1==$23 && $1==$25 && $1==$27 && $1==$29 && $1==$31 && $1==$33 && $1==$35 && $1==$37 && $1==$39 && $1==$41 && $1==$43 && $1==$45 && $1==$47 && $1==$49 && $1==$51 && $1==$53 && $1==$55 && $1==$57 && $1==$59 && $1==$61 && $1==$63 && $1==$65 && $1==$67 && $1==$69 && $1==$71){print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36,$38,$40,$42,$44,$46,$48,$50,$52,$54,$56,$58,$60,$62,$64,$66,$68,$70,$72}}' >all_samples_chr_mapped_uniq_reads_stat.summary.txt

```

```bash
# 绘制reads 分布热图时要用到 rep123合并的bam文件
# 创建 bam merge 脚本 4_make_all_WN_siNURD_CTCF_ChIP_Seq_rep123_merged_bam_and_index.sh 内容如下：
# # 合并各重复的 bam
# samtools merge --threads 72 HEK293T_siCHD3_input_rep123.merged.sorted.bam HEK293T_siCHD3_input_rep1.bowtie2.sorted.bam HEK293T_siCHD3_input_rep2.bowtie2.sorted.bam HEK293T_siCHD3_input_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_siCHD3_IP_CTCF_rep123.merged.sorted.bam HEK293T_siCHD3_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD3_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD3_IP_CTCF_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_siCHD4_input_rep123.merged.sorted.bam HEK293T_siCHD4_input_rep1.bowtie2.sorted.bam HEK293T_siCHD4_input_rep2.bowtie2.sorted.bam HEK293T_siCHD4_input_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_siCHD4_IP_CTCF_rep123.merged.sorted.bam HEK293T_siCHD4_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD4_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD4_IP_CTCF_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_siCHD5_input_rep123.merged.sorted.bam HEK293T_siCHD5_input_rep1.bowtie2.sorted.bam HEK293T_siCHD5_input_rep2.bowtie2.sorted.bam HEK293T_siCHD5_input_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_siCHD5_IP_CTCF_rep123.merged.sorted.bam HEK293T_siCHD5_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD5_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD5_IP_CTCF_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_siMBD2_input_rep123.merged.sorted.bam HEK293T_siMBD2_input_rep1.bowtie2.sorted.bam HEK293T_siMBD2_input_rep2.bowtie2.sorted.bam HEK293T_siMBD2_input_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_siMBD2_IP_CTCF_rep123.merged.sorted.bam HEK293T_siMBD2_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siMBD2_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siMBD2_IP_CTCF_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_siMBD3_input_rep123.merged.sorted.bam HEK293T_siMBD3_input_rep1.bowtie2.sorted.bam HEK293T_siMBD3_input_rep2.bowtie2.sorted.bam HEK293T_siMBD3_input_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_siMBD3_IP_CTCF_rep123.merged.sorted.bam HEK293T_siMBD3_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siMBD3_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siMBD3_IP_CTCF_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_siNC_input_rep123.merged.sorted.bam HEK293T_siNC_input_rep1.bowtie2.sorted.bam HEK293T_siNC_input_rep2.bowtie2.sorted.bam HEK293T_siNC_input_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam
# # 为所有 merged bam 建立索引
# for i in `ls -1 *_rep123.merged.sorted.bam`; do echo $i; samtools index -@ 72 $i;done
nohup bash 4_make_all_WN_siNURD_CTCF_ChIP_Seq_rep123_merged_bam_and_index.sh >4_make_all_WN_siNURD_CTCF_ChIP_Seq_rep123_merged_bam_and_index.sh.runlog 2>&1 &
tail -f 4_make_all_WN_siNURD_CTCF_ChIP_Seq_rep123_merged_bam_and_index.sh.runlog

```

# Detection of enrichment and repeatability

```bash
# 参考： https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html
# 可以判断实验的成败、富集程度的好坏等
# siCHD3
plotFingerprint -b HEK293T_siCHD3_input_rep1.bowtie2.sorted.bam HEK293T_siCHD3_input_rep2.bowtie2.sorted.bam HEK293T_siCHD3_input_rep3.bowtie2.sorted.bam HEK293T_siCHD3_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD3_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD3_IP_CTCF_rep3.bowtie2.sorted.bam --labels siCHD3_input_rep1 siCHD3_input_rep2 siCHD3_input_rep3 siCHD3_IP_CTCF_rep1 siCHD3_IP_CTCF_rep2 siCHD3_IP_CTCF_rep3 -o plotFingerprint_for_siCHD3_CTCF_ChIP_Seq_samples.pdf --extendReads --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max"
# siCHD4
plotFingerprint -b HEK293T_siCHD4_input_rep1.bowtie2.sorted.bam HEK293T_siCHD4_input_rep2.bowtie2.sorted.bam HEK293T_siCHD4_input_rep3.bowtie2.sorted.bam HEK293T_siCHD4_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD4_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD4_IP_CTCF_rep3.bowtie2.sorted.bam --labels siCHD4_input_rep1 siCHD4_input_rep2 siCHD4_input_rep3 siCHD4_IP_CTCF_rep1 siCHD4_IP_CTCF_rep2 siCHD4_IP_CTCF_rep3 -o plotFingerprint_for_siCHD4_CTCF_ChIP_Seq_samples.pdf --extendReads --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max"
# siCHD5
plotFingerprint -b HEK293T_siCHD5_input_rep1.bowtie2.sorted.bam HEK293T_siCHD5_input_rep2.bowtie2.sorted.bam HEK293T_siCHD5_input_rep3.bowtie2.sorted.bam HEK293T_siCHD5_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD5_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD5_IP_CTCF_rep3.bowtie2.sorted.bam --labels siCHD5_input_rep1 siCHD5_input_rep2 siCHD5_input_rep3 siCHD5_IP_CTCF_rep1 siCHD5_IP_CTCF_rep2 siCHD5_IP_CTCF_rep3 -o plotFingerprint_for_siCHD5_CTCF_ChIP_Seq_samples.pdf --extendReads --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max"
# siMBD2
plotFingerprint -b HEK293T_siMBD2_input_rep1.bowtie2.sorted.bam HEK293T_siMBD2_input_rep2.bowtie2.sorted.bam HEK293T_siMBD2_input_rep3.bowtie2.sorted.bam HEK293T_siMBD2_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siMBD2_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siMBD2_IP_CTCF_rep3.bowtie2.sorted.bam --labels siMBD2_input_rep1 siMBD2_input_rep2 siMBD2_input_rep3 siMBD2_IP_CTCF_rep1 siMBD2_IP_CTCF_rep2 siMBD2_IP_CTCF_rep3 -o plotFingerprint_for_siMBD2_CTCF_ChIP_Seq_samples.pdf --extendReads --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max"
# siMBD3
plotFingerprint -b HEK293T_siMBD3_input_rep1.bowtie2.sorted.bam HEK293T_siMBD3_input_rep2.bowtie2.sorted.bam HEK293T_siMBD3_input_rep3.bowtie2.sorted.bam HEK293T_siMBD3_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siMBD3_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siMBD3_IP_CTCF_rep3.bowtie2.sorted.bam --labels siMBD3_input_rep1 siMBD3_input_rep2 siMBD3_input_rep3 siMBD3_IP_CTCF_rep1 siMBD3_IP_CTCF_rep2 siMBD3_IP_CTCF_rep3 -o plotFingerprint_for_siMBD3_CTCF_ChIP_Seq_samples.pdf --extendReads --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max"
# siNC
plotFingerprint -b HEK293T_siNC_input_rep1.bowtie2.sorted.bam HEK293T_siNC_input_rep2.bowtie2.sorted.bam HEK293T_siNC_input_rep3.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam --labels siNC_input_rep1 siNC_input_rep2 siNC_input_rep3 siNC_IP_CTCF_rep1 siNC_IP_CTCF_rep2 siNC_IP_CTCF_rep3 -o plotFingerprint_for_siNC_CTCF_ChIP_Seq_samples.pdf --extendReads --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max"

```

```bash
# 使用 deeptools 分析 样本相似度
# 先对全基因组分bin统计 uniq reads数量
multiBamSummary bins --bamfiles HEK293T_siCHD3_input_rep1.bowtie2.sorted.bam HEK293T_siCHD3_input_rep2.bowtie2.sorted.bam HEK293T_siCHD3_input_rep3.bowtie2.sorted.bam HEK293T_siCHD3_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD3_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD3_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siCHD4_input_rep1.bowtie2.sorted.bam HEK293T_siCHD4_input_rep2.bowtie2.sorted.bam HEK293T_siCHD4_input_rep3.bowtie2.sorted.bam HEK293T_siCHD4_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD4_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD4_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siCHD5_input_rep1.bowtie2.sorted.bam HEK293T_siCHD5_input_rep2.bowtie2.sorted.bam HEK293T_siCHD5_input_rep3.bowtie2.sorted.bam HEK293T_siCHD5_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD5_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD5_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siMBD2_input_rep1.bowtie2.sorted.bam HEK293T_siMBD2_input_rep2.bowtie2.sorted.bam HEK293T_siMBD2_input_rep3.bowtie2.sorted.bam HEK293T_siMBD2_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siMBD2_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siMBD2_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siMBD3_input_rep1.bowtie2.sorted.bam HEK293T_siMBD3_input_rep2.bowtie2.sorted.bam HEK293T_siMBD3_input_rep3.bowtie2.sorted.bam HEK293T_siMBD3_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siMBD3_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siMBD3_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siNC_input_rep1.bowtie2.sorted.bam HEK293T_siNC_input_rep2.bowtie2.sorted.bam HEK293T_siNC_input_rep3.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam -out multiBamSummary_results_for_all_bam_files.npz --labels siCHD3_input_rep1 siCHD3_input_rep2 siCHD3_input_rep3 siCHD3_IP_CTCF_rep1 siCHD3_IP_CTCF_rep2 siCHD3_IP_CTCF_rep3 siCHD4_input_rep1 siCHD4_input_rep2 siCHD4_input_rep3 siCHD4_IP_CTCF_rep1 siCHD4_IP_CTCF_rep2 siCHD4_IP_CTCF_rep3 siCHD5_input_rep1 siCHD5_input_rep2 siCHD5_input_rep3 siCHD5_IP_CTCF_rep1 siCHD5_IP_CTCF_rep2 siCHD5_IP_CTCF_rep3 siMBD2_input_rep1 siMBD2_input_rep2 siMBD2_input_rep3 siMBD2_IP_CTCF_rep1 siMBD2_IP_CTCF_rep2 siMBD2_IP_CTCF_rep3 siMBD3_input_rep1 siMBD3_input_rep2 siMBD3_input_rep3 siMBD3_IP_CTCF_rep1 siMBD3_IP_CTCF_rep2 siMBD3_IP_CTCF_rep3 siNC_input_rep1 siNC_input_rep2 siNC_input_rep3 siNC_IP_CTCF_rep1 siNC_IP_CTCF_rep2 siNC_IP_CTCF_rep3 --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# 用 plotCorrelation 和 plotPCA 分析样本的相似度
plotCorrelation --corData multiBamSummary_results_for_all_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_all_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_all_samples.txt --colorMap coolwarm --plotNumbers
plotPCA --corData multiBamSummary_results_for_all_bam_files.npz --plotFile plotPCA_plot_for_all_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_all_samples.txt
plotPCA --rowCenter --corData multiBamSummary_results_for_all_bam_files.npz --plotFile plotPCA_plot_for_all_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_all_samples.rowCenter.txt

# siCHD3 ChIP_Seq 样本的相似度
multiBamSummary bins --bamfiles HEK293T_siCHD3_input_rep1.bowtie2.sorted.bam HEK293T_siCHD3_input_rep2.bowtie2.sorted.bam HEK293T_siCHD3_input_rep3.bowtie2.sorted.bam HEK293T_siCHD3_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD3_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD3_IP_CTCF_rep3.bowtie2.sorted.bam -out multiBamSummary_results_for_siCHD3_CTCF_ChIP_Seq_bam_files.npz --labels siCHD3_input_rep1 siCHD3_input_rep2 siCHD3_input_rep3 siCHD3_IP_CTCF_rep1 siCHD3_IP_CTCF_rep2 siCHD3_IP_CTCF_rep3 --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
plotCorrelation --corData multiBamSummary_results_for_siCHD3_CTCF_ChIP_Seq_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_siCHD3_CTCF_ChIP_Seq_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_siCHD3_CTCF_ChIP_Seq_samples.txt --colorMap coolwarm --plotNumbers
plotPCA --corData multiBamSummary_results_for_siCHD3_CTCF_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_siCHD3_CTCF_ChIP_Seq_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_siCHD3_CTCF_ChIP_Seq_samples.txt
plotPCA --rowCenter --corData multiBamSummary_results_for_siCHD3_CTCF_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_siCHD3_CTCF_ChIP_Seq_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_siCHD3_CTCF_ChIP_Seq_samples.rowCenter.txt
# siCHD4 ChIP_Seq 样本的相似度
multiBamSummary bins --bamfiles HEK293T_siCHD4_input_rep1.bowtie2.sorted.bam HEK293T_siCHD4_input_rep2.bowtie2.sorted.bam HEK293T_siCHD4_input_rep3.bowtie2.sorted.bam HEK293T_siCHD4_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD4_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD4_IP_CTCF_rep3.bowtie2.sorted.bam -out multiBamSummary_results_for_siCHD4_CTCF_ChIP_Seq_bam_files.npz --labels siCHD4_input_rep1 siCHD4_input_rep2 siCHD4_input_rep3 siCHD4_IP_CTCF_rep1 siCHD4_IP_CTCF_rep2 siCHD4_IP_CTCF_rep3 --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
plotCorrelation --corData multiBamSummary_results_for_siCHD4_CTCF_ChIP_Seq_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_siCHD4_CTCF_ChIP_Seq_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_siCHD4_CTCF_ChIP_Seq_samples.txt --colorMap coolwarm --plotNumbers
plotPCA --corData multiBamSummary_results_for_siCHD4_CTCF_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_siCHD4_CTCF_ChIP_Seq_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_siCHD4_CTCF_ChIP_Seq_samples.txt
plotPCA --rowCenter --corData multiBamSummary_results_for_siCHD4_CTCF_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_siCHD4_CTCF_ChIP_Seq_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_siCHD4_CTCF_ChIP_Seq_samples.rowCenter.txt
# siCHD5 ChIP_Seq 样本的相似度
multiBamSummary bins --bamfiles HEK293T_siCHD5_input_rep1.bowtie2.sorted.bam HEK293T_siCHD5_input_rep2.bowtie2.sorted.bam HEK293T_siCHD5_input_rep3.bowtie2.sorted.bam HEK293T_siCHD5_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD5_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD5_IP_CTCF_rep3.bowtie2.sorted.bam -out multiBamSummary_results_for_siCHD5_CTCF_ChIP_Seq_bam_files.npz --labels siCHD5_input_rep1 siCHD5_input_rep2 siCHD5_input_rep3 siCHD5_IP_CTCF_rep1 siCHD5_IP_CTCF_rep2 siCHD5_IP_CTCF_rep3 --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
plotCorrelation --corData multiBamSummary_results_for_siCHD5_CTCF_ChIP_Seq_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_siCHD5_CTCF_ChIP_Seq_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_siCHD5_CTCF_ChIP_Seq_samples.txt --colorMap coolwarm --plotNumbers
plotPCA --corData multiBamSummary_results_for_siCHD5_CTCF_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_siCHD5_CTCF_ChIP_Seq_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_siCHD5_CTCF_ChIP_Seq_samples.txt
plotPCA --rowCenter --corData multiBamSummary_results_for_siCHD5_CTCF_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_siCHD5_CTCF_ChIP_Seq_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_siCHD5_CTCF_ChIP_Seq_samples.rowCenter.txt
# siMBD2 ChIP_Seq 样本的相似度
multiBamSummary bins --bamfiles HEK293T_siMBD2_input_rep1.bowtie2.sorted.bam HEK293T_siMBD2_input_rep2.bowtie2.sorted.bam HEK293T_siMBD2_input_rep3.bowtie2.sorted.bam HEK293T_siMBD2_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siMBD2_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siMBD2_IP_CTCF_rep3.bowtie2.sorted.bam -out multiBamSummary_results_for_siMBD2_CTCF_ChIP_Seq_bam_files.npz --labels siMBD2_input_rep1 siMBD2_input_rep2 siMBD2_input_rep3 siMBD2_IP_CTCF_rep1 siMBD2_IP_CTCF_rep2 siMBD2_IP_CTCF_rep3 --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
plotCorrelation --corData multiBamSummary_results_for_siMBD2_CTCF_ChIP_Seq_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_siMBD2_CTCF_ChIP_Seq_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_siMBD2_CTCF_ChIP_Seq_samples.txt --colorMap coolwarm --plotNumbers
plotPCA --corData multiBamSummary_results_for_siMBD2_CTCF_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_siMBD2_CTCF_ChIP_Seq_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_siMBD2_CTCF_ChIP_Seq_samples.txt
plotPCA --rowCenter --corData multiBamSummary_results_for_siMBD2_CTCF_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_siMBD2_CTCF_ChIP_Seq_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_siMBD2_CTCF_ChIP_Seq_samples.rowCenter.txt
# siMBD3 ChIP_Seq 样本的相似度
multiBamSummary bins --bamfiles HEK293T_siMBD3_input_rep1.bowtie2.sorted.bam HEK293T_siMBD3_input_rep2.bowtie2.sorted.bam HEK293T_siMBD3_input_rep3.bowtie2.sorted.bam HEK293T_siMBD3_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siMBD3_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siMBD3_IP_CTCF_rep3.bowtie2.sorted.bam -out multiBamSummary_results_for_siMBD3_CTCF_ChIP_Seq_bam_files.npz --labels siMBD3_input_rep1 siMBD3_input_rep2 siMBD3_input_rep3 siMBD3_IP_CTCF_rep1 siMBD3_IP_CTCF_rep2 siMBD3_IP_CTCF_rep3 --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
plotCorrelation --corData multiBamSummary_results_for_siMBD3_CTCF_ChIP_Seq_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_siMBD3_CTCF_ChIP_Seq_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_siMBD3_CTCF_ChIP_Seq_samples.txt --colorMap coolwarm --plotNumbers
plotPCA --corData multiBamSummary_results_for_siMBD3_CTCF_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_siMBD3_CTCF_ChIP_Seq_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_siMBD3_CTCF_ChIP_Seq_samples.txt
plotPCA --rowCenter --corData multiBamSummary_results_for_siMBD3_CTCF_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_siMBD3_CTCF_ChIP_Seq_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_siMBD3_CTCF_ChIP_Seq_samples.rowCenter.txt
# siNC ChIP_Seq 样本的相似度
multiBamSummary bins --bamfiles HEK293T_siNC_input_rep1.bowtie2.sorted.bam HEK293T_siNC_input_rep2.bowtie2.sorted.bam HEK293T_siNC_input_rep3.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam -out multiBamSummary_results_for_siNC_CTCF_ChIP_Seq_bam_files.npz --labels siNC_input_rep1 siNC_input_rep2 siNC_input_rep3 siNC_IP_CTCF_rep1 siNC_IP_CTCF_rep2 siNC_IP_CTCF_rep3 --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
plotCorrelation --corData multiBamSummary_results_for_siNC_CTCF_ChIP_Seq_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_siNC_CTCF_ChIP_Seq_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_siNC_CTCF_ChIP_Seq_samples.txt --colorMap coolwarm --plotNumbers
plotPCA --corData multiBamSummary_results_for_siNC_CTCF_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_siNC_CTCF_ChIP_Seq_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_siNC_CTCF_ChIP_Seq_samples.txt
plotPCA --rowCenter --corData multiBamSummary_results_for_siNC_CTCF_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_siNC_CTCF_ChIP_Seq_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_siNC_CTCF_ChIP_Seq_samples.rowCenter.txt

```

# macs2 peak calling

## peak calling for each repetition separately 

### siCHD3 CTCF ChIP-Seq peak calling

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq
mkdir macs2_peak_calling_1_for_siCHD3_CTCF_each_repeat
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_1_for_siCHD3_CTCF_each_repeat
# 创建所需的bam文件连接
# bam
ln -s ../HEK293T_siCHD3_input_rep1.bowtie2.sorted.bam HEK293T_siCHD3_input_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD3_input_rep2.bowtie2.sorted.bam HEK293T_siCHD3_input_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD3_input_rep3.bowtie2.sorted.bam HEK293T_siCHD3_input_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD3_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD3_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD3_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD3_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD3_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siCHD3_IP_CTCF_rep3.bowtie2.sorted.bam
# bam.bai
ln -s ../HEK293T_siCHD3_input_rep1.bowtie2.sorted.bam.bai HEK293T_siCHD3_input_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD3_input_rep2.bowtie2.sorted.bam.bai HEK293T_siCHD3_input_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD3_input_rep3.bowtie2.sorted.bam.bai HEK293T_siCHD3_input_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD3_IP_CTCF_rep1.bowtie2.sorted.bam.bai HEK293T_siCHD3_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD3_IP_CTCF_rep2.bowtie2.sorted.bam.bai HEK293T_siCHD3_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD3_IP_CTCF_rep3.bowtie2.sorted.bam.bai HEK293T_siCHD3_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *.bam

# siCHD3
nohup macs2 callpeak -t HEK293T_siCHD3_IP_CTCF_rep1.bowtie2.sorted.bam -c HEK293T_siCHD3_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siCHD3_IP_CTCF_rep1_TvsC -n siCHD3_IP_CTCF_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_siCHD3_IP_CTCF_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_siCHD3_IP_CTCF_rep2.bowtie2.sorted.bam -c HEK293T_siCHD3_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siCHD3_IP_CTCF_rep2_TvsC -n siCHD3_IP_CTCF_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_siCHD3_IP_CTCF_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_siCHD3_IP_CTCF_rep3.bowtie2.sorted.bam -c HEK293T_siCHD3_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siCHD3_IP_CTCF_rep3_TvsC -n siCHD3_IP_CTCF_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_siCHD3_IP_CTCF_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_siCHD3_IP_CTCF_rep1_TvsC.runlog
tail -f macs2_callpeak_for_siCHD3_IP_CTCF_rep2_TvsC.runlog
tail -f macs2_callpeak_for_siCHD3_IP_CTCF_rep3_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_1_for_siCHD3_CTCF_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siCHD3_IP_CTCF_rep1_TvsC/siCHD3_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# 14323
# ./siCHD3_IP_CTCF_rep2_TvsC/siCHD3_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# 13775
# ./siCHD3_IP_CTCF_rep3_TvsC/siCHD3_IP_CTCF_rep3_TvsC_peaks.narrowPeak
# 17246

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siCHD3_IP_CTCF_rep1_TvsC/siCHD3_IP_CTCF_rep1_TvsC_summits.bed
# 14988
# ./siCHD3_IP_CTCF_rep2_TvsC/siCHD3_IP_CTCF_rep2_TvsC_summits.bed
# 14456
# ./siCHD3_IP_CTCF_rep3_TvsC/siCHD3_IP_CTCF_rep3_TvsC_summits.bed
# 17979

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_siCHD3_IP_CTCF_ChIP_Seq_samples_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bigwig|sort`;do echo $i;done

```



### siCHD4 CTCF ChIP-Seq peak calling

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq
mkdir macs2_peak_calling_1_for_siCHD4_CTCF_each_repeat
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_1_for_siCHD4_CTCF_each_repeat
# 创建所需的bam文件连接
# bam
ln -s ../HEK293T_siCHD4_input_rep1.bowtie2.sorted.bam HEK293T_siCHD4_input_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD4_input_rep2.bowtie2.sorted.bam HEK293T_siCHD4_input_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD4_input_rep3.bowtie2.sorted.bam HEK293T_siCHD4_input_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD4_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD4_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD4_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD4_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD4_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siCHD4_IP_CTCF_rep3.bowtie2.sorted.bam
# bam.bai
ln -s ../HEK293T_siCHD4_input_rep1.bowtie2.sorted.bam.bai HEK293T_siCHD4_input_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD4_input_rep2.bowtie2.sorted.bam.bai HEK293T_siCHD4_input_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD4_input_rep3.bowtie2.sorted.bam.bai HEK293T_siCHD4_input_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD4_IP_CTCF_rep1.bowtie2.sorted.bam.bai HEK293T_siCHD4_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD4_IP_CTCF_rep2.bowtie2.sorted.bam.bai HEK293T_siCHD4_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD4_IP_CTCF_rep3.bowtie2.sorted.bam.bai HEK293T_siCHD4_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *.bam

# siCHD4
nohup macs2 callpeak -t HEK293T_siCHD4_IP_CTCF_rep1.bowtie2.sorted.bam -c HEK293T_siCHD4_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siCHD4_IP_CTCF_rep1_TvsC -n siCHD4_IP_CTCF_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_siCHD4_IP_CTCF_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_siCHD4_IP_CTCF_rep2.bowtie2.sorted.bam -c HEK293T_siCHD4_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siCHD4_IP_CTCF_rep2_TvsC -n siCHD4_IP_CTCF_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_siCHD4_IP_CTCF_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_siCHD4_IP_CTCF_rep3.bowtie2.sorted.bam -c HEK293T_siCHD4_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siCHD4_IP_CTCF_rep3_TvsC -n siCHD4_IP_CTCF_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_siCHD4_IP_CTCF_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_siCHD4_IP_CTCF_rep1_TvsC.runlog
tail -f macs2_callpeak_for_siCHD4_IP_CTCF_rep2_TvsC.runlog
tail -f macs2_callpeak_for_siCHD4_IP_CTCF_rep3_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_1_for_siCHD4_CTCF_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siCHD4_IP_CTCF_rep1_TvsC/siCHD4_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# 6093
# ./siCHD4_IP_CTCF_rep2_TvsC/siCHD4_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# 8888
# ./siCHD4_IP_CTCF_rep3_TvsC/siCHD4_IP_CTCF_rep3_TvsC_peaks.narrowPeak
# 9480

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siCHD4_IP_CTCF_rep1_TvsC/siCHD4_IP_CTCF_rep1_TvsC_summits.bed
# 6405
# ./siCHD4_IP_CTCF_rep2_TvsC/siCHD4_IP_CTCF_rep2_TvsC_summits.bed
# 9317
# ./siCHD4_IP_CTCF_rep3_TvsC/siCHD4_IP_CTCF_rep3_TvsC_summits.bed
# 9995

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_siCHD4_IP_CTCF_ChIP_Seq_samples_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bigwig|sort`;do echo $i;done

```



### siCHD5 CTCF ChIP-Seq peak calling

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq
mkdir macs2_peak_calling_1_for_siCHD5_CTCF_each_repeat
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_1_for_siCHD5_CTCF_each_repeat
# 创建所需的bam文件连接
# bam
ln -s ../HEK293T_siCHD5_input_rep1.bowtie2.sorted.bam HEK293T_siCHD5_input_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD5_input_rep2.bowtie2.sorted.bam HEK293T_siCHD5_input_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD5_input_rep3.bowtie2.sorted.bam HEK293T_siCHD5_input_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD5_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD5_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD5_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD5_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD5_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siCHD5_IP_CTCF_rep3.bowtie2.sorted.bam
# bam.bai
ln -s ../HEK293T_siCHD5_input_rep1.bowtie2.sorted.bam.bai HEK293T_siCHD5_input_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD5_input_rep2.bowtie2.sorted.bam.bai HEK293T_siCHD5_input_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD5_input_rep3.bowtie2.sorted.bam.bai HEK293T_siCHD5_input_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD5_IP_CTCF_rep1.bowtie2.sorted.bam.bai HEK293T_siCHD5_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD5_IP_CTCF_rep2.bowtie2.sorted.bam.bai HEK293T_siCHD5_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD5_IP_CTCF_rep3.bowtie2.sorted.bam.bai HEK293T_siCHD5_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *.bam

# siCHD5
nohup macs2 callpeak -t HEK293T_siCHD5_IP_CTCF_rep1.bowtie2.sorted.bam -c HEK293T_siCHD5_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siCHD5_IP_CTCF_rep1_TvsC -n siCHD5_IP_CTCF_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_siCHD5_IP_CTCF_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_siCHD5_IP_CTCF_rep2.bowtie2.sorted.bam -c HEK293T_siCHD5_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siCHD5_IP_CTCF_rep2_TvsC -n siCHD5_IP_CTCF_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_siCHD5_IP_CTCF_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_siCHD5_IP_CTCF_rep3.bowtie2.sorted.bam -c HEK293T_siCHD5_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siCHD5_IP_CTCF_rep3_TvsC -n siCHD5_IP_CTCF_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_siCHD5_IP_CTCF_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_siCHD5_IP_CTCF_rep1_TvsC.runlog
tail -f macs2_callpeak_for_siCHD5_IP_CTCF_rep2_TvsC.runlog
tail -f macs2_callpeak_for_siCHD5_IP_CTCF_rep3_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_1_for_siCHD5_CTCF_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siCHD5_IP_CTCF_rep1_TvsC/siCHD5_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# 18290
# ./siCHD5_IP_CTCF_rep2_TvsC/siCHD5_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# 19296
# ./siCHD5_IP_CTCF_rep3_TvsC/siCHD5_IP_CTCF_rep3_TvsC_peaks.narrowPeak
# 8274

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siCHD5_IP_CTCF_rep1_TvsC/siCHD5_IP_CTCF_rep1_TvsC_summits.bed
# 18452
# ./siCHD5_IP_CTCF_rep2_TvsC/siCHD5_IP_CTCF_rep2_TvsC_summits.bed
# 19491
# ./siCHD5_IP_CTCF_rep3_TvsC/siCHD5_IP_CTCF_rep3_TvsC_summits.bed
# 8429

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_siCHD5_IP_CTCF_ChIP_Seq_samples_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bigwig|sort`;do echo $i;done

```


### siMBD2 CTCF ChIP-Seq peak calling

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq
mkdir macs2_peak_calling_1_for_siMBD2_CTCF_each_repeat
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_1_for_siMBD2_CTCF_each_repeat
# 创建所需的bam文件连接
# bam
ln -s ../HEK293T_siMBD2_input_rep1.bowtie2.sorted.bam HEK293T_siMBD2_input_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD2_input_rep2.bowtie2.sorted.bam HEK293T_siMBD2_input_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD2_input_rep3.bowtie2.sorted.bam HEK293T_siMBD2_input_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD2_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siMBD2_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD2_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siMBD2_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD2_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siMBD2_IP_CTCF_rep3.bowtie2.sorted.bam
# bam.bai
ln -s ../HEK293T_siMBD2_input_rep1.bowtie2.sorted.bam.bai HEK293T_siMBD2_input_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD2_input_rep2.bowtie2.sorted.bam.bai HEK293T_siMBD2_input_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD2_input_rep3.bowtie2.sorted.bam.bai HEK293T_siMBD2_input_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD2_IP_CTCF_rep1.bowtie2.sorted.bam.bai HEK293T_siMBD2_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD2_IP_CTCF_rep2.bowtie2.sorted.bam.bai HEK293T_siMBD2_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD2_IP_CTCF_rep3.bowtie2.sorted.bam.bai HEK293T_siMBD2_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *.bam

# siMBD2
nohup macs2 callpeak -t HEK293T_siMBD2_IP_CTCF_rep1.bowtie2.sorted.bam -c HEK293T_siMBD2_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siMBD2_IP_CTCF_rep1_TvsC -n siMBD2_IP_CTCF_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_siMBD2_IP_CTCF_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_siMBD2_IP_CTCF_rep2.bowtie2.sorted.bam -c HEK293T_siMBD2_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siMBD2_IP_CTCF_rep2_TvsC -n siMBD2_IP_CTCF_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_siMBD2_IP_CTCF_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_siMBD2_IP_CTCF_rep3.bowtie2.sorted.bam -c HEK293T_siMBD2_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siMBD2_IP_CTCF_rep3_TvsC -n siMBD2_IP_CTCF_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_siMBD2_IP_CTCF_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_siMBD2_IP_CTCF_rep1_TvsC.runlog
tail -f macs2_callpeak_for_siMBD2_IP_CTCF_rep2_TvsC.runlog
tail -f macs2_callpeak_for_siMBD2_IP_CTCF_rep3_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_1_for_siMBD2_CTCF_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siMBD2_IP_CTCF_rep1_TvsC/siMBD2_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# 9433
# ./siMBD2_IP_CTCF_rep2_TvsC/siMBD2_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# 7748
# ./siMBD2_IP_CTCF_rep3_TvsC/siMBD2_IP_CTCF_rep3_TvsC_peaks.narrowPeak
# 12693

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siMBD2_IP_CTCF_rep1_TvsC/siMBD2_IP_CTCF_rep1_TvsC_summits.bed
# 9567
# ./siMBD2_IP_CTCF_rep2_TvsC/siMBD2_IP_CTCF_rep2_TvsC_summits.bed
# 7808
# ./siMBD2_IP_CTCF_rep3_TvsC/siMBD2_IP_CTCF_rep3_TvsC_summits.bed
# 12837

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_siMBD2_IP_CTCF_ChIP_Seq_samples_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done

```



### siMBD3 CTCF ChIP-Seq peak calling

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq
mkdir macs2_peak_calling_1_for_siMBD3_CTCF_each_repeat
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_1_for_siMBD3_CTCF_each_repeat
# 创建所需的bam文件连接
# bam
ln -s ../HEK293T_siMBD3_input_rep1.bowtie2.sorted.bam HEK293T_siMBD3_input_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD3_input_rep2.bowtie2.sorted.bam HEK293T_siMBD3_input_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD3_input_rep3.bowtie2.sorted.bam HEK293T_siMBD3_input_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD3_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siMBD3_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD3_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siMBD3_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD3_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siMBD3_IP_CTCF_rep3.bowtie2.sorted.bam
# bam.bai
ln -s ../HEK293T_siMBD3_input_rep1.bowtie2.sorted.bam.bai HEK293T_siMBD3_input_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD3_input_rep2.bowtie2.sorted.bam.bai HEK293T_siMBD3_input_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD3_input_rep3.bowtie2.sorted.bam.bai HEK293T_siMBD3_input_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD3_IP_CTCF_rep1.bowtie2.sorted.bam.bai HEK293T_siMBD3_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD3_IP_CTCF_rep2.bowtie2.sorted.bam.bai HEK293T_siMBD3_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD3_IP_CTCF_rep3.bowtie2.sorted.bam.bai HEK293T_siMBD3_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *.bam

# siMBD3
nohup macs2 callpeak -t HEK293T_siMBD3_IP_CTCF_rep1.bowtie2.sorted.bam -c HEK293T_siMBD3_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siMBD3_IP_CTCF_rep1_TvsC -n siMBD3_IP_CTCF_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_siMBD3_IP_CTCF_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_siMBD3_IP_CTCF_rep2.bowtie2.sorted.bam -c HEK293T_siMBD3_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siMBD3_IP_CTCF_rep2_TvsC -n siMBD3_IP_CTCF_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_siMBD3_IP_CTCF_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_siMBD3_IP_CTCF_rep3.bowtie2.sorted.bam -c HEK293T_siMBD3_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siMBD3_IP_CTCF_rep3_TvsC -n siMBD3_IP_CTCF_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_siMBD3_IP_CTCF_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_siMBD3_IP_CTCF_rep1_TvsC.runlog
tail -f macs2_callpeak_for_siMBD3_IP_CTCF_rep2_TvsC.runlog
tail -f macs2_callpeak_for_siMBD3_IP_CTCF_rep3_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_1_for_siMBD3_CTCF_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siMBD3_IP_CTCF_rep1_TvsC/siMBD3_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# 10513
# ./siMBD3_IP_CTCF_rep2_TvsC/siMBD3_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# 5170
# ./siMBD3_IP_CTCF_rep3_TvsC/siMBD3_IP_CTCF_rep3_TvsC_peaks.narrowPeak
# 11468

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siMBD3_IP_CTCF_rep1_TvsC/siMBD3_IP_CTCF_rep1_TvsC_summits.bed
# 10631
# ./siMBD3_IP_CTCF_rep2_TvsC/siMBD3_IP_CTCF_rep2_TvsC_summits.bed
# 5231
# ./siMBD3_IP_CTCF_rep3_TvsC/siMBD3_IP_CTCF_rep3_TvsC_summits.bed
# 11605

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_siMBD3_IP_CTCF_ChIP_Seq_samples_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bigwig|sort`;do echo $i;done

```



### siNC CTCF ChIP-Seq peak calling

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq
mkdir macs2_peak_calling_1_for_siNC_CTCF_each_repeat
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_1_for_siNC_CTCF_each_repeat
# 创建所需的bam文件连接
# bam
ln -s ../HEK293T_siNC_input_rep1.bowtie2.sorted.bam HEK293T_siNC_input_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_input_rep2.bowtie2.sorted.bam HEK293T_siNC_input_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_input_rep3.bowtie2.sorted.bam HEK293T_siNC_input_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam
# bam.bai
ln -s ../HEK293T_siNC_input_rep1.bowtie2.sorted.bam.bai HEK293T_siNC_input_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_input_rep2.bowtie2.sorted.bam.bai HEK293T_siNC_input_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_input_rep3.bowtie2.sorted.bam.bai HEK293T_siNC_input_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *.bam

# siNC
nohup macs2 callpeak -t HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam -c HEK293T_siNC_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siNC_IP_CTCF_rep1_TvsC -n siNC_IP_CTCF_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_siNC_IP_CTCF_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam -c HEK293T_siNC_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siNC_IP_CTCF_rep2_TvsC -n siNC_IP_CTCF_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_siNC_IP_CTCF_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam -c HEK293T_siNC_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siNC_IP_CTCF_rep3_TvsC -n siNC_IP_CTCF_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_siNC_IP_CTCF_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_siNC_IP_CTCF_rep1_TvsC.runlog
tail -f macs2_callpeak_for_siNC_IP_CTCF_rep2_TvsC.runlog
tail -f macs2_callpeak_for_siNC_IP_CTCF_rep3_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_1_for_siNC_CTCF_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siNC_IP_CTCF_rep1_TvsC/siNC_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# 20819
# ./siNC_IP_CTCF_rep2_TvsC/siNC_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# 23515
# ./siNC_IP_CTCF_rep3_TvsC/siNC_IP_CTCF_rep3_TvsC_peaks.narrowPeak
# 22466

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siNC_IP_CTCF_rep1_TvsC/siNC_IP_CTCF_rep1_TvsC_summits.bed
# 21231
# ./siNC_IP_CTCF_rep2_TvsC/siNC_IP_CTCF_rep2_TvsC_summits.bed
# 24005
# ./siNC_IP_CTCF_rep3_TvsC/siNC_IP_CTCF_rep3_TvsC_summits.bed
# 22879

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_siNC_IP_CTCF_ChIP_Seq_samples_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bigwig|sort`;do echo $i;done

```



## peak calling for rep123 merged

### siCHD3 rep123 merged peak calling

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq
mkdir macs2_peak_calling_2_for_siCHD3_rep123_together
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siCHD3_rep123_together
# 创建所需的bam文件连接
# bam
ln -s ../HEK293T_siCHD3_input_rep1.bowtie2.sorted.bam HEK293T_siCHD3_input_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD3_input_rep2.bowtie2.sorted.bam HEK293T_siCHD3_input_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD3_input_rep3.bowtie2.sorted.bam HEK293T_siCHD3_input_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD3_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD3_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD3_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD3_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD3_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siCHD3_IP_CTCF_rep3.bowtie2.sorted.bam
# bam.bai
ln -s ../HEK293T_siCHD3_input_rep1.bowtie2.sorted.bam.bai HEK293T_siCHD3_input_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD3_input_rep2.bowtie2.sorted.bam.bai HEK293T_siCHD3_input_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD3_input_rep3.bowtie2.sorted.bam.bai HEK293T_siCHD3_input_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD3_IP_CTCF_rep1.bowtie2.sorted.bam.bai HEK293T_siCHD3_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD3_IP_CTCF_rep2.bowtie2.sorted.bam.bai HEK293T_siCHD3_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD3_IP_CTCF_rep3.bowtie2.sorted.bam.bai HEK293T_siCHD3_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *.bam

# siCHD3_IP_CTCF_rep123_TvsC
nohup macs2 callpeak -t HEK293T_siCHD3_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD3_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD3_IP_CTCF_rep3.bowtie2.sorted.bam -c HEK293T_siCHD3_input_rep1.bowtie2.sorted.bam HEK293T_siCHD3_input_rep2.bowtie2.sorted.bam HEK293T_siCHD3_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siCHD3_IP_CTCF_rep123_TvsC -n siCHD3_IP_CTCF_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_siCHD3_IP_CTCF_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_siCHD3_IP_CTCF_rep123_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siCHD3_rep123_together
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siCHD3_IP_CTCF_rep123_TvsC/siCHD3_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# 35356

# 统计 peak summit 数量
for i in `find . -name *summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siCHD3_IP_CTCF_rep123_TvsC/siCHD3_IP_CTCF_rep123_TvsC_summits.bed
# 37785

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_siCHD3_IP_CTCF_ChIP_Seq_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siCHD3_rep123_together
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

```
#### Intersection with the siNC

```bash
# 两组结果交集
# cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together
# for i in `find . -name *narrowPeak`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
# wc -l `find . -name *_rep123_TvsC_peaks.narrowPeak.bed`
# # 36855 ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
# bedtools slop -b 50 -i ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
# wc -l `find . -name *summits*.bed`
# # 38169 ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.bed
# # 38169 ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed

cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siCHD3_rep123_together
for i in `find . -name *narrowPeak`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
wc -l `find . -name *_rep123_TvsC_peaks.narrowPeak.bed`
# 35356 ./siCHD3_IP_CTCF_rep123_TvsC/siCHD3_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
bedtools slop -b 50 -i ./siCHD3_IP_CTCF_rep123_TvsC/siCHD3_IP_CTCF_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./siCHD3_IP_CTCF_rep123_TvsC/siCHD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
wc -l `find . -name *summits*.bed`
# 37785 ./siCHD3_IP_CTCF_rep123_TvsC/siCHD3_IP_CTCF_rep123_TvsC_summits.bed
# 37785 ./siCHD3_IP_CTCF_rep123_TvsC/siCHD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed

# 创建 siNC narrowPeak bed 以及 summits extB50 bed 文件连接
ln -s ../macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
ln -s ../macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
# 创建 siCHD3 narrowPeak bed 以及 summits extB50 bed 文件连接
ln -s ./siCHD3_IP_CTCF_rep123_TvsC/siCHD3_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siCHD3_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
ln -s ./siCHD3_IP_CTCF_rep123_TvsC/siCHD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siCHD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed

# 两组 narrowPeak 交集数量
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -b siCHD3_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 23922
bedtools intersect -a siCHD3_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -b siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 23607

# 两组 summit extB50 交集数量
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siCHD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 24107
bedtools intersect -a siCHD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 24101

# 构造 三类 summit
# common summit
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siCHD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siCHD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >IP_CTCF_rep123_TvsC_siNC_and_siCHD3.common.summit.bed
# siCHD3 only summit
bedtools intersect -a siCHD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a siCHD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >IP_CTCF_rep123_TvsC_siCHD3.only.summit.bed
# siNC only summit
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siCHD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siCHD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >IP_CTCF_rep123_TvsC_siNC.only.summit.bed
wc -l *summit.bed
# 13684 IP_CTCF_rep123_TvsC_siCHD3.only.summit.bed
# 24110 IP_CTCF_rep123_TvsC_siNC_and_siCHD3.common.summit.bed
# 14062 IP_CTCF_rep123_TvsC_siNC.only.summit.bed


# intervene plot
# siNC and siMBD3 peaks venn
intervene venn -i siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siCHD3_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed --names=siNC_CTCF_peaks,siCHD3_CTCF_peaks -o 0_siNC_siMBD3_CTCF_peaks_venn --save-overlaps
# siNC and siMBD3 summit extB50 venn
intervene venn -i siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siCHD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed --names=siNC_CTCF_summits_eB50,siCHD3_CTCF_summits_eB50 -o 0_siNC_siMBD3_CTCF_summits_extB50_venn --save-overlaps

# jupyter make scaled venn plot based on number 
# 0_make_siNC_siCHD3_scaled_venn_plot.based_on_number.ipynb

```



### siCHD4 rep123 merged peak calling

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq
mkdir macs2_peak_calling_2_for_siCHD4_rep123_together
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siCHD4_rep123_together
# 创建所需的bam文件连接
# bam
ln -s ../HEK293T_siCHD4_input_rep1.bowtie2.sorted.bam HEK293T_siCHD4_input_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD4_input_rep2.bowtie2.sorted.bam HEK293T_siCHD4_input_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD4_input_rep3.bowtie2.sorted.bam HEK293T_siCHD4_input_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD4_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD4_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD4_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD4_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD4_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siCHD4_IP_CTCF_rep3.bowtie2.sorted.bam
# bam.bai
ln -s ../HEK293T_siCHD4_input_rep1.bowtie2.sorted.bam.bai HEK293T_siCHD4_input_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD4_input_rep2.bowtie2.sorted.bam.bai HEK293T_siCHD4_input_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD4_input_rep3.bowtie2.sorted.bam.bai HEK293T_siCHD4_input_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD4_IP_CTCF_rep1.bowtie2.sorted.bam.bai HEK293T_siCHD4_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD4_IP_CTCF_rep2.bowtie2.sorted.bam.bai HEK293T_siCHD4_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD4_IP_CTCF_rep3.bowtie2.sorted.bam.bai HEK293T_siCHD4_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *.bam

# siCHD4_IP_CTCF_rep123_TvsC
nohup macs2 callpeak -t HEK293T_siCHD4_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD4_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD4_IP_CTCF_rep3.bowtie2.sorted.bam -c HEK293T_siCHD4_input_rep1.bowtie2.sorted.bam HEK293T_siCHD4_input_rep2.bowtie2.sorted.bam HEK293T_siCHD4_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siCHD4_IP_CTCF_rep123_TvsC -n siCHD4_IP_CTCF_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_siCHD4_IP_CTCF_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_siCHD4_IP_CTCF_rep123_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siCHD4_rep123_together
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siCHD4_IP_CTCF_rep123_TvsC/siCHD4_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# 26017

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siCHD4_IP_CTCF_rep123_TvsC/siCHD4_IP_CTCF_rep123_TvsC_summits.bed
# 27547

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_siCHD4_IP_CTCF_ChIP_Seq_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siCHD4_rep123_together
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

```

#### Intersection with the siNC

```bash
# 两组结果交集
# cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together
# for i in `find . -name *narrowPeak`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
# wc -l `find . -name *_rep123_TvsC_peaks.narrowPeak.bed`
# # 36855 ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
# bedtools slop -b 50 -i ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
# wc -l `find . -name *summits*.bed`
# # 38169 ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.bed
# # 38169 ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed

cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siCHD4_rep123_together
for i in `find . -name *narrowPeak`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
wc -l `find . -name *_rep123_TvsC_peaks.narrowPeak.bed`
# 26017 ./siCHD4_IP_CTCF_rep123_TvsC/siCHD4_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
bedtools slop -b 50 -i ./siCHD4_IP_CTCF_rep123_TvsC/siCHD4_IP_CTCF_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./siCHD4_IP_CTCF_rep123_TvsC/siCHD4_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
wc -l `find . -name *summits*.bed`
# 27547 ./siCHD4_IP_CTCF_rep123_TvsC/siCHD4_IP_CTCF_rep123_TvsC_summits.bed
# 27547 ./siCHD4_IP_CTCF_rep123_TvsC/siCHD4_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed

# 创建 siNC narrowPeak bed 以及 summits extB50 bed 文件连接
ln -s ../macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
ln -s ../macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
# 创建 siCHD4 narrowPeak bed 以及 summits extB50 bed 文件连接
ln -s ./siCHD4_IP_CTCF_rep123_TvsC/siCHD4_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siCHD4_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
ln -s ./siCHD4_IP_CTCF_rep123_TvsC/siCHD4_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siCHD4_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed

# 两组 narrowPeak 交集数量
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -b siCHD4_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 21666
bedtools intersect -a siCHD4_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -b siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 21528

# 两组 summit extB50 交集数量
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siCHD4_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 21984
bedtools intersect -a siCHD4_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 21987

# 构造 三类 summit
# common summit
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siCHD4_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siCHD4_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >IP_CTCF_rep123_TvsC_siNC_and_siCHD4.common.summit.bed
# siCHD4 only summit
bedtools intersect -a siCHD4_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a siCHD4_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >IP_CTCF_rep123_TvsC_siCHD4.only.summit.bed
# siNC only summit
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siCHD4_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siCHD4_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >IP_CTCF_rep123_TvsC_siNC.only.summit.bed
wc -l *summit.bed
#  5560 IP_CTCF_rep123_TvsC_siCHD4.only.summit.bed
# 21998 IP_CTCF_rep123_TvsC_siNC_and_siCHD4.common.summit.bed
# 16185 IP_CTCF_rep123_TvsC_siNC.only.summit.bed


# intervene plot
# siNC and siCHD4 peaks venn
intervene venn -i siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siCHD4_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed --names=siNC_CTCF_peaks,siCHD4_CTCF_peaks -o 0_siNC_siCHD4_CTCF_peaks_venn --save-overlaps
# siNC and siCHD4 summit extB50 venn
intervene venn -i siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siCHD4_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed --names=siNC_CTCF_summits_eB50,siCHD4_CTCF_summits_eB50 -o 0_siNC_siCHD4_CTCF_summits_extB50_venn --save-overlaps

# jupyter make scaled venn plot based on number 
# 0_make_siNC_siCHD4_scaled_venn_plot.based_on_number.ipynb

```

### siCHD5 rep123 merged peak calling

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq
mkdir macs2_peak_calling_2_for_siCHD5_rep123_together
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siCHD5_rep123_together
# 创建所需的bam文件连接
# bam
ln -s ../HEK293T_siCHD5_input_rep1.bowtie2.sorted.bam HEK293T_siCHD5_input_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD5_input_rep2.bowtie2.sorted.bam HEK293T_siCHD5_input_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD5_input_rep3.bowtie2.sorted.bam HEK293T_siCHD5_input_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD5_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD5_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD5_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD5_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD5_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siCHD5_IP_CTCF_rep3.bowtie2.sorted.bam
# bam.bai
ln -s ../HEK293T_siCHD5_input_rep1.bowtie2.sorted.bam.bai HEK293T_siCHD5_input_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD5_input_rep2.bowtie2.sorted.bam.bai HEK293T_siCHD5_input_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD5_input_rep3.bowtie2.sorted.bam.bai HEK293T_siCHD5_input_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD5_IP_CTCF_rep1.bowtie2.sorted.bam.bai HEK293T_siCHD5_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD5_IP_CTCF_rep2.bowtie2.sorted.bam.bai HEK293T_siCHD5_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD5_IP_CTCF_rep3.bowtie2.sorted.bam.bai HEK293T_siCHD5_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *.bam

# siCHD5_IP_CTCF_rep123_TvsC
nohup macs2 callpeak -t HEK293T_siCHD5_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD5_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD5_IP_CTCF_rep3.bowtie2.sorted.bam -c HEK293T_siCHD5_input_rep1.bowtie2.sorted.bam HEK293T_siCHD5_input_rep2.bowtie2.sorted.bam HEK293T_siCHD5_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siCHD5_IP_CTCF_rep123_TvsC -n siCHD5_IP_CTCF_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_siCHD5_IP_CTCF_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_siCHD5_IP_CTCF_rep123_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siCHD5_rep123_together
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siCHD5_IP_CTCF_rep123_TvsC/siCHD5_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# 32037

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siCHD5_IP_CTCF_rep123_TvsC/siCHD5_IP_CTCF_rep123_TvsC_summits.bed
# 32626

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_siCHD5_IP_CTCF_ChIP_Seq_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siCHD5_rep123_together
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

```
#### Intersection with the siNC

```bash
# 两组结果交集
# cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together
# for i in `find . -name *narrowPeak`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
# wc -l `find . -name *_rep123_TvsC_peaks.narrowPeak.bed`
# # 36855 ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
# bedtools slop -b 50 -i ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
# wc -l `find . -name *summits*.bed`
# # 38169 ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.bed
# # 38169 ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed

cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siCHD5_rep123_together
for i in `find . -name *narrowPeak`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
wc -l `find . -name *_rep123_TvsC_peaks.narrowPeak.bed`
# 32037 ./siCHD5_IP_CTCF_rep123_TvsC/siCHD5_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
bedtools slop -b 50 -i ./siCHD5_IP_CTCF_rep123_TvsC/siCHD5_IP_CTCF_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./siCHD5_IP_CTCF_rep123_TvsC/siCHD5_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
wc -l `find . -name *summits*.bed`
# 32626 ./siCHD5_IP_CTCF_rep123_TvsC/siCHD5_IP_CTCF_rep123_TvsC_summits.bed
# 32626 ./siCHD5_IP_CTCF_rep123_TvsC/siCHD5_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed

# 创建 siNC narrowPeak bed 以及 summits extB50 bed 文件连接
ln -s ../macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
ln -s ../macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
# 创建 siCHD5 narrowPeak bed 以及 summits extB50 bed 文件连接
ln -s ./siCHD5_IP_CTCF_rep123_TvsC/siCHD5_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siCHD5_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
ln -s ./siCHD5_IP_CTCF_rep123_TvsC/siCHD5_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siCHD5_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed

# 两组 narrowPeak 交集数量
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -b siCHD5_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 28929
bedtools intersect -a siCHD5_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -b siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 28955

# 两组 summit extB50 交集数量
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siCHD5_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 29223
bedtools intersect -a siCHD5_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 29219

# 构造 三类 summit
# common summit
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siCHD5_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siCHD5_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >IP_CTCF_rep123_TvsC_siNC_and_siCHD5.common.summit.bed
# siCHD5 only summit
bedtools intersect -a siCHD5_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a siCHD5_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >IP_CTCF_rep123_TvsC_siCHD5.only.summit.bed
# siNC only summit
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siCHD5_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siCHD5_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >IP_CTCF_rep123_TvsC_siNC.only.summit.bed
wc -l *summit.bed
#  3407 IP_CTCF_rep123_TvsC_siCHD5.only.summit.bed
# 29224 IP_CTCF_rep123_TvsC_siNC_and_siCHD5.common.summit.bed
#  8946 IP_CTCF_rep123_TvsC_siNC.only.summit.bed

# intervene plot
# siNC and siCHD5 peaks venn
intervene venn -i siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siCHD5_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed --names=siNC_CTCF_peaks,siCHD5_CTCF_peaks -o 0_siNC_siCHD5_CTCF_peaks_venn --save-overlaps
# siNC and siCHD5 summit extB50 venn
intervene venn -i siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siCHD5_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed --names=siNC_CTCF_summits_eB50,siCHD5_CTCF_summits_eB50 -o 0_siNC_siCHD5_CTCF_summits_extB50_venn --save-overlaps

# jupyter make scaled venn plot based on number 
# 0_make_siNC_siCHD5_scaled_venn_plot.based_on_number.ipynb

```

### siMBD2 rep123 merged peak calling

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq
mkdir macs2_peak_calling_2_for_siMBD2_rep123_together
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siMBD2_rep123_together
# 创建所需的bam文件连接
# bam
ln -s ../HEK293T_siMBD2_input_rep1.bowtie2.sorted.bam HEK293T_siMBD2_input_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD2_input_rep2.bowtie2.sorted.bam HEK293T_siMBD2_input_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD2_input_rep3.bowtie2.sorted.bam HEK293T_siMBD2_input_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD2_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siMBD2_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD2_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siMBD2_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD2_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siMBD2_IP_CTCF_rep3.bowtie2.sorted.bam
# bam.bai
ln -s ../HEK293T_siMBD2_input_rep1.bowtie2.sorted.bam.bai HEK293T_siMBD2_input_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD2_input_rep2.bowtie2.sorted.bam.bai HEK293T_siMBD2_input_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD2_input_rep3.bowtie2.sorted.bam.bai HEK293T_siMBD2_input_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD2_IP_CTCF_rep1.bowtie2.sorted.bam.bai HEK293T_siMBD2_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD2_IP_CTCF_rep2.bowtie2.sorted.bam.bai HEK293T_siMBD2_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD2_IP_CTCF_rep3.bowtie2.sorted.bam.bai HEK293T_siMBD2_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *.bam

# siMBD2_IP_CTCF_rep123_TvsC
nohup macs2 callpeak -t HEK293T_siMBD2_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siMBD2_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siMBD2_IP_CTCF_rep3.bowtie2.sorted.bam -c HEK293T_siMBD2_input_rep1.bowtie2.sorted.bam HEK293T_siMBD2_input_rep2.bowtie2.sorted.bam HEK293T_siMBD2_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siMBD2_IP_CTCF_rep123_TvsC -n siMBD2_IP_CTCF_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_siMBD2_IP_CTCF_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_siMBD2_IP_CTCF_rep123_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siMBD2_rep123_together
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siMBD2_IP_CTCF_rep123_TvsC/siMBD2_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# 26690

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siMBD2_IP_CTCF_rep123_TvsC/siMBD2_IP_CTCF_rep123_TvsC_summits.bed
# 27187

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_siMBD2_IP_CTCF_ChIP_Seq_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siMBD2_rep123_together
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

```
#### Intersection with the siNC

```bash
# 两组结果交集
# cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together
# for i in `find . -name *narrowPeak`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
# wc -l `find . -name *_rep123_TvsC_peaks.narrowPeak.bed`
# # 36855 ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
# bedtools slop -b 50 -i ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
# wc -l `find . -name *summits*.bed`
# # 38169 ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.bed
# # 38169 ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed

cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siMBD2_rep123_together
for i in `find . -name *narrowPeak`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
wc -l `find . -name *_rep123_TvsC_peaks.narrowPeak.bed`
# 26690 ./siMBD2_IP_CTCF_rep123_TvsC/siMBD2_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
bedtools slop -b 50 -i ./siMBD2_IP_CTCF_rep123_TvsC/siMBD2_IP_CTCF_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./siMBD2_IP_CTCF_rep123_TvsC/siMBD2_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
wc -l `find . -name *summits*.bed`
# 27187 ./siMBD2_IP_CTCF_rep123_TvsC/siMBD2_IP_CTCF_rep123_TvsC_summits.bed
# 27187 ./siMBD2_IP_CTCF_rep123_TvsC/siMBD2_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed

# 创建 siNC narrowPeak bed 以及 summits extB50 bed 文件连接
ln -s ../macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
ln -s ../macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
# 创建 siMBD2 narrowPeak bed 以及 summits extB50 bed 文件连接
ln -s ./siMBD2_IP_CTCF_rep123_TvsC/siMBD2_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siMBD2_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
ln -s ./siMBD2_IP_CTCF_rep123_TvsC/siMBD2_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siMBD2_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed

# 两组 narrowPeak 交集数量
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -b siMBD2_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 24055
bedtools intersect -a siMBD2_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -b siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 24112

# 两组 summit extB50 交集数量
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siMBD2_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 24233
bedtools intersect -a siMBD2_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 24228

# 构造 三类 summit
# common summit
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siMBD2_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siMBD2_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >IP_CTCF_rep123_TvsC_siNC_and_siMBD2.common.summit.bed
# siMBD2 only summit
bedtools intersect -a siMBD2_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a siMBD2_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >IP_CTCF_rep123_TvsC_siMBD2.only.summit.bed
# siNC only summit
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siMBD2_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siMBD2_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >IP_CTCF_rep123_TvsC_siNC.only.summit.bed
wc -l *summit.bed
#  2959 IP_CTCF_rep123_TvsC_siMBD2.only.summit.bed
# 24236 IP_CTCF_rep123_TvsC_siNC_and_siMBD2.common.summit.bed
# 13936 IP_CTCF_rep123_TvsC_siNC.only.summit.bed

# intervene plot
# siNC and siMBD2 peaks venn
intervene venn -i siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siMBD2_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed --names=siNC_CTCF_peaks,siMBD2_CTCF_peaks -o 0_siNC_siMBD2_CTCF_peaks_venn --save-overlaps
# siNC and siMBD2 summit extB50 venn
intervene venn -i siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siMBD2_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed --names=siNC_CTCF_summits_eB50,siMBD2_CTCF_summits_eB50 -o 0_siNC_siMBD2_CTCF_summits_extB50_venn --save-overlaps

# jupyter make scaled venn plot based on number 
# 0_make_siNC_siMBD2_scaled_venn_plot.based_on_number.ipynb

```

### siMBD3 rep123 merged peak calling

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq
mkdir macs2_peak_calling_2_for_siMBD3_rep123_together
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siMBD3_rep123_together
# 创建所需的bam文件连接
# bam
ln -s ../HEK293T_siMBD3_input_rep1.bowtie2.sorted.bam HEK293T_siMBD3_input_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD3_input_rep2.bowtie2.sorted.bam HEK293T_siMBD3_input_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD3_input_rep3.bowtie2.sorted.bam HEK293T_siMBD3_input_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD3_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siMBD3_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD3_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siMBD3_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD3_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siMBD3_IP_CTCF_rep3.bowtie2.sorted.bam
# bam.bai
ln -s ../HEK293T_siMBD3_input_rep1.bowtie2.sorted.bam.bai HEK293T_siMBD3_input_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD3_input_rep2.bowtie2.sorted.bam.bai HEK293T_siMBD3_input_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD3_input_rep3.bowtie2.sorted.bam.bai HEK293T_siMBD3_input_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD3_IP_CTCF_rep1.bowtie2.sorted.bam.bai HEK293T_siMBD3_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD3_IP_CTCF_rep2.bowtie2.sorted.bam.bai HEK293T_siMBD3_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD3_IP_CTCF_rep3.bowtie2.sorted.bam.bai HEK293T_siMBD3_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *.bam

# siMBD3_IP_CTCF_rep123_TvsC
nohup macs2 callpeak -t HEK293T_siMBD3_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siMBD3_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siMBD3_IP_CTCF_rep3.bowtie2.sorted.bam -c HEK293T_siMBD3_input_rep1.bowtie2.sorted.bam HEK293T_siMBD3_input_rep2.bowtie2.sorted.bam HEK293T_siMBD3_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siMBD3_IP_CTCF_rep123_TvsC -n siMBD3_IP_CTCF_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_siMBD3_IP_CTCF_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_siMBD3_IP_CTCF_rep123_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siMBD3_rep123_together
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siMBD3_IP_CTCF_rep123_TvsC/siMBD3_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# 26459

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siMBD3_IP_CTCF_rep123_TvsC/siMBD3_IP_CTCF_rep123_TvsC_summits.bed
# 26973

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_siMBD3_IP_CTCF_ChIP_Seq_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siMBD3_rep123_together
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

```
#### Intersection with the siNC

```bash
# 两组结果交集
# cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together
# for i in `find . -name *narrowPeak`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
# wc -l `find . -name *_rep123_TvsC_peaks.narrowPeak.bed`
# # 36855 ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
# bedtools slop -b 50 -i ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
# wc -l `find . -name *summits*.bed`
# # 38169 ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.bed
# # 38169 ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed

cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siMBD3_rep123_together
for i in `find . -name *narrowPeak`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
wc -l `find . -name *_rep123_TvsC_peaks.narrowPeak.bed`
# 26459 ./siMBD3_IP_CTCF_rep123_TvsC/siMBD3_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
bedtools slop -b 50 -i ./siMBD3_IP_CTCF_rep123_TvsC/siMBD3_IP_CTCF_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./siMBD3_IP_CTCF_rep123_TvsC/siMBD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
wc -l `find . -name *summits*.bed`
# 26973 ./siMBD3_IP_CTCF_rep123_TvsC/siMBD3_IP_CTCF_rep123_TvsC_summits.bed
# 26973 ./siMBD3_IP_CTCF_rep123_TvsC/siMBD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed

# 创建 siNC narrowPeak bed 以及 summits extB50 bed 文件连接
ln -s ../macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
ln -s ../macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
# 创建 siMBD3 narrowPeak bed 以及 summits extB50 bed 文件连接
ln -s ./siMBD3_IP_CTCF_rep123_TvsC/siMBD3_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siMBD3_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
ln -s ./siMBD3_IP_CTCF_rep123_TvsC/siMBD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siMBD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed

# 两组 narrowPeak 交集数量
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -b siMBD3_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 22737
bedtools intersect -a siMBD3_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -b siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 22806

# 两组 summit extB50 交集数量
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siMBD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 22850
bedtools intersect -a siMBD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 22848

# 构造 三类 summit
# common summit
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siMBD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siMBD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >IP_CTCF_rep123_TvsC_siNC_and_siMBD3.common.summit.bed
# siMBD3 only summit
bedtools intersect -a siMBD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a siMBD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >IP_CTCF_rep123_TvsC_siMBD3.only.summit.bed
# siNC only summit
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siMBD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b siMBD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >IP_CTCF_rep123_TvsC_siNC.only.summit.bed
wc -l *summit.bed
#  4125 IP_CTCF_rep123_TvsC_siMBD3.only.summit.bed
# 22856 IP_CTCF_rep123_TvsC_siNC_and_siMBD3.common.summit.bed
# 15319 IP_CTCF_rep123_TvsC_siNC.only.summit.bed

# intervene plot
# siNC and siMBD3 peaks venn
intervene venn -i siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siMBD3_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed --names=siNC_CTCF_peaks,siMBD3_CTCF_peaks -o 0_siNC_siMBD3_CTCF_peaks_venn --save-overlaps
# siNC and siMBD3 summit extB50 venn
intervene venn -i siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siMBD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed --names=siNC_CTCF_summits_eB50,siMBD3_CTCF_summits_eB50 -o 0_siNC_siMBD3_CTCF_summits_extB50_venn --save-overlaps

# jupyter make scaled venn plot based on number 
# 0_make_siNC_siMBD3_scaled_venn_plot.based_on_number.ipynb

```

### siNC rep123 merged peak calling

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq
mkdir macs2_peak_calling_2_for_siNC_rep123_together
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together
du -hs
rm *.bam
rm *.bam.bai
du -hs
# 创建所需的bam文件连接
# bam
ln -s ../HEK293T_siNC_input_rep1.bowtie2.sorted.bam HEK293T_siNC_input_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_input_rep2.bowtie2.sorted.bam HEK293T_siNC_input_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_input_rep3.bowtie2.sorted.bam HEK293T_siNC_input_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam
# bam.bai
ln -s ../HEK293T_siNC_input_rep1.bowtie2.sorted.bam.bai HEK293T_siNC_input_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_input_rep2.bowtie2.sorted.bam.bai HEK293T_siNC_input_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_input_rep3.bowtie2.sorted.bam.bai HEK293T_siNC_input_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *.bam

# siNC_IP_CTCF_rep123_TvsC
nohup macs2 callpeak -t HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam -c HEK293T_siNC_input_rep1.bowtie2.sorted.bam HEK293T_siNC_input_rep2.bowtie2.sorted.bam HEK293T_siNC_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir siNC_IP_CTCF_rep123_TvsC -n siNC_IP_CTCF_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_siNC_IP_CTCF_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_siNC_IP_CTCF_rep123_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# 36855

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.bed
# 38169

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_siNC_IP_CTCF_ChIP_Seq_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

```

#### Intersection with the siCHD3/4/5 and siMBD2/3

```bash
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together
# 创建 siNC narrowPeak bed 以及 summits extB50 bed 文件连接
ln -s ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
ln -s ./siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
# 创建 siCHD3 narrowPeak bed 以及 summits extB50 bed 文件连接
ln -s ../macs2_peak_calling_2_for_siCHD3_rep123_together/siCHD3_IP_CTCF_rep123_TvsC/siCHD3_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siCHD3_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
ln -s ../macs2_peak_calling_2_for_siCHD3_rep123_together/siCHD3_IP_CTCF_rep123_TvsC/siCHD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siCHD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
# 创建 siCHD4 narrowPeak bed 以及 summits extB50 bed 文件连接
ln -s ../macs2_peak_calling_2_for_siCHD4_rep123_together/siCHD4_IP_CTCF_rep123_TvsC/siCHD4_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siCHD4_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
ln -s ../macs2_peak_calling_2_for_siCHD4_rep123_together/siCHD4_IP_CTCF_rep123_TvsC/siCHD4_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siCHD4_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
# 创建 siCHD5 narrowPeak bed 以及 summits extB50 bed 文件连接
ln -s ../macs2_peak_calling_2_for_siCHD5_rep123_together/siCHD5_IP_CTCF_rep123_TvsC/siCHD5_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siCHD5_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
ln -s ../macs2_peak_calling_2_for_siCHD5_rep123_together/siCHD5_IP_CTCF_rep123_TvsC/siCHD5_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siCHD5_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
# 创建 siMBD2 narrowPeak bed 以及 summits extB50 bed 文件连接
ln -s ../macs2_peak_calling_2_for_siMBD2_rep123_together/siMBD2_IP_CTCF_rep123_TvsC/siMBD2_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siMBD2_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
ln -s ../macs2_peak_calling_2_for_siMBD2_rep123_together/siMBD2_IP_CTCF_rep123_TvsC/siMBD2_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siMBD2_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
# 创建 siMBD3 narrowPeak bed 以及 summits extB50 bed 文件连接
ln -s ../macs2_peak_calling_2_for_siMBD3_rep123_together/siMBD3_IP_CTCF_rep123_TvsC/siMBD3_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siMBD3_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
ln -s ../macs2_peak_calling_2_for_siMBD3_rep123_together/siMBD3_IP_CTCF_rep123_TvsC/siMBD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siMBD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed

# intervene plot
# all samples peaks venn
intervene venn -i siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siCHD3_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siCHD4_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siCHD5_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siMBD2_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siMBD3_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed --names=siNC_CTCF_peaks,siCHD3_CTCF_peaks,siCHD4_CTCF_peaks,siCHD5_CTCF_peaks,siMBD2_CTCF_peaks,siMBD3_CTCF_peaks -o 0_all_samples_CTCF_peaks_venn --project 0_all_samples_CTCF_peaks --save-overlaps
# all samples summit extB50 venn
intervene venn -i siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siCHD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siCHD4_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siCHD5_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siMBD2_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siMBD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed --names=siNC_CTCF_summits_eB50,siCHD3_CTCF_summits_eB50,siCHD4_CTCF_summits_eB50,siCHD5_CTCF_summits_eB50,siMBD2_CTCF_summits_eB50,siMBD3_CTCF_summits_eB50 -o 0_all_samples_CTCF_summits_extB50_venn --project 0_all_samples_CTCF_summits_extB50 --save-overlaps

# venn 绘图结果没法看，改用upset
# intervene plot
# all samples peaks upset
intervene upset -i siNC_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siCHD3_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siCHD4_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siCHD5_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siMBD2_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed siMBD3_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed --names=siNC_CTCF_peaks,siCHD3_CTCF_peaks,siCHD4_CTCF_peaks,siCHD5_CTCF_peaks,siMBD2_CTCF_peaks,siMBD3_CTCF_peaks -o 0_all_samples_CTCF_peaks_upset --project 0_all_samples_CTCF_peaks --save-overlaps
# all samples summit extB50 upset
intervene upset -i siNC_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siCHD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siCHD4_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siCHD5_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siMBD2_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed siMBD3_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed --names=siNC_CTCF_summits_eB50,siCHD3_CTCF_summits_eB50,siCHD4_CTCF_summits_eB50,siCHD5_CTCF_summits_eB50,siMBD2_CTCF_summits_eB50,siMBD3_CTCF_summits_eB50 -o 0_all_samples_CTCF_summits_extB50_upset --project 0_all_samples_CTCF_summits_extB50 --save-overlaps

```



#  peak summit reads enrichment heatmap

## siNC & siCHD3 IP_VS_input

### Use each individual bam separately

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/
mkdir plotHeatMap_for_siNC_and_siCHD3
cd plotHeatMap_for_siNC_and_siCHD3
# 创建 bam 和 bai 文件链接
## bam
ln -s ../HEK293T_siNC_input_rep1.bowtie2.sorted.bam HEK293T_siNC_input_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_input_rep2.bowtie2.sorted.bam HEK293T_siNC_input_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_input_rep3.bowtie2.sorted.bam HEK293T_siNC_input_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD3_input_rep1.bowtie2.sorted.bam HEK293T_siCHD3_input_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD3_input_rep2.bowtie2.sorted.bam HEK293T_siCHD3_input_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD3_input_rep3.bowtie2.sorted.bam HEK293T_siCHD3_input_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD3_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD3_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD3_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD3_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD3_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siCHD3_IP_CTCF_rep3.bowtie2.sorted.bam
## bam.bai
ln -s ../HEK293T_siNC_input_rep1.bowtie2.sorted.bam.bai HEK293T_siNC_input_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_input_rep2.bowtie2.sorted.bam.bai HEK293T_siNC_input_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_input_rep3.bowtie2.sorted.bam.bai HEK293T_siNC_input_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD3_input_rep1.bowtie2.sorted.bam.bai HEK293T_siCHD3_input_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD3_input_rep2.bowtie2.sorted.bam.bai HEK293T_siCHD3_input_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD3_input_rep3.bowtie2.sorted.bam.bai HEK293T_siCHD3_input_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD3_IP_CTCF_rep1.bowtie2.sorted.bam.bai HEK293T_siCHD3_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD3_IP_CTCF_rep2.bowtie2.sorted.bam.bai HEK293T_siCHD3_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD3_IP_CTCF_rep3.bowtie2.sorted.bam.bai HEK293T_siCHD3_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *bowtie2.sorted.bam

# 创建 siNC CTCF peak summit 文件连接
ln -s ../macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_siCHD3_rep123_together/IP_CTCF_rep123_TvsC_siNC.only.summit.bed siNC_only.CTCF.summit.c3su.bed
ln -s ../macs2_peak_calling_2_for_siCHD3_rep123_together/IP_CTCF_rep123_TvsC_siNC_and_siCHD3.common.summit.bed siNC_siCHD3_common.CTCF.summit.c3su.bed
ln -s ../macs2_peak_calling_2_for_siCHD3_rep123_together/IP_CTCF_rep123_TvsC_siCHD3.only.summit.bed siCHD3_only.CTCF.summit.c3su.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_siNC_vs_siCHD3_bdgdiff/ROI_siNC.CTCF.summit.sorted.bed ROI_siNC.CTCF.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_siNC_vs_siCHD3_bdgdiff/RIB_siNC_and_siCHD3.CTCF.summit.sorted.bed RIB_siNC_and_siCHD3.CTCF.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_siNC_vs_siCHD3_bdgdiff/ROI_siCHD3.CTCF.summit.sorted.bed ROI_siCHD3.CTCF.summit.sorted.bed
wc -l *.bed
#     0 ROI_siNC.CTCF.summit.sorted.bed
#  1832 RIB_siNC_and_siCHD3.CTCF.summit.sorted.bed
#  7928 ROI_siCHD3.CTCF.summit.sorted.bed
# 14062 siNC_only.CTCF.summit.c3su.bed
# 24110 siNC_siCHD3_common.CTCF.summit.c3su.bed
# 13684 siCHD3_only.CTCF.summit.c3su.bed
# 38169 siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed

```

#### Create log2 transformed bigwig files

```bash
# 构造 7_run_bamCompare_for_siNC_siCHD3_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# bamCompare -b1 HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam -b2 HEK293T_siNC_input_rep1.bowtie2.sorted.bam -o siNC_IP_CTCF_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam -b2 HEK293T_siNC_input_rep2.bowtie2.sorted.bam -o siNC_IP_CTCF_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam -b2 HEK293T_siNC_input_rep3.bowtie2.sorted.bam -o siNC_IP_CTCF_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siCHD3_IP_CTCF_rep1.bowtie2.sorted.bam -b2 HEK293T_siCHD3_input_rep1.bowtie2.sorted.bam -o siCHD3_IP_CTCF_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siCHD3_IP_CTCF_rep2.bowtie2.sorted.bam -b2 HEK293T_siCHD3_input_rep2.bowtie2.sorted.bam -o siCHD3_IP_CTCF_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siCHD3_IP_CTCF_rep3.bowtie2.sorted.bam -b2 HEK293T_siCHD3_input_rep3.bowtie2.sorted.bam -o siCHD3_IP_CTCF_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
nohup bash 7_run_bamCompare_for_siNC_siCHD3_each_repeat.sh >7_run_bamCompare_for_siNC_siCHD3_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bamCompare_for_siNC_siCHD3_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD3
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siCHD3_IP_CTCF_rep1_TvsC.bigwig siCHD3_IP_CTCF_rep2_TvsC.bigwig siCHD3_IP_CTCF_rep3_TvsC.bigwig -R siNC_only.CTCF.summit.c3su.bed siNC_siCHD3_common.CTCF.summit.c3su.bed siCHD3_only.CTCF.summit.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD3_all_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD3_all_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_siNC_siCHD3_all_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_only(14062)" "common_summits(24110)" "siCHD3_only(13684)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siCHD3_CTCF_rep1_TvsC" "siCHD3_CTCF_rep2_TvsC" "siCHD3_CTCF_rep3_TvsC"

# use summits generated by bdgdiff
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD3
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siCHD3_IP_CTCF_rep1_TvsC.bigwig siCHD3_IP_CTCF_rep2_TvsC.bigwig siCHD3_IP_CTCF_rep3_TvsC.bigwig -R RIB_siNC_and_siCHD3.CTCF.summit.sorted.bed ROI_siCHD3.CTCF.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD3_all_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD3_all_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_siNC_siCHD3_all_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "RIB_siNC_and_siCHD3(1832)" "ROI_siCHD3(7928)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siCHD3_CTCF_rep1_TvsC" "siCHD3_CTCF_rep2_TvsC" "siCHD3_CTCF_rep3_TvsC"

# use siNC CTCF summits
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD3
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siCHD3_IP_CTCF_rep1_TvsC.bigwig siCHD3_IP_CTCF_rep2_TvsC.bigwig siCHD3_IP_CTCF_rep3_TvsC.bigwig -R siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD3_all_samples_for_plotHeatMap.siNC_CTCF_summit.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD3_all_samples_for_plotHeatMap.siNC_CTCF_summit.gz -o bamCompare_matrix_of_siNC_siCHD3_all_samples_for_plotHeatMap.siNC_CTCF_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_CTCF(38169)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siCHD3_CTCF_rep1_TvsC" "siCHD3_CTCF_rep2_TvsC" "siCHD3_CTCF_rep3_TvsC"

# use 293T 2w CTCF summit
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD3
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siCHD3_IP_CTCF_rep1_TvsC.bigwig siCHD3_IP_CTCF_rep2_TvsC.bigwig siCHD3_IP_CTCF_rep3_TvsC.bigwig -R ../CTCF.summit.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD3_all_samples_for_plotHeatMap.293T2wCTCF.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD3_all_samples_for_plotHeatMap.293T2wCTCF.gz -o bamCompare_matrix_of_siNC_siCHD3_all_samples_for_plotHeatMap.293T2wCTCF.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_summits(20000)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siCHD3_CTCF_rep1_TvsC" "siCHD3_CTCF_rep2_TvsC" "siCHD3_CTCF_rep3_TvsC"


# use siNC TET1 summits (28450)
# wc -l /home/sunwenju/projects_on_940xa/2023122801_WN_siNuRD_TET1_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_TET1_rep123_TvsC/siNC_IP_TET1_rep123_TvsC_summits.c3su.bed
# 28450 /home/sunwenju/projects_on_940xa/2023122801_WN_siNuRD_TET1_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_TET1_rep123_TvsC/siNC_IP_TET1_rep123_TvsC_summits.c3su.bed
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD3
# 创建 siNC TET1 peak summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2023122801_WN_siNuRD_TET1_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_TET1_rep123_TvsC/siNC_IP_TET1_rep123_TvsC_summits.c3su.bed siNC_IP_TET1_rep123_TvsC_summits.c3su.bed
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siCHD3_IP_CTCF_rep1_TvsC.bigwig siCHD3_IP_CTCF_rep2_TvsC.bigwig siCHD3_IP_CTCF_rep3_TvsC.bigwig -R siNC_IP_TET1_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD3_all_samples_for_plotHeatMap.siNC_TET1_summit.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD3_all_samples_for_plotHeatMap.siNC_TET1_summit.gz -o bamCompare_matrix_of_siNC_siCHD3_all_samples_for_plotHeatMap.siNC_TET1_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_TET1_summits(28450)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siCHD3_CTCF_rep1_TvsC" "siCHD3_CTCF_rep2_TvsC" "siCHD3_CTCF_rep3_TvsC"

```

### Use rep123 merged bam

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD3/
# 创建 rep123 merged bam 和 bai 文件链接
## bam
ln -s ../HEK293T_siNC_input_rep123.merged.sorted.bam HEK293T_siNC_input_rep123.merged.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam
ln -s ../HEK293T_siCHD3_input_rep123.merged.sorted.bam HEK293T_siCHD3_input_rep123.merged.sorted.bam
ln -s ../HEK293T_siCHD3_IP_CTCF_rep123.merged.sorted.bam HEK293T_siCHD3_IP_CTCF_rep123.merged.sorted.bam
## bam.bai
ln -s ../HEK293T_siNC_input_rep123.merged.sorted.bam.bai HEK293T_siNC_input_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_siCHD3_input_rep123.merged.sorted.bam.bai HEK293T_siCHD3_input_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_siCHD3_IP_CTCF_rep123.merged.sorted.bam.bai HEK293T_siCHD3_IP_CTCF_rep123.merged.sorted.bam.bai
ls -1 *rep123.merged.sorted.bam

```

#### Create log2 transformed bigwig files

```bash
# 构造 7_run_bamCompare_for_siNC_siCHD3_rep123_merged_bam.sh 内容如下
# # Compare for rep123 merged bam
# bamCompare -b1 HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam -b2 HEK293T_siNC_input_rep123.merged.sorted.bam -o siNC_IP_CTCF_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siCHD3_IP_CTCF_rep123.merged.sorted.bam -b2 HEK293T_siCHD3_input_rep123.merged.sorted.bam -o siCHD3_IP_CTCF_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
nohup bash 7_run_bamCompare_for_siNC_siCHD3_rep123_merged_bam.sh >7_run_bamCompare_for_siNC_siCHD3_rep123_merged_bam.sh.runlog 2>&1 &
tail -f 7_run_bamCompare_for_siNC_siCHD3_rep123_merged_bam.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD3
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siCHD3_IP_CTCF_rep123_merged_TvsC.bigwig -R siNC_only.CTCF.summit.c3su.bed siNC_siCHD3_common.CTCF.summit.c3su.bed siCHD3_only.CTCF.summit.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD3_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD3_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_siNC_siCHD3_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_only(14062)" "common_summits(24110)" "siCHD3_only(13684)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siCHD3_CTCF_rep123_TvsC"

# use summits generated by bdgdiff
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD3
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siCHD3_IP_CTCF_rep123_merged_TvsC.bigwig -R RIB_siNC_and_siCHD3.CTCF.summit.sorted.bed ROI_siCHD3.CTCF.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD3_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD3_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_siNC_siCHD3_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "RIB_siNC_and_siCHD3(1832)" "ROI_siCHD3(7928)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siCHD3_CTCF_rep123_TvsC"

# use siNC CTCF summits
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD3
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siCHD3_IP_CTCF_rep123_merged_TvsC.bigwig -R siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD3_rep123_merged_for_plotHeatMap.siNC_CTCF_summit.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD3_rep123_merged_for_plotHeatMap.siNC_CTCF_summit.gz -o bamCompare_matrix_of_siNC_siCHD3_rep123_merged_for_plotHeatMap.siNC_CTCF_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_CTCF(38169)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siCHD3_CTCF_rep123_TvsC"

# use 293T 2w CTCF summit
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD3
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siCHD3_IP_CTCF_rep123_merged_TvsC.bigwig -R ../CTCF.summit.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD3_rep123_merged_for_plotHeatMap.293T2wCTCF.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD3_rep123_merged_for_plotHeatMap.293T2wCTCF.gz -o bamCompare_matrix_of_siNC_siCHD3_rep123_merged_for_plotHeatMap.293T2wCTCF.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_summits(20000)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siCHD3_CTCF_rep123_TvsC"

# use siNC TET1 summits (28450)
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD3
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siCHD3_IP_CTCF_rep123_merged_TvsC.bigwig -R siNC_IP_TET1_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD3_rep123_merged_for_plotHeatMap.siNC_TET1_summit.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD3_rep123_merged_for_plotHeatMap.siNC_TET1_summit.gz -o bamCompare_matrix_of_siNC_siCHD3_rep123_merged_for_plotHeatMap.siNC_TET1_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_TET1_summits(28450)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siCHD3_CTCF_rep123_TvsC"

```



## siNC & siCHD4 IP_VS_input
### Use each individual bam separately

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/
mkdir plotHeatMap_for_siNC_and_siCHD4
cd plotHeatMap_for_siNC_and_siCHD4
# 创建 bam 和 bai 文件链接
## bam
ln -s ../HEK293T_siNC_input_rep1.bowtie2.sorted.bam HEK293T_siNC_input_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_input_rep2.bowtie2.sorted.bam HEK293T_siNC_input_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_input_rep3.bowtie2.sorted.bam HEK293T_siNC_input_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD4_input_rep1.bowtie2.sorted.bam HEK293T_siCHD4_input_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD4_input_rep2.bowtie2.sorted.bam HEK293T_siCHD4_input_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD4_input_rep3.bowtie2.sorted.bam HEK293T_siCHD4_input_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD4_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD4_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD4_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD4_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD4_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siCHD4_IP_CTCF_rep3.bowtie2.sorted.bam
## bam.bai
ln -s ../HEK293T_siNC_input_rep1.bowtie2.sorted.bam.bai HEK293T_siNC_input_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_input_rep2.bowtie2.sorted.bam.bai HEK293T_siNC_input_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_input_rep3.bowtie2.sorted.bam.bai HEK293T_siNC_input_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD4_input_rep1.bowtie2.sorted.bam.bai HEK293T_siCHD4_input_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD4_input_rep2.bowtie2.sorted.bam.bai HEK293T_siCHD4_input_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD4_input_rep3.bowtie2.sorted.bam.bai HEK293T_siCHD4_input_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD4_IP_CTCF_rep1.bowtie2.sorted.bam.bai HEK293T_siCHD4_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD4_IP_CTCF_rep2.bowtie2.sorted.bam.bai HEK293T_siCHD4_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD4_IP_CTCF_rep3.bowtie2.sorted.bam.bai HEK293T_siCHD4_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *bowtie2.sorted.bam

# 创建 siNC CTCF peak summit 文件连接
ln -s ../macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_siCHD4_rep123_together/IP_CTCF_rep123_TvsC_siNC.only.summit.bed siNC_only.CTCF.summit.c3su.bed
ln -s ../macs2_peak_calling_2_for_siCHD4_rep123_together/IP_CTCF_rep123_TvsC_siNC_and_siCHD4.common.summit.bed siNC_siCHD4_common.CTCF.summit.c3su.bed
ln -s ../macs2_peak_calling_2_for_siCHD4_rep123_together/IP_CTCF_rep123_TvsC_siCHD4.only.summit.bed siCHD4_only.CTCF.summit.c3su.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_siNC_vs_siCHD4_bdgdiff/ROI_siNC.CTCF.summit.sorted.bed ROI_siNC.CTCF.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_siNC_vs_siCHD4_bdgdiff/RIB_siNC_and_siCHD4.CTCF.summit.sorted.bed RIB_siNC_and_siCHD4.CTCF.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_siNC_vs_siCHD4_bdgdiff/ROI_siCHD4.CTCF.summit.sorted.bed ROI_siCHD4.CTCF.summit.sorted.bed
wc -l *.bed
#  8318 RIB_siNC_and_siCHD4.CTCF.summit.sorted.bed
#  2171 ROI_siCHD4.CTCF.summit.sorted.bed
#    26 ROI_siNC.CTCF.summit.sorted.bed
# 16185 siNC_only.CTCF.summit.c3su.bed
# 21998 siNC_siCHD4_common.CTCF.summit.c3su.bed
#  5560 siCHD4_only.CTCF.summit.c3su.bed
# 38169 siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed

```

#### Create log2 transformed bigwig files

```bash
# 构造 7_run_bamCompare_for_siNC_siCHD4_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# bamCompare -b1 HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam -b2 HEK293T_siNC_input_rep1.bowtie2.sorted.bam -o siNC_IP_CTCF_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam -b2 HEK293T_siNC_input_rep2.bowtie2.sorted.bam -o siNC_IP_CTCF_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam -b2 HEK293T_siNC_input_rep3.bowtie2.sorted.bam -o siNC_IP_CTCF_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siCHD4_IP_CTCF_rep1.bowtie2.sorted.bam -b2 HEK293T_siCHD4_input_rep1.bowtie2.sorted.bam -o siCHD4_IP_CTCF_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siCHD4_IP_CTCF_rep2.bowtie2.sorted.bam -b2 HEK293T_siCHD4_input_rep2.bowtie2.sorted.bam -o siCHD4_IP_CTCF_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siCHD4_IP_CTCF_rep3.bowtie2.sorted.bam -b2 HEK293T_siCHD4_input_rep3.bowtie2.sorted.bam -o siCHD4_IP_CTCF_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
nohup bash 7_run_bamCompare_for_siNC_siCHD4_each_repeat.sh >7_run_bamCompare_for_siNC_siCHD4_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bamCompare_for_siNC_siCHD4_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD4
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siCHD4_IP_CTCF_rep1_TvsC.bigwig siCHD4_IP_CTCF_rep2_TvsC.bigwig siCHD4_IP_CTCF_rep3_TvsC.bigwig -R siNC_only.CTCF.summit.c3su.bed siNC_siCHD4_common.CTCF.summit.c3su.bed siCHD4_only.CTCF.summit.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD4_all_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD4_all_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_siNC_siCHD4_all_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_only(16185)" "common_summits(21998)" "siCHD4_only(5560)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siCHD4_CTCF_rep1_TvsC" "siCHD4_CTCF_rep2_TvsC" "siCHD4_CTCF_rep3_TvsC"

# use summits generated by bdgdiff
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD4
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siCHD4_IP_CTCF_rep1_TvsC.bigwig siCHD4_IP_CTCF_rep2_TvsC.bigwig siCHD4_IP_CTCF_rep3_TvsC.bigwig -R ROI_siNC.CTCF.summit.sorted.bed RIB_siNC_and_siCHD4.CTCF.summit.sorted.bed ROI_siCHD4.CTCF.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD4_all_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD4_all_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_siNC_siCHD4_all_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "ROI_siNC(26)" "RIB_siNC_and_siCHD4(8318)" "ROI_siCHD4(2171)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siCHD4_CTCF_rep1_TvsC" "siCHD4_CTCF_rep2_TvsC" "siCHD4_CTCF_rep3_TvsC"

# use siNC CTCF summits
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD4
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siCHD4_IP_CTCF_rep1_TvsC.bigwig siCHD4_IP_CTCF_rep2_TvsC.bigwig siCHD4_IP_CTCF_rep3_TvsC.bigwig -R siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD4_all_samples_for_plotHeatMap.siNC_CTCF_summit.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD4_all_samples_for_plotHeatMap.siNC_CTCF_summit.gz -o bamCompare_matrix_of_siNC_siCHD4_all_samples_for_plotHeatMap.siNC_CTCF_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_CTCF(38169)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siCHD4_CTCF_rep1_TvsC" "siCHD4_CTCF_rep2_TvsC" "siCHD4_CTCF_rep3_TvsC"

# use 293T 2w CTCF summit
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD4
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siCHD4_IP_CTCF_rep1_TvsC.bigwig siCHD4_IP_CTCF_rep2_TvsC.bigwig siCHD4_IP_CTCF_rep3_TvsC.bigwig -R ../CTCF.summit.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD4_all_samples_for_plotHeatMap.293T2wCTCF.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD4_all_samples_for_plotHeatMap.293T2wCTCF.gz -o bamCompare_matrix_of_siNC_siCHD4_all_samples_for_plotHeatMap.293T2wCTCF.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_summits(20000)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siCHD4_CTCF_rep1_TvsC" "siCHD4_CTCF_rep2_TvsC" "siCHD4_CTCF_rep3_TvsC"

# use siNC TET1 summits (28450)
# wc -l /home/sunwenju/projects_on_940xa/2023122801_WN_siNuRD_TET1_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_TET1_rep123_TvsC/siNC_IP_TET1_rep123_TvsC_summits.c3su.bed
# 28450 /home/sunwenju/projects_on_940xa/2023122801_WN_siNuRD_TET1_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_TET1_rep123_TvsC/siNC_IP_TET1_rep123_TvsC_summits.c3su.bed
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD4
# 创建 siNC TET1 peak summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2023122801_WN_siNuRD_TET1_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_TET1_rep123_TvsC/siNC_IP_TET1_rep123_TvsC_summits.c3su.bed siNC_IP_TET1_rep123_TvsC_summits.c3su.bed
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siCHD4_IP_CTCF_rep1_TvsC.bigwig siCHD4_IP_CTCF_rep2_TvsC.bigwig siCHD4_IP_CTCF_rep3_TvsC.bigwig -R siNC_IP_TET1_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD4_all_samples_for_plotHeatMap.siNC_TET1_summit.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD4_all_samples_for_plotHeatMap.siNC_TET1_summit.gz -o bamCompare_matrix_of_siNC_siCHD4_all_samples_for_plotHeatMap.siNC_TET1_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_TET1_summits(28450)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siCHD4_CTCF_rep1_TvsC" "siCHD4_CTCF_rep2_TvsC" "siCHD4_CTCF_rep3_TvsC"

```

### Use rep123 merged bam

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD4/
# 创建 rep123 merged bam 和 bai 文件链接
## bam
ln -s ../HEK293T_siNC_input_rep123.merged.sorted.bam HEK293T_siNC_input_rep123.merged.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam
ln -s ../HEK293T_siCHD4_input_rep123.merged.sorted.bam HEK293T_siCHD4_input_rep123.merged.sorted.bam
ln -s ../HEK293T_siCHD4_IP_CTCF_rep123.merged.sorted.bam HEK293T_siCHD4_IP_CTCF_rep123.merged.sorted.bam
## bam.bai
ln -s ../HEK293T_siNC_input_rep123.merged.sorted.bam.bai HEK293T_siNC_input_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_siCHD4_input_rep123.merged.sorted.bam.bai HEK293T_siCHD4_input_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_siCHD4_IP_CTCF_rep123.merged.sorted.bam.bai HEK293T_siCHD4_IP_CTCF_rep123.merged.sorted.bam.bai
ls -1 *rep123.merged.sorted.bam

```

#### Create log2 transformed bigwig files

```bash
# 构造 7_run_bamCompare_for_siNC_siCHD4_rep123_merged_bam.sh 内容如下
# # Compare for rep123 merged bam
# bamCompare -b1 HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam -b2 HEK293T_siNC_input_rep123.merged.sorted.bam -o siNC_IP_CTCF_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siCHD4_IP_CTCF_rep123.merged.sorted.bam -b2 HEK293T_siCHD4_input_rep123.merged.sorted.bam -o siCHD4_IP_CTCF_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
nohup bash 7_run_bamCompare_for_siNC_siCHD4_rep123_merged_bam.sh >7_run_bamCompare_for_siNC_siCHD4_rep123_merged_bam.sh.runlog 2>&1 &
tail -f 7_run_bamCompare_for_siNC_siCHD4_rep123_merged_bam.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD4
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siCHD4_IP_CTCF_rep123_merged_TvsC.bigwig -R siNC_only.CTCF.summit.c3su.bed siNC_siCHD4_common.CTCF.summit.c3su.bed siCHD4_only.CTCF.summit.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD4_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD4_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_siNC_siCHD4_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_only(16185)" "common_summits(21998)" "siCHD4_only(5560)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siCHD4_CTCF_rep123_TvsC"

# use summits generated by bdgdiff
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD4
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siCHD4_IP_CTCF_rep123_merged_TvsC.bigwig -R ROI_siNC.CTCF.summit.sorted.bed RIB_siNC_and_siCHD4.CTCF.summit.sorted.bed ROI_siCHD4.CTCF.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD4_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD4_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_siNC_siCHD4_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "ROI_siNC(26)" "RIB_siNC_and_siCHD4(8318)" "ROI_siCHD4(2171)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siCHD4_CTCF_rep123_TvsC"

# use siNC CTCF summits
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD4
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siCHD4_IP_CTCF_rep123_merged_TvsC.bigwig -R siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD4_rep123_merged_for_plotHeatMap.siNC_CTCF_summit.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD4_rep123_merged_for_plotHeatMap.siNC_CTCF_summit.gz -o bamCompare_matrix_of_siNC_siCHD4_rep123_merged_for_plotHeatMap.siNC_CTCF_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_CTCF(38169)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siCHD4_CTCF_rep123_TvsC"

# use 293T 2w CTCF summit
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD4
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siCHD4_IP_CTCF_rep123_merged_TvsC.bigwig -R ../CTCF.summit.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD4_rep123_merged_for_plotHeatMap.293T2wCTCF.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD4_rep123_merged_for_plotHeatMap.293T2wCTCF.gz -o bamCompare_matrix_of_siNC_siCHD4_rep123_merged_for_plotHeatMap.293T2wCTCF.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_summits(20000)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siCHD4_CTCF_rep123_TvsC"

# use siNC TET1 summits (28450)
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD4
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siCHD4_IP_CTCF_rep123_merged_TvsC.bigwig -R siNC_IP_TET1_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD4_rep123_merged_for_plotHeatMap.siNC_TET1_summit.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD4_rep123_merged_for_plotHeatMap.siNC_TET1_summit.gz -o bamCompare_matrix_of_siNC_siCHD4_rep123_merged_for_plotHeatMap.siNC_TET1_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_TET1_summits(28450)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siCHD4_CTCF_rep123_TvsC"


```




## siNC & siCHD5 IP_VS_input
### Use each individual bam separately

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/
mkdir plotHeatMap_for_siNC_and_siCHD5
cd plotHeatMap_for_siNC_and_siCHD5
# 创建 bam 和 bai 文件链接
## bam
ln -s ../HEK293T_siNC_input_rep1.bowtie2.sorted.bam HEK293T_siNC_input_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_input_rep2.bowtie2.sorted.bam HEK293T_siNC_input_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_input_rep3.bowtie2.sorted.bam HEK293T_siNC_input_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD5_input_rep1.bowtie2.sorted.bam HEK293T_siCHD5_input_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD5_input_rep2.bowtie2.sorted.bam HEK293T_siCHD5_input_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD5_input_rep3.bowtie2.sorted.bam HEK293T_siCHD5_input_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD5_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siCHD5_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD5_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siCHD5_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siCHD5_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siCHD5_IP_CTCF_rep3.bowtie2.sorted.bam
## bam.bai
ln -s ../HEK293T_siNC_input_rep1.bowtie2.sorted.bam.bai HEK293T_siNC_input_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_input_rep2.bowtie2.sorted.bam.bai HEK293T_siNC_input_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_input_rep3.bowtie2.sorted.bam.bai HEK293T_siNC_input_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD5_input_rep1.bowtie2.sorted.bam.bai HEK293T_siCHD5_input_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD5_input_rep2.bowtie2.sorted.bam.bai HEK293T_siCHD5_input_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD5_input_rep3.bowtie2.sorted.bam.bai HEK293T_siCHD5_input_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD5_IP_CTCF_rep1.bowtie2.sorted.bam.bai HEK293T_siCHD5_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD5_IP_CTCF_rep2.bowtie2.sorted.bam.bai HEK293T_siCHD5_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siCHD5_IP_CTCF_rep3.bowtie2.sorted.bam.bai HEK293T_siCHD5_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *bowtie2.sorted.bam

# 创建 siNC CTCF peak summit 文件连接
ln -s ../macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_siCHD5_rep123_together/IP_CTCF_rep123_TvsC_siNC.only.summit.bed siNC_only.CTCF.summit.c3su.bed
ln -s ../macs2_peak_calling_2_for_siCHD5_rep123_together/IP_CTCF_rep123_TvsC_siNC_and_siCHD5.common.summit.bed siNC_siCHD5_common.CTCF.summit.c3su.bed
ln -s ../macs2_peak_calling_2_for_siCHD5_rep123_together/IP_CTCF_rep123_TvsC_siCHD5.only.summit.bed siCHD5_only.CTCF.summit.c3su.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_siNC_vs_siCHD5_bdgdiff/ROI_siNC.CTCF.summit.sorted.bed ROI_siNC.CTCF.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_siNC_vs_siCHD5_bdgdiff/RIB_siNC_and_siCHD5.CTCF.summit.sorted.bed RIB_siNC_and_siCHD5.CTCF.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_siNC_vs_siCHD5_bdgdiff/ROI_siCHD5.CTCF.summit.sorted.bed ROI_siCHD5.CTCF.summit.sorted.bed
wc -l *.bed
# 20744 RIB_siNC_and_siCHD5.CTCF.summit.sorted.bed
#   257 ROI_siCHD5.CTCF.summit.sorted.bed
#  1027 ROI_siNC.CTCF.summit.sorted.bed
#  8946 siNC_only.CTCF.summit.c3su.bed
# 29224 siNC_siCHD5_common.CTCF.summit.c3su.bed
#  3407 siCHD5_only.CTCF.summit.c3su.bed
# 38169 siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed

```

#### Create log2 transformed bigwig files

```bash
# 构造 7_run_bamCompare_for_siNC_siCHD5_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# bamCompare -b1 HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam -b2 HEK293T_siNC_input_rep1.bowtie2.sorted.bam -o siNC_IP_CTCF_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam -b2 HEK293T_siNC_input_rep2.bowtie2.sorted.bam -o siNC_IP_CTCF_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam -b2 HEK293T_siNC_input_rep3.bowtie2.sorted.bam -o siNC_IP_CTCF_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siCHD5_IP_CTCF_rep1.bowtie2.sorted.bam -b2 HEK293T_siCHD5_input_rep1.bowtie2.sorted.bam -o siCHD5_IP_CTCF_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siCHD5_IP_CTCF_rep2.bowtie2.sorted.bam -b2 HEK293T_siCHD5_input_rep2.bowtie2.sorted.bam -o siCHD5_IP_CTCF_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siCHD5_IP_CTCF_rep3.bowtie2.sorted.bam -b2 HEK293T_siCHD5_input_rep3.bowtie2.sorted.bam -o siCHD5_IP_CTCF_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
nohup bash 7_run_bamCompare_for_siNC_siCHD5_each_repeat.sh >7_run_bamCompare_for_siNC_siCHD5_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bamCompare_for_siNC_siCHD5_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD5
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siCHD5_IP_CTCF_rep1_TvsC.bigwig siCHD5_IP_CTCF_rep2_TvsC.bigwig siCHD5_IP_CTCF_rep3_TvsC.bigwig -R siNC_only.CTCF.summit.c3su.bed siNC_siCHD5_common.CTCF.summit.c3su.bed siCHD5_only.CTCF.summit.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD5_all_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD5_all_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_siNC_siCHD5_all_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_only(8946)" "common_summits(29224)" "siCHD5_only(257)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siCHD5_CTCF_rep1_TvsC" "siCHD5_CTCF_rep2_TvsC" "siCHD5_CTCF_rep3_TvsC"

# use summits generated by bdgdiff
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD5
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siCHD5_IP_CTCF_rep1_TvsC.bigwig siCHD5_IP_CTCF_rep2_TvsC.bigwig siCHD5_IP_CTCF_rep3_TvsC.bigwig -R ROI_siNC.CTCF.summit.sorted.bed RIB_siNC_and_siCHD5.CTCF.summit.sorted.bed ROI_siCHD5.CTCF.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD5_all_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD5_all_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_siNC_siCHD5_all_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "ROI_siNC(1027)" "RIB_siNC_and_siCHD5(20744)" "ROI_siCHD5(2171)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siCHD5_CTCF_rep1_TvsC" "siCHD5_CTCF_rep2_TvsC" "siCHD5_CTCF_rep3_TvsC"

# use siNC CTCF summits
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD5
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siCHD5_IP_CTCF_rep1_TvsC.bigwig siCHD5_IP_CTCF_rep2_TvsC.bigwig siCHD5_IP_CTCF_rep3_TvsC.bigwig -R siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD5_all_samples_for_plotHeatMap.siNC_CTCF_summit.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD5_all_samples_for_plotHeatMap.siNC_CTCF_summit.gz -o bamCompare_matrix_of_siNC_siCHD5_all_samples_for_plotHeatMap.siNC_CTCF_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_CTCF(38169)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siCHD5_CTCF_rep1_TvsC" "siCHD5_CTCF_rep2_TvsC" "siCHD5_CTCF_rep3_TvsC"

# use 293T 2w CTCF summit
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD5
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siCHD5_IP_CTCF_rep1_TvsC.bigwig siCHD5_IP_CTCF_rep2_TvsC.bigwig siCHD5_IP_CTCF_rep3_TvsC.bigwig -R ../CTCF.summit.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD5_all_samples_for_plotHeatMap.293T2wCTCF.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD5_all_samples_for_plotHeatMap.293T2wCTCF.gz -o bamCompare_matrix_of_siNC_siCHD5_all_samples_for_plotHeatMap.293T2wCTCF.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_summits(20000)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siCHD5_CTCF_rep1_TvsC" "siCHD5_CTCF_rep2_TvsC" "siCHD5_CTCF_rep3_TvsC"

# use siNC TET1 summits (28450)
# wc -l /home/sunwenju/projects_on_940xa/2023122801_WN_siNuRD_TET1_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_TET1_rep123_TvsC/siNC_IP_TET1_rep123_TvsC_summits.c3su.bed
# 28450 /home/sunwenju/projects_on_940xa/2023122801_WN_siNuRD_TET1_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_TET1_rep123_TvsC/siNC_IP_TET1_rep123_TvsC_summits.c3su.bed
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD5
# 创建 siNC TET1 peak summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2023122801_WN_siNuRD_TET1_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_TET1_rep123_TvsC/siNC_IP_TET1_rep123_TvsC_summits.c3su.bed siNC_IP_TET1_rep123_TvsC_summits.c3su.bed
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siCHD5_IP_CTCF_rep1_TvsC.bigwig siCHD5_IP_CTCF_rep2_TvsC.bigwig siCHD5_IP_CTCF_rep3_TvsC.bigwig -R siNC_IP_TET1_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD5_all_samples_for_plotHeatMap.siNC_TET1_summit.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD5_all_samples_for_plotHeatMap.siNC_TET1_summit.gz -o bamCompare_matrix_of_siNC_siCHD5_all_samples_for_plotHeatMap.siNC_TET1_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_TET1_summits(28450)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siCHD5_CTCF_rep1_TvsC" "siCHD5_CTCF_rep2_TvsC" "siCHD5_CTCF_rep3_TvsC"


```

### Use rep123 merged bam

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD5/
# 创建 rep123 merged bam 和 bai 文件链接
## bam
ln -s ../HEK293T_siNC_input_rep123.merged.sorted.bam HEK293T_siNC_input_rep123.merged.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam
ln -s ../HEK293T_siCHD5_input_rep123.merged.sorted.bam HEK293T_siCHD5_input_rep123.merged.sorted.bam
ln -s ../HEK293T_siCHD5_IP_CTCF_rep123.merged.sorted.bam HEK293T_siCHD5_IP_CTCF_rep123.merged.sorted.bam
## bam.bai
ln -s ../HEK293T_siNC_input_rep123.merged.sorted.bam.bai HEK293T_siNC_input_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_siCHD5_input_rep123.merged.sorted.bam.bai HEK293T_siCHD5_input_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_siCHD5_IP_CTCF_rep123.merged.sorted.bam.bai HEK293T_siCHD5_IP_CTCF_rep123.merged.sorted.bam.bai
ls -1 *rep123.merged.sorted.bam

```

#### Create log2 transformed bigwig files

```bash
# 构造 7_run_bamCompare_for_siNC_siCHD5_rep123_merged_bam.sh 内容如下
# # Compare for rep123 merged bam
# bamCompare -b1 HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam -b2 HEK293T_siNC_input_rep123.merged.sorted.bam -o siNC_IP_CTCF_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siCHD5_IP_CTCF_rep123.merged.sorted.bam -b2 HEK293T_siCHD5_input_rep123.merged.sorted.bam -o siCHD5_IP_CTCF_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
nohup bash 7_run_bamCompare_for_siNC_siCHD5_rep123_merged_bam.sh >7_run_bamCompare_for_siNC_siCHD5_rep123_merged_bam.sh.runlog 2>&1 &
tail -f 7_run_bamCompare_for_siNC_siCHD5_rep123_merged_bam.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD5
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siCHD5_IP_CTCF_rep123_merged_TvsC.bigwig -R siNC_only.CTCF.summit.c3su.bed siNC_siCHD5_common.CTCF.summit.c3su.bed siCHD5_only.CTCF.summit.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD5_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD5_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_siNC_siCHD5_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_only(8946)" "common_summits(29224)" "siCHD5_only(257)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siCHD5_CTCF_rep123_TvsC"

# use summits generated by bdgdiff
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD5
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siCHD5_IP_CTCF_rep123_merged_TvsC.bigwig -R ROI_siNC.CTCF.summit.sorted.bed RIB_siNC_and_siCHD5.CTCF.summit.sorted.bed ROI_siCHD5.CTCF.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD5_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD5_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_siNC_siCHD5_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "ROI_siNC(1027)" "RIB_siNC_and_siCHD5(20744)" "ROI_siCHD5(2171)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siCHD5_CTCF_rep123_TvsC"

# use siNC CTCF summits
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD5
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siCHD5_IP_CTCF_rep123_merged_TvsC.bigwig -R siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD5_rep123_merged_for_plotHeatMap.siNC_CTCF_summit.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD5_rep123_merged_for_plotHeatMap.siNC_CTCF_summit.gz -o bamCompare_matrix_of_siNC_siCHD5_rep123_merged_for_plotHeatMap.siNC_CTCF_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_CTCF(38169)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siCHD5_CTCF_rep123_TvsC"

# use 293T 2w CTCF summit
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD5
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siCHD5_IP_CTCF_rep123_merged_TvsC.bigwig -R ../CTCF.summit.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD5_rep123_merged_for_plotHeatMap.293T2wCTCF.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD5_rep123_merged_for_plotHeatMap.293T2wCTCF.gz -o bamCompare_matrix_of_siNC_siCHD5_rep123_merged_for_plotHeatMap.293T2wCTCF.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_summits(20000)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siCHD5_CTCF_rep123_TvsC"

# use siNC TET1 summits (28450)
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siCHD5
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siCHD5_IP_CTCF_rep123_merged_TvsC.bigwig -R siNC_IP_TET1_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siCHD5_rep123_merged_for_plotHeatMap.siNC_TET1_summit.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siCHD5_rep123_merged_for_plotHeatMap.siNC_TET1_summit.gz -o bamCompare_matrix_of_siNC_siCHD5_rep123_merged_for_plotHeatMap.siNC_TET1_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_TET1_summits(28450)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siCHD5_CTCF_rep123_TvsC"

```



## siNC & siMBD2 IP_VS_input
### Use each individual bam separately

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/
mkdir plotHeatMap_for_siNC_and_siMBD2
cd plotHeatMap_for_siNC_and_siMBD2
# 创建 bam 和 bai 文件链接
## bam
ln -s ../HEK293T_siNC_input_rep1.bowtie2.sorted.bam HEK293T_siNC_input_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_input_rep2.bowtie2.sorted.bam HEK293T_siNC_input_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_input_rep3.bowtie2.sorted.bam HEK293T_siNC_input_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD2_input_rep1.bowtie2.sorted.bam HEK293T_siMBD2_input_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD2_input_rep2.bowtie2.sorted.bam HEK293T_siMBD2_input_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD2_input_rep3.bowtie2.sorted.bam HEK293T_siMBD2_input_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD2_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siMBD2_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD2_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siMBD2_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD2_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siMBD2_IP_CTCF_rep3.bowtie2.sorted.bam
## bam.bai
ln -s ../HEK293T_siNC_input_rep1.bowtie2.sorted.bam.bai HEK293T_siNC_input_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_input_rep2.bowtie2.sorted.bam.bai HEK293T_siNC_input_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_input_rep3.bowtie2.sorted.bam.bai HEK293T_siNC_input_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD2_input_rep1.bowtie2.sorted.bam.bai HEK293T_siMBD2_input_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD2_input_rep2.bowtie2.sorted.bam.bai HEK293T_siMBD2_input_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD2_input_rep3.bowtie2.sorted.bam.bai HEK293T_siMBD2_input_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD2_IP_CTCF_rep1.bowtie2.sorted.bam.bai HEK293T_siMBD2_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD2_IP_CTCF_rep2.bowtie2.sorted.bam.bai HEK293T_siMBD2_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD2_IP_CTCF_rep3.bowtie2.sorted.bam.bai HEK293T_siMBD2_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *bowtie2.sorted.bam

# 创建 siNC CTCF peak summit 文件连接
ln -s ../macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_siMBD2_rep123_together/IP_CTCF_rep123_TvsC_siNC.only.summit.bed siNC_only.CTCF.summit.c3su.bed
ln -s ../macs2_peak_calling_2_for_siMBD2_rep123_together/IP_CTCF_rep123_TvsC_siNC_and_siMBD2.common.summit.bed siNC_siMBD2_common.CTCF.summit.c3su.bed
ln -s ../macs2_peak_calling_2_for_siMBD2_rep123_together/IP_CTCF_rep123_TvsC_siMBD2.only.summit.bed siMBD2_only.CTCF.summit.c3su.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_siNC_vs_siMBD2_bdgdiff/ROI_siNC.CTCF.summit.sorted.bed ROI_siNC.CTCF.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_siNC_vs_siMBD2_bdgdiff/RIB_siNC_and_siMBD2.CTCF.summit.sorted.bed RIB_siNC_and_siMBD2.CTCF.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_siNC_vs_siMBD2_bdgdiff/ROI_siMBD2.CTCF.summit.sorted.bed ROI_siMBD2.CTCF.summit.sorted.bed
wc -l *.bed
#     0 ROI_siNC.CTCF.summit.sorted.bed
#  6180 RIB_siNC_and_siMBD2.CTCF.summit.sorted.bed
#   663 ROI_siMBD2.CTCF.summit.sorted.bed
# 13936 siNC_only.CTCF.summit.c3su.bed
# 24236 siNC_siMBD2_common.CTCF.summit.c3su.bed
#  2959 siMBD2_only.CTCF.summit.c3su.bed
# 38169 siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed

```

#### Create log2 transformed bigwig files

```bash
# 构造 7_run_bamCompare_for_siNC_siMBD2_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# bamCompare -b1 HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam -b2 HEK293T_siNC_input_rep1.bowtie2.sorted.bam -o siNC_IP_CTCF_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam -b2 HEK293T_siNC_input_rep2.bowtie2.sorted.bam -o siNC_IP_CTCF_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam -b2 HEK293T_siNC_input_rep3.bowtie2.sorted.bam -o siNC_IP_CTCF_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siMBD2_IP_CTCF_rep1.bowtie2.sorted.bam -b2 HEK293T_siMBD2_input_rep1.bowtie2.sorted.bam -o siMBD2_IP_CTCF_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siMBD2_IP_CTCF_rep2.bowtie2.sorted.bam -b2 HEK293T_siMBD2_input_rep2.bowtie2.sorted.bam -o siMBD2_IP_CTCF_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siMBD2_IP_CTCF_rep3.bowtie2.sorted.bam -b2 HEK293T_siMBD2_input_rep3.bowtie2.sorted.bam -o siMBD2_IP_CTCF_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
nohup bash 7_run_bamCompare_for_siNC_siMBD2_each_repeat.sh >7_run_bamCompare_for_siNC_siMBD2_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bamCompare_for_siNC_siMBD2_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siMBD2
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siMBD2_IP_CTCF_rep1_TvsC.bigwig siMBD2_IP_CTCF_rep2_TvsC.bigwig siMBD2_IP_CTCF_rep3_TvsC.bigwig -R siNC_only.CTCF.summit.c3su.bed siNC_siMBD2_common.CTCF.summit.c3su.bed siMBD2_only.CTCF.summit.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siMBD2_all_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siMBD2_all_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_siNC_siMBD2_all_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_only(13936)" "common_summits(24236)" "siMBD2_only(2959)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siMBD2_CTCF_rep1_TvsC" "siMBD2_CTCF_rep2_TvsC" "siMBD2_CTCF_rep3_TvsC"

# use summits generated by bdgdiff
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siMBD2
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siMBD2_IP_CTCF_rep1_TvsC.bigwig siMBD2_IP_CTCF_rep2_TvsC.bigwig siMBD2_IP_CTCF_rep3_TvsC.bigwig -R RIB_siNC_and_siMBD2.CTCF.summit.sorted.bed ROI_siMBD2.CTCF.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siMBD2_all_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siMBD2_all_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_siNC_siMBD2_all_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "RIB_siNC_and_siMBD2(6180)" "ROI_siMBD2(663)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siMBD2_CTCF_rep1_TvsC" "siMBD2_CTCF_rep2_TvsC" "siMBD2_CTCF_rep3_TvsC"

# use siNC CTCF summits
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siMBD2
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siMBD2_IP_CTCF_rep1_TvsC.bigwig siMBD2_IP_CTCF_rep2_TvsC.bigwig siMBD2_IP_CTCF_rep3_TvsC.bigwig -R siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siMBD2_all_samples_for_plotHeatMap.siNC_CTCF_summit.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siMBD2_all_samples_for_plotHeatMap.siNC_CTCF_summit.gz -o bamCompare_matrix_of_siNC_siMBD2_all_samples_for_plotHeatMap.siNC_CTCF_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_CTCF(38169)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siMBD2_CTCF_rep1_TvsC" "siMBD2_CTCF_rep2_TvsC" "siMBD2_CTCF_rep3_TvsC"

# use 293T 2w CTCF summit
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siMBD2
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siMBD2_IP_CTCF_rep1_TvsC.bigwig siMBD2_IP_CTCF_rep2_TvsC.bigwig siMBD2_IP_CTCF_rep3_TvsC.bigwig -R ../CTCF.summit.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siMBD2_all_samples_for_plotHeatMap.293T2wCTCF.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siMBD2_all_samples_for_plotHeatMap.293T2wCTCF.gz -o bamCompare_matrix_of_siNC_siMBD2_all_samples_for_plotHeatMap.293T2wCTCF.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_summits(20000)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siMBD2_CTCF_rep1_TvsC" "siMBD2_CTCF_rep2_TvsC" "siMBD2_CTCF_rep3_TvsC"

# use siNC TET1 summits (28450)
# wc -l /home/sunwenju/projects_on_940xa/2023122801_WN_siNuRD_TET1_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_TET1_rep123_TvsC/siNC_IP_TET1_rep123_TvsC_summits.c3su.bed
# 28450 /home/sunwenju/projects_on_940xa/2023122801_WN_siNuRD_TET1_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_TET1_rep123_TvsC/siNC_IP_TET1_rep123_TvsC_summits.c3su.bed
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siMBD2
# 创建 siNC TET1 peak summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2023122801_WN_siNuRD_TET1_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_TET1_rep123_TvsC/siNC_IP_TET1_rep123_TvsC_summits.c3su.bed siNC_IP_TET1_rep123_TvsC_summits.c3su.bed
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siMBD2_IP_CTCF_rep1_TvsC.bigwig siMBD2_IP_CTCF_rep2_TvsC.bigwig siMBD2_IP_CTCF_rep3_TvsC.bigwig -R siNC_IP_TET1_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siMBD2_all_samples_for_plotHeatMap.siNC_TET1_summit.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siMBD2_all_samples_for_plotHeatMap.siNC_TET1_summit.gz -o bamCompare_matrix_of_siNC_siMBD2_all_samples_for_plotHeatMap.siNC_TET1_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_TET1_summits(28450)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siMBD2_CTCF_rep1_TvsC" "siMBD2_CTCF_rep2_TvsC" "siMBD2_CTCF_rep3_TvsC"

```

### Use rep123 merged bam

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siMBD2/
# 创建 rep123 merged bam 和 bai 文件链接
## bam
ln -s ../HEK293T_siNC_input_rep123.merged.sorted.bam HEK293T_siNC_input_rep123.merged.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam
ln -s ../HEK293T_siMBD2_input_rep123.merged.sorted.bam HEK293T_siMBD2_input_rep123.merged.sorted.bam
ln -s ../HEK293T_siMBD2_IP_CTCF_rep123.merged.sorted.bam HEK293T_siMBD2_IP_CTCF_rep123.merged.sorted.bam
## bam.bai
ln -s ../HEK293T_siNC_input_rep123.merged.sorted.bam.bai HEK293T_siNC_input_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_siMBD2_input_rep123.merged.sorted.bam.bai HEK293T_siMBD2_input_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_siMBD2_IP_CTCF_rep123.merged.sorted.bam.bai HEK293T_siMBD2_IP_CTCF_rep123.merged.sorted.bam.bai
ls -1 *rep123.merged.sorted.bam

```

#### Create log2 transformed bigwig files

```bash
# 构造 7_run_bamCompare_for_siNC_siMBD2_rep123_merged_bam.sh 内容如下
# # Compare for rep123 merged bam
# bamCompare -b1 HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam -b2 HEK293T_siNC_input_rep123.merged.sorted.bam -o siNC_IP_CTCF_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siMBD2_IP_CTCF_rep123.merged.sorted.bam -b2 HEK293T_siMBD2_input_rep123.merged.sorted.bam -o siMBD2_IP_CTCF_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
nohup bash 7_run_bamCompare_for_siNC_siMBD2_rep123_merged_bam.sh >7_run_bamCompare_for_siNC_siMBD2_rep123_merged_bam.sh.runlog 2>&1 &
tail -f 7_run_bamCompare_for_siNC_siMBD2_rep123_merged_bam.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siMBD2
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siMBD2_IP_CTCF_rep123_merged_TvsC.bigwig -R siNC_only.CTCF.summit.c3su.bed siNC_siMBD2_common.CTCF.summit.c3su.bed siMBD2_only.CTCF.summit.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siMBD2_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siMBD2_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_siNC_siMBD2_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_only(13936)" "common_summits(24236)" "siMBD2_only(2959)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siMBD2_CTCF_rep123_TvsC"

# use summits generated by bdgdiff
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siMBD2
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siMBD2_IP_CTCF_rep123_merged_TvsC.bigwig -R RIB_siNC_and_siMBD2.CTCF.summit.sorted.bed ROI_siMBD2.CTCF.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siMBD2_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siMBD2_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_siNC_siMBD2_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "RIB_siNC_and_siMBD2(6180)" "ROI_siMBD2(663)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siMBD2_CTCF_rep123_TvsC"

# use siNC CTCF summits
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siMBD2
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siMBD2_IP_CTCF_rep123_merged_TvsC.bigwig -R siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siMBD2_rep123_merged_for_plotHeatMap.siNC_CTCF_summit.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siMBD2_rep123_merged_for_plotHeatMap.siNC_CTCF_summit.gz -o bamCompare_matrix_of_siNC_siMBD2_rep123_merged_for_plotHeatMap.siNC_CTCF_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_CTCF(38169)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siMBD2_CTCF_rep123_TvsC"

# use 293T 2w CTCF summit
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siMBD2
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siMBD2_IP_CTCF_rep123_merged_TvsC.bigwig -R ../CTCF.summit.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siMBD2_rep123_merged_for_plotHeatMap.293T2wCTCF.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siMBD2_rep123_merged_for_plotHeatMap.293T2wCTCF.gz -o bamCompare_matrix_of_siNC_siMBD2_rep123_merged_for_plotHeatMap.293T2wCTCF.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_summits(20000)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siMBD2_CTCF_rep123_TvsC"

# use siNC TET1 summits (28450)
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siMBD2
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siMBD2_IP_CTCF_rep123_merged_TvsC.bigwig -R siNC_IP_TET1_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siMBD2_rep123_merged_for_plotHeatMap.siNC_TET1_summit.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siMBD2_rep123_merged_for_plotHeatMap.siNC_TET1_summit.gz -o bamCompare_matrix_of_siNC_siMBD2_rep123_merged_for_plotHeatMap.siNC_TET1_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_TET1_summits(28450)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siMBD2_CTCF_rep123_TvsC"


```



## siNC & siMBD3 IP_VS_input
### Use each individual bam separately

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/
mkdir plotHeatMap_for_siNC_and_siMBD3
cd plotHeatMap_for_siNC_and_siMBD3
# 创建 bam 和 bai 文件链接
## bam
ln -s ../HEK293T_siNC_input_rep1.bowtie2.sorted.bam HEK293T_siNC_input_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_input_rep2.bowtie2.sorted.bam HEK293T_siNC_input_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_input_rep3.bowtie2.sorted.bam HEK293T_siNC_input_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD3_input_rep1.bowtie2.sorted.bam HEK293T_siMBD3_input_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD3_input_rep2.bowtie2.sorted.bam HEK293T_siMBD3_input_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD3_input_rep3.bowtie2.sorted.bam HEK293T_siMBD3_input_rep3.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD3_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_siMBD3_IP_CTCF_rep1.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD3_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_siMBD3_IP_CTCF_rep2.bowtie2.sorted.bam
ln -s ../HEK293T_siMBD3_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_siMBD3_IP_CTCF_rep3.bowtie2.sorted.bam
## bam.bai
ln -s ../HEK293T_siNC_input_rep1.bowtie2.sorted.bam.bai HEK293T_siNC_input_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_input_rep2.bowtie2.sorted.bam.bai HEK293T_siNC_input_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_input_rep3.bowtie2.sorted.bam.bai HEK293T_siNC_input_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD3_input_rep1.bowtie2.sorted.bam.bai HEK293T_siMBD3_input_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD3_input_rep2.bowtie2.sorted.bam.bai HEK293T_siMBD3_input_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD3_input_rep3.bowtie2.sorted.bam.bai HEK293T_siMBD3_input_rep3.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD3_IP_CTCF_rep1.bowtie2.sorted.bam.bai HEK293T_siMBD3_IP_CTCF_rep1.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD3_IP_CTCF_rep2.bowtie2.sorted.bam.bai HEK293T_siMBD3_IP_CTCF_rep2.bowtie2.sorted.bam.bai
ln -s ../HEK293T_siMBD3_IP_CTCF_rep3.bowtie2.sorted.bam.bai HEK293T_siMBD3_IP_CTCF_rep3.bowtie2.sorted.bam.bai
ls -1 *bowtie2.sorted.bam

# 创建 siNC CTCF peak summit 文件连接
ln -s ../macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_siMBD3_rep123_together/IP_CTCF_rep123_TvsC_siNC.only.summit.bed siNC_only.CTCF.summit.c3su.bed
ln -s ../macs2_peak_calling_2_for_siMBD3_rep123_together/IP_CTCF_rep123_TvsC_siNC_and_siMBD3.common.summit.bed siNC_siMBD3_common.CTCF.summit.c3su.bed
ln -s ../macs2_peak_calling_2_for_siMBD3_rep123_together/IP_CTCF_rep123_TvsC_siMBD3.only.summit.bed siMBD3_only.CTCF.summit.c3su.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_siNC_vs_siMBD3_bdgdiff/ROI_siNC.CTCF.summit.sorted.bed ROI_siNC.CTCF.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_siNC_vs_siMBD3_bdgdiff/RIB_siNC_and_siMBD3.CTCF.summit.sorted.bed RIB_siNC_and_siMBD3.CTCF.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_siNC_vs_siMBD3_bdgdiff/ROI_siMBD3.CTCF.summit.sorted.bed ROI_siMBD3.CTCF.summit.sorted.bed
wc -l *.bed
#     0 ROI_siNC.CTCF.summit.sorted.bed
#  1578 RIB_siNC_and_siMBD3.CTCF.summit.sorted.bed
#   983 ROI_siMBD3.CTCF.summit.sorted.bed
# 15319 siNC_only.CTCF.summit.c3su.bed
# 22856 siNC_siMBD3_common.CTCF.summit.c3su.bed
#  4125 siMBD3_only.CTCF.summit.c3su.bed
# 38169 siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed

```

#### Create log2 transformed bigwig files

```bash
# 构造 7_run_bamCompare_for_siNC_siMBD3_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# bamCompare -b1 HEK293T_siNC_IP_CTCF_rep1.bowtie2.sorted.bam -b2 HEK293T_siNC_input_rep1.bowtie2.sorted.bam -o siNC_IP_CTCF_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siNC_IP_CTCF_rep2.bowtie2.sorted.bam -b2 HEK293T_siNC_input_rep2.bowtie2.sorted.bam -o siNC_IP_CTCF_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siNC_IP_CTCF_rep3.bowtie2.sorted.bam -b2 HEK293T_siNC_input_rep3.bowtie2.sorted.bam -o siNC_IP_CTCF_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siMBD3_IP_CTCF_rep1.bowtie2.sorted.bam -b2 HEK293T_siMBD3_input_rep1.bowtie2.sorted.bam -o siMBD3_IP_CTCF_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siMBD3_IP_CTCF_rep2.bowtie2.sorted.bam -b2 HEK293T_siMBD3_input_rep2.bowtie2.sorted.bam -o siMBD3_IP_CTCF_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siMBD3_IP_CTCF_rep3.bowtie2.sorted.bam -b2 HEK293T_siMBD3_input_rep3.bowtie2.sorted.bam -o siMBD3_IP_CTCF_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
nohup bash 7_run_bamCompare_for_siNC_siMBD3_each_repeat.sh >7_run_bamCompare_for_siNC_siMBD3_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bamCompare_for_siNC_siMBD3_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siMBD3
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siMBD3_IP_CTCF_rep1_TvsC.bigwig siMBD3_IP_CTCF_rep2_TvsC.bigwig siMBD3_IP_CTCF_rep3_TvsC.bigwig -R siNC_only.CTCF.summit.c3su.bed siNC_siMBD3_common.CTCF.summit.c3su.bed siMBD3_only.CTCF.summit.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siMBD3_all_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siMBD3_all_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_siNC_siMBD3_all_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_only(15319)" "common_summits(22856)" "siMBD3_only(4125)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siMBD3_CTCF_rep1_TvsC" "siMBD3_CTCF_rep2_TvsC" "siMBD3_CTCF_rep3_TvsC"

# use summits generated by bdgdiff
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siMBD3
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siMBD3_IP_CTCF_rep1_TvsC.bigwig siMBD3_IP_CTCF_rep2_TvsC.bigwig siMBD3_IP_CTCF_rep3_TvsC.bigwig -R RIB_siNC_and_siMBD3.CTCF.summit.sorted.bed ROI_siMBD3.CTCF.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siMBD3_all_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siMBD3_all_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_siNC_siMBD3_all_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "RIB_siNC_and_siMBD3(1578)" "ROI_siMBD3(983)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siMBD3_CTCF_rep1_TvsC" "siMBD3_CTCF_rep2_TvsC" "siMBD3_CTCF_rep3_TvsC"

# use siNC CTCF summits
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siMBD3
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siMBD3_IP_CTCF_rep1_TvsC.bigwig siMBD3_IP_CTCF_rep2_TvsC.bigwig siMBD3_IP_CTCF_rep3_TvsC.bigwig -R siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siMBD3_all_samples_for_plotHeatMap.siNC_CTCF_summit.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siMBD3_all_samples_for_plotHeatMap.siNC_CTCF_summit.gz -o bamCompare_matrix_of_siNC_siMBD3_all_samples_for_plotHeatMap.siNC_CTCF_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_CTCF(38169)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siMBD3_CTCF_rep1_TvsC" "siMBD3_CTCF_rep2_TvsC" "siMBD3_CTCF_rep3_TvsC"

# use 293T 2w CTCF summit
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siMBD3
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siMBD3_IP_CTCF_rep1_TvsC.bigwig siMBD3_IP_CTCF_rep2_TvsC.bigwig siMBD3_IP_CTCF_rep3_TvsC.bigwig -R ../CTCF.summit.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siMBD3_all_samples_for_plotHeatMap.293T2wCTCF.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siMBD3_all_samples_for_plotHeatMap.293T2wCTCF.gz -o bamCompare_matrix_of_siNC_siMBD3_all_samples_for_plotHeatMap.293T2wCTCF.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_summits(20000)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siMBD3_CTCF_rep1_TvsC" "siMBD3_CTCF_rep2_TvsC" "siMBD3_CTCF_rep3_TvsC"

# use siNC TET1 summits (28450)
# wc -l /home/sunwenju/projects_on_940xa/2023122801_WN_siNuRD_TET1_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_TET1_rep123_TvsC/siNC_IP_TET1_rep123_TvsC_summits.c3su.bed
# 28450 /home/sunwenju/projects_on_940xa/2023122801_WN_siNuRD_TET1_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_TET1_rep123_TvsC/siNC_IP_TET1_rep123_TvsC_summits.c3su.bed
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siMBD3
# 创建 siNC TET1 peak summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2023122801_WN_siNuRD_TET1_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_TET1_rep123_TvsC/siNC_IP_TET1_rep123_TvsC_summits.c3su.bed siNC_IP_TET1_rep123_TvsC_summits.c3su.bed
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siMBD3_IP_CTCF_rep1_TvsC.bigwig siMBD3_IP_CTCF_rep2_TvsC.bigwig siMBD3_IP_CTCF_rep3_TvsC.bigwig -R siNC_IP_TET1_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siMBD3_all_samples_for_plotHeatMap.siNC_TET1_summit.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siMBD3_all_samples_for_plotHeatMap.siNC_TET1_summit.gz -o bamCompare_matrix_of_siNC_siMBD3_all_samples_for_plotHeatMap.siNC_TET1_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_TET1_summits(28450)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siMBD3_CTCF_rep1_TvsC" "siMBD3_CTCF_rep2_TvsC" "siMBD3_CTCF_rep3_TvsC"

```

### Use rep123 merged bam

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siMBD3/
# 创建 rep123 merged bam 和 bai 文件链接
## bam
ln -s ../HEK293T_siNC_input_rep123.merged.sorted.bam HEK293T_siNC_input_rep123.merged.sorted.bam
ln -s ../HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam
ln -s ../HEK293T_siMBD3_input_rep123.merged.sorted.bam HEK293T_siMBD3_input_rep123.merged.sorted.bam
ln -s ../HEK293T_siMBD3_IP_CTCF_rep123.merged.sorted.bam HEK293T_siMBD3_IP_CTCF_rep123.merged.sorted.bam
## bam.bai
ln -s ../HEK293T_siNC_input_rep123.merged.sorted.bam.bai HEK293T_siNC_input_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam.bai HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_siMBD3_input_rep123.merged.sorted.bam.bai HEK293T_siMBD3_input_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_siMBD3_IP_CTCF_rep123.merged.sorted.bam.bai HEK293T_siMBD3_IP_CTCF_rep123.merged.sorted.bam.bai
ls -1 *rep123.merged.sorted.bam

```

#### Create log2 transformed bigwig files

```bash
# 构造 7_run_bamCompare_for_siNC_siMBD3_rep123_merged_bam.sh 内容如下
# # Compare for rep123 merged bam
# bamCompare -b1 HEK293T_siNC_IP_CTCF_rep123.merged.sorted.bam -b2 HEK293T_siNC_input_rep123.merged.sorted.bam -o siNC_IP_CTCF_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_siMBD3_IP_CTCF_rep123.merged.sorted.bam -b2 HEK293T_siMBD3_input_rep123.merged.sorted.bam -o siMBD3_IP_CTCF_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
nohup bash 7_run_bamCompare_for_siNC_siMBD3_rep123_merged_bam.sh >7_run_bamCompare_for_siNC_siMBD3_rep123_merged_bam.sh.runlog 2>&1 &
tail -f 7_run_bamCompare_for_siNC_siMBD3_rep123_merged_bam.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siMBD3
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siMBD3_IP_CTCF_rep123_merged_TvsC.bigwig -R siNC_only.CTCF.summit.c3su.bed siNC_siMBD3_common.CTCF.summit.c3su.bed siMBD3_only.CTCF.summit.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siMBD3_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siMBD3_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_siNC_siMBD3_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_only(15319)" "common_summits(22856)" "siMBD3_only(4125)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siMBD3_CTCF_rep123_TvsC"

# use summits generated by bdgdiff
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siMBD3
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siMBD3_IP_CTCF_rep123_merged_TvsC.bigwig -R RIB_siNC_and_siMBD3.CTCF.summit.sorted.bed ROI_siMBD3.CTCF.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siMBD3_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siMBD3_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_siNC_siMBD3_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "RIB_siNC_and_siMBD3(1578)" "ROI_siMBD3(983)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siMBD3_CTCF_rep123_TvsC"


# use siNC CTCF summits
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siMBD3
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siMBD3_IP_CTCF_rep123_merged_TvsC.bigwig -R siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siMBD3_rep123_merged_for_plotHeatMap.siNC_CTCF_summit.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siMBD3_rep123_merged_for_plotHeatMap.siNC_CTCF_summit.gz -o bamCompare_matrix_of_siNC_siMBD3_rep123_merged_for_plotHeatMap.siNC_CTCF_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_CTCF(38169)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siMBD3_CTCF_rep123_TvsC"

# use 293T 2w CTCF summit
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siMBD3
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siMBD3_IP_CTCF_rep123_merged_TvsC.bigwig -R ../CTCF.summit.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siMBD3_rep123_merged_for_plotHeatMap.293T2wCTCF.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siMBD3_rep123_merged_for_plotHeatMap.293T2wCTCF.gz -o bamCompare_matrix_of_siNC_siMBD3_rep123_merged_for_plotHeatMap.293T2wCTCF.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_summits(20000)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siMBD3_CTCF_rep123_TvsC"

# use siNC TET1 summits (28450)
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_siNC_and_siMBD3
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siMBD3_IP_CTCF_rep123_merged_TvsC.bigwig -R siNC_IP_TET1_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNC_siMBD3_rep123_merged_for_plotHeatMap.siNC_TET1_summit.gz
plotHeatmap -m bamCompare_matrix_of_siNC_siMBD3_rep123_merged_for_plotHeatMap.siNC_TET1_summit.gz -o bamCompare_matrix_of_siNC_siMBD3_rep123_merged_for_plotHeatMap.siNC_TET1_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_TET1_summits(28450)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siMBD3_CTCF_rep123_TvsC"


```



## all samples IP_VS_input
### Use each individual bam separately

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq
mkdir plotHeatMap_for_all
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_all
ln -s ../plotHeatMap_for_siNC_and_siCHD3/siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD3/siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD3/siNC_IP_CTCF_rep3_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD3/siCHD3_IP_CTCF_rep1_TvsC.bigwig siCHD3_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD3/siCHD3_IP_CTCF_rep2_TvsC.bigwig siCHD3_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD3/siCHD3_IP_CTCF_rep3_TvsC.bigwig siCHD3_IP_CTCF_rep3_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD4/siCHD4_IP_CTCF_rep1_TvsC.bigwig siCHD4_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD4/siCHD4_IP_CTCF_rep2_TvsC.bigwig siCHD4_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD4/siCHD4_IP_CTCF_rep3_TvsC.bigwig siCHD4_IP_CTCF_rep3_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD5/siCHD5_IP_CTCF_rep1_TvsC.bigwig siCHD5_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD5/siCHD5_IP_CTCF_rep2_TvsC.bigwig siCHD5_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD5/siCHD5_IP_CTCF_rep3_TvsC.bigwig siCHD5_IP_CTCF_rep3_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siMBD2/siMBD2_IP_CTCF_rep1_TvsC.bigwig siMBD2_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siMBD2/siMBD2_IP_CTCF_rep2_TvsC.bigwig siMBD2_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siMBD2/siMBD2_IP_CTCF_rep3_TvsC.bigwig siMBD2_IP_CTCF_rep3_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siMBD3/siMBD3_IP_CTCF_rep1_TvsC.bigwig siMBD3_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siMBD3/siMBD3_IP_CTCF_rep2_TvsC.bigwig siMBD3_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siMBD3/siMBD3_IP_CTCF_rep3_TvsC.bigwig siMBD3_IP_CTCF_rep3_TvsC.bigwig
ls -1 *rep[123]_TvsC.bigwig

# 创建 siNC CTCF peak summit 文件连接
ln -s ../macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_CTCF_rep123_TvsC/siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed

```

#### computeMatrix and plotHeatmap

```bash
# use siNC CTCF summits
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_all/
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siCHD3_IP_CTCF_rep1_TvsC.bigwig siCHD3_IP_CTCF_rep2_TvsC.bigwig siCHD3_IP_CTCF_rep3_TvsC.bigwig siCHD4_IP_CTCF_rep1_TvsC.bigwig siCHD4_IP_CTCF_rep2_TvsC.bigwig siCHD4_IP_CTCF_rep3_TvsC.bigwig siCHD5_IP_CTCF_rep1_TvsC.bigwig siCHD5_IP_CTCF_rep2_TvsC.bigwig siCHD5_IP_CTCF_rep3_TvsC.bigwig siMBD2_IP_CTCF_rep1_TvsC.bigwig siMBD2_IP_CTCF_rep2_TvsC.bigwig siMBD2_IP_CTCF_rep3_TvsC.bigwig siMBD3_IP_CTCF_rep1_TvsC.bigwig siMBD3_IP_CTCF_rep2_TvsC.bigwig siMBD3_IP_CTCF_rep3_TvsC.bigwig -R siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_all_samples_for_plotHeatMap.siNC_CTCF_summit.gz
plotHeatmap -m bamCompare_matrix_of_all_samples_for_plotHeatMap.siNC_CTCF_summit.gz -o bamCompare_matrix_of_all_samples_for_plotHeatMap.siNC_CTCF_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_CTCF(38169)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siCHD3_CTCF_rep1_TvsC" "siCHD3_CTCF_rep2_TvsC" "siCHD3_CTCF_rep3_TvsC" "siCHD4_CTCF_rep1_TvsC" "siCHD4_CTCF_rep2_TvsC" "siCHD4_CTCF_rep3_TvsC" "siCHD5_CTCF_rep1_TvsC" "siCHD5_CTCF_rep2_TvsC" "siCHD5_CTCF_rep3_TvsC" "siMBD2_CTCF_rep1_TvsC" "siMBD2_CTCF_rep2_TvsC" "siMBD2_CTCF_rep3_TvsC" "siMBD3_CTCF_rep1_TvsC" "siMBD3_CTCF_rep2_TvsC" "siMBD3_CTCF_rep3_TvsC"

# use 293T 2w CTCF summit
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_all/
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siCHD3_IP_CTCF_rep1_TvsC.bigwig siCHD3_IP_CTCF_rep2_TvsC.bigwig siCHD3_IP_CTCF_rep3_TvsC.bigwig siCHD4_IP_CTCF_rep1_TvsC.bigwig siCHD4_IP_CTCF_rep2_TvsC.bigwig siCHD4_IP_CTCF_rep3_TvsC.bigwig siCHD5_IP_CTCF_rep1_TvsC.bigwig siCHD5_IP_CTCF_rep2_TvsC.bigwig siCHD5_IP_CTCF_rep3_TvsC.bigwig siMBD2_IP_CTCF_rep1_TvsC.bigwig siMBD2_IP_CTCF_rep2_TvsC.bigwig siMBD2_IP_CTCF_rep3_TvsC.bigwig siMBD3_IP_CTCF_rep1_TvsC.bigwig siMBD3_IP_CTCF_rep2_TvsC.bigwig siMBD3_IP_CTCF_rep3_TvsC.bigwig -R ../CTCF.summit.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_all_samples_for_plotHeatMap.293T2wCTCF.gz
plotHeatmap -m bamCompare_matrix_of_all_samples_for_plotHeatMap.293T2wCTCF.gz -o bamCompare_matrix_of_all_samples_for_plotHeatMap.293T2wCTCF.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_summits(20000)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siCHD3_CTCF_rep1_TvsC" "siCHD3_CTCF_rep2_TvsC" "siCHD3_CTCF_rep3_TvsC" "siCHD4_CTCF_rep1_TvsC" "siCHD4_CTCF_rep2_TvsC" "siCHD4_CTCF_rep3_TvsC" "siCHD5_CTCF_rep1_TvsC" "siCHD5_CTCF_rep2_TvsC" "siCHD5_CTCF_rep3_TvsC" "siMBD2_CTCF_rep1_TvsC" "siMBD2_CTCF_rep2_TvsC" "siMBD2_CTCF_rep3_TvsC" "siMBD3_CTCF_rep1_TvsC" "siMBD3_CTCF_rep2_TvsC" "siMBD3_CTCF_rep3_TvsC"

# use siNC TET1 summits (28450)
# wc -l /home/sunwenju/projects_on_940xa/2023122801_WN_siNuRD_TET1_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_TET1_rep123_TvsC/siNC_IP_TET1_rep123_TvsC_summits.c3su.bed
# 28450 /home/sunwenju/projects_on_940xa/2023122801_WN_siNuRD_TET1_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_TET1_rep123_TvsC/siNC_IP_TET1_rep123_TvsC_summits.c3su.bed
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_all/
# 创建 siNC TET1 peak summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2023122801_WN_siNuRD_TET1_ChIP_Seq/macs2_peak_calling_2_for_siNC_rep123_together/siNC_IP_TET1_rep123_TvsC/siNC_IP_TET1_rep123_TvsC_summits.c3su.bed siNC_IP_TET1_rep123_TvsC_summits.c3su.bed
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siCHD3_IP_CTCF_rep1_TvsC.bigwig siCHD3_IP_CTCF_rep2_TvsC.bigwig siCHD3_IP_CTCF_rep3_TvsC.bigwig siCHD4_IP_CTCF_rep1_TvsC.bigwig siCHD4_IP_CTCF_rep2_TvsC.bigwig siCHD4_IP_CTCF_rep3_TvsC.bigwig siCHD5_IP_CTCF_rep1_TvsC.bigwig siCHD5_IP_CTCF_rep2_TvsC.bigwig siCHD5_IP_CTCF_rep3_TvsC.bigwig siMBD2_IP_CTCF_rep1_TvsC.bigwig siMBD2_IP_CTCF_rep2_TvsC.bigwig siMBD2_IP_CTCF_rep3_TvsC.bigwig siMBD3_IP_CTCF_rep1_TvsC.bigwig siMBD3_IP_CTCF_rep2_TvsC.bigwig siMBD3_IP_CTCF_rep3_TvsC.bigwig -R siNC_IP_TET1_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_all_samples_for_plotHeatMap.siNC_TET1_summit.gz
plotHeatmap -m bamCompare_matrix_of_all_samples_for_plotHeatMap.siNC_TET1_summit.gz -o bamCompare_matrix_of_all_samples_for_plotHeatMap.siNC_TET1_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_TET1_summit(28450)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siCHD3_CTCF_rep1_TvsC" "siCHD3_CTCF_rep2_TvsC" "siCHD3_CTCF_rep3_TvsC" "siCHD4_CTCF_rep1_TvsC" "siCHD4_CTCF_rep2_TvsC" "siCHD4_CTCF_rep3_TvsC" "siCHD5_CTCF_rep1_TvsC" "siCHD5_CTCF_rep2_TvsC" "siCHD5_CTCF_rep3_TvsC" "siMBD2_CTCF_rep1_TvsC" "siMBD2_CTCF_rep2_TvsC" "siMBD2_CTCF_rep3_TvsC" "siMBD3_CTCF_rep1_TvsC" "siMBD3_CTCF_rep2_TvsC" "siMBD3_CTCF_rep3_TvsC"

```

### Use rep123 merged bam

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_all
ln -s ../plotHeatMap_for_siNC_and_siCHD3/siNC_IP_CTCF_rep123_merged_TvsC.bigwig siNC_IP_CTCF_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD3/siCHD3_IP_CTCF_rep123_merged_TvsC.bigwig siCHD3_IP_CTCF_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD4/siCHD4_IP_CTCF_rep123_merged_TvsC.bigwig siCHD4_IP_CTCF_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD5/siCHD5_IP_CTCF_rep123_merged_TvsC.bigwig siCHD5_IP_CTCF_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siMBD2/siMBD2_IP_CTCF_rep123_merged_TvsC.bigwig siMBD2_IP_CTCF_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siMBD3/siMBD3_IP_CTCF_rep123_merged_TvsC.bigwig siMBD3_IP_CTCF_rep123_merged_TvsC.bigwig
ls -1 *rep123_merged_TvsC.bigwig

```

#### computeMatrix and plotHeatmap

```bash
# use siNC CTCF summits
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_all/
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siCHD3_IP_CTCF_rep123_merged_TvsC.bigwig siCHD4_IP_CTCF_rep123_merged_TvsC.bigwig siCHD5_IP_CTCF_rep123_merged_TvsC.bigwig siMBD2_IP_CTCF_rep123_merged_TvsC.bigwig siMBD3_IP_CTCF_rep123_merged_TvsC.bigwig -R siNC_IP_CTCF_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.siNC_CTCF_summit.gz
plotHeatmap -m bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.siNC_CTCF_summit.gz -o bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.siNC_CTCF_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_CTCF(38169)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siCHD3_CTCF_rep123_TvsC" "siCHD4_CTCF_rep123_TvsC" "siCHD5_CTCF_rep123_TvsC" "siMBD2_CTCF_rep123_TvsC" "siMBD3_CTCF_rep123_TvsC"

# use 293T 2w CTCF summit
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_all/
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siCHD3_IP_CTCF_rep123_merged_TvsC.bigwig siCHD4_IP_CTCF_rep123_merged_TvsC.bigwig siCHD5_IP_CTCF_rep123_merged_TvsC.bigwig siMBD2_IP_CTCF_rep123_merged_TvsC.bigwig siMBD3_IP_CTCF_rep123_merged_TvsC.bigwig -R ../CTCF.summit.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.293T2wCTCF.gz
plotHeatmap -m bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.293T2wCTCF.gz -o bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.293T2wCTCF.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_summits(20000)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siCHD3_CTCF_rep123_TvsC" "siCHD4_CTCF_rep123_TvsC" "siCHD5_CTCF_rep123_TvsC" "siMBD2_CTCF_rep123_TvsC" "siMBD3_CTCF_rep123_TvsC"

# use siNC TET1 summits (28450)
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_all/
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siCHD3_IP_CTCF_rep123_merged_TvsC.bigwig siCHD4_IP_CTCF_rep123_merged_TvsC.bigwig siCHD5_IP_CTCF_rep123_merged_TvsC.bigwig siMBD2_IP_CTCF_rep123_merged_TvsC.bigwig siMBD3_IP_CTCF_rep123_merged_TvsC.bigwig -R siNC_IP_TET1_rep123_TvsC_summits.c3su.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.siNC_TET1_summit.gz
plotHeatmap -m bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.siNC_TET1_summit.gz -o bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.siNC_TET1_summit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "siNC_TET1_summit(28450)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siCHD3_CTCF_rep123_TvsC" "siCHD4_CTCF_rep123_TvsC" "siCHD5_CTCF_rep123_TvsC" "siMBD2_CTCF_rep123_TvsC" "siMBD3_CTCF_rep123_TvsC"

```



# plot reads enrichment use public CTCF rep12 union summit

##  data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/
mkdir plotHeatMap_for_all_use_public_CTCF_rep12_union_summit
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_all_use_public_CTCF_rep12_union_summit

# 所需的几类 public CTCF rep12 union summit bed 文件在如下分析文件夹中已经构造过了
# /home/sunwenju/projects_on_940xa/2022072201_CEF_FSR_WN_BGI_data/04_WN_BPK25_IP_CTCF_MBD3_ChIP_new_input_downsp/plotHeatMap_use_bamCompare_bigwig_for_CTCF_public_rep12_union
# 为所需的bed文件创建软连接创建 
ln -s /home/sunwenju/projects_on_940xa/2022072201_CEF_FSR_WN_BGI_data/04_WN_BPK25_IP_CTCF_MBD3_ChIP_new_input_downsp/plotHeatMap_use_bamCompare_bigwig_for_CTCF_public_rep12_union/CTCF_rep12_narrowPeak.summit.union.all.sorted.bed CTCF_rep12_narrowPeak.summit.union.all.sorted.bed
ln -s /home/sunwenju/projects_on_940xa/2022072201_CEF_FSR_WN_BGI_data/04_WN_BPK25_IP_CTCF_MBD3_ChIP_new_input_downsp/plotHeatMap_use_bamCompare_bigwig_for_CTCF_public_rep12_union/CTCF_rep12_union.rep12_common.summit.sorted.bed CTCF_rep12_union.rep12_common.summit.sorted.bed
ln -s /home/sunwenju/projects_on_940xa/2022072201_CEF_FSR_WN_BGI_data/04_WN_BPK25_IP_CTCF_MBD3_ChIP_new_input_downsp/plotHeatMap_use_bamCompare_bigwig_for_CTCF_public_rep12_union/CTCF_rep12_union.rep1_only.summit.sorted.bed CTCF_rep12_union.rep1_only.summit.sorted.bed
ln -s /home/sunwenju/projects_on_940xa/2022072201_CEF_FSR_WN_BGI_data/04_WN_BPK25_IP_CTCF_MBD3_ChIP_new_input_downsp/plotHeatMap_use_bamCompare_bigwig_for_CTCF_public_rep12_union/CTCF_rep12_union.rep2_only.summit.sorted.bed CTCF_rep12_union.rep2_only.summit.sorted.bed
ln -s /home/sunwenju/projects_on_940xa/2022072201_CEF_FSR_WN_BGI_data/04_WN_BPK25_IP_CTCF_MBD3_ChIP_new_input_downsp/plotHeatMap_use_bamCompare_bigwig_for_CTCF_public_rep12_union/CTCF_rep12_union.phq1.summit.sorted.bed CTCF_rep12_union.phq1.summit.sorted.bed
ln -s /home/sunwenju/projects_on_940xa/2022072201_CEF_FSR_WN_BGI_data/04_WN_BPK25_IP_CTCF_MBD3_ChIP_new_input_downsp/plotHeatMap_use_bamCompare_bigwig_for_CTCF_public_rep12_union/CTCF_rep12_union.phq2.summit.sorted.bed CTCF_rep12_union.phq2.summit.sorted.bed
ln -s /home/sunwenju/projects_on_940xa/2022072201_CEF_FSR_WN_BGI_data/04_WN_BPK25_IP_CTCF_MBD3_ChIP_new_input_downsp/plotHeatMap_use_bamCompare_bigwig_for_CTCF_public_rep12_union/CTCF_rep12_union.phq3.summit.sorted.bed CTCF_rep12_union.phq3.summit.sorted.bed
ln -s /home/sunwenju/projects_on_940xa/2022072201_CEF_FSR_WN_BGI_data/04_WN_BPK25_IP_CTCF_MBD3_ChIP_new_input_downsp/plotHeatMap_use_bamCompare_bigwig_for_CTCF_public_rep12_union/CTCF_rep12_union.phq4.summit.sorted.bed CTCF_rep12_union.phq4.summit.sorted.bed
ls -1 *.bed
wc -l *.bed
# 97424 CTCF_rep12_narrowPeak.summit.union.all.sorted.bed
# 31737 CTCF_rep12_union.rep12_common.summit.sorted.bed
# 47355 CTCF_rep12_union.rep1_only.summit.sorted.bed
# 18332 CTCF_rep12_union.rep2_only.summit.sorted.bed
# 24356 CTCF_rep12_union.phq1.summit.sorted.bed
# 24356 CTCF_rep12_union.phq2.summit.sorted.bed
# 24356 CTCF_rep12_union.phq3.summit.sorted.bed
# 24356 CTCF_rep12_union.phq4.summit.sorted.bed

# 创建全部 bigwig 文件连接
ln -s ../plotHeatMap_for_siNC_and_siCHD3/siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD3/siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD3/siNC_IP_CTCF_rep3_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD3/siCHD3_IP_CTCF_rep1_TvsC.bigwig siCHD3_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD3/siCHD3_IP_CTCF_rep2_TvsC.bigwig siCHD3_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD3/siCHD3_IP_CTCF_rep3_TvsC.bigwig siCHD3_IP_CTCF_rep3_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD4/siCHD4_IP_CTCF_rep1_TvsC.bigwig siCHD4_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD4/siCHD4_IP_CTCF_rep2_TvsC.bigwig siCHD4_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD4/siCHD4_IP_CTCF_rep3_TvsC.bigwig siCHD4_IP_CTCF_rep3_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD5/siCHD5_IP_CTCF_rep1_TvsC.bigwig siCHD5_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD5/siCHD5_IP_CTCF_rep2_TvsC.bigwig siCHD5_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD5/siCHD5_IP_CTCF_rep3_TvsC.bigwig siCHD5_IP_CTCF_rep3_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siMBD2/siMBD2_IP_CTCF_rep1_TvsC.bigwig siMBD2_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siMBD2/siMBD2_IP_CTCF_rep2_TvsC.bigwig siMBD2_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siMBD2/siMBD2_IP_CTCF_rep3_TvsC.bigwig siMBD2_IP_CTCF_rep3_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siMBD3/siMBD3_IP_CTCF_rep1_TvsC.bigwig siMBD3_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siMBD3/siMBD3_IP_CTCF_rep2_TvsC.bigwig siMBD3_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siMBD3/siMBD3_IP_CTCF_rep3_TvsC.bigwig siMBD3_IP_CTCF_rep3_TvsC.bigwig
ls -1 *rep[123]_TvsC.bigwig

ln -s ../plotHeatMap_for_siNC_and_siCHD3/siNC_IP_CTCF_rep123_merged_TvsC.bigwig siNC_IP_CTCF_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD3/siCHD3_IP_CTCF_rep123_merged_TvsC.bigwig siCHD3_IP_CTCF_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD4/siCHD4_IP_CTCF_rep123_merged_TvsC.bigwig siCHD4_IP_CTCF_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siCHD5/siCHD5_IP_CTCF_rep123_merged_TvsC.bigwig siCHD5_IP_CTCF_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siMBD2/siMBD2_IP_CTCF_rep123_merged_TvsC.bigwig siMBD2_IP_CTCF_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_for_siNC_and_siMBD3/siMBD3_IP_CTCF_rep123_merged_TvsC.bigwig siMBD3_IP_CTCF_rep123_merged_TvsC.bigwig
ls -1 *rep123_merged_TvsC.bigwig

```

## plot use CTCF rep12 union all summits

```bash
# plot use CTCF rep12 union all summits
## compute and plot use rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siCHD3_IP_CTCF_rep1_TvsC.bigwig siCHD3_IP_CTCF_rep2_TvsC.bigwig siCHD3_IP_CTCF_rep3_TvsC.bigwig siCHD4_IP_CTCF_rep1_TvsC.bigwig siCHD4_IP_CTCF_rep2_TvsC.bigwig siCHD4_IP_CTCF_rep3_TvsC.bigwig siCHD5_IP_CTCF_rep1_TvsC.bigwig siCHD5_IP_CTCF_rep2_TvsC.bigwig siCHD5_IP_CTCF_rep3_TvsC.bigwig siMBD2_IP_CTCF_rep1_TvsC.bigwig siMBD2_IP_CTCF_rep2_TvsC.bigwig siMBD2_IP_CTCF_rep3_TvsC.bigwig siMBD3_IP_CTCF_rep1_TvsC.bigwig siMBD3_IP_CTCF_rep2_TvsC.bigwig siMBD3_IP_CTCF_rep3_TvsC.bigwig -R CTCF_rep12_narrowPeak.summit.union.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNuRD_all_samples_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz
plotHeatmap -m bamCompare_matrix_of_siNuRD_all_samples_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz -o bamCompare_matrix_of_siNuRD_all_samples_for_plotHeatMap.public_CTCF_rep12_union_all_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_rep12_union_all_summits(97424)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siCHD3_CTCF_rep1_TvsC" "siCHD3_CTCF_rep2_TvsC" "siCHD3_CTCF_rep3_TvsC" "siCHD4_CTCF_rep1_TvsC" "siCHD4_CTCF_rep2_TvsC" "siCHD4_CTCF_rep3_TvsC" "siCHD5_CTCF_rep1_TvsC" "siCHD5_CTCF_rep2_TvsC" "siCHD5_CTCF_rep3_TvsC" "siMBD2_CTCF_rep1_TvsC" "siMBD2_CTCF_rep2_TvsC" "siMBD2_CTCF_rep3_TvsC" "siMBD3_CTCF_rep1_TvsC" "siMBD3_CTCF_rep2_TvsC" "siMBD3_CTCF_rep3_TvsC"

## compute and plot use rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siCHD3_IP_CTCF_rep123_merged_TvsC.bigwig siCHD4_IP_CTCF_rep123_merged_TvsC.bigwig siCHD5_IP_CTCF_rep123_merged_TvsC.bigwig siMBD2_IP_CTCF_rep123_merged_TvsC.bigwig siMBD3_IP_CTCF_rep123_merged_TvsC.bigwig -R CTCF_rep12_narrowPeak.summit.union.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNuRD_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz
plotHeatmap -m bamCompare_matrix_of_siNuRD_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz -o bamCompare_matrix_of_siNuRD_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_rep12_union_all_summits(97424)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siCHD3_CTCF_rep123_TvsC" "siCHD4_CTCF_rep123_TvsC" "siCHD5_CTCF_rep123_TvsC" "siMBD2_CTCF_rep123_TvsC" "siMBD3_CTCF_rep123_TvsC"

```

## plot use union  source 3 class summit

```bash
# use CTCF_rep12_union rep1_only、rep12_common、rep2_only
## compute and plot use rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siCHD3_IP_CTCF_rep1_TvsC.bigwig siCHD3_IP_CTCF_rep2_TvsC.bigwig siCHD3_IP_CTCF_rep3_TvsC.bigwig siCHD4_IP_CTCF_rep1_TvsC.bigwig siCHD4_IP_CTCF_rep2_TvsC.bigwig siCHD4_IP_CTCF_rep3_TvsC.bigwig siCHD5_IP_CTCF_rep1_TvsC.bigwig siCHD5_IP_CTCF_rep2_TvsC.bigwig siCHD5_IP_CTCF_rep3_TvsC.bigwig siMBD2_IP_CTCF_rep1_TvsC.bigwig siMBD2_IP_CTCF_rep2_TvsC.bigwig siMBD2_IP_CTCF_rep3_TvsC.bigwig siMBD3_IP_CTCF_rep1_TvsC.bigwig siMBD3_IP_CTCF_rep2_TvsC.bigwig siMBD3_IP_CTCF_rep3_TvsC.bigwig -R CTCF_rep12_union.rep1_only.summit.sorted.bed CTCF_rep12_union.rep12_common.summit.sorted.bed CTCF_rep12_union.rep2_only.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNuRD_all_samples_for_plotHeatMap.public_CTCF_rep12_union_3source.gz
plotHeatmap -m bamCompare_matrix_of_siNuRD_all_samples_for_plotHeatMap.public_CTCF_rep12_union_3source.gz -o bamCompare_matrix_of_siNuRD_all_samples_for_plotHeatMap.public_CTCF_rep12_union_3source.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "rep1_only(47355)_summit" "rep12_common_summit(31737)" "rep2_only_summit(18332)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siCHD3_CTCF_rep1_TvsC" "siCHD3_CTCF_rep2_TvsC" "siCHD3_CTCF_rep3_TvsC" "siCHD4_CTCF_rep1_TvsC" "siCHD4_CTCF_rep2_TvsC" "siCHD4_CTCF_rep3_TvsC" "siCHD5_CTCF_rep1_TvsC" "siCHD5_CTCF_rep2_TvsC" "siCHD5_CTCF_rep3_TvsC" "siMBD2_CTCF_rep1_TvsC" "siMBD2_CTCF_rep2_TvsC" "siMBD2_CTCF_rep3_TvsC" "siMBD3_CTCF_rep1_TvsC" "siMBD3_CTCF_rep2_TvsC" "siMBD3_CTCF_rep3_TvsC"

## compute and plot use rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siCHD3_IP_CTCF_rep123_merged_TvsC.bigwig siCHD4_IP_CTCF_rep123_merged_TvsC.bigwig siCHD5_IP_CTCF_rep123_merged_TvsC.bigwig siMBD2_IP_CTCF_rep123_merged_TvsC.bigwig siMBD3_IP_CTCF_rep123_merged_TvsC.bigwig -R CTCF_rep12_union.rep1_only.summit.sorted.bed CTCF_rep12_union.rep12_common.summit.sorted.bed CTCF_rep12_union.rep2_only.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNuRD_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.gz
plotHeatmap -m bamCompare_matrix_of_siNuRD_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.gz -o bamCompare_matrix_of_siNuRD_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "rep1_only(47355)_summit" "rep12_common_summit(31737)" "rep2_only_summit(18332)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siCHD3_CTCF_rep123_TvsC" "siCHD4_CTCF_rep123_TvsC" "siCHD5_CTCF_rep123_TvsC" "siMBD2_CTCF_rep123_TvsC" "siMBD3_CTCF_rep123_TvsC"

```

## plot use union  4 class summit ( peak height quartile)

```bash
# use TCF_rep12_union.phq1234.summit.sorted.bed
## compute and plot use rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep1_TvsC.bigwig siNC_IP_CTCF_rep2_TvsC.bigwig siNC_IP_CTCF_rep3_TvsC.bigwig siCHD3_IP_CTCF_rep1_TvsC.bigwig siCHD3_IP_CTCF_rep2_TvsC.bigwig siCHD3_IP_CTCF_rep3_TvsC.bigwig siCHD4_IP_CTCF_rep1_TvsC.bigwig siCHD4_IP_CTCF_rep2_TvsC.bigwig siCHD4_IP_CTCF_rep3_TvsC.bigwig siCHD5_IP_CTCF_rep1_TvsC.bigwig siCHD5_IP_CTCF_rep2_TvsC.bigwig siCHD5_IP_CTCF_rep3_TvsC.bigwig siMBD2_IP_CTCF_rep1_TvsC.bigwig siMBD2_IP_CTCF_rep2_TvsC.bigwig siMBD2_IP_CTCF_rep3_TvsC.bigwig siMBD3_IP_CTCF_rep1_TvsC.bigwig siMBD3_IP_CTCF_rep2_TvsC.bigwig siMBD3_IP_CTCF_rep3_TvsC.bigwig -R CTCF_rep12_union.phq4.summit.sorted.bed CTCF_rep12_union.phq3.summit.sorted.bed CTCF_rep12_union.phq2.summit.sorted.bed CTCF_rep12_union.phq1.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNuRD_all_samples_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz
plotHeatmap -m bamCompare_matrix_of_siNuRD_all_samples_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz -o bamCompare_matrix_of_siNuRD_all_samples_for_plotHeatMap.public_CTCF_rep12_union.phq1234.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "CTCF_rep12_union_phq4_summit(24356)" "CTCF_rep12_union_phq3_summit(24356)" "CTCF_rep12_union_phq2_summit(24356)" "CTCF_rep12_union_phq1_summit(24356)" --samplesLabel "siNC_CTCF_rep1_TvsC" "siNC_CTCF_rep2_TvsC" "siNC_CTCF_rep3_TvsC" "siCHD3_CTCF_rep1_TvsC" "siCHD3_CTCF_rep2_TvsC" "siCHD3_CTCF_rep3_TvsC" "siCHD4_CTCF_rep1_TvsC" "siCHD4_CTCF_rep2_TvsC" "siCHD4_CTCF_rep3_TvsC" "siCHD5_CTCF_rep1_TvsC" "siCHD5_CTCF_rep2_TvsC" "siCHD5_CTCF_rep3_TvsC" "siMBD2_CTCF_rep1_TvsC" "siMBD2_CTCF_rep2_TvsC" "siMBD2_CTCF_rep3_TvsC" "siMBD3_CTCF_rep1_TvsC" "siMBD3_CTCF_rep2_TvsC" "siMBD3_CTCF_rep3_TvsC"

## compute and plot use rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023081601_WN_siNURD_CTCF_ChIP_Seq/plotHeatMap_for_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S siNC_IP_CTCF_rep123_merged_TvsC.bigwig siCHD3_IP_CTCF_rep123_merged_TvsC.bigwig siCHD4_IP_CTCF_rep123_merged_TvsC.bigwig siCHD5_IP_CTCF_rep123_merged_TvsC.bigwig siMBD2_IP_CTCF_rep123_merged_TvsC.bigwig siMBD3_IP_CTCF_rep123_merged_TvsC.bigwig -R CTCF_rep12_union.phq4.summit.sorted.bed CTCF_rep12_union.phq3.summit.sorted.bed CTCF_rep12_union.phq2.summit.sorted.bed CTCF_rep12_union.phq1.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_siNuRD_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz
plotHeatmap -m bamCompare_matrix_of_siNuRD_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz -o bamCompare_matrix_of_siNuRD_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "CTCF_rep12_union_phq4_summit(24356)" "CTCF_rep12_union_phq3_summit(24356)" "CTCF_rep12_union_phq2_summit(24356)" "CTCF_rep12_union_phq1_summit(24356)" --samplesLabel "siNC_CTCF_rep123_TvsC" "siCHD3_CTCF_rep123_TvsC" "siCHD4_CTCF_rep123_TvsC" "siCHD5_CTCF_rep123_TvsC" "siMBD2_CTCF_rep123_TvsC" "siMBD3_CTCF_rep123_TvsC"

```

