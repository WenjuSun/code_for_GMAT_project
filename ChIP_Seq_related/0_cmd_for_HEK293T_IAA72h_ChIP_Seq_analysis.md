

# QC

## raw data

```bash
# 创建项目分析目录
cd /home/sunwenju/projects_on_940xa/
mkdir 2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/

# 创建 new raw data 质量评估目录
mkdir 00_WN_293T_IAA72h_ChIP_Seq_new_raw_data_qc
# 先对新加测的数据的 raw data 进行质量评估
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/00_WN_293T_IAA72h_ChIP_Seq_new_raw_data_qc
# 根据样本信息表构造 new raw data 数据连接脚本，上传至该目录并运行
bash 00_make_WN_293T_IAA72h_ChIP_Seq_new_raw_data_fq_gz_file_links.sh
ls -1
#
mkdir 1_all_293T_IAA72h_ChIP_Seq_new_raw_data_fastqc_result
nohup fastqc -outdir 1_all_293T_IAA72h_ChIP_Seq_new_raw_data_fastqc_result --threads 72 *.fq.gz >1_all_293T_IAA72h_ChIP_Seq_new_raw_data_fastqc.runlog 2>&1 &
tail -f 1_all_293T_IAA72h_ChIP_Seq_new_raw_data_fastqc.runlog

# 并用 multiqc 进行结果汇总
# 提取原始数据reads数量
cd 1_all_293T_IAA72h_ChIP_Seq_new_raw_data_fastqc_result
multiqc ./*fastqc.zip --ignore *.html --no-megaqc-upload
cd ../

```

## clean data

```bash
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/
# 创建 new raw data 质量评估目录
mkdir 01_WN_293T_IAA72h_ChIP_Seq_new_clean_data_qc
# 先对新加测的数据的 clean data 进行质量评估
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/01_WN_293T_IAA72h_ChIP_Seq_new_clean_data_qc
# 根据样本信息表构造 new clean data 数据连接脚本，上传至该目录并运行
bash 0_make_WN_293T_IAA72h_ChIP_Seq_new_clean_data_fq_gz_file_links.sh
ls -1
#
mkdir 1_all_293T_IAA72h_ChIP_Seq_new_clean_data_fastqc_result
nohup fastqc -outdir 1_all_293T_IAA72h_ChIP_Seq_new_clean_data_fastqc_result --threads 72 *.fq.gz >1_all_293T_IAA72h_ChIP_Seq_new_clean_data_fastqc.runlog 2>&1 &
tail -f 1_all_293T_IAA72h_ChIP_Seq_new_clean_data_fastqc.runlog

# 并用 multiqc 进行结果汇总
# 提取原始数据reads数量
cd 1_all_293T_IAA72h_ChIP_Seq_new_clean_data_fastqc_result
multiqc ./*fastqc.zip --ignore *.html --no-megaqc-upload
cd ../

```

## merged 

```bash
# 构造两批数据 raw data 的文件连接，上传并运行
# 0_make_WN_293T_IAA72h_ChIP_Seq_raw_data_part1_fq_gz_file_links.sh （基于上次构造的脚本改造，F:\projects_on_940xa\2022081701_WN_ChIP_Seq_data\01_WN_293T_IAA72h_ChIP_Seq）
# 0_make_WN_293T_IAA72h_ChIP_Seq_raw_data_part2_fq_gz_file_links.sh
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data
bash 0_make_WN_293T_IAA72h_ChIP_Seq_raw_data_part1_fq_gz_file_links.sh
ls -1
bash 0_make_WN_293T_IAA72h_ChIP_Seq_raw_data_part2_fq_gz_file_links.sh
ls -1

# 编写合并两批数据的脚本并运行 (使用cat比 zcat|gzip快很多)
bash 01_WN_293T_IAA72h_ChIP_Seq_raw_data_part12_merged.sh
# 合并完成后删除连接文件
rm *raw.R[12].fq.gz

# 对合并的数据进行质量评估
mkdir 1_all_293T_IAA72h_ChIP_Seq_data_fastqc_result
nohup fastqc -outdir 1_all_293T_IAA72h_ChIP_Seq_data_fastqc_result --threads 72 *.fq.gz >1_all_293T_IAA72h_ChIP_Seq_data_fastqc.runlog 2>&1 &
tail -f 1_all_293T_IAA72h_ChIP_Seq_data_fastqc.runlog

# 并用 multiqc 进行结果汇总
# 提取原始数据reads数量
cd 1_all_293T_IAA72h_ChIP_Seq_data_fastqc_result
multiqc ./*fastqc.zip --ignore *.html --no-megaqc-upload
cd ../

```



# cutadapt and QC

```bash
# 根据 fastqc 的结果提示，在 illumina 的接头序列文件中找到如下序列，并对其中一个文件进行了grep统计，确认是接头
# Read 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
# Read 2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
# 使用如下命令生成去接头脚本并运行
for i in `ls -1 *R1.fq.gz`;do echo 'cutadapt -j 72 -q 20,20 --trim-n -m 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o '${i%.R1.fq.gz}.R1.cutadapt.fq.gz' -p '${i%.R1.fq.gz}.R2.cutadapt.fq.gz' '$i' '${i%.R1.fq.gz}.R2.fq.gz;done >2.1_all_293T_IAA72h_ChIP_Seq_raw_data_run_cutadapt.sh
nohup bash 2.1_all_293T_IAA72h_ChIP_Seq_raw_data_run_cutadapt.sh >2.1_all_293T_IAA72h_ChIP_Seq_raw_data_run_cutadapt.sh.runlog 2>&1 &
tail -f 2.1_all_293T_IAA72h_ChIP_Seq_raw_data_run_cutadapt.sh.runlog

# 去接头后再次进行质量评估
mkdir 2_all_293T_IAA72h_ChIP_Seq_cutadapt_fastqc_result
nohup fastqc -outdir 2_all_293T_IAA72h_ChIP_Seq_cutadapt_fastqc_result --threads 72 *.cutadapt.fq.gz >2.2_all_293T_IAA72h_ChIP_Seq_cutadapt_fastqc.runlog 2>&1 &
tail -f 2.2_all_293T_IAA72h_ChIP_Seq_cutadapt_fastqc.runlog

# 并用 multiqc 进行结果汇总
# 提取原始数据reads数量
cd 2_all_293T_IAA72h_ChIP_Seq_cutadapt_fastqc_result
multiqc ./*fastqc.zip --ignore *.html --no-megaqc-upload
cd ../

# 删除原始合并的 fq.gz 文件, 后续只使用去接头后的 fq.gz
rm *rep[123].R[12].fq.gz

```



# mapping

```bash
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data
mkdir 3_merge_clean_data_and_mapping_stat
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/3_merge_clean_data_and_mapping_stat
# 基于 raw data 软连接创建脚本，构造连接文件脚本，上传并运行
bash 0_make_WN_293T_IAA72h_ChIP_Seq_clean_data_part1_fq_gz_file_links.sh
ls -1 *part1*.fq.gz
bash 0_make_WN_293T_IAA72h_ChIP_Seq_clean_data_part2_fq_gz_file_links.sh
ls -1 *part2*.fq.gz
# 合并 part1和part2; 编写合并两批数据的脚本并运行 (使用cat比 zcat|gzip快很多)
bash 01_WN_293T_IAA72h_ChIP_Seq_clean_data_part12_merged.sh
ls -1
# 合并完成后删除连接文件
rm *clean.R[12].fq.gz

# 对合并的数据进行质量评估
mkdir 1_all_293T_IAA72h_ChIP_Seq_data_fastqc_result
nohup fastqc -outdir 1_all_293T_IAA72h_ChIP_Seq_data_fastqc_result --threads 72 *.fq.gz >1_all_293T_IAA72h_ChIP_Seq_data_fastqc.runlog 2>&1 &
tail -f 1_all_293T_IAA72h_ChIP_Seq_data_fastqc.runlog

# 并用 multiqc 进行结果汇总
# 提取原始数据reads数量
cd 1_all_293T_IAA72h_ChIP_Seq_data_fastqc_result
multiqc ./*fastqc.zip --ignore *.html --no-megaqc-upload
cd ../

# 比对汇总查看比对效率和reads数量
# 可使用如下命令构造比对脚本， 然后运行
for i in `ls -1 *rep[123].R1.fq.gz`;do echo 'bowtie2 -t -p 72 --no-unal -x /home/sunwenju/3_genomeData/human/hg19/bowtie2Index/hg19bt2idx -1 '$i' -2 '${i%.R1.fq.gz}.R2.fq.gz' -S '${i%.R1.fq.gz}.bowtie2.sam' --un-gz '${i%.R1.fq.gz}.notAlign.unpairedReads.fq.gz' --un-conc-gz '${i%.R1.fq.gz}.notAlign.pairedReads.fq.gz' >'${i%.R1.fq.gz}.bowtie2.runlog' 2>&1';done >3_all_293T_IAA72h_ChIP_Seq_bowtie2_map.sh
nohup bash 3_all_293T_IAA72h_ChIP_Seq_bowtie2_map.sh >3_all_293T_IAA72h_ChIP_Seq_bowtie2_map.sh.runlog 2>&1 &
# tail -f 3_all_293T_IAA72h_ChIP_Seq_bowtie2_map.sh.runlog
tail -f `ls -1 *bowtie2.runlog|tail -n 1`

# 使用 multiqc 对 bowtie2 比对情况进行汇总
multiqc --no-megaqc-upload -m bowtie2 -o 3_all_293T_IAA72h_ChIP_Seq_bowtie2_map_results_summary ./*.bowtie2.runlog

# 从结果来看，合并clean data 比对效率仅提高不到4%，可用的read 绝对数量反而变少了，所以就用 raw data 合并后去接头的数据比对的结果
# 删除无用结果
rm *.bowtie2.sam
rm *notAlign.pairedReads.fq.1.gz
rm *notAlign.pairedReads.fq.2.gz
rm *notAlign.unpairedReads.fq.gz
rm *.fq.gz
# 同时与之前的比对结果数量相比较

```



# convert sam to bam

```bash
# 将比对sam结果转换为 按座位排序的 bam 并建立索引
# for i in `ls -1 *.bowtie2.sam|sort`; do echo $i; samtools view --threads 36 -bS $i | samtools sort --threads 36 -o ${i%.sam}.sorted.bam; samtools index ${i%.sam}.sorted.bam; done
nohup bash 4_WN_IAA72h_all_ChIP_Seq_data_bowtie2_sam_to_sorted_bam.sh >4_WN_IAA72h_all_ChIP_Seq_data_bowtie2_sam_to_sorted_bam.sh.runlog 2>&1 &
tail -f 4_WN_IAA72h_all_ChIP_Seq_data_bowtie2_sam_to_sorted_bam.sh.runlog
# rm *.bowtie2.sam

# 使用samtools flagstat 统计比对情况
# 每个bam文件不到一分钟
for i in `ls -1 *bowtie2.sorted.bam`;do echo $i; samtools flagstat --threads 16 $i;done >3_293T_IAA72h_ChIP_Seq_bowtie2_map_samtools_flagstat.txt

# 统计比对到各个染色体上的 总read数量 
for i in `ls -1 *sorted.bam|sort`;do echo $i; samtools view $i |grep -v -e '^@' |cut -f 3 |sort |uniq -c |awk 'BEGIN{OFS="\t"}{if($0 !~ /_/){print $2,$1}}' |sort -k1,1 >${i%.sam}.chr_mapped_reads_stat.txt;done
# 统计比对到各个染色体上的 uniq read数量 
for i in `ls -1 *sorted.bam|sort`;do echo $i; samtools view $i |grep -v -e '^@' |cut -f 3,4 |sort |uniq |cut -f 1 |sort |uniq -c |awk 'BEGIN{OFS="\t"}{if($0 !~ /_/){print $2,$1}}' |sort -k1,1 >${i%.sam}.chr_mapped_uniq_reads_stat.txt;done

# 比对到各个染色体上的 read数量 统计结果合并
paste `ls -1 *chr_mapped_reads_stat.txt|sort` |head
# 检测各个统计结果的染色体数量和顺序是否一致
paste `ls -1 *chr_mapped_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"}{if($1!=$3 || $1!=$5 || $1!=$7 || $1!=$9 || $1!=$11 || $1!=$13 || $1!=$15 || $1!=$17 || $1!=$19 || $1!=$21 || $1!=$23 || $1!=$25 || $1!=$27 || $1!=$29 || $1!=$31 || $1!=$33 || $1!=$35 || $1!=$37 || $1!=$39 || $1!=$41 || $1!=$43 || $1!=$45 || $1!=$47 || $1!=$49 || $1!=$51 || $1!=$53 || $1!=$55 || $1!=$57 || $1!=$59 ){print $1,$2}}' |head
# 染色体数量和顺序一致，则合并数据列
paste `ls -1 *chr_mapped_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"; print "chrom","IAA_minus_INPUT_rep1","IAA_minus_INPUT_rep2","IAA_minus_INPUT_rep3","IAA_minus_IP_CTCF_rep1","IAA_minus_IP_CTCF_rep2","IAA_minus_IP_CTCF_rep3","IAA_minus_IP_MBD3_rep1","IAA_minus_IP_MBD3_rep2","IAA_minus_IP_MBD3_rep3","IAA_minus_IP_TET1_rep1","IAA_minus_IP_TET1_rep2","IAA_minus_IP_TET1_rep3","IAA_minus_IP_TET2_rep1","IAA_minus_IP_TET2_rep2","IAA_minus_IP_TET2_rep3","IAA_plus_INPUT_rep1","IAA_plus_INPUT_rep2","IAA_plus_INPUT_rep3","IAA_plus_IP_CTCF_rep1","IAA_plus_IP_CTCF_rep2","IAA_plus_IP_CTCF_rep3","IAA_plus_IP_MBD3_rep1","IAA_plus_IP_MBD3_rep2","IAA_plus_IP_MBD3_rep3","IAA_plus_IP_TET1_rep1","IAA_plus_IP_TET1_rep2","IAA_plus_IP_TET1_rep3","IAA_plus_IP_TET2_rep1","IAA_plus_IP_TET2_rep2","IAA_plus_IP_TET2_rep3"}{if($1==$3 && $1==$5 && $1==$7 && $1==$9 && $1==$11 && $1==$13 && $1==$15 && $1==$17 && $1==$19 && $1==$21 && $1==$23 && $1==$25 && $1==$27 && $1==$29 && $1==$31 && $1==$33 && $1==$35 && $1==$37 && $1==$39 && $1==$41 && $1==$43 && $1==$45 && $1==$47 && $1==$49 && $1==$51 && $1==$53 && $1==$55 && $1==$57 && $1==$59){print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36,$38,$40,$42,$44,$46,$48,$50,$52,$54,$56,$58,$60}}' |head
paste `ls -1 *chr_mapped_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"; print "chrom","IAA_minus_INPUT_rep1","IAA_minus_INPUT_rep2","IAA_minus_INPUT_rep3","IAA_minus_IP_CTCF_rep1","IAA_minus_IP_CTCF_rep2","IAA_minus_IP_CTCF_rep3","IAA_minus_IP_MBD3_rep1","IAA_minus_IP_MBD3_rep2","IAA_minus_IP_MBD3_rep3","IAA_minus_IP_TET1_rep1","IAA_minus_IP_TET1_rep2","IAA_minus_IP_TET1_rep3","IAA_minus_IP_TET2_rep1","IAA_minus_IP_TET2_rep2","IAA_minus_IP_TET2_rep3","IAA_plus_INPUT_rep1","IAA_plus_INPUT_rep2","IAA_plus_INPUT_rep3","IAA_plus_IP_CTCF_rep1","IAA_plus_IP_CTCF_rep2","IAA_plus_IP_CTCF_rep3","IAA_plus_IP_MBD3_rep1","IAA_plus_IP_MBD3_rep2","IAA_plus_IP_MBD3_rep3","IAA_plus_IP_TET1_rep1","IAA_plus_IP_TET1_rep2","IAA_plus_IP_TET1_rep3","IAA_plus_IP_TET2_rep1","IAA_plus_IP_TET2_rep2","IAA_plus_IP_TET2_rep3"}{if($1==$3 && $1==$5 && $1==$7 && $1==$9 && $1==$11 && $1==$13 && $1==$15 && $1==$17 && $1==$19 && $1==$21 && $1==$23 && $1==$25 && $1==$27 && $1==$29 && $1==$31 && $1==$33 && $1==$35 && $1==$37 && $1==$39 && $1==$41 && $1==$43 && $1==$45 && $1==$47 && $1==$49 && $1==$51 && $1==$53 && $1==$55 && $1==$57 && $1==$59){print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36,$38,$40,$42,$44,$46,$48,$50,$52,$54,$56,$58,$60}}' >all_samples_chr_mapped_reads_stat.summary.txt

# 比对到各个染色体上的 uniq read数量 统计结果合并
paste `ls -1 *chr_mapped_uniq_reads_stat.txt|sort` |head
# 检测各个统计结果的染色体数量和顺序是否一致
paste `ls -1 *chr_mapped_uniq_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"}{if($1!=$3 || $1!=$5 || $1!=$7 || $1!=$9 || $1!=$11 || $1!=$13 || $1!=$15 || $1!=$17 || $1!=$19 || $1!=$21 || $1!=$23 || $1!=$25 || $1!=$27 || $1!=$29 || $1!=$31 || $1!=$33 || $1!=$35 || $1!=$37 || $1!=$39 || $1!=$41 || $1!=$43 || $1!=$45 || $1!=$47 || $1!=$49 || $1!=$51 || $1!=$53 || $1!=$55 || $1!=$57 || $1!=$59){print $1,$2}}' |head
# 染色体数量和顺序一致，则合并数据列
paste `ls -1 *chr_mapped_uniq_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"; print "chrom","IAA_minus_INPUT_rep1","IAA_minus_INPUT_rep2","IAA_minus_INPUT_rep3","IAA_minus_IP_CTCF_rep1","IAA_minus_IP_CTCF_rep2","IAA_minus_IP_CTCF_rep3","IAA_minus_IP_MBD3_rep1","IAA_minus_IP_MBD3_rep2","IAA_minus_IP_MBD3_rep3","IAA_minus_IP_TET1_rep1","IAA_minus_IP_TET1_rep2","IAA_minus_IP_TET1_rep3","IAA_minus_IP_TET2_rep1","IAA_minus_IP_TET2_rep2","IAA_minus_IP_TET2_rep3","IAA_plus_INPUT_rep1","IAA_plus_INPUT_rep2","IAA_plus_INPUT_rep3","IAA_plus_IP_CTCF_rep1","IAA_plus_IP_CTCF_rep2","IAA_plus_IP_CTCF_rep3","IAA_plus_IP_MBD3_rep1","IAA_plus_IP_MBD3_rep2","IAA_plus_IP_MBD3_rep3","IAA_plus_IP_TET1_rep1","IAA_plus_IP_TET1_rep2","IAA_plus_IP_TET1_rep3","IAA_plus_IP_TET2_rep1","IAA_plus_IP_TET2_rep2","IAA_plus_IP_TET2_rep3"}{if($1==$3 && $1==$5 && $1==$7 && $1==$9 && $1==$11 && $1==$13 && $1==$15 && $1==$17 && $1==$19 && $1==$21 && $1==$23 && $1==$25 && $1==$27 && $1==$29 && $1==$31 && $1==$33 && $1==$35 && $1==$37 && $1==$39 && $1==$41 && $1==$43 && $1==$45 && $1==$47 && $1==$49 && $1==$51 && $1==$53 && $1==$55 && $1==$57 && $1==$59){print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36,$38,$40,$42,$44,$46,$48,$50,$52,$54,$56,$58,$60}}' |head
paste `ls -1 *chr_mapped_uniq_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"; print "chrom","IAA_minus_INPUT_rep1","IAA_minus_INPUT_rep2","IAA_minus_INPUT_rep3","IAA_minus_IP_CTCF_rep1","IAA_minus_IP_CTCF_rep2","IAA_minus_IP_CTCF_rep3","IAA_minus_IP_MBD3_rep1","IAA_minus_IP_MBD3_rep2","IAA_minus_IP_MBD3_rep3","IAA_minus_IP_TET1_rep1","IAA_minus_IP_TET1_rep2","IAA_minus_IP_TET1_rep3","IAA_minus_IP_TET2_rep1","IAA_minus_IP_TET2_rep2","IAA_minus_IP_TET2_rep3","IAA_plus_INPUT_rep1","IAA_plus_INPUT_rep2","IAA_plus_INPUT_rep3","IAA_plus_IP_CTCF_rep1","IAA_plus_IP_CTCF_rep2","IAA_plus_IP_CTCF_rep3","IAA_plus_IP_MBD3_rep1","IAA_plus_IP_MBD3_rep2","IAA_plus_IP_MBD3_rep3","IAA_plus_IP_TET1_rep1","IAA_plus_IP_TET1_rep2","IAA_plus_IP_TET1_rep3","IAA_plus_IP_TET2_rep1","IAA_plus_IP_TET2_rep2","IAA_plus_IP_TET2_rep3"}{if($1==$3 && $1==$5 && $1==$7 && $1==$9 && $1==$11 && $1==$13 && $1==$15 && $1==$17 && $1==$19 && $1==$21 && $1==$23 && $1==$25 && $1==$27 && $1==$29 && $1==$31 && $1==$33 && $1==$35 && $1==$37 && $1==$39 && $1==$41 && $1==$43 && $1==$45 && $1==$47 && $1==$49 && $1==$51 && $1==$53 && $1==$55 && $1==$57 && $1==$59){print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36,$38,$40,$42,$44,$46,$48,$50,$52,$54,$56,$58,$60}}' >all_samples_chr_mapped_uniq_reads_stat.summary.txt

```

```bash
# 绘制reads 分布热图时要用到 rep123合并的bam文件
# 创建 bam merge 脚本 5_make_all_293T_IAA72h_ChIP_Seq_rep123_merged_bam_and_index.sh 内容如下：
# # 合并各重复的 bam
# samtools merge --threads 72 HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_IAA_minus_IP_CTCF_rep123.merged.sorted.bam HEK293T_IAA_minus_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_CTCF_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_IAA_minus_IP_MBD3_rep123.merged.sorted.bam HEK293T_IAA_minus_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_MBD3_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_IAA_minus_IP_TET1_rep123.merged.sorted.bam HEK293T_IAA_minus_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET1_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_IAA_minus_IP_TET2_rep123.merged.sorted.bam HEK293T_IAA_minus_IP_TET2_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET2_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET2_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_IAA_plus_IP_CTCF_rep123.merged.sorted.bam HEK293T_IAA_plus_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_CTCF_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_IAA_plus_IP_MBD3_rep123.merged.sorted.bam HEK293T_IAA_plus_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_MBD3_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_IAA_plus_IP_TET1_rep123.merged.sorted.bam HEK293T_IAA_plus_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET1_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_IAA_plus_IP_TET2_rep123.merged.sorted.bam HEK293T_IAA_plus_IP_TET2_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET2_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET2_rep3.bowtie2.sorted.bam
# # 为所有 merged bam 建立索引
# for i in `ls -1 *_rep123.merged.sorted.bam`; do echo $i; samtools index -@ 72 $i;done
nohup bash 5_make_all_293T_IAA72h_ChIP_Seq_rep123_merged_bam_and_index.sh >5_make_all_293T_IAA72h_ChIP_Seq_rep123_merged_bam_and_index.sh.runlog 2>&1 &
tail -f 5_make_all_293T_IAA72h_ChIP_Seq_rep123_merged_bam_and_index.sh.runlog

```



# Detection of enrichment and repeatability

```bash
# 参考： https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html
# 可以判断实验的成败、富集程度的好坏等
# IAA_minus samples
plotFingerprint -b HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam HEK293T_IAA_minus_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_IAA_minus_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_MBD3_rep3.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET1_rep3.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET2_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET2_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET2_rep3.bowtie2.sorted.bam --labels IAA_minus_INPUT_rep1 IAA_minus_INPUT_rep2 IAA_minus_INPUT_rep3 IAA_minus_IP_CTCF_rep1 IAA_minus_IP_CTCF_rep2 IAA_minus_IP_CTCF_rep3 IAA_minus_IP_MBD3_rep1 IAA_minus_IP_MBD3_rep2 IAA_minus_IP_MBD3_rep3 IAA_minus_IP_TET1_rep1 IAA_minus_IP_TET1_rep2 IAA_minus_IP_TET1_rep3 IAA_minus_IP_TET2_rep1 IAA_minus_IP_TET2_rep2 IAA_minus_IP_TET2_rep3 -o plotFingerprint_for_IAA_minus_samples.pdf --extendReads --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max"
# IAA_plus samples
plotFingerprint -b HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_MBD3_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET1_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET2_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET2_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET2_rep3.bowtie2.sorted.bam --labels IAA_plus_INPUT_rep1 IAA_plus_INPUT_rep2 IAA_plus_INPUT_rep3 IAA_plus_IP_CTCF_rep1 IAA_plus_IP_CTCF_rep2 IAA_plus_IP_CTCF_rep3 IAA_plus_IP_MBD3_rep1 IAA_plus_IP_MBD3_rep2 IAA_plus_IP_MBD3_rep3 IAA_plus_IP_TET1_rep1 IAA_plus_IP_TET1_rep2 IAA_plus_IP_TET1_rep3 IAA_plus_IP_TET2_rep1 IAA_plus_IP_TET2_rep2 IAA_plus_IP_TET2_rep3 -o plotFingerprint_for_IAA_plus_samples.pdf --extendReads --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max"

# CTCF
plotFingerprint -b HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam HEK293T_IAA_minus_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_CTCF_rep3.bowtie2.sorted.bam --labels IAA_minus_INPUT_rep1 IAA_minus_INPUT_rep2 IAA_minus_INPUT_rep3 IAA_minus_IP_CTCF_rep1 IAA_minus_IP_CTCF_rep2 IAA_minus_IP_CTCF_rep3 IAA_plus_INPUT_rep1 IAA_plus_INPUT_rep2 IAA_plus_INPUT_rep3 IAA_plus_IP_CTCF_rep1 IAA_plus_IP_CTCF_rep2 IAA_plus_IP_CTCF_rep3 -o plotFingerprint_for_IAA72h_CTCF_ChIP_Seq_samples.pdf --extendReads --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max"
# MBD3
plotFingerprint -b HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam HEK293T_IAA_minus_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_MBD3_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_MBD3_rep3.bowtie2.sorted.bam --labels IAA_minus_INPUT_rep1 IAA_minus_INPUT_rep2 IAA_minus_INPUT_rep3 IAA_minus_IP_MBD3_rep1 IAA_minus_IP_MBD3_rep2 IAA_minus_IP_MBD3_rep3 IAA_plus_INPUT_rep1 IAA_plus_INPUT_rep2 IAA_plus_INPUT_rep3 IAA_plus_IP_MBD3_rep1 IAA_plus_IP_MBD3_rep2 IAA_plus_IP_MBD3_rep3 -o plotFingerprint_for_IAA72h_MBD3_ChIP_Seq_samples.pdf --extendReads --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max"
# TET1
plotFingerprint -b HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET1_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET1_rep3.bowtie2.sorted.bam --labels IAA_minus_INPUT_rep1 IAA_minus_INPUT_rep2 IAA_minus_INPUT_rep3 IAA_minus_IP_TET1_rep1 IAA_minus_IP_TET1_rep2 IAA_minus_IP_TET1_rep3 IAA_plus_INPUT_rep1 IAA_plus_INPUT_rep2 IAA_plus_INPUT_rep3 IAA_plus_IP_TET1_rep1 IAA_plus_IP_TET1_rep2 IAA_plus_IP_TET1_rep3 -o plotFingerprint_for_IAA72h_TET1_ChIP_Seq_samples.pdf --extendReads --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max"
# TET2
plotFingerprint -b HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET2_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET2_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET2_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET2_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET2_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET2_rep3.bowtie2.sorted.bam --labels IAA_minus_INPUT_rep1 IAA_minus_INPUT_rep2 IAA_minus_INPUT_rep3 IAA_minus_IP_TET2_rep1 IAA_minus_IP_TET2_rep2 IAA_minus_IP_TET2_rep3 IAA_plus_INPUT_rep1 IAA_plus_INPUT_rep2 IAA_plus_INPUT_rep3 IAA_plus_IP_TET2_rep1 IAA_plus_IP_TET2_rep2 IAA_plus_IP_TET2_rep3 -o plotFingerprint_for_IAA72h_TET2_ChIP_Seq_samples.pdf --extendReads --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max"

# 使用 deeptools 分析 样本相似度
# 先对全基因组分bin统计 uniq reads数量
multiBamSummary bins --bamfiles HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam HEK293T_IAA_minus_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_IAA_minus_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_MBD3_rep3.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET1_rep3.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET2_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET2_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET2_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_MBD3_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET1_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET2_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET2_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET2_rep3.bowtie2.sorted.bam -out multiBamSummary_results_for_all_bam_files.npz --labels IAA_minus_INPUT_rep1 IAA_minus_INPUT_rep2 IAA_minus_INPUT_rep3 IAA_minus_IP_CTCF_rep1 IAA_minus_IP_CTCF_rep2 IAA_minus_IP_CTCF_rep3 IAA_minus_IP_MBD3_rep1 IAA_minus_IP_MBD3_rep2 IAA_minus_IP_MBD3_rep3 IAA_minus_IP_TET1_rep1 IAA_minus_IP_TET1_rep2 IAA_minus_IP_TET1_rep3 IAA_minus_IP_TET2_rep1 IAA_minus_IP_TET2_rep2 IAA_minus_IP_TET2_rep3 IAA_plus_INPUT_rep1 IAA_plus_INPUT_rep2 IAA_plus_INPUT_rep3 IAA_plus_IP_CTCF_rep1 IAA_plus_IP_CTCF_rep2 IAA_plus_IP_CTCF_rep3 IAA_plus_IP_MBD3_rep1 IAA_plus_IP_MBD3_rep2 IAA_plus_IP_MBD3_rep3 IAA_plus_IP_TET1_rep1 IAA_plus_IP_TET1_rep2 IAA_plus_IP_TET1_rep3 IAA_plus_IP_TET2_rep1 IAA_plus_IP_TET2_rep2 IAA_plus_IP_TET2_rep3 --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# 用 plotCorrelation 和 plotPCA 分析样本的相似度
plotCorrelation --corData multiBamSummary_results_for_all_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_all_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_all_samples.txt --colorMap coolwarm --plotNumbers
plotPCA --corData multiBamSummary_results_for_all_bam_files.npz --plotFile plotPCA_plot_for_all_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_all_samples.txt
plotPCA --rowCenter --corData multiBamSummary_results_for_all_bam_files.npz --plotFile plotPCA_plot_for_all_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_all_samples.rowCenter.txt

# IAA72h CTCF ChIP_Seq 样本的相似度
multiBamSummary bins --bamfiles HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam HEK293T_IAA_minus_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_CTCF_rep3.bowtie2.sorted.bam -out multiBamSummary_results_for_IAA72h_CTCF_ChIP_Seq_bam_files.npz --labels IAA_minus_INPUT_rep1 IAA_minus_INPUT_rep2 IAA_minus_INPUT_rep3 IAA_minus_IP_CTCF_rep1 IAA_minus_IP_CTCF_rep2 IAA_minus_IP_CTCF_rep3 IAA_plus_INPUT_rep1 IAA_plus_INPUT_rep2 IAA_plus_INPUT_rep3 IAA_plus_IP_CTCF_rep1 IAA_plus_IP_CTCF_rep2 IAA_plus_IP_CTCF_rep3 --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
plotCorrelation --corData multiBamSummary_results_for_IAA72h_CTCF_ChIP_Seq_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_IAA72h_CTCF_ChIP_Seq_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_IAA72h_CTCF_ChIP_Seq_samples.txt --colorMap coolwarm --plotNumbers
plotPCA --corData multiBamSummary_results_for_IAA72h_CTCF_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_IAA72h_CTCF_ChIP_Seq_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_IAA72h_CTCF_ChIP_Seq_samples.txt
plotPCA --rowCenter --corData multiBamSummary_results_for_IAA72h_CTCF_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_IAA72h_CTCF_ChIP_Seq_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_IAA72h_CTCF_ChIP_Seq_samples.rowCenter.txt
# IAA72h MBD3 ChIP_Seq 样本的相似度
multiBamSummary bins --bamfiles HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam HEK293T_IAA_minus_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_MBD3_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_MBD3_rep3.bowtie2.sorted.bam -out multiBamSummary_results_for_IAA72h_MBD3_ChIP_Seq_bam_files.npz --labels IAA_minus_INPUT_rep1 IAA_minus_INPUT_rep2 IAA_minus_INPUT_rep3 IAA_minus_IP_MBD3_rep1 IAA_minus_IP_MBD3_rep2 IAA_minus_IP_MBD3_rep3 IAA_plus_INPUT_rep1 IAA_plus_INPUT_rep2 IAA_plus_INPUT_rep3 IAA_plus_IP_MBD3_rep1 IAA_plus_IP_MBD3_rep2 IAA_plus_IP_MBD3_rep3 --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
plotCorrelation --corData multiBamSummary_results_for_IAA72h_MBD3_ChIP_Seq_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_IAA72h_MBD3_ChIP_Seq_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_IAA72h_MBD3_ChIP_Seq_samples.txt --colorMap coolwarm --plotNumbers
plotPCA --corData multiBamSummary_results_for_IAA72h_MBD3_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_IAA72h_MBD3_ChIP_Seq_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_IAA72h_MBD3_ChIP_Seq_samples.txt
plotPCA --rowCenter --corData multiBamSummary_results_for_IAA72h_MBD3_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_IAA72h_MBD3_ChIP_Seq_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_IAA72h_MBD3_ChIP_Seq_samples.rowCenter.txt
# IAA72h TET1 ChIP_Seq 样本的相似度
multiBamSummary bins --bamfiles HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET1_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET1_rep3.bowtie2.sorted.bam -out multiBamSummary_results_for_IAA72h_TET1_ChIP_Seq_bam_files.npz --labels IAA_minus_INPUT_rep1 IAA_minus_INPUT_rep2 IAA_minus_INPUT_rep3 IAA_minus_IP_TET1_rep1 IAA_minus_IP_TET1_rep2 IAA_minus_IP_TET1_rep3 IAA_plus_INPUT_rep1 IAA_plus_INPUT_rep2 IAA_plus_INPUT_rep3 IAA_plus_IP_TET1_rep1 IAA_plus_IP_TET1_rep2 IAA_plus_IP_TET1_rep3 --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
plotCorrelation --corData multiBamSummary_results_for_IAA72h_TET1_ChIP_Seq_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_IAA72h_TET1_ChIP_Seq_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_IAA72h_TET1_ChIP_Seq_samples.txt --colorMap coolwarm --plotNumbers
plotPCA --corData multiBamSummary_results_for_IAA72h_TET1_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_IAA72h_TET1_ChIP_Seq_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_IAA72h_TET1_ChIP_Seq_samples.txt
plotPCA --rowCenter --corData multiBamSummary_results_for_IAA72h_TET1_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_IAA72h_TET1_ChIP_Seq_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_IAA72h_TET1_ChIP_Seq_samples.rowCenter.txt
# IAA72h TET2 ChIP_Seq 样本的相似度
multiBamSummary bins --bamfiles HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET2_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET2_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET2_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET2_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET2_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET2_rep3.bowtie2.sorted.bam -out multiBamSummary_results_for_IAA72h_TET2_ChIP_Seq_bam_files.npz --labels IAA_minus_INPUT_rep1 IAA_minus_INPUT_rep2 IAA_minus_INPUT_rep3 IAA_minus_IP_TET2_rep1 IAA_minus_IP_TET2_rep2 IAA_minus_IP_TET2_rep3 IAA_plus_INPUT_rep1 IAA_plus_INPUT_rep2 IAA_plus_INPUT_rep3 IAA_plus_IP_TET2_rep1 IAA_plus_IP_TET2_rep2 IAA_plus_IP_TET2_rep3 --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
plotCorrelation --corData multiBamSummary_results_for_IAA72h_TET2_ChIP_Seq_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_IAA72h_TET2_ChIP_Seq_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_IAA72h_TET2_ChIP_Seq_samples.txt --colorMap coolwarm --plotNumbers
plotPCA --corData multiBamSummary_results_for_IAA72h_TET2_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_IAA72h_TET2_ChIP_Seq_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_IAA72h_TET2_ChIP_Seq_samples.txt
plotPCA --rowCenter --corData multiBamSummary_results_for_IAA72h_TET2_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_IAA72h_TET2_ChIP_Seq_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_IAA72h_TET2_ChIP_Seq_samples.rowCenter.txt

```



# macs2 peak calling

## peak calling for each repetition separately 

### CTCF ChIP-Seq

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data
mkdir macs2_peak_calling_1_for_CTCF_each_repeat
cd macs2_peak_calling_1_for_CTCF_each_repeat
# 上传 bam 连接构造脚本并运行
bash 6_make_IAA72h_CTCF_ChIP_Seq_sorted_bam_file_links.sh
ls -1 *.bam

# IAA_minus
nohup macs2 callpeak -t HEK293T_IAA_minus_IP_CTCF_rep1.bowtie2.sorted.bam -c HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_minus_IP_CTCF_rep1_TvsC -n IAA_minus_IP_CTCF_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_minus_IP_CTCF_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_IAA_minus_IP_CTCF_rep2.bowtie2.sorted.bam -c HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_minus_IP_CTCF_rep2_TvsC -n IAA_minus_IP_CTCF_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_minus_IP_CTCF_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_IAA_minus_IP_CTCF_rep3.bowtie2.sorted.bam -c HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_minus_IP_CTCF_rep3_TvsC -n IAA_minus_IP_CTCF_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_minus_IP_CTCF_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_IAA_minus_IP_CTCF_rep1_TvsC.runlog
tail -f macs2_callpeak_for_IAA_minus_IP_CTCF_rep2_TvsC.runlog
tail -f macs2_callpeak_for_IAA_minus_IP_CTCF_rep3_TvsC.runlog

# IAA_plus
nohup macs2 callpeak -t HEK293T_IAA_plus_IP_CTCF_rep1.bowtie2.sorted.bam -c HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_plus_IP_CTCF_rep1_TvsC -n IAA_plus_IP_CTCF_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_plus_IP_CTCF_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t HEK293T_IAA_plus_IP_CTCF_rep2.bowtie2.sorted.bam -c HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_plus_IP_CTCF_rep2_TvsC -n IAA_plus_IP_CTCF_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_plus_IP_CTCF_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t HEK293T_IAA_plus_IP_CTCF_rep3.bowtie2.sorted.bam -c HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_plus_IP_CTCF_rep3_TvsC -n IAA_plus_IP_CTCF_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_plus_IP_CTCF_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_IAA_plus_IP_CTCF_rep1_TvsC.runlog
tail -f macs2_callpeak_for_IAA_plus_IP_CTCF_rep2_TvsC.runlog
tail -f macs2_callpeak_for_IAA_plus_IP_CTCF_rep3_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/macs2_peak_calling_1_for_CTCF_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./IAA_minus_IP_CTCF_rep1_TvsC/IAA_minus_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# 38480
# ./IAA_minus_IP_CTCF_rep2_TvsC/IAA_minus_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# 38564
# ./IAA_minus_IP_CTCF_rep3_TvsC/IAA_minus_IP_CTCF_rep3_TvsC_peaks.narrowPeak
# 39886
# ./IAA_plus_IP_CTCF_rep1_TvsC/IAA_plus_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# 41
# ./IAA_plus_IP_CTCF_rep2_TvsC/IAA_plus_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# 50
# ./IAA_plus_IP_CTCF_rep3_TvsC/IAA_plus_IP_CTCF_rep3_TvsC_peaks.narrowPeak
# 29

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./IAA_minus_IP_CTCF_rep1_TvsC/IAA_minus_IP_CTCF_rep1_TvsC_summits.bed
# 39111
# ./IAA_minus_IP_CTCF_rep2_TvsC/IAA_minus_IP_CTCF_rep2_TvsC_summits.bed
# 39179
# ./IAA_minus_IP_CTCF_rep3_TvsC/IAA_minus_IP_CTCF_rep3_TvsC_summits.bed
# 40561
# ./IAA_plus_IP_CTCF_rep1_TvsC/IAA_plus_IP_CTCF_rep1_TvsC_summits.bed
# 42
# ./IAA_plus_IP_CTCF_rep2_TvsC/IAA_plus_IP_CTCF_rep2_TvsC_summits.bed
# 54
# ./IAA_plus_IP_CTCF_rep3_TvsC/IAA_plus_IP_CTCF_rep3_TvsC_summits.bed
# 32

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_IAA72h_IP_CTCF_ChIP_Seq_samples_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bigwig|sort`;do echo $i;done

```



### MBD3 ChIP-Seq

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data
mkdir macs2_peak_calling_1_for_MBD3_each_repeat
cd macs2_peak_calling_1_for_MBD3_each_repeat
# 上传 bam 连接构造脚本并运行
bash 6_make_IAA72h_MBD3_ChIP_Seq_sorted_bam_file_links.sh
ls -1 *.bam

# IAA_minus
nohup macs2 callpeak -t HEK293T_IAA_minus_IP_MBD3_rep1.bowtie2.sorted.bam -c HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_minus_IP_MBD3_rep1_TvsC -n IAA_minus_IP_MBD3_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_minus_IP_MBD3_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_IAA_minus_IP_MBD3_rep2.bowtie2.sorted.bam -c HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_minus_IP_MBD3_rep2_TvsC -n IAA_minus_IP_MBD3_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_minus_IP_MBD3_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_IAA_minus_IP_MBD3_rep3.bowtie2.sorted.bam -c HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_minus_IP_MBD3_rep3_TvsC -n IAA_minus_IP_MBD3_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_minus_IP_MBD3_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_IAA_minus_IP_MBD3_rep1_TvsC.runlog
tail -f macs2_callpeak_for_IAA_minus_IP_MBD3_rep2_TvsC.runlog
tail -f macs2_callpeak_for_IAA_minus_IP_MBD3_rep3_TvsC.runlog

# IAA_plus
nohup macs2 callpeak -t HEK293T_IAA_plus_IP_MBD3_rep1.bowtie2.sorted.bam -c HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_plus_IP_MBD3_rep1_TvsC -n IAA_plus_IP_MBD3_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_plus_IP_MBD3_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t HEK293T_IAA_plus_IP_MBD3_rep2.bowtie2.sorted.bam -c HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_plus_IP_MBD3_rep2_TvsC -n IAA_plus_IP_MBD3_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_plus_IP_MBD3_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t HEK293T_IAA_plus_IP_MBD3_rep3.bowtie2.sorted.bam -c HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_plus_IP_MBD3_rep3_TvsC -n IAA_plus_IP_MBD3_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_plus_IP_MBD3_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_IAA_plus_IP_MBD3_rep1_TvsC.runlog
tail -f macs2_callpeak_for_IAA_plus_IP_MBD3_rep2_TvsC.runlog
tail -f macs2_callpeak_for_IAA_plus_IP_MBD3_rep3_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/macs2_peak_calling_1_for_MBD3_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./IAA_minus_IP_MBD3_rep1_TvsC/IAA_minus_IP_MBD3_rep1_TvsC_peaks.narrowPeak
# 3002
# ./IAA_minus_IP_MBD3_rep2_TvsC/IAA_minus_IP_MBD3_rep2_TvsC_peaks.narrowPeak
# 2821
# ./IAA_minus_IP_MBD3_rep3_TvsC/IAA_minus_IP_MBD3_rep3_TvsC_peaks.narrowPeak
# 2788
# ./IAA_plus_IP_MBD3_rep1_TvsC/IAA_plus_IP_MBD3_rep1_TvsC_peaks.narrowPeak
# 29
# ./IAA_plus_IP_MBD3_rep2_TvsC/IAA_plus_IP_MBD3_rep2_TvsC_peaks.narrowPeak
# 29
# ./IAA_plus_IP_MBD3_rep3_TvsC/IAA_plus_IP_MBD3_rep3_TvsC_peaks.narrowPeak
# 56

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./IAA_minus_IP_MBD3_rep1_TvsC/IAA_minus_IP_MBD3_rep1_TvsC_summits.bed
# 3048
# ./IAA_minus_IP_MBD3_rep2_TvsC/IAA_minus_IP_MBD3_rep2_TvsC_summits.bed
# 2857
# ./IAA_minus_IP_MBD3_rep3_TvsC/IAA_minus_IP_MBD3_rep3_TvsC_summits.bed
# 2837
# ./IAA_plus_IP_MBD3_rep1_TvsC/IAA_plus_IP_MBD3_rep1_TvsC_summits.bed
# 33
# ./IAA_plus_IP_MBD3_rep2_TvsC/IAA_plus_IP_MBD3_rep2_TvsC_summits.bed
# 34
# ./IAA_plus_IP_MBD3_rep3_TvsC/IAA_plus_IP_MBD3_rep3_TvsC_summits.bed
# 71

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_IAA72h_IP_MBD3_ChIP_Seq_samples_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/macs2_peak_calling_1_for_MBD3_each_repeat
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bigwig|sort`;do echo $i;done

```


### TET1 ChIP-Seq

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data
mkdir macs2_peak_calling_1_for_TET1_each_repeat
cd macs2_peak_calling_1_for_TET1_each_repeat
# 上传 bam 连接构造脚本并运行
bash 6_make_IAA72h_TET1_ChIP_Seq_sorted_bam_file_links.sh
ls -1 *.bam

# IAA_minus
nohup macs2 callpeak -t HEK293T_IAA_minus_IP_TET1_rep1.bowtie2.sorted.bam -c HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_minus_IP_TET1_rep1_TvsC -n IAA_minus_IP_TET1_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_minus_IP_TET1_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_IAA_minus_IP_TET1_rep2.bowtie2.sorted.bam -c HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_minus_IP_TET1_rep2_TvsC -n IAA_minus_IP_TET1_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_minus_IP_TET1_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_IAA_minus_IP_TET1_rep3.bowtie2.sorted.bam -c HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_minus_IP_TET1_rep3_TvsC -n IAA_minus_IP_TET1_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_minus_IP_TET1_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_IAA_minus_IP_TET1_rep1_TvsC.runlog
tail -f macs2_callpeak_for_IAA_minus_IP_TET1_rep2_TvsC.runlog
tail -f macs2_callpeak_for_IAA_minus_IP_TET1_rep3_TvsC.runlog

# IAA_plus
nohup macs2 callpeak -t HEK293T_IAA_plus_IP_TET1_rep1.bowtie2.sorted.bam -c HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_plus_IP_TET1_rep1_TvsC -n IAA_plus_IP_TET1_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_plus_IP_TET1_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t HEK293T_IAA_plus_IP_TET1_rep2.bowtie2.sorted.bam -c HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_plus_IP_TET1_rep2_TvsC -n IAA_plus_IP_TET1_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_plus_IP_TET1_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t HEK293T_IAA_plus_IP_TET1_rep3.bowtie2.sorted.bam -c HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_plus_IP_TET1_rep3_TvsC -n IAA_plus_IP_TET1_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_plus_IP_TET1_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_IAA_plus_IP_TET1_rep1_TvsC.runlog
tail -f macs2_callpeak_for_IAA_plus_IP_TET1_rep2_TvsC.runlog
tail -f macs2_callpeak_for_IAA_plus_IP_TET1_rep3_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/macs2_peak_calling_1_for_TET1_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./IAA_minus_IP_TET1_rep1_TvsC/IAA_minus_IP_TET1_rep1_TvsC_peaks.narrowPeak
# 11287
# ./IAA_minus_IP_TET1_rep2_TvsC/IAA_minus_IP_TET1_rep2_TvsC_peaks.narrowPeak
# 11584
# ./IAA_minus_IP_TET1_rep3_TvsC/IAA_minus_IP_TET1_rep3_TvsC_peaks.narrowPeak
# 11156
# ./IAA_plus_IP_TET1_rep1_TvsC/IAA_plus_IP_TET1_rep1_TvsC_peaks.narrowPeak
# 121
# ./IAA_plus_IP_TET1_rep2_TvsC/IAA_plus_IP_TET1_rep2_TvsC_peaks.narrowPeak
# 136
# ./IAA_plus_IP_TET1_rep3_TvsC/IAA_plus_IP_TET1_rep3_TvsC_peaks.narrowPeak
# 124

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./IAA_minus_IP_TET1_rep1_TvsC/IAA_minus_IP_TET1_rep1_TvsC_summits.bed
# 11417
# ./IAA_minus_IP_TET1_rep2_TvsC/IAA_minus_IP_TET1_rep2_TvsC_summits.bed
# 11712
# ./IAA_minus_IP_TET1_rep3_TvsC/IAA_minus_IP_TET1_rep3_TvsC_summits.bed
# 11272
# ./IAA_plus_IP_TET1_rep1_TvsC/IAA_plus_IP_TET1_rep1_TvsC_summits.bed
# 137
# ./IAA_plus_IP_TET1_rep2_TvsC/IAA_plus_IP_TET1_rep2_TvsC_summits.bed
# 160
# ./IAA_plus_IP_TET1_rep3_TvsC/IAA_plus_IP_TET1_rep3_TvsC_summits.bed
# 142

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_IAA72h_IP_TET1_ChIP_Seq_samples_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/macs2_peak_calling_1_for_TET1_each_repeat
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bigwig|sort`;do echo $i;done

```


### TET2 ChIP-Seq

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data
mkdir macs2_peak_calling_1_for_TET2_each_repeat
cd macs2_peak_calling_1_for_TET2_each_repeat
# 上传 bam 连接构造脚本并运行
bash 6_make_IAA72h_TET2_ChIP_Seq_sorted_bam_file_links.sh
ls -1 *.bam

# IAA_minus
nohup macs2 callpeak -t HEK293T_IAA_minus_IP_TET2_rep1.bowtie2.sorted.bam -c HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_minus_IP_TET2_rep1_TvsC -n IAA_minus_IP_TET2_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_minus_IP_TET2_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_IAA_minus_IP_TET2_rep2.bowtie2.sorted.bam -c HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_minus_IP_TET2_rep2_TvsC -n IAA_minus_IP_TET2_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_minus_IP_TET2_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_IAA_minus_IP_TET2_rep3.bowtie2.sorted.bam -c HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_minus_IP_TET2_rep3_TvsC -n IAA_minus_IP_TET2_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_minus_IP_TET2_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_IAA_minus_IP_TET2_rep1_TvsC.runlog
tail -f macs2_callpeak_for_IAA_minus_IP_TET2_rep2_TvsC.runlog
tail -f macs2_callpeak_for_IAA_minus_IP_TET2_rep3_TvsC.runlog

# IAA_plus
nohup macs2 callpeak -t HEK293T_IAA_plus_IP_TET2_rep1.bowtie2.sorted.bam -c HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_plus_IP_TET2_rep1_TvsC -n IAA_plus_IP_TET2_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_plus_IP_TET2_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t HEK293T_IAA_plus_IP_TET2_rep2.bowtie2.sorted.bam -c HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_plus_IP_TET2_rep2_TvsC -n IAA_plus_IP_TET2_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_plus_IP_TET2_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t HEK293T_IAA_plus_IP_TET2_rep3.bowtie2.sorted.bam -c HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_plus_IP_TET2_rep3_TvsC -n IAA_plus_IP_TET2_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_plus_IP_TET2_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_IAA_plus_IP_TET2_rep1_TvsC.runlog
tail -f macs2_callpeak_for_IAA_plus_IP_TET2_rep2_TvsC.runlog
tail -f macs2_callpeak_for_IAA_plus_IP_TET2_rep3_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/macs2_peak_calling_1_for_TET2_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./IAA_minus_IP_TET2_rep1_TvsC/IAA_minus_IP_TET2_rep1_TvsC_peaks.narrowPeak
# 473
# ./IAA_minus_IP_TET2_rep2_TvsC/IAA_minus_IP_TET2_rep2_TvsC_peaks.narrowPeak
# 651
# ./IAA_minus_IP_TET2_rep3_TvsC/IAA_minus_IP_TET2_rep3_TvsC_peaks.narrowPeak
# 708
# ./IAA_plus_IP_TET2_rep1_TvsC/IAA_plus_IP_TET2_rep1_TvsC_peaks.narrowPeak
# 1566
# ./IAA_plus_IP_TET2_rep2_TvsC/IAA_plus_IP_TET2_rep2_TvsC_peaks.narrowPeak
# 1729
# ./IAA_plus_IP_TET2_rep3_TvsC/IAA_plus_IP_TET2_rep3_TvsC_peaks.narrowPeak
# 1986

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./IAA_minus_IP_TET2_rep1_TvsC/IAA_minus_IP_TET2_rep1_TvsC_summits.bed
# 500
# ./IAA_minus_IP_TET2_rep2_TvsC/IAA_minus_IP_TET2_rep2_TvsC_summits.bed
# 682
# ./IAA_minus_IP_TET2_rep3_TvsC/IAA_minus_IP_TET2_rep3_TvsC_summits.bed
# 750
# ./IAA_plus_IP_TET2_rep1_TvsC/IAA_plus_IP_TET2_rep1_TvsC_summits.bed
# 1750
# ./IAA_plus_IP_TET2_rep2_TvsC/IAA_plus_IP_TET2_rep2_TvsC_summits.bed
# 1885
# ./IAA_plus_IP_TET2_rep3_TvsC/IAA_plus_IP_TET2_rep3_TvsC_summits.bed
# 2187

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_IAA72h_IP_TET2_ChIP_Seq_samples_macs2_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/macs2_peak_calling_1_for_TET2_each_repeat
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bigwig|sort`;do echo $i;done

```



## peak calling for rep123 merged

### CTCF ChIP-Seq

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data
mkdir macs2_peak_calling_2_for_CTCF_rep123_together
cd macs2_peak_calling_2_for_CTCF_rep123_together
# 上传 bam 连接构造脚本并运行
bash 6_make_IAA72h_CTCF_ChIP_Seq_sorted_bam_file_links.sh
ls -1 *.bam

# IAA_minus_IP_CTCF_rep123_TvsC
nohup macs2 callpeak -t HEK293T_IAA_minus_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_CTCF_rep3.bowtie2.sorted.bam -c HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_minus_IP_CTCF_rep123_TvsC -n IAA_minus_IP_CTCF_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_minus_IP_CTCF_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_IAA_minus_IP_CTCF_rep123_TvsC.runlog

# IAA_plus_IP_CTCF_rep123_TvsC
nohup macs2 callpeak -t HEK293T_IAA_plus_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_CTCF_rep3.bowtie2.sorted.bam -c HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_plus_IP_CTCF_rep123_TvsC -n IAA_plus_IP_CTCF_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_plus_IP_CTCF_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_IAA_plus_IP_CTCF_rep123_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/macs2_peak_calling_2_for_CTCF_rep123_together
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./IAA_minus_IP_CTCF_rep123_TvsC/IAA_minus_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# 53727
# ./IAA_plus_IP_CTCF_rep123_TvsC/IAA_plus_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# 384

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./IAA_minus_IP_CTCF_rep123_TvsC/IAA_minus_IP_CTCF_rep123_TvsC_summits.bed
# 55645
# ./IAA_plus_IP_CTCF_rep123_TvsC/IAA_plus_IP_CTCF_rep123_TvsC_summits.bed
# 400

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_IAA72h_CTCF_ChIP_Seq_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/macs2_peak_calling_2_for_CTCF_rep123_together
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

```

```bash
# 两组结果交集
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/macs2_peak_calling_2_for_CTCF_rep123_together
for i in `find . -name *narrowPeak`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
wc -l `find . -name *_rep123_TvsC_peaks.narrowPeak.bed`

# 两组 narrowPeak 交集数量
bedtools intersect -a ./IAA_minus_IP_CTCF_rep123_TvsC/IAA_minus_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -b ./IAA_plus_IP_CTCF_rep123_TvsC/IAA_plus_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 320
bedtools intersect -a ./IAA_plus_IP_CTCF_rep123_TvsC/IAA_plus_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -b ./IAA_minus_IP_CTCF_rep123_TvsC/IAA_minus_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 322

# 两组 summit extB50 交集数量
bedtools slop -b 50 -i ./IAA_plus_IP_CTCF_rep123_TvsC/IAA_plus_IP_CTCF_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./IAA_plus_IP_CTCF_rep123_TvsC/IAA_plus_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
bedtools slop -b 50 -i ./IAA_minus_IP_CTCF_rep123_TvsC/IAA_minus_IP_CTCF_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./IAA_minus_IP_CTCF_rep123_TvsC/IAA_minus_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
# 求交集
bedtools intersect -a ./IAA_minus_IP_CTCF_rep123_TvsC/IAA_minus_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_plus_IP_CTCF_rep123_TvsC/IAA_plus_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 317
bedtools intersect -a ./IAA_plus_IP_CTCF_rep123_TvsC/IAA_plus_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_minus_IP_CTCF_rep123_TvsC/IAA_minus_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 317

# 构造 三类 summit
# common summit
bedtools intersect -a ./IAA_plus_IP_CTCF_rep123_TvsC/IAA_plus_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_minus_IP_CTCF_rep123_TvsC/IAA_minus_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a ./IAA_plus_IP_CTCF_rep123_TvsC/IAA_plus_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_minus_IP_CTCF_rep123_TvsC/IAA_minus_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >IP_CTCF_rep123_TvsC_IAA_minus_and_plus.common.summit.bed
# IAA_minus only summit
bedtools intersect -a ./IAA_minus_IP_CTCF_rep123_TvsC/IAA_minus_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_plus_IP_CTCF_rep123_TvsC/IAA_plus_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./IAA_minus_IP_CTCF_rep123_TvsC/IAA_minus_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_plus_IP_CTCF_rep123_TvsC/IAA_plus_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >IP_CTCF_rep123_TvsC_IAA_minus.only.summit.bed
# IAA_plus only summit
bedtools intersect -a ./IAA_plus_IP_CTCF_rep123_TvsC/IAA_plus_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_minus_IP_CTCF_rep123_TvsC/IAA_minus_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./IAA_plus_IP_CTCF_rep123_TvsC/IAA_plus_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_minus_IP_CTCF_rep123_TvsC/IAA_minus_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >IP_CTCF_rep123_TvsC_IAA_plus.only.summit.bed

wc -l *summit.bed
# 55328 IP_CTCF_rep123_TvsC_IAA_minus.only.summit.bed
#   317 IP_CTCF_rep123_TvsC_IAA_minus_and_plus.common.summit.bed
#    83 IP_CTCF_rep123_TvsC_IAA_plus.only.summit.bed

```



### MBD3 ChIP-Seq

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data
mkdir macs2_peak_calling_2_for_MBD3_rep123_together
cd macs2_peak_calling_2_for_MBD3_rep123_together
# 上传 bam 连接构造脚本并运行
bash 6_make_IAA72h_MBD3_ChIP_Seq_sorted_bam_file_links.sh
ls -1 *.bam

# IAA_minus_IP_MBD3_rep123_TvsC
nohup macs2 callpeak -t HEK293T_IAA_minus_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_MBD3_rep3.bowtie2.sorted.bam -c HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_minus_IP_MBD3_rep123_TvsC -n IAA_minus_IP_MBD3_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_minus_IP_MBD3_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_IAA_minus_IP_MBD3_rep123_TvsC.runlog

# IAA_plus_IP_MBD3_rep123_TvsC
nohup macs2 callpeak -t HEK293T_IAA_plus_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_MBD3_rep3.bowtie2.sorted.bam -c HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_plus_IP_MBD3_rep123_TvsC -n IAA_plus_IP_MBD3_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_plus_IP_MBD3_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_IAA_plus_IP_MBD3_rep123_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/macs2_peak_calling_2_for_MBD3_rep123_together
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./IAA_minus_IP_MBD3_rep123_TvsC/IAA_minus_IP_MBD3_rep123_TvsC_peaks.narrowPeak
# 16279
# ./IAA_plus_IP_MBD3_rep123_TvsC/IAA_plus_IP_MBD3_rep123_TvsC_peaks.narrowPeak
# 194

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./IAA_minus_IP_MBD3_rep123_TvsC/IAA_minus_IP_MBD3_rep123_TvsC_summits.bed
# 16985
# ./IAA_plus_IP_MBD3_rep123_TvsC/IAA_plus_IP_MBD3_rep123_TvsC_summits.bed
# 224

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_IAA72h_MBD3_ChIP_Seq_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/macs2_peak_calling_2_for_MBD3_rep123_together
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

```
```bash
# 两组结果交集
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/macs2_peak_calling_2_for_MBD3_rep123_together
for i in `find . -name *narrowPeak`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
wc -l `find . -name *_rep123_TvsC_peaks.narrowPeak.bed`

# 两组 narrowPeak 交集数量
bedtools intersect -a ./IAA_minus_IP_MBD3_rep123_TvsC/IAA_minus_IP_MBD3_rep123_TvsC_peaks.narrowPeak.bed -b ./IAA_plus_IP_MBD3_rep123_TvsC/IAA_plus_IP_MBD3_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 47
bedtools intersect -a ./IAA_plus_IP_MBD3_rep123_TvsC/IAA_plus_IP_MBD3_rep123_TvsC_peaks.narrowPeak.bed -b ./IAA_minus_IP_MBD3_rep123_TvsC/IAA_minus_IP_MBD3_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 50

# 两组 summit extB50 交集数量
bedtools slop -b 50 -i ./IAA_plus_IP_MBD3_rep123_TvsC/IAA_plus_IP_MBD3_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./IAA_plus_IP_MBD3_rep123_TvsC/IAA_plus_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed
bedtools slop -b 50 -i ./IAA_minus_IP_MBD3_rep123_TvsC/IAA_minus_IP_MBD3_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./IAA_minus_IP_MBD3_rep123_TvsC/IAA_minus_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed
# 求交集
bedtools intersect -a ./IAA_minus_IP_MBD3_rep123_TvsC/IAA_minus_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_plus_IP_MBD3_rep123_TvsC/IAA_plus_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 56
bedtools intersect -a ./IAA_plus_IP_MBD3_rep123_TvsC/IAA_plus_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_minus_IP_MBD3_rep123_TvsC/IAA_minus_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 56

# 构造 三类 summit
# common summit
bedtools intersect -a ./IAA_plus_IP_MBD3_rep123_TvsC/IAA_plus_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_minus_IP_MBD3_rep123_TvsC/IAA_minus_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a ./IAA_plus_IP_MBD3_rep123_TvsC/IAA_plus_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_minus_IP_MBD3_rep123_TvsC/IAA_minus_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >IP_MBD3_rep123_TvsC_IAA_minus_and_plus.common.summit.bed
# IAA_minus only summit
bedtools intersect -a ./IAA_minus_IP_MBD3_rep123_TvsC/IAA_minus_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_plus_IP_MBD3_rep123_TvsC/IAA_plus_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./IAA_minus_IP_MBD3_rep123_TvsC/IAA_minus_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_plus_IP_MBD3_rep123_TvsC/IAA_plus_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >IP_MBD3_rep123_TvsC_IAA_minus.only.summit.bed
# IAA_plus only summit
bedtools intersect -a ./IAA_plus_IP_MBD3_rep123_TvsC/IAA_plus_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_minus_IP_MBD3_rep123_TvsC/IAA_minus_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./IAA_plus_IP_MBD3_rep123_TvsC/IAA_plus_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_minus_IP_MBD3_rep123_TvsC/IAA_minus_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >IP_MBD3_rep123_TvsC_IAA_plus.only.summit.bed

wc -l *summit.bed
# 16929 IP_MBD3_rep123_TvsC_IAA_minus.only.summit.bed
#    56 IP_MBD3_rep123_TvsC_IAA_minus_and_plus.common.summit.bed
#   168 IP_MBD3_rep123_TvsC_IAA_plus.only.summit.bed

```



### TET1 ChIP-Seq

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data
mkdir macs2_peak_calling_2_for_TET1_rep123_together
cd macs2_peak_calling_2_for_TET1_rep123_together
# 上传 bam 连接构造脚本并运行
bash 6_make_IAA72h_TET1_ChIP_Seq_sorted_bam_file_links.sh
ls -1 *.bam

# IAA_minus_IP_TET1_rep123_TvsC
nohup macs2 callpeak -t HEK293T_IAA_minus_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET1_rep3.bowtie2.sorted.bam -c HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_minus_IP_TET1_rep123_TvsC -n IAA_minus_IP_TET1_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_minus_IP_TET1_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_IAA_minus_IP_TET1_rep123_TvsC.runlog

# IAA_plus_IP_TET1_rep123_TvsC
nohup macs2 callpeak -t HEK293T_IAA_plus_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET1_rep3.bowtie2.sorted.bam -c HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_plus_IP_TET1_rep123_TvsC -n IAA_plus_IP_TET1_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_plus_IP_TET1_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_IAA_plus_IP_TET1_rep123_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/macs2_peak_calling_2_for_TET1_rep123_together
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./IAA_minus_IP_TET1_rep123_TvsC/IAA_minus_IP_TET1_rep123_TvsC_peaks.narrowPeak
# 25161
# ./IAA_plus_IP_TET1_rep123_TvsC/IAA_plus_IP_TET1_rep123_TvsC_peaks.narrowPeak
# 1241

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./IAA_minus_IP_TET1_rep123_TvsC/IAA_minus_IP_TET1_rep123_TvsC_summits.bed
# 25892
# ./IAA_plus_IP_TET1_rep123_TvsC/IAA_plus_IP_TET1_rep123_TvsC_summits.bed
# 1343

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_IAA72h_TET1_ChIP_Seq_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/macs2_peak_calling_2_for_TET1_rep123_together
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

```
```bash
# 两组结果交集
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/macs2_peak_calling_2_for_TET1_rep123_together
for i in `find . -name *narrowPeak`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
wc -l `find . -name *_rep123_TvsC_peaks.narrowPeak.bed`

# 两组 narrowPeak 交集数量
bedtools intersect -a ./IAA_minus_IP_TET1_rep123_TvsC/IAA_minus_IP_TET1_rep123_TvsC_peaks.narrowPeak.bed -b ./IAA_plus_IP_TET1_rep123_TvsC/IAA_plus_IP_TET1_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 158
bedtools intersect -a ./IAA_plus_IP_TET1_rep123_TvsC/IAA_plus_IP_TET1_rep123_TvsC_peaks.narrowPeak.bed -b ./IAA_minus_IP_TET1_rep123_TvsC/IAA_minus_IP_TET1_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 160

# 两组 summit extB50 交集数量
bedtools slop -b 50 -i ./IAA_plus_IP_TET1_rep123_TvsC/IAA_plus_IP_TET1_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./IAA_plus_IP_TET1_rep123_TvsC/IAA_plus_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed
bedtools slop -b 50 -i ./IAA_minus_IP_TET1_rep123_TvsC/IAA_minus_IP_TET1_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./IAA_minus_IP_TET1_rep123_TvsC/IAA_minus_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed
# 求交集
bedtools intersect -a ./IAA_minus_IP_TET1_rep123_TvsC/IAA_minus_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_plus_IP_TET1_rep123_TvsC/IAA_plus_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 138
bedtools intersect -a ./IAA_plus_IP_TET1_rep123_TvsC/IAA_plus_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_minus_IP_TET1_rep123_TvsC/IAA_minus_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 138

# 构造 三类 summit
# common summit
bedtools intersect -a ./IAA_plus_IP_TET1_rep123_TvsC/IAA_plus_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_minus_IP_TET1_rep123_TvsC/IAA_minus_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a ./IAA_plus_IP_TET1_rep123_TvsC/IAA_plus_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_minus_IP_TET1_rep123_TvsC/IAA_minus_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >IP_TET1_rep123_TvsC_IAA_minus_and_plus.common.summit.bed
# IAA_minus only summit
bedtools intersect -a ./IAA_minus_IP_TET1_rep123_TvsC/IAA_minus_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_plus_IP_TET1_rep123_TvsC/IAA_plus_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./IAA_minus_IP_TET1_rep123_TvsC/IAA_minus_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_plus_IP_TET1_rep123_TvsC/IAA_plus_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >IP_TET1_rep123_TvsC_IAA_minus.only.summit.bed
# IAA_plus only summit
bedtools intersect -a ./IAA_plus_IP_TET1_rep123_TvsC/IAA_plus_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_minus_IP_TET1_rep123_TvsC/IAA_minus_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./IAA_plus_IP_TET1_rep123_TvsC/IAA_plus_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_minus_IP_TET1_rep123_TvsC/IAA_minus_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >IP_TET1_rep123_TvsC_IAA_plus.only.summit.bed

wc -l *summit.bed
# 25754 IP_TET1_rep123_TvsC_IAA_minus.only.summit.bed
#   138 IP_TET1_rep123_TvsC_IAA_minus_and_plus.common.summit.bed
#  1205 IP_TET1_rep123_TvsC_IAA_plus.only.summit.bed

```



### TET2 ChIP-Seq

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data
mkdir macs2_peak_calling_2_for_TET2_rep123_together
cd macs2_peak_calling_2_for_TET2_rep123_together
# 上传 bam 连接构造脚本并运行
bash 6_make_IAA72h_TET2_ChIP_Seq_sorted_bam_file_links.sh
ls -1 *.bam

# IAA_minus_IP_TET2_rep123_TvsC
nohup macs2 callpeak -t HEK293T_IAA_minus_IP_TET2_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET2_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_IP_TET2_rep3.bowtie2.sorted.bam -c HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_minus_IP_TET2_rep123_TvsC -n IAA_minus_IP_TET2_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_minus_IP_TET2_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_IAA_minus_IP_TET2_rep123_TvsC.runlog

# IAA_plus_IP_TET2_rep123_TvsC
nohup macs2 callpeak -t HEK293T_IAA_plus_IP_TET2_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET2_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_IP_TET2_rep3.bowtie2.sorted.bam -c HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_plus_IP_TET2_rep123_TvsC -n IAA_plus_IP_TET2_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_plus_IP_TET2_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_IAA_plus_IP_TET2_rep123_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/macs2_peak_calling_2_for_TET2_rep123_together
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./IAA_minus_IP_TET2_rep123_TvsC/IAA_minus_IP_TET2_rep123_TvsC_peaks.narrowPeak
# 6209
# ./IAA_plus_IP_TET2_rep123_TvsC/IAA_plus_IP_TET2_rep123_TvsC_peaks.narrowPeak
# 7763

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./IAA_minus_IP_TET2_rep123_TvsC/IAA_minus_IP_TET2_rep123_TvsC_summits.bed
# 6732
# ./IAA_plus_IP_TET2_rep123_TvsC/IAA_plus_IP_TET2_rep123_TvsC_summits.bed
# 9202

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_IAA72h_TET2_ChIP_Seq_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/macs2_peak_calling_2_for_TET2_rep123_together
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

```

```bash
# 两组结果交集
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/macs2_peak_calling_2_for_TET2_rep123_together
for i in `find . -name *narrowPeak`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
wc -l `find . -name *_rep123_TvsC_peaks.narrowPeak.bed`

# 两组 narrowPeak 交集数量
bedtools intersect -a ./IAA_minus_IP_TET2_rep123_TvsC/IAA_minus_IP_TET2_rep123_TvsC_peaks.narrowPeak.bed -b ./IAA_plus_IP_TET2_rep123_TvsC/IAA_plus_IP_TET2_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 35
bedtools intersect -a ./IAA_plus_IP_TET2_rep123_TvsC/IAA_plus_IP_TET2_rep123_TvsC_peaks.narrowPeak.bed -b ./IAA_minus_IP_TET2_rep123_TvsC/IAA_minus_IP_TET2_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 35

# 两组 summit extB50 交集数量
bedtools slop -b 50 -i ./IAA_plus_IP_TET2_rep123_TvsC/IAA_plus_IP_TET2_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./IAA_plus_IP_TET2_rep123_TvsC/IAA_plus_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed
bedtools slop -b 50 -i ./IAA_minus_IP_TET2_rep123_TvsC/IAA_minus_IP_TET2_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./IAA_minus_IP_TET2_rep123_TvsC/IAA_minus_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed
# 求交集
bedtools intersect -a ./IAA_minus_IP_TET2_rep123_TvsC/IAA_minus_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_plus_IP_TET2_rep123_TvsC/IAA_plus_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 37
bedtools intersect -a ./IAA_plus_IP_TET2_rep123_TvsC/IAA_plus_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_minus_IP_TET2_rep123_TvsC/IAA_minus_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 37

# 构造 三类 summit
# common summit
bedtools intersect -a ./IAA_plus_IP_TET2_rep123_TvsC/IAA_plus_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_minus_IP_TET2_rep123_TvsC/IAA_minus_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a ./IAA_plus_IP_TET2_rep123_TvsC/IAA_plus_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_minus_IP_TET2_rep123_TvsC/IAA_minus_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >IP_TET2_rep123_TvsC_IAA_minus_and_plus.common.summit.bed
# IAA_minus only summit
bedtools intersect -a ./IAA_minus_IP_TET2_rep123_TvsC/IAA_minus_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_plus_IP_TET2_rep123_TvsC/IAA_plus_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./IAA_minus_IP_TET2_rep123_TvsC/IAA_minus_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_plus_IP_TET2_rep123_TvsC/IAA_plus_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >IP_TET2_rep123_TvsC_IAA_minus.only.summit.bed
# IAA_plus only summit
bedtools intersect -a ./IAA_plus_IP_TET2_rep123_TvsC/IAA_plus_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_minus_IP_TET2_rep123_TvsC/IAA_minus_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./IAA_plus_IP_TET2_rep123_TvsC/IAA_plus_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./IAA_minus_IP_TET2_rep123_TvsC/IAA_minus_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >IP_TET2_rep123_TvsC_IAA_plus.only.summit.bed

wc -l *summit.bed
# 6695 IP_TET2_rep123_TvsC_IAA_minus.only.summit.bed
#   37 IP_TET2_rep123_TvsC_IAA_minus_and_plus.common.summit.bed
# 9165 IP_TET2_rep123_TvsC_IAA_plus.only.summit.bed

```



#  summit f1k reads enrichment heatmap

## CTCF ChIP-Seq

### use individual bam of each repeat

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data
mkdir plotHeatMap_use_bamCompare_bigwig_for_CTCF
cd plotHeatMap_use_bamCompare_bigwig_for_CTCF
# 上传 6_make_IAA72h_CTCF_ChIP_Seq_sorted_bam_file_links.sh
# 创建 bam 和 bai 文件链接
bash 6_make_IAA72h_CTCF_ChIP_Seq_sorted_bam_file_links.sh
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_CTCF_rep123_together/IP_CTCF_rep123_TvsC_IAA_minus.only.summit.bed ./CTCF_IAA_minus.only.summit.bed
ln -s ../macs2_peak_calling_2_for_CTCF_rep123_together/IP_CTCF_rep123_TvsC_IAA_minus_and_plus.common.summit.bed ./CTCF_IAA_minus_and_plus.common.summit.bed
ln -s ../macs2_peak_calling_2_for_CTCF_rep123_together/IP_CTCF_rep123_TvsC_IAA_plus.only.summit.bed ./CTCF_IAA_plus.only.summit.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_CTCF_bdgdiff/RIB_IAA_plus_and_minus.CTCF.summit.sorted.bed ./RIB_IAA_plus_and_minus.CTCF.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_CTCF_bdgdiff/ROI_IAA_minus.CTCF.summit.sorted.bed ./ROI_IAA_minus.CTCF.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_CTCF_bdgdiff/ROI_IAA_plus.CTCF.summit.sorted.bed ./ROI_IAA_plus.CTCF.summit.sorted.bed
ls -1

```

#### creat log2 trans bigwig by using bamCompare 

```bash
# 构造 7_run_bam_comCompare_for_CTCF_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# bamCompare -b1 HEK293T_IAA_minus_IP_CTCF_rep1.bowtie2.sorted.bam -b2 HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam -o IAA_minus_IP_CTCF_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_IAA_minus_IP_CTCF_rep2.bowtie2.sorted.bam -b2 HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam -o IAA_minus_IP_CTCF_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_IAA_minus_IP_CTCF_rep3.bowtie2.sorted.bam -b2 HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam -o IAA_minus_IP_CTCF_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_IAA_plus_IP_CTCF_rep1.bowtie2.sorted.bam -b2 HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam -o IAA_plus_IP_CTCF_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_IAA_plus_IP_CTCF_rep2.bowtie2.sorted.bam -b2 HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam -o IAA_plus_IP_CTCF_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_IAA_plus_IP_CTCF_rep3.bowtie2.sorted.bam -b2 HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam -o IAA_plus_IP_CTCF_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
nohup bash 7_run_bam_comCompare_for_CTCF_each_repeat.sh >7_run_bam_comCompare_for_CTCF_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bam_comCompare_for_CTCF_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_CTCF_rep1_TvsC.bigwig IAA_minus_IP_CTCF_rep2_TvsC.bigwig IAA_minus_IP_CTCF_rep3_TvsC.bigwig IAA_plus_IP_CTCF_rep1_TvsC.bigwig IAA_plus_IP_CTCF_rep2_TvsC.bigwig IAA_plus_IP_CTCF_rep3_TvsC.bigwig -R CTCF_IAA_minus.only.summit.bed CTCF_IAA_minus_and_plus.common.summit.bed CTCF_IAA_plus.only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_CTCF_ChIP_Seq_all_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_CTCF_ChIP_Seq_all_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_IAA72h_CTCF_ChIP_Seq_all_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_IAA_minus_only(55328)" "common_summits(317)" "CTCF_IAA_plus_only(83)" --samplesLabel "IAA_minus_CTCF_rep1_TvsC" "IAA_minus_CTCF_rep2_TvsC" "IAA_minus_CTCF_rep3_TvsC" "IAA_plus_CTCF_rep1_TvsC" "IAA_plus_CTCF_rep2_TvsC" "IAA_plus_CTCF_rep3_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_CTCF_rep1_TvsC.bigwig IAA_minus_IP_CTCF_rep2_TvsC.bigwig IAA_minus_IP_CTCF_rep3_TvsC.bigwig IAA_plus_IP_CTCF_rep1_TvsC.bigwig IAA_plus_IP_CTCF_rep2_TvsC.bigwig IAA_plus_IP_CTCF_rep3_TvsC.bigwig -R ROI_IAA_minus.CTCF.summit.sorted.bed RIB_IAA_plus_and_minus.CTCF.summit.sorted.bed ROI_IAA_plus.CTCF.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_CTCF_ChIP_Seq_all_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_CTCF_ChIP_Seq_all_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_IAA72h_CTCF_ChIP_Seq_all_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_ROI_IAA_minus(32366)" "CTCF_RIB_IAA_plus_and_minus(4)" "CTCF_ROI_IAA_plus(131)" --samplesLabel "IAA_minus_CTCF_rep1_TvsC" "IAA_minus_CTCF_rep2_TvsC" "IAA_minus_CTCF_rep3_TvsC" "IAA_plus_CTCF_rep1_TvsC" "IAA_plus_CTCF_rep2_TvsC" "IAA_plus_CTCF_rep3_TvsC"

# use 293T 2w CTCF summit
# 上传 CTCF.summit.all.sorted.bed 到该目录
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/plotHeatMap_use_bamCompare_bigwig_for_CTCF
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_CTCF_rep1_TvsC.bigwig IAA_minus_IP_CTCF_rep2_TvsC.bigwig IAA_minus_IP_CTCF_rep3_TvsC.bigwig IAA_plus_IP_CTCF_rep1_TvsC.bigwig IAA_plus_IP_CTCF_rep2_TvsC.bigwig IAA_plus_IP_CTCF_rep3_TvsC.bigwig -R CTCF.summit.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_CTCF_ChIP_Seq_all_samples_for_plotHeatMap.293T2wCTCF.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_CTCF_ChIP_Seq_all_samples_for_plotHeatMap.293T2wCTCF.gz -o bamCompare_matrix_of_IAA72h_CTCF_ChIP_Seq_all_samples_for_plotHeatMap.293T2wCTCF.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_summits(20000)" --samplesLabel "IAA_minus_CTCF_rep1_TvsC" "IAA_minus_CTCF_rep2_TvsC" "IAA_minus_CTCF_rep3_TvsC" "IAA_plus_CTCF_rep1_TvsC" "IAA_plus_CTCF_rep2_TvsC" "IAA_plus_CTCF_rep3_TvsC"

```

### use rep123 bam together in one

#### make rep123 merge bam file links

```bash
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/plotHeatMap_use_bamCompare_bigwig_for_CTCF
# bam
ln -s ../HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam
ln -s ../HEK293T_IAA_minus_IP_CTCF_rep123.merged.sorted.bam HEK293T_IAA_minus_IP_CTCF_rep123.merged.sorted.bam
ln -s ../HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam
ln -s ../HEK293T_IAA_plus_IP_CTCF_rep123.merged.sorted.bam HEK293T_IAA_plus_IP_CTCF_rep123.merged.sorted.bam
# bai
ln -s ../HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam.bai HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_IAA_minus_IP_CTCF_rep123.merged.sorted.bam.bai HEK293T_IAA_minus_IP_CTCF_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam.bai HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_IAA_plus_IP_CTCF_rep123.merged.sorted.bam.bai HEK293T_IAA_plus_IP_CTCF_rep123.merged.sorted.bam.bai
ls -1 *_rep123.merged.sorted.bam
# HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam
# HEK293T_IAA_minus_IP_CTCF_rep123.merged.sorted.bam
# HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam
# HEK293T_IAA_plus_IP_CTCF_rep123.merged.sorted.bam

```

#### make log2 trans bigwig by bamCompare

```bash
bamCompare -b1 HEK293T_IAA_minus_IP_CTCF_rep123.merged.sorted.bam -b2 HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam -o IAA_minus_IP_CTCF_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
bamCompare -b1 HEK293T_IAA_plus_IP_CTCF_rep123.merged.sorted.bam -b2 HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam -o IAA_plus_IP_CTCF_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_CTCF_rep123_merged_TvsC.bigwig IAA_plus_IP_CTCF_rep123_merged_TvsC.bigwig -R CTCF_IAA_minus.only.summit.bed CTCF_IAA_minus_and_plus.common.summit.bed CTCF_IAA_plus.only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_CTCF_ChIP_Seq_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_CTCF_ChIP_Seq_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_IAA72h_CTCF_ChIP_Seq_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_IAA_minus_only(55328)" "common_summits(317)" "CTCF_IAA_plus_only(83)" --samplesLabel "IAA_minus_CTCF_rep123_TvsC" "IAA_plus_CTCF_rep123_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_CTCF_rep123_merged_TvsC.bigwig IAA_plus_IP_CTCF_rep123_merged_TvsC.bigwig -R ROI_IAA_minus.CTCF.summit.sorted.bed RIB_IAA_plus_and_minus.CTCF.summit.sorted.bed ROI_IAA_plus.CTCF.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_CTCF_ChIP_Seq_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_CTCF_ChIP_Seq_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_IAA72h_CTCF_ChIP_Seq_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_ROI_IAA_minus(32366)" "CTCF_RIB_IAA_plus_and_minus(4)" "CTCF_ROI_IAA_plus(131)" --samplesLabel "IAA_minus_CTCF_rep123_TvsC" "IAA_plus_CTCF_rep123_TvsC"

# use 293T 2w CTCF summit
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_CTCF_rep123_merged_TvsC.bigwig IAA_plus_IP_CTCF_rep123_merged_TvsC.bigwig -R CTCF.summit.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_CTCF_ChIP_Seq_rep123_merged_for_plotHeatMap.293T2wCTCF.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_CTCF_ChIP_Seq_rep123_merged_for_plotHeatMap.293T2wCTCF.gz -o bamCompare_matrix_of_IAA72h_CTCF_ChIP_Seq_rep123_merged_for_plotHeatMap.293T2wCTCF.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_summits(20000)" --samplesLabel "IAA_minus_CTCF_rep123_TvsC" "IAA_plus_CTCF_rep123_TvsC"

```

## MBD3 ChIP-Seq

### use individual bam of each repeat

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data
mkdir plotHeatMap_use_bamCompare_bigwig_for_MBD3
cd plotHeatMap_use_bamCompare_bigwig_for_MBD3
# 上传 6_make_IAA72h_MBD3_ChIP_Seq_sorted_bam_file_links.sh
# 创建 bam 和 bai 文件链接
bash 6_make_IAA72h_MBD3_ChIP_Seq_sorted_bam_file_links.sh
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_MBD3_rep123_together/IP_MBD3_rep123_TvsC_IAA_minus.only.summit.bed ./MBD3_IAA_minus.only.summit.bed
ln -s ../macs2_peak_calling_2_for_MBD3_rep123_together/IP_MBD3_rep123_TvsC_IAA_minus_and_plus.common.summit.bed ./MBD3_IAA_minus_and_plus.common.summit.bed
ln -s ../macs2_peak_calling_2_for_MBD3_rep123_together/IP_MBD3_rep123_TvsC_IAA_plus.only.summit.bed ./MBD3_IAA_plus.only.summit.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_MBD3_bdgdiff/RIB_IAA_plus_and_minus.MBD3.summit.sorted.bed ./RIB_IAA_plus_and_minus.MBD3.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_MBD3_bdgdiff/ROI_IAA_minus.MBD3.summit.sorted.bed ./ROI_IAA_minus.MBD3.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_MBD3_bdgdiff/ROI_IAA_plus.MBD3.summit.sorted.bed ./ROI_IAA_plus.MBD3.summit.sorted.bed
ls -1

```

#### creat log2 trans bigwig by using bamCompare 

```bash
# 构造 7_run_bam_comCompare_for_MBD3_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# bamCompare -b1 HEK293T_IAA_minus_IP_MBD3_rep1.bowtie2.sorted.bam -b2 HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam -o IAA_minus_IP_MBD3_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_IAA_minus_IP_MBD3_rep2.bowtie2.sorted.bam -b2 HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam -o IAA_minus_IP_MBD3_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_IAA_minus_IP_MBD3_rep3.bowtie2.sorted.bam -b2 HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam -o IAA_minus_IP_MBD3_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_IAA_plus_IP_MBD3_rep1.bowtie2.sorted.bam -b2 HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam -o IAA_plus_IP_MBD3_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_IAA_plus_IP_MBD3_rep2.bowtie2.sorted.bam -b2 HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam -o IAA_plus_IP_MBD3_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_IAA_plus_IP_MBD3_rep3.bowtie2.sorted.bam -b2 HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam -o IAA_plus_IP_MBD3_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
nohup bash 7_run_bam_comCompare_for_MBD3_each_repeat.sh >7_run_bam_comCompare_for_MBD3_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bam_comCompare_for_MBD3_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_MBD3_rep1_TvsC.bigwig IAA_minus_IP_MBD3_rep2_TvsC.bigwig IAA_minus_IP_MBD3_rep3_TvsC.bigwig IAA_plus_IP_MBD3_rep1_TvsC.bigwig IAA_plus_IP_MBD3_rep2_TvsC.bigwig IAA_plus_IP_MBD3_rep3_TvsC.bigwig -R MBD3_IAA_minus.only.summit.bed MBD3_IAA_minus_and_plus.common.summit.bed MBD3_IAA_plus.only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_MBD3_ChIP_Seq_all_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_MBD3_ChIP_Seq_all_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_IAA72h_MBD3_ChIP_Seq_all_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_IAA_minus_only(16929)" "common_summits(56)" "MBD3_IAA_plus_only(168)" --samplesLabel "IAA_minus_MBD3_rep1_TvsC" "IAA_minus_MBD3_rep2_TvsC" "IAA_minus_MBD3_rep3_TvsC" "IAA_plus_MBD3_rep1_TvsC" "IAA_plus_MBD3_rep2_TvsC" "IAA_plus_MBD3_rep3_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_MBD3_rep1_TvsC.bigwig IAA_minus_IP_MBD3_rep2_TvsC.bigwig IAA_minus_IP_MBD3_rep3_TvsC.bigwig IAA_plus_IP_MBD3_rep1_TvsC.bigwig IAA_plus_IP_MBD3_rep2_TvsC.bigwig IAA_plus_IP_MBD3_rep3_TvsC.bigwig -R ROI_IAA_minus.MBD3.summit.sorted.bed RIB_IAA_plus_and_minus.MBD3.summit.sorted.bed ROI_IAA_plus.MBD3.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_MBD3_ChIP_Seq_all_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_MBD3_ChIP_Seq_all_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_IAA72h_MBD3_ChIP_Seq_all_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_ROI_IAA_minus(2292)" "MBD3_RIB_IAA_plus_and_minus(6)" "MBD3_ROI_IAA_plus(323)" --samplesLabel "IAA_minus_MBD3_rep1_TvsC" "IAA_minus_MBD3_rep2_TvsC" "IAA_minus_MBD3_rep3_TvsC" "IAA_plus_MBD3_rep1_TvsC" "IAA_plus_MBD3_rep2_TvsC" "IAA_plus_MBD3_rep3_TvsC"

```

### use rep123 bam together in one

#### make rep123 merge bam file links

```bash
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/plotHeatMap_use_bamCompare_bigwig_for_MBD3
# bam
ln -s ../HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam
ln -s ../HEK293T_IAA_minus_IP_MBD3_rep123.merged.sorted.bam HEK293T_IAA_minus_IP_MBD3_rep123.merged.sorted.bam
ln -s ../HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam
ln -s ../HEK293T_IAA_plus_IP_MBD3_rep123.merged.sorted.bam HEK293T_IAA_plus_IP_MBD3_rep123.merged.sorted.bam
# bai
ln -s ../HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam.bai HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_IAA_minus_IP_MBD3_rep123.merged.sorted.bam.bai HEK293T_IAA_minus_IP_MBD3_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam.bai HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_IAA_plus_IP_MBD3_rep123.merged.sorted.bam.bai HEK293T_IAA_plus_IP_MBD3_rep123.merged.sorted.bam.bai
ls -1 *_rep123.merged.sorted.bam
# HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam
# HEK293T_IAA_minus_IP_MBD3_rep123.merged.sorted.bam
# HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam
# HEK293T_IAA_plus_IP_MBD3_rep123.merged.sorted.bam

```

#### make log2 trans bigwig by bamCompare

```bash
bamCompare -b1 HEK293T_IAA_minus_IP_MBD3_rep123.merged.sorted.bam -b2 HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam -o IAA_minus_IP_MBD3_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
bamCompare -b1 HEK293T_IAA_plus_IP_MBD3_rep123.merged.sorted.bam -b2 HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam -o IAA_plus_IP_MBD3_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_MBD3_rep123_merged_TvsC.bigwig IAA_plus_IP_MBD3_rep123_merged_TvsC.bigwig -R MBD3_IAA_minus.only.summit.bed MBD3_IAA_minus_and_plus.common.summit.bed MBD3_IAA_plus.only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_MBD3_ChIP_Seq_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_MBD3_ChIP_Seq_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_IAA72h_MBD3_ChIP_Seq_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_IAA_minus_only(16929)" "common_summits(56)" "MBD3_IAA_plus_only(168)" --samplesLabel "IAA_minus_MBD3_rep123_TvsC" "IAA_plus_MBD3_rep123_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_MBD3_rep123_merged_TvsC.bigwig IAA_plus_IP_MBD3_rep123_merged_TvsC.bigwig -R ROI_IAA_minus.MBD3.summit.sorted.bed RIB_IAA_plus_and_minus.MBD3.summit.sorted.bed ROI_IAA_plus.MBD3.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_MBD3_ChIP_Seq_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_MBD3_ChIP_Seq_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_IAA72h_MBD3_ChIP_Seq_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_ROI_IAA_minus(2292)" "MBD3_RIB_IAA_plus_and_minus(6)" "MBD3_ROI_IAA_plus(323)" --samplesLabel "IAA_minus_MBD3_rep123_TvsC" "IAA_plus_MBD3_rep123_TvsC"

```
#### use 2w CTCF summit

```bash
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data
mkdir reads_heatmap_of_MBD3_around_2w_CTCF_summits
# 上传 2w CTCF summit, 创建 bigwig 文件连接，绘图
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/reads_heatmap_of_MBD3_around_2w_CTCF_summits
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_MBD3/IAA_minus_IP_MBD3_rep123_merged_TvsC.bigwig IAA_minus_IP_MBD3_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_MBD3/IAA_minus_IP_MBD3_rep1_TvsC.bigwig IAA_minus_IP_MBD3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_MBD3/IAA_minus_IP_MBD3_rep2_TvsC.bigwig IAA_minus_IP_MBD3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_MBD3/IAA_minus_IP_MBD3_rep3_TvsC.bigwig IAA_minus_IP_MBD3_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_MBD3/IAA_plus_IP_MBD3_rep123_merged_TvsC.bigwig IAA_plus_IP_MBD3_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_MBD3/IAA_plus_IP_MBD3_rep1_TvsC.bigwig IAA_plus_IP_MBD3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_MBD3/IAA_plus_IP_MBD3_rep2_TvsC.bigwig IAA_plus_IP_MBD3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_MBD3/IAA_plus_IP_MBD3_rep3_TvsC.bigwig IAA_plus_IP_MBD3_rep3_TvsC.bigwig

# use 2w CTCF summits for all samples
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_MBD3_rep1_TvsC.bigwig IAA_minus_IP_MBD3_rep2_TvsC.bigwig IAA_minus_IP_MBD3_rep3_TvsC.bigwig IAA_plus_IP_MBD3_rep1_TvsC.bigwig IAA_plus_IP_MBD3_rep2_TvsC.bigwig IAA_plus_IP_MBD3_rep3_TvsC.bigwig -R CTCF.summit.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_MBD3_ChIP_Seq_all_samples_for_plotHeatMap.2wCTCFsummit.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_MBD3_ChIP_Seq_all_samples_for_plotHeatMap.2wCTCFsummit.gz -o bamCompare_matrix_of_IAA72h_MBD3_ChIP_Seq_all_samples_for_plotHeatMap.2wCTCFsummit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF peak summits(20,000)" --samplesLabel "IAA_minus_MBD3_rep1_TvsC" "IAA_minus_MBD3_rep2_TvsC" "IAA_minus_MBD3_rep3_TvsC" "IAA_plus_MBD3_rep1_TvsC" "IAA_plus_MBD3_rep2_TvsC" "IAA_plus_MBD3_rep3_TvsC"

# use 2w CTCF summits for rep123 merged
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_MBD3_rep123_merged_TvsC.bigwig IAA_plus_IP_MBD3_rep123_merged_TvsC.bigwig -R CTCF.summit.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_MBD3_ChIP_Seq_rep123_merged_for_plotHeatMap.2wCTCFsummit.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_MBD3_ChIP_Seq_rep123_merged_for_plotHeatMap.2wCTCFsummit.gz -o bamCompare_matrix_of_IAA72h_MBD3_ChIP_Seq_rep123_merged_for_plotHeatMap.2wCTCFsummit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF peak summits(20,000)" --samplesLabel "IAA_minus_MBD3_rep123_TvsC" "IAA_plus_MBD3_rep123_TvsC"

```



## TET1 ChIP-Seq

### use individual bam of each repeat

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data
mkdir plotHeatMap_use_bamCompare_bigwig_for_TET1
cd plotHeatMap_use_bamCompare_bigwig_for_TET1
# 上传 6_make_IAA72h_TET1_ChIP_Seq_sorted_bam_file_links.sh
# 创建 bam 和 bai 文件链接
bash 6_make_IAA72h_TET1_ChIP_Seq_sorted_bam_file_links.sh
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_TET1_rep123_together/IP_TET1_rep123_TvsC_IAA_minus.only.summit.bed ./TET1_IAA_minus.only.summit.bed
ln -s ../macs2_peak_calling_2_for_TET1_rep123_together/IP_TET1_rep123_TvsC_IAA_minus_and_plus.common.summit.bed ./TET1_IAA_minus_and_plus.common.summit.bed
ln -s ../macs2_peak_calling_2_for_TET1_rep123_together/IP_TET1_rep123_TvsC_IAA_plus.only.summit.bed ./TET1_IAA_plus.only.summit.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_TET1_bdgdiff/RIB_IAA_plus_and_minus.TET1.summit.sorted.bed ./RIB_IAA_plus_and_minus.TET1.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_TET1_bdgdiff/ROI_IAA_minus.TET1.summit.sorted.bed ./ROI_IAA_minus.TET1.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_TET1_bdgdiff/ROI_IAA_plus.TET1.summit.sorted.bed ./ROI_IAA_plus.TET1.summit.sorted.bed
ls -1

```

#### creat log2 trans bigwig by using bamCompare 

```bash
# 构造 7_run_bam_comCompare_for_TET1_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# bamCompare -b1 HEK293T_IAA_minus_IP_TET1_rep1.bowtie2.sorted.bam -b2 HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam -o IAA_minus_IP_TET1_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_IAA_minus_IP_TET1_rep2.bowtie2.sorted.bam -b2 HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam -o IAA_minus_IP_TET1_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_IAA_minus_IP_TET1_rep3.bowtie2.sorted.bam -b2 HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam -o IAA_minus_IP_TET1_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_IAA_plus_IP_TET1_rep1.bowtie2.sorted.bam -b2 HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam -o IAA_plus_IP_TET1_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_IAA_plus_IP_TET1_rep2.bowtie2.sorted.bam -b2 HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam -o IAA_plus_IP_TET1_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_IAA_plus_IP_TET1_rep3.bowtie2.sorted.bam -b2 HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam -o IAA_plus_IP_TET1_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
nohup bash 7_run_bam_comCompare_for_TET1_each_repeat.sh >7_run_bam_comCompare_for_TET1_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bam_comCompare_for_TET1_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_TET1_rep1_TvsC.bigwig IAA_minus_IP_TET1_rep2_TvsC.bigwig IAA_minus_IP_TET1_rep3_TvsC.bigwig IAA_plus_IP_TET1_rep1_TvsC.bigwig IAA_plus_IP_TET1_rep2_TvsC.bigwig IAA_plus_IP_TET1_rep3_TvsC.bigwig -R TET1_IAA_minus.only.summit.bed TET1_IAA_minus_and_plus.common.summit.bed TET1_IAA_plus.only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_TET1_ChIP_Seq_all_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_TET1_ChIP_Seq_all_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_IAA72h_TET1_ChIP_Seq_all_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "TET1_IAA_minus_only(25754)" "common_summits(138)" "TET1_IAA_plus_only(1205)" --samplesLabel "IAA_minus_TET1_rep1_TvsC" "IAA_minus_TET1_rep2_TvsC" "IAA_minus_TET1_rep3_TvsC" "IAA_plus_TET1_rep1_TvsC" "IAA_plus_TET1_rep2_TvsC" "IAA_plus_TET1_rep3_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_TET1_rep1_TvsC.bigwig IAA_minus_IP_TET1_rep2_TvsC.bigwig IAA_minus_IP_TET1_rep3_TvsC.bigwig IAA_plus_IP_TET1_rep1_TvsC.bigwig IAA_plus_IP_TET1_rep2_TvsC.bigwig IAA_plus_IP_TET1_rep3_TvsC.bigwig -R ROI_IAA_minus.TET1.summit.sorted.bed RIB_IAA_plus_and_minus.TET1.summit.sorted.bed ROI_IAA_plus.TET1.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_TET1_ChIP_Seq_all_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_TET1_ChIP_Seq_all_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_IAA72h_TET1_ChIP_Seq_all_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "TET1_ROI_IAA_minus(6684)" "TET1_RIB_IAA_plus_and_minus(15)" "TET1_ROI_IAA_plus(879)" --samplesLabel "IAA_minus_TET1_rep1_TvsC" "IAA_minus_TET1_rep2_TvsC" "IAA_minus_TET1_rep3_TvsC" "IAA_plus_TET1_rep1_TvsC" "IAA_plus_TET1_rep2_TvsC" "IAA_plus_TET1_rep3_TvsC"

```

### use rep123 bam together in one

#### make rep123 merge bam file links

```bash
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/plotHeatMap_use_bamCompare_bigwig_for_TET1
# bam
ln -s ../HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam
ln -s ../HEK293T_IAA_minus_IP_TET1_rep123.merged.sorted.bam HEK293T_IAA_minus_IP_TET1_rep123.merged.sorted.bam
ln -s ../HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam
ln -s ../HEK293T_IAA_plus_IP_TET1_rep123.merged.sorted.bam HEK293T_IAA_plus_IP_TET1_rep123.merged.sorted.bam
# bai
ln -s ../HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam.bai HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_IAA_minus_IP_TET1_rep123.merged.sorted.bam.bai HEK293T_IAA_minus_IP_TET1_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam.bai HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_IAA_plus_IP_TET1_rep123.merged.sorted.bam.bai HEK293T_IAA_plus_IP_TET1_rep123.merged.sorted.bam.bai
ls -1 *_rep123.merged.sorted.bam
# HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam
# HEK293T_IAA_minus_IP_TET1_rep123.merged.sorted.bam
# HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam
# HEK293T_IAA_plus_IP_TET1_rep123.merged.sorted.bam

```

#### make log2 trans bigwig by bamCompare

```bash
bamCompare -b1 HEK293T_IAA_minus_IP_TET1_rep123.merged.sorted.bam -b2 HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam -o IAA_minus_IP_TET1_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
bamCompare -b1 HEK293T_IAA_plus_IP_TET1_rep123.merged.sorted.bam -b2 HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam -o IAA_plus_IP_TET1_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_TET1_rep123_merged_TvsC.bigwig IAA_plus_IP_TET1_rep123_merged_TvsC.bigwig -R TET1_IAA_minus.only.summit.bed TET1_IAA_minus_and_plus.common.summit.bed TET1_IAA_plus.only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_TET1_ChIP_Seq_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_TET1_ChIP_Seq_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_IAA72h_TET1_ChIP_Seq_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "TET1_IAA_minus_only(25754)" "common_summits(138)" "TET1_IAA_plus_only(1205)" --samplesLabel "IAA_minus_TET1_rep123_TvsC" "IAA_plus_TET1_rep123_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_TET1_rep123_merged_TvsC.bigwig IAA_plus_IP_TET1_rep123_merged_TvsC.bigwig -R ROI_IAA_minus.TET1.summit.sorted.bed RIB_IAA_plus_and_minus.TET1.summit.sorted.bed ROI_IAA_plus.TET1.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_TET1_ChIP_Seq_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_TET1_ChIP_Seq_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_IAA72h_TET1_ChIP_Seq_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "TET1_ROI_IAA_minus(6684)" "TET1_RIB_IAA_plus_and_minus(15)" "TET1_ROI_IAA_plus(879)" --samplesLabel "IAA_minus_TET1_rep123_TvsC" "IAA_plus_TET1_rep123_TvsC"

```
#### use 2w CTCF summit

```bash
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data
mkdir reads_heatmap_of_TET1_around_2w_CTCF_summits
# 上传 2w CTCF summit, 创建 bigwig 文件连接，绘图
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/reads_heatmap_of_TET1_around_2w_CTCF_summits
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET1/IAA_minus_IP_TET1_rep123_merged_TvsC.bigwig IAA_minus_IP_TET1_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET1/IAA_minus_IP_TET1_rep1_TvsC.bigwig IAA_minus_IP_TET1_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET1/IAA_minus_IP_TET1_rep2_TvsC.bigwig IAA_minus_IP_TET1_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET1/IAA_minus_IP_TET1_rep3_TvsC.bigwig IAA_minus_IP_TET1_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET1/IAA_plus_IP_TET1_rep123_merged_TvsC.bigwig IAA_plus_IP_TET1_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET1/IAA_plus_IP_TET1_rep1_TvsC.bigwig IAA_plus_IP_TET1_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET1/IAA_plus_IP_TET1_rep2_TvsC.bigwig IAA_plus_IP_TET1_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET1/IAA_plus_IP_TET1_rep3_TvsC.bigwig IAA_plus_IP_TET1_rep3_TvsC.bigwig

# use 2w CTCF summits for all samples
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_TET1_rep1_TvsC.bigwig IAA_minus_IP_TET1_rep2_TvsC.bigwig IAA_minus_IP_TET1_rep3_TvsC.bigwig IAA_plus_IP_TET1_rep1_TvsC.bigwig IAA_plus_IP_TET1_rep2_TvsC.bigwig IAA_plus_IP_TET1_rep3_TvsC.bigwig -R CTCF.summit.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_TET1_ChIP_Seq_all_samples_for_plotHeatMap.2wCTCFsummit.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_TET1_ChIP_Seq_all_samples_for_plotHeatMap.2wCTCFsummit.gz -o bamCompare_matrix_of_IAA72h_TET1_ChIP_Seq_all_samples_for_plotHeatMap.2wCTCFsummit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF peak summits(20,000)" --samplesLabel "IAA_minus_TET1_rep1_TvsC" "IAA_minus_TET1_rep2_TvsC" "IAA_minus_TET1_rep3_TvsC" "IAA_plus_TET1_rep1_TvsC" "IAA_plus_TET1_rep2_TvsC" "IAA_plus_TET1_rep3_TvsC"

# use 2w CTCF summits for rep123 merged
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_TET1_rep123_merged_TvsC.bigwig IAA_plus_IP_TET1_rep123_merged_TvsC.bigwig -R CTCF.summit.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_TET1_ChIP_Seq_rep123_merged_for_plotHeatMap.2wCTCFsummit.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_TET1_ChIP_Seq_rep123_merged_for_plotHeatMap.2wCTCFsummit.gz -o bamCompare_matrix_of_IAA72h_TET1_ChIP_Seq_rep123_merged_for_plotHeatMap.2wCTCFsummit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF peak summits(20,000)" --samplesLabel "IAA_minus_TET1_rep123_TvsC" "IAA_plus_TET1_rep123_TvsC"

```



## TET2 ChIP-Seq

### use individual bam of each repeat

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data
mkdir plotHeatMap_use_bamCompare_bigwig_for_TET2
cd plotHeatMap_use_bamCompare_bigwig_for_TET2
# 上传 6_make_IAA72h_TET2_ChIP_Seq_sorted_bam_file_links.sh
# 创建 bam 和 bai 文件链接
bash 6_make_IAA72h_TET2_ChIP_Seq_sorted_bam_file_links.sh
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_TET2_rep123_together/IP_TET2_rep123_TvsC_IAA_minus.only.summit.bed ./TET2_IAA_minus.only.summit.bed
ln -s ../macs2_peak_calling_2_for_TET2_rep123_together/IP_TET2_rep123_TvsC_IAA_minus_and_plus.common.summit.bed ./TET2_IAA_minus_and_plus.common.summit.bed
ln -s ../macs2_peak_calling_2_for_TET2_rep123_together/IP_TET2_rep123_TvsC_IAA_plus.only.summit.bed ./TET2_IAA_plus.only.summit.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_TET2_bdgdiff/RIB_IAA_plus_and_minus.TET2.summit.sorted.bed ./RIB_IAA_plus_and_minus.TET2.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_TET2_bdgdiff/ROI_IAA_minus.TET2.summit.sorted.bed ./ROI_IAA_minus.TET2.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_TET2_bdgdiff/ROI_IAA_plus.TET2.summit.sorted.bed ./ROI_IAA_plus.TET2.summit.sorted.bed
ls -1

```

#### creat log2 trans bigwig by using bamCompare 

```bash
# 构造 7_run_bam_comCompare_for_TET2_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# bamCompare -b1 HEK293T_IAA_minus_IP_TET2_rep1.bowtie2.sorted.bam -b2 HEK293T_IAA_minus_INPUT_rep1.bowtie2.sorted.bam -o IAA_minus_IP_TET2_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_IAA_minus_IP_TET2_rep2.bowtie2.sorted.bam -b2 HEK293T_IAA_minus_INPUT_rep2.bowtie2.sorted.bam -o IAA_minus_IP_TET2_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_IAA_minus_IP_TET2_rep3.bowtie2.sorted.bam -b2 HEK293T_IAA_minus_INPUT_rep3.bowtie2.sorted.bam -o IAA_minus_IP_TET2_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_IAA_plus_IP_TET2_rep1.bowtie2.sorted.bam -b2 HEK293T_IAA_plus_INPUT_rep1.bowtie2.sorted.bam -o IAA_plus_IP_TET2_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_IAA_plus_IP_TET2_rep2.bowtie2.sorted.bam -b2 HEK293T_IAA_plus_INPUT_rep2.bowtie2.sorted.bam -o IAA_plus_IP_TET2_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
# bamCompare -b1 HEK293T_IAA_plus_IP_TET2_rep3.bowtie2.sorted.bam -b2 HEK293T_IAA_plus_INPUT_rep3.bowtie2.sorted.bam -o IAA_plus_IP_TET2_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
nohup bash 7_run_bam_comCompare_for_TET2_each_repeat.sh >7_run_bam_comCompare_for_TET2_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bam_comCompare_for_TET2_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_TET2_rep1_TvsC.bigwig IAA_minus_IP_TET2_rep2_TvsC.bigwig IAA_minus_IP_TET2_rep3_TvsC.bigwig IAA_plus_IP_TET2_rep1_TvsC.bigwig IAA_plus_IP_TET2_rep2_TvsC.bigwig IAA_plus_IP_TET2_rep3_TvsC.bigwig -R TET2_IAA_minus.only.summit.bed TET2_IAA_minus_and_plus.common.summit.bed TET2_IAA_plus.only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_TET2_ChIP_Seq_all_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_TET2_ChIP_Seq_all_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_IAA72h_TET2_ChIP_Seq_all_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "TET2_IAA_minus_only(6695)" "common_summits(37)" "TET2_IAA_plus_only(9165)" --samplesLabel "IAA_minus_TET2_rep1_TvsC" "IAA_minus_TET2_rep2_TvsC" "IAA_minus_TET2_rep3_TvsC" "IAA_plus_TET2_rep1_TvsC" "IAA_plus_TET2_rep2_TvsC" "IAA_plus_TET2_rep3_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_TET2_rep1_TvsC.bigwig IAA_minus_IP_TET2_rep2_TvsC.bigwig IAA_minus_IP_TET2_rep3_TvsC.bigwig IAA_plus_IP_TET2_rep1_TvsC.bigwig IAA_plus_IP_TET2_rep2_TvsC.bigwig IAA_plus_IP_TET2_rep3_TvsC.bigwig -R ROI_IAA_minus.TET2.summit.sorted.bed RIB_IAA_plus_and_minus.TET2.summit.sorted.bed ROI_IAA_plus.TET2.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_TET2_ChIP_Seq_all_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_TET2_ChIP_Seq_all_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_IAA72h_TET2_ChIP_Seq_all_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "TET2_ROI_IAA_minus(261)" "TET2_RIB_IAA_plus_and_minus(3)" "TET2_ROI_IAA_plus(3676)" --samplesLabel "IAA_minus_TET2_rep1_TvsC" "IAA_minus_TET2_rep2_TvsC" "IAA_minus_TET2_rep3_TvsC" "IAA_plus_TET2_rep1_TvsC" "IAA_plus_TET2_rep2_TvsC" "IAA_plus_TET2_rep3_TvsC"

```

### use rep123 bam together in one

#### make rep123 merge bam file links

```bash
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/plotHeatMap_use_bamCompare_bigwig_for_TET2
# bam
ln -s ../HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam
ln -s ../HEK293T_IAA_minus_IP_TET2_rep123.merged.sorted.bam HEK293T_IAA_minus_IP_TET2_rep123.merged.sorted.bam
ln -s ../HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam
ln -s ../HEK293T_IAA_plus_IP_TET2_rep123.merged.sorted.bam HEK293T_IAA_plus_IP_TET2_rep123.merged.sorted.bam
# bai
ln -s ../HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam.bai HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_IAA_minus_IP_TET2_rep123.merged.sorted.bam.bai HEK293T_IAA_minus_IP_TET2_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam.bai HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam.bai
ln -s ../HEK293T_IAA_plus_IP_TET2_rep123.merged.sorted.bam.bai HEK293T_IAA_plus_IP_TET2_rep123.merged.sorted.bam.bai
ls -1 *_rep123.merged.sorted.bam
# HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam
# HEK293T_IAA_minus_IP_TET2_rep123.merged.sorted.bam
# HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam
# HEK293T_IAA_plus_IP_TET2_rep123.merged.sorted.bam

```

#### make log2 trans bigwig by bamCompare

```bash
bamCompare -b1 HEK293T_IAA_minus_IP_TET2_rep123.merged.sorted.bam -b2 HEK293T_IAA_minus_INPUT_rep123.merged.sorted.bam -o IAA_minus_IP_TET2_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads
bamCompare -b1 HEK293T_IAA_plus_IP_TET2_rep123.merged.sorted.bam -b2 HEK293T_IAA_plus_INPUT_rep123.merged.sorted.bam -o IAA_plus_IP_TET2_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --ignoreDuplicates --centerReads

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_TET2_rep123_merged_TvsC.bigwig IAA_plus_IP_TET2_rep123_merged_TvsC.bigwig -R TET2_IAA_minus.only.summit.bed TET2_IAA_minus_and_plus.common.summit.bed TET2_IAA_plus.only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_TET2_ChIP_Seq_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_TET2_ChIP_Seq_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_IAA72h_TET2_ChIP_Seq_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "TET2_IAA_minus_only(6695)" "common_summits(37)" "TET2_IAA_plus_only(9165)" --samplesLabel "IAA_minus_TET2_rep123_TvsC" "IAA_plus_TET2_rep123_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_TET2_rep123_merged_TvsC.bigwig IAA_plus_IP_TET2_rep123_merged_TvsC.bigwig -R ROI_IAA_minus.TET2.summit.sorted.bed RIB_IAA_plus_and_minus.TET2.summit.sorted.bed ROI_IAA_plus.TET2.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_TET2_ChIP_Seq_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_TET2_ChIP_Seq_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_IAA72h_TET2_ChIP_Seq_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "TET2_ROI_IAA_minus(261)" "TET2_RIB_IAA_plus_and_minus(3)" "TET2_ROI_IAA_plus(3676)" --samplesLabel "IAA_minus_TET2_rep123_TvsC" "IAA_plus_TET2_rep123_TvsC"

```

#### use 2w CTCF summit

```bash
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data
mkdir reads_heatmap_of_TET2_around_2w_CTCF_summits
# 上传 2w CTCF summit, 创建 bigwig 文件连接，绘图
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/reads_heatmap_of_TET2_around_2w_CTCF_summits
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET2/IAA_minus_IP_TET2_rep123_merged_TvsC.bigwig IAA_minus_IP_TET2_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET2/IAA_minus_IP_TET2_rep1_TvsC.bigwig IAA_minus_IP_TET2_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET2/IAA_minus_IP_TET2_rep2_TvsC.bigwig IAA_minus_IP_TET2_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET2/IAA_minus_IP_TET2_rep3_TvsC.bigwig IAA_minus_IP_TET2_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET2/IAA_plus_IP_TET2_rep123_merged_TvsC.bigwig IAA_plus_IP_TET2_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET2/IAA_plus_IP_TET2_rep1_TvsC.bigwig IAA_plus_IP_TET2_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET2/IAA_plus_IP_TET2_rep2_TvsC.bigwig IAA_plus_IP_TET2_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET2/IAA_plus_IP_TET2_rep3_TvsC.bigwig IAA_plus_IP_TET2_rep3_TvsC.bigwig

# use 2w CTCF summits for all samples
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_TET2_rep1_TvsC.bigwig IAA_minus_IP_TET2_rep2_TvsC.bigwig IAA_minus_IP_TET2_rep3_TvsC.bigwig IAA_plus_IP_TET2_rep1_TvsC.bigwig IAA_plus_IP_TET2_rep2_TvsC.bigwig IAA_plus_IP_TET2_rep3_TvsC.bigwig -R CTCF.summit.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_TET2_ChIP_Seq_all_samples_for_plotHeatMap.2wCTCFsummit.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_TET2_ChIP_Seq_all_samples_for_plotHeatMap.2wCTCFsummit.gz -o bamCompare_matrix_of_IAA72h_TET2_ChIP_Seq_all_samples_for_plotHeatMap.2wCTCFsummit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF peak summits(20,000)" --samplesLabel "IAA_minus_TET2_rep1_TvsC" "IAA_minus_TET2_rep2_TvsC" "IAA_minus_TET2_rep3_TvsC" "IAA_plus_TET2_rep1_TvsC" "IAA_plus_TET2_rep2_TvsC" "IAA_plus_TET2_rep3_TvsC"

# use 2w CTCF summits for rep123 merged
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_TET2_rep123_merged_TvsC.bigwig IAA_plus_IP_TET2_rep123_merged_TvsC.bigwig -R CTCF.summit.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_TET2_ChIP_Seq_rep123_merged_for_plotHeatMap.2wCTCFsummit.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_TET2_ChIP_Seq_rep123_merged_for_plotHeatMap.2wCTCFsummit.gz -o bamCompare_matrix_of_IAA72h_TET2_ChIP_Seq_rep123_merged_for_plotHeatMap.2wCTCFsummit.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF peak summits(20,000)" --samplesLabel "IAA_minus_TET2_rep123_TvsC" "IAA_plus_TET2_rep123_TvsC"

```



# plot reads enrichment heatmap use public CTCF rep12 union summit

## data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data
mkdir plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit

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
# 24356 CTCF_rep12_union.phq1.summit.sorted.bed
# 24356 CTCF_rep12_union.phq2.summit.sorted.bed
# 24356 CTCF_rep12_union.phq3.summit.sorted.bed
# 24356 CTCF_rep12_union.phq4.summit.sorted.bed
# 31737 CTCF_rep12_union.rep12_common.summit.sorted.bed
# 47355 CTCF_rep12_union.rep1_only.summit.sorted.bed
# 18332 CTCF_rep12_union.rep2_only.summit.sorted.bed

# 创建全部所需的 bigwig 文件的软连接
## CTCF
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_CTCF/IAA_minus_IP_CTCF_rep123_merged_TvsC.bigwig IAA_minus_IP_CTCF_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_CTCF/IAA_plus_IP_CTCF_rep123_merged_TvsC.bigwig IAA_plus_IP_CTCF_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_CTCF/IAA_minus_IP_CTCF_rep1_TvsC.bigwig IAA_minus_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_CTCF/IAA_minus_IP_CTCF_rep2_TvsC.bigwig IAA_minus_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_CTCF/IAA_minus_IP_CTCF_rep3_TvsC.bigwig IAA_minus_IP_CTCF_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_CTCF/IAA_plus_IP_CTCF_rep1_TvsC.bigwig IAA_plus_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_CTCF/IAA_plus_IP_CTCF_rep2_TvsC.bigwig IAA_plus_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_CTCF/IAA_plus_IP_CTCF_rep3_TvsC.bigwig IAA_plus_IP_CTCF_rep3_TvsC.bigwig
## MBD3
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_MBD3/IAA_minus_IP_MBD3_rep123_merged_TvsC.bigwig IAA_minus_IP_MBD3_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_MBD3/IAA_plus_IP_MBD3_rep123_merged_TvsC.bigwig IAA_plus_IP_MBD3_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_MBD3/IAA_minus_IP_MBD3_rep1_TvsC.bigwig IAA_minus_IP_MBD3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_MBD3/IAA_minus_IP_MBD3_rep2_TvsC.bigwig IAA_minus_IP_MBD3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_MBD3/IAA_minus_IP_MBD3_rep3_TvsC.bigwig IAA_minus_IP_MBD3_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_MBD3/IAA_plus_IP_MBD3_rep1_TvsC.bigwig IAA_plus_IP_MBD3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_MBD3/IAA_plus_IP_MBD3_rep2_TvsC.bigwig IAA_plus_IP_MBD3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_MBD3/IAA_plus_IP_MBD3_rep3_TvsC.bigwig IAA_plus_IP_MBD3_rep3_TvsC.bigwig
## TET1
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET1/IAA_minus_IP_TET1_rep123_merged_TvsC.bigwig IAA_minus_IP_TET1_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET1/IAA_plus_IP_TET1_rep123_merged_TvsC.bigwig IAA_plus_IP_TET1_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET1/IAA_minus_IP_TET1_rep1_TvsC.bigwig IAA_minus_IP_TET1_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET1/IAA_minus_IP_TET1_rep2_TvsC.bigwig IAA_minus_IP_TET1_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET1/IAA_minus_IP_TET1_rep3_TvsC.bigwig IAA_minus_IP_TET1_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET1/IAA_plus_IP_TET1_rep1_TvsC.bigwig IAA_plus_IP_TET1_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET1/IAA_plus_IP_TET1_rep2_TvsC.bigwig IAA_plus_IP_TET1_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET1/IAA_plus_IP_TET1_rep3_TvsC.bigwig IAA_plus_IP_TET1_rep3_TvsC.bigwig
## TET2
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET2/IAA_minus_IP_TET2_rep123_merged_TvsC.bigwig IAA_minus_IP_TET2_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET2/IAA_plus_IP_TET2_rep123_merged_TvsC.bigwig IAA_plus_IP_TET2_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET2/IAA_minus_IP_TET2_rep1_TvsC.bigwig IAA_minus_IP_TET2_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET2/IAA_minus_IP_TET2_rep2_TvsC.bigwig IAA_minus_IP_TET2_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET2/IAA_minus_IP_TET2_rep3_TvsC.bigwig IAA_minus_IP_TET2_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET2/IAA_plus_IP_TET2_rep1_TvsC.bigwig IAA_plus_IP_TET2_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET2/IAA_plus_IP_TET2_rep2_TvsC.bigwig IAA_plus_IP_TET2_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_bamCompare_bigwig_for_TET2/IAA_plus_IP_TET2_rep3_TvsC.bigwig IAA_plus_IP_TET2_rep3_TvsC.bigwig
ls -1 *.bigwig

```

## plot use CTCF rep12 union all summits

```bash
# plot use CTCF rep12 union all summits
## compute and plot use rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_CTCF_rep1_TvsC.bigwig IAA_minus_IP_CTCF_rep2_TvsC.bigwig IAA_minus_IP_CTCF_rep3_TvsC.bigwig IAA_plus_IP_CTCF_rep1_TvsC.bigwig IAA_plus_IP_CTCF_rep2_TvsC.bigwig IAA_plus_IP_CTCF_rep3_TvsC.bigwig IAA_minus_IP_MBD3_rep1_TvsC.bigwig IAA_minus_IP_MBD3_rep2_TvsC.bigwig IAA_minus_IP_MBD3_rep3_TvsC.bigwig IAA_plus_IP_MBD3_rep1_TvsC.bigwig IAA_plus_IP_MBD3_rep2_TvsC.bigwig IAA_plus_IP_MBD3_rep3_TvsC.bigwig IAA_minus_IP_TET1_rep1_TvsC.bigwig IAA_minus_IP_TET1_rep2_TvsC.bigwig IAA_minus_IP_TET1_rep3_TvsC.bigwig IAA_plus_IP_TET1_rep1_TvsC.bigwig IAA_plus_IP_TET1_rep2_TvsC.bigwig IAA_plus_IP_TET1_rep3_TvsC.bigwig IAA_minus_IP_TET2_rep1_TvsC.bigwig IAA_minus_IP_TET2_rep2_TvsC.bigwig IAA_minus_IP_TET2_rep3_TvsC.bigwig IAA_plus_IP_TET2_rep1_TvsC.bigwig IAA_plus_IP_TET2_rep2_TvsC.bigwig IAA_plus_IP_TET2_rep3_TvsC.bigwig -R CTCF_rep12_narrowPeak.summit.union.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_ChIP_all_samples_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_ChIP_all_samples_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz -o bamCompare_matrix_of_IAA72h_ChIP_all_samples_for_plotHeatMap.public_CTCF_rep12_union_all_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_rep12_union_all_summits(97424)" --samplesLabel "IAA_minus_CTCF_rep1" "IAA_minus_CTCF_rep2" "IAA_minus_CTCF_rep3" "IAA_plus_CTCF_rep1" "IAA_plus_CTCF_rep2" "IAA_plus_CTCF_rep3" "IAA_minus_MBD3_rep1" "IAA_minus_MBD3_rep2" "IAA_minus_MBD3_rep3" "IAA_plus_MBD3_rep1" "IAA_plus_MBD3_rep2" "IAA_plus_MBD3_rep3" "IAA_minus_TET1_rep1" "IAA_minus_TET1_rep2" "IAA_minus_TET1_rep3" "IAA_plus_TET1_rep1" "IAA_plus_TET1_rep2" "IAA_plus_TET1_rep3" "IAA_minus_TET2_rep1" "IAA_minus_TET2_rep2" "IAA_minus_TET2_rep3" "IAA_plus_TET2_rep1" "IAA_plus_TET2_rep2" "IAA_plus_TET2_rep3"

## compute and plot use rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_CTCF_rep123_merged_TvsC.bigwig IAA_plus_IP_CTCF_rep123_merged_TvsC.bigwig IAA_minus_IP_MBD3_rep123_merged_TvsC.bigwig IAA_plus_IP_MBD3_rep123_merged_TvsC.bigwig IAA_minus_IP_TET1_rep123_merged_TvsC.bigwig IAA_plus_IP_TET1_rep123_merged_TvsC.bigwig IAA_minus_IP_TET2_rep123_merged_TvsC.bigwig IAA_plus_IP_TET2_rep123_merged_TvsC.bigwig -R CTCF_rep12_narrowPeak.summit.union.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_ChIP_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_ChIP_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz -o bamCompare_matrix_of_IAA72h_ChIP_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_rep12_union_all_summits(97424)" --samplesLabel "IAA_minus_CTCF_rep123" "IAA_plus_CTCF_rep123" "IAA_minus_MBD3_rep123" "IAA_plus_MBD3_rep123" "IAA_minus_TET1_rep123" "IAA_plus_TET1_rep123" "IAA_minus_TET2_rep123" "IAA_plus_TET2_rep123"

## compute and plot use only CTCF rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_CTCF_rep123_merged_TvsC.bigwig IAA_plus_IP_CTCF_rep123_merged_TvsC.bigwig -R CTCF_rep12_narrowPeak.summit.union.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_ChIP_CTCF_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_ChIP_CTCF_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz -o bamCompare_matrix_of_IAA72h_ChIP_CTCF_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_rep12_union_all_summits(97424)" --samplesLabel "DMSO_CTCF_rep123_TvsC" "IAA_CTCF_rep123_TvsC"

## compute and plot use only MBD3 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_MBD3_rep123_merged_TvsC.bigwig IAA_plus_IP_MBD3_rep123_merged_TvsC.bigwig -R CTCF_rep12_narrowPeak.summit.union.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_ChIP_MBD3_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_ChIP_MBD3_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz -o bamCompare_matrix_of_IAA72h_ChIP_MBD3_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_rep12_union_all_summits(97424)" --samplesLabel "DMSO_MBD3_rep123_TvsC" "IAA_MBD3_rep123_TvsC"

# 调整颜色尺度 --zMin=-1.3 --zMax=1.9
plotHeatmap -m bamCompare_matrix_of_IAA72h_ChIP_MBD3_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz -o bamCompare_matrix_of_IAA72h_ChIP_MBD3_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.plot.changeColorScale.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --zMin=-1.3 --zMax=1.9 --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_rep12_union_all_summits(97424)" --samplesLabel "DMSO_MBD3_rep123_TvsC" "IAA_MBD3_rep123_TvsC"

```

## plot use union  source 3 class summit

```bash
# plot use CTCF rep12 union all summits
## compute and plot use rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_CTCF_rep1_TvsC.bigwig IAA_minus_IP_CTCF_rep2_TvsC.bigwig IAA_minus_IP_CTCF_rep3_TvsC.bigwig IAA_plus_IP_CTCF_rep1_TvsC.bigwig IAA_plus_IP_CTCF_rep2_TvsC.bigwig IAA_plus_IP_CTCF_rep3_TvsC.bigwig IAA_minus_IP_MBD3_rep1_TvsC.bigwig IAA_minus_IP_MBD3_rep2_TvsC.bigwig IAA_minus_IP_MBD3_rep3_TvsC.bigwig IAA_plus_IP_MBD3_rep1_TvsC.bigwig IAA_plus_IP_MBD3_rep2_TvsC.bigwig IAA_plus_IP_MBD3_rep3_TvsC.bigwig IAA_minus_IP_TET1_rep1_TvsC.bigwig IAA_minus_IP_TET1_rep2_TvsC.bigwig IAA_minus_IP_TET1_rep3_TvsC.bigwig IAA_plus_IP_TET1_rep1_TvsC.bigwig IAA_plus_IP_TET1_rep2_TvsC.bigwig IAA_plus_IP_TET1_rep3_TvsC.bigwig IAA_minus_IP_TET2_rep1_TvsC.bigwig IAA_minus_IP_TET2_rep2_TvsC.bigwig IAA_minus_IP_TET2_rep3_TvsC.bigwig IAA_plus_IP_TET2_rep1_TvsC.bigwig IAA_plus_IP_TET2_rep2_TvsC.bigwig IAA_plus_IP_TET2_rep3_TvsC.bigwig -R CTCF_rep12_union.rep1_only.summit.sorted.bed CTCF_rep12_union.rep12_common.summit.sorted.bed CTCF_rep12_union.rep2_only.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_ChIP_all_samples_for_plotHeatMap.public_CTCF_rep12_union_3source.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_ChIP_all_samples_for_plotHeatMap.public_CTCF_rep12_union_3source.gz -o bamCompare_matrix_of_IAA72h_ChIP_all_samples_for_plotHeatMap.public_CTCF_rep12_union_3source.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "rep1_only(47355)_summit" "rep12_common_summit(31737)" "rep2_only_summit(18332)" --samplesLabel "IAA_minus_CTCF_rep1" "IAA_minus_CTCF_rep2" "IAA_minus_CTCF_rep3" "IAA_plus_CTCF_rep1" "IAA_plus_CTCF_rep2" "IAA_plus_CTCF_rep3" "IAA_minus_MBD3_rep1" "IAA_minus_MBD3_rep2" "IAA_minus_MBD3_rep3" "IAA_plus_MBD3_rep1" "IAA_plus_MBD3_rep2" "IAA_plus_MBD3_rep3" "IAA_minus_TET1_rep1" "IAA_minus_TET1_rep2" "IAA_minus_TET1_rep3" "IAA_plus_TET1_rep1" "IAA_plus_TET1_rep2" "IAA_plus_TET1_rep3" "IAA_minus_TET2_rep1" "IAA_minus_TET2_rep2" "IAA_minus_TET2_rep3" "IAA_plus_TET2_rep1" "IAA_plus_TET2_rep2" "IAA_plus_TET2_rep3"

## compute and plot use rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_CTCF_rep123_merged_TvsC.bigwig IAA_plus_IP_CTCF_rep123_merged_TvsC.bigwig IAA_minus_IP_MBD3_rep123_merged_TvsC.bigwig IAA_plus_IP_MBD3_rep123_merged_TvsC.bigwig IAA_minus_IP_TET1_rep123_merged_TvsC.bigwig IAA_plus_IP_TET1_rep123_merged_TvsC.bigwig IAA_minus_IP_TET2_rep123_merged_TvsC.bigwig IAA_plus_IP_TET2_rep123_merged_TvsC.bigwig -R CTCF_rep12_union.rep1_only.summit.sorted.bed CTCF_rep12_union.rep12_common.summit.sorted.bed CTCF_rep12_union.rep2_only.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_ChIP_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_ChIP_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.gz -o bamCompare_matrix_of_IAA72h_ChIP_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "rep1_only(47355)_summit" "rep12_common_summit(31737)" "rep2_only_summit(18332)" --samplesLabel "IAA_minus_CTCF_rep123" "IAA_plus_CTCF_rep123" "IAA_minus_MBD3_rep123" "IAA_plus_MBD3_rep123" "IAA_minus_TET1_rep123" "IAA_plus_TET1_rep123" "IAA_minus_TET2_rep123" "IAA_plus_TET2_rep123"

## compute and plot use only CTCF rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_CTCF_rep123_merged_TvsC.bigwig IAA_plus_IP_CTCF_rep123_merged_TvsC.bigwig -R CTCF_rep12_union.rep1_only.summit.sorted.bed CTCF_rep12_union.rep12_common.summit.sorted.bed CTCF_rep12_union.rep2_only.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_ChIP_CTCF_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_ChIP_CTCF_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.gz -o bamCompare_matrix_of_IAA72h_ChIP_CTCF_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "rep1_only(47355)_summit" "rep12_common_summit(31737)" "rep2_only_summit(18332)" --samplesLabel "DMSO_CTCF_rep123_TvsC" "IAA_CTCF_rep123_TvsC"

## compute and plot use only MBD3 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_MBD3_rep123_merged_TvsC.bigwig IAA_plus_IP_MBD3_rep123_merged_TvsC.bigwig -R CTCF_rep12_union.rep1_only.summit.sorted.bed CTCF_rep12_union.rep12_common.summit.sorted.bed CTCF_rep12_union.rep2_only.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_ChIP_MBD3_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_ChIP_MBD3_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.gz -o bamCompare_matrix_of_IAA72h_ChIP_MBD3_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "rep1_only(47355)_summit" "rep12_common_summit(31737)" "rep2_only_summit(18332)" --samplesLabel "DMSO_MBD3_rep123_TvsC" "IAA_MBD3_rep123_TvsC"

```

## plot use union  4 class summit ( peak height quartile)

```bash
# plot use CTCF rep12 union all summits
## compute and plot use rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_CTCF_rep1_TvsC.bigwig IAA_minus_IP_CTCF_rep2_TvsC.bigwig IAA_minus_IP_CTCF_rep3_TvsC.bigwig IAA_plus_IP_CTCF_rep1_TvsC.bigwig IAA_plus_IP_CTCF_rep2_TvsC.bigwig IAA_plus_IP_CTCF_rep3_TvsC.bigwig IAA_minus_IP_MBD3_rep1_TvsC.bigwig IAA_minus_IP_MBD3_rep2_TvsC.bigwig IAA_minus_IP_MBD3_rep3_TvsC.bigwig IAA_plus_IP_MBD3_rep1_TvsC.bigwig IAA_plus_IP_MBD3_rep2_TvsC.bigwig IAA_plus_IP_MBD3_rep3_TvsC.bigwig IAA_minus_IP_TET1_rep1_TvsC.bigwig IAA_minus_IP_TET1_rep2_TvsC.bigwig IAA_minus_IP_TET1_rep3_TvsC.bigwig IAA_plus_IP_TET1_rep1_TvsC.bigwig IAA_plus_IP_TET1_rep2_TvsC.bigwig IAA_plus_IP_TET1_rep3_TvsC.bigwig IAA_minus_IP_TET2_rep1_TvsC.bigwig IAA_minus_IP_TET2_rep2_TvsC.bigwig IAA_minus_IP_TET2_rep3_TvsC.bigwig IAA_plus_IP_TET2_rep1_TvsC.bigwig IAA_plus_IP_TET2_rep2_TvsC.bigwig IAA_plus_IP_TET2_rep3_TvsC.bigwig -R CTCF_rep12_union.phq4.summit.sorted.bed CTCF_rep12_union.phq3.summit.sorted.bed CTCF_rep12_union.phq2.summit.sorted.bed CTCF_rep12_union.phq1.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_ChIP_all_samples_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_ChIP_all_samples_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz -o bamCompare_matrix_of_IAA72h_ChIP_all_samples_for_plotHeatMap.public_CTCF_rep12_union.phq1234.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "CTCF_rep12_union_phq4_summit(24356)" "CTCF_rep12_union_phq3_summit(24356)" "CTCF_rep12_union_phq2_summit(24356)" "CTCF_rep12_union_phq1_summit(24356)" --samplesLabel "IAA_minus_CTCF_rep1" "IAA_minus_CTCF_rep2" "IAA_minus_CTCF_rep3" "IAA_plus_CTCF_rep1" "IAA_plus_CTCF_rep2" "IAA_plus_CTCF_rep3" "IAA_minus_MBD3_rep1" "IAA_minus_MBD3_rep2" "IAA_minus_MBD3_rep3" "IAA_plus_MBD3_rep1" "IAA_plus_MBD3_rep2" "IAA_plus_MBD3_rep3" "IAA_minus_TET1_rep1" "IAA_minus_TET1_rep2" "IAA_minus_TET1_rep3" "IAA_plus_TET1_rep1" "IAA_plus_TET1_rep2" "IAA_plus_TET1_rep3" "IAA_minus_TET2_rep1" "IAA_minus_TET2_rep2" "IAA_minus_TET2_rep3" "IAA_plus_TET2_rep1" "IAA_plus_TET2_rep2" "IAA_plus_TET2_rep3"

## compute and plot use rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_CTCF_rep123_merged_TvsC.bigwig IAA_plus_IP_CTCF_rep123_merged_TvsC.bigwig IAA_minus_IP_MBD3_rep123_merged_TvsC.bigwig IAA_plus_IP_MBD3_rep123_merged_TvsC.bigwig IAA_minus_IP_TET1_rep123_merged_TvsC.bigwig IAA_plus_IP_TET1_rep123_merged_TvsC.bigwig IAA_minus_IP_TET2_rep123_merged_TvsC.bigwig IAA_plus_IP_TET2_rep123_merged_TvsC.bigwig -R CTCF_rep12_union.phq4.summit.sorted.bed CTCF_rep12_union.phq3.summit.sorted.bed CTCF_rep12_union.phq2.summit.sorted.bed CTCF_rep12_union.phq1.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_ChIP_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_ChIP_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz -o bamCompare_matrix_of_IAA72h_ChIP_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "CTCF_rep12_union_phq4_summit(24356)" "CTCF_rep12_union_phq3_summit(24356)" "CTCF_rep12_union_phq2_summit(24356)" "CTCF_rep12_union_phq1_summit(24356)" --samplesLabel "IAA_minus_CTCF_rep123" "IAA_plus_CTCF_rep123" "IAA_minus_MBD3_rep123" "IAA_plus_MBD3_rep123" "IAA_minus_TET1_rep123" "IAA_plus_TET1_rep123" "IAA_minus_TET2_rep123" "IAA_plus_TET2_rep123"

## compute and plot use only CTCF rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_CTCF_rep123_merged_TvsC.bigwig IAA_plus_IP_CTCF_rep123_merged_TvsC.bigwig -R CTCF_rep12_union.phq4.summit.sorted.bed CTCF_rep12_union.phq3.summit.sorted.bed CTCF_rep12_union.phq2.summit.sorted.bed CTCF_rep12_union.phq1.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_ChIP_CTCF_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_ChIP_CTCF_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz -o bamCompare_matrix_of_IAA72h_ChIP_CTCF_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "CTCF_rep12_union_phq4_summit(24356)" "CTCF_rep12_union_phq3_summit(24356)" "CTCF_rep12_union_phq2_summit(24356)" "CTCF_rep12_union_phq1_summit(24356)" --samplesLabel "DMSO_CTCF_rep123_TvsC" "IAA_CTCF_rep123_TvsC"

## compute and plot use only MBD3 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2022103101_WN_293T_IAA72h_ChIP_Seq_with_add_data/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S IAA_minus_IP_MBD3_rep123_merged_TvsC.bigwig IAA_plus_IP_MBD3_rep123_merged_TvsC.bigwig -R CTCF_rep12_union.phq4.summit.sorted.bed CTCF_rep12_union.phq3.summit.sorted.bed CTCF_rep12_union.phq2.summit.sorted.bed CTCF_rep12_union.phq1.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_IAA72h_ChIP_MBD3_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz
plotHeatmap -m bamCompare_matrix_of_IAA72h_ChIP_MBD3_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz -o bamCompare_matrix_of_IAA72h_ChIP_MBD3_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "CTCF_rep12_union_phq4_summit(24356)" "CTCF_rep12_union_phq3_summit(24356)" "CTCF_rep12_union_phq2_summit(24356)" "CTCF_rep12_union_phq1_summit(24356)" --samplesLabel "DMSO_MBD3_rep123_TvsC" "IAA_MBD3_rep123_TvsC"

```

