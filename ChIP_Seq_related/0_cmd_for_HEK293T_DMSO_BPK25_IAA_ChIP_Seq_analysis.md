

# QC

```bash
# 创建项目分析目录
cd /home/sunwenju/projects_on_940xa/
mkdir 2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged
# 构造并上传数据文件连接脚本，然后运行脚本
bash 0_make_and_rename_HEK293T_DMSO_IAA_BPK25_ChIP_Seq_raw_fq_gz_file.sh
ls -1
# 编写合并两批数据的脚本并运行 (使用cat比 zcat|gzip快很多)
bash 0_merge_DMSO_part1_and_part2_fq_file.sh
# 合并完成后删除连接文件
rm *DMSO*part[12]*.R[12].fq.gz
ls -1
# 
ls -1 HEK293T*DMSO*R1.fq.gz |sort |wc -l
# 15
ls -1 HEK293T*DMSO*R2.fq.gz |sort |wc -l
# 15
ls -1 HEK293T*BPK25*R1.fq.gz |sort |wc -l
# 15
ls -1 HEK293T*BPK25*R2.fq.gz |sort |wc -l
# 15
ls -1 HEK293T*IAA*R1.fq.gz |sort |wc -l
# 15
ls -1 HEK293T*IAA*R2.fq.gz |sort |wc -l
# 15

mkdir 1_all_HEK293T_ChIP_Seq_raw_data_fastqc_result
nohup fastqc -outdir 1_all_HEK293T_ChIP_Seq_raw_data_fastqc_result --threads 72 *.fq.gz >1_all_HEK293T_ChIP_Seq_raw_data_fastqc.runlog 2>&1 &
tail -f 1_all_HEK293T_ChIP_Seq_raw_data_fastqc.runlog

# 并用 multiqc 进行结果汇总
# 提取原始数据reads数量
cd 1_all_HEK293T_ChIP_Seq_raw_data_fastqc_result
multiqc ./*fastqc.zip --ignore *.html --no-megaqc-upload --interactive
cd ../

```



# mapping

```bash
# 使用如下命令构造比对脚本， 然后运行
for i in `ls -1 *.R1.fq.gz`;do echo 'bowtie2 -t -p 72 --no-unal -x /home/sunwenju/3_genomeData/human/hg19/bowtie2Index/hg19bt2idx -1 '$i' -2 '${i%.R1.fq.gz}.R2.fq.gz' -S '${i%.R1.fq.gz}.bowtie2.sam' --un-gz '${i%.R1.fq.gz}.notAlign.unpairedReads.fq.gz' --un-conc-gz '${i%.R1.fq.gz}.notAlign.pairedReads.fq.gz' >'${i%.R1.fq.gz}.bowtie2.runlog' 2>&1';done >2_all_HEK293T_ChIP_Seq_bowtie2_map.sh
nohup bash 2_all_HEK293T_ChIP_Seq_bowtie2_map.sh >2_all_HEK293T_ChIP_Seq_bowtie2_map.sh.runlog 2>&1 &
# tail -f 2_all_HEK293T_ChIP_Seq_bowtie2_map.sh.runlog
tail -f `ls -1 *bowtie2.runlog|tail -n 1`

# 使用 multiqc 对 bowtie2 比对情况进行汇总
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged
multiqc --no-megaqc-upload -m bowtie2 -o 2_all_HEK293T_ChIP_Seq_bowtie2_map_results_summary ./*.bowtie2.runlog

```



# Convert sam to bam

```bash
# 将比对sam结果转换为 按座位排序的 bam 并建立索引
# for i in `ls -1 *.bowtie2.sam|sort`; do echo $i; samtools view --threads 36 -bS $i | samtools sort --threads 36 -o ${i%.sam}.sorted.bam; samtools index ${i%.sam}.sorted.bam; done
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged
nohup bash 3_all_HEK293T_ChIP_Seq_data_bowtie2_sam_to_sorted_bam.sh >3_all_HEK293T_ChIP_Seq_data_bowtie2_sam_to_sorted_bam.sh.runlog 2>&1 &
tail -f 3_all_HEK293T_ChIP_Seq_data_bowtie2_sam_to_sorted_bam.sh.runlog
# rm *.bowtie2.sam

# 使用samtools flagstat 统计比对情况
# 每个bam文件不到一分钟
for i in `ls -1 *rep[123].bowtie2.sorted.bam`;do echo $i; samtools flagstat --threads 16 $i;done >3_all_HEK293T_ChIP_Seq_bowtie2_map_samtools_flagstat.txt

# 统计比对到各个染色体上的 总read数量 
for i in `ls -1 *rep[123].bowtie2.sorted.bam|sort`;do echo $i; samtools view $i |grep -v -e '^@' |cut -f 3 |sort |uniq -c |awk 'BEGIN{OFS="\t"}{if($0 !~ /_/){print $2,$1}}' |sort -k1,1 >${i%.sorted.bam}.chr_mapped_reads_stat.txt;done
# 统计比对到各个染色体上的 uniq read数量 
for i in `ls -1 *rep[123].bowtie2.sorted.bam|sort`;do echo $i; samtools view $i |grep -v -e '^@' |cut -f 3,4 |sort |uniq |cut -f 1 |sort |uniq -c |awk 'BEGIN{OFS="\t"}{if($0 !~ /_/){print $2,$1}}' |sort -k1,1 >${i%.sorted.bam}.chr_mapped_uniq_reads_stat.txt;done

# 比对到各个染色体上的 read数量 统计结果合并
paste `ls -1 *chr_mapped_reads_stat.txt|sort` |head
# 检测各个统计结果的染色体数量和顺序是否一致
paste `ls -1 *chr_mapped_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"}{if($1!=$3 || $1!=$5 || $1!=$7 || $1!=$9 || $1!=$11 || $1!=$13 || $1!=$15 || $1!=$17 || $1!=$19 || $1!=$21 || $1!=$23 || $1!=$25 || $1!=$27 || $1!=$29 || $1!=$31 || $1!=$33 || $1!=$35 || $1!=$37 || $1!=$39 || $1!=$41 || $1!=$43 || $1!=$45 || $1!=$47 || $1!=$49 || $1!=$51 || $1!=$53 || $1!=$55 || $1!=$57 || $1!=$59 || $1!=$61 || $1!=$63 || $1!=$65 || $1!=$67 || $1!=$69 || $1!=$71 || $1!=$73 || $1!=$75 || $1!=$77 || $1!=$79 || $1!=$81 || $1!=$83 || $1!=$85 || $1!=$87 || $1!=$89 || $1!=$91 || $1!=$93 || $1!=$95 || $1!=$97 || $1!=$99 || $1!=$101 || $1!=$103 || $1!=$105 || $1!=$107 || $1!=$109 || $1!=$111 || $1!=$113 || $1!=$115 || $1!=$117 || $1!=$119){print $1,$2}}' |head
# 染色体数量和顺序一致，则合并数据列
paste `ls -1 *chr_mapped_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"; print "chrom","BPK25_input_rep1","BPK25_input_rep2","BPK25_input_rep3","BPK25_IP_CTCF_rep1","BPK25_IP_CTCF_rep2","BPK25_IP_CTCF_rep3","BPK25_IP_MBD2_rep1","BPK25_IP_MBD2_rep2","BPK25_IP_MBD2_rep3","BPK25_IP_MBD3_rep1","BPK25_IP_MBD3_rep2","BPK25_IP_MBD3_rep3","BPK25_IP_TET1_rep1","BPK25_IP_TET1_rep2","BPK25_IP_TET1_rep3","DMSO_input_rep1","DMSO_input_rep2","DMSO_input_rep3","DMSO_IP_CTCF_rep1","DMSO_IP_CTCF_rep2","DMSO_IP_CTCF_rep3","DMSO_IP_MBD2_rep1","DMSO_IP_MBD2_rep2","DMSO_IP_MBD2_rep3","DMSO_IP_MBD3_rep1","DMSO_IP_MBD3_rep2","DMSO_IP_MBD3_rep3","DMSO_IP_TET1_rep1","DMSO_IP_TET1_rep2","DMSO_IP_TET1_rep3","IAA_input_rep1","IAA_input_rep2","IAA_input_rep3","IAA_IP_CTCF_rep1","IAA_IP_CTCF_rep2","IAA_IP_CTCF_rep3","IAA_IP_MBD2_rep1","IAA_IP_MBD2_rep2","IAA_IP_MBD2_rep3","IAA_IP_MBD3_rep1","IAA_IP_MBD3_rep2","IAA_IP_MBD3_rep3","IAA_IP_TET1_rep1","IAA_IP_TET1_rep2","IAA_IP_TET1_rep3"}{if($1==$3 && $1==$5 && $1==$7 && $1==$9 && $1==$11 && $1==$13 && $1==$15 && $1==$17 && $1==$19 && $1==$21 && $1==$23 && $1==$25 && $1==$27 && $1==$29 && $1==$31 && $1==$33 && $1==$35 && $1==$37 && $1==$39 && $1==$41 && $1==$43 && $1==$45 && $1==$47 && $1==$49 && $1==$51 && $1==$53 && $1==$55 && $1==$57 && $1==$59 && $1==$61 && $1==$63 && $1==$65 && $1==$67 && $1==$69 && $1==$71 && $1==$73 && $1==$75 && $1==$77 && $1==$79 && $1==$81 && $1==$83 && $1==$85 && $1==$87 && $1==$89){print $1, $2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36,$38,$40,$42,$44,$46,$48,$50,$52,$54,$56,$58,$60,$62,$64,$66,$68,$70,$72,$74,$76,$78,$80,$82,$84,$86,$88,$90}}' |head
paste `ls -1 *chr_mapped_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"; print "chrom","BPK25_input_rep1","BPK25_input_rep2","BPK25_input_rep3","BPK25_IP_CTCF_rep1","BPK25_IP_CTCF_rep2","BPK25_IP_CTCF_rep3","BPK25_IP_MBD2_rep1","BPK25_IP_MBD2_rep2","BPK25_IP_MBD2_rep3","BPK25_IP_MBD3_rep1","BPK25_IP_MBD3_rep2","BPK25_IP_MBD3_rep3","BPK25_IP_TET1_rep1","BPK25_IP_TET1_rep2","BPK25_IP_TET1_rep3","DMSO_input_rep1","DMSO_input_rep2","DMSO_input_rep3","DMSO_IP_CTCF_rep1","DMSO_IP_CTCF_rep2","DMSO_IP_CTCF_rep3","DMSO_IP_MBD2_rep1","DMSO_IP_MBD2_rep2","DMSO_IP_MBD2_rep3","DMSO_IP_MBD3_rep1","DMSO_IP_MBD3_rep2","DMSO_IP_MBD3_rep3","DMSO_IP_TET1_rep1","DMSO_IP_TET1_rep2","DMSO_IP_TET1_rep3","IAA_input_rep1","IAA_input_rep2","IAA_input_rep3","IAA_IP_CTCF_rep1","IAA_IP_CTCF_rep2","IAA_IP_CTCF_rep3","IAA_IP_MBD2_rep1","IAA_IP_MBD2_rep2","IAA_IP_MBD2_rep3","IAA_IP_MBD3_rep1","IAA_IP_MBD3_rep2","IAA_IP_MBD3_rep3","IAA_IP_TET1_rep1","IAA_IP_TET1_rep2","IAA_IP_TET1_rep3"}{if($1==$3 && $1==$5 && $1==$7 && $1==$9 && $1==$11 && $1==$13 && $1==$15 && $1==$17 && $1==$19 && $1==$21 && $1==$23 && $1==$25 && $1==$27 && $1==$29 && $1==$31 && $1==$33 && $1==$35 && $1==$37 && $1==$39 && $1==$41 && $1==$43 && $1==$45 && $1==$47 && $1==$49 && $1==$51 && $1==$53 && $1==$55 && $1==$57 && $1==$59 && $1==$61 && $1==$63 && $1==$65 && $1==$67 && $1==$69 && $1==$71 && $1==$73 && $1==$75 && $1==$77 && $1==$79 && $1==$81 && $1==$83 && $1==$85 && $1==$87 && $1==$89){print $1, $2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36,$38,$40,$42,$44,$46,$48,$50,$52,$54,$56,$58,$60,$62,$64,$66,$68,$70,$72,$74,$76,$78,$80,$82,$84,$86,$88,$90}}' >all_samples_chr_mapped_reads_stat.summary.txt

# 比对到各个染色体上的 uniq read数量 统计结果合并
paste `ls -1 *chr_mapped_uniq_reads_stat.txt|sort` |head
# 检测各个统计结果的染色体数量和顺序是否一致
paste `ls -1 *chr_mapped_uniq_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"}{if($1!=$3 || $1!=$5 || $1!=$7 || $1!=$9 || $1!=$11 || $1!=$13 || $1!=$15 || $1!=$17 || $1!=$19 || $1!=$21 || $1!=$23 || $1!=$25 || $1!=$27 || $1!=$29 || $1!=$31 || $1!=$33 || $1!=$35 || $1!=$37 || $1!=$39 || $1!=$41 || $1!=$43 || $1!=$45 || $1!=$47 || $1!=$49 || $1!=$51 || $1!=$53 || $1!=$55 || $1!=$57 || $1!=$59 || $1!=$61 || $1!=$63 || $1!=$65 || $1!=$67 || $1!=$69 || $1!=$71 || $1!=$73 || $1!=$75 || $1!=$77 || $1!=$79 || $1!=$81 || $1!=$83 || $1!=$85 || $1!=$87 || $1!=$89){print $1,$2}}' |head
# 染色体数量和顺序一致，则合并数据列
paste `ls -1 *chr_mapped_uniq_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"; print "chrom","BPK25_input_rep1","BPK25_input_rep2","BPK25_input_rep3","BPK25_IP_CTCF_rep1","BPK25_IP_CTCF_rep2","BPK25_IP_CTCF_rep3","BPK25_IP_MBD2_rep1","BPK25_IP_MBD2_rep2","BPK25_IP_MBD2_rep3","BPK25_IP_MBD3_rep1","BPK25_IP_MBD3_rep2","BPK25_IP_MBD3_rep3","BPK25_IP_TET1_rep1","BPK25_IP_TET1_rep2","BPK25_IP_TET1_rep3","DMSO_input_rep1","DMSO_input_rep2","DMSO_input_rep3","DMSO_IP_CTCF_rep1","DMSO_IP_CTCF_rep2","DMSO_IP_CTCF_rep3","DMSO_IP_MBD2_rep1","DMSO_IP_MBD2_rep2","DMSO_IP_MBD2_rep3","DMSO_IP_MBD3_rep1","DMSO_IP_MBD3_rep2","DMSO_IP_MBD3_rep3","DMSO_IP_TET1_rep1","DMSO_IP_TET1_rep2","DMSO_IP_TET1_rep3","IAA_input_rep1","IAA_input_rep2","IAA_input_rep3","IAA_IP_CTCF_rep1","IAA_IP_CTCF_rep2","IAA_IP_CTCF_rep3","IAA_IP_MBD2_rep1","IAA_IP_MBD2_rep2","IAA_IP_MBD2_rep3","IAA_IP_MBD3_rep1","IAA_IP_MBD3_rep2","IAA_IP_MBD3_rep3","IAA_IP_TET1_rep1","IAA_IP_TET1_rep2","IAA_IP_TET1_rep3"}{if($1==$3 && $1==$5 && $1==$7 && $1==$9 && $1==$11 && $1==$13 && $1==$15 && $1==$17 && $1==$19 && $1==$21 && $1==$23 && $1==$25 && $1==$27 && $1==$29 && $1==$31 && $1==$33 && $1==$35 && $1==$37 && $1==$39 && $1==$41 && $1==$43 && $1==$45 && $1==$47 && $1==$49 && $1==$51 && $1==$53 && $1==$55 && $1==$57 && $1==$59 && $1==$61 && $1==$63 && $1==$65 && $1==$67 && $1==$69 && $1==$71 && $1==$73 && $1==$75 && $1==$77 && $1==$79 && $1==$81 && $1==$83 && $1==$85 && $1==$87 && $1==$89){print $1, $2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36,$38,$40,$42,$44,$46,$48,$50,$52,$54,$56,$58,$60,$62,$64,$66,$68,$70,$72,$74,$76,$78,$80,$82,$84,$86,$88,$90}}' |head
paste `ls -1 *chr_mapped_uniq_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"; print "chrom","BPK25_input_rep1","BPK25_input_rep2","BPK25_input_rep3","BPK25_IP_CTCF_rep1","BPK25_IP_CTCF_rep2","BPK25_IP_CTCF_rep3","BPK25_IP_MBD2_rep1","BPK25_IP_MBD2_rep2","BPK25_IP_MBD2_rep3","BPK25_IP_MBD3_rep1","BPK25_IP_MBD3_rep2","BPK25_IP_MBD3_rep3","BPK25_IP_TET1_rep1","BPK25_IP_TET1_rep2","BPK25_IP_TET1_rep3","DMSO_input_rep1","DMSO_input_rep2","DMSO_input_rep3","DMSO_IP_CTCF_rep1","DMSO_IP_CTCF_rep2","DMSO_IP_CTCF_rep3","DMSO_IP_MBD2_rep1","DMSO_IP_MBD2_rep2","DMSO_IP_MBD2_rep3","DMSO_IP_MBD3_rep1","DMSO_IP_MBD3_rep2","DMSO_IP_MBD3_rep3","DMSO_IP_TET1_rep1","DMSO_IP_TET1_rep2","DMSO_IP_TET1_rep3","IAA_input_rep1","IAA_input_rep2","IAA_input_rep3","IAA_IP_CTCF_rep1","IAA_IP_CTCF_rep2","IAA_IP_CTCF_rep3","IAA_IP_MBD2_rep1","IAA_IP_MBD2_rep2","IAA_IP_MBD2_rep3","IAA_IP_MBD3_rep1","IAA_IP_MBD3_rep2","IAA_IP_MBD3_rep3","IAA_IP_TET1_rep1","IAA_IP_TET1_rep2","IAA_IP_TET1_rep3"}{if($1==$3 && $1==$5 && $1==$7 && $1==$9 && $1==$11 && $1==$13 && $1==$15 && $1==$17 && $1==$19 && $1==$21 && $1==$23 && $1==$25 && $1==$27 && $1==$29 && $1==$31 && $1==$33 && $1==$35 && $1==$37 && $1==$39 && $1==$41 && $1==$43 && $1==$45 && $1==$47 && $1==$49 && $1==$51 && $1==$53 && $1==$55 && $1==$57 && $1==$59 && $1==$61 && $1==$63 && $1==$65 && $1==$67 && $1==$69 && $1==$71 && $1==$73 && $1==$75 && $1==$77 && $1==$79 && $1==$81 && $1==$83 && $1==$85 && $1==$87 && $1==$89){print $1, $2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36,$38,$40,$42,$44,$46,$48,$50,$52,$54,$56,$58,$60,$62,$64,$66,$68,$70,$72,$74,$76,$78,$80,$82,$84,$86,$88,$90}}' >all_samples_chr_mapped_uniq_reads_stat.summary.txt

```

```bash
# 绘制reads 分布热图时要用到 rep123合并的bam文件
# 创建 bam merge 脚本 4_make_all_HEK293T_ChIP_Seq_rep123_merged_bam_and_index.sh 内容如下：
# # 合并各重复的 bam
# ## DMSO
# samtools merge --threads 72 HEK293T_DMSO_input_rep123.merged.sorted.bam HEK293T_DMSO_input_rep1.bowtie2.sorted.bam HEK293T_DMSO_input_rep2.bowtie2.sorted.bam HEK293T_DMSO_input_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_DMSO_IP_CTCF_rep123.merged.sorted.bam HEK293T_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_DMSO_IP_MBD2_rep123.merged.sorted.bam HEK293T_DMSO_IP_MBD2_rep1.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD2_rep2.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD2_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_DMSO_IP_MBD3_rep123.merged.sorted.bam HEK293T_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_DMSO_IP_TET1_rep123.merged.sorted.bam HEK293T_DMSO_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_DMSO_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_DMSO_IP_TET1_rep3.bowtie2.sorted.bam
# ## BPK25
# samtools merge --threads 72 HEK293T_BPK25_input_rep123.merged.sorted.bam HEK293T_BPK25_input_rep1.bowtie2.sorted.bam HEK293T_BPK25_input_rep2.bowtie2.sorted.bam HEK293T_BPK25_input_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_BPK25_IP_CTCF_rep123.merged.sorted.bam HEK293T_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_BPK25_IP_MBD2_rep123.merged.sorted.bam HEK293T_BPK25_IP_MBD2_rep1.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD2_rep2.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD2_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_BPK25_IP_MBD3_rep123.merged.sorted.bam HEK293T_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_BPK25_IP_TET1_rep123.merged.sorted.bam HEK293T_BPK25_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_BPK25_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_BPK25_IP_TET1_rep3.bowtie2.sorted.bam
# ## IAA
# samtools merge --threads 72 HEK293T_IAA_input_rep123.merged.sorted.bam HEK293T_IAA_input_rep1.bowtie2.sorted.bam HEK293T_IAA_input_rep2.bowtie2.sorted.bam HEK293T_IAA_input_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_IAA_IP_CTCF_rep123.merged.sorted.bam HEK293T_IAA_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_IAA_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_IAA_IP_CTCF_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_IAA_IP_MBD2_rep123.merged.sorted.bam HEK293T_IAA_IP_MBD2_rep1.bowtie2.sorted.bam HEK293T_IAA_IP_MBD2_rep2.bowtie2.sorted.bam HEK293T_IAA_IP_MBD2_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_IAA_IP_MBD3_rep123.merged.sorted.bam HEK293T_IAA_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_IAA_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_IAA_IP_MBD3_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 HEK293T_IAA_IP_TET1_rep123.merged.sorted.bam HEK293T_IAA_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_IAA_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_IAA_IP_TET1_rep3.bowtie2.sorted.bam
# # 为所有 merged bam 建立索引
# for i in `ls -1 *rep123.merged.sorted.bam`; do echo $i; samtools index -@ 72 $i;done

nohup bash 4_make_all_HEK293T_ChIP_Seq_rep123_merged_bam_and_index.sh >4_make_all_HEK293T_ChIP_Seq_rep123_merged_bam_and_index.sh.runlog 2>&1 &
tail -f 4_make_all_HEK293T_ChIP_Seq_rep123_merged_bam_and_index.sh.runlog

```



# Detection of enrichment and repeatability

```bash
# 参考： https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html
# 可以判断实验的成败、富集程度的好坏等
# DMSO all samples
plotFingerprint -b HEK293T_DMSO_input_rep1.bowtie2.sorted.bam HEK293T_DMSO_input_rep2.bowtie2.sorted.bam HEK293T_DMSO_input_rep3.bowtie2.sorted.bam HEK293T_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD2_rep1.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD2_rep2.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD2_rep3.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam HEK293T_DMSO_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_DMSO_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_DMSO_IP_TET1_rep3.bowtie2.sorted.bam --labels HEK293T_DMSO_input_rep1 HEK293T_DMSO_input_rep2 HEK293T_DMSO_input_rep3 HEK293T_DMSO_IP_CTCF_rep1 HEK293T_DMSO_IP_CTCF_rep2 HEK293T_DMSO_IP_CTCF_rep3 HEK293T_DMSO_IP_MBD2_rep1 HEK293T_DMSO_IP_MBD2_rep2 HEK293T_DMSO_IP_MBD2_rep3 HEK293T_DMSO_IP_MBD3_rep1 HEK293T_DMSO_IP_MBD3_rep2 HEK293T_DMSO_IP_MBD3_rep3 HEK293T_DMSO_IP_TET1_rep1 HEK293T_DMSO_IP_TET1_rep2 HEK293T_DMSO_IP_TET1_rep3 -o plotFingerprint_for_DMSO_all_samples.pdf --plotTitle "fingerprint_plot_for_DMSO_all_samples" --extendReads --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max"
# BPK25 all samples
plotFingerprint -b HEK293T_BPK25_input_rep1.bowtie2.sorted.bam HEK293T_BPK25_input_rep2.bowtie2.sorted.bam HEK293T_BPK25_input_rep3.bowtie2.sorted.bam HEK293T_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD2_rep1.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD2_rep2.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD2_rep3.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam HEK293T_BPK25_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_BPK25_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_BPK25_IP_TET1_rep3.bowtie2.sorted.bam --labels HEK293T_BPK25_input_rep1 HEK293T_BPK25_input_rep2 HEK293T_BPK25_input_rep3 HEK293T_BPK25_IP_CTCF_rep1 HEK293T_BPK25_IP_CTCF_rep2 HEK293T_BPK25_IP_CTCF_rep3 HEK293T_BPK25_IP_MBD2_rep1 HEK293T_BPK25_IP_MBD2_rep2 HEK293T_BPK25_IP_MBD2_rep3 HEK293T_BPK25_IP_MBD3_rep1 HEK293T_BPK25_IP_MBD3_rep2 HEK293T_BPK25_IP_MBD3_rep3 HEK293T_BPK25_IP_TET1_rep1 HEK293T_BPK25_IP_TET1_rep2 HEK293T_BPK25_IP_TET1_rep3 -o plotFingerprint_for_BPK25_all_samples.pdf --plotTitle "fingerprint_plot_for_BPK25_all_samples" --extendReads --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max"
# IAA all samples
plotFingerprint -b HEK293T_IAA_input_rep1.bowtie2.sorted.bam HEK293T_IAA_input_rep2.bowtie2.sorted.bam HEK293T_IAA_input_rep3.bowtie2.sorted.bam HEK293T_IAA_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_IAA_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_IAA_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_IAA_IP_MBD2_rep1.bowtie2.sorted.bam HEK293T_IAA_IP_MBD2_rep2.bowtie2.sorted.bam HEK293T_IAA_IP_MBD2_rep3.bowtie2.sorted.bam HEK293T_IAA_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_IAA_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_IAA_IP_MBD3_rep3.bowtie2.sorted.bam HEK293T_IAA_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_IAA_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_IAA_IP_TET1_rep3.bowtie2.sorted.bam --labels HEK293T_IAA_input_rep1 HEK293T_IAA_input_rep2 HEK293T_IAA_input_rep3 HEK293T_IAA_IP_CTCF_rep1 HEK293T_IAA_IP_CTCF_rep2 HEK293T_IAA_IP_CTCF_rep3 HEK293T_IAA_IP_MBD2_rep1 HEK293T_IAA_IP_MBD2_rep2 HEK293T_IAA_IP_MBD2_rep3 HEK293T_IAA_IP_MBD3_rep1 HEK293T_IAA_IP_MBD3_rep2 HEK293T_IAA_IP_MBD3_rep3 HEK293T_IAA_IP_TET1_rep1 HEK293T_IAA_IP_TET1_rep2 HEK293T_IAA_IP_TET1_rep3 -o plotFingerprint_for_IAA_all_samples.pdf --plotTitle "fingerprint_plot_for_IAA_all_samples" --extendReads --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max"

```

```bash
# CTCF (DMSO/BPK25/IAA)
plotFingerprint -b HEK293T_DMSO_input_rep1.bowtie2.sorted.bam HEK293T_DMSO_input_rep2.bowtie2.sorted.bam HEK293T_DMSO_input_rep3.bowtie2.sorted.bam HEK293T_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_BPK25_input_rep1.bowtie2.sorted.bam HEK293T_BPK25_input_rep2.bowtie2.sorted.bam HEK293T_BPK25_input_rep3.bowtie2.sorted.bam HEK293T_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_IAA_input_rep1.bowtie2.sorted.bam HEK293T_IAA_input_rep2.bowtie2.sorted.bam HEK293T_IAA_input_rep3.bowtie2.sorted.bam HEK293T_IAA_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_IAA_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_IAA_IP_CTCF_rep3.bowtie2.sorted.bam --labels HEK293T_DMSO_input_rep1 HEK293T_DMSO_input_rep2 HEK293T_DMSO_input_rep3 HEK293T_DMSO_IP_CTCF_rep1 HEK293T_DMSO_IP_CTCF_rep2 HEK293T_DMSO_IP_CTCF_rep3 HEK293T_BPK25_input_rep1 HEK293T_BPK25_input_rep2 HEK293T_BPK25_input_rep3 HEK293T_BPK25_IP_CTCF_rep1 HEK293T_BPK25_IP_CTCF_rep2 HEK293T_BPK25_IP_CTCF_rep3 HEK293T_IAA_input_rep1 HEK293T_IAA_input_rep2 HEK293T_IAA_input_rep3 HEK293T_IAA_IP_CTCF_rep1 HEK293T_IAA_IP_CTCF_rep2 HEK293T_IAA_IP_CTCF_rep3 -o plotFingerprint_for_all_treatment_CTCF_ChIP_Seq_samples.pdf --plotTitle "fingerprint_plot_for_all_treatment_CTCF_ChIP_Seq_samples" --extendReads --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max"

# MBD2 (DMSO/BPK25/IAA)
plotFingerprint -b HEK293T_DMSO_input_rep1.bowtie2.sorted.bam HEK293T_DMSO_input_rep2.bowtie2.sorted.bam HEK293T_DMSO_input_rep3.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD2_rep1.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD2_rep2.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD2_rep3.bowtie2.sorted.bam HEK293T_BPK25_input_rep1.bowtie2.sorted.bam HEK293T_BPK25_input_rep2.bowtie2.sorted.bam HEK293T_BPK25_input_rep3.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD2_rep1.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD2_rep2.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD2_rep3.bowtie2.sorted.bam HEK293T_IAA_input_rep1.bowtie2.sorted.bam HEK293T_IAA_input_rep2.bowtie2.sorted.bam HEK293T_IAA_input_rep3.bowtie2.sorted.bam HEK293T_IAA_IP_MBD2_rep1.bowtie2.sorted.bam HEK293T_IAA_IP_MBD2_rep2.bowtie2.sorted.bam HEK293T_IAA_IP_MBD2_rep3.bowtie2.sorted.bam --labels HEK293T_DMSO_input_rep1 HEK293T_DMSO_input_rep2 HEK293T_DMSO_input_rep3 HEK293T_DMSO_IP_MBD2_rep1 HEK293T_DMSO_IP_MBD2_rep2 HEK293T_DMSO_IP_MBD2_rep3 HEK293T_BPK25_input_rep1 HEK293T_BPK25_input_rep2 HEK293T_BPK25_input_rep3 HEK293T_BPK25_IP_MBD2_rep1 HEK293T_BPK25_IP_MBD2_rep2 HEK293T_BPK25_IP_MBD2_rep3 HEK293T_IAA_input_rep1 HEK293T_IAA_input_rep2 HEK293T_IAA_input_rep3 HEK293T_IAA_IP_MBD2_rep1 HEK293T_IAA_IP_MBD2_rep2 HEK293T_IAA_IP_MBD2_rep3 -o plotFingerprint_for_all_treatment_MBD2_ChIP_Seq_samples.pdf --plotTitle "fingerprint_plot_for_all_treatment_MBD2_ChIP_Seq_samples" --extendReads --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max"

# MBD3 (DMSO/BPK25/IAA)
plotFingerprint -b HEK293T_DMSO_input_rep1.bowtie2.sorted.bam HEK293T_DMSO_input_rep2.bowtie2.sorted.bam HEK293T_DMSO_input_rep3.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam HEK293T_BPK25_input_rep1.bowtie2.sorted.bam HEK293T_BPK25_input_rep2.bowtie2.sorted.bam HEK293T_BPK25_input_rep3.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam HEK293T_IAA_input_rep1.bowtie2.sorted.bam HEK293T_IAA_input_rep2.bowtie2.sorted.bam HEK293T_IAA_input_rep3.bowtie2.sorted.bam HEK293T_IAA_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_IAA_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_IAA_IP_MBD3_rep3.bowtie2.sorted.bam --labels HEK293T_DMSO_input_rep1 HEK293T_DMSO_input_rep2 HEK293T_DMSO_input_rep3 HEK293T_DMSO_IP_MBD3_rep1 HEK293T_DMSO_IP_MBD3_rep2 HEK293T_DMSO_IP_MBD3_rep3 HEK293T_BPK25_input_rep1 HEK293T_BPK25_input_rep2 HEK293T_BPK25_input_rep3 HEK293T_BPK25_IP_MBD3_rep1 HEK293T_BPK25_IP_MBD3_rep2 HEK293T_BPK25_IP_MBD3_rep3 HEK293T_IAA_input_rep1 HEK293T_IAA_input_rep2 HEK293T_IAA_input_rep3 HEK293T_IAA_IP_MBD3_rep1 HEK293T_IAA_IP_MBD3_rep2 HEK293T_IAA_IP_MBD3_rep3 -o plotFingerprint_for_all_treatment_MBD3_ChIP_Seq_samples.pdf --plotTitle "fingerprint_plot_for_all_treatment_MBD3_ChIP_Seq_samples" --extendReads --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max"

# TET1 (DMSO/BPK25/IAA)
plotFingerprint -b HEK293T_DMSO_input_rep1.bowtie2.sorted.bam HEK293T_DMSO_input_rep2.bowtie2.sorted.bam HEK293T_DMSO_input_rep3.bowtie2.sorted.bam HEK293T_DMSO_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_DMSO_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_DMSO_IP_TET1_rep3.bowtie2.sorted.bam HEK293T_BPK25_input_rep1.bowtie2.sorted.bam HEK293T_BPK25_input_rep2.bowtie2.sorted.bam HEK293T_BPK25_input_rep3.bowtie2.sorted.bam HEK293T_BPK25_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_BPK25_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_BPK25_IP_TET1_rep3.bowtie2.sorted.bam HEK293T_IAA_input_rep1.bowtie2.sorted.bam HEK293T_IAA_input_rep2.bowtie2.sorted.bam HEK293T_IAA_input_rep3.bowtie2.sorted.bam HEK293T_IAA_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_IAA_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_IAA_IP_TET1_rep3.bowtie2.sorted.bam --labels HEK293T_DMSO_input_rep1 HEK293T_DMSO_input_rep2 HEK293T_DMSO_input_rep3 HEK293T_DMSO_IP_TET1_rep1 HEK293T_DMSO_IP_TET1_rep2 HEK293T_DMSO_IP_TET1_rep3 HEK293T_BPK25_input_rep1 HEK293T_BPK25_input_rep2 HEK293T_BPK25_input_rep3 HEK293T_BPK25_IP_TET1_rep1 HEK293T_BPK25_IP_TET1_rep2 HEK293T_BPK25_IP_TET1_rep3 HEK293T_IAA_input_rep1 HEK293T_IAA_input_rep2 HEK293T_IAA_input_rep3 HEK293T_IAA_IP_TET1_rep1 HEK293T_IAA_IP_TET1_rep2 HEK293T_IAA_IP_TET1_rep3 -o plotFingerprint_for_all_treatment_TET1_ChIP_Seq_samples.pdf --plotTitle "fingerprint_plot_for_all_treatment_TET1_ChIP_Seq_samples" --extendReads --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max"

```

```bash
# 使用 deeptools 分析 样本相似度
# 先对全基因组分bin统计 uniq reads数量
# use DMSO
multiBamSummary bins --bamfiles HEK293T_DMSO_input_rep1.bowtie2.sorted.bam HEK293T_DMSO_input_rep2.bowtie2.sorted.bam HEK293T_DMSO_input_rep3.bowtie2.sorted.bam HEK293T_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD2_rep1.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD2_rep2.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD2_rep3.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam HEK293T_DMSO_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_DMSO_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_DMSO_IP_TET1_rep3.bowtie2.sorted.bam HEK293T_BPK25_input_rep1.bowtie2.sorted.bam HEK293T_BPK25_input_rep2.bowtie2.sorted.bam HEK293T_BPK25_input_rep3.bowtie2.sorted.bam HEK293T_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD2_rep1.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD2_rep2.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD2_rep3.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam HEK293T_BPK25_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_BPK25_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_BPK25_IP_TET1_rep3.bowtie2.sorted.bam HEK293T_IAA_input_rep1.bowtie2.sorted.bam HEK293T_IAA_input_rep2.bowtie2.sorted.bam HEK293T_IAA_input_rep3.bowtie2.sorted.bam HEK293T_IAA_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_IAA_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_IAA_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_IAA_IP_MBD2_rep1.bowtie2.sorted.bam HEK293T_IAA_IP_MBD2_rep2.bowtie2.sorted.bam HEK293T_IAA_IP_MBD2_rep3.bowtie2.sorted.bam HEK293T_IAA_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_IAA_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_IAA_IP_MBD3_rep3.bowtie2.sorted.bam HEK293T_IAA_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_IAA_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_IAA_IP_TET1_rep3.bowtie2.sorted.bam -out multiBamSummary_results_for_all_samples_bam_files.npz --labels HEK293T_DMSO_input_rep1 HEK293T_DMSO_input_rep2 HEK293T_DMSO_input_rep3 HEK293T_DMSO_IP_CTCF_rep1 HEK293T_DMSO_IP_CTCF_rep2 HEK293T_DMSO_IP_CTCF_rep3 HEK293T_DMSO_IP_MBD2_rep1 HEK293T_DMSO_IP_MBD2_rep2 HEK293T_DMSO_IP_MBD2_rep3 HEK293T_DMSO_IP_MBD3_rep1 HEK293T_DMSO_IP_MBD3_rep2 HEK293T_DMSO_IP_MBD3_rep3 HEK293T_DMSO_IP_TET1_rep1 HEK293T_DMSO_IP_TET1_rep2 HEK293T_DMSO_IP_TET1_rep3 HEK293T_BPK25_input_rep1 HEK293T_BPK25_input_rep2 HEK293T_BPK25_input_rep3 HEK293T_BPK25_IP_CTCF_rep1 HEK293T_BPK25_IP_CTCF_rep2 HEK293T_BPK25_IP_CTCF_rep3 HEK293T_BPK25_IP_MBD2_rep1 HEK293T_BPK25_IP_MBD2_rep2 HEK293T_BPK25_IP_MBD2_rep3 HEK293T_BPK25_IP_MBD3_rep1 HEK293T_BPK25_IP_MBD3_rep2 HEK293T_BPK25_IP_MBD3_rep3 HEK293T_BPK25_IP_TET1_rep1 HEK293T_BPK25_IP_TET1_rep2 HEK293T_BPK25_IP_TET1_rep3 HEK293T_IAA_input_rep1 HEK293T_IAA_input_rep2 HEK293T_IAA_input_rep3 HEK293T_IAA_IP_CTCF_rep1 HEK293T_IAA_IP_CTCF_rep2 HEK293T_IAA_IP_CTCF_rep3 HEK293T_IAA_IP_MBD2_rep1 HEK293T_IAA_IP_MBD2_rep2 HEK293T_IAA_IP_MBD2_rep3 HEK293T_IAA_IP_MBD3_rep1 HEK293T_IAA_IP_MBD3_rep2 HEK293T_IAA_IP_MBD3_rep3 HEK293T_IAA_IP_TET1_rep1 HEK293T_IAA_IP_TET1_rep2 HEK293T_IAA_IP_TET1_rep3 --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# 用 plotCorrelation 和 plotPCA 分析样本的相似度
plotCorrelation --corData multiBamSummary_results_for_all_samples_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_all_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_all_samples.txt --colorMap coolwarm --plotNumbers
plotPCA --corData multiBamSummary_results_for_all_samples_bam_files.npz --plotFile plotPCA_plot_for_all_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_all_samples.txt
plotPCA --rowCenter --corData multiBamSummary_results_for_all_samples_bam_files.npz --plotFile plotPCA_plot_for_all_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_all_samples.rowCenter.txt

```

```bash
# CTCF ChIP_Seq 样本的相似度
multiBamSummary bins --bamfiles HEK293T_DMSO_input_rep1.bowtie2.sorted.bam HEK293T_DMSO_input_rep2.bowtie2.sorted.bam HEK293T_DMSO_input_rep3.bowtie2.sorted.bam HEK293T_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_BPK25_input_rep1.bowtie2.sorted.bam HEK293T_BPK25_input_rep2.bowtie2.sorted.bam HEK293T_BPK25_input_rep3.bowtie2.sorted.bam HEK293T_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam HEK293T_IAA_input_rep1.bowtie2.sorted.bam HEK293T_IAA_input_rep2.bowtie2.sorted.bam HEK293T_IAA_input_rep3.bowtie2.sorted.bam HEK293T_IAA_IP_CTCF_rep1.bowtie2.sorted.bam HEK293T_IAA_IP_CTCF_rep2.bowtie2.sorted.bam HEK293T_IAA_IP_CTCF_rep3.bowtie2.sorted.bam -out multiBamSummary_results_for_all_treatment_CTCF_ChIP_Seq_bam_files.npz --labels HEK293T_DMSO_input_rep1 HEK293T_DMSO_input_rep2 HEK293T_DMSO_input_rep3 HEK293T_DMSO_IP_CTCF_rep1 HEK293T_DMSO_IP_CTCF_rep2 HEK293T_DMSO_IP_CTCF_rep3 HEK293T_BPK25_input_rep1 HEK293T_BPK25_input_rep2 HEK293T_BPK25_input_rep3 HEK293T_BPK25_IP_CTCF_rep1 HEK293T_BPK25_IP_CTCF_rep2 HEK293T_BPK25_IP_CTCF_rep3 HEK293T_IAA_input_rep1 HEK293T_IAA_input_rep2 HEK293T_IAA_input_rep3 HEK293T_IAA_IP_CTCF_rep1 HEK293T_IAA_IP_CTCF_rep2 HEK293T_IAA_IP_CTCF_rep3 --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
plotCorrelation --corData multiBamSummary_results_for_all_treatment_CTCF_ChIP_Seq_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_all_treatment_CTCF_ChIP_Seq_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_all_treatment_CTCF_ChIP_Seq_samples.txt --colorMap coolwarm --plotNumbers
plotPCA --corData multiBamSummary_results_for_all_treatment_CTCF_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_all_treatment_CTCF_ChIP_Seq_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_all_treatment_CTCF_ChIP_Seq_samples.txt
plotPCA --rowCenter --corData multiBamSummary_results_for_all_treatment_CTCF_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_all_treatment_CTCF_ChIP_Seq_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_all_treatment_CTCF_ChIP_Seq_samples.rowCenter.txt

# MBD3 ChIP_Seq 样本的相似度
multiBamSummary bins --bamfiles HEK293T_DMSO_input_rep1.bowtie2.sorted.bam HEK293T_DMSO_input_rep2.bowtie2.sorted.bam HEK293T_DMSO_input_rep3.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam HEK293T_BPK25_input_rep1.bowtie2.sorted.bam HEK293T_BPK25_input_rep2.bowtie2.sorted.bam HEK293T_BPK25_input_rep3.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam HEK293T_IAA_input_rep1.bowtie2.sorted.bam HEK293T_IAA_input_rep2.bowtie2.sorted.bam HEK293T_IAA_input_rep3.bowtie2.sorted.bam HEK293T_IAA_IP_MBD3_rep1.bowtie2.sorted.bam HEK293T_IAA_IP_MBD3_rep2.bowtie2.sorted.bam HEK293T_IAA_IP_MBD3_rep3.bowtie2.sorted.bam -out multiBamSummary_results_for_all_treatment_MBD3_ChIP_Seq_bam_files.npz --labels HEK293T_DMSO_input_rep1 HEK293T_DMSO_input_rep2 HEK293T_DMSO_input_rep3 HEK293T_DMSO_IP_MBD3_rep1 HEK293T_DMSO_IP_MBD3_rep2 HEK293T_DMSO_IP_MBD3_rep3 HEK293T_BPK25_input_rep1 HEK293T_BPK25_input_rep2 HEK293T_BPK25_input_rep3 HEK293T_BPK25_IP_MBD3_rep1 HEK293T_BPK25_IP_MBD3_rep2 HEK293T_BPK25_IP_MBD3_rep3 HEK293T_IAA_input_rep1 HEK293T_IAA_input_rep2 HEK293T_IAA_input_rep3 HEK293T_IAA_IP_MBD3_rep1 HEK293T_IAA_IP_MBD3_rep2 HEK293T_IAA_IP_MBD3_rep3 --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
plotCorrelation --corData multiBamSummary_results_for_all_treatment_MBD3_ChIP_Seq_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_all_treatment_MBD3_ChIP_Seq_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_all_treatment_MBD3_ChIP_Seq_samples.txt --colorMap coolwarm --plotNumbers
plotPCA --corData multiBamSummary_results_for_all_treatment_MBD3_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_all_treatment_MBD3_ChIP_Seq_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_all_treatment_MBD3_ChIP_Seq_samples.txt
plotPCA --rowCenter --corData multiBamSummary_results_for_all_treatment_MBD3_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_all_treatment_MBD3_ChIP_Seq_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_all_treatment_MBD3_ChIP_Seq_samples.rowCenter.txt

# MBD2 ChIP_Seq 样本的相似度
multiBamSummary bins --bamfiles HEK293T_DMSO_input_rep1.bowtie2.sorted.bam HEK293T_DMSO_input_rep2.bowtie2.sorted.bam HEK293T_DMSO_input_rep3.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD2_rep1.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD2_rep2.bowtie2.sorted.bam HEK293T_DMSO_IP_MBD2_rep3.bowtie2.sorted.bam HEK293T_BPK25_input_rep1.bowtie2.sorted.bam HEK293T_BPK25_input_rep2.bowtie2.sorted.bam HEK293T_BPK25_input_rep3.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD2_rep1.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD2_rep2.bowtie2.sorted.bam HEK293T_BPK25_IP_MBD2_rep3.bowtie2.sorted.bam HEK293T_IAA_input_rep1.bowtie2.sorted.bam HEK293T_IAA_input_rep2.bowtie2.sorted.bam HEK293T_IAA_input_rep3.bowtie2.sorted.bam HEK293T_IAA_IP_MBD2_rep1.bowtie2.sorted.bam HEK293T_IAA_IP_MBD2_rep2.bowtie2.sorted.bam HEK293T_IAA_IP_MBD2_rep3.bowtie2.sorted.bam -out multiBamSummary_results_for_all_treatment_MBD2_ChIP_Seq_bam_files.npz --labels HEK293T_DMSO_input_rep1 HEK293T_DMSO_input_rep2 HEK293T_DMSO_input_rep3 HEK293T_DMSO_IP_MBD2_rep1 HEK293T_DMSO_IP_MBD2_rep2 HEK293T_DMSO_IP_MBD2_rep3 HEK293T_BPK25_input_rep1 HEK293T_BPK25_input_rep2 HEK293T_BPK25_input_rep3 HEK293T_BPK25_IP_MBD2_rep1 HEK293T_BPK25_IP_MBD2_rep2 HEK293T_BPK25_IP_MBD2_rep3 HEK293T_IAA_input_rep1 HEK293T_IAA_input_rep2 HEK293T_IAA_input_rep3 HEK293T_IAA_IP_MBD2_rep1 HEK293T_IAA_IP_MBD2_rep2 HEK293T_IAA_IP_MBD2_rep3 --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
plotCorrelation --corData multiBamSummary_results_for_all_treatment_MBD2_ChIP_Seq_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_all_treatment_MBD2_ChIP_Seq_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_all_treatment_MBD2_ChIP_Seq_samples.txt --colorMap coolwarm --plotNumbers
plotPCA --corData multiBamSummary_results_for_all_treatment_MBD2_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_all_treatment_MBD2_ChIP_Seq_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_all_treatment_MBD2_ChIP_Seq_samples.txt
plotPCA --rowCenter --corData multiBamSummary_results_for_all_treatment_MBD2_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_all_treatment_MBD2_ChIP_Seq_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_all_treatment_MBD2_ChIP_Seq_samples.rowCenter.txt

# TET1 ChIP_Seq 样本的相似度
multiBamSummary bins --bamfiles HEK293T_DMSO_input_rep1.bowtie2.sorted.bam HEK293T_DMSO_input_rep2.bowtie2.sorted.bam HEK293T_DMSO_input_rep3.bowtie2.sorted.bam HEK293T_DMSO_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_DMSO_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_DMSO_IP_TET1_rep3.bowtie2.sorted.bam HEK293T_BPK25_input_rep1.bowtie2.sorted.bam HEK293T_BPK25_input_rep2.bowtie2.sorted.bam HEK293T_BPK25_input_rep3.bowtie2.sorted.bam HEK293T_BPK25_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_BPK25_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_BPK25_IP_TET1_rep3.bowtie2.sorted.bam HEK293T_IAA_input_rep1.bowtie2.sorted.bam HEK293T_IAA_input_rep2.bowtie2.sorted.bam HEK293T_IAA_input_rep3.bowtie2.sorted.bam HEK293T_IAA_IP_TET1_rep1.bowtie2.sorted.bam HEK293T_IAA_IP_TET1_rep2.bowtie2.sorted.bam HEK293T_IAA_IP_TET1_rep3.bowtie2.sorted.bam -out multiBamSummary_results_for_all_treatment_TET1_ChIP_Seq_bam_files.npz --labels HEK293T_DMSO_input_rep1 HEK293T_DMSO_input_rep2 HEK293T_DMSO_input_rep3 HEK293T_DMSO_IP_TET1_rep1 HEK293T_DMSO_IP_TET1_rep2 HEK293T_DMSO_IP_TET1_rep3 HEK293T_BPK25_input_rep1 HEK293T_BPK25_input_rep2 HEK293T_BPK25_input_rep3 HEK293T_BPK25_IP_TET1_rep1 HEK293T_BPK25_IP_TET1_rep2 HEK293T_BPK25_IP_TET1_rep3 HEK293T_IAA_input_rep1 HEK293T_IAA_input_rep2 HEK293T_IAA_input_rep3 HEK293T_IAA_IP_TET1_rep1 HEK293T_IAA_IP_TET1_rep2 HEK293T_IAA_IP_TET1_rep3 --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
plotCorrelation --corData multiBamSummary_results_for_all_treatment_TET1_ChIP_Seq_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_all_treatment_TET1_ChIP_Seq_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_all_treatment_TET1_ChIP_Seq_samples.txt --colorMap coolwarm --plotNumbers
plotPCA --corData multiBamSummary_results_for_all_treatment_TET1_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_all_treatment_TET1_ChIP_Seq_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_all_treatment_TET1_ChIP_Seq_samples.txt
plotPCA --rowCenter --corData multiBamSummary_results_for_all_treatment_TET1_ChIP_Seq_bam_files.npz --plotFile plotPCA_plot_for_all_treatment_TET1_ChIP_Seq_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_all_treatment_TET1_ChIP_Seq_samples.rowCenter.txt

```



# macs2 peak calling

## peak calling for each repetition separately 

### DMSO ChIP-Seq

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged
mkdir macs2_peak_calling_1_for_DMSO_each_repeat
cd macs2_peak_calling_1_for_DMSO_each_repeat
# 上传 bam 连接构造脚本并运行
bash 5_make_DMSO_all_ChIP_Seq_sorted_bam_file_links.sh
ls -1 *.bam

# DMSO_IP_CTCF
nohup macs2 callpeak -t HEK293T_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam -c HEK293T_DMSO_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir DMSO_IP_CTCF_rep1_TvsC -n DMSO_IP_CTCF_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_DMSO_IP_CTCF_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam -c HEK293T_DMSO_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir DMSO_IP_CTCF_rep2_TvsC -n DMSO_IP_CTCF_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_DMSO_IP_CTCF_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam -c HEK293T_DMSO_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir DMSO_IP_CTCF_rep3_TvsC -n DMSO_IP_CTCF_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_DMSO_IP_CTCF_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_DMSO_IP_CTCF_rep1_TvsC.runlog
tail -f macs2_callpeak_for_DMSO_IP_CTCF_rep2_TvsC.runlog
tail -f macs2_callpeak_for_DMSO_IP_CTCF_rep3_TvsC.runlog

# DMSO_IP_MBD2
nohup macs2 callpeak -t HEK293T_DMSO_IP_MBD2_rep1.bowtie2.sorted.bam -c HEK293T_DMSO_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir DMSO_IP_MBD2_rep1_TvsC -n DMSO_IP_MBD2_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_DMSO_IP_MBD2_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_DMSO_IP_MBD2_rep2.bowtie2.sorted.bam -c HEK293T_DMSO_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir DMSO_IP_MBD2_rep2_TvsC -n DMSO_IP_MBD2_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_DMSO_IP_MBD2_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_DMSO_IP_MBD2_rep3.bowtie2.sorted.bam -c HEK293T_DMSO_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir DMSO_IP_MBD2_rep3_TvsC -n DMSO_IP_MBD2_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_DMSO_IP_MBD2_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_DMSO_IP_MBD2_rep1_TvsC.runlog
tail -f macs2_callpeak_for_DMSO_IP_MBD2_rep2_TvsC.runlog
tail -f macs2_callpeak_for_DMSO_IP_MBD2_rep3_TvsC.runlog

# DMSO_IP_MBD3
nohup macs2 callpeak -t HEK293T_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam -c HEK293T_DMSO_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir DMSO_IP_MBD3_rep1_TvsC -n DMSO_IP_MBD3_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_DMSO_IP_MBD3_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam -c HEK293T_DMSO_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir DMSO_IP_MBD3_rep2_TvsC -n DMSO_IP_MBD3_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_DMSO_IP_MBD3_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam -c HEK293T_DMSO_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir DMSO_IP_MBD3_rep3_TvsC -n DMSO_IP_MBD3_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_DMSO_IP_MBD3_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_DMSO_IP_MBD3_rep1_TvsC.runlog
tail -f macs2_callpeak_for_DMSO_IP_MBD3_rep2_TvsC.runlog
tail -f macs2_callpeak_for_DMSO_IP_MBD3_rep3_TvsC.runlog

# DMSO_IP_TET1
nohup macs2 callpeak -t HEK293T_DMSO_IP_TET1_rep1.bowtie2.sorted.bam -c HEK293T_DMSO_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir DMSO_IP_TET1_rep1_TvsC -n DMSO_IP_TET1_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_DMSO_IP_TET1_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_DMSO_IP_TET1_rep2.bowtie2.sorted.bam -c HEK293T_DMSO_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir DMSO_IP_TET1_rep2_TvsC -n DMSO_IP_TET1_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_DMSO_IP_TET1_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_DMSO_IP_TET1_rep3.bowtie2.sorted.bam -c HEK293T_DMSO_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir DMSO_IP_TET1_rep3_TvsC -n DMSO_IP_TET1_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_DMSO_IP_TET1_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_DMSO_IP_TET1_rep1_TvsC.runlog
tail -f macs2_callpeak_for_DMSO_IP_TET1_rep2_TvsC.runlog
tail -f macs2_callpeak_for_DMSO_IP_TET1_rep3_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_1_for_DMSO_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./DMSO_IP_CTCF_rep1_TvsC/DMSO_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# 3769
# ./DMSO_IP_CTCF_rep2_TvsC/DMSO_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# 5006
# ./DMSO_IP_CTCF_rep3_TvsC/DMSO_IP_CTCF_rep3_TvsC_peaks.narrowPeak
# 3636
# ./DMSO_IP_MBD2_rep1_TvsC/DMSO_IP_MBD2_rep1_TvsC_peaks.narrowPeak
# 6350
# ./DMSO_IP_MBD2_rep2_TvsC/DMSO_IP_MBD2_rep2_TvsC_peaks.narrowPeak
# 1389
# ./DMSO_IP_MBD2_rep3_TvsC/DMSO_IP_MBD2_rep3_TvsC_peaks.narrowPeak
# 3581
# ./DMSO_IP_MBD3_rep1_TvsC/DMSO_IP_MBD3_rep1_TvsC_peaks.narrowPeak
# 1078
# ./DMSO_IP_MBD3_rep2_TvsC/DMSO_IP_MBD3_rep2_TvsC_peaks.narrowPeak
# 2221
# ./DMSO_IP_MBD3_rep3_TvsC/DMSO_IP_MBD3_rep3_TvsC_peaks.narrowPeak
# 1886
# ./DMSO_IP_TET1_rep1_TvsC/DMSO_IP_TET1_rep1_TvsC_peaks.narrowPeak
# 873
# ./DMSO_IP_TET1_rep2_TvsC/DMSO_IP_TET1_rep2_TvsC_peaks.narrowPeak
# 706
# ./DMSO_IP_TET1_rep3_TvsC/DMSO_IP_TET1_rep3_TvsC_peaks.narrowPeak
# 1321

# 统计 peak summit 数量
for i in `find . -name *summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./DMSO_IP_CTCF_rep1_TvsC/DMSO_IP_CTCF_rep1_TvsC_summits.bed
# 3851
# ./DMSO_IP_CTCF_rep2_TvsC/DMSO_IP_CTCF_rep2_TvsC_summits.bed
# 5156
# ./DMSO_IP_CTCF_rep3_TvsC/DMSO_IP_CTCF_rep3_TvsC_summits.bed
# 3688
# ./DMSO_IP_MBD2_rep1_TvsC/DMSO_IP_MBD2_rep1_TvsC_summits.bed
# 6479
# ./DMSO_IP_MBD2_rep2_TvsC/DMSO_IP_MBD2_rep2_TvsC_summits.bed
# 1404
# ./DMSO_IP_MBD2_rep3_TvsC/DMSO_IP_MBD2_rep3_TvsC_summits.bed
# 3628
# ./DMSO_IP_MBD3_rep1_TvsC/DMSO_IP_MBD3_rep1_TvsC_summits.bed
# 1108
# ./DMSO_IP_MBD3_rep2_TvsC/DMSO_IP_MBD3_rep2_TvsC_summits.bed
# 2321
# ./DMSO_IP_MBD3_rep3_TvsC/DMSO_IP_MBD3_rep3_TvsC_summits.bed
# 1950
# ./DMSO_IP_TET1_rep1_TvsC/DMSO_IP_TET1_rep1_TvsC_summits.bed
# 896
# ./DMSO_IP_TET1_rep2_TvsC/DMSO_IP_TET1_rep2_TvsC_summits.bed
# 746
# ./DMSO_IP_TET1_rep3_TvsC/DMSO_IP_TET1_rep3_TvsC_summits.bed
# 1367

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_DMSO_ChIP_Seq_samples_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_1_for_DMSO_each_repeat
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bigwig|sort`;do echo $i;done

```



### BPK25 ChIP-Seq

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged
mkdir macs2_peak_calling_1_for_BPK25_each_repeat
cd macs2_peak_calling_1_for_BPK25_each_repeat
# 上传 bam 连接构造脚本并运行
bash 5_make_BPK25_all_ChIP_Seq_sorted_bam_file_links.sh
ls -1 *.bam

# BPK25_IP_CTCF
nohup macs2 callpeak -t HEK293T_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam -c HEK293T_BPK25_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir BPK25_IP_CTCF_rep1_TvsC -n BPK25_IP_CTCF_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_BPK25_IP_CTCF_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam -c HEK293T_BPK25_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir BPK25_IP_CTCF_rep2_TvsC -n BPK25_IP_CTCF_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_BPK25_IP_CTCF_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam -c HEK293T_BPK25_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir BPK25_IP_CTCF_rep3_TvsC -n BPK25_IP_CTCF_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_BPK25_IP_CTCF_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_BPK25_IP_CTCF_rep1_TvsC.runlog
tail -f macs2_callpeak_for_BPK25_IP_CTCF_rep2_TvsC.runlog
tail -f macs2_callpeak_for_BPK25_IP_CTCF_rep3_TvsC.runlog

# BPK25_IP_MBD2
nohup macs2 callpeak -t HEK293T_BPK25_IP_MBD2_rep1.bowtie2.sorted.bam -c HEK293T_BPK25_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir BPK25_IP_MBD2_rep1_TvsC -n BPK25_IP_MBD2_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_BPK25_IP_MBD2_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_BPK25_IP_MBD2_rep2.bowtie2.sorted.bam -c HEK293T_BPK25_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir BPK25_IP_MBD2_rep2_TvsC -n BPK25_IP_MBD2_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_BPK25_IP_MBD2_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_BPK25_IP_MBD2_rep3.bowtie2.sorted.bam -c HEK293T_BPK25_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir BPK25_IP_MBD2_rep3_TvsC -n BPK25_IP_MBD2_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_BPK25_IP_MBD2_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_BPK25_IP_MBD2_rep1_TvsC.runlog
tail -f macs2_callpeak_for_BPK25_IP_MBD2_rep2_TvsC.runlog
tail -f macs2_callpeak_for_BPK25_IP_MBD2_rep3_TvsC.runlog

# BPK25_IP_MBD3
nohup macs2 callpeak -t HEK293T_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam -c HEK293T_BPK25_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir BPK25_IP_MBD3_rep1_TvsC -n BPK25_IP_MBD3_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_BPK25_IP_MBD3_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam -c HEK293T_BPK25_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir BPK25_IP_MBD3_rep2_TvsC -n BPK25_IP_MBD3_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_BPK25_IP_MBD3_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam -c HEK293T_BPK25_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir BPK25_IP_MBD3_rep3_TvsC -n BPK25_IP_MBD3_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_BPK25_IP_MBD3_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_BPK25_IP_MBD3_rep1_TvsC.runlog
tail -f macs2_callpeak_for_BPK25_IP_MBD3_rep2_TvsC.runlog
tail -f macs2_callpeak_for_BPK25_IP_MBD3_rep3_TvsC.runlog

# BPK25_IP_TET1
nohup macs2 callpeak -t HEK293T_BPK25_IP_TET1_rep1.bowtie2.sorted.bam -c HEK293T_BPK25_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir BPK25_IP_TET1_rep1_TvsC -n BPK25_IP_TET1_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_BPK25_IP_TET1_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_BPK25_IP_TET1_rep2.bowtie2.sorted.bam -c HEK293T_BPK25_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir BPK25_IP_TET1_rep2_TvsC -n BPK25_IP_TET1_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_BPK25_IP_TET1_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_BPK25_IP_TET1_rep3.bowtie2.sorted.bam -c HEK293T_BPK25_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir BPK25_IP_TET1_rep3_TvsC -n BPK25_IP_TET1_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_BPK25_IP_TET1_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_BPK25_IP_TET1_rep1_TvsC.runlog
tail -f macs2_callpeak_for_BPK25_IP_TET1_rep2_TvsC.runlog
tail -f macs2_callpeak_for_BPK25_IP_TET1_rep3_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_1_for_BPK25_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./BPK25_IP_CTCF_rep1_TvsC/BPK25_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# 144
# ./BPK25_IP_CTCF_rep2_TvsC/BPK25_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# 274
# ./BPK25_IP_CTCF_rep3_TvsC/BPK25_IP_CTCF_rep3_TvsC_peaks.narrowPeak
# 601
# ./BPK25_IP_MBD2_rep1_TvsC/BPK25_IP_MBD2_rep1_TvsC_peaks.narrowPeak
# 121
# ./BPK25_IP_MBD2_rep2_TvsC/BPK25_IP_MBD2_rep2_TvsC_peaks.narrowPeak
# 225
# ./BPK25_IP_MBD2_rep3_TvsC/BPK25_IP_MBD2_rep3_TvsC_peaks.narrowPeak
# 84
# ./BPK25_IP_MBD3_rep1_TvsC/BPK25_IP_MBD3_rep1_TvsC_peaks.narrowPeak
# 31
# ./BPK25_IP_MBD3_rep2_TvsC/BPK25_IP_MBD3_rep2_TvsC_peaks.narrowPeak
# 311
# ./BPK25_IP_MBD3_rep3_TvsC/BPK25_IP_MBD3_rep3_TvsC_peaks.narrowPeak
# 756
# ./BPK25_IP_TET1_rep1_TvsC/BPK25_IP_TET1_rep1_TvsC_peaks.narrowPeak
# 91
# ./BPK25_IP_TET1_rep2_TvsC/BPK25_IP_TET1_rep2_TvsC_peaks.narrowPeak
# 119
# ./BPK25_IP_TET1_rep3_TvsC/BPK25_IP_TET1_rep3_TvsC_peaks.narrowPeak
# 18

# 统计 peak summit 数量
for i in `find . -name *summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./BPK25_IP_CTCF_rep1_TvsC/BPK25_IP_CTCF_rep1_TvsC_summits.bed
# 145
# ./BPK25_IP_CTCF_rep2_TvsC/BPK25_IP_CTCF_rep2_TvsC_summits.bed
# 287
# ./BPK25_IP_CTCF_rep3_TvsC/BPK25_IP_CTCF_rep3_TvsC_summits.bed
# 604
# ./BPK25_IP_MBD2_rep1_TvsC/BPK25_IP_MBD2_rep1_TvsC_summits.bed
# 132
# ./BPK25_IP_MBD2_rep2_TvsC/BPK25_IP_MBD2_rep2_TvsC_summits.bed
# 254
# ./BPK25_IP_MBD2_rep3_TvsC/BPK25_IP_MBD2_rep3_TvsC_summits.bed
# 86
# ./BPK25_IP_MBD3_rep1_TvsC/BPK25_IP_MBD3_rep1_TvsC_summits.bed
# 31
# ./BPK25_IP_MBD3_rep2_TvsC/BPK25_IP_MBD3_rep2_TvsC_summits.bed
# 345
# ./BPK25_IP_MBD3_rep3_TvsC/BPK25_IP_MBD3_rep3_TvsC_summits.bed
# 777
# ./BPK25_IP_TET1_rep1_TvsC/BPK25_IP_TET1_rep1_TvsC_summits.bed
# 95
# ./BPK25_IP_TET1_rep2_TvsC/BPK25_IP_TET1_rep2_TvsC_summits.bed
# 122
# ./BPK25_IP_TET1_rep3_TvsC/BPK25_IP_TET1_rep3_TvsC_summits.bed
# 19

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_BPK25_ChIP_Seq_samples_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_1_for_BPK25_each_repeat
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bigwig|sort`;do echo $i;done

```


### IAA ChIP-Seq

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged
mkdir macs2_peak_calling_1_for_IAA_each_repeat
cd macs2_peak_calling_1_for_IAA_each_repeat
# 上传 bam 连接构造脚本并运行
bash 5_make_IAA_all_ChIP_Seq_sorted_bam_file_links.sh
ls -1 *.bam

# IAA_IP_CTCF
nohup macs2 callpeak -t HEK293T_IAA_IP_CTCF_rep1.bowtie2.sorted.bam -c HEK293T_IAA_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_IP_CTCF_rep1_TvsC -n IAA_IP_CTCF_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_IP_CTCF_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_IAA_IP_CTCF_rep2.bowtie2.sorted.bam -c HEK293T_IAA_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_IP_CTCF_rep2_TvsC -n IAA_IP_CTCF_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_IP_CTCF_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_IAA_IP_CTCF_rep3.bowtie2.sorted.bam -c HEK293T_IAA_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_IP_CTCF_rep3_TvsC -n IAA_IP_CTCF_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_IP_CTCF_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_IAA_IP_CTCF_rep1_TvsC.runlog
tail -f macs2_callpeak_for_IAA_IP_CTCF_rep2_TvsC.runlog
tail -f macs2_callpeak_for_IAA_IP_CTCF_rep3_TvsC.runlog

# IAA_IP_MBD2
nohup macs2 callpeak -t HEK293T_IAA_IP_MBD2_rep1.bowtie2.sorted.bam -c HEK293T_IAA_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_IP_MBD2_rep1_TvsC -n IAA_IP_MBD2_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_IP_MBD2_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_IAA_IP_MBD2_rep2.bowtie2.sorted.bam -c HEK293T_IAA_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_IP_MBD2_rep2_TvsC -n IAA_IP_MBD2_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_IP_MBD2_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_IAA_IP_MBD2_rep3.bowtie2.sorted.bam -c HEK293T_IAA_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_IP_MBD2_rep3_TvsC -n IAA_IP_MBD2_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_IP_MBD2_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_IAA_IP_MBD2_rep1_TvsC.runlog
tail -f macs2_callpeak_for_IAA_IP_MBD2_rep2_TvsC.runlog
tail -f macs2_callpeak_for_IAA_IP_MBD2_rep3_TvsC.runlog

# IAA_IP_MBD3
nohup macs2 callpeak -t HEK293T_IAA_IP_MBD3_rep1.bowtie2.sorted.bam -c HEK293T_IAA_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_IP_MBD3_rep1_TvsC -n IAA_IP_MBD3_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_IP_MBD3_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_IAA_IP_MBD3_rep2.bowtie2.sorted.bam -c HEK293T_IAA_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_IP_MBD3_rep2_TvsC -n IAA_IP_MBD3_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_IP_MBD3_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_IAA_IP_MBD3_rep3.bowtie2.sorted.bam -c HEK293T_IAA_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_IP_MBD3_rep3_TvsC -n IAA_IP_MBD3_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_IP_MBD3_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_IAA_IP_MBD3_rep1_TvsC.runlog
tail -f macs2_callpeak_for_IAA_IP_MBD3_rep2_TvsC.runlog
tail -f macs2_callpeak_for_IAA_IP_MBD3_rep3_TvsC.runlog

# IAA_IP_TET1
nohup macs2 callpeak -t HEK293T_IAA_IP_TET1_rep1.bowtie2.sorted.bam -c HEK293T_IAA_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_IP_TET1_rep1_TvsC -n IAA_IP_TET1_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_IP_TET1_rep1_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_IAA_IP_TET1_rep2.bowtie2.sorted.bam -c HEK293T_IAA_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_IP_TET1_rep2_TvsC -n IAA_IP_TET1_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_IP_TET1_rep2_TvsC.runlog 2>&1 & 
nohup macs2 callpeak -t HEK293T_IAA_IP_TET1_rep3.bowtie2.sorted.bam -c HEK293T_IAA_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_IP_TET1_rep3_TvsC -n IAA_IP_TET1_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_IP_TET1_rep3_TvsC.runlog 2>&1 & 
tail -f macs2_callpeak_for_IAA_IP_TET1_rep1_TvsC.runlog
tail -f macs2_callpeak_for_IAA_IP_TET1_rep2_TvsC.runlog
tail -f macs2_callpeak_for_IAA_IP_TET1_rep3_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_1_for_IAA_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./IAA_IP_CTCF_rep1_TvsC/IAA_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# 296
# ./IAA_IP_CTCF_rep2_TvsC/IAA_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# 169
# ./IAA_IP_CTCF_rep3_TvsC/IAA_IP_CTCF_rep3_TvsC_peaks.narrowPeak
# 291
# ./IAA_IP_MBD2_rep1_TvsC/IAA_IP_MBD2_rep1_TvsC_peaks.narrowPeak
# 31
# ./IAA_IP_MBD2_rep2_TvsC/IAA_IP_MBD2_rep2_TvsC_peaks.narrowPeak
# 13
# ./IAA_IP_MBD2_rep3_TvsC/IAA_IP_MBD2_rep3_TvsC_peaks.narrowPeak
# 83
# ./IAA_IP_MBD3_rep1_TvsC/IAA_IP_MBD3_rep1_TvsC_peaks.narrowPeak
# 20
# ./IAA_IP_MBD3_rep2_TvsC/IAA_IP_MBD3_rep2_TvsC_peaks.narrowPeak
# 22
# ./IAA_IP_MBD3_rep3_TvsC/IAA_IP_MBD3_rep3_TvsC_peaks.narrowPeak
# 8
# ./IAA_IP_TET1_rep1_TvsC/IAA_IP_TET1_rep1_TvsC_peaks.narrowPeak
# 25
# ./IAA_IP_TET1_rep2_TvsC/IAA_IP_TET1_rep2_TvsC_peaks.narrowPeak
# 16
# ./IAA_IP_TET1_rep3_TvsC/IAA_IP_TET1_rep3_TvsC_peaks.narrowPeak
# 18

# 统计 peak summit 数量
for i in `find . -name *summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./IAA_IP_CTCF_rep1_TvsC/IAA_IP_CTCF_rep1_TvsC_summits.bed
# 325
# ./IAA_IP_CTCF_rep2_TvsC/IAA_IP_CTCF_rep2_TvsC_summits.bed
# 178
# ./IAA_IP_CTCF_rep3_TvsC/IAA_IP_CTCF_rep3_TvsC_summits.bed
# 315
# ./IAA_IP_MBD2_rep1_TvsC/IAA_IP_MBD2_rep1_TvsC_summits.bed
# 37
# ./IAA_IP_MBD2_rep2_TvsC/IAA_IP_MBD2_rep2_TvsC_summits.bed
# 15
# ./IAA_IP_MBD2_rep3_TvsC/IAA_IP_MBD2_rep3_TvsC_summits.bed
# 107
# ./IAA_IP_MBD3_rep1_TvsC/IAA_IP_MBD3_rep1_TvsC_summits.bed
# 24
# ./IAA_IP_MBD3_rep2_TvsC/IAA_IP_MBD3_rep2_TvsC_summits.bed
# 28
# ./IAA_IP_MBD3_rep3_TvsC/IAA_IP_MBD3_rep3_TvsC_summits.bed
# 9
# ./IAA_IP_TET1_rep1_TvsC/IAA_IP_TET1_rep1_TvsC_summits.bed
# 30
# ./IAA_IP_TET1_rep2_TvsC/IAA_IP_TET1_rep2_TvsC_summits.bed
# 21
# ./IAA_IP_TET1_rep3_TvsC/IAA_IP_TET1_rep3_TvsC_summits.bed
# 20

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_IAA_ChIP_Seq_samples_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_1_for_IAA_each_repeat
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bigwig|sort`;do echo $i;done

```



##  peak calling for rep123 merged

### DMSO  ChIP-Seq

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged
mkdir macs2_peak_calling_2_for_DMSO_rep123_together
cd macs2_peak_calling_2_for_DMSO_rep123_together
# 上传 bam 连接构造脚本并运行
bash 5_make_DMSO_all_rep123_merged_bam_file_links.sh
ls -1 *.bam

# DMSO_IP_CTCF_rep123_TvsC
nohup macs2 callpeak -t HEK293T_DMSO_IP_CTCF_rep123.merged.sorted.bam -c HEK293T_DMSO_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir DMSO_IP_CTCF_rep123_TvsC -n DMSO_IP_CTCF_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_DMSO_IP_CTCF_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_DMSO_IP_CTCF_rep123_TvsC.runlog

# DMSO_IP_MBD2_rep123_TvsC
nohup macs2 callpeak -t HEK293T_DMSO_IP_MBD2_rep123.merged.sorted.bam -c HEK293T_DMSO_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir DMSO_IP_MBD2_rep123_TvsC -n DMSO_IP_MBD2_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_DMSO_IP_MBD2_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_DMSO_IP_MBD2_rep123_TvsC.runlog

# DMSO_IP_MBD3_rep123_TvsC
nohup macs2 callpeak -t HEK293T_DMSO_IP_MBD3_rep123.merged.sorted.bam -c HEK293T_DMSO_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir DMSO_IP_MBD3_rep123_TvsC -n DMSO_IP_MBD3_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_DMSO_IP_MBD3_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_DMSO_IP_MBD3_rep123_TvsC.runlog

# DMSO_IP_TET1_rep123_TvsC
nohup macs2 callpeak -t HEK293T_DMSO_IP_TET1_rep123.merged.sorted.bam -c HEK293T_DMSO_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir DMSO_IP_TET1_rep123_TvsC -n DMSO_IP_TET1_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_DMSO_IP_TET1_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_DMSO_IP_TET1_rep123_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_2_for_DMSO_rep123_together
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./DMSO_IP_CTCF_rep123_TvsC/DMSO_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# 17440
# ./DMSO_IP_MBD2_rep123_TvsC/DMSO_IP_MBD2_rep123_TvsC_peaks.narrowPeak
# 15991
# ./DMSO_IP_MBD3_rep123_TvsC/DMSO_IP_MBD3_rep123_TvsC_peaks.narrowPeak
# 10137
# ./DMSO_IP_TET1_rep123_TvsC/DMSO_IP_TET1_rep123_TvsC_peaks.narrowPeak
# 7211

# 将 narrowPeak 转为uniq sorted bed3 用于后续交集分析
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_2_for_DMSO_rep123_together
for i in `find . -name *narrowPeak`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
find . -name *narrowPeak.bed
# ./DMSO_IP_CTCF_rep123_TvsC/DMSO_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
# ./DMSO_IP_MBD2_rep123_TvsC/DMSO_IP_MBD2_rep123_TvsC_peaks.narrowPeak.bed
# ./DMSO_IP_MBD3_rep123_TvsC/DMSO_IP_MBD3_rep123_TvsC_peaks.narrowPeak.bed
# ./DMSO_IP_TET1_rep123_TvsC/DMSO_IP_TET1_rep123_TvsC_peaks.narrowPeak.bed
head `find . -name *narrowPeak.bed`

# 统计 peak summit 数量
for i in `find . -name *summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./DMSO_IP_CTCF_rep123_TvsC/DMSO_IP_CTCF_rep123_TvsC_summits.bed
# 17945
# ./DMSO_IP_MBD2_rep123_TvsC/DMSO_IP_MBD2_rep123_TvsC_summits.bed
# 16415
# ./DMSO_IP_MBD3_rep123_TvsC/DMSO_IP_MBD3_rep123_TvsC_summits.bed
# 10498
# ./DMSO_IP_TET1_rep123_TvsC/DMSO_IP_TET1_rep123_TvsC_summits.bed
# 7515

# 将 peak summits 转为 extB50 uniq sorted bed3 用于后续交集分析
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_2_for_DMSO_rep123_together
for i in `find . -name *summits.bed|sort`;do echo $i; bedtools slop -b 50 -i $i -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >${i%.bed}.extB50.sorted.bed ;done
find . -name *summits.extB50.sorted.bed
# ./DMSO_IP_CTCF_rep123_TvsC/DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
# ./DMSO_IP_MBD2_rep123_TvsC/DMSO_IP_MBD2_rep123_TvsC_summits.extB50.sorted.bed
# ./DMSO_IP_MBD3_rep123_TvsC/DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed
# ./DMSO_IP_TET1_rep123_TvsC/DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed
head `find . -name *summits.extB50.sorted.bed`


# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_IAA72h_CTCF_ChIP_Seq_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_2_for_DMSO_rep123_together
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 构造 peak summits bed3 (uniq sorted)
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_2_for_DMSO_rep123_together
for i in `find . -name *summits.bed |sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i%.bed}.bed3_uniq.sorted.bed; done
wc -l `find . -name *sorted.bed |sort`
# 17945 ./DMSO_IP_CTCF_rep123_TvsC/DMSO_IP_CTCF_rep123_TvsC_summits.bed3_uniq.sorted.bed
# 16415 ./DMSO_IP_MBD2_rep123_TvsC/DMSO_IP_MBD2_rep123_TvsC_summits.bed3_uniq.sorted.bed
# 10498 ./DMSO_IP_MBD3_rep123_TvsC/DMSO_IP_MBD3_rep123_TvsC_summits.bed3_uniq.sorted.bed
#  7515 ./DMSO_IP_TET1_rep123_TvsC/DMSO_IP_TET1_rep123_TvsC_summits.bed3_uniq.sorted.bed

# 构造 narrow peak bed3  (uniq sorted)
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_2_for_DMSO_rep123_together
for i in `find . -name *narrowPeak`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
wc -l `find . -name *_rep123_TvsC_peaks.narrowPeak.bed`
# 17440 ./DMSO_IP_CTCF_rep123_TvsC/DMSO_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
# 15991 ./DMSO_IP_MBD2_rep123_TvsC/DMSO_IP_MBD2_rep123_TvsC_peaks.narrowPeak.bed
# 10137 ./DMSO_IP_MBD3_rep123_TvsC/DMSO_IP_MBD3_rep123_TvsC_peaks.narrowPeak.bed
#  7211 ./DMSO_IP_TET1_rep123_TvsC/DMSO_IP_TET1_rep123_TvsC_peaks.narrowPeak.bed

```



### BPK25 ChIP-Seq

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged
mkdir macs2_peak_calling_2_for_BPK25_rep123_together
cd macs2_peak_calling_2_for_BPK25_rep123_together
# 上传 bam 连接构造脚本并运行
bash 5_make_BPK25_all_rep123_merged_bam_file_links.sh
ls -1 *.bam

# BPK25_IP_CTCF_rep123_TvsC
nohup macs2 callpeak -t HEK293T_BPK25_IP_CTCF_rep123.merged.sorted.bam -c HEK293T_BPK25_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir BPK25_IP_CTCF_rep123_TvsC -n BPK25_IP_CTCF_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_BPK25_IP_CTCF_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_BPK25_IP_CTCF_rep123_TvsC.runlog

# BPK25_IP_MBD2_rep123_TvsC
nohup macs2 callpeak -t HEK293T_BPK25_IP_MBD2_rep123.merged.sorted.bam -c HEK293T_BPK25_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir BPK25_IP_MBD2_rep123_TvsC -n BPK25_IP_MBD2_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_BPK25_IP_MBD2_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_BPK25_IP_MBD2_rep123_TvsC.runlog

# BPK25_IP_MBD3_rep123_TvsC
nohup macs2 callpeak -t HEK293T_BPK25_IP_MBD3_rep123.merged.sorted.bam -c HEK293T_BPK25_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir BPK25_IP_MBD3_rep123_TvsC -n BPK25_IP_MBD3_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_BPK25_IP_MBD3_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_BPK25_IP_MBD3_rep123_TvsC.runlog

# BPK25_IP_TET1_rep123_TvsC
nohup macs2 callpeak -t HEK293T_BPK25_IP_TET1_rep123.merged.sorted.bam -c HEK293T_BPK25_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir BPK25_IP_TET1_rep123_TvsC -n BPK25_IP_TET1_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_BPK25_IP_TET1_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_BPK25_IP_TET1_rep123_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_2_for_BPK25_rep123_together
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./BPK25_IP_CTCF_rep123_TvsC/BPK25_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# 375
# ./BPK25_IP_MBD2_rep123_TvsC/BPK25_IP_MBD2_rep123_TvsC_peaks.narrowPeak
# 510
# ./BPK25_IP_MBD3_rep123_TvsC/BPK25_IP_MBD3_rep123_TvsC_peaks.narrowPeak
# 374
# ./BPK25_IP_TET1_rep123_TvsC/BPK25_IP_TET1_rep123_TvsC_peaks.narrowPeak
# 777

# 将 narrowPeak 转为uniq sorted bed3 用于后续交集分析
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_2_for_BPK25_rep123_together
for i in `find . -name *narrowPeak`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
find . -name *narrowPeak.bed
head `find . -name *narrowPeak.bed`
# ./BPK25_IP_CTCF_rep123_TvsC/BPK25_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
# ./BPK25_IP_MBD2_rep123_TvsC/BPK25_IP_MBD2_rep123_TvsC_peaks.narrowPeak.bed
# ./BPK25_IP_MBD3_rep123_TvsC/BPK25_IP_MBD3_rep123_TvsC_peaks.narrowPeak.bed
# ./BPK25_IP_TET1_rep123_TvsC/BPK25_IP_TET1_rep123_TvsC_peaks.narrowPeak.bed


# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./BPK25_IP_CTCF_rep123_TvsC/BPK25_IP_CTCF_rep123_TvsC_summits.bed
# 414
# ./BPK25_IP_MBD2_rep123_TvsC/BPK25_IP_MBD2_rep123_TvsC_summits.bed
# 590
# ./BPK25_IP_MBD3_rep123_TvsC/BPK25_IP_MBD3_rep123_TvsC_summits.bed
# 421
# ./BPK25_IP_TET1_rep123_TvsC/BPK25_IP_TET1_rep123_TvsC_summits.bed
# 803

# 将 peak summits 转为 extB50 uniq sorted bed3 用于后续交集分析
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_2_for_BPK25_rep123_together
for i in `find . -name *summits.bed|sort`;do echo $i; bedtools slop -b 50 -i $i -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >${i%.bed}.extB50.sorted.bed ;done
find . -name *summits.extB50.sorted.bed
# ./BPK25_IP_CTCF_rep123_TvsC/BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
# ./BPK25_IP_MBD2_rep123_TvsC/BPK25_IP_MBD2_rep123_TvsC_summits.extB50.sorted.bed
# ./BPK25_IP_MBD3_rep123_TvsC/BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed
# ./BPK25_IP_TET1_rep123_TvsC/BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed
head `find . -name *summits.extB50.sorted.bed`


# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_IAA72h_CTCF_ChIP_Seq_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_2_for_BPK25_rep123_together
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 构造 narrow peak bed3  (uniq sorted)
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_2_for_BPK25_rep123_together
for i in `find . -name *narrowPeak`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
wc -l `find . -name *_rep123_TvsC_peaks.narrowPeak.bed`
# 375 ./BPK25_IP_CTCF_rep123_TvsC/BPK25_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
# 510 ./BPK25_IP_MBD2_rep123_TvsC/BPK25_IP_MBD2_rep123_TvsC_peaks.narrowPeak.bed
# 374 ./BPK25_IP_MBD3_rep123_TvsC/BPK25_IP_MBD3_rep123_TvsC_peaks.narrowPeak.bed
# 777 ./BPK25_IP_TET1_rep123_TvsC/BPK25_IP_TET1_rep123_TvsC_peaks.narrowPeak.bed

```


### IAA ChIP-Seq

```bash
# 数据准备
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged
mkdir macs2_peak_calling_2_for_IAA_rep123_together
cd macs2_peak_calling_2_for_IAA_rep123_together
# 上传 bam 连接构造脚本并运行
bash 5_make_IAA_all_rep123_merged_bam_file_links.sh
ls -1 *.bam

# IAA_IP_CTCF_rep123_TvsC
nohup macs2 callpeak -t HEK293T_IAA_IP_CTCF_rep123.merged.sorted.bam -c HEK293T_IAA_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_IP_CTCF_rep123_TvsC -n IAA_IP_CTCF_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_IP_CTCF_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_IAA_IP_CTCF_rep123_TvsC.runlog

# IAA_IP_MBD2_rep123_TvsC
nohup macs2 callpeak -t HEK293T_IAA_IP_MBD2_rep123.merged.sorted.bam -c HEK293T_IAA_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_IP_MBD2_rep123_TvsC -n IAA_IP_MBD2_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_IP_MBD2_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_IAA_IP_MBD2_rep123_TvsC.runlog

# IAA_IP_MBD3_rep123_TvsC
nohup macs2 callpeak -t HEK293T_IAA_IP_MBD3_rep123.merged.sorted.bam -c HEK293T_IAA_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_IP_MBD3_rep123_TvsC -n IAA_IP_MBD3_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_IP_MBD3_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_IAA_IP_MBD3_rep123_TvsC.runlog

# IAA_IP_TET1_rep123_TvsC
nohup macs2 callpeak -t HEK293T_IAA_IP_TET1_rep123.merged.sorted.bam -c HEK293T_IAA_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir IAA_IP_TET1_rep123_TvsC -n IAA_IP_TET1_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_IAA_IP_TET1_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_IAA_IP_TET1_rep123_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_2_for_IAA_rep123_together
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./IAA_IP_CTCF_rep123_TvsC/IAA_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# 595
# ./IAA_IP_MBD2_rep123_TvsC/IAA_IP_MBD2_rep123_TvsC_peaks.narrowPeak
# 138
# ./IAA_IP_MBD3_rep123_TvsC/IAA_IP_MBD3_rep123_TvsC_peaks.narrowPeak
# 95
# ./IAA_IP_TET1_rep123_TvsC/IAA_IP_TET1_rep123_TvsC_peaks.narrowPeak
# 100

# 将 narrowPeak 转为uniq sorted bed3 用于后续交集分析
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_2_for_IAA_rep123_together
for i in `find . -name *narrowPeak`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
find . -name *narrowPeak.bed
# ./IAA_IP_CTCF_rep123_TvsC/IAA_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
# ./IAA_IP_MBD2_rep123_TvsC/IAA_IP_MBD2_rep123_TvsC_peaks.narrowPeak.bed
# ./IAA_IP_MBD3_rep123_TvsC/IAA_IP_MBD3_rep123_TvsC_peaks.narrowPeak.bed
# ./IAA_IP_TET1_rep123_TvsC/IAA_IP_TET1_rep123_TvsC_peaks.narrowPeak.bed
head `find . -name *narrowPeak.bed`


# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./IAA_IP_CTCF_rep123_TvsC/IAA_IP_CTCF_rep123_TvsC_summits.bed
# 698
# ./IAA_IP_MBD2_rep123_TvsC/IAA_IP_MBD2_rep123_TvsC_summits.bed
# 203
# ./IAA_IP_MBD3_rep123_TvsC/IAA_IP_MBD3_rep123_TvsC_summits.bed
# 136
# ./IAA_IP_TET1_rep123_TvsC/IAA_IP_TET1_rep123_TvsC_summits.bed
# 135

# 将 peak summits 转为 extB50 uniq sorted bed3 用于后续交集分析
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_2_for_IAA_rep123_together
for i in `find . -name *summits.bed|sort`;do echo $i; bedtools slop -b 50 -i $i -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >${i%.bed}.extB50.sorted.bed ;done
find . -name *summits.extB50.sorted.bed
# ./IAA_IP_CTCF_rep123_TvsC/IAA_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
# ./IAA_IP_MBD2_rep123_TvsC/IAA_IP_MBD2_rep123_TvsC_summits.extB50.sorted.bed
# ./IAA_IP_MBD3_rep123_TvsC/IAA_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed
# ./IAA_IP_TET1_rep123_TvsC/IAA_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed
head `find . -name *summits.extB50.sorted.bed`


# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_IAA72h_CTCF_ChIP_Seq_macs2_TvsC_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_2_for_IAA_rep123_together
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名


# 构造 narrow peak bed3  (uniq sorted)
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_2_for_IAA_rep123_together
for i in `find . -name *narrowPeak`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
wc -l `find . -name *_rep123_TvsC_peaks.narrowPeak.bed`
# 595 ./IAA_IP_CTCF_rep123_TvsC/IAA_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed
# 138 ./IAA_IP_MBD2_rep123_TvsC/IAA_IP_MBD2_rep123_TvsC_peaks.narrowPeak.bed
#  95 ./IAA_IP_MBD3_rep123_TvsC/IAA_IP_MBD3_rep123_TvsC_peaks.narrowPeak.bed
# 100 ./IAA_IP_TET1_rep123_TvsC/IAA_IP_TET1_rep123_TvsC_peaks.narrowPeak.bed

```



# Create log2 transformed bigwig files

## DMSO ChIP-Seq

```bash
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged
mkdir make_bamCompare_bigwig_for_DMSO_all_samples
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/make_bamCompare_bigwig_for_DMSO_all_samples
# 上传前面构造好的创建rep1/2/3单独的 bam 文件连接脚本
bash 5_make_DMSO_all_ChIP_Seq_sorted_bam_file_links.sh
ls -1 *rep[123].*bam
# 补充创建 rep123 merged bam / bai 文件连接
bash 5_make_DMSO_all_rep123_merged_bam_file_links.sh
ls -1 *rep123.*bam

# 构造 6_run_bam_comCompare_for_DMSO_all_samples.sh 内容如下
# # Compare for bam of each repeat
# bamCompare -b1 HEK293T_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam -b2 HEK293T_DMSO_input_rep1.bowtie2.sorted.bam -o HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam -b2 HEK293T_DMSO_input_rep2.bowtie2.sorted.bam -o HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam -b2 HEK293T_DMSO_input_rep3.bowtie2.sorted.bam -o HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_DMSO_IP_MBD2_rep1.bowtie2.sorted.bam -b2 HEK293T_DMSO_input_rep1.bowtie2.sorted.bam -o HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_DMSO_IP_MBD2_rep2.bowtie2.sorted.bam -b2 HEK293T_DMSO_input_rep2.bowtie2.sorted.bam -o HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_DMSO_IP_MBD2_rep3.bowtie2.sorted.bam -b2 HEK293T_DMSO_input_rep3.bowtie2.sorted.bam -o HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam -b2 HEK293T_DMSO_input_rep1.bowtie2.sorted.bam -o HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam -b2 HEK293T_DMSO_input_rep2.bowtie2.sorted.bam -o HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam -b2 HEK293T_DMSO_input_rep3.bowtie2.sorted.bam -o HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_DMSO_IP_TET1_rep1.bowtie2.sorted.bam -b2 HEK293T_DMSO_input_rep1.bowtie2.sorted.bam -o HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_DMSO_IP_TET1_rep2.bowtie2.sorted.bam -b2 HEK293T_DMSO_input_rep2.bowtie2.sorted.bam -o HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_DMSO_IP_TET1_rep3.bowtie2.sorted.bam -b2 HEK293T_DMSO_input_rep3.bowtie2.sorted.bam -o HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# # Compare for rep123 merged bam
# bamCompare -b1 HEK293T_DMSO_IP_CTCF_rep123.merged.sorted.bam -b2 HEK293T_DMSO_input_rep123.merged.sorted.bam -o HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_DMSO_IP_MBD2_rep123.merged.sorted.bam -b2 HEK293T_DMSO_input_rep123.merged.sorted.bam -o HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_DMSO_IP_MBD3_rep123.merged.sorted.bam -b2 HEK293T_DMSO_input_rep123.merged.sorted.bam -o HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_DMSO_IP_TET1_rep123.merged.sorted.bam -b2 HEK293T_DMSO_input_rep123.merged.sorted.bam -o HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads

# 上传并运行
nohup bash 6_run_bam_comCompare_for_DMSO_all_samples.sh >6_run_bam_comCompare_for_DMSO_all_samples.sh.runlog 2>&1 &
tail -f 6_run_bam_comCompare_for_DMSO_all_samples.sh.runlog

```

## BPK25 ChIP-Seq

```bash
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged
mkdir make_bamCompare_bigwig_for_BPK25_all_samples
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/make_bamCompare_bigwig_for_BPK25_all_samples
# 上传前面构造好的创建rep1/2/3单独的 bam 文件连接脚本
bash 5_make_BPK25_all_ChIP_Seq_sorted_bam_file_links.sh
ls -1 *rep[123].*bam
# 补充创建 rep123 merged bam / bai 文件连接
bash 5_make_BPK25_all_rep123_merged_bam_file_links.sh
ls -1 *rep123.*bam

# 构造 6_run_bam_comCompare_for_BPK25_all_samples.sh 内容如下
# # Compare for bam of each repeat
# bamCompare -b1 HEK293T_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam -b2 HEK293T_BPK25_input_rep1.bowtie2.sorted.bam -o HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam -b2 HEK293T_BPK25_input_rep2.bowtie2.sorted.bam -o HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam -b2 HEK293T_BPK25_input_rep3.bowtie2.sorted.bam -o HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_BPK25_IP_MBD2_rep1.bowtie2.sorted.bam -b2 HEK293T_BPK25_input_rep1.bowtie2.sorted.bam -o HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_BPK25_IP_MBD2_rep2.bowtie2.sorted.bam -b2 HEK293T_BPK25_input_rep2.bowtie2.sorted.bam -o HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_BPK25_IP_MBD2_rep3.bowtie2.sorted.bam -b2 HEK293T_BPK25_input_rep3.bowtie2.sorted.bam -o HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam -b2 HEK293T_BPK25_input_rep1.bowtie2.sorted.bam -o HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam -b2 HEK293T_BPK25_input_rep2.bowtie2.sorted.bam -o HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam -b2 HEK293T_BPK25_input_rep3.bowtie2.sorted.bam -o HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_BPK25_IP_TET1_rep1.bowtie2.sorted.bam -b2 HEK293T_BPK25_input_rep1.bowtie2.sorted.bam -o HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_BPK25_IP_TET1_rep2.bowtie2.sorted.bam -b2 HEK293T_BPK25_input_rep2.bowtie2.sorted.bam -o HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_BPK25_IP_TET1_rep3.bowtie2.sorted.bam -b2 HEK293T_BPK25_input_rep3.bowtie2.sorted.bam -o HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# # Compare for rep123 merged bam
# bamCompare -b1 HEK293T_BPK25_IP_CTCF_rep123.merged.sorted.bam -b2 HEK293T_BPK25_input_rep123.merged.sorted.bam -o HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_BPK25_IP_MBD2_rep123.merged.sorted.bam -b2 HEK293T_BPK25_input_rep123.merged.sorted.bam -o HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_BPK25_IP_MBD3_rep123.merged.sorted.bam -b2 HEK293T_BPK25_input_rep123.merged.sorted.bam -o HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_BPK25_IP_TET1_rep123.merged.sorted.bam -b2 HEK293T_BPK25_input_rep123.merged.sorted.bam -o HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads

# 上传并运行
nohup bash 6_run_bam_comCompare_for_BPK25_all_samples.sh >6_run_bam_comCompare_for_BPK25_all_samples.sh.runlog 2>&1 &
tail -f 6_run_bam_comCompare_for_BPK25_all_samples.sh.runlog

```

## IAA ChIP-Seq

```bash
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged
mkdir make_bamCompare_bigwig_for_IAA_all_samples
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/make_bamCompare_bigwig_for_IAA_all_samples
# 上传前面构造好的创建rep1/2/3单独的 bam 文件连接脚本
bash 5_make_IAA_all_ChIP_Seq_sorted_bam_file_links.sh
ls -1 *rep[123].*bam
# 补充创建 rep123 merged bam / bai 文件连接
bash 5_make_IAA_all_rep123_merged_bam_file_links.sh
ls -1 *rep123.*bam

# 构造 6_run_bam_comCompare_for_IAA_all_samples.sh 内容如下
# # Compare for bam of each repeat
# bamCompare -b1 HEK293T_IAA_IP_CTCF_rep1.bowtie2.sorted.bam -b2 HEK293T_IAA_input_rep1.bowtie2.sorted.bam -o HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_IAA_IP_CTCF_rep2.bowtie2.sorted.bam -b2 HEK293T_IAA_input_rep2.bowtie2.sorted.bam -o HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_IAA_IP_CTCF_rep3.bowtie2.sorted.bam -b2 HEK293T_IAA_input_rep3.bowtie2.sorted.bam -o HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_IAA_IP_MBD2_rep1.bowtie2.sorted.bam -b2 HEK293T_IAA_input_rep1.bowtie2.sorted.bam -o HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_IAA_IP_MBD2_rep2.bowtie2.sorted.bam -b2 HEK293T_IAA_input_rep2.bowtie2.sorted.bam -o HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_IAA_IP_MBD2_rep3.bowtie2.sorted.bam -b2 HEK293T_IAA_input_rep3.bowtie2.sorted.bam -o HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_IAA_IP_MBD3_rep1.bowtie2.sorted.bam -b2 HEK293T_IAA_input_rep1.bowtie2.sorted.bam -o HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_IAA_IP_MBD3_rep2.bowtie2.sorted.bam -b2 HEK293T_IAA_input_rep2.bowtie2.sorted.bam -o HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_IAA_IP_MBD3_rep3.bowtie2.sorted.bam -b2 HEK293T_IAA_input_rep3.bowtie2.sorted.bam -o HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_IAA_IP_TET1_rep1.bowtie2.sorted.bam -b2 HEK293T_IAA_input_rep1.bowtie2.sorted.bam -o HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_IAA_IP_TET1_rep2.bowtie2.sorted.bam -b2 HEK293T_IAA_input_rep2.bowtie2.sorted.bam -o HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_IAA_IP_TET1_rep3.bowtie2.sorted.bam -b2 HEK293T_IAA_input_rep3.bowtie2.sorted.bam -o HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# # Compare for rep123 merged bam
# bamCompare -b1 HEK293T_IAA_IP_CTCF_rep123.merged.sorted.bam -b2 HEK293T_IAA_input_rep123.merged.sorted.bam -o HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_IAA_IP_MBD2_rep123.merged.sorted.bam -b2 HEK293T_IAA_input_rep123.merged.sorted.bam -o HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_IAA_IP_MBD3_rep123.merged.sorted.bam -b2 HEK293T_IAA_input_rep123.merged.sorted.bam -o HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# bamCompare -b1 HEK293T_IAA_IP_TET1_rep123.merged.sorted.bam -b2 HEK293T_IAA_input_rep123.merged.sorted.bam -o HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig -of bigwig --skipZeroOverZero --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads

# 上传并运行
nohup bash 6_run_bam_comCompare_for_IAA_all_samples.sh >6_run_bam_comCompare_for_IAA_all_samples.sh.runlog 2>&1 &
tail -f 6_run_bam_comCompare_for_IAA_all_samples.sh.runlog

```

# public CTCF rep12 union summit f1k reads enrichment heatmap

## data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged
mkdir plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit

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
# DMSO rep1/2/3 individual
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig
# DMSO rep123 merged
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig
ls -1 HEK293T_DMSO*.bigwig
ls -1 HEK293T_DMSO*.bigwig |wc -l
# 16

# BPK25 rep1/2/3 individual
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig
# BPK25 rep123 merged
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig
ls -1 HEK293T_BPK25*.bigwig
ls -1 HEK293T_BPK25*.bigwig |wc -l
# 16

# IAA rep1/2/3 individual
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig
# IAA rep123 merged
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig
ls -1 HEK293T_IAA*.bigwig
ls -1 HEK293T_IAA*.bigwig |wc -l
# 16

```

## plot use CTCF rep12 union all summits

```bash
# plot use CTCF rep12 union all summits
## compute and plot use all 4 TFs rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig -R CTCF_rep12_narrowPeak.summit.union.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz
plotHeatmap -m bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz -o bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_rep12_union_all_summits(97424)" --samplesLabel "DMSO_IP_CTCF_rep123" "DMSO_IP_MBD2_rep123" "DMSO_IP_MBD3_rep123" "DMSO_IP_TET1_rep123" "BPK25_IP_CTCF_rep123" "BPK25_IP_MBD2_rep123" "BPK25_IP_MBD3_rep123" "BPK25_IP_TET1_rep123" "IAA_IP_CTCF_rep123" "IAA_IP_MBD2_rep123" "IAA_IP_MBD3_rep123" "IAA_IP_TET1_rep123"

## compute and plot use all 4 TFs rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig -R CTCF_rep12_narrowPeak.summit.union.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_all_samples_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz
plotHeatmap -m bamCompare_matrix_of_all_samples_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz -o bamCompare_matrix_of_all_samples_for_plotHeatMap.public_CTCF_rep12_union_all_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_rep12_union_all_summits(97424)" --samplesLabel "DMSO_IP_CTCF_rep1" "DMSO_IP_CTCF_rep2" "DMSO_IP_CTCF_rep3" "DMSO_IP_MBD2_rep1" "DMSO_IP_MBD2_rep2" "DMSO_IP_MBD2_rep3" "DMSO_IP_MBD3_rep1" "DMSO_IP_MBD3_rep2" "DMSO_IP_MBD3_rep3" "DMSO_IP_TET1_rep1" "DMSO_IP_TET1_rep2" "DMSO_IP_TET1_rep3" "BPK25_IP_CTCF_rep1" "BPK25_IP_CTCF_rep2" "BPK25_IP_CTCF_rep3" "BPK25_IP_MBD2_rep1" "BPK25_IP_MBD2_rep2" "BPK25_IP_MBD2_rep3" "BPK25_IP_MBD3_rep1" "BPK25_IP_MBD3_rep2" "BPK25_IP_MBD3_rep3" "BPK25_IP_TET1_rep1" "BPK25_IP_TET1_rep2" "BPK25_IP_TET1_rep3" "IAA_IP_CTCF_rep1" "IAA_IP_CTCF_rep2" "IAA_IP_CTCF_rep3" "IAA_IP_MBD2_rep1" "IAA_IP_MBD2_rep2" "IAA_IP_MBD2_rep3" "IAA_IP_MBD3_rep1" "IAA_IP_MBD3_rep2" "IAA_IP_MBD3_rep3" "IAA_IP_TET1_rep1" "IAA_IP_TET1_rep2" "IAA_IP_TET1_rep3"


## compute and plot use only CTCF rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig -R CTCF_rep12_narrowPeak.summit.union.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_CTCF_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz
plotHeatmap -m bamCompare_matrix_of_CTCF_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz -o bamCompare_matrix_of_CTCF_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_rep12_union_all_summits(97424)" --samplesLabel "DMSO_IP_CTCF_rep123" "BPK25_IP_CTCF_rep123" "IAA_IP_CTCF_rep123"

## compute and plot use only CTCF rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig -R CTCF_rep12_narrowPeak.summit.union.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_CTCF_samples_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz
plotHeatmap -m bamCompare_matrix_of_CTCF_samples_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz -o bamCompare_matrix_of_CTCF_samples_for_plotHeatMap.public_CTCF_rep12_union_all_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_rep12_union_all_summits(97424)" --samplesLabel "DMSO_IP_CTCF_rep1" "DMSO_IP_CTCF_rep2" "DMSO_IP_CTCF_rep3" "BPK25_IP_CTCF_rep1" "BPK25_IP_CTCF_rep2" "BPK25_IP_CTCF_rep3" "IAA_IP_CTCF_rep1" "IAA_IP_CTCF_rep2" "IAA_IP_CTCF_rep3"


## compute and plot use only MBD2 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig -R CTCF_rep12_narrowPeak.summit.union.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD2_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz
plotHeatmap -m bamCompare_matrix_of_MBD2_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz -o bamCompare_matrix_of_MBD2_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_rep12_union_all_summits(97424)" --samplesLabel "DMSO_IP_MBD2_rep123" "BPK25_IP_MBD2_rep123" "IAA_IP_MBD2_rep123"

## compute and plot use only MBD2 rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig -R CTCF_rep12_narrowPeak.summit.union.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD2_samples_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz
plotHeatmap -m bamCompare_matrix_of_MBD2_samples_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz -o bamCompare_matrix_of_MBD2_samples_for_plotHeatMap.public_CTCF_rep12_union_all_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_rep12_union_all_summits(97424)" --samplesLabel "DMSO_IP_MBD2_rep1" "DMSO_IP_MBD2_rep2" "DMSO_IP_MBD2_rep3" "BPK25_IP_MBD2_rep1" "BPK25_IP_MBD2_rep2" "BPK25_IP_MBD2_rep3" "IAA_IP_MBD2_rep1" "IAA_IP_MBD2_rep2" "IAA_IP_MBD2_rep3"


## compute and plot use only MBD3 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig -R CTCF_rep12_narrowPeak.summit.union.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD3_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz
plotHeatmap -m bamCompare_matrix_of_MBD3_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz -o bamCompare_matrix_of_MBD3_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_rep12_union_all_summits(97424)" --samplesLabel "DMSO_IP_MBD3_rep123" "BPK25_IP_MBD3_rep123" "IAA_IP_MBD3_rep123"

## compute and plot use only MBD3 rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig -R CTCF_rep12_narrowPeak.summit.union.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD3_samples_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz
plotHeatmap -m bamCompare_matrix_of_MBD3_samples_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz -o bamCompare_matrix_of_MBD3_samples_for_plotHeatMap.public_CTCF_rep12_union_all_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_rep12_union_all_summits(97424)" --samplesLabel "DMSO_IP_MBD3_rep1" "DMSO_IP_MBD3_rep2" "DMSO_IP_MBD3_rep3" "BPK25_IP_MBD3_rep1" "BPK25_IP_MBD3_rep2" "BPK25_IP_MBD3_rep3" "IAA_IP_MBD3_rep1" "IAA_IP_MBD3_rep2" "IAA_IP_MBD3_rep3"


## compute and plot use only TET1 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig -R CTCF_rep12_narrowPeak.summit.union.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_TET1_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz
plotHeatmap -m bamCompare_matrix_of_TET1_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz -o bamCompare_matrix_of_TET1_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_all_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_rep12_union_all_summits(97424)" --samplesLabel "DMSO_IP_TET1_rep123" "BPK25_IP_TET1_rep123" "IAA_IP_TET1_rep123"

## compute and plot use only TET1 rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig -R CTCF_rep12_narrowPeak.summit.union.all.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_TET1_samples_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz
plotHeatmap -m bamCompare_matrix_of_TET1_samples_for_plotHeatMap.public_CTCF_rep12_union_all_summits.gz -o bamCompare_matrix_of_TET1_samples_for_plotHeatMap.public_CTCF_rep12_union_all_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "HEK293T_CTCF_rep12_union_all_summits(97424)" --samplesLabel "DMSO_IP_TET1_rep1" "DMSO_IP_TET1_rep2" "DMSO_IP_TET1_rep3" "BPK25_IP_TET1_rep1" "BPK25_IP_TET1_rep2" "BPK25_IP_TET1_rep3" "IAA_IP_TET1_rep1" "IAA_IP_TET1_rep2" "IAA_IP_TET1_rep3"

```

## plot use union  source 3 class summit

```bash
# plot use CTCF rep12 union all summits
## compute and plot use all 4 TFs rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig -R CTCF_rep12_union.rep1_only.summit.sorted.bed CTCF_rep12_union.rep12_common.summit.sorted.bed CTCF_rep12_union.rep2_only.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.gz
plotHeatmap -m bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.gz -o bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "rep1_only(47355)_summit" "rep12_common_summit(31737)" "rep2_only_summit(18332)" --samplesLabel "DMSO_IP_CTCF_rep123" "DMSO_IP_MBD2_rep123" "DMSO_IP_MBD3_rep123" "DMSO_IP_TET1_rep123" "BPK25_IP_CTCF_rep123" "BPK25_IP_MBD2_rep123" "BPK25_IP_MBD3_rep123" "BPK25_IP_TET1_rep123" "IAA_IP_CTCF_rep123" "IAA_IP_MBD2_rep123" "IAA_IP_MBD3_rep123" "IAA_IP_TET1_rep123"

## compute and plot use all 4 TFs rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig -R CTCF_rep12_union.rep1_only.summit.sorted.bed CTCF_rep12_union.rep12_common.summit.sorted.bed CTCF_rep12_union.rep2_only.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_all_samples_for_plotHeatMap.public_CTCF_rep12_union_3source.gz
plotHeatmap -m bamCompare_matrix_of_all_samples_for_plotHeatMap.public_CTCF_rep12_union_3source.gz -o bamCompare_matrix_of_all_samples_for_plotHeatMap.public_CTCF_rep12_union_3source.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "rep1_only(47355)_summit" "rep12_common_summit(31737)" "rep2_only_summit(18332)" --samplesLabel "DMSO_IP_CTCF_rep1" "DMSO_IP_CTCF_rep2" "DMSO_IP_CTCF_rep3" "DMSO_IP_MBD2_rep1" "DMSO_IP_MBD2_rep2" "DMSO_IP_MBD2_rep3" "DMSO_IP_MBD3_rep1" "DMSO_IP_MBD3_rep2" "DMSO_IP_MBD3_rep3" "DMSO_IP_TET1_rep1" "DMSO_IP_TET1_rep2" "DMSO_IP_TET1_rep3" "BPK25_IP_CTCF_rep1" "BPK25_IP_CTCF_rep2" "BPK25_IP_CTCF_rep3" "BPK25_IP_MBD2_rep1" "BPK25_IP_MBD2_rep2" "BPK25_IP_MBD2_rep3" "BPK25_IP_MBD3_rep1" "BPK25_IP_MBD3_rep2" "BPK25_IP_MBD3_rep3" "BPK25_IP_TET1_rep1" "BPK25_IP_TET1_rep2" "BPK25_IP_TET1_rep3" "IAA_IP_CTCF_rep1" "IAA_IP_CTCF_rep2" "IAA_IP_CTCF_rep3" "IAA_IP_MBD2_rep1" "IAA_IP_MBD2_rep2" "IAA_IP_MBD2_rep3" "IAA_IP_MBD3_rep1" "IAA_IP_MBD3_rep2" "IAA_IP_MBD3_rep3" "IAA_IP_TET1_rep1" "IAA_IP_TET1_rep2" "IAA_IP_TET1_rep3"


## compute and plot use only CTCF rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig -R CTCF_rep12_union.rep1_only.summit.sorted.bed CTCF_rep12_union.rep12_common.summit.sorted.bed CTCF_rep12_union.rep2_only.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_CTCF_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.gz
plotHeatmap -m bamCompare_matrix_of_CTCF_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.gz -o bamCompare_matrix_of_CTCF_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "rep1_only(47355)_summit" "rep12_common_summit(31737)" "rep2_only_summit(18332)" --samplesLabel "DMSO_IP_CTCF_rep123" "BPK25_IP_CTCF_rep123" "IAA_IP_CTCF_rep123"

## compute and plot use only CTCF rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig -R CTCF_rep12_union.rep1_only.summit.sorted.bed CTCF_rep12_union.rep12_common.summit.sorted.bed CTCF_rep12_union.rep2_only.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_CTCF_samples_for_plotHeatMap.public_CTCF_rep12_union_3source.gz
plotHeatmap -m bamCompare_matrix_of_CTCF_samples_for_plotHeatMap.public_CTCF_rep12_union_3source.gz -o bamCompare_matrix_of_CTCF_samples_for_plotHeatMap.public_CTCF_rep12_union_3source.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "rep1_only(47355)_summit" "rep12_common_summit(31737)" "rep2_only_summit(18332)" --samplesLabel "DMSO_IP_CTCF_rep1" "DMSO_IP_CTCF_rep2" "DMSO_IP_CTCF_rep3" "BPK25_IP_CTCF_rep1" "BPK25_IP_CTCF_rep2" "BPK25_IP_CTCF_rep3" "IAA_IP_CTCF_rep1" "IAA_IP_CTCF_rep2" "IAA_IP_CTCF_rep3"


## compute and plot use only MBD2 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig -R CTCF_rep12_union.rep1_only.summit.sorted.bed CTCF_rep12_union.rep12_common.summit.sorted.bed CTCF_rep12_union.rep2_only.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD2_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.gz
plotHeatmap -m bamCompare_matrix_of_MBD2_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.gz -o bamCompare_matrix_of_MBD2_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "rep1_only(47355)_summit" "rep12_common_summit(31737)" "rep2_only_summit(18332)" --samplesLabel "DMSO_IP_MBD2_rep123" "BPK25_IP_MBD2_rep123" "IAA_IP_MBD2_rep123"

## compute and plot use only MBD2 rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig -R CTCF_rep12_union.rep1_only.summit.sorted.bed CTCF_rep12_union.rep12_common.summit.sorted.bed CTCF_rep12_union.rep2_only.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD2_samples_for_plotHeatMap.public_CTCF_rep12_union_3source.gz
plotHeatmap -m bamCompare_matrix_of_MBD2_samples_for_plotHeatMap.public_CTCF_rep12_union_3source.gz -o bamCompare_matrix_of_MBD2_samples_for_plotHeatMap.public_CTCF_rep12_union_3source.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "rep1_only(47355)_summit" "rep12_common_summit(31737)" "rep2_only_summit(18332)" --samplesLabel "DMSO_IP_MBD2_rep1" "DMSO_IP_MBD2_rep2" "DMSO_IP_MBD2_rep3" "BPK25_IP_MBD2_rep1" "BPK25_IP_MBD2_rep2" "BPK25_IP_MBD2_rep3" "IAA_IP_MBD2_rep1" "IAA_IP_MBD2_rep2" "IAA_IP_MBD2_rep3"


## compute and plot use only MBD3 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig -R CTCF_rep12_union.rep1_only.summit.sorted.bed CTCF_rep12_union.rep12_common.summit.sorted.bed CTCF_rep12_union.rep2_only.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD3_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.gz
plotHeatmap -m bamCompare_matrix_of_MBD3_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.gz -o bamCompare_matrix_of_MBD3_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "rep1_only(47355)_summit" "rep12_common_summit(31737)" "rep2_only_summit(18332)" --samplesLabel "DMSO_IP_MBD3_rep123" "BPK25_IP_MBD3_rep123" "IAA_IP_MBD3_rep123"

## compute and plot use only MBD3 rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig -R CTCF_rep12_union.rep1_only.summit.sorted.bed CTCF_rep12_union.rep12_common.summit.sorted.bed CTCF_rep12_union.rep2_only.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD3_samples_for_plotHeatMap.public_CTCF_rep12_union_3source.gz
plotHeatmap -m bamCompare_matrix_of_MBD3_samples_for_plotHeatMap.public_CTCF_rep12_union_3source.gz -o bamCompare_matrix_of_MBD3_samples_for_plotHeatMap.public_CTCF_rep12_union_3source.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "rep1_only(47355)_summit" "rep12_common_summit(31737)" "rep2_only_summit(18332)" --samplesLabel "DMSO_IP_MBD3_rep1" "DMSO_IP_MBD3_rep2" "DMSO_IP_MBD3_rep3" "BPK25_IP_MBD3_rep1" "BPK25_IP_MBD3_rep2" "BPK25_IP_MBD3_rep3" "IAA_IP_MBD3_rep1" "IAA_IP_MBD3_rep2" "IAA_IP_MBD3_rep3"


## compute and plot use only TET1 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig -R CTCF_rep12_union.rep1_only.summit.sorted.bed CTCF_rep12_union.rep12_common.summit.sorted.bed CTCF_rep12_union.rep2_only.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_TET1_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.gz
plotHeatmap -m bamCompare_matrix_of_TET1_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.gz -o bamCompare_matrix_of_TET1_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union_3source.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "rep1_only(47355)_summit" "rep12_common_summit(31737)" "rep2_only_summit(18332)" --samplesLabel "DMSO_IP_TET1_rep123" "BPK25_IP_TET1_rep123" "IAA_IP_TET1_rep123"

## compute and plot use only TET1 rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig -R CTCF_rep12_union.rep1_only.summit.sorted.bed CTCF_rep12_union.rep12_common.summit.sorted.bed CTCF_rep12_union.rep2_only.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_TET1_samples_for_plotHeatMap.public_CTCF_rep12_union_3source.gz
plotHeatmap -m bamCompare_matrix_of_TET1_samples_for_plotHeatMap.public_CTCF_rep12_union_3source.gz -o bamCompare_matrix_of_TET1_samples_for_plotHeatMap.public_CTCF_rep12_union_3source.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "rep1_only(47355)_summit" "rep12_common_summit(31737)" "rep2_only_summit(18332)" --samplesLabel "DMSO_IP_TET1_rep1" "DMSO_IP_TET1_rep2" "DMSO_IP_TET1_rep3" "BPK25_IP_TET1_rep1" "BPK25_IP_TET1_rep2" "BPK25_IP_TET1_rep3" "IAA_IP_TET1_rep1" "IAA_IP_TET1_rep2" "IAA_IP_TET1_rep3"

```

## plot use union  4 class summit ( peak height quartile)

```bash
# plot use CTCF rep12 union all summits
## compute and plot use all 4 TFs rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig -R CTCF_rep12_union.phq4.summit.sorted.bed CTCF_rep12_union.phq3.summit.sorted.bed CTCF_rep12_union.phq2.summit.sorted.bed CTCF_rep12_union.phq1.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz
plotHeatmap -m bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz -o bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "CTCF_rep12_union_phq4_summit(24356)" "CTCF_rep12_union_phq3_summit(24356)" "CTCF_rep12_union_phq2_summit(24356)" "CTCF_rep12_union_phq1_summit(24356)" --samplesLabel "DMSO_IP_CTCF_rep123" "DMSO_IP_MBD2_rep123" "DMSO_IP_MBD3_rep123" "DMSO_IP_TET1_rep123" "BPK25_IP_CTCF_rep123" "BPK25_IP_MBD2_rep123" "BPK25_IP_MBD3_rep123" "BPK25_IP_TET1_rep123" "IAA_IP_CTCF_rep123" "IAA_IP_MBD2_rep123" "IAA_IP_MBD3_rep123" "IAA_IP_TET1_rep123"

## compute and plot use all 4 TFs rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig -R CTCF_rep12_union.phq4.summit.sorted.bed CTCF_rep12_union.phq3.summit.sorted.bed CTCF_rep12_union.phq2.summit.sorted.bed CTCF_rep12_union.phq1.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_all_samples_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz
plotHeatmap -m bamCompare_matrix_of_all_samples_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz -o bamCompare_matrix_of_all_samples_for_plotHeatMap.public_CTCF_rep12_union.phq1234.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "CTCF_rep12_union_phq4_summit(24356)" "CTCF_rep12_union_phq3_summit(24356)" "CTCF_rep12_union_phq2_summit(24356)" "CTCF_rep12_union_phq1_summit(24356)" --samplesLabel "DMSO_IP_CTCF_rep1" "DMSO_IP_CTCF_rep2" "DMSO_IP_CTCF_rep3" "DMSO_IP_MBD2_rep1" "DMSO_IP_MBD2_rep2" "DMSO_IP_MBD2_rep3" "DMSO_IP_MBD3_rep1" "DMSO_IP_MBD3_rep2" "DMSO_IP_MBD3_rep3" "DMSO_IP_TET1_rep1" "DMSO_IP_TET1_rep2" "DMSO_IP_TET1_rep3" "BPK25_IP_CTCF_rep1" "BPK25_IP_CTCF_rep2" "BPK25_IP_CTCF_rep3" "BPK25_IP_MBD2_rep1" "BPK25_IP_MBD2_rep2" "BPK25_IP_MBD2_rep3" "BPK25_IP_MBD3_rep1" "BPK25_IP_MBD3_rep2" "BPK25_IP_MBD3_rep3" "BPK25_IP_TET1_rep1" "BPK25_IP_TET1_rep2" "BPK25_IP_TET1_rep3" "IAA_IP_CTCF_rep1" "IAA_IP_CTCF_rep2" "IAA_IP_CTCF_rep3" "IAA_IP_MBD2_rep1" "IAA_IP_MBD2_rep2" "IAA_IP_MBD2_rep3" "IAA_IP_MBD3_rep1" "IAA_IP_MBD3_rep2" "IAA_IP_MBD3_rep3" "IAA_IP_TET1_rep1" "IAA_IP_TET1_rep2" "IAA_IP_TET1_rep3"

## compute and plot use only CTCF rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig -R CTCF_rep12_union.phq4.summit.sorted.bed CTCF_rep12_union.phq3.summit.sorted.bed CTCF_rep12_union.phq2.summit.sorted.bed CTCF_rep12_union.phq1.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_CTCF_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz
plotHeatmap -m bamCompare_matrix_of_CTCF_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz -o bamCompare_matrix_of_CTCF_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "CTCF_rep12_union_phq4_summit(24356)" "CTCF_rep12_union_phq3_summit(24356)" "CTCF_rep12_union_phq2_summit(24356)" "CTCF_rep12_union_phq1_summit(24356)" --samplesLabel "DMSO_IP_CTCF_rep123" "BPK25_IP_CTCF_rep123" "IAA_IP_CTCF_rep123"

## compute and plot use only CTCF rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig -R CTCF_rep12_union.phq4.summit.sorted.bed CTCF_rep12_union.phq3.summit.sorted.bed CTCF_rep12_union.phq2.summit.sorted.bed CTCF_rep12_union.phq1.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_CTCF_samples_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz
plotHeatmap -m bamCompare_matrix_of_CTCF_samples_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz -o bamCompare_matrix_of_CTCF_samples_for_plotHeatMap.public_CTCF_rep12_union.phq1234.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "CTCF_rep12_union_phq4_summit(24356)" "CTCF_rep12_union_phq3_summit(24356)" "CTCF_rep12_union_phq2_summit(24356)" "CTCF_rep12_union_phq1_summit(24356)" --samplesLabel "DMSO_IP_CTCF_rep1" "DMSO_IP_CTCF_rep2" "DMSO_IP_CTCF_rep3" "BPK25_IP_CTCF_rep1" "BPK25_IP_CTCF_rep2" "BPK25_IP_CTCF_rep3" "IAA_IP_CTCF_rep1" "IAA_IP_CTCF_rep2" "IAA_IP_CTCF_rep3"


## compute and plot use only MBD2 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig -R CTCF_rep12_union.phq4.summit.sorted.bed CTCF_rep12_union.phq3.summit.sorted.bed CTCF_rep12_union.phq2.summit.sorted.bed CTCF_rep12_union.phq1.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD2_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz
plotHeatmap -m bamCompare_matrix_of_MBD2_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz -o bamCompare_matrix_of_MBD2_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "CTCF_rep12_union_phq4_summit(24356)" "CTCF_rep12_union_phq3_summit(24356)" "CTCF_rep12_union_phq2_summit(24356)" "CTCF_rep12_union_phq1_summit(24356)" --samplesLabel "DMSO_IP_MBD2_rep123" "BPK25_IP_MBD2_rep123" "IAA_IP_MBD2_rep123"

## compute and plot use only MBD2 rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig -R CTCF_rep12_union.phq4.summit.sorted.bed CTCF_rep12_union.phq3.summit.sorted.bed CTCF_rep12_union.phq2.summit.sorted.bed CTCF_rep12_union.phq1.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD2_samples_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz
plotHeatmap -m bamCompare_matrix_of_MBD2_samples_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz -o bamCompare_matrix_of_MBD2_samples_for_plotHeatMap.public_CTCF_rep12_union.phq1234.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "CTCF_rep12_union_phq4_summit(24356)" "CTCF_rep12_union_phq3_summit(24356)" "CTCF_rep12_union_phq2_summit(24356)" "CTCF_rep12_union_phq1_summit(24356)" --samplesLabel "DMSO_IP_MBD2_rep1" "DMSO_IP_MBD2_rep2" "DMSO_IP_MBD2_rep3" "BPK25_IP_MBD2_rep1" "BPK25_IP_MBD2_rep2" "BPK25_IP_MBD2_rep3" "IAA_IP_MBD2_rep1" "IAA_IP_MBD2_rep2" "IAA_IP_MBD2_rep3"


## compute and plot use only MBD3 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig -R CTCF_rep12_union.phq4.summit.sorted.bed CTCF_rep12_union.phq3.summit.sorted.bed CTCF_rep12_union.phq2.summit.sorted.bed CTCF_rep12_union.phq1.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD3_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz
plotHeatmap -m bamCompare_matrix_of_MBD3_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz -o bamCompare_matrix_of_MBD3_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "CTCF_rep12_union_phq4_summit(24356)" "CTCF_rep12_union_phq3_summit(24356)" "CTCF_rep12_union_phq2_summit(24356)" "CTCF_rep12_union_phq1_summit(24356)" --samplesLabel "DMSO_IP_MBD3_rep123" "BPK25_IP_MBD3_rep123" "IAA_IP_MBD3_rep123"

## compute and plot use only MBD3 rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig -R CTCF_rep12_union.phq4.summit.sorted.bed CTCF_rep12_union.phq3.summit.sorted.bed CTCF_rep12_union.phq2.summit.sorted.bed CTCF_rep12_union.phq1.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD3_samples_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz
plotHeatmap -m bamCompare_matrix_of_MBD3_samples_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz -o bamCompare_matrix_of_MBD3_samples_for_plotHeatMap.public_CTCF_rep12_union.phq1234.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "CTCF_rep12_union_phq4_summit(24356)" "CTCF_rep12_union_phq3_summit(24356)" "CTCF_rep12_union_phq2_summit(24356)" "CTCF_rep12_union_phq1_summit(24356)" --samplesLabel "DMSO_IP_MBD3_rep1" "DMSO_IP_MBD3_rep2" "DMSO_IP_MBD3_rep3" "BPK25_IP_MBD3_rep1" "BPK25_IP_MBD3_rep2" "BPK25_IP_MBD3_rep3" "IAA_IP_MBD3_rep1" "IAA_IP_MBD3_rep2" "IAA_IP_MBD3_rep3"


## compute and plot use only TET1 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig -R CTCF_rep12_union.phq4.summit.sorted.bed CTCF_rep12_union.phq3.summit.sorted.bed CTCF_rep12_union.phq2.summit.sorted.bed CTCF_rep12_union.phq1.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_TET1_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz
plotHeatmap -m bamCompare_matrix_of_TET1_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz -o bamCompare_matrix_of_TET1_rep123_merged_for_plotHeatMap.public_CTCF_rep12_union.phq1234.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "CTCF_rep12_union_phq4_summit(24356)" "CTCF_rep12_union_phq3_summit(24356)" "CTCF_rep12_union_phq2_summit(24356)" "CTCF_rep12_union_phq1_summit(24356)" --samplesLabel "DMSO_IP_TET1_rep123" "BPK25_IP_TET1_rep123" "IAA_IP_TET1_rep123"

## compute and plot use only TET1 rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_public_CTCF_rep12_union_summit
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig -R CTCF_rep12_union.phq4.summit.sorted.bed CTCF_rep12_union.phq3.summit.sorted.bed CTCF_rep12_union.phq2.summit.sorted.bed CTCF_rep12_union.phq1.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_TET1_samples_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz
plotHeatmap -m bamCompare_matrix_of_TET1_samples_for_plotHeatMap.public_CTCF_rep12_union.phq1234.gz -o bamCompare_matrix_of_TET1_samples_for_plotHeatMap.public_CTCF_rep12_union.phq1234.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to CTCF summits" --regionsLabel "CTCF_rep12_union_phq4_summit(24356)" "CTCF_rep12_union_phq3_summit(24356)" "CTCF_rep12_union_phq2_summit(24356)" "CTCF_rep12_union_phq1_summit(24356)" --samplesLabel "DMSO_IP_TET1_rep1" "DMSO_IP_TET1_rep2" "DMSO_IP_TET1_rep3" "BPK25_IP_TET1_rep1" "BPK25_IP_TET1_rep2" "BPK25_IP_TET1_rep3" "IAA_IP_TET1_rep1" "IAA_IP_TET1_rep2" "IAA_IP_TET1_rep3"

```

# DMSO CTCF summits f1k reads reads enrichment heatmap

## data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged
mkdir plotHeatMap_of_all_use_DMSO_CTCF_summits
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_CTCF_summits
# 所需的 summit bed 文件在如下分析文件夹中已经构造过了
# /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_2_for_DMSO_rep123_together
# 为所需的bed文件创建软连接创建 
ln -s ../macs2_peak_calling_2_for_DMSO_rep123_together/DMSO_IP_CTCF_rep123_TvsC/DMSO_IP_CTCF_rep123_TvsC_summits.bed3_uniq.sorted.bed DMSO_IP_CTCF_rep123_TvsC_summits.bed3_uniq.sorted.bed
wc -l *.sorted.bed
# 17945 DMSO_IP_CTCF_rep123_TvsC_summits.bed3_uniq.sorted.bed

# 创建全部所需的 bigwig 文件的软连接
# DMSO rep1/2/3 individual
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig
# DMSO rep123 merged
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig
ls -1 HEK293T_DMSO*.bigwig
ls -1 HEK293T_DMSO*.bigwig |wc -l
# 16

# BPK25 rep1/2/3 individual
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig
# BPK25 rep123 merged
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig
ls -1 HEK293T_BPK25*.bigwig
ls -1 HEK293T_BPK25*.bigwig |wc -l
# 16

# IAA rep1/2/3 individual
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig
# IAA rep123 merged
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig
ls -1 HEK293T_IAA*.bigwig
ls -1 HEK293T_IAA*.bigwig |wc -l
# 16

```

## computeMatrix and plotHeatmap use DMSO CTCF summits

```bash
# plot use DMSO_CTCF_summits
## compute and plot use all 4 TFs rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_CTCF_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig -R DMSO_IP_CTCF_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.DMSO_CTCF_summits.gz
plotHeatmap -m bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.DMSO_CTCF_summits.gz -o bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.DMSO_CTCF_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_CTCF_summits" --regionsLabel "DMSO_CTCF_summits(17945)" --samplesLabel "DMSO_IP_CTCF_rep123" "DMSO_IP_MBD2_rep123" "DMSO_IP_MBD3_rep123" "DMSO_IP_TET1_rep123" "BPK25_IP_CTCF_rep123" "BPK25_IP_MBD2_rep123" "BPK25_IP_MBD3_rep123" "BPK25_IP_TET1_rep123" "IAA_IP_CTCF_rep123" "IAA_IP_MBD2_rep123" "IAA_IP_MBD3_rep123" "IAA_IP_TET1_rep123"

## compute and plot use all 4 TFs rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_CTCF_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig -R DMSO_IP_CTCF_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_all_samples_for_plotHeatMap.DMSO_CTCF_summits.gz
plotHeatmap -m bamCompare_matrix_of_all_samples_for_plotHeatMap.DMSO_CTCF_summits.gz -o bamCompare_matrix_of_all_samples_for_plotHeatMap.DMSO_CTCF_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_CTCF_summits" --regionsLabel "DMSO_CTCF_summits(17945)" --samplesLabel "DMSO_IP_CTCF_rep1" "DMSO_IP_CTCF_rep2" "DMSO_IP_CTCF_rep3" "DMSO_IP_MBD2_rep1" "DMSO_IP_MBD2_rep2" "DMSO_IP_MBD2_rep3" "DMSO_IP_MBD3_rep1" "DMSO_IP_MBD3_rep2" "DMSO_IP_MBD3_rep3" "DMSO_IP_TET1_rep1" "DMSO_IP_TET1_rep2" "DMSO_IP_TET1_rep3" "BPK25_IP_CTCF_rep1" "BPK25_IP_CTCF_rep2" "BPK25_IP_CTCF_rep3" "BPK25_IP_MBD2_rep1" "BPK25_IP_MBD2_rep2" "BPK25_IP_MBD2_rep3" "BPK25_IP_MBD3_rep1" "BPK25_IP_MBD3_rep2" "BPK25_IP_MBD3_rep3" "BPK25_IP_TET1_rep1" "BPK25_IP_TET1_rep2" "BPK25_IP_TET1_rep3" "IAA_IP_CTCF_rep1" "IAA_IP_CTCF_rep2" "IAA_IP_CTCF_rep3" "IAA_IP_MBD2_rep1" "IAA_IP_MBD2_rep2" "IAA_IP_MBD2_rep3" "IAA_IP_MBD3_rep1" "IAA_IP_MBD3_rep2" "IAA_IP_MBD3_rep3" "IAA_IP_TET1_rep1" "IAA_IP_TET1_rep2" "IAA_IP_TET1_rep3"


## compute and plot use only CTCF rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_CTCF_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig -R DMSO_IP_CTCF_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_CTCF_rep123_merged_for_plotHeatMap.DMSO_CTCF_summits.gz
plotHeatmap -m bamCompare_matrix_of_CTCF_rep123_merged_for_plotHeatMap.DMSO_CTCF_summits.gz -o bamCompare_matrix_of_CTCF_rep123_merged_for_plotHeatMap.DMSO_CTCF_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_CTCF_summits" --regionsLabel "DMSO_CTCF_summits(17945)" --samplesLabel "DMSO_IP_CTCF_rep123" "BPK25_IP_CTCF_rep123" "IAA_IP_CTCF_rep123"

## compute and plot use only CTCF rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_CTCF_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig -R DMSO_IP_CTCF_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_CTCF_samples_for_plotHeatMap.DMSO_CTCF_summits.gz
plotHeatmap -m bamCompare_matrix_of_CTCF_samples_for_plotHeatMap.DMSO_CTCF_summits.gz -o bamCompare_matrix_of_CTCF_samples_for_plotHeatMap.DMSO_CTCF_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_CTCF_summits" --regionsLabel "DMSO_CTCF_summits(17945)" --samplesLabel "DMSO_IP_CTCF_rep1" "DMSO_IP_CTCF_rep2" "DMSO_IP_CTCF_rep3" "BPK25_IP_CTCF_rep1" "BPK25_IP_CTCF_rep2" "BPK25_IP_CTCF_rep3" "IAA_IP_CTCF_rep1" "IAA_IP_CTCF_rep2" "IAA_IP_CTCF_rep3"


## compute and plot use only MBD2 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_CTCF_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig -R DMSO_IP_CTCF_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD2_rep123_merged_for_plotHeatMap.DMSO_CTCF_summits.gz
plotHeatmap -m bamCompare_matrix_of_MBD2_rep123_merged_for_plotHeatMap.DMSO_CTCF_summits.gz -o bamCompare_matrix_of_MBD2_rep123_merged_for_plotHeatMap.DMSO_CTCF_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_CTCF_summits" --regionsLabel "DMSO_CTCF_summits(17945)" --samplesLabel "DMSO_IP_MBD2_rep123" "BPK25_IP_MBD2_rep123" "IAA_IP_MBD2_rep123"

## compute and plot use only MBD2 rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_CTCF_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig -R DMSO_IP_CTCF_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD2_samples_for_plotHeatMap.DMSO_CTCF_summits.gz
plotHeatmap -m bamCompare_matrix_of_MBD2_samples_for_plotHeatMap.DMSO_CTCF_summits.gz -o bamCompare_matrix_of_MBD2_samples_for_plotHeatMap.DMSO_CTCF_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_CTCF_summits" --regionsLabel "DMSO_CTCF_summits(17945)" --samplesLabel "DMSO_IP_MBD2_rep1" "DMSO_IP_MBD2_rep2" "DMSO_IP_MBD2_rep3" "BPK25_IP_MBD2_rep1" "BPK25_IP_MBD2_rep2" "BPK25_IP_MBD2_rep3" "IAA_IP_MBD2_rep1" "IAA_IP_MBD2_rep2" "IAA_IP_MBD2_rep3"


## compute and plot use only MBD3 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_CTCF_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig -R DMSO_IP_CTCF_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD3_rep123_merged_for_plotHeatMap.DMSO_CTCF_summits.gz
plotHeatmap -m bamCompare_matrix_of_MBD3_rep123_merged_for_plotHeatMap.DMSO_CTCF_summits.gz -o bamCompare_matrix_of_MBD3_rep123_merged_for_plotHeatMap.DMSO_CTCF_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_CTCF_summits" --regionsLabel "DMSO_CTCF_summits(17945)" --samplesLabel "DMSO_IP_MBD3_rep123" "BPK25_IP_MBD3_rep123" "IAA_IP_MBD3_rep123"

## compute and plot use only MBD3 rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_CTCF_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig -R DMSO_IP_CTCF_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD3_samples_for_plotHeatMap.DMSO_CTCF_summits.gz
plotHeatmap -m bamCompare_matrix_of_MBD3_samples_for_plotHeatMap.DMSO_CTCF_summits.gz -o bamCompare_matrix_of_MBD3_samples_for_plotHeatMap.DMSO_CTCF_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_CTCF_summits" --regionsLabel "DMSO_CTCF_summits(17945)" --samplesLabel "DMSO_IP_MBD3_rep1" "DMSO_IP_MBD3_rep2" "DMSO_IP_MBD3_rep3" "BPK25_IP_MBD3_rep1" "BPK25_IP_MBD3_rep2" "BPK25_IP_MBD3_rep3" "IAA_IP_MBD3_rep1" "IAA_IP_MBD3_rep2" "IAA_IP_MBD3_rep3"


## compute and plot use only TET1 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_CTCF_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig -R DMSO_IP_CTCF_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_TET1_rep123_merged_for_plotHeatMap.DMSO_CTCF_summits.gz
plotHeatmap -m bamCompare_matrix_of_TET1_rep123_merged_for_plotHeatMap.DMSO_CTCF_summits.gz -o bamCompare_matrix_of_TET1_rep123_merged_for_plotHeatMap.DMSO_CTCF_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_CTCF_summits" --regionsLabel "DMSO_CTCF_summits(17945)" --samplesLabel "DMSO_IP_TET1_rep123" "BPK25_IP_TET1_rep123" "IAA_IP_TET1_rep123"

## compute and plot use only TET1 rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_CTCF_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig -R DMSO_IP_CTCF_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_TET1_samples_for_plotHeatMap.DMSO_CTCF_summits.gz
plotHeatmap -m bamCompare_matrix_of_TET1_samples_for_plotHeatMap.DMSO_CTCF_summits.gz -o bamCompare_matrix_of_TET1_samples_for_plotHeatMap.DMSO_CTCF_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_CTCF_summits" --regionsLabel "DMSO_CTCF_summits(17945)" --samplesLabel "DMSO_IP_TET1_rep1" "DMSO_IP_TET1_rep2" "DMSO_IP_TET1_rep3" "BPK25_IP_TET1_rep1" "BPK25_IP_TET1_rep2" "BPK25_IP_TET1_rep3" "IAA_IP_TET1_rep1" "IAA_IP_TET1_rep2" "IAA_IP_TET1_rep3"

```



# DMSO MBD2 summits f1k reads enrichment heatmap

## data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged
mkdir plotHeatMap_of_all_use_DMSO_MBD2_summits
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_MBD2_summits
# 所需的 summit bed 文件在如下分析文件夹中已经构造过了
# /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_2_for_DMSO_rep123_together
# 为所需的bed文件创建软连接创建 
ln -s ../macs2_peak_calling_2_for_DMSO_rep123_together/DMSO_IP_MBD2_rep123_TvsC/DMSO_IP_MBD2_rep123_TvsC_summits.bed3_uniq.sorted.bed DMSO_IP_MBD2_rep123_TvsC_summits.bed3_uniq.sorted.bed
wc -l *.sorted.bed
# 16415 DMSO_IP_MBD2_rep123_TvsC_summits.bed3_uniq.sorted.bed

# 创建全部所需的 bigwig 文件的软连接
# DMSO rep1/2/3 individual
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig
# DMSO rep123 merged
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig
ls -1 HEK293T_DMSO*.bigwig
ls -1 HEK293T_DMSO*.bigwig |wc -l
# 16

# BPK25 rep1/2/3 individual
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig
# BPK25 rep123 merged
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig
ls -1 HEK293T_BPK25*.bigwig
ls -1 HEK293T_BPK25*.bigwig |wc -l
# 16

# IAA rep1/2/3 individual
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig
# IAA rep123 merged
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig
ls -1 HEK293T_IAA*.bigwig
ls -1 HEK293T_IAA*.bigwig |wc -l
# 16

```

## computeMatrix and plotHeatmap use DMSO MBD2 summits

```bash
# plot use DMSO_MBD2_summits
## compute and plot use all 4 TFs rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_MBD2_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig -R DMSO_IP_MBD2_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.DMSO_MBD2_summits.gz
plotHeatmap -m bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.DMSO_MBD2_summits.gz -o bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.DMSO_MBD2_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_MBD2_summits" --regionsLabel "DMSO_MBD2_summits(16415)" --samplesLabel "DMSO_IP_CTCF_rep123" "DMSO_IP_MBD2_rep123" "DMSO_IP_MBD3_rep123" "DMSO_IP_TET1_rep123" "BPK25_IP_CTCF_rep123" "BPK25_IP_MBD2_rep123" "BPK25_IP_MBD3_rep123" "BPK25_IP_TET1_rep123" "IAA_IP_CTCF_rep123" "IAA_IP_MBD2_rep123" "IAA_IP_MBD3_rep123" "IAA_IP_TET1_rep123"

## compute and plot use all 4 TFs rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_MBD2_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig -R DMSO_IP_MBD2_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_all_samples_for_plotHeatMap.DMSO_MBD2_summits.gz
plotHeatmap -m bamCompare_matrix_of_all_samples_for_plotHeatMap.DMSO_MBD2_summits.gz -o bamCompare_matrix_of_all_samples_for_plotHeatMap.DMSO_MBD2_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_MBD2_summits" --regionsLabel "DMSO_MBD2_summits(16415)" --samplesLabel "DMSO_IP_CTCF_rep1" "DMSO_IP_CTCF_rep2" "DMSO_IP_CTCF_rep3" "DMSO_IP_MBD2_rep1" "DMSO_IP_MBD2_rep2" "DMSO_IP_MBD2_rep3" "DMSO_IP_MBD3_rep1" "DMSO_IP_MBD3_rep2" "DMSO_IP_MBD3_rep3" "DMSO_IP_TET1_rep1" "DMSO_IP_TET1_rep2" "DMSO_IP_TET1_rep3" "BPK25_IP_CTCF_rep1" "BPK25_IP_CTCF_rep2" "BPK25_IP_CTCF_rep3" "BPK25_IP_MBD2_rep1" "BPK25_IP_MBD2_rep2" "BPK25_IP_MBD2_rep3" "BPK25_IP_MBD3_rep1" "BPK25_IP_MBD3_rep2" "BPK25_IP_MBD3_rep3" "BPK25_IP_TET1_rep1" "BPK25_IP_TET1_rep2" "BPK25_IP_TET1_rep3" "IAA_IP_CTCF_rep1" "IAA_IP_CTCF_rep2" "IAA_IP_CTCF_rep3" "IAA_IP_MBD2_rep1" "IAA_IP_MBD2_rep2" "IAA_IP_MBD2_rep3" "IAA_IP_MBD3_rep1" "IAA_IP_MBD3_rep2" "IAA_IP_MBD3_rep3" "IAA_IP_TET1_rep1" "IAA_IP_TET1_rep2" "IAA_IP_TET1_rep3"


## compute and plot use only CTCF rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_MBD2_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig -R DMSO_IP_MBD2_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_CTCF_rep123_merged_for_plotHeatMap.DMSO_MBD2_summits.gz
plotHeatmap -m bamCompare_matrix_of_CTCF_rep123_merged_for_plotHeatMap.DMSO_MBD2_summits.gz -o bamCompare_matrix_of_CTCF_rep123_merged_for_plotHeatMap.DMSO_MBD2_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_MBD2_summits" --regionsLabel "DMSO_MBD2_summits(16415)" --samplesLabel "DMSO_IP_CTCF_rep123" "BPK25_IP_CTCF_rep123" "IAA_IP_CTCF_rep123"

## compute and plot use only CTCF rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_MBD2_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig -R DMSO_IP_MBD2_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_CTCF_samples_for_plotHeatMap.DMSO_MBD2_summits.gz
plotHeatmap -m bamCompare_matrix_of_CTCF_samples_for_plotHeatMap.DMSO_MBD2_summits.gz -o bamCompare_matrix_of_CTCF_samples_for_plotHeatMap.DMSO_MBD2_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_MBD2_summits" --regionsLabel "DMSO_MBD2_summits(16415)" --samplesLabel "DMSO_IP_CTCF_rep1" "DMSO_IP_CTCF_rep2" "DMSO_IP_CTCF_rep3" "BPK25_IP_CTCF_rep1" "BPK25_IP_CTCF_rep2" "BPK25_IP_CTCF_rep3" "IAA_IP_CTCF_rep1" "IAA_IP_CTCF_rep2" "IAA_IP_CTCF_rep3"


## compute and plot use only MBD2 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_MBD2_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig -R DMSO_IP_MBD2_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD2_rep123_merged_for_plotHeatMap.DMSO_MBD2_summits.gz
plotHeatmap -m bamCompare_matrix_of_MBD2_rep123_merged_for_plotHeatMap.DMSO_MBD2_summits.gz -o bamCompare_matrix_of_MBD2_rep123_merged_for_plotHeatMap.DMSO_MBD2_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_MBD2_summits" --regionsLabel "DMSO_MBD2_summits(16415)" --samplesLabel "DMSO_IP_MBD2_rep123" "BPK25_IP_MBD2_rep123" "IAA_IP_MBD2_rep123"

## compute and plot use only MBD2 rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_MBD2_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig -R DMSO_IP_MBD2_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD2_samples_for_plotHeatMap.DMSO_MBD2_summits.gz
plotHeatmap -m bamCompare_matrix_of_MBD2_samples_for_plotHeatMap.DMSO_MBD2_summits.gz -o bamCompare_matrix_of_MBD2_samples_for_plotHeatMap.DMSO_MBD2_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_MBD2_summits" --regionsLabel "DMSO_MBD2_summits(16415)" --samplesLabel "DMSO_IP_MBD2_rep1" "DMSO_IP_MBD2_rep2" "DMSO_IP_MBD2_rep3" "BPK25_IP_MBD2_rep1" "BPK25_IP_MBD2_rep2" "BPK25_IP_MBD2_rep3" "IAA_IP_MBD2_rep1" "IAA_IP_MBD2_rep2" "IAA_IP_MBD2_rep3"


## compute and plot use only MBD3 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_MBD2_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig -R DMSO_IP_MBD2_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD3_rep123_merged_for_plotHeatMap.DMSO_MBD2_summits.gz
plotHeatmap -m bamCompare_matrix_of_MBD3_rep123_merged_for_plotHeatMap.DMSO_MBD2_summits.gz -o bamCompare_matrix_of_MBD3_rep123_merged_for_plotHeatMap.DMSO_MBD2_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_MBD2_summits" --regionsLabel "DMSO_MBD2_summits(16415)" --samplesLabel "DMSO_IP_MBD3_rep123" "BPK25_IP_MBD3_rep123" "IAA_IP_MBD3_rep123"

## compute and plot use only MBD3 rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_MBD2_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig -R DMSO_IP_MBD2_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD3_samples_for_plotHeatMap.DMSO_MBD2_summits.gz
plotHeatmap -m bamCompare_matrix_of_MBD3_samples_for_plotHeatMap.DMSO_MBD2_summits.gz -o bamCompare_matrix_of_MBD3_samples_for_plotHeatMap.DMSO_MBD2_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_MBD2_summits" --regionsLabel "DMSO_MBD2_summits(16415)" --samplesLabel "DMSO_IP_MBD3_rep1" "DMSO_IP_MBD3_rep2" "DMSO_IP_MBD3_rep3" "BPK25_IP_MBD3_rep1" "BPK25_IP_MBD3_rep2" "BPK25_IP_MBD3_rep3" "IAA_IP_MBD3_rep1" "IAA_IP_MBD3_rep2" "IAA_IP_MBD3_rep3"


## compute and plot use only TET1 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_MBD2_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig -R DMSO_IP_MBD2_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_TET1_rep123_merged_for_plotHeatMap.DMSO_MBD2_summits.gz
plotHeatmap -m bamCompare_matrix_of_TET1_rep123_merged_for_plotHeatMap.DMSO_MBD2_summits.gz -o bamCompare_matrix_of_TET1_rep123_merged_for_plotHeatMap.DMSO_MBD2_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_MBD2_summits" --regionsLabel "DMSO_MBD2_summits(16415)" --samplesLabel "DMSO_IP_TET1_rep123" "BPK25_IP_TET1_rep123" "IAA_IP_TET1_rep123"

## compute and plot use only TET1 rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_MBD2_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig -R DMSO_IP_MBD2_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_TET1_samples_for_plotHeatMap.DMSO_MBD2_summits.gz
plotHeatmap -m bamCompare_matrix_of_TET1_samples_for_plotHeatMap.DMSO_MBD2_summits.gz -o bamCompare_matrix_of_TET1_samples_for_plotHeatMap.DMSO_MBD2_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_MBD2_summits" --regionsLabel "DMSO_MBD2_summits(16415)" --samplesLabel "DMSO_IP_TET1_rep1" "DMSO_IP_TET1_rep2" "DMSO_IP_TET1_rep3" "BPK25_IP_TET1_rep1" "BPK25_IP_TET1_rep2" "BPK25_IP_TET1_rep3" "IAA_IP_TET1_rep1" "IAA_IP_TET1_rep2" "IAA_IP_TET1_rep3"

```




# DMSO MBD3 summits f1k reads enrichment heatmap

## data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged
mkdir plotHeatMap_of_all_use_DMSO_MBD3_summits
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_MBD3_summits
# 所需的 summit bed 文件在如下分析文件夹中已经构造过了
# /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_2_for_DMSO_rep123_together
# 为所需的bed文件创建软连接创建 
ln -s ../macs2_peak_calling_2_for_DMSO_rep123_together/DMSO_IP_MBD3_rep123_TvsC/DMSO_IP_MBD3_rep123_TvsC_summits.bed3_uniq.sorted.bed DMSO_IP_MBD3_rep123_TvsC_summits.bed3_uniq.sorted.bed
wc -l *.sorted.bed
# 10498 DMSO_IP_MBD3_rep123_TvsC_summits.bed3_uniq.sorted.bed

# 创建全部所需的 bigwig 文件的软连接
# DMSO rep1/2/3 individual
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig
# DMSO rep123 merged
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig
ls -1 HEK293T_DMSO*.bigwig
ls -1 HEK293T_DMSO*.bigwig |wc -l
# 16

# BPK25 rep1/2/3 individual
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig
# BPK25 rep123 merged
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig
ls -1 HEK293T_BPK25*.bigwig
ls -1 HEK293T_BPK25*.bigwig |wc -l
# 16

# IAA rep1/2/3 individual
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig
# IAA rep123 merged
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig
ls -1 HEK293T_IAA*.bigwig
ls -1 HEK293T_IAA*.bigwig |wc -l
# 16

```

## computeMatrix and plotHeatmap use DMSO MBD3 summits

```bash
# plot use DMSO_MBD3_summits
## compute and plot use all 4 TFs rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_MBD3_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig -R DMSO_IP_MBD3_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.DMSO_MBD3_summits.gz
plotHeatmap -m bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.DMSO_MBD3_summits.gz -o bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.DMSO_MBD3_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_MBD3_summits" --regionsLabel "DMSO_MBD3_summits(10498)" --samplesLabel "DMSO_IP_CTCF_rep123" "DMSO_IP_MBD2_rep123" "DMSO_IP_MBD3_rep123" "DMSO_IP_TET1_rep123" "BPK25_IP_CTCF_rep123" "BPK25_IP_MBD2_rep123" "BPK25_IP_MBD3_rep123" "BPK25_IP_TET1_rep123" "IAA_IP_CTCF_rep123" "IAA_IP_MBD2_rep123" "IAA_IP_MBD3_rep123" "IAA_IP_TET1_rep123"

## compute and plot use all 4 TFs rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_MBD3_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig -R DMSO_IP_MBD3_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_all_samples_for_plotHeatMap.DMSO_MBD3_summits.gz
plotHeatmap -m bamCompare_matrix_of_all_samples_for_plotHeatMap.DMSO_MBD3_summits.gz -o bamCompare_matrix_of_all_samples_for_plotHeatMap.DMSO_MBD3_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_MBD3_summits" --regionsLabel "DMSO_MBD3_summits(10498)" --samplesLabel "DMSO_IP_CTCF_rep1" "DMSO_IP_CTCF_rep2" "DMSO_IP_CTCF_rep3" "DMSO_IP_MBD2_rep1" "DMSO_IP_MBD2_rep2" "DMSO_IP_MBD2_rep3" "DMSO_IP_MBD3_rep1" "DMSO_IP_MBD3_rep2" "DMSO_IP_MBD3_rep3" "DMSO_IP_TET1_rep1" "DMSO_IP_TET1_rep2" "DMSO_IP_TET1_rep3" "BPK25_IP_CTCF_rep1" "BPK25_IP_CTCF_rep2" "BPK25_IP_CTCF_rep3" "BPK25_IP_MBD2_rep1" "BPK25_IP_MBD2_rep2" "BPK25_IP_MBD2_rep3" "BPK25_IP_MBD3_rep1" "BPK25_IP_MBD3_rep2" "BPK25_IP_MBD3_rep3" "BPK25_IP_TET1_rep1" "BPK25_IP_TET1_rep2" "BPK25_IP_TET1_rep3" "IAA_IP_CTCF_rep1" "IAA_IP_CTCF_rep2" "IAA_IP_CTCF_rep3" "IAA_IP_MBD2_rep1" "IAA_IP_MBD2_rep2" "IAA_IP_MBD2_rep3" "IAA_IP_MBD3_rep1" "IAA_IP_MBD3_rep2" "IAA_IP_MBD3_rep3" "IAA_IP_TET1_rep1" "IAA_IP_TET1_rep2" "IAA_IP_TET1_rep3"


## compute and plot use only CTCF rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_MBD3_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig -R DMSO_IP_MBD3_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_CTCF_rep123_merged_for_plotHeatMap.DMSO_MBD3_summits.gz
plotHeatmap -m bamCompare_matrix_of_CTCF_rep123_merged_for_plotHeatMap.DMSO_MBD3_summits.gz -o bamCompare_matrix_of_CTCF_rep123_merged_for_plotHeatMap.DMSO_MBD3_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_MBD3_summits" --regionsLabel "DMSO_MBD3_summits(10498)" --samplesLabel "DMSO_IP_CTCF_rep123" "BPK25_IP_CTCF_rep123" "IAA_IP_CTCF_rep123"

## compute and plot use only CTCF rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_MBD3_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig -R DMSO_IP_MBD3_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_CTCF_samples_for_plotHeatMap.DMSO_MBD3_summits.gz
plotHeatmap -m bamCompare_matrix_of_CTCF_samples_for_plotHeatMap.DMSO_MBD3_summits.gz -o bamCompare_matrix_of_CTCF_samples_for_plotHeatMap.DMSO_MBD3_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_MBD3_summits" --regionsLabel "DMSO_MBD3_summits(10498)" --samplesLabel "DMSO_IP_CTCF_rep1" "DMSO_IP_CTCF_rep2" "DMSO_IP_CTCF_rep3" "BPK25_IP_CTCF_rep1" "BPK25_IP_CTCF_rep2" "BPK25_IP_CTCF_rep3" "IAA_IP_CTCF_rep1" "IAA_IP_CTCF_rep2" "IAA_IP_CTCF_rep3"


## compute and plot use only MBD2 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_MBD3_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig -R DMSO_IP_MBD3_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD2_rep123_merged_for_plotHeatMap.DMSO_MBD3_summits.gz
plotHeatmap -m bamCompare_matrix_of_MBD2_rep123_merged_for_plotHeatMap.DMSO_MBD3_summits.gz -o bamCompare_matrix_of_MBD2_rep123_merged_for_plotHeatMap.DMSO_MBD3_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_MBD3_summits" --regionsLabel "DMSO_MBD3_summits(10498)" --samplesLabel "DMSO_IP_MBD2_rep123" "BPK25_IP_MBD2_rep123" "IAA_IP_MBD2_rep123"

## compute and plot use only MBD2 rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_MBD3_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig -R DMSO_IP_MBD3_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD2_samples_for_plotHeatMap.DMSO_MBD3_summits.gz
plotHeatmap -m bamCompare_matrix_of_MBD2_samples_for_plotHeatMap.DMSO_MBD3_summits.gz -o bamCompare_matrix_of_MBD2_samples_for_plotHeatMap.DMSO_MBD3_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_MBD3_summits" --regionsLabel "DMSO_MBD3_summits(10498)" --samplesLabel "DMSO_IP_MBD2_rep1" "DMSO_IP_MBD2_rep2" "DMSO_IP_MBD2_rep3" "BPK25_IP_MBD2_rep1" "BPK25_IP_MBD2_rep2" "BPK25_IP_MBD2_rep3" "IAA_IP_MBD2_rep1" "IAA_IP_MBD2_rep2" "IAA_IP_MBD2_rep3"


## compute and plot use only MBD3 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_MBD3_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig -R DMSO_IP_MBD3_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD3_rep123_merged_for_plotHeatMap.DMSO_MBD3_summits.gz
plotHeatmap -m bamCompare_matrix_of_MBD3_rep123_merged_for_plotHeatMap.DMSO_MBD3_summits.gz -o bamCompare_matrix_of_MBD3_rep123_merged_for_plotHeatMap.DMSO_MBD3_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_MBD3_summits" --regionsLabel "DMSO_MBD3_summits(10498)" --samplesLabel "DMSO_IP_MBD3_rep123" "BPK25_IP_MBD3_rep123" "IAA_IP_MBD3_rep123"

## compute and plot use only MBD3 rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_MBD3_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig -R DMSO_IP_MBD3_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD3_samples_for_plotHeatMap.DMSO_MBD3_summits.gz
plotHeatmap -m bamCompare_matrix_of_MBD3_samples_for_plotHeatMap.DMSO_MBD3_summits.gz -o bamCompare_matrix_of_MBD3_samples_for_plotHeatMap.DMSO_MBD3_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_MBD3_summits" --regionsLabel "DMSO_MBD3_summits(10498)" --samplesLabel "DMSO_IP_MBD3_rep1" "DMSO_IP_MBD3_rep2" "DMSO_IP_MBD3_rep3" "BPK25_IP_MBD3_rep1" "BPK25_IP_MBD3_rep2" "BPK25_IP_MBD3_rep3" "IAA_IP_MBD3_rep1" "IAA_IP_MBD3_rep2" "IAA_IP_MBD3_rep3"


## compute and plot use only TET1 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_MBD3_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig -R DMSO_IP_MBD3_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_TET1_rep123_merged_for_plotHeatMap.DMSO_MBD3_summits.gz
plotHeatmap -m bamCompare_matrix_of_TET1_rep123_merged_for_plotHeatMap.DMSO_MBD3_summits.gz -o bamCompare_matrix_of_TET1_rep123_merged_for_plotHeatMap.DMSO_MBD3_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_MBD3_summits" --regionsLabel "DMSO_MBD3_summits(10498)" --samplesLabel "DMSO_IP_TET1_rep123" "BPK25_IP_TET1_rep123" "IAA_IP_TET1_rep123"

## compute and plot use only TET1 rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_MBD3_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig -R DMSO_IP_MBD3_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_TET1_samples_for_plotHeatMap.DMSO_MBD3_summits.gz
plotHeatmap -m bamCompare_matrix_of_TET1_samples_for_plotHeatMap.DMSO_MBD3_summits.gz -o bamCompare_matrix_of_TET1_samples_for_plotHeatMap.DMSO_MBD3_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_MBD3_summits" --regionsLabel "DMSO_MBD3_summits(10498)" --samplesLabel "DMSO_IP_TET1_rep1" "DMSO_IP_TET1_rep2" "DMSO_IP_TET1_rep3" "BPK25_IP_TET1_rep1" "BPK25_IP_TET1_rep2" "BPK25_IP_TET1_rep3" "IAA_IP_TET1_rep1" "IAA_IP_TET1_rep2" "IAA_IP_TET1_rep3"

```




# DMSO TET1 summits f1k reads enrichment heatmap

## data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged
mkdir plotHeatMap_of_all_use_DMSO_TET1_summits
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_TET1_summits
# 所需的 summit bed 文件在如下分析文件夹中已经构造过了
# /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/macs2_peak_calling_2_for_DMSO_rep123_together
# 为所需的bed文件创建软连接创建 
ln -s ../macs2_peak_calling_2_for_DMSO_rep123_together/DMSO_IP_TET1_rep123_TvsC/DMSO_IP_TET1_rep123_TvsC_summits.bed3_uniq.sorted.bed DMSO_IP_TET1_rep123_TvsC_summits.bed3_uniq.sorted.bed
wc -l *.sorted.bed
# 7515 DMSO_IP_TET1_rep123_TvsC_summits.bed3_uniq.sorted.bed

# 创建全部所需的 bigwig 文件的软连接
# DMSO rep1/2/3 individual
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig
# DMSO rep123 merged
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_DMSO_all_samples/HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig
ls -1 HEK293T_DMSO*.bigwig
ls -1 HEK293T_DMSO*.bigwig |wc -l
# 16

# BPK25 rep1/2/3 individual
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig
# BPK25 rep123 merged
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_BPK25_all_samples/HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig
ls -1 HEK293T_BPK25*.bigwig
ls -1 HEK293T_BPK25*.bigwig |wc -l
# 16

# IAA rep1/2/3 individual
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig
# IAA rep123 merged
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig
ln -s ../make_bamCompare_bigwig_for_IAA_all_samples/HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig
ls -1 HEK293T_IAA*.bigwig
ls -1 HEK293T_IAA*.bigwig |wc -l
# 16

```

## computeMatrix and plotHeatmap use DMSO TET1 summits

```bash
# plot use DMSO_TET1_summits
## compute and plot use all 4 TFs rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_TET1_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig -R DMSO_IP_TET1_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.DMSO_TET1_summits.gz
plotHeatmap -m bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.DMSO_TET1_summits.gz -o bamCompare_matrix_of_all_samples_rep123_merged_for_plotHeatMap.DMSO_TET1_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_TET1_summits" --regionsLabel "DMSO_TET1_summits(7515)" --samplesLabel "DMSO_IP_CTCF_rep123" "DMSO_IP_MBD2_rep123" "DMSO_IP_MBD3_rep123" "DMSO_IP_TET1_rep123" "BPK25_IP_CTCF_rep123" "BPK25_IP_MBD2_rep123" "BPK25_IP_MBD3_rep123" "BPK25_IP_TET1_rep123" "IAA_IP_CTCF_rep123" "IAA_IP_MBD2_rep123" "IAA_IP_MBD3_rep123" "IAA_IP_TET1_rep123"

## compute and plot use all 4 TFs rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_TET1_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig -R DMSO_IP_TET1_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_all_samples_for_plotHeatMap.DMSO_TET1_summits.gz
plotHeatmap -m bamCompare_matrix_of_all_samples_for_plotHeatMap.DMSO_TET1_summits.gz -o bamCompare_matrix_of_all_samples_for_plotHeatMap.DMSO_TET1_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_TET1_summits" --regionsLabel "DMSO_TET1_summits(7515)" --samplesLabel "DMSO_IP_CTCF_rep1" "DMSO_IP_CTCF_rep2" "DMSO_IP_CTCF_rep3" "DMSO_IP_MBD2_rep1" "DMSO_IP_MBD2_rep2" "DMSO_IP_MBD2_rep3" "DMSO_IP_MBD3_rep1" "DMSO_IP_MBD3_rep2" "DMSO_IP_MBD3_rep3" "DMSO_IP_TET1_rep1" "DMSO_IP_TET1_rep2" "DMSO_IP_TET1_rep3" "BPK25_IP_CTCF_rep1" "BPK25_IP_CTCF_rep2" "BPK25_IP_CTCF_rep3" "BPK25_IP_MBD2_rep1" "BPK25_IP_MBD2_rep2" "BPK25_IP_MBD2_rep3" "BPK25_IP_MBD3_rep1" "BPK25_IP_MBD3_rep2" "BPK25_IP_MBD3_rep3" "BPK25_IP_TET1_rep1" "BPK25_IP_TET1_rep2" "BPK25_IP_TET1_rep3" "IAA_IP_CTCF_rep1" "IAA_IP_CTCF_rep2" "IAA_IP_CTCF_rep3" "IAA_IP_MBD2_rep1" "IAA_IP_MBD2_rep2" "IAA_IP_MBD2_rep3" "IAA_IP_MBD3_rep1" "IAA_IP_MBD3_rep2" "IAA_IP_MBD3_rep3" "IAA_IP_TET1_rep1" "IAA_IP_TET1_rep2" "IAA_IP_TET1_rep3"


## compute and plot use only CTCF rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_TET1_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep123.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep123.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep123.TvsC.bigwig -R DMSO_IP_TET1_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_CTCF_rep123_merged_for_plotHeatMap.DMSO_TET1_summits.gz
plotHeatmap -m bamCompare_matrix_of_CTCF_rep123_merged_for_plotHeatMap.DMSO_TET1_summits.gz -o bamCompare_matrix_of_CTCF_rep123_merged_for_plotHeatMap.DMSO_TET1_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_TET1_summits" --regionsLabel "DMSO_TET1_summits(7515)" --samplesLabel "DMSO_IP_CTCF_rep123" "BPK25_IP_CTCF_rep123" "IAA_IP_CTCF_rep123"

## compute and plot use only CTCF rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_TET1_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_CTCF_rep1.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep2.TvsC.bigwig HEK293T_DMSO_IP_CTCF_rep3.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep1.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep2.TvsC.bigwig HEK293T_BPK25_IP_CTCF_rep3.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep1.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep2.TvsC.bigwig HEK293T_IAA_IP_CTCF_rep3.TvsC.bigwig -R DMSO_IP_TET1_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_CTCF_samples_for_plotHeatMap.DMSO_TET1_summits.gz
plotHeatmap -m bamCompare_matrix_of_CTCF_samples_for_plotHeatMap.DMSO_TET1_summits.gz -o bamCompare_matrix_of_CTCF_samples_for_plotHeatMap.DMSO_TET1_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_TET1_summits" --regionsLabel "DMSO_TET1_summits(7515)" --samplesLabel "DMSO_IP_CTCF_rep1" "DMSO_IP_CTCF_rep2" "DMSO_IP_CTCF_rep3" "BPK25_IP_CTCF_rep1" "BPK25_IP_CTCF_rep2" "BPK25_IP_CTCF_rep3" "IAA_IP_CTCF_rep1" "IAA_IP_CTCF_rep2" "IAA_IP_CTCF_rep3"


## compute and plot use only MBD2 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_TET1_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD2_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep123.TvsC.bigwig -R DMSO_IP_TET1_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD2_rep123_merged_for_plotHeatMap.DMSO_TET1_summits.gz
plotHeatmap -m bamCompare_matrix_of_MBD2_rep123_merged_for_plotHeatMap.DMSO_TET1_summits.gz -o bamCompare_matrix_of_MBD2_rep123_merged_for_plotHeatMap.DMSO_TET1_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_TET1_summits" --regionsLabel "DMSO_TET1_summits(7515)" --samplesLabel "DMSO_IP_MBD2_rep123" "BPK25_IP_MBD2_rep123" "IAA_IP_MBD2_rep123"

## compute and plot use only MBD2 rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_TET1_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD2_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD2_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD2_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD2_rep3.TvsC.bigwig -R DMSO_IP_TET1_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD2_samples_for_plotHeatMap.DMSO_TET1_summits.gz
plotHeatmap -m bamCompare_matrix_of_MBD2_samples_for_plotHeatMap.DMSO_TET1_summits.gz -o bamCompare_matrix_of_MBD2_samples_for_plotHeatMap.DMSO_TET1_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_TET1_summits" --regionsLabel "DMSO_TET1_summits(7515)" --samplesLabel "DMSO_IP_MBD2_rep1" "DMSO_IP_MBD2_rep2" "DMSO_IP_MBD2_rep3" "BPK25_IP_MBD2_rep1" "BPK25_IP_MBD2_rep2" "BPK25_IP_MBD2_rep3" "IAA_IP_MBD2_rep1" "IAA_IP_MBD2_rep2" "IAA_IP_MBD2_rep3"


## compute and plot use only MBD3 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_TET1_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD3_rep123.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep123.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep123.TvsC.bigwig -R DMSO_IP_TET1_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD3_rep123_merged_for_plotHeatMap.DMSO_TET1_summits.gz
plotHeatmap -m bamCompare_matrix_of_MBD3_rep123_merged_for_plotHeatMap.DMSO_TET1_summits.gz -o bamCompare_matrix_of_MBD3_rep123_merged_for_plotHeatMap.DMSO_TET1_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_TET1_summits" --regionsLabel "DMSO_TET1_summits(7515)" --samplesLabel "DMSO_IP_MBD3_rep123" "BPK25_IP_MBD3_rep123" "IAA_IP_MBD3_rep123"

## compute and plot use only MBD3 rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_TET1_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_MBD3_rep1.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep2.TvsC.bigwig HEK293T_DMSO_IP_MBD3_rep3.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep1.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep2.TvsC.bigwig HEK293T_BPK25_IP_MBD3_rep3.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep1.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep2.TvsC.bigwig HEK293T_IAA_IP_MBD3_rep3.TvsC.bigwig -R DMSO_IP_TET1_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_MBD3_samples_for_plotHeatMap.DMSO_TET1_summits.gz
plotHeatmap -m bamCompare_matrix_of_MBD3_samples_for_plotHeatMap.DMSO_TET1_summits.gz -o bamCompare_matrix_of_MBD3_samples_for_plotHeatMap.DMSO_TET1_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_TET1_summits" --regionsLabel "DMSO_TET1_summits(7515)" --samplesLabel "DMSO_IP_MBD3_rep1" "DMSO_IP_MBD3_rep2" "DMSO_IP_MBD3_rep3" "BPK25_IP_MBD3_rep1" "BPK25_IP_MBD3_rep2" "BPK25_IP_MBD3_rep3" "IAA_IP_MBD3_rep1" "IAA_IP_MBD3_rep2" "IAA_IP_MBD3_rep3"


## compute and plot use only TET1 rep123 merged bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_TET1_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_TET1_rep123.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep123.TvsC.bigwig HEK293T_IAA_IP_TET1_rep123.TvsC.bigwig -R DMSO_IP_TET1_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_TET1_rep123_merged_for_plotHeatMap.DMSO_TET1_summits.gz
plotHeatmap -m bamCompare_matrix_of_TET1_rep123_merged_for_plotHeatMap.DMSO_TET1_summits.gz -o bamCompare_matrix_of_TET1_rep123_merged_for_plotHeatMap.DMSO_TET1_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_TET1_summits" --regionsLabel "DMSO_TET1_summits(7515)" --samplesLabel "DMSO_IP_TET1_rep123" "BPK25_IP_TET1_rep123" "IAA_IP_TET1_rep123"

## compute and plot use only TET1 rep1/2/3 respective bigwig
cd /home/sunwenju/projects_on_940xa/2023112201_WN_IAA_BPK25_ChIP_Seq_DMSO_part12_merged/plotHeatMap_of_all_use_DMSO_TET1_summits
computeMatrix reference-point -p 32 --referencePoint center -S HEK293T_DMSO_IP_TET1_rep1.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep2.TvsC.bigwig HEK293T_DMSO_IP_TET1_rep3.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep1.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep2.TvsC.bigwig HEK293T_BPK25_IP_TET1_rep3.TvsC.bigwig HEK293T_IAA_IP_TET1_rep1.TvsC.bigwig HEK293T_IAA_IP_TET1_rep2.TvsC.bigwig HEK293T_IAA_IP_TET1_rep3.TvsC.bigwig -R DMSO_IP_TET1_rep123_TvsC_summits.bed3_uniq.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_TET1_samples_for_plotHeatMap.DMSO_TET1_summits.gz
plotHeatmap -m bamCompare_matrix_of_TET1_samples_for_plotHeatMap.DMSO_TET1_summits.gz -o bamCompare_matrix_of_TET1_samples_for_plotHeatMap.DMSO_TET1_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to DMSO_TET1_summits" --regionsLabel "DMSO_TET1_summits(7515)" --samplesLabel "DMSO_IP_TET1_rep1" "DMSO_IP_TET1_rep2" "DMSO_IP_TET1_rep3" "BPK25_IP_TET1_rep1" "BPK25_IP_TET1_rep2" "BPK25_IP_TET1_rep3" "IAA_IP_TET1_rep1" "IAA_IP_TET1_rep2" "IAA_IP_TET1_rep3"

```

