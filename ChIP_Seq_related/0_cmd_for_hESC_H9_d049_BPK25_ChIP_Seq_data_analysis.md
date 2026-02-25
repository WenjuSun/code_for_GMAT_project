
# QC

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
# 
bash 1_WN_H9_d049_BPK25_input_raw_fq_part12_merge.sh
#
mkdir 1_all_H9_d049_BPK25_ChIP_Seq_raw_data_fastqc_result
nohup fastqc -outdir 1_all_H9_d049_BPK25_ChIP_Seq_raw_data_fastqc_result --threads 72 *rep[123].R[12].fq.gz >1_all_H9_d049_BPK25_ChIP_Seq_raw_data_fastqc.runlog 2>&1 &
tail -f 1_all_H9_d049_BPK25_ChIP_Seq_raw_data_fastqc.runlog

# 并用 multiqc 进行结果汇总
cd 1_all_H9_d049_BPK25_ChIP_Seq_raw_data_fastqc_result
multiqc ./*fastqc.zip --ignore *.html --no-megaqc-upload --interactive --cl_config "max_table_rows: 3000"
cd ../

```

# cutadapt and QC

```bash
# 根据 fastqc 的结果提示，在 illumina 的接头序列文件中找到如下序列，并对其中一个文件进行了grep统计，确认是接头
# Read 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
# Read 2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
# 使用如下命令生成去接头脚本并运行
for i in `ls -1 *rep[123].R1.fq.gz`;do echo 'cutadapt -j 72 -q 20,20 --trim-n -m 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o '${i%.R1.fq.gz}.R1.cutadapt.fq.gz' -p '${i%.R1.fq.gz}.R2.cutadapt.fq.gz' '$i' '${i%.R1.fq.gz}.R2.fq.gz;done >2_all_H9_d049_BPK25_ChIP_Seq_raw_data_run_cutadapt.sh
nohup bash 2_all_H9_d049_BPK25_ChIP_Seq_raw_data_run_cutadapt.sh >2_all_H9_d049_BPK25_ChIP_Seq_raw_data_run_cutadapt.sh.runlog 2>&1 &
tail -f 2_all_H9_d049_BPK25_ChIP_Seq_raw_data_run_cutadapt.sh.runlog

# 去接头后再次进行质量评估
mkdir 2_all_H9_d049_BPK25_ChIP_Seq_cutadapt_fastqc_result
nohup fastqc -outdir 2_all_H9_d049_BPK25_ChIP_Seq_cutadapt_fastqc_result --threads 72 *.cutadapt.fq.gz >2_all_H9_d049_BPK25_ChIP_Seq_cutadapt_fastqc.runlog 2>&1 &
tail -f 2_all_H9_d049_BPK25_ChIP_Seq_cutadapt_fastqc.runlog

# 并用 multiqc 进行结果汇总
# 提取原始数据reads数量
cd 2_all_H9_d049_BPK25_ChIP_Seq_cutadapt_fastqc_result
multiqc ./*fastqc.zip --ignore *.html --no-megaqc-upload --interactive --cl_config "max_table_rows: 3000"
cd ../

```



# Genome mapping

```bash
# 由于测试后发现 将全部 bowtie2 日志给 multiqc 汇总，只汇总了最后一个
# 将 bowtie2 日志拆分，发现可以汇总全部，所以需要给每条bowtie2的比对命令生成一个日志，使用通配符给multiqc，才能将全部比对日志汇总到一起
# 同时将没有比对上的 reads 分别输出到文件，如果比对效率低，可对这些文件进行 collapser 用top序列blast 分析原因，如果比对效率很好90%，可将这些文件删除。
# 可使用如下命令构造比对脚本， 然后运行
for i in `ls -1 *rep[123].R1.cutadapt.fq.gz`;do echo 'bowtie2 -t -p 72 --no-unal -x /home/sunwenju/3_genomeData/human/hg19/bowtie2Index/hg19bt2idx -1 '$i' -2 '${i%.R1.cutadapt.fq.gz}.R2.cutadapt.fq.gz' -S '${i%.R1.cutadapt.fq.gz}.bowtie2.sam' --un-gz '${i%.R1.cutadapt.fq.gz}.notAlign.unpairedReads.fq.gz' --un-conc-gz '${i%.R1.cutadapt.fq.gz}.notAlign.pairedReads.fq.gz' >'${i%.R1.cutadapt.fq.gz}.bowtie2.runlog' 2>&1';done >3_all_H9_d049_BPK25_ChIP_Seq_cutadapt_data_run_bowtie2.sh
nohup bash 3_all_H9_d049_BPK25_ChIP_Seq_cutadapt_data_run_bowtie2.sh >3_all_H9_d049_BPK25_ChIP_Seq_cutadapt_data_run_bowtie2.sh.runlog 2>&1 &
tail -f 3_all_H9_d049_BPK25_ChIP_Seq_cutadapt_data_run_bowtie2.sh.runlog

# 使用 multiqc 对 bowtie2 比对情况进行汇总
multiqc -m bowtie2 -o 3_all_H9_d049_BPK25_ChIP_Seq_bowtie2_map_results_summary ./*.bowtie2.runlog --no-megaqc-upload --interactive

# 比对效率过低时，统计汇总没有比对上的reads，数量最多的在前面，可用来查看是什么序列
for i in `ls -1 *notAlign.pairedReads.fq.1.gz`;do echo $i; zcat $i | fastx_collapser -Q 33 -v -o ${i%.fq.1.gz}.1.fa;done
for i in `ls -1 *notAlign.pairedReads.fq.2.gz`;do echo $i; zcat $i | fastx_collapser -Q 33 -v -o ${i%.fq.2.gz}.2.fa;done
for i in `ls -1 *notAlign.unpairedReads.fq.gz`;do echo $i; zcat $i | fastx_collapser -Q 33 -v -o ${i%.cutadapt.fq.gz}.fa;done

# 比对效率很好时，collapser 之后可以删除这些文件
rm *notAlign.pairedReads.fq.1.gz
rm *notAlign.pairedReads.fq.2.gz
rm *notAlign.unpairedReads.fq.gz

```

# Convert sam to sorted bam

```bash
# 将比对sam结果转换为 按座位排序的 bam 并建立索引
# for i in `ls -1 *.bowtie2.sam|sort`; do echo $i; samtools view --threads 36 -bS $i | samtools sort --threads 36 -o ${i%.sam}.sorted.bam; samtools index ${i%.sam}.sorted.bam; done
nohup bash 4_all_H9_d049_BPK25_ChIP_Seq_bowtie2_sam_to_sorted_bam.sh >4_all_H9_d049_BPK25_ChIP_Seq_bowtie2_sam_to_sorted_bam.sh.runlog 2>&1 &
tail -f 4_all_H9_d049_BPK25_ChIP_Seq_bowtie2_sam_to_sorted_bam.sh.runlog
# rm *.bowtie2.sam

# 使用samtools flagstat 统计比对情况
# 每个bam文件不到一分钟
for i in `ls -1 *bowtie2.sorted.bam`;do echo $i; samtools flagstat --threads 16 $i;done >4_all_H9_d049_BPK25_ChIP_Seq_bowtie2_map_samtools_flagstat.txt
# 统计比对到各个染色体上的 总read数量 
for i in `ls -1 *sorted.bam|sort`;do echo $i; samtools view $i |grep -v -e '^@' |cut -f 3 |sort |uniq -c |awk 'BEGIN{OFS="\t"}{if($0 !~ /_/){print $2,$1}}' |sort -k1,1 >${i%.sam}.chr_mapped_reads_stat.txt;done
# 统计比对到各个染色体上的 uniq read数量 
for i in `ls -1 *sorted.bam|sort`;do echo $i; samtools view $i |grep -v -e '^@' |cut -f 3,4 |sort |uniq |cut -f 1 |sort |uniq -c |awk 'BEGIN{OFS="\t"}{if($0 !~ /_/){print $2,$1}}' |sort -k1,1 >${i%.sam}.chr_mapped_uniq_reads_stat.txt;done

# 比对到各个染色体上的 read数量 统计结果合并
paste `ls -1 *chr_mapped_reads_stat.txt|sort` |head
# 检测各个统计结果的染色体数量和顺序是否一致
paste `ls -1 *chr_mapped_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"}{if($1!=$3 || $1!=$5 || $1!=$7 || $1!=$9 || $1!=$11 || $1!=$13 || $1!=$15 || $1!=$17 || $1!=$19 || $1!=$21 || $1!=$23 || $1!=$25 || $1!=$27 || $1!=$29 || $1!=$31 || $1!=$33 || $1!=$35 || $1!=$37 || $1!=$39 || $1!=$41 || $1!=$43 || $1!=$45 || $1!=$47 || $1!=$49 || $1!=$51 || $1!=$53 || $1!=$55 || $1!=$57 || $1!=$59 || $1!=$61 || $1!=$63 || $1!=$65 || $1!=$67 || $1!=$69 || $1!=$71 || $1!=$73 || $1!=$75 || $1!=$77 || $1!=$79 || $1!=$81 || $1!=$83 || $1!=$85 || $1!=$87 || $1!=$89 || $1!=$91 || $1!=$93 || $1!=$95 || $1!=$97 || $1!=$99 || $1!=$101 || $1!=$103 || $1!=$105 || $1!=$107 || $1!=$109 || $1!=$111 || $1!=$113 || $1!=$115 || $1!=$117 || $1!=$119 || $1!=$121 || $1!=$123 || $1!=$125 || $1!=$127 || $1!=$129 || $1!=$131 || $1!=$133 || $1!=$135 || $1!=$137 || $1!=$139 || $1!=$141 || $1!=$143 || $1!=$145 || $1!=$147 || $1!=$149 || $1!=$151 || $1!=$153 || $1!=$155 || $1!=$157 || $1!=$159 || $1!=$161 || $1!=$163 || $1!=$165 || $1!=$167 || $1!=$169 || $1!=$171 || $1!=$173 || $1!=$175 || $1!=$177 || $1!=$179 || $1!=$181 || $1!=$183 || $1!=$185 || $1!=$187 || $1!=$189 || $1!=$191 || $1!=$193 || $1!=$195 || $1!=$197 || $1!=$199 || $1!=$201 || $1!=$203 || $1!=$205 || $1!=$207 || $1!=$209 || $1!=$211 || $1!=$213 || $1!=$215 ){print $1,$2}}' |head
# 染色体数量和顺序一致，则合并数据列
paste `ls -1 *chr_mapped_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"; print "chrom","H9_d0_BPK25_input_rep1","H9_d0_BPK25_input_rep2","H9_d0_BPK25_input_rep3","H9_d0_BPK25_IP_CTCF_rep1","H9_d0_BPK25_IP_CTCF_rep2","H9_d0_BPK25_IP_CTCF_rep3","H9_d0_BPK25_IP_MBD3_rep1","H9_d0_BPK25_IP_MBD3_rep2","H9_d0_BPK25_IP_MBD3_rep3","H9_d0_BPK25_IP_TET1_rep1","H9_d0_BPK25_IP_TET1_rep2","H9_d0_BPK25_IP_TET1_rep3","H9_d0_BPK25_IP_TET2_rep1","H9_d0_BPK25_IP_TET2_rep2","H9_d0_BPK25_IP_TET2_rep3","H9_d0_BPK25_IP_TET3_rep1","H9_d0_BPK25_IP_TET3_rep2","H9_d0_BPK25_IP_TET3_rep3","H9_d0_DMSO_input_rep1","H9_d0_DMSO_input_rep2","H9_d0_DMSO_input_rep3","H9_d0_DMSO_IP_CTCF_rep1","H9_d0_DMSO_IP_CTCF_rep2","H9_d0_DMSO_IP_CTCF_rep3","H9_d0_DMSO_IP_MBD3_rep1","H9_d0_DMSO_IP_MBD3_rep2","H9_d0_DMSO_IP_MBD3_rep3","H9_d0_DMSO_IP_TET1_rep1","H9_d0_DMSO_IP_TET1_rep2","H9_d0_DMSO_IP_TET1_rep3","H9_d0_DMSO_IP_TET2_rep1","H9_d0_DMSO_IP_TET2_rep2","H9_d0_DMSO_IP_TET2_rep3","H9_d0_DMSO_IP_TET3_rep1","H9_d0_DMSO_IP_TET3_rep2","H9_d0_DMSO_IP_TET3_rep3","H9_d4_BPK25_input_rep1","H9_d4_BPK25_input_rep2","H9_d4_BPK25_input_rep3","H9_d4_BPK25_IP_CTCF_rep1","H9_d4_BPK25_IP_CTCF_rep2","H9_d4_BPK25_IP_CTCF_rep3","H9_d4_BPK25_IP_MBD3_rep1","H9_d4_BPK25_IP_MBD3_rep2","H9_d4_BPK25_IP_MBD3_rep3","H9_d4_BPK25_IP_TET1_rep1","H9_d4_BPK25_IP_TET1_rep2","H9_d4_BPK25_IP_TET1_rep3","H9_d4_BPK25_IP_TET2_rep1","H9_d4_BPK25_IP_TET2_rep2","H9_d4_BPK25_IP_TET2_rep3","H9_d4_BPK25_IP_TET3_rep1","H9_d4_BPK25_IP_TET3_rep2","H9_d4_BPK25_IP_TET3_rep3","H9_d4_DMSO_input_rep1","H9_d4_DMSO_input_rep2","H9_d4_DMSO_input_rep3","H9_d4_DMSO_IP_CTCF_rep1","H9_d4_DMSO_IP_CTCF_rep2","H9_d4_DMSO_IP_CTCF_rep3","H9_d4_DMSO_IP_MBD3_rep1","H9_d4_DMSO_IP_MBD3_rep2","H9_d4_DMSO_IP_MBD3_rep3","H9_d4_DMSO_IP_TET1_rep1","H9_d4_DMSO_IP_TET1_rep2","H9_d4_DMSO_IP_TET1_rep3","H9_d4_DMSO_IP_TET2_rep1","H9_d4_DMSO_IP_TET2_rep2","H9_d4_DMSO_IP_TET2_rep3","H9_d4_DMSO_IP_TET3_rep1","H9_d4_DMSO_IP_TET3_rep2","H9_d4_DMSO_IP_TET3_rep3","H9_d9_BPK25_input_rep1","H9_d9_BPK25_input_rep2","H9_d9_BPK25_input_rep3","H9_d9_BPK25_IP_CTCF_rep1","H9_d9_BPK25_IP_CTCF_rep2","H9_d9_BPK25_IP_CTCF_rep3","H9_d9_BPK25_IP_MBD3_rep1","H9_d9_BPK25_IP_MBD3_rep2","H9_d9_BPK25_IP_MBD3_rep3","H9_d9_BPK25_IP_TET1_rep1","H9_d9_BPK25_IP_TET1_rep2","H9_d9_BPK25_IP_TET1_rep3","H9_d9_BPK25_IP_TET2_rep1","H9_d9_BPK25_IP_TET2_rep2","H9_d9_BPK25_IP_TET2_rep3","H9_d9_BPK25_IP_TET3_rep1","H9_d9_BPK25_IP_TET3_rep2","H9_d9_BPK25_IP_TET3_rep3","H9_d9_DMSO_input_rep1","H9_d9_DMSO_input_rep2","H9_d9_DMSO_input_rep3","H9_d9_DMSO_IP_CTCF_rep1","H9_d9_DMSO_IP_CTCF_rep2","H9_d9_DMSO_IP_CTCF_rep3","H9_d9_DMSO_IP_MBD3_rep1","H9_d9_DMSO_IP_MBD3_rep2","H9_d9_DMSO_IP_MBD3_rep3","H9_d9_DMSO_IP_TET1_rep1","H9_d9_DMSO_IP_TET1_rep2","H9_d9_DMSO_IP_TET1_rep3","H9_d9_DMSO_IP_TET2_rep1","H9_d9_DMSO_IP_TET2_rep2","H9_d9_DMSO_IP_TET2_rep3","H9_d9_DMSO_IP_TET3_rep1","H9_d9_DMSO_IP_TET3_rep2","H9_d9_DMSO_IP_TET3_rep3"}{if($1!=$3 && $1!=$5 && $1!=$7 && $1!=$9 && $1!=$11 && $1!=$13 && $1!=$15 && $1!=$17 && $1!=$19 && $1!=$21 && $1!=$23 && $1!=$25 && $1!=$27 && $1!=$29 && $1!=$31 && $1!=$33 && $1!=$35 && $1!=$37 && $1!=$39 && $1!=$41 && $1!=$43 && $1!=$45 && $1!=$47 && $1!=$49 && $1!=$51 && $1!=$53 && $1!=$55 && $1!=$57 && $1!=$59 && $1!=$61 && $1!=$63 && $1!=$65 && $1!=$67 && $1!=$69 && $1!=$71 && $1!=$73 && $1!=$75 && $1!=$77 && $1!=$79 && $1!=$81 && $1!=$83 && $1!=$85 && $1!=$87 && $1!=$89 && $1!=$91 && $1!=$93 && $1!=$95 && $1!=$97 && $1!=$99 && $1!=$101 && $1!=$103 && $1!=$105 && $1!=$107 && $1!=$109 && $1!=$111 && $1!=$113 && $1!=$115 && $1!=$117 && $1!=$119 && $1!=$121 && $1!=$123 && $1!=$125 && $1!=$127 && $1!=$129 && $1!=$131 && $1!=$133 && $1!=$135 && $1!=$137 && $1!=$139 && $1!=$141 && $1!=$143 && $1!=$145 && $1!=$147 && $1!=$149 && $1!=$151 && $1!=$153 && $1!=$155 && $1!=$157 && $1!=$159 && $1!=$161 && $1!=$163 && $1!=$165 && $1!=$167 && $1!=$169 && $1!=$171 && $1!=$173 && $1!=$175 && $1!=$177 && $1!=$179 && $1!=$181 && $1!=$183 && $1!=$185 && $1!=$187 && $1!=$189 && $1!=$191 && $1!=$193 && $1!=$195 && $1!=$197 && $1!=$199 && $1!=$201 && $1!=$203 && $1!=$205 && $1!=$207 && $1!=$209 && $1!=$211 && $1!=$213 && $1!=$215){print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36,$38,$40,$42,$44,$46,$48,$50,$52,$54,$56,$58,$60,$62,$64,$66,$68,$70,$72,$74,$76,$78,$80,$82,$84,$86,$88,$90,$92,$94,$96,$98,$100,$102,$104,$106,$108,$110,$112,$114,$116,$118,$120,$122,$124,$126,$128,$130,$132,$134,$136,$138,$140,$142,$144,$146,$148,$150,$152,$154,$156,$158,$160,$162,$164,$166,$168,$170,$172,$174,$176,$178,$180,$182,$184,$186,$188,$190,$192,$194,$196,$198,$200,$202,$204,$206,$208,$210,$212,$214,$216}}' |head
paste `ls -1 *chr_mapped_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"; print "chrom","H9_d0_BPK25_input_rep1","H9_d0_BPK25_input_rep2","H9_d0_BPK25_input_rep3","H9_d0_BPK25_IP_CTCF_rep1","H9_d0_BPK25_IP_CTCF_rep2","H9_d0_BPK25_IP_CTCF_rep3","H9_d0_BPK25_IP_MBD3_rep1","H9_d0_BPK25_IP_MBD3_rep2","H9_d0_BPK25_IP_MBD3_rep3","H9_d0_BPK25_IP_TET1_rep1","H9_d0_BPK25_IP_TET1_rep2","H9_d0_BPK25_IP_TET1_rep3","H9_d0_BPK25_IP_TET2_rep1","H9_d0_BPK25_IP_TET2_rep2","H9_d0_BPK25_IP_TET2_rep3","H9_d0_BPK25_IP_TET3_rep1","H9_d0_BPK25_IP_TET3_rep2","H9_d0_BPK25_IP_TET3_rep3","H9_d0_DMSO_input_rep1","H9_d0_DMSO_input_rep2","H9_d0_DMSO_input_rep3","H9_d0_DMSO_IP_CTCF_rep1","H9_d0_DMSO_IP_CTCF_rep2","H9_d0_DMSO_IP_CTCF_rep3","H9_d0_DMSO_IP_MBD3_rep1","H9_d0_DMSO_IP_MBD3_rep2","H9_d0_DMSO_IP_MBD3_rep3","H9_d0_DMSO_IP_TET1_rep1","H9_d0_DMSO_IP_TET1_rep2","H9_d0_DMSO_IP_TET1_rep3","H9_d0_DMSO_IP_TET2_rep1","H9_d0_DMSO_IP_TET2_rep2","H9_d0_DMSO_IP_TET2_rep3","H9_d0_DMSO_IP_TET3_rep1","H9_d0_DMSO_IP_TET3_rep2","H9_d0_DMSO_IP_TET3_rep3","H9_d4_BPK25_input_rep1","H9_d4_BPK25_input_rep2","H9_d4_BPK25_input_rep3","H9_d4_BPK25_IP_CTCF_rep1","H9_d4_BPK25_IP_CTCF_rep2","H9_d4_BPK25_IP_CTCF_rep3","H9_d4_BPK25_IP_MBD3_rep1","H9_d4_BPK25_IP_MBD3_rep2","H9_d4_BPK25_IP_MBD3_rep3","H9_d4_BPK25_IP_TET1_rep1","H9_d4_BPK25_IP_TET1_rep2","H9_d4_BPK25_IP_TET1_rep3","H9_d4_BPK25_IP_TET2_rep1","H9_d4_BPK25_IP_TET2_rep2","H9_d4_BPK25_IP_TET2_rep3","H9_d4_BPK25_IP_TET3_rep1","H9_d4_BPK25_IP_TET3_rep2","H9_d4_BPK25_IP_TET3_rep3","H9_d4_DMSO_input_rep1","H9_d4_DMSO_input_rep2","H9_d4_DMSO_input_rep3","H9_d4_DMSO_IP_CTCF_rep1","H9_d4_DMSO_IP_CTCF_rep2","H9_d4_DMSO_IP_CTCF_rep3","H9_d4_DMSO_IP_MBD3_rep1","H9_d4_DMSO_IP_MBD3_rep2","H9_d4_DMSO_IP_MBD3_rep3","H9_d4_DMSO_IP_TET1_rep1","H9_d4_DMSO_IP_TET1_rep2","H9_d4_DMSO_IP_TET1_rep3","H9_d4_DMSO_IP_TET2_rep1","H9_d4_DMSO_IP_TET2_rep2","H9_d4_DMSO_IP_TET2_rep3","H9_d4_DMSO_IP_TET3_rep1","H9_d4_DMSO_IP_TET3_rep2","H9_d4_DMSO_IP_TET3_rep3","H9_d9_BPK25_input_rep1","H9_d9_BPK25_input_rep2","H9_d9_BPK25_input_rep3","H9_d9_BPK25_IP_CTCF_rep1","H9_d9_BPK25_IP_CTCF_rep2","H9_d9_BPK25_IP_CTCF_rep3","H9_d9_BPK25_IP_MBD3_rep1","H9_d9_BPK25_IP_MBD3_rep2","H9_d9_BPK25_IP_MBD3_rep3","H9_d9_BPK25_IP_TET1_rep1","H9_d9_BPK25_IP_TET1_rep2","H9_d9_BPK25_IP_TET1_rep3","H9_d9_BPK25_IP_TET2_rep1","H9_d9_BPK25_IP_TET2_rep2","H9_d9_BPK25_IP_TET2_rep3","H9_d9_BPK25_IP_TET3_rep1","H9_d9_BPK25_IP_TET3_rep2","H9_d9_BPK25_IP_TET3_rep3","H9_d9_DMSO_input_rep1","H9_d9_DMSO_input_rep2","H9_d9_DMSO_input_rep3","H9_d9_DMSO_IP_CTCF_rep1","H9_d9_DMSO_IP_CTCF_rep2","H9_d9_DMSO_IP_CTCF_rep3","H9_d9_DMSO_IP_MBD3_rep1","H9_d9_DMSO_IP_MBD3_rep2","H9_d9_DMSO_IP_MBD3_rep3","H9_d9_DMSO_IP_TET1_rep1","H9_d9_DMSO_IP_TET1_rep2","H9_d9_DMSO_IP_TET1_rep3","H9_d9_DMSO_IP_TET2_rep1","H9_d9_DMSO_IP_TET2_rep2","H9_d9_DMSO_IP_TET2_rep3","H9_d9_DMSO_IP_TET3_rep1","H9_d9_DMSO_IP_TET3_rep2","H9_d9_DMSO_IP_TET3_rep3"}{if($1!=$3 && $1!=$5 && $1!=$7 && $1!=$9 && $1!=$11 && $1!=$13 && $1!=$15 && $1!=$17 && $1!=$19 && $1!=$21 && $1!=$23 && $1!=$25 && $1!=$27 && $1!=$29 && $1!=$31 && $1!=$33 && $1!=$35 && $1!=$37 && $1!=$39 && $1!=$41 && $1!=$43 && $1!=$45 && $1!=$47 && $1!=$49 && $1!=$51 && $1!=$53 && $1!=$55 && $1!=$57 && $1!=$59 && $1!=$61 && $1!=$63 && $1!=$65 && $1!=$67 && $1!=$69 && $1!=$71 && $1!=$73 && $1!=$75 && $1!=$77 && $1!=$79 && $1!=$81 && $1!=$83 && $1!=$85 && $1!=$87 && $1!=$89 && $1!=$91 && $1!=$93 && $1!=$95 && $1!=$97 && $1!=$99 && $1!=$101 && $1!=$103 && $1!=$105 && $1!=$107 && $1!=$109 && $1!=$111 && $1!=$113 && $1!=$115 && $1!=$117 && $1!=$119 && $1!=$121 && $1!=$123 && $1!=$125 && $1!=$127 && $1!=$129 && $1!=$131 && $1!=$133 && $1!=$135 && $1!=$137 && $1!=$139 && $1!=$141 && $1!=$143 && $1!=$145 && $1!=$147 && $1!=$149 && $1!=$151 && $1!=$153 && $1!=$155 && $1!=$157 && $1!=$159 && $1!=$161 && $1!=$163 && $1!=$165 && $1!=$167 && $1!=$169 && $1!=$171 && $1!=$173 && $1!=$175 && $1!=$177 && $1!=$179 && $1!=$181 && $1!=$183 && $1!=$185 && $1!=$187 && $1!=$189 && $1!=$191 && $1!=$193 && $1!=$195 && $1!=$197 && $1!=$199 && $1!=$201 && $1!=$203 && $1!=$205 && $1!=$207 && $1!=$209 && $1!=$211 && $1!=$213 && $1!=$215){print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36,$38,$40,$42,$44,$46,$48,$50,$52,$54,$56,$58,$60,$62,$64,$66,$68,$70,$72,$74,$76,$78,$80,$82,$84,$86,$88,$90,$92,$94,$96,$98,$100,$102,$104,$106,$108,$110,$112,$114,$116,$118,$120,$122,$124,$126,$128,$130,$132,$134,$136,$138,$140,$142,$144,$146,$148,$150,$152,$154,$156,$158,$160,$162,$164,$166,$168,$170,$172,$174,$176,$178,$180,$182,$184,$186,$188,$190,$192,$194,$196,$198,$200,$202,$204,$206,$208,$210,$212,$214,$216}}' >all_samples_chr_mapped_reads_stat.summary.txt

# 比对到各个染色体上的 uniq read数量 统计结果合并
paste `ls -1 *chr_mapped_uniq_reads_stat.txt|sort` |head
# 检测各个统计结果的染色体数量和顺序是否一致
paste `ls -1 *chr_mapped_uniq_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"}{if($1!=$3 || $1!=$5 || $1!=$7 || $1!=$9 || $1!=$11 || $1!=$13 || $1!=$15 || $1!=$17 || $1!=$19 || $1!=$21 || $1!=$23 || $1!=$25 || $1!=$27 || $1!=$29 || $1!=$31 || $1!=$33 || $1!=$35 || $1!=$37 || $1!=$39 || $1!=$41 || $1!=$43 || $1!=$45 || $1!=$47 || $1!=$49 || $1!=$51 || $1!=$53 || $1!=$55 || $1!=$57 || $1!=$59 || $1!=$61 || $1!=$63 || $1!=$65 || $1!=$67 || $1!=$69 || $1!=$71 || $1!=$73 || $1!=$75 || $1!=$77 || $1!=$79 || $1!=$81 || $1!=$83 || $1!=$85 || $1!=$87 || $1!=$89 || $1!=$91 || $1!=$93 || $1!=$95 || $1!=$97 || $1!=$99 || $1!=$101 || $1!=$103 || $1!=$105 || $1!=$107 || $1!=$109 || $1!=$111 || $1!=$113 || $1!=$115 || $1!=$117 || $1!=$119 || $1!=$121 || $1!=$123 || $1!=$125 || $1!=$127 || $1!=$129 || $1!=$131 || $1!=$133 || $1!=$135 || $1!=$137 || $1!=$139 || $1!=$141 || $1!=$143 || $1!=$145 || $1!=$147 || $1!=$149 || $1!=$151 || $1!=$153 || $1!=$155 || $1!=$157 || $1!=$159 || $1!=$161 || $1!=$163 || $1!=$165 || $1!=$167 || $1!=$169 || $1!=$171 || $1!=$173 || $1!=$175 || $1!=$177 || $1!=$179 || $1!=$181 || $1!=$183 || $1!=$185 || $1!=$187 || $1!=$189 || $1!=$191 || $1!=$193 || $1!=$195 || $1!=$197 || $1!=$199 || $1!=$201 || $1!=$203 || $1!=$205 || $1!=$207 || $1!=$209 || $1!=$211 || $1!=$213 || $1!=$215 ){print $1,$2}}' |head
# 染色体数量和顺序一致，则合并数据列
paste `ls -1 *chr_mapped_uniq_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"; print "chrom","H9_d0_BPK25_input_rep1","H9_d0_BPK25_input_rep2","H9_d0_BPK25_input_rep3","H9_d0_BPK25_IP_CTCF_rep1","H9_d0_BPK25_IP_CTCF_rep2","H9_d0_BPK25_IP_CTCF_rep3","H9_d0_BPK25_IP_MBD3_rep1","H9_d0_BPK25_IP_MBD3_rep2","H9_d0_BPK25_IP_MBD3_rep3","H9_d0_BPK25_IP_TET1_rep1","H9_d0_BPK25_IP_TET1_rep2","H9_d0_BPK25_IP_TET1_rep3","H9_d0_BPK25_IP_TET2_rep1","H9_d0_BPK25_IP_TET2_rep2","H9_d0_BPK25_IP_TET2_rep3","H9_d0_BPK25_IP_TET3_rep1","H9_d0_BPK25_IP_TET3_rep2","H9_d0_BPK25_IP_TET3_rep3","H9_d0_DMSO_input_rep1","H9_d0_DMSO_input_rep2","H9_d0_DMSO_input_rep3","H9_d0_DMSO_IP_CTCF_rep1","H9_d0_DMSO_IP_CTCF_rep2","H9_d0_DMSO_IP_CTCF_rep3","H9_d0_DMSO_IP_MBD3_rep1","H9_d0_DMSO_IP_MBD3_rep2","H9_d0_DMSO_IP_MBD3_rep3","H9_d0_DMSO_IP_TET1_rep1","H9_d0_DMSO_IP_TET1_rep2","H9_d0_DMSO_IP_TET1_rep3","H9_d0_DMSO_IP_TET2_rep1","H9_d0_DMSO_IP_TET2_rep2","H9_d0_DMSO_IP_TET2_rep3","H9_d0_DMSO_IP_TET3_rep1","H9_d0_DMSO_IP_TET3_rep2","H9_d0_DMSO_IP_TET3_rep3","H9_d4_BPK25_input_rep1","H9_d4_BPK25_input_rep2","H9_d4_BPK25_input_rep3","H9_d4_BPK25_IP_CTCF_rep1","H9_d4_BPK25_IP_CTCF_rep2","H9_d4_BPK25_IP_CTCF_rep3","H9_d4_BPK25_IP_MBD3_rep1","H9_d4_BPK25_IP_MBD3_rep2","H9_d4_BPK25_IP_MBD3_rep3","H9_d4_BPK25_IP_TET1_rep1","H9_d4_BPK25_IP_TET1_rep2","H9_d4_BPK25_IP_TET1_rep3","H9_d4_BPK25_IP_TET2_rep1","H9_d4_BPK25_IP_TET2_rep2","H9_d4_BPK25_IP_TET2_rep3","H9_d4_BPK25_IP_TET3_rep1","H9_d4_BPK25_IP_TET3_rep2","H9_d4_BPK25_IP_TET3_rep3","H9_d4_DMSO_input_rep1","H9_d4_DMSO_input_rep2","H9_d4_DMSO_input_rep3","H9_d4_DMSO_IP_CTCF_rep1","H9_d4_DMSO_IP_CTCF_rep2","H9_d4_DMSO_IP_CTCF_rep3","H9_d4_DMSO_IP_MBD3_rep1","H9_d4_DMSO_IP_MBD3_rep2","H9_d4_DMSO_IP_MBD3_rep3","H9_d4_DMSO_IP_TET1_rep1","H9_d4_DMSO_IP_TET1_rep2","H9_d4_DMSO_IP_TET1_rep3","H9_d4_DMSO_IP_TET2_rep1","H9_d4_DMSO_IP_TET2_rep2","H9_d4_DMSO_IP_TET2_rep3","H9_d4_DMSO_IP_TET3_rep1","H9_d4_DMSO_IP_TET3_rep2","H9_d4_DMSO_IP_TET3_rep3","H9_d9_BPK25_input_rep1","H9_d9_BPK25_input_rep2","H9_d9_BPK25_input_rep3","H9_d9_BPK25_IP_CTCF_rep1","H9_d9_BPK25_IP_CTCF_rep2","H9_d9_BPK25_IP_CTCF_rep3","H9_d9_BPK25_IP_MBD3_rep1","H9_d9_BPK25_IP_MBD3_rep2","H9_d9_BPK25_IP_MBD3_rep3","H9_d9_BPK25_IP_TET1_rep1","H9_d9_BPK25_IP_TET1_rep2","H9_d9_BPK25_IP_TET1_rep3","H9_d9_BPK25_IP_TET2_rep1","H9_d9_BPK25_IP_TET2_rep2","H9_d9_BPK25_IP_TET2_rep3","H9_d9_BPK25_IP_TET3_rep1","H9_d9_BPK25_IP_TET3_rep2","H9_d9_BPK25_IP_TET3_rep3","H9_d9_DMSO_input_rep1","H9_d9_DMSO_input_rep2","H9_d9_DMSO_input_rep3","H9_d9_DMSO_IP_CTCF_rep1","H9_d9_DMSO_IP_CTCF_rep2","H9_d9_DMSO_IP_CTCF_rep3","H9_d9_DMSO_IP_MBD3_rep1","H9_d9_DMSO_IP_MBD3_rep2","H9_d9_DMSO_IP_MBD3_rep3","H9_d9_DMSO_IP_TET1_rep1","H9_d9_DMSO_IP_TET1_rep2","H9_d9_DMSO_IP_TET1_rep3","H9_d9_DMSO_IP_TET2_rep1","H9_d9_DMSO_IP_TET2_rep2","H9_d9_DMSO_IP_TET2_rep3","H9_d9_DMSO_IP_TET3_rep1","H9_d9_DMSO_IP_TET3_rep2","H9_d9_DMSO_IP_TET3_rep3"}{if($1!=$3 && $1!=$5 && $1!=$7 && $1!=$9 && $1!=$11 && $1!=$13 && $1!=$15 && $1!=$17 && $1!=$19 && $1!=$21 && $1!=$23 && $1!=$25 && $1!=$27 && $1!=$29 && $1!=$31 && $1!=$33 && $1!=$35 && $1!=$37 && $1!=$39 && $1!=$41 && $1!=$43 && $1!=$45 && $1!=$47 && $1!=$49 && $1!=$51 && $1!=$53 && $1!=$55 && $1!=$57 && $1!=$59 && $1!=$61 && $1!=$63 && $1!=$65 && $1!=$67 && $1!=$69 && $1!=$71 && $1!=$73 && $1!=$75 && $1!=$77 && $1!=$79 && $1!=$81 && $1!=$83 && $1!=$85 && $1!=$87 && $1!=$89 && $1!=$91 && $1!=$93 && $1!=$95 && $1!=$97 && $1!=$99 && $1!=$101 && $1!=$103 && $1!=$105 && $1!=$107 && $1!=$109 && $1!=$111 && $1!=$113 && $1!=$115 && $1!=$117 && $1!=$119 && $1!=$121 && $1!=$123 && $1!=$125 && $1!=$127 && $1!=$129 && $1!=$131 && $1!=$133 && $1!=$135 && $1!=$137 && $1!=$139 && $1!=$141 && $1!=$143 && $1!=$145 && $1!=$147 && $1!=$149 && $1!=$151 && $1!=$153 && $1!=$155 && $1!=$157 && $1!=$159 && $1!=$161 && $1!=$163 && $1!=$165 && $1!=$167 && $1!=$169 && $1!=$171 && $1!=$173 && $1!=$175 && $1!=$177 && $1!=$179 && $1!=$181 && $1!=$183 && $1!=$185 && $1!=$187 && $1!=$189 && $1!=$191 && $1!=$193 && $1!=$195 && $1!=$197 && $1!=$199 && $1!=$201 && $1!=$203 && $1!=$205 && $1!=$207 && $1!=$209 && $1!=$211 && $1!=$213 && $1!=$215){print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36,$38,$40,$42,$44,$46,$48,$50,$52,$54,$56,$58,$60,$62,$64,$66,$68,$70,$72,$74,$76,$78,$80,$82,$84,$86,$88,$90,$92,$94,$96,$98,$100,$102,$104,$106,$108,$110,$112,$114,$116,$118,$120,$122,$124,$126,$128,$130,$132,$134,$136,$138,$140,$142,$144,$146,$148,$150,$152,$154,$156,$158,$160,$162,$164,$166,$168,$170,$172,$174,$176,$178,$180,$182,$184,$186,$188,$190,$192,$194,$196,$198,$200,$202,$204,$206,$208,$210,$212,$214,$216}}' |head
paste `ls -1 *chr_mapped_uniq_reads_stat.txt|sort` |awk 'BEGIN{FS="\t";OFS="\t"; print "chrom","H9_d0_BPK25_input_rep1","H9_d0_BPK25_input_rep2","H9_d0_BPK25_input_rep3","H9_d0_BPK25_IP_CTCF_rep1","H9_d0_BPK25_IP_CTCF_rep2","H9_d0_BPK25_IP_CTCF_rep3","H9_d0_BPK25_IP_MBD3_rep1","H9_d0_BPK25_IP_MBD3_rep2","H9_d0_BPK25_IP_MBD3_rep3","H9_d0_BPK25_IP_TET1_rep1","H9_d0_BPK25_IP_TET1_rep2","H9_d0_BPK25_IP_TET1_rep3","H9_d0_BPK25_IP_TET2_rep1","H9_d0_BPK25_IP_TET2_rep2","H9_d0_BPK25_IP_TET2_rep3","H9_d0_BPK25_IP_TET3_rep1","H9_d0_BPK25_IP_TET3_rep2","H9_d0_BPK25_IP_TET3_rep3","H9_d0_DMSO_input_rep1","H9_d0_DMSO_input_rep2","H9_d0_DMSO_input_rep3","H9_d0_DMSO_IP_CTCF_rep1","H9_d0_DMSO_IP_CTCF_rep2","H9_d0_DMSO_IP_CTCF_rep3","H9_d0_DMSO_IP_MBD3_rep1","H9_d0_DMSO_IP_MBD3_rep2","H9_d0_DMSO_IP_MBD3_rep3","H9_d0_DMSO_IP_TET1_rep1","H9_d0_DMSO_IP_TET1_rep2","H9_d0_DMSO_IP_TET1_rep3","H9_d0_DMSO_IP_TET2_rep1","H9_d0_DMSO_IP_TET2_rep2","H9_d0_DMSO_IP_TET2_rep3","H9_d0_DMSO_IP_TET3_rep1","H9_d0_DMSO_IP_TET3_rep2","H9_d0_DMSO_IP_TET3_rep3","H9_d4_BPK25_input_rep1","H9_d4_BPK25_input_rep2","H9_d4_BPK25_input_rep3","H9_d4_BPK25_IP_CTCF_rep1","H9_d4_BPK25_IP_CTCF_rep2","H9_d4_BPK25_IP_CTCF_rep3","H9_d4_BPK25_IP_MBD3_rep1","H9_d4_BPK25_IP_MBD3_rep2","H9_d4_BPK25_IP_MBD3_rep3","H9_d4_BPK25_IP_TET1_rep1","H9_d4_BPK25_IP_TET1_rep2","H9_d4_BPK25_IP_TET1_rep3","H9_d4_BPK25_IP_TET2_rep1","H9_d4_BPK25_IP_TET2_rep2","H9_d4_BPK25_IP_TET2_rep3","H9_d4_BPK25_IP_TET3_rep1","H9_d4_BPK25_IP_TET3_rep2","H9_d4_BPK25_IP_TET3_rep3","H9_d4_DMSO_input_rep1","H9_d4_DMSO_input_rep2","H9_d4_DMSO_input_rep3","H9_d4_DMSO_IP_CTCF_rep1","H9_d4_DMSO_IP_CTCF_rep2","H9_d4_DMSO_IP_CTCF_rep3","H9_d4_DMSO_IP_MBD3_rep1","H9_d4_DMSO_IP_MBD3_rep2","H9_d4_DMSO_IP_MBD3_rep3","H9_d4_DMSO_IP_TET1_rep1","H9_d4_DMSO_IP_TET1_rep2","H9_d4_DMSO_IP_TET1_rep3","H9_d4_DMSO_IP_TET2_rep1","H9_d4_DMSO_IP_TET2_rep2","H9_d4_DMSO_IP_TET2_rep3","H9_d4_DMSO_IP_TET3_rep1","H9_d4_DMSO_IP_TET3_rep2","H9_d4_DMSO_IP_TET3_rep3","H9_d9_BPK25_input_rep1","H9_d9_BPK25_input_rep2","H9_d9_BPK25_input_rep3","H9_d9_BPK25_IP_CTCF_rep1","H9_d9_BPK25_IP_CTCF_rep2","H9_d9_BPK25_IP_CTCF_rep3","H9_d9_BPK25_IP_MBD3_rep1","H9_d9_BPK25_IP_MBD3_rep2","H9_d9_BPK25_IP_MBD3_rep3","H9_d9_BPK25_IP_TET1_rep1","H9_d9_BPK25_IP_TET1_rep2","H9_d9_BPK25_IP_TET1_rep3","H9_d9_BPK25_IP_TET2_rep1","H9_d9_BPK25_IP_TET2_rep2","H9_d9_BPK25_IP_TET2_rep3","H9_d9_BPK25_IP_TET3_rep1","H9_d9_BPK25_IP_TET3_rep2","H9_d9_BPK25_IP_TET3_rep3","H9_d9_DMSO_input_rep1","H9_d9_DMSO_input_rep2","H9_d9_DMSO_input_rep3","H9_d9_DMSO_IP_CTCF_rep1","H9_d9_DMSO_IP_CTCF_rep2","H9_d9_DMSO_IP_CTCF_rep3","H9_d9_DMSO_IP_MBD3_rep1","H9_d9_DMSO_IP_MBD3_rep2","H9_d9_DMSO_IP_MBD3_rep3","H9_d9_DMSO_IP_TET1_rep1","H9_d9_DMSO_IP_TET1_rep2","H9_d9_DMSO_IP_TET1_rep3","H9_d9_DMSO_IP_TET2_rep1","H9_d9_DMSO_IP_TET2_rep2","H9_d9_DMSO_IP_TET2_rep3","H9_d9_DMSO_IP_TET3_rep1","H9_d9_DMSO_IP_TET3_rep2","H9_d9_DMSO_IP_TET3_rep3"}{if($1!=$3 && $1!=$5 && $1!=$7 && $1!=$9 && $1!=$11 && $1!=$13 && $1!=$15 && $1!=$17 && $1!=$19 && $1!=$21 && $1!=$23 && $1!=$25 && $1!=$27 && $1!=$29 && $1!=$31 && $1!=$33 && $1!=$35 && $1!=$37 && $1!=$39 && $1!=$41 && $1!=$43 && $1!=$45 && $1!=$47 && $1!=$49 && $1!=$51 && $1!=$53 && $1!=$55 && $1!=$57 && $1!=$59 && $1!=$61 && $1!=$63 && $1!=$65 && $1!=$67 && $1!=$69 && $1!=$71 && $1!=$73 && $1!=$75 && $1!=$77 && $1!=$79 && $1!=$81 && $1!=$83 && $1!=$85 && $1!=$87 && $1!=$89 && $1!=$91 && $1!=$93 && $1!=$95 && $1!=$97 && $1!=$99 && $1!=$101 && $1!=$103 && $1!=$105 && $1!=$107 && $1!=$109 && $1!=$111 && $1!=$113 && $1!=$115 && $1!=$117 && $1!=$119 && $1!=$121 && $1!=$123 && $1!=$125 && $1!=$127 && $1!=$129 && $1!=$131 && $1!=$133 && $1!=$135 && $1!=$137 && $1!=$139 && $1!=$141 && $1!=$143 && $1!=$145 && $1!=$147 && $1!=$149 && $1!=$151 && $1!=$153 && $1!=$155 && $1!=$157 && $1!=$159 && $1!=$161 && $1!=$163 && $1!=$165 && $1!=$167 && $1!=$169 && $1!=$171 && $1!=$173 && $1!=$175 && $1!=$177 && $1!=$179 && $1!=$181 && $1!=$183 && $1!=$185 && $1!=$187 && $1!=$189 && $1!=$191 && $1!=$193 && $1!=$195 && $1!=$197 && $1!=$199 && $1!=$201 && $1!=$203 && $1!=$205 && $1!=$207 && $1!=$209 && $1!=$211 && $1!=$213 && $1!=$215){print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32,$34,$36,$38,$40,$42,$44,$46,$48,$50,$52,$54,$56,$58,$60,$62,$64,$66,$68,$70,$72,$74,$76,$78,$80,$82,$84,$86,$88,$90,$92,$94,$96,$98,$100,$102,$104,$106,$108,$110,$112,$114,$116,$118,$120,$122,$124,$126,$128,$130,$132,$134,$136,$138,$140,$142,$144,$146,$148,$150,$152,$154,$156,$158,$160,$162,$164,$166,$168,$170,$172,$174,$176,$178,$180,$182,$184,$186,$188,$190,$192,$194,$196,$198,$200,$202,$204,$206,$208,$210,$212,$214,$216}}' >all_samples_chr_mapped_uniq_reads_stat.summary.txt

```

```bash
# 绘制reads 分布热图时要用到 rep123合并的bam文件
# 创建 bam merge 脚本 5_make_all_H9_d049_BPK25_ChIP_Seq_rep123_merged_bam_and_index.sh 内容如下：
# # 合并各重复的 bam
# samtools merge --threads 72 hESC_H9_d0_BPK25_input_rep123.merged.sorted.bam hESC_H9_d0_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_input_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d0_BPK25_IP_CTCF_rep123.merged.sorted.bam hESC_H9_d0_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d0_BPK25_IP_MBD3_rep123.merged.sorted.bam hESC_H9_d0_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d0_BPK25_IP_TET1_rep123.merged.sorted.bam hESC_H9_d0_BPK25_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET1_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d0_BPK25_IP_TET2_rep123.merged.sorted.bam hESC_H9_d0_BPK25_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET2_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d0_BPK25_IP_TET3_rep123.merged.sorted.bam hESC_H9_d0_BPK25_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET3_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d0_DMSO_input_rep123.merged.sorted.bam hESC_H9_d0_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_input_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d0_DMSO_IP_CTCF_rep123.merged.sorted.bam hESC_H9_d0_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d0_DMSO_IP_MBD3_rep123.merged.sorted.bam hESC_H9_d0_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d0_DMSO_IP_TET1_rep123.merged.sorted.bam hESC_H9_d0_DMSO_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET1_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d0_DMSO_IP_TET2_rep123.merged.sorted.bam hESC_H9_d0_DMSO_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET2_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d0_DMSO_IP_TET3_rep123.merged.sorted.bam hESC_H9_d0_DMSO_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET3_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d4_BPK25_input_rep123.merged.sorted.bam hESC_H9_d4_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_input_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d4_BPK25_IP_CTCF_rep123.merged.sorted.bam hESC_H9_d4_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d4_BPK25_IP_MBD3_rep123.merged.sorted.bam hESC_H9_d4_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d4_BPK25_IP_TET1_rep123.merged.sorted.bam hESC_H9_d4_BPK25_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET1_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d4_BPK25_IP_TET2_rep123.merged.sorted.bam hESC_H9_d4_BPK25_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET2_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d4_BPK25_IP_TET3_rep123.merged.sorted.bam hESC_H9_d4_BPK25_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET3_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d4_DMSO_input_rep123.merged.sorted.bam hESC_H9_d4_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_input_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d4_DMSO_IP_CTCF_rep123.merged.sorted.bam hESC_H9_d4_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d4_DMSO_IP_MBD3_rep123.merged.sorted.bam hESC_H9_d4_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d4_DMSO_IP_TET1_rep123.merged.sorted.bam hESC_H9_d4_DMSO_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET1_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d4_DMSO_IP_TET2_rep123.merged.sorted.bam hESC_H9_d4_DMSO_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET2_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d4_DMSO_IP_TET3_rep123.merged.sorted.bam hESC_H9_d4_DMSO_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET3_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d9_BPK25_input_rep123.merged.sorted.bam hESC_H9_d9_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_input_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d9_BPK25_IP_CTCF_rep123.merged.sorted.bam hESC_H9_d9_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d9_BPK25_IP_MBD3_rep123.merged.sorted.bam hESC_H9_d9_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d9_BPK25_IP_TET1_rep123.merged.sorted.bam hESC_H9_d9_BPK25_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET1_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d9_BPK25_IP_TET2_rep123.merged.sorted.bam hESC_H9_d9_BPK25_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET2_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d9_BPK25_IP_TET3_rep123.merged.sorted.bam hESC_H9_d9_BPK25_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET3_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d9_DMSO_input_rep123.merged.sorted.bam hESC_H9_d9_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_input_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d9_DMSO_IP_CTCF_rep123.merged.sorted.bam hESC_H9_d9_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d9_DMSO_IP_MBD3_rep123.merged.sorted.bam hESC_H9_d9_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d9_DMSO_IP_TET1_rep123.merged.sorted.bam hESC_H9_d9_DMSO_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET1_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d9_DMSO_IP_TET2_rep123.merged.sorted.bam hESC_H9_d9_DMSO_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET2_rep3.bowtie2.sorted.bam
# samtools merge --threads 72 hESC_H9_d9_DMSO_IP_TET3_rep123.merged.sorted.bam hESC_H9_d9_DMSO_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET3_rep3.bowtie2.sorted.bam
# # 为所有 merged bam 建立索引
# for i in `ls -1 *_rep123.merged.sorted.bam`; do echo $i; samtools index -@ 72 $i;done

nohup bash 5_make_all_H9_d049_BPK25_ChIP_Seq_rep123_merged_bam_and_index.sh >5_make_all_H9_d049_BPK25_ChIP_Seq_rep123_merged_bam_and_index.sh.runlog 2>&1 &
tail -f 5_make_all_H9_d049_BPK25_ChIP_Seq_rep123_merged_bam_and_index.sh.runlog

```

# Detection of enrichment and repeatability

## heatmap

### d049 all

```bash
# 由于文件多，命令长，构造 5_make_cor_heatmap_and_PCA_plot_for_d049_all_samples.sh 内容如下：
# # 使用 deeptools 分析全部样本相似度
# ## 先对全基因组分 bin (default 10kb) 统计 uniq reads数量
# multiBamSummary bins --bamfiles hESC_H9_d0_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_input_rep3.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET1_rep3.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET2_rep3.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET3_rep3.bowtie2.sorted.bam hESC_H9_d0_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_input_rep3.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET1_rep3.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET2_rep3.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET3_rep3.bowtie2.sorted.bam hESC_H9_d4_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_input_rep3.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET1_rep3.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET2_rep3.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET3_rep3.bowtie2.sorted.bam hESC_H9_d4_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_input_rep3.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET1_rep3.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET2_rep3.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET3_rep3.bowtie2.sorted.bam hESC_H9_d9_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_input_rep3.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET1_rep3.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET2_rep3.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET3_rep3.bowtie2.sorted.bam hESC_H9_d9_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_input_rep3.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET1_rep3.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET2_rep3.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET3_rep3.bowtie2.sorted.bam --labels H9_d0_DMSO_input_rep1 H9_d0_DMSO_input_rep2 H9_d0_DMSO_input_rep3 H9_d0_DMSO_IP_CTCF_rep1 H9_d0_DMSO_IP_CTCF_rep2 H9_d0_DMSO_IP_CTCF_rep3 H9_d0_DMSO_IP_MBD3_rep1 H9_d0_DMSO_IP_MBD3_rep2 H9_d0_DMSO_IP_MBD3_rep3 H9_d0_DMSO_IP_TET1_rep1 H9_d0_DMSO_IP_TET1_rep2 H9_d0_DMSO_IP_TET1_rep3 H9_d0_DMSO_IP_TET2_rep1 H9_d0_DMSO_IP_TET2_rep2 H9_d0_DMSO_IP_TET2_rep3 H9_d0_DMSO_IP_TET3_rep1 H9_d0_DMSO_IP_TET3_rep2 H9_d0_DMSO_IP_TET3_rep3 H9_d0_BPK25_input_rep1 H9_d0_BPK25_input_rep2 H9_d0_BPK25_input_rep3 H9_d0_BPK25_IP_CTCF_rep1 H9_d0_BPK25_IP_CTCF_rep2 H9_d0_BPK25_IP_CTCF_rep3 H9_d0_BPK25_IP_MBD3_rep1 H9_d0_BPK25_IP_MBD3_rep2 H9_d0_BPK25_IP_MBD3_rep3 H9_d0_BPK25_IP_TET1_rep1 H9_d0_BPK25_IP_TET1_rep2 H9_d0_BPK25_IP_TET1_rep3 H9_d0_BPK25_IP_TET2_rep1 H9_d0_BPK25_IP_TET2_rep2 H9_d0_BPK25_IP_TET2_rep3 H9_d0_BPK25_IP_TET3_rep1 H9_d0_BPK25_IP_TET3_rep2 H9_d0_BPK25_IP_TET3_rep3 H9_d4_DMSO_input_rep1 H9_d4_DMSO_input_rep2 H9_d4_DMSO_input_rep3 H9_d4_DMSO_IP_CTCF_rep1 H9_d4_DMSO_IP_CTCF_rep2 H9_d4_DMSO_IP_CTCF_rep3 H9_d4_DMSO_IP_MBD3_rep1 H9_d4_DMSO_IP_MBD3_rep2 H9_d4_DMSO_IP_MBD3_rep3 H9_d4_DMSO_IP_TET1_rep1 H9_d4_DMSO_IP_TET1_rep2 H9_d4_DMSO_IP_TET1_rep3 H9_d4_DMSO_IP_TET2_rep1 H9_d4_DMSO_IP_TET2_rep2 H9_d4_DMSO_IP_TET2_rep3 H9_d4_DMSO_IP_TET3_rep1 H9_d4_DMSO_IP_TET3_rep2 H9_d4_DMSO_IP_TET3_rep3 H9_d4_BPK25_input_rep1 H9_d4_BPK25_input_rep2 H9_d4_BPK25_input_rep3 H9_d4_BPK25_IP_CTCF_rep1 H9_d4_BPK25_IP_CTCF_rep2 H9_d4_BPK25_IP_CTCF_rep3 H9_d4_BPK25_IP_MBD3_rep1 H9_d4_BPK25_IP_MBD3_rep2 H9_d4_BPK25_IP_MBD3_rep3 H9_d4_BPK25_IP_TET1_rep1 H9_d4_BPK25_IP_TET1_rep2 H9_d4_BPK25_IP_TET1_rep3 H9_d4_BPK25_IP_TET2_rep1 H9_d4_BPK25_IP_TET2_rep2 H9_d4_BPK25_IP_TET2_rep3 H9_d4_BPK25_IP_TET3_rep1 H9_d4_BPK25_IP_TET3_rep2 H9_d4_BPK25_IP_TET3_rep3 H9_d9_DMSO_input_rep1 H9_d9_DMSO_input_rep2 H9_d9_DMSO_input_rep3 H9_d9_DMSO_IP_CTCF_rep1 H9_d9_DMSO_IP_CTCF_rep2 H9_d9_DMSO_IP_CTCF_rep3 H9_d9_DMSO_IP_MBD3_rep1 H9_d9_DMSO_IP_MBD3_rep2 H9_d9_DMSO_IP_MBD3_rep3 H9_d9_DMSO_IP_TET1_rep1 H9_d9_DMSO_IP_TET1_rep2 H9_d9_DMSO_IP_TET1_rep3 H9_d9_DMSO_IP_TET2_rep1 H9_d9_DMSO_IP_TET2_rep2 H9_d9_DMSO_IP_TET2_rep3 H9_d9_DMSO_IP_TET3_rep1 H9_d9_DMSO_IP_TET3_rep2 H9_d9_DMSO_IP_TET3_rep3 H9_d9_BPK25_input_rep1 H9_d9_BPK25_input_rep2 H9_d9_BPK25_input_rep3 H9_d9_BPK25_IP_CTCF_rep1 H9_d9_BPK25_IP_CTCF_rep2 H9_d9_BPK25_IP_CTCF_rep3 H9_d9_BPK25_IP_MBD3_rep1 H9_d9_BPK25_IP_MBD3_rep2 H9_d9_BPK25_IP_MBD3_rep3 H9_d9_BPK25_IP_TET1_rep1 H9_d9_BPK25_IP_TET1_rep2 H9_d9_BPK25_IP_TET1_rep3 H9_d9_BPK25_IP_TET2_rep1 H9_d9_BPK25_IP_TET2_rep2 H9_d9_BPK25_IP_TET2_rep3 H9_d9_BPK25_IP_TET3_rep1 H9_d9_BPK25_IP_TET3_rep2 H9_d9_BPK25_IP_TET3_rep3 -out multiBamSummary_results_for_H9_d049_DMSO_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# ## 用 plotCorrelation 和 plotPCA 分析样本的相似度
# plotCorrelation --corData multiBamSummary_results_for_H9_d049_DMSO_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_H9_d049_DMSO_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_H9_d049_DMSO_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.txt --colorMap coolwarm --plotNumbers
# plotPCA --corData multiBamSummary_results_for_H9_d049_DMSO_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --plotFile plotPCA_plot_for_H9_d049_DMSO_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_H9_d049_DMSO_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.txt
# plotPCA --rowCenter --corData multiBamSummary_results_for_H9_d049_DMSO_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --plotFile plotPCA_plot_for_H9_d049_DMSO_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_H9_d049_DMSO_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.rowCenter.txt

# 上传并运行
nohup bash 5_make_cor_heatmap_and_PCA_plot_for_d049_all_samples.sh >5_make_cor_heatmap_and_PCA_plot_for_d049_all_samples.sh.runlog 2>&1 &
tail -f 5_make_cor_heatmap_and_PCA_plot_for_d049_all_samples.sh.runlog

```

### d0_DMSO_all

```bash
# 由于文件多，命令长，构造 5_make_cor_heatmap_and_PCA_plot_for_d0_DMSO_all_samples.sh 内容如下：
# # 使用 deeptools 分析全部样本相似度
# ## 先对全基因组分 bin (default 10kb) 统计 uniq reads数量
# multiBamSummary bins --bamfiles hESC_H9_d0_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_input_rep3.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET1_rep3.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET2_rep3.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET3_rep3.bowtie2.sorted.bam --labels H9_d0_DMSO_input_rep1 H9_d0_DMSO_input_rep2 H9_d0_DMSO_input_rep3 H9_d0_DMSO_IP_CTCF_rep1 H9_d0_DMSO_IP_CTCF_rep2 H9_d0_DMSO_IP_CTCF_rep3 H9_d0_DMSO_IP_MBD3_rep1 H9_d0_DMSO_IP_MBD3_rep2 H9_d0_DMSO_IP_MBD3_rep3 H9_d0_DMSO_IP_TET1_rep1 H9_d0_DMSO_IP_TET1_rep2 H9_d0_DMSO_IP_TET1_rep3 H9_d0_DMSO_IP_TET2_rep1 H9_d0_DMSO_IP_TET2_rep2 H9_d0_DMSO_IP_TET2_rep3 H9_d0_DMSO_IP_TET3_rep1 H9_d0_DMSO_IP_TET3_rep2 H9_d0_DMSO_IP_TET3_rep3 -out multiBamSummary_results_for_H9_d0_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# ## 用 plotCorrelation 和 plotPCA 分析样本的相似度
# plotCorrelation --corData multiBamSummary_results_for_H9_d0_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_H9_d0_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_H9_d0_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.txt --colorMap coolwarm --plotNumbers
# plotPCA --corData multiBamSummary_results_for_H9_d0_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --plotFile plotPCA_plot_for_H9_d0_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_H9_d0_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.txt
# plotPCA --rowCenter --corData multiBamSummary_results_for_H9_d0_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --plotFile plotPCA_plot_for_H9_d0_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_H9_d0_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.rowCenter.txt

# 上传并运行
nohup bash 5_make_cor_heatmap_and_PCA_plot_for_d0_DMSO_all_samples.sh >5_make_cor_heatmap_and_PCA_plot_for_d0_DMSO_all_samples.sh.runlog 2>&1 &
tail -f 5_make_cor_heatmap_and_PCA_plot_for_d0_DMSO_all_samples.sh.runlog

```

### d0_BPK25_all

```bash
# 由于文件多，命令长，构造 5_make_cor_heatmap_and_PCA_plot_for_d0_BPK25_all_samples.sh 内容如下：
# # 使用 deeptools 分析全部样本相似度
# ## 先对全基因组分 bin (default 10kb) 统计 uniq reads数量
# multiBamSummary bins --bamfiles hESC_H9_d0_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_input_rep3.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET1_rep3.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET2_rep3.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET3_rep3.bowtie2.sorted.bam --labels H9_d0_BPK25_input_rep1 H9_d0_BPK25_input_rep2 H9_d0_BPK25_input_rep3 H9_d0_BPK25_IP_CTCF_rep1 H9_d0_BPK25_IP_CTCF_rep2 H9_d0_BPK25_IP_CTCF_rep3 H9_d0_BPK25_IP_MBD3_rep1 H9_d0_BPK25_IP_MBD3_rep2 H9_d0_BPK25_IP_MBD3_rep3 H9_d0_BPK25_IP_TET1_rep1 H9_d0_BPK25_IP_TET1_rep2 H9_d0_BPK25_IP_TET1_rep3 H9_d0_BPK25_IP_TET2_rep1 H9_d0_BPK25_IP_TET2_rep2 H9_d0_BPK25_IP_TET2_rep3 H9_d0_BPK25_IP_TET3_rep1 H9_d0_BPK25_IP_TET3_rep2 H9_d0_BPK25_IP_TET3_rep3 -out multiBamSummary_results_for_H9_d0_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# ## 用 plotCorrelation 和 plotPCA 分析样本的相似度
# plotCorrelation --corData multiBamSummary_results_for_H9_d0_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_H9_d0_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_H9_d0_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.txt --colorMap coolwarm --plotNumbers
# plotPCA --corData multiBamSummary_results_for_H9_d0_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --plotFile plotPCA_plot_for_H9_d0_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_H9_d0_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.txt
# plotPCA --rowCenter --corData multiBamSummary_results_for_H9_d0_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --plotFile plotPCA_plot_for_H9_d0_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_H9_d0_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.rowCenter.txt

# 上传并运行
nohup bash 5_make_cor_heatmap_and_PCA_plot_for_d0_BPK25_all_samples.sh >5_make_cor_heatmap_and_PCA_plot_for_d0_BPK25_all_samples.sh.runlog 2>&1 &
tail -f 5_make_cor_heatmap_and_PCA_plot_for_d0_BPK25_all_samples.sh.runlog

```

### d4_DMSO_all

```bash
# 由于文件多，命令长，构造 5_make_cor_heatmap_and_PCA_plot_for_d4_DMSO_all_samples.sh 内容如下：
# # 使用 deeptools 分析全部样本相似度
# ## 先对全基因组分 bin (default 10kb) 统计 uniq reads数量
# multiBamSummary bins --bamfiles hESC_H9_d4_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_input_rep3.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET1_rep3.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET2_rep3.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET3_rep3.bowtie2.sorted.bam --labels H9_d4_DMSO_input_rep1 H9_d4_DMSO_input_rep2 H9_d4_DMSO_input_rep3 H9_d4_DMSO_IP_CTCF_rep1 H9_d4_DMSO_IP_CTCF_rep2 H9_d4_DMSO_IP_CTCF_rep3 H9_d4_DMSO_IP_MBD3_rep1 H9_d4_DMSO_IP_MBD3_rep2 H9_d4_DMSO_IP_MBD3_rep3 H9_d4_DMSO_IP_TET1_rep1 H9_d4_DMSO_IP_TET1_rep2 H9_d4_DMSO_IP_TET1_rep3 H9_d4_DMSO_IP_TET2_rep1 H9_d4_DMSO_IP_TET2_rep2 H9_d4_DMSO_IP_TET2_rep3 H9_d4_DMSO_IP_TET3_rep1 H9_d4_DMSO_IP_TET3_rep2 H9_d4_DMSO_IP_TET3_rep3 -out multiBamSummary_results_for_H9_d4_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# ## 用 plotCorrelation 和 plotPCA 分析样本的相似度
# plotCorrelation --corData multiBamSummary_results_for_H9_d4_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_H9_d4_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_H9_d4_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.txt --colorMap coolwarm --plotNumbers
# plotPCA --corData multiBamSummary_results_for_H9_d4_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --plotFile plotPCA_plot_for_H9_d4_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_H9_d4_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.txt
# plotPCA --rowCenter --corData multiBamSummary_results_for_H9_d4_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --plotFile plotPCA_plot_for_H9_d4_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_H9_d4_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.rowCenter.txt

# 上传并运行
nohup bash 5_make_cor_heatmap_and_PCA_plot_for_d4_DMSO_all_samples.sh >5_make_cor_heatmap_and_PCA_plot_for_d4_DMSO_all_samples.sh.runlog 2>&1 &
tail -f 5_make_cor_heatmap_and_PCA_plot_for_d4_DMSO_all_samples.sh.runlog

```

### d4_BPK25_all

```bash
# 由于文件多，命令长，构造 5_make_cor_heatmap_and_PCA_plot_for_d4_BPK25_all_samples.sh 内容如下：
# # 使用 deeptools 分析全部样本相似度
# ## 先对全基因组分 bin (default 10kb) 统计 uniq reads数量
# multiBamSummary bins --bamfiles hESC_H9_d4_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_input_rep3.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET1_rep3.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET2_rep3.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET3_rep3.bowtie2.sorted.bam --labels H9_d4_BPK25_input_rep1 H9_d4_BPK25_input_rep2 H9_d4_BPK25_input_rep3 H9_d4_BPK25_IP_CTCF_rep1 H9_d4_BPK25_IP_CTCF_rep2 H9_d4_BPK25_IP_CTCF_rep3 H9_d4_BPK25_IP_MBD3_rep1 H9_d4_BPK25_IP_MBD3_rep2 H9_d4_BPK25_IP_MBD3_rep3 H9_d4_BPK25_IP_TET1_rep1 H9_d4_BPK25_IP_TET1_rep2 H9_d4_BPK25_IP_TET1_rep3 H9_d4_BPK25_IP_TET2_rep1 H9_d4_BPK25_IP_TET2_rep2 H9_d4_BPK25_IP_TET2_rep3 H9_d4_BPK25_IP_TET3_rep1 H9_d4_BPK25_IP_TET3_rep2 H9_d4_BPK25_IP_TET3_rep3 -out multiBamSummary_results_for_H9_d4_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# ## 用 plotCorrelation 和 plotPCA 分析样本的相似度
# plotCorrelation --corData multiBamSummary_results_for_H9_d4_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_H9_d4_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_H9_d4_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.txt --colorMap coolwarm --plotNumbers
# plotPCA --corData multiBamSummary_results_for_H9_d4_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --plotFile plotPCA_plot_for_H9_d4_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_H9_d4_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.txt
# plotPCA --rowCenter --corData multiBamSummary_results_for_H9_d4_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --plotFile plotPCA_plot_for_H9_d4_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_H9_d4_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.rowCenter.txt

# 上传并运行
nohup bash 5_make_cor_heatmap_and_PCA_plot_for_d4_BPK25_all_samples.sh >5_make_cor_heatmap_and_PCA_plot_for_d4_BPK25_all_samples.sh.runlog 2>&1 &
tail -f 5_make_cor_heatmap_and_PCA_plot_for_d4_BPK25_all_samples.sh.runlog

```

### d9_DMSO_all

```bash
# 由于文件多，命令长，构造 5_make_cor_heatmap_and_PCA_plot_for_d9_DMSO_all_samples.sh 内容如下：
# # 使用 deeptools 分析全部样本相似度
# ## 先对全基因组分 bin (default 10kb) 统计 uniq reads数量
# multiBamSummary bins --bamfiles hESC_H9_d9_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_input_rep3.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET1_rep3.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET2_rep3.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET3_rep3.bowtie2.sorted.bam --labels H9_d9_DMSO_input_rep1 H9_d9_DMSO_input_rep2 H9_d9_DMSO_input_rep3 H9_d9_DMSO_IP_CTCF_rep1 H9_d9_DMSO_IP_CTCF_rep2 H9_d9_DMSO_IP_CTCF_rep3 H9_d9_DMSO_IP_MBD3_rep1 H9_d9_DMSO_IP_MBD3_rep2 H9_d9_DMSO_IP_MBD3_rep3 H9_d9_DMSO_IP_TET1_rep1 H9_d9_DMSO_IP_TET1_rep2 H9_d9_DMSO_IP_TET1_rep3 H9_d9_DMSO_IP_TET2_rep1 H9_d9_DMSO_IP_TET2_rep2 H9_d9_DMSO_IP_TET2_rep3 H9_d9_DMSO_IP_TET3_rep1 H9_d9_DMSO_IP_TET3_rep2 H9_d9_DMSO_IP_TET3_rep3 -out multiBamSummary_results_for_H9_d9_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# ## 用 plotCorrelation 和 plotPCA 分析样本的相似度
# plotCorrelation --corData multiBamSummary_results_for_H9_d9_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_H9_d9_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_H9_d9_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.txt --colorMap coolwarm --plotNumbers
# plotPCA --corData multiBamSummary_results_for_H9_d9_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --plotFile plotPCA_plot_for_H9_d9_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_H9_d9_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.txt
# plotPCA --rowCenter --corData multiBamSummary_results_for_H9_d9_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --plotFile plotPCA_plot_for_H9_d9_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_H9_d9_DMSO_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.rowCenter.txt

# 上传并运行
nohup bash 5_make_cor_heatmap_and_PCA_plot_for_d9_DMSO_all_samples.sh >5_make_cor_heatmap_and_PCA_plot_for_d9_DMSO_all_samples.sh.runlog 2>&1 &
tail -f 5_make_cor_heatmap_and_PCA_plot_for_d9_DMSO_all_samples.sh.runlog

```

### d9_BPK25_all

```bash
# 由于文件多，命令长，构造 5_make_cor_heatmap_and_PCA_plot_for_d9_BPK25_all_samples.sh 内容如下：
# # 使用 deeptools 分析全部样本相似度
# ## 先对全基因组分 bin (default 10kb) 统计 uniq reads数量
# multiBamSummary bins --bamfiles hESC_H9_d9_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_input_rep3.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET1_rep3.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET2_rep3.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET3_rep3.bowtie2.sorted.bam --labels H9_d9_BPK25_input_rep1 H9_d9_BPK25_input_rep2 H9_d9_BPK25_input_rep3 H9_d9_BPK25_IP_CTCF_rep1 H9_d9_BPK25_IP_CTCF_rep2 H9_d9_BPK25_IP_CTCF_rep3 H9_d9_BPK25_IP_MBD3_rep1 H9_d9_BPK25_IP_MBD3_rep2 H9_d9_BPK25_IP_MBD3_rep3 H9_d9_BPK25_IP_TET1_rep1 H9_d9_BPK25_IP_TET1_rep2 H9_d9_BPK25_IP_TET1_rep3 H9_d9_BPK25_IP_TET2_rep1 H9_d9_BPK25_IP_TET2_rep2 H9_d9_BPK25_IP_TET2_rep3 H9_d9_BPK25_IP_TET3_rep1 H9_d9_BPK25_IP_TET3_rep2 H9_d9_BPK25_IP_TET3_rep3 -out multiBamSummary_results_for_H9_d9_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -p "max" --extendReads --ignoreDuplicates --minMappingQuality 30 --centerReads
# ## 用 plotCorrelation 和 plotPCA 分析样本的相似度
# plotCorrelation --corData multiBamSummary_results_for_H9_d9_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --corMethod pearson --whatToPlot heatmap --plotFile plotCorrelation_heatmap_plot_for_H9_d9_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.pdf --skipZeros --removeOutliers --outFileCorMatrix plotCorrelation_CorMatrix_for_H9_d9_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.txt --colorMap coolwarm --plotNumbers
# plotPCA --corData multiBamSummary_results_for_H9_d9_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --plotFile plotPCA_plot_for_H9_d9_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_H9_d9_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.txt
# plotPCA --rowCenter --corData multiBamSummary_results_for_H9_d9_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_bam_files.npz --plotFile plotPCA_plot_for_H9_d9_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.rowCenter.pdf --ntop 0 --outFileNameData plotPCA_stat_data_for_H9_d9_BPK25_IP_CTCF_MBD3_TET123_ChIP_Seq_all_samples.rowCenter.txt

# 上传并运行
nohup bash 5_make_cor_heatmap_and_PCA_plot_for_d9_BPK25_all_samples.sh >5_make_cor_heatmap_and_PCA_plot_for_d9_BPK25_all_samples.sh.runlog 2>&1 &
tail -f 5_make_cor_heatmap_and_PCA_plot_for_d9_BPK25_all_samples.sh.runlog

```

## fingerprint plot

### d049_DMSO_BPK25_CTCF_ChIP_Seq

```bash
# 参考： https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html
# 可以判断实验的成败、富集程度的好坏等
# 构造 5_make_fingerprint_plot_for_d049_DMSO_BPK25_CTCF_ChIP_Seq_samples.sh 内容如下：
# # d0_DMSO_IP_CTCF
# plotFingerprint -b hESC_H9_d0_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_input_rep3.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam --labels H9_d0_DMSO_input_rep1 H9_d0_DMSO_input_rep2 H9_d0_DMSO_input_rep3 H9_d0_DMSO_IP_CTCF_rep1 H9_d0_DMSO_IP_CTCF_rep2 H9_d0_DMSO_IP_CTCF_rep3 -o plotFingerprint_for_d0_DMSO_IP_CTCF_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d0_BPK25_IP_CTCF
# plotFingerprint -b hESC_H9_d0_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_input_rep3.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam --labels H9_d0_BPK25_input_rep1 H9_d0_BPK25_input_rep2 H9_d0_BPK25_input_rep3 H9_d0_BPK25_IP_CTCF_rep1 H9_d0_BPK25_IP_CTCF_rep2 H9_d0_BPK25_IP_CTCF_rep3 -o plotFingerprint_for_d0_BPK25_IP_CTCF_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d4_DMSO_IP_CTCF
# plotFingerprint -b hESC_H9_d4_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_input_rep3.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam --labels H9_d4_DMSO_input_rep1 H9_d4_DMSO_input_rep2 H9_d4_DMSO_input_rep3 H9_d4_DMSO_IP_CTCF_rep1 H9_d4_DMSO_IP_CTCF_rep2 H9_d4_DMSO_IP_CTCF_rep3 -o plotFingerprint_for_d4_DMSO_IP_CTCF_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d4_BPK25_IP_CTCF
# plotFingerprint -b hESC_H9_d4_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_input_rep3.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam --labels H9_d4_BPK25_input_rep1 H9_d4_BPK25_input_rep2 H9_d4_BPK25_input_rep3 H9_d4_BPK25_IP_CTCF_rep1 H9_d4_BPK25_IP_CTCF_rep2 H9_d4_BPK25_IP_CTCF_rep3 -o plotFingerprint_for_d4_BPK25_IP_CTCF_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d9_DMSO_IP_CTCF
# plotFingerprint -b hESC_H9_d9_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_input_rep3.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam --labels H9_d9_DMSO_input_rep1 H9_d9_DMSO_input_rep2 H9_d9_DMSO_input_rep3 H9_d9_DMSO_IP_CTCF_rep1 H9_d9_DMSO_IP_CTCF_rep2 H9_d9_DMSO_IP_CTCF_rep3 -o plotFingerprint_for_d9_DMSO_IP_CTCF_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d9_BPK25_IP_CTCF
# plotFingerprint -b hESC_H9_d9_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_input_rep3.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam --labels H9_d9_BPK25_input_rep1 H9_d9_BPK25_input_rep2 H9_d9_BPK25_input_rep3 H9_d9_BPK25_IP_CTCF_rep1 H9_d9_BPK25_IP_CTCF_rep2 H9_d9_BPK25_IP_CTCF_rep3 -o plotFingerprint_for_d9_BPK25_IP_CTCF_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"

# 上传并运行
nohup bash 5_make_fingerprint_plot_for_d049_DMSO_BPK25_CTCF_ChIP_Seq_samples.sh >5_make_fingerprint_plot_for_d049_DMSO_BPK25_CTCF_ChIP_Seq_samples.sh.runlog 2>&1 &
tail -f 5_make_fingerprint_plot_for_d049_DMSO_BPK25_CTCF_ChIP_Seq_samples.sh.runlog

```

### d049_DMSO_BPK25_MDB3_ChIP_Seq

```bash
# 参考： https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html
# 可以判断实验的成败、富集程度的好坏等
# 构造 5_make_fingerprint_plot_for_d049_DMSO_BPK25_MBD3_ChIP_Seq_samples.sh 内容如下：
# # d0_DMSO_IP_MBD3
# plotFingerprint -b hESC_H9_d0_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_input_rep3.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam --labels H9_d0_DMSO_input_rep1 H9_d0_DMSO_input_rep2 H9_d0_DMSO_input_rep3 H9_d0_DMSO_IP_MBD3_rep1 H9_d0_DMSO_IP_MBD3_rep2 H9_d0_DMSO_IP_MBD3_rep3 -o plotFingerprint_for_d0_DMSO_IP_MBD3_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d0_BPK25_IP_MBD3
# plotFingerprint -b hESC_H9_d0_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_input_rep3.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam --labels H9_d0_BPK25_input_rep1 H9_d0_BPK25_input_rep2 H9_d0_BPK25_input_rep3 H9_d0_BPK25_IP_MBD3_rep1 H9_d0_BPK25_IP_MBD3_rep2 H9_d0_BPK25_IP_MBD3_rep3 -o plotFingerprint_for_d0_BPK25_IP_MBD3_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d4_DMSO_IP_MBD3
# plotFingerprint -b hESC_H9_d4_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_input_rep3.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam --labels H9_d4_DMSO_input_rep1 H9_d4_DMSO_input_rep2 H9_d4_DMSO_input_rep3 H9_d4_DMSO_IP_MBD3_rep1 H9_d4_DMSO_IP_MBD3_rep2 H9_d4_DMSO_IP_MBD3_rep3 -o plotFingerprint_for_d4_DMSO_IP_MBD3_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d4_BPK25_IP_MBD3
# plotFingerprint -b hESC_H9_d4_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_input_rep3.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam --labels H9_d4_BPK25_input_rep1 H9_d4_BPK25_input_rep2 H9_d4_BPK25_input_rep3 H9_d4_BPK25_IP_MBD3_rep1 H9_d4_BPK25_IP_MBD3_rep2 H9_d4_BPK25_IP_MBD3_rep3 -o plotFingerprint_for_d4_BPK25_IP_MBD3_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d9_DMSO_IP_MBD3
# plotFingerprint -b hESC_H9_d9_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_input_rep3.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam --labels H9_d9_DMSO_input_rep1 H9_d9_DMSO_input_rep2 H9_d9_DMSO_input_rep3 H9_d9_DMSO_IP_MBD3_rep1 H9_d9_DMSO_IP_MBD3_rep2 H9_d9_DMSO_IP_MBD3_rep3 -o plotFingerprint_for_d9_DMSO_IP_MBD3_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d9_BPK25_IP_MBD3
# plotFingerprint -b hESC_H9_d9_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_input_rep3.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam --labels H9_d9_BPK25_input_rep1 H9_d9_BPK25_input_rep2 H9_d9_BPK25_input_rep3 H9_d9_BPK25_IP_MBD3_rep1 H9_d9_BPK25_IP_MBD3_rep2 H9_d9_BPK25_IP_MBD3_rep3 -o plotFingerprint_for_d9_BPK25_IP_MBD3_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"

# 上传并运行
nohup bash 5_make_fingerprint_plot_for_d049_DMSO_BPK25_MBD3_ChIP_Seq_samples.sh >5_make_fingerprint_plot_for_d049_DMSO_BPK25_MBD3_ChIP_Seq_samples.sh.runlog 2>&1 &
tail -f 5_make_fingerprint_plot_for_d049_DMSO_BPK25_MBD3_ChIP_Seq_samples.sh.runlog

```
### d049_DMSO_BPK25_TET1_ChIP_Seq

```bash
# 参考： https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html
# 可以判断实验的成败、富集程度的好坏等
# 构造 5_make_fingerprint_plot_for_d049_DMSO_BPK25_TET1_ChIP_Seq_samples.sh 内容如下：
# # d0_DMSO_IP_TET1
# plotFingerprint -b hESC_H9_d0_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_input_rep3.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET1_rep3.bowtie2.sorted.bam --labels H9_d0_DMSO_input_rep1 H9_d0_DMSO_input_rep2 H9_d0_DMSO_input_rep3 H9_d0_DMSO_IP_TET1_rep1 H9_d0_DMSO_IP_TET1_rep2 H9_d0_DMSO_IP_TET1_rep3 -o plotFingerprint_for_d0_DMSO_IP_TET1_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d0_BPK25_IP_TET1
# plotFingerprint -b hESC_H9_d0_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_input_rep3.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET1_rep3.bowtie2.sorted.bam --labels H9_d0_BPK25_input_rep1 H9_d0_BPK25_input_rep2 H9_d0_BPK25_input_rep3 H9_d0_BPK25_IP_TET1_rep1 H9_d0_BPK25_IP_TET1_rep2 H9_d0_BPK25_IP_TET1_rep3 -o plotFingerprint_for_d0_BPK25_IP_TET1_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d4_DMSO_IP_TET1
# plotFingerprint -b hESC_H9_d4_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_input_rep3.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET1_rep3.bowtie2.sorted.bam --labels H9_d4_DMSO_input_rep1 H9_d4_DMSO_input_rep2 H9_d4_DMSO_input_rep3 H9_d4_DMSO_IP_TET1_rep1 H9_d4_DMSO_IP_TET1_rep2 H9_d4_DMSO_IP_TET1_rep3 -o plotFingerprint_for_d4_DMSO_IP_TET1_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d4_BPK25_IP_TET1
# plotFingerprint -b hESC_H9_d4_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_input_rep3.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET1_rep3.bowtie2.sorted.bam --labels H9_d4_BPK25_input_rep1 H9_d4_BPK25_input_rep2 H9_d4_BPK25_input_rep3 H9_d4_BPK25_IP_TET1_rep1 H9_d4_BPK25_IP_TET1_rep2 H9_d4_BPK25_IP_TET1_rep3 -o plotFingerprint_for_d4_BPK25_IP_TET1_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d9_DMSO_IP_TET1
# plotFingerprint -b hESC_H9_d9_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_input_rep3.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET1_rep3.bowtie2.sorted.bam --labels H9_d9_DMSO_input_rep1 H9_d9_DMSO_input_rep2 H9_d9_DMSO_input_rep3 H9_d9_DMSO_IP_TET1_rep1 H9_d9_DMSO_IP_TET1_rep2 H9_d9_DMSO_IP_TET1_rep3 -o plotFingerprint_for_d9_DMSO_IP_TET1_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d9_BPK25_IP_TET1
# plotFingerprint -b hESC_H9_d9_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_input_rep3.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET1_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET1_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET1_rep3.bowtie2.sorted.bam --labels H9_d9_BPK25_input_rep1 H9_d9_BPK25_input_rep2 H9_d9_BPK25_input_rep3 H9_d9_BPK25_IP_TET1_rep1 H9_d9_BPK25_IP_TET1_rep2 H9_d9_BPK25_IP_TET1_rep3 -o plotFingerprint_for_d9_BPK25_IP_TET1_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"

# 上传并运行
nohup bash 5_make_fingerprint_plot_for_d049_DMSO_BPK25_TET1_ChIP_Seq_samples.sh >5_make_fingerprint_plot_for_d049_DMSO_BPK25_TET1_ChIP_Seq_samples.sh.runlog 2>&1 &
tail -f 5_make_fingerprint_plot_for_d049_DMSO_BPK25_TET1_ChIP_Seq_samples.sh.runlog

```
### d049_DMSO_BPK25_TET2_ChIP_Seq

```bash
# 参考： https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html
# 可以判断实验的成败、富集程度的好坏等
# 构造 5_make_fingerprint_plot_for_d049_DMSO_BPK25_TET2_ChIP_Seq_samples.sh 内容如下：
# # d0_DMSO_IP_TET2
# plotFingerprint -b hESC_H9_d0_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_input_rep3.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET2_rep3.bowtie2.sorted.bam --labels H9_d0_DMSO_input_rep1 H9_d0_DMSO_input_rep2 H9_d0_DMSO_input_rep3 H9_d0_DMSO_IP_TET2_rep1 H9_d0_DMSO_IP_TET2_rep2 H9_d0_DMSO_IP_TET2_rep3 -o plotFingerprint_for_d0_DMSO_IP_TET2_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d0_BPK25_IP_TET2
# plotFingerprint -b hESC_H9_d0_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_input_rep3.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET2_rep3.bowtie2.sorted.bam --labels H9_d0_BPK25_input_rep1 H9_d0_BPK25_input_rep2 H9_d0_BPK25_input_rep3 H9_d0_BPK25_IP_TET2_rep1 H9_d0_BPK25_IP_TET2_rep2 H9_d0_BPK25_IP_TET2_rep3 -o plotFingerprint_for_d0_BPK25_IP_TET2_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d4_DMSO_IP_TET2
# plotFingerprint -b hESC_H9_d4_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_input_rep3.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET2_rep3.bowtie2.sorted.bam --labels H9_d4_DMSO_input_rep1 H9_d4_DMSO_input_rep2 H9_d4_DMSO_input_rep3 H9_d4_DMSO_IP_TET2_rep1 H9_d4_DMSO_IP_TET2_rep2 H9_d4_DMSO_IP_TET2_rep3 -o plotFingerprint_for_d4_DMSO_IP_TET2_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d4_BPK25_IP_TET2
# plotFingerprint -b hESC_H9_d4_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_input_rep3.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET2_rep3.bowtie2.sorted.bam --labels H9_d4_BPK25_input_rep1 H9_d4_BPK25_input_rep2 H9_d4_BPK25_input_rep3 H9_d4_BPK25_IP_TET2_rep1 H9_d4_BPK25_IP_TET2_rep2 H9_d4_BPK25_IP_TET2_rep3 -o plotFingerprint_for_d4_BPK25_IP_TET2_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d9_DMSO_IP_TET2
# plotFingerprint -b hESC_H9_d9_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_input_rep3.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET2_rep3.bowtie2.sorted.bam --labels H9_d9_DMSO_input_rep1 H9_d9_DMSO_input_rep2 H9_d9_DMSO_input_rep3 H9_d9_DMSO_IP_TET2_rep1 H9_d9_DMSO_IP_TET2_rep2 H9_d9_DMSO_IP_TET2_rep3 -o plotFingerprint_for_d9_DMSO_IP_TET2_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d9_BPK25_IP_TET2
# plotFingerprint -b hESC_H9_d9_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_input_rep3.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET2_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET2_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET2_rep3.bowtie2.sorted.bam --labels H9_d9_BPK25_input_rep1 H9_d9_BPK25_input_rep2 H9_d9_BPK25_input_rep3 H9_d9_BPK25_IP_TET2_rep1 H9_d9_BPK25_IP_TET2_rep2 H9_d9_BPK25_IP_TET2_rep3 -o plotFingerprint_for_d9_BPK25_IP_TET2_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"

# 上传并运行
nohup bash 5_make_fingerprint_plot_for_d049_DMSO_BPK25_TET2_ChIP_Seq_samples.sh >5_make_fingerprint_plot_for_d049_DMSO_BPK25_TET2_ChIP_Seq_samples.sh.runlog 2>&1 &
tail -f 5_make_fingerprint_plot_for_d049_DMSO_BPK25_TET2_ChIP_Seq_samples.sh.runlog

```
### d049_DMSO_BPK25_TET3_ChIP_Seq

```bash
# 参考： https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html
# 可以判断实验的成败、富集程度的好坏等
# 构造 5_make_fingerprint_plot_for_d049_DMSO_BPK25_TET3_ChIP_Seq_samples.sh 内容如下：
# # d0_DMSO_IP_TET3
# plotFingerprint -b hESC_H9_d0_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_input_rep3.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d0_DMSO_IP_TET3_rep3.bowtie2.sorted.bam --labels H9_d0_DMSO_input_rep1 H9_d0_DMSO_input_rep2 H9_d0_DMSO_input_rep3 H9_d0_DMSO_IP_TET3_rep1 H9_d0_DMSO_IP_TET3_rep2 H9_d0_DMSO_IP_TET3_rep3 -o plotFingerprint_for_d0_DMSO_IP_TET3_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d0_BPK25_IP_TET3
# plotFingerprint -b hESC_H9_d0_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_input_rep3.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d0_BPK25_IP_TET3_rep3.bowtie2.sorted.bam --labels H9_d0_BPK25_input_rep1 H9_d0_BPK25_input_rep2 H9_d0_BPK25_input_rep3 H9_d0_BPK25_IP_TET3_rep1 H9_d0_BPK25_IP_TET3_rep2 H9_d0_BPK25_IP_TET3_rep3 -o plotFingerprint_for_d0_BPK25_IP_TET3_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d4_DMSO_IP_TET3
# plotFingerprint -b hESC_H9_d4_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_input_rep3.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d4_DMSO_IP_TET3_rep3.bowtie2.sorted.bam --labels H9_d4_DMSO_input_rep1 H9_d4_DMSO_input_rep2 H9_d4_DMSO_input_rep3 H9_d4_DMSO_IP_TET3_rep1 H9_d4_DMSO_IP_TET3_rep2 H9_d4_DMSO_IP_TET3_rep3 -o plotFingerprint_for_d4_DMSO_IP_TET3_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d4_BPK25_IP_TET3
# plotFingerprint -b hESC_H9_d4_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_input_rep3.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d4_BPK25_IP_TET3_rep3.bowtie2.sorted.bam --labels H9_d4_BPK25_input_rep1 H9_d4_BPK25_input_rep2 H9_d4_BPK25_input_rep3 H9_d4_BPK25_IP_TET3_rep1 H9_d4_BPK25_IP_TET3_rep2 H9_d4_BPK25_IP_TET3_rep3 -o plotFingerprint_for_d4_BPK25_IP_TET3_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d9_DMSO_IP_TET3
# plotFingerprint -b hESC_H9_d9_DMSO_input_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_input_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_input_rep3.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d9_DMSO_IP_TET3_rep3.bowtie2.sorted.bam --labels H9_d9_DMSO_input_rep1 H9_d9_DMSO_input_rep2 H9_d9_DMSO_input_rep3 H9_d9_DMSO_IP_TET3_rep1 H9_d9_DMSO_IP_TET3_rep2 H9_d9_DMSO_IP_TET3_rep3 -o plotFingerprint_for_d9_DMSO_IP_TET3_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"
# # d9_BPK25_IP_TET3
# plotFingerprint -b hESC_H9_d9_BPK25_input_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_input_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_input_rep3.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET3_rep1.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET3_rep2.bowtie2.sorted.bam hESC_H9_d9_BPK25_IP_TET3_rep3.bowtie2.sorted.bam --labels H9_d9_BPK25_input_rep1 H9_d9_BPK25_input_rep2 H9_d9_BPK25_input_rep3 H9_d9_BPK25_IP_TET3_rep1 H9_d9_BPK25_IP_TET3_rep2 H9_d9_BPK25_IP_TET3_rep3 -o plotFingerprint_for_d9_BPK25_IP_TET3_ChIP_Seq_samples.pdf --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed --ignoreDuplicates --centerReads --minMappingQuality 30 --skipZeros -p "max"

# 上传并运行
nohup bash 5_make_fingerprint_plot_for_d049_DMSO_BPK25_TET3_ChIP_Seq_samples.sh >5_make_fingerprint_plot_for_d049_DMSO_BPK25_TET3_ChIP_Seq_samples.sh.runlog 2>&1 &
tail -f 5_make_fingerprint_plot_for_d049_DMSO_BPK25_TET3_ChIP_Seq_samples.sh.runlog

```



# macs2 peak calling

## peak calling for each repetition separately 

### d049_DMSO_BPK25_CTCF_ChIP_Seq

```bash
# data preperation
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir macs2_peak_calling_1_for_CTCF_each_repeat
cd macs2_peak_calling_1_for_CTCF_each_repeat
# 上传 bam 连接构造脚本并运行
bash 6_make_H9_d0_DMSO_BPK25_CTCF_ChIP_Seq_sorted_bam_file_links.sh
bash 6_make_H9_d4_DMSO_BPK25_CTCF_ChIP_Seq_sorted_bam_file_links.sh
bash 6_make_H9_d9_DMSO_BPK25_CTCF_ChIP_Seq_sorted_bam_file_links.sh
ls -1

# d0_DMSO
nohup macs2 callpeak -t hESC_H9_d0_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam -c hESC_H9_d0_DMSO_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_DMSO_IP_CTCF_rep1_TvsC -n d0_DMSO_IP_CTCF_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_DMSO_IP_CTCF_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d0_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam -c hESC_H9_d0_DMSO_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_DMSO_IP_CTCF_rep2_TvsC -n d0_DMSO_IP_CTCF_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_DMSO_IP_CTCF_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d0_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam -c hESC_H9_d0_DMSO_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_DMSO_IP_CTCF_rep3_TvsC -n d0_DMSO_IP_CTCF_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_DMSO_IP_CTCF_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d0_DMSO_IP_CTCF_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d0_DMSO_IP_CTCF_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d0_DMSO_IP_CTCF_rep3_TvsC.runlog
# d0_BPK25
nohup macs2 callpeak -t hESC_H9_d0_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam -c hESC_H9_d0_BPK25_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_BPK25_IP_CTCF_rep1_TvsC -n d0_BPK25_IP_CTCF_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_BPK25_IP_CTCF_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d0_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam -c hESC_H9_d0_BPK25_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_BPK25_IP_CTCF_rep2_TvsC -n d0_BPK25_IP_CTCF_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_BPK25_IP_CTCF_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d0_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam -c hESC_H9_d0_BPK25_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_BPK25_IP_CTCF_rep3_TvsC -n d0_BPK25_IP_CTCF_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_BPK25_IP_CTCF_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d0_BPK25_IP_CTCF_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d0_BPK25_IP_CTCF_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d0_BPK25_IP_CTCF_rep3_TvsC.runlog
# d4_DMSO
nohup macs2 callpeak -t hESC_H9_d4_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam -c hESC_H9_d4_DMSO_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_DMSO_IP_CTCF_rep1_TvsC -n d4_DMSO_IP_CTCF_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_DMSO_IP_CTCF_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d4_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam -c hESC_H9_d4_DMSO_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_DMSO_IP_CTCF_rep2_TvsC -n d4_DMSO_IP_CTCF_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_DMSO_IP_CTCF_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d4_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam -c hESC_H9_d4_DMSO_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_DMSO_IP_CTCF_rep3_TvsC -n d4_DMSO_IP_CTCF_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_DMSO_IP_CTCF_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d4_DMSO_IP_CTCF_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d4_DMSO_IP_CTCF_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d4_DMSO_IP_CTCF_rep3_TvsC.runlog
# d4_BPK25
nohup macs2 callpeak -t hESC_H9_d4_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam -c hESC_H9_d4_BPK25_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_BPK25_IP_CTCF_rep1_TvsC -n d4_BPK25_IP_CTCF_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_BPK25_IP_CTCF_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d4_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam -c hESC_H9_d4_BPK25_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_BPK25_IP_CTCF_rep2_TvsC -n d4_BPK25_IP_CTCF_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_BPK25_IP_CTCF_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d4_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam -c hESC_H9_d4_BPK25_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_BPK25_IP_CTCF_rep3_TvsC -n d4_BPK25_IP_CTCF_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_BPK25_IP_CTCF_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d4_BPK25_IP_CTCF_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d4_BPK25_IP_CTCF_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d4_BPK25_IP_CTCF_rep3_TvsC.runlog
# d9_DMSO
nohup macs2 callpeak -t hESC_H9_d9_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam -c hESC_H9_d9_DMSO_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_DMSO_IP_CTCF_rep1_TvsC -n d9_DMSO_IP_CTCF_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_DMSO_IP_CTCF_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d9_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam -c hESC_H9_d9_DMSO_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_DMSO_IP_CTCF_rep2_TvsC -n d9_DMSO_IP_CTCF_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_DMSO_IP_CTCF_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d9_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam -c hESC_H9_d9_DMSO_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_DMSO_IP_CTCF_rep3_TvsC -n d9_DMSO_IP_CTCF_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_DMSO_IP_CTCF_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d9_DMSO_IP_CTCF_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d9_DMSO_IP_CTCF_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d9_DMSO_IP_CTCF_rep3_TvsC.runlog
# d9_BPK25
nohup macs2 callpeak -t hESC_H9_d9_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam -c hESC_H9_d9_BPK25_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_BPK25_IP_CTCF_rep1_TvsC -n d9_BPK25_IP_CTCF_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_BPK25_IP_CTCF_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d9_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam -c hESC_H9_d9_BPK25_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_BPK25_IP_CTCF_rep2_TvsC -n d9_BPK25_IP_CTCF_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_BPK25_IP_CTCF_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d9_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam -c hESC_H9_d9_BPK25_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_BPK25_IP_CTCF_rep3_TvsC -n d9_BPK25_IP_CTCF_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_BPK25_IP_CTCF_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d9_BPK25_IP_CTCF_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d9_BPK25_IP_CTCF_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d9_BPK25_IP_CTCF_rep3_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_1_for_CTCF_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./d0_BPK25_IP_CTCF_rep1_TvsC/d0_BPK25_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# 206
# ./d0_BPK25_IP_CTCF_rep2_TvsC/d0_BPK25_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# 171
# ./d0_BPK25_IP_CTCF_rep3_TvsC/d0_BPK25_IP_CTCF_rep3_TvsC_peaks.narrowPeak
# 312
# ./d0_DMSO_IP_CTCF_rep1_TvsC/d0_DMSO_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# 501
# ./d0_DMSO_IP_CTCF_rep2_TvsC/d0_DMSO_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# 553
# ./d0_DMSO_IP_CTCF_rep3_TvsC/d0_DMSO_IP_CTCF_rep3_TvsC_peaks.narrowPeak
# 562
# ./d4_BPK25_IP_CTCF_rep1_TvsC/d4_BPK25_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# 4005
# ./d4_BPK25_IP_CTCF_rep2_TvsC/d4_BPK25_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# 4322
# ./d4_BPK25_IP_CTCF_rep3_TvsC/d4_BPK25_IP_CTCF_rep3_TvsC_peaks.narrowPeak
# 4019
# ./d4_DMSO_IP_CTCF_rep1_TvsC/d4_DMSO_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# 4291
# ./d4_DMSO_IP_CTCF_rep2_TvsC/d4_DMSO_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# 4994
# ./d4_DMSO_IP_CTCF_rep3_TvsC/d4_DMSO_IP_CTCF_rep3_TvsC_peaks.narrowPeak
# 5277
# ./d9_BPK25_IP_CTCF_rep1_TvsC/d9_BPK25_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# 806
# ./d9_BPK25_IP_CTCF_rep2_TvsC/d9_BPK25_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# 946
# ./d9_BPK25_IP_CTCF_rep3_TvsC/d9_BPK25_IP_CTCF_rep3_TvsC_peaks.narrowPeak
# 1091
# ./d9_DMSO_IP_CTCF_rep1_TvsC/d9_DMSO_IP_CTCF_rep1_TvsC_peaks.narrowPeak
# 923
# ./d9_DMSO_IP_CTCF_rep2_TvsC/d9_DMSO_IP_CTCF_rep2_TvsC_peaks.narrowPeak
# 561
# ./d9_DMSO_IP_CTCF_rep3_TvsC/d9_DMSO_IP_CTCF_rep3_TvsC_peaks.narrowPeak
# 711

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./d0_BPK25_IP_CTCF_rep1_TvsC/d0_BPK25_IP_CTCF_rep1_TvsC_summits.bed
# 247
# ./d0_BPK25_IP_CTCF_rep2_TvsC/d0_BPK25_IP_CTCF_rep2_TvsC_summits.bed
# 215
# ./d0_BPK25_IP_CTCF_rep3_TvsC/d0_BPK25_IP_CTCF_rep3_TvsC_summits.bed
# 377
# ./d0_DMSO_IP_CTCF_rep1_TvsC/d0_DMSO_IP_CTCF_rep1_TvsC_summits.bed
# 636
# ./d0_DMSO_IP_CTCF_rep2_TvsC/d0_DMSO_IP_CTCF_rep2_TvsC_summits.bed
# 686
# ./d0_DMSO_IP_CTCF_rep3_TvsC/d0_DMSO_IP_CTCF_rep3_TvsC_summits.bed
# 694
# ./d4_BPK25_IP_CTCF_rep1_TvsC/d4_BPK25_IP_CTCF_rep1_TvsC_summits.bed
# 4774
# ./d4_BPK25_IP_CTCF_rep2_TvsC/d4_BPK25_IP_CTCF_rep2_TvsC_summits.bed
# 5178
# ./d4_BPK25_IP_CTCF_rep3_TvsC/d4_BPK25_IP_CTCF_rep3_TvsC_summits.bed
# 4815
# ./d4_DMSO_IP_CTCF_rep1_TvsC/d4_DMSO_IP_CTCF_rep1_TvsC_summits.bed
# 4834
# ./d4_DMSO_IP_CTCF_rep2_TvsC/d4_DMSO_IP_CTCF_rep2_TvsC_summits.bed
# 5694
# ./d4_DMSO_IP_CTCF_rep3_TvsC/d4_DMSO_IP_CTCF_rep3_TvsC_summits.bed
# 5909
# ./d9_BPK25_IP_CTCF_rep1_TvsC/d9_BPK25_IP_CTCF_rep1_TvsC_summits.bed
# 990
# ./d9_BPK25_IP_CTCF_rep2_TvsC/d9_BPK25_IP_CTCF_rep2_TvsC_summits.bed
# 1154
# ./d9_BPK25_IP_CTCF_rep3_TvsC/d9_BPK25_IP_CTCF_rep3_TvsC_summits.bed
# 1297
# ./d9_DMSO_IP_CTCF_rep1_TvsC/d9_DMSO_IP_CTCF_rep1_TvsC_summits.bed
# 979
# ./d9_DMSO_IP_CTCF_rep2_TvsC/d9_DMSO_IP_CTCF_rep2_TvsC_summits.bed
# 619
# ./d9_DMSO_IP_CTCF_rep3_TvsC/d9_DMSO_IP_CTCF_rep3_TvsC_summits.bed
# 765

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_BPK25_DMSO_IP_CTCF_ChIP_Seq_samples_macs2_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done

```



### d049_DMSO_BPK25_MBD3_ChIP_Seq

```bash
# data preperation
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir macs2_peak_calling_1_for_MBD3_each_repeat
cd macs2_peak_calling_1_for_MBD3_each_repeat
# 上传 bam 连接构造脚本并运行
bash 6_make_H9_d0_DMSO_BPK25_MBD3_ChIP_Seq_sorted_bam_file_links.sh
bash 6_make_H9_d4_DMSO_BPK25_MBD3_ChIP_Seq_sorted_bam_file_links.sh
bash 6_make_H9_d9_DMSO_BPK25_MBD3_ChIP_Seq_sorted_bam_file_links.sh
ls -1

# d0_DMSO
nohup macs2 callpeak -t hESC_H9_d0_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam -c hESC_H9_d0_DMSO_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_DMSO_IP_MBD3_rep1_TvsC -n d0_DMSO_IP_MBD3_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_DMSO_IP_MBD3_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d0_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam -c hESC_H9_d0_DMSO_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_DMSO_IP_MBD3_rep2_TvsC -n d0_DMSO_IP_MBD3_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_DMSO_IP_MBD3_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d0_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam -c hESC_H9_d0_DMSO_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_DMSO_IP_MBD3_rep3_TvsC -n d0_DMSO_IP_MBD3_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_DMSO_IP_MBD3_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d0_DMSO_IP_MBD3_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d0_DMSO_IP_MBD3_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d0_DMSO_IP_MBD3_rep3_TvsC.runlog
# d0_BPK25
nohup macs2 callpeak -t hESC_H9_d0_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam -c hESC_H9_d0_BPK25_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_BPK25_IP_MBD3_rep1_TvsC -n d0_BPK25_IP_MBD3_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_BPK25_IP_MBD3_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d0_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam -c hESC_H9_d0_BPK25_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_BPK25_IP_MBD3_rep2_TvsC -n d0_BPK25_IP_MBD3_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_BPK25_IP_MBD3_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d0_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam -c hESC_H9_d0_BPK25_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_BPK25_IP_MBD3_rep3_TvsC -n d0_BPK25_IP_MBD3_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_BPK25_IP_MBD3_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d0_BPK25_IP_MBD3_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d0_BPK25_IP_MBD3_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d0_BPK25_IP_MBD3_rep3_TvsC.runlog
# d4_DMSO
nohup macs2 callpeak -t hESC_H9_d4_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam -c hESC_H9_d4_DMSO_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_DMSO_IP_MBD3_rep1_TvsC -n d4_DMSO_IP_MBD3_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_DMSO_IP_MBD3_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d4_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam -c hESC_H9_d4_DMSO_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_DMSO_IP_MBD3_rep2_TvsC -n d4_DMSO_IP_MBD3_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_DMSO_IP_MBD3_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d4_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam -c hESC_H9_d4_DMSO_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_DMSO_IP_MBD3_rep3_TvsC -n d4_DMSO_IP_MBD3_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_DMSO_IP_MBD3_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d4_DMSO_IP_MBD3_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d4_DMSO_IP_MBD3_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d4_DMSO_IP_MBD3_rep3_TvsC.runlog
# d4_BPK25
nohup macs2 callpeak -t hESC_H9_d4_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam -c hESC_H9_d4_BPK25_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_BPK25_IP_MBD3_rep1_TvsC -n d4_BPK25_IP_MBD3_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_BPK25_IP_MBD3_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d4_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam -c hESC_H9_d4_BPK25_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_BPK25_IP_MBD3_rep2_TvsC -n d4_BPK25_IP_MBD3_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_BPK25_IP_MBD3_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d4_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam -c hESC_H9_d4_BPK25_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_BPK25_IP_MBD3_rep3_TvsC -n d4_BPK25_IP_MBD3_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_BPK25_IP_MBD3_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d4_BPK25_IP_MBD3_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d4_BPK25_IP_MBD3_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d4_BPK25_IP_MBD3_rep3_TvsC.runlog
# d9_DMSO
nohup macs2 callpeak -t hESC_H9_d9_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam -c hESC_H9_d9_DMSO_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_DMSO_IP_MBD3_rep1_TvsC -n d9_DMSO_IP_MBD3_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_DMSO_IP_MBD3_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d9_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam -c hESC_H9_d9_DMSO_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_DMSO_IP_MBD3_rep2_TvsC -n d9_DMSO_IP_MBD3_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_DMSO_IP_MBD3_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d9_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam -c hESC_H9_d9_DMSO_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_DMSO_IP_MBD3_rep3_TvsC -n d9_DMSO_IP_MBD3_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_DMSO_IP_MBD3_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d9_DMSO_IP_MBD3_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d9_DMSO_IP_MBD3_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d9_DMSO_IP_MBD3_rep3_TvsC.runlog
# d9_BPK25
nohup macs2 callpeak -t hESC_H9_d9_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam -c hESC_H9_d9_BPK25_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_BPK25_IP_MBD3_rep1_TvsC -n d9_BPK25_IP_MBD3_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_BPK25_IP_MBD3_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d9_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam -c hESC_H9_d9_BPK25_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_BPK25_IP_MBD3_rep2_TvsC -n d9_BPK25_IP_MBD3_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_BPK25_IP_MBD3_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d9_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam -c hESC_H9_d9_BPK25_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_BPK25_IP_MBD3_rep3_TvsC -n d9_BPK25_IP_MBD3_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_BPK25_IP_MBD3_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d9_BPK25_IP_MBD3_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d9_BPK25_IP_MBD3_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d9_BPK25_IP_MBD3_rep3_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_1_for_MBD3_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./d0_BPK25_IP_MBD3_rep1_TvsC/d0_BPK25_IP_MBD3_rep1_TvsC_peaks.narrowPeak
# 143
# ./d0_BPK25_IP_MBD3_rep2_TvsC/d0_BPK25_IP_MBD3_rep2_TvsC_peaks.narrowPeak
# 178
# ./d0_BPK25_IP_MBD3_rep3_TvsC/d0_BPK25_IP_MBD3_rep3_TvsC_peaks.narrowPeak
# 118
# ./d0_DMSO_IP_MBD3_rep1_TvsC/d0_DMSO_IP_MBD3_rep1_TvsC_peaks.narrowPeak
# 1817
# ./d0_DMSO_IP_MBD3_rep2_TvsC/d0_DMSO_IP_MBD3_rep2_TvsC_peaks.narrowPeak
# 1885
# ./d0_DMSO_IP_MBD3_rep3_TvsC/d0_DMSO_IP_MBD3_rep3_TvsC_peaks.narrowPeak
# 2139
# ./d4_BPK25_IP_MBD3_rep1_TvsC/d4_BPK25_IP_MBD3_rep1_TvsC_peaks.narrowPeak
# 3346
# ./d4_BPK25_IP_MBD3_rep2_TvsC/d4_BPK25_IP_MBD3_rep2_TvsC_peaks.narrowPeak
# 4296
# ./d4_BPK25_IP_MBD3_rep3_TvsC/d4_BPK25_IP_MBD3_rep3_TvsC_peaks.narrowPeak
# 2543
# ./d4_DMSO_IP_MBD3_rep1_TvsC/d4_DMSO_IP_MBD3_rep1_TvsC_peaks.narrowPeak
# 6578
# ./d4_DMSO_IP_MBD3_rep2_TvsC/d4_DMSO_IP_MBD3_rep2_TvsC_peaks.narrowPeak
# 12160
# ./d4_DMSO_IP_MBD3_rep3_TvsC/d4_DMSO_IP_MBD3_rep3_TvsC_peaks.narrowPeak
# 8042
# ./d9_BPK25_IP_MBD3_rep1_TvsC/d9_BPK25_IP_MBD3_rep1_TvsC_peaks.narrowPeak
# 3560
# ./d9_BPK25_IP_MBD3_rep2_TvsC/d9_BPK25_IP_MBD3_rep2_TvsC_peaks.narrowPeak
# 5091
# ./d9_BPK25_IP_MBD3_rep3_TvsC/d9_BPK25_IP_MBD3_rep3_TvsC_peaks.narrowPeak
# 6900
# ./d9_DMSO_IP_MBD3_rep1_TvsC/d9_DMSO_IP_MBD3_rep1_TvsC_peaks.narrowPeak
# 5143
# ./d9_DMSO_IP_MBD3_rep2_TvsC/d9_DMSO_IP_MBD3_rep2_TvsC_peaks.narrowPeak
# 4580
# ./d9_DMSO_IP_MBD3_rep3_TvsC/d9_DMSO_IP_MBD3_rep3_TvsC_peaks.narrowPeak
# 5793

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./d0_BPK25_IP_MBD3_rep1_TvsC/d0_BPK25_IP_MBD3_rep1_TvsC_summits.bed
# 204
# ./d0_BPK25_IP_MBD3_rep2_TvsC/d0_BPK25_IP_MBD3_rep2_TvsC_summits.bed
# 242
# ./d0_BPK25_IP_MBD3_rep3_TvsC/d0_BPK25_IP_MBD3_rep3_TvsC_summits.bed
# 178
# ./d0_DMSO_IP_MBD3_rep1_TvsC/d0_DMSO_IP_MBD3_rep1_TvsC_summits.bed
# 2091
# ./d0_DMSO_IP_MBD3_rep2_TvsC/d0_DMSO_IP_MBD3_rep2_TvsC_summits.bed
# 2132
# ./d0_DMSO_IP_MBD3_rep3_TvsC/d0_DMSO_IP_MBD3_rep3_TvsC_summits.bed
# 2507
# ./d4_BPK25_IP_MBD3_rep1_TvsC/d4_BPK25_IP_MBD3_rep1_TvsC_summits.bed
# 3969
# ./d4_BPK25_IP_MBD3_rep2_TvsC/d4_BPK25_IP_MBD3_rep2_TvsC_summits.bed
# 5021
# ./d4_BPK25_IP_MBD3_rep3_TvsC/d4_BPK25_IP_MBD3_rep3_TvsC_summits.bed
# 3096
# ./d4_DMSO_IP_MBD3_rep1_TvsC/d4_DMSO_IP_MBD3_rep1_TvsC_summits.bed
# 6894
# ./d4_DMSO_IP_MBD3_rep2_TvsC/d4_DMSO_IP_MBD3_rep2_TvsC_summits.bed
# 12630
# ./d4_DMSO_IP_MBD3_rep3_TvsC/d4_DMSO_IP_MBD3_rep3_TvsC_summits.bed
# 8265
# ./d9_BPK25_IP_MBD3_rep1_TvsC/d9_BPK25_IP_MBD3_rep1_TvsC_summits.bed
# 3918
# ./d9_BPK25_IP_MBD3_rep2_TvsC/d9_BPK25_IP_MBD3_rep2_TvsC_summits.bed
# 5472
# ./d9_BPK25_IP_MBD3_rep3_TvsC/d9_BPK25_IP_MBD3_rep3_TvsC_summits.bed
# 7399
# ./d9_DMSO_IP_MBD3_rep1_TvsC/d9_DMSO_IP_MBD3_rep1_TvsC_summits.bed
# 5576
# ./d9_DMSO_IP_MBD3_rep2_TvsC/d9_DMSO_IP_MBD3_rep2_TvsC_summits.bed
# 5034
# ./d9_DMSO_IP_MBD3_rep3_TvsC/d9_DMSO_IP_MBD3_rep3_TvsC_summits.bed
# 6330

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_BPK25_DMSO_IP_MBD3_ChIP_Seq_samples_macs2_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done

```

### d049_DMSO_BPK25_TET1_ChIP_Seq

```bash
# data preperation
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir macs2_peak_calling_1_for_TET1_each_repeat
cd macs2_peak_calling_1_for_TET1_each_repeat
# 上传 bam 连接构造脚本并运行
bash 6_make_H9_d0_DMSO_BPK25_TET1_ChIP_Seq_sorted_bam_file_links.sh
bash 6_make_H9_d4_DMSO_BPK25_TET1_ChIP_Seq_sorted_bam_file_links.sh
bash 6_make_H9_d9_DMSO_BPK25_TET1_ChIP_Seq_sorted_bam_file_links.sh
ls -1

# d0_DMSO
nohup macs2 callpeak -t hESC_H9_d0_DMSO_IP_TET1_rep1.bowtie2.sorted.bam -c hESC_H9_d0_DMSO_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_DMSO_IP_TET1_rep1_TvsC -n d0_DMSO_IP_TET1_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_DMSO_IP_TET1_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d0_DMSO_IP_TET1_rep2.bowtie2.sorted.bam -c hESC_H9_d0_DMSO_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_DMSO_IP_TET1_rep2_TvsC -n d0_DMSO_IP_TET1_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_DMSO_IP_TET1_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d0_DMSO_IP_TET1_rep3.bowtie2.sorted.bam -c hESC_H9_d0_DMSO_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_DMSO_IP_TET1_rep3_TvsC -n d0_DMSO_IP_TET1_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_DMSO_IP_TET1_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d0_DMSO_IP_TET1_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d0_DMSO_IP_TET1_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d0_DMSO_IP_TET1_rep3_TvsC.runlog
# d0_BPK25
nohup macs2 callpeak -t hESC_H9_d0_BPK25_IP_TET1_rep1.bowtie2.sorted.bam -c hESC_H9_d0_BPK25_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_BPK25_IP_TET1_rep1_TvsC -n d0_BPK25_IP_TET1_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_BPK25_IP_TET1_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d0_BPK25_IP_TET1_rep2.bowtie2.sorted.bam -c hESC_H9_d0_BPK25_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_BPK25_IP_TET1_rep2_TvsC -n d0_BPK25_IP_TET1_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_BPK25_IP_TET1_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d0_BPK25_IP_TET1_rep3.bowtie2.sorted.bam -c hESC_H9_d0_BPK25_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_BPK25_IP_TET1_rep3_TvsC -n d0_BPK25_IP_TET1_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_BPK25_IP_TET1_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d0_BPK25_IP_TET1_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d0_BPK25_IP_TET1_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d0_BPK25_IP_TET1_rep3_TvsC.runlog
# d4_DMSO
nohup macs2 callpeak -t hESC_H9_d4_DMSO_IP_TET1_rep1.bowtie2.sorted.bam -c hESC_H9_d4_DMSO_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_DMSO_IP_TET1_rep1_TvsC -n d4_DMSO_IP_TET1_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_DMSO_IP_TET1_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d4_DMSO_IP_TET1_rep2.bowtie2.sorted.bam -c hESC_H9_d4_DMSO_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_DMSO_IP_TET1_rep2_TvsC -n d4_DMSO_IP_TET1_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_DMSO_IP_TET1_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d4_DMSO_IP_TET1_rep3.bowtie2.sorted.bam -c hESC_H9_d4_DMSO_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_DMSO_IP_TET1_rep3_TvsC -n d4_DMSO_IP_TET1_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_DMSO_IP_TET1_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d4_DMSO_IP_TET1_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d4_DMSO_IP_TET1_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d4_DMSO_IP_TET1_rep3_TvsC.runlog
# d4_BPK25
nohup macs2 callpeak -t hESC_H9_d4_BPK25_IP_TET1_rep1.bowtie2.sorted.bam -c hESC_H9_d4_BPK25_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_BPK25_IP_TET1_rep1_TvsC -n d4_BPK25_IP_TET1_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_BPK25_IP_TET1_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d4_BPK25_IP_TET1_rep2.bowtie2.sorted.bam -c hESC_H9_d4_BPK25_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_BPK25_IP_TET1_rep2_TvsC -n d4_BPK25_IP_TET1_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_BPK25_IP_TET1_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d4_BPK25_IP_TET1_rep3.bowtie2.sorted.bam -c hESC_H9_d4_BPK25_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_BPK25_IP_TET1_rep3_TvsC -n d4_BPK25_IP_TET1_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_BPK25_IP_TET1_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d4_BPK25_IP_TET1_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d4_BPK25_IP_TET1_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d4_BPK25_IP_TET1_rep3_TvsC.runlog
# d9_DMSO
nohup macs2 callpeak -t hESC_H9_d9_DMSO_IP_TET1_rep1.bowtie2.sorted.bam -c hESC_H9_d9_DMSO_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_DMSO_IP_TET1_rep1_TvsC -n d9_DMSO_IP_TET1_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_DMSO_IP_TET1_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d9_DMSO_IP_TET1_rep2.bowtie2.sorted.bam -c hESC_H9_d9_DMSO_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_DMSO_IP_TET1_rep2_TvsC -n d9_DMSO_IP_TET1_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_DMSO_IP_TET1_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d9_DMSO_IP_TET1_rep3.bowtie2.sorted.bam -c hESC_H9_d9_DMSO_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_DMSO_IP_TET1_rep3_TvsC -n d9_DMSO_IP_TET1_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_DMSO_IP_TET1_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d9_DMSO_IP_TET1_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d9_DMSO_IP_TET1_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d9_DMSO_IP_TET1_rep3_TvsC.runlog
# d9_BPK25
nohup macs2 callpeak -t hESC_H9_d9_BPK25_IP_TET1_rep1.bowtie2.sorted.bam -c hESC_H9_d9_BPK25_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_BPK25_IP_TET1_rep1_TvsC -n d9_BPK25_IP_TET1_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_BPK25_IP_TET1_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d9_BPK25_IP_TET1_rep2.bowtie2.sorted.bam -c hESC_H9_d9_BPK25_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_BPK25_IP_TET1_rep2_TvsC -n d9_BPK25_IP_TET1_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_BPK25_IP_TET1_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d9_BPK25_IP_TET1_rep3.bowtie2.sorted.bam -c hESC_H9_d9_BPK25_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_BPK25_IP_TET1_rep3_TvsC -n d9_BPK25_IP_TET1_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_BPK25_IP_TET1_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d9_BPK25_IP_TET1_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d9_BPK25_IP_TET1_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d9_BPK25_IP_TET1_rep3_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_1_for_TET1_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./d0_BPK25_IP_TET1_rep1_TvsC/d0_BPK25_IP_TET1_rep1_TvsC_peaks.narrowPeak
# 15105
# ./d0_BPK25_IP_TET1_rep2_TvsC/d0_BPK25_IP_TET1_rep2_TvsC_peaks.narrowPeak
# 3951
# ./d0_BPK25_IP_TET1_rep3_TvsC/d0_BPK25_IP_TET1_rep3_TvsC_peaks.narrowPeak
# 12116
# ./d0_DMSO_IP_TET1_rep1_TvsC/d0_DMSO_IP_TET1_rep1_TvsC_peaks.narrowPeak
# 1593
# ./d0_DMSO_IP_TET1_rep2_TvsC/d0_DMSO_IP_TET1_rep2_TvsC_peaks.narrowPeak
# 1432
# ./d0_DMSO_IP_TET1_rep3_TvsC/d0_DMSO_IP_TET1_rep3_TvsC_peaks.narrowPeak
# 2418
# ./d4_BPK25_IP_TET1_rep1_TvsC/d4_BPK25_IP_TET1_rep1_TvsC_peaks.narrowPeak
# 13106
# ./d4_BPK25_IP_TET1_rep2_TvsC/d4_BPK25_IP_TET1_rep2_TvsC_peaks.narrowPeak
# 9035
# ./d4_BPK25_IP_TET1_rep3_TvsC/d4_BPK25_IP_TET1_rep3_TvsC_peaks.narrowPeak
# 19716
# ./d4_DMSO_IP_TET1_rep1_TvsC/d4_DMSO_IP_TET1_rep1_TvsC_peaks.narrowPeak
# 3666
# ./d4_DMSO_IP_TET1_rep2_TvsC/d4_DMSO_IP_TET1_rep2_TvsC_peaks.narrowPeak
# 3765
# ./d4_DMSO_IP_TET1_rep3_TvsC/d4_DMSO_IP_TET1_rep3_TvsC_peaks.narrowPeak
# 3719
# ./d9_BPK25_IP_TET1_rep1_TvsC/d9_BPK25_IP_TET1_rep1_TvsC_peaks.narrowPeak
# 3394
# ./d9_BPK25_IP_TET1_rep2_TvsC/d9_BPK25_IP_TET1_rep2_TvsC_peaks.narrowPeak
# 4561
# ./d9_BPK25_IP_TET1_rep3_TvsC/d9_BPK25_IP_TET1_rep3_TvsC_peaks.narrowPeak
# 6076
# ./d9_DMSO_IP_TET1_rep1_TvsC/d9_DMSO_IP_TET1_rep1_TvsC_peaks.narrowPeak
# 1894
# ./d9_DMSO_IP_TET1_rep2_TvsC/d9_DMSO_IP_TET1_rep2_TvsC_peaks.narrowPeak
# 1367
# ./d9_DMSO_IP_TET1_rep3_TvsC/d9_DMSO_IP_TET1_rep3_TvsC_peaks.narrowPeak
# 531

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./d0_BPK25_IP_TET1_rep1_TvsC/d0_BPK25_IP_TET1_rep1_TvsC_summits.bed
# 16110
# ./d0_BPK25_IP_TET1_rep2_TvsC/d0_BPK25_IP_TET1_rep2_TvsC_summits.bed
# 4290
# ./d0_BPK25_IP_TET1_rep3_TvsC/d0_BPK25_IP_TET1_rep3_TvsC_summits.bed
# 12788
# ./d0_DMSO_IP_TET1_rep1_TvsC/d0_DMSO_IP_TET1_rep1_TvsC_summits.bed
# 1780
# ./d0_DMSO_IP_TET1_rep2_TvsC/d0_DMSO_IP_TET1_rep2_TvsC_summits.bed
# 1600
# ./d0_DMSO_IP_TET1_rep3_TvsC/d0_DMSO_IP_TET1_rep3_TvsC_summits.bed
# 2680
# ./d4_BPK25_IP_TET1_rep1_TvsC/d4_BPK25_IP_TET1_rep1_TvsC_summits.bed
# 13625
# ./d4_BPK25_IP_TET1_rep2_TvsC/d4_BPK25_IP_TET1_rep2_TvsC_summits.bed
# 9393
# ./d4_BPK25_IP_TET1_rep3_TvsC/d4_BPK25_IP_TET1_rep3_TvsC_summits.bed
# 20320
# ./d4_DMSO_IP_TET1_rep1_TvsC/d4_DMSO_IP_TET1_rep1_TvsC_summits.bed
# 3732
# ./d4_DMSO_IP_TET1_rep2_TvsC/d4_DMSO_IP_TET1_rep2_TvsC_summits.bed
# 3856
# ./d4_DMSO_IP_TET1_rep3_TvsC/d4_DMSO_IP_TET1_rep3_TvsC_summits.bed
# 3796
# ./d9_BPK25_IP_TET1_rep1_TvsC/d9_BPK25_IP_TET1_rep1_TvsC_summits.bed
# 3575
# ./d9_BPK25_IP_TET1_rep2_TvsC/d9_BPK25_IP_TET1_rep2_TvsC_summits.bed
# 4748
# ./d9_BPK25_IP_TET1_rep3_TvsC/d9_BPK25_IP_TET1_rep3_TvsC_summits.bed
# 6356
# ./d9_DMSO_IP_TET1_rep1_TvsC/d9_DMSO_IP_TET1_rep1_TvsC_summits.bed
# 1955
# ./d9_DMSO_IP_TET1_rep2_TvsC/d9_DMSO_IP_TET1_rep2_TvsC_summits.bed
# 1422
# ./d9_DMSO_IP_TET1_rep3_TvsC/d9_DMSO_IP_TET1_rep3_TvsC_summits.bed
# 563

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_BPK25_DMSO_IP_TET1_ChIP_Seq_samples_macs2_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done

```



### d049_DMSO_BPK25_TET2_ChIP_Seq

```bash
# data preperation
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir macs2_peak_calling_1_for_TET2_each_repeat
cd macs2_peak_calling_1_for_TET2_each_repeat
# 上传 bam 连接构造脚本并运行
bash 6_make_H9_d0_DMSO_BPK25_TET2_ChIP_Seq_sorted_bam_file_links.sh
bash 6_make_H9_d4_DMSO_BPK25_TET2_ChIP_Seq_sorted_bam_file_links.sh
bash 6_make_H9_d9_DMSO_BPK25_TET2_ChIP_Seq_sorted_bam_file_links.sh
ls -1

# d0_DMSO
nohup macs2 callpeak -t hESC_H9_d0_DMSO_IP_TET2_rep1.bowtie2.sorted.bam -c hESC_H9_d0_DMSO_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_DMSO_IP_TET2_rep1_TvsC -n d0_DMSO_IP_TET2_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_DMSO_IP_TET2_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d0_DMSO_IP_TET2_rep2.bowtie2.sorted.bam -c hESC_H9_d0_DMSO_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_DMSO_IP_TET2_rep2_TvsC -n d0_DMSO_IP_TET2_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_DMSO_IP_TET2_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d0_DMSO_IP_TET2_rep3.bowtie2.sorted.bam -c hESC_H9_d0_DMSO_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_DMSO_IP_TET2_rep3_TvsC -n d0_DMSO_IP_TET2_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_DMSO_IP_TET2_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d0_DMSO_IP_TET2_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d0_DMSO_IP_TET2_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d0_DMSO_IP_TET2_rep3_TvsC.runlog
# d0_BPK25
nohup macs2 callpeak -t hESC_H9_d0_BPK25_IP_TET2_rep1.bowtie2.sorted.bam -c hESC_H9_d0_BPK25_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_BPK25_IP_TET2_rep1_TvsC -n d0_BPK25_IP_TET2_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_BPK25_IP_TET2_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d0_BPK25_IP_TET2_rep2.bowtie2.sorted.bam -c hESC_H9_d0_BPK25_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_BPK25_IP_TET2_rep2_TvsC -n d0_BPK25_IP_TET2_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_BPK25_IP_TET2_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d0_BPK25_IP_TET2_rep3.bowtie2.sorted.bam -c hESC_H9_d0_BPK25_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_BPK25_IP_TET2_rep3_TvsC -n d0_BPK25_IP_TET2_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_BPK25_IP_TET2_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d0_BPK25_IP_TET2_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d0_BPK25_IP_TET2_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d0_BPK25_IP_TET2_rep3_TvsC.runlog
# d4_DMSO
nohup macs2 callpeak -t hESC_H9_d4_DMSO_IP_TET2_rep1.bowtie2.sorted.bam -c hESC_H9_d4_DMSO_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_DMSO_IP_TET2_rep1_TvsC -n d4_DMSO_IP_TET2_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_DMSO_IP_TET2_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d4_DMSO_IP_TET2_rep2.bowtie2.sorted.bam -c hESC_H9_d4_DMSO_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_DMSO_IP_TET2_rep2_TvsC -n d4_DMSO_IP_TET2_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_DMSO_IP_TET2_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d4_DMSO_IP_TET2_rep3.bowtie2.sorted.bam -c hESC_H9_d4_DMSO_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_DMSO_IP_TET2_rep3_TvsC -n d4_DMSO_IP_TET2_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_DMSO_IP_TET2_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d4_DMSO_IP_TET2_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d4_DMSO_IP_TET2_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d4_DMSO_IP_TET2_rep3_TvsC.runlog
# d4_BPK25
nohup macs2 callpeak -t hESC_H9_d4_BPK25_IP_TET2_rep1.bowtie2.sorted.bam -c hESC_H9_d4_BPK25_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_BPK25_IP_TET2_rep1_TvsC -n d4_BPK25_IP_TET2_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_BPK25_IP_TET2_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d4_BPK25_IP_TET2_rep2.bowtie2.sorted.bam -c hESC_H9_d4_BPK25_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_BPK25_IP_TET2_rep2_TvsC -n d4_BPK25_IP_TET2_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_BPK25_IP_TET2_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d4_BPK25_IP_TET2_rep3.bowtie2.sorted.bam -c hESC_H9_d4_BPK25_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_BPK25_IP_TET2_rep3_TvsC -n d4_BPK25_IP_TET2_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_BPK25_IP_TET2_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d4_BPK25_IP_TET2_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d4_BPK25_IP_TET2_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d4_BPK25_IP_TET2_rep3_TvsC.runlog
# d9_DMSO
nohup macs2 callpeak -t hESC_H9_d9_DMSO_IP_TET2_rep1.bowtie2.sorted.bam -c hESC_H9_d9_DMSO_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_DMSO_IP_TET2_rep1_TvsC -n d9_DMSO_IP_TET2_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_DMSO_IP_TET2_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d9_DMSO_IP_TET2_rep2.bowtie2.sorted.bam -c hESC_H9_d9_DMSO_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_DMSO_IP_TET2_rep2_TvsC -n d9_DMSO_IP_TET2_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_DMSO_IP_TET2_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d9_DMSO_IP_TET2_rep3.bowtie2.sorted.bam -c hESC_H9_d9_DMSO_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_DMSO_IP_TET2_rep3_TvsC -n d9_DMSO_IP_TET2_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_DMSO_IP_TET2_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d9_DMSO_IP_TET2_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d9_DMSO_IP_TET2_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d9_DMSO_IP_TET2_rep3_TvsC.runlog
# d9_BPK25
nohup macs2 callpeak -t hESC_H9_d9_BPK25_IP_TET2_rep1.bowtie2.sorted.bam -c hESC_H9_d9_BPK25_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_BPK25_IP_TET2_rep1_TvsC -n d9_BPK25_IP_TET2_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_BPK25_IP_TET2_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d9_BPK25_IP_TET2_rep2.bowtie2.sorted.bam -c hESC_H9_d9_BPK25_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_BPK25_IP_TET2_rep2_TvsC -n d9_BPK25_IP_TET2_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_BPK25_IP_TET2_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d9_BPK25_IP_TET2_rep3.bowtie2.sorted.bam -c hESC_H9_d9_BPK25_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_BPK25_IP_TET2_rep3_TvsC -n d9_BPK25_IP_TET2_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_BPK25_IP_TET2_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d9_BPK25_IP_TET2_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d9_BPK25_IP_TET2_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d9_BPK25_IP_TET2_rep3_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_1_for_TET2_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./d0_BPK25_IP_TET2_rep1_TvsC/d0_BPK25_IP_TET2_rep1_TvsC_peaks.narrowPeak
# 2039
# ./d0_BPK25_IP_TET2_rep2_TvsC/d0_BPK25_IP_TET2_rep2_TvsC_peaks.narrowPeak
# 30507
# ./d0_BPK25_IP_TET2_rep3_TvsC/d0_BPK25_IP_TET2_rep3_TvsC_peaks.narrowPeak
# 5975
# ./d0_DMSO_IP_TET2_rep1_TvsC/d0_DMSO_IP_TET2_rep1_TvsC_peaks.narrowPeak
# 910
# ./d0_DMSO_IP_TET2_rep2_TvsC/d0_DMSO_IP_TET2_rep2_TvsC_peaks.narrowPeak
# 448
# ./d0_DMSO_IP_TET2_rep3_TvsC/d0_DMSO_IP_TET2_rep3_TvsC_peaks.narrowPeak
# 440
# ./d4_BPK25_IP_TET2_rep1_TvsC/d4_BPK25_IP_TET2_rep1_TvsC_peaks.narrowPeak
# 4198
# ./d4_BPK25_IP_TET2_rep2_TvsC/d4_BPK25_IP_TET2_rep2_TvsC_peaks.narrowPeak
# 965
# ./d4_BPK25_IP_TET2_rep3_TvsC/d4_BPK25_IP_TET2_rep3_TvsC_peaks.narrowPeak
# 1075
# ./d4_DMSO_IP_TET2_rep1_TvsC/d4_DMSO_IP_TET2_rep1_TvsC_peaks.narrowPeak
# 3713
# ./d4_DMSO_IP_TET2_rep2_TvsC/d4_DMSO_IP_TET2_rep2_TvsC_peaks.narrowPeak
# 4430
# ./d4_DMSO_IP_TET2_rep3_TvsC/d4_DMSO_IP_TET2_rep3_TvsC_peaks.narrowPeak
# 3424
# ./d9_BPK25_IP_TET2_rep1_TvsC/d9_BPK25_IP_TET2_rep1_TvsC_peaks.narrowPeak
# 64
# ./d9_BPK25_IP_TET2_rep2_TvsC/d9_BPK25_IP_TET2_rep2_TvsC_peaks.narrowPeak
# 60
# ./d9_BPK25_IP_TET2_rep3_TvsC/d9_BPK25_IP_TET2_rep3_TvsC_peaks.narrowPeak
# 168
# ./d9_DMSO_IP_TET2_rep1_TvsC/d9_DMSO_IP_TET2_rep1_TvsC_peaks.narrowPeak
# 2199
# ./d9_DMSO_IP_TET2_rep2_TvsC/d9_DMSO_IP_TET2_rep2_TvsC_peaks.narrowPeak
# 1090
# ./d9_DMSO_IP_TET2_rep3_TvsC/d9_DMSO_IP_TET2_rep3_TvsC_peaks.narrowPeak
# 296

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./d0_BPK25_IP_TET2_rep1_TvsC/d0_BPK25_IP_TET2_rep1_TvsC_summits.bed
# 2277
# ./d0_BPK25_IP_TET2_rep2_TvsC/d0_BPK25_IP_TET2_rep2_TvsC_summits.bed
# 31514
# ./d0_BPK25_IP_TET2_rep3_TvsC/d0_BPK25_IP_TET2_rep3_TvsC_summits.bed
# 6332
# ./d0_DMSO_IP_TET2_rep1_TvsC/d0_DMSO_IP_TET2_rep1_TvsC_summits.bed
# 1028
# ./d0_DMSO_IP_TET2_rep2_TvsC/d0_DMSO_IP_TET2_rep2_TvsC_summits.bed
# 547
# ./d0_DMSO_IP_TET2_rep3_TvsC/d0_DMSO_IP_TET2_rep3_TvsC_summits.bed
# 541
# ./d4_BPK25_IP_TET2_rep1_TvsC/d4_BPK25_IP_TET2_rep1_TvsC_summits.bed
# 4566
# ./d4_BPK25_IP_TET2_rep2_TvsC/d4_BPK25_IP_TET2_rep2_TvsC_summits.bed
# 1054
# ./d4_BPK25_IP_TET2_rep3_TvsC/d4_BPK25_IP_TET2_rep3_TvsC_summits.bed
# 1169
# ./d4_DMSO_IP_TET2_rep1_TvsC/d4_DMSO_IP_TET2_rep1_TvsC_summits.bed
# 3787
# ./d4_DMSO_IP_TET2_rep2_TvsC/d4_DMSO_IP_TET2_rep2_TvsC_summits.bed
# 4496
# ./d4_DMSO_IP_TET2_rep3_TvsC/d4_DMSO_IP_TET2_rep3_TvsC_summits.bed
# 3481
# ./d9_BPK25_IP_TET2_rep1_TvsC/d9_BPK25_IP_TET2_rep1_TvsC_summits.bed
# 78
# ./d9_BPK25_IP_TET2_rep2_TvsC/d9_BPK25_IP_TET2_rep2_TvsC_summits.bed
# 76
# ./d9_BPK25_IP_TET2_rep3_TvsC/d9_BPK25_IP_TET2_rep3_TvsC_summits.bed
# 182
# ./d9_DMSO_IP_TET2_rep1_TvsC/d9_DMSO_IP_TET2_rep1_TvsC_summits.bed
# 2234
# ./d9_DMSO_IP_TET2_rep2_TvsC/d9_DMSO_IP_TET2_rep2_TvsC_summits.bed
# 1125
# ./d9_DMSO_IP_TET2_rep3_TvsC/d9_DMSO_IP_TET2_rep3_TvsC_summits.bed
# 316

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_BPK25_DMSO_IP_TET2_ChIP_Seq_samples_macs2_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done

```



### d049_DMSO_BPK25_TET3_ChIP_Seq

```bash
# data preperation
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir macs2_peak_calling_1_for_TET3_each_repeat
cd macs2_peak_calling_1_for_TET3_each_repeat
# 上传 bam 连接构造脚本并运行
bash 6_make_H9_d0_DMSO_BPK25_TET3_ChIP_Seq_sorted_bam_file_links.sh
bash 6_make_H9_d4_DMSO_BPK25_TET3_ChIP_Seq_sorted_bam_file_links.sh
bash 6_make_H9_d9_DMSO_BPK25_TET3_ChIP_Seq_sorted_bam_file_links.sh
ls -1

# d0_DMSO
nohup macs2 callpeak -t hESC_H9_d0_DMSO_IP_TET3_rep1.bowtie2.sorted.bam -c hESC_H9_d0_DMSO_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_DMSO_IP_TET3_rep1_TvsC -n d0_DMSO_IP_TET3_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_DMSO_IP_TET3_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d0_DMSO_IP_TET3_rep2.bowtie2.sorted.bam -c hESC_H9_d0_DMSO_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_DMSO_IP_TET3_rep2_TvsC -n d0_DMSO_IP_TET3_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_DMSO_IP_TET3_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d0_DMSO_IP_TET3_rep3.bowtie2.sorted.bam -c hESC_H9_d0_DMSO_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_DMSO_IP_TET3_rep3_TvsC -n d0_DMSO_IP_TET3_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_DMSO_IP_TET3_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d0_DMSO_IP_TET3_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d0_DMSO_IP_TET3_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d0_DMSO_IP_TET3_rep3_TvsC.runlog
# d0_BPK25
nohup macs2 callpeak -t hESC_H9_d0_BPK25_IP_TET3_rep1.bowtie2.sorted.bam -c hESC_H9_d0_BPK25_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_BPK25_IP_TET3_rep1_TvsC -n d0_BPK25_IP_TET3_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_BPK25_IP_TET3_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d0_BPK25_IP_TET3_rep2.bowtie2.sorted.bam -c hESC_H9_d0_BPK25_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_BPK25_IP_TET3_rep2_TvsC -n d0_BPK25_IP_TET3_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_BPK25_IP_TET3_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d0_BPK25_IP_TET3_rep3.bowtie2.sorted.bam -c hESC_H9_d0_BPK25_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_BPK25_IP_TET3_rep3_TvsC -n d0_BPK25_IP_TET3_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d0_BPK25_IP_TET3_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d0_BPK25_IP_TET3_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d0_BPK25_IP_TET3_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d0_BPK25_IP_TET3_rep3_TvsC.runlog
# d4_DMSO
nohup macs2 callpeak -t hESC_H9_d4_DMSO_IP_TET3_rep1.bowtie2.sorted.bam -c hESC_H9_d4_DMSO_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_DMSO_IP_TET3_rep1_TvsC -n d4_DMSO_IP_TET3_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_DMSO_IP_TET3_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d4_DMSO_IP_TET3_rep2.bowtie2.sorted.bam -c hESC_H9_d4_DMSO_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_DMSO_IP_TET3_rep2_TvsC -n d4_DMSO_IP_TET3_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_DMSO_IP_TET3_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d4_DMSO_IP_TET3_rep3.bowtie2.sorted.bam -c hESC_H9_d4_DMSO_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_DMSO_IP_TET3_rep3_TvsC -n d4_DMSO_IP_TET3_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_DMSO_IP_TET3_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d4_DMSO_IP_TET3_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d4_DMSO_IP_TET3_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d4_DMSO_IP_TET3_rep3_TvsC.runlog
# d4_BPK25
nohup macs2 callpeak -t hESC_H9_d4_BPK25_IP_TET3_rep1.bowtie2.sorted.bam -c hESC_H9_d4_BPK25_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_BPK25_IP_TET3_rep1_TvsC -n d4_BPK25_IP_TET3_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_BPK25_IP_TET3_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d4_BPK25_IP_TET3_rep2.bowtie2.sorted.bam -c hESC_H9_d4_BPK25_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_BPK25_IP_TET3_rep2_TvsC -n d4_BPK25_IP_TET3_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_BPK25_IP_TET3_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d4_BPK25_IP_TET3_rep3.bowtie2.sorted.bam -c hESC_H9_d4_BPK25_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_BPK25_IP_TET3_rep3_TvsC -n d4_BPK25_IP_TET3_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d4_BPK25_IP_TET3_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d4_BPK25_IP_TET3_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d4_BPK25_IP_TET3_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d4_BPK25_IP_TET3_rep3_TvsC.runlog
# d9_DMSO
nohup macs2 callpeak -t hESC_H9_d9_DMSO_IP_TET3_rep1.bowtie2.sorted.bam -c hESC_H9_d9_DMSO_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_DMSO_IP_TET3_rep1_TvsC -n d9_DMSO_IP_TET3_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_DMSO_IP_TET3_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d9_DMSO_IP_TET3_rep2.bowtie2.sorted.bam -c hESC_H9_d9_DMSO_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_DMSO_IP_TET3_rep2_TvsC -n d9_DMSO_IP_TET3_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_DMSO_IP_TET3_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d9_DMSO_IP_TET3_rep3.bowtie2.sorted.bam -c hESC_H9_d9_DMSO_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_DMSO_IP_TET3_rep3_TvsC -n d9_DMSO_IP_TET3_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_DMSO_IP_TET3_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d9_DMSO_IP_TET3_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d9_DMSO_IP_TET3_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d9_DMSO_IP_TET3_rep3_TvsC.runlog
# d9_BPK25
nohup macs2 callpeak -t hESC_H9_d9_BPK25_IP_TET3_rep1.bowtie2.sorted.bam -c hESC_H9_d9_BPK25_input_rep1.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_BPK25_IP_TET3_rep1_TvsC -n d9_BPK25_IP_TET3_rep1_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_BPK25_IP_TET3_rep1_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d9_BPK25_IP_TET3_rep2.bowtie2.sorted.bam -c hESC_H9_d9_BPK25_input_rep2.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_BPK25_IP_TET3_rep2_TvsC -n d9_BPK25_IP_TET3_rep2_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_BPK25_IP_TET3_rep2_TvsC.runlog 2>&1 &
nohup macs2 callpeak -t hESC_H9_d9_BPK25_IP_TET3_rep3.bowtie2.sorted.bam -c hESC_H9_d9_BPK25_input_rep3.bowtie2.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_BPK25_IP_TET3_rep3_TvsC -n d9_BPK25_IP_TET3_rep3_TvsC -B --SPMR --call-summits --cutoff-analysis >macs2_callpeak_for_d9_BPK25_IP_TET3_rep3_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d9_BPK25_IP_TET3_rep1_TvsC.runlog
tail -f macs2_callpeak_for_d9_BPK25_IP_TET3_rep2_TvsC.runlog
tail -f macs2_callpeak_for_d9_BPK25_IP_TET3_rep3_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_1_for_TET3_each_repeat
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./d0_BPK25_IP_TET3_rep1_TvsC/d0_BPK25_IP_TET3_rep1_TvsC_peaks.narrowPeak
# 3390
# ./d0_BPK25_IP_TET3_rep2_TvsC/d0_BPK25_IP_TET3_rep2_TvsC_peaks.narrowPeak
# 1861
# ./d0_BPK25_IP_TET3_rep3_TvsC/d0_BPK25_IP_TET3_rep3_TvsC_peaks.narrowPeak
# 1035
# ./d0_DMSO_IP_TET3_rep1_TvsC/d0_DMSO_IP_TET3_rep1_TvsC_peaks.narrowPeak
# 9513
# ./d0_DMSO_IP_TET3_rep2_TvsC/d0_DMSO_IP_TET3_rep2_TvsC_peaks.narrowPeak
# 7307
# ./d0_DMSO_IP_TET3_rep3_TvsC/d0_DMSO_IP_TET3_rep3_TvsC_peaks.narrowPeak
# 11088
# ./d4_BPK25_IP_TET3_rep1_TvsC/d4_BPK25_IP_TET3_rep1_TvsC_peaks.narrowPeak
# 493
# ./d4_BPK25_IP_TET3_rep2_TvsC/d4_BPK25_IP_TET3_rep2_TvsC_peaks.narrowPeak
# 620
# ./d4_BPK25_IP_TET3_rep3_TvsC/d4_BPK25_IP_TET3_rep3_TvsC_peaks.narrowPeak
# 343
# ./d4_DMSO_IP_TET3_rep1_TvsC/d4_DMSO_IP_TET3_rep1_TvsC_peaks.narrowPeak
# 2616
# ./d4_DMSO_IP_TET3_rep2_TvsC/d4_DMSO_IP_TET3_rep2_TvsC_peaks.narrowPeak
# 2755
# ./d4_DMSO_IP_TET3_rep3_TvsC/d4_DMSO_IP_TET3_rep3_TvsC_peaks.narrowPeak
# 5207
# ./d9_BPK25_IP_TET3_rep1_TvsC/d9_BPK25_IP_TET3_rep1_TvsC_peaks.narrowPeak
# 36
# ./d9_BPK25_IP_TET3_rep2_TvsC/d9_BPK25_IP_TET3_rep2_TvsC_peaks.narrowPeak
# 29
# ./d9_BPK25_IP_TET3_rep3_TvsC/d9_BPK25_IP_TET3_rep3_TvsC_peaks.narrowPeak
# 36
# ./d9_DMSO_IP_TET3_rep1_TvsC/d9_DMSO_IP_TET3_rep1_TvsC_peaks.narrowPeak
# 2085
# ./d9_DMSO_IP_TET3_rep2_TvsC/d9_DMSO_IP_TET3_rep2_TvsC_peaks.narrowPeak
# 2118
# ./d9_DMSO_IP_TET3_rep3_TvsC/d9_DMSO_IP_TET3_rep3_TvsC_peaks.narrowPeak
# 1010

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./d0_BPK25_IP_TET3_rep1_TvsC/d0_BPK25_IP_TET3_rep1_TvsC_summits.bed
# 3639
# ./d0_BPK25_IP_TET3_rep2_TvsC/d0_BPK25_IP_TET3_rep2_TvsC_summits.bed
# 2088
# ./d0_BPK25_IP_TET3_rep3_TvsC/d0_BPK25_IP_TET3_rep3_TvsC_summits.bed
# 1185
# ./d0_DMSO_IP_TET3_rep1_TvsC/d0_DMSO_IP_TET3_rep1_TvsC_summits.bed
# 10293
# ./d0_DMSO_IP_TET3_rep2_TvsC/d0_DMSO_IP_TET3_rep2_TvsC_summits.bed
# 7885
# ./d0_DMSO_IP_TET3_rep3_TvsC/d0_DMSO_IP_TET3_rep3_TvsC_summits.bed
# 11857
# ./d4_BPK25_IP_TET3_rep1_TvsC/d4_BPK25_IP_TET3_rep1_TvsC_summits.bed
# 543
# ./d4_BPK25_IP_TET3_rep2_TvsC/d4_BPK25_IP_TET3_rep2_TvsC_summits.bed
# 668
# ./d4_BPK25_IP_TET3_rep3_TvsC/d4_BPK25_IP_TET3_rep3_TvsC_summits.bed
# 392
# ./d4_DMSO_IP_TET3_rep1_TvsC/d4_DMSO_IP_TET3_rep1_TvsC_summits.bed
# 2684
# ./d4_DMSO_IP_TET3_rep2_TvsC/d4_DMSO_IP_TET3_rep2_TvsC_summits.bed
# 2838
# ./d4_DMSO_IP_TET3_rep3_TvsC/d4_DMSO_IP_TET3_rep3_TvsC_summits.bed
# 5302
# ./d9_BPK25_IP_TET3_rep1_TvsC/d9_BPK25_IP_TET3_rep1_TvsC_summits.bed
# 50
# ./d9_BPK25_IP_TET3_rep2_TvsC/d9_BPK25_IP_TET3_rep2_TvsC_summits.bed
# 50
# ./d9_BPK25_IP_TET3_rep3_TvsC/d9_BPK25_IP_TET3_rep3_TvsC_summits.bed
# 51
# ./d9_DMSO_IP_TET3_rep1_TvsC/d9_DMSO_IP_TET3_rep1_TvsC_summits.bed
# 2147
# ./d9_DMSO_IP_TET3_rep2_TvsC/d9_DMSO_IP_TET3_rep2_TvsC_summits.bed
# 2172
# ./d9_DMSO_IP_TET3_rep3_TvsC/d9_DMSO_IP_TET3_rep3_TvsC_summits.bed
# 1053

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_BPK25_DMSO_IP_TET3_ChIP_Seq_samples_macs2_results_summary `find . -name *peaks.xls`

# 将bdg转为bigwig用于可视化
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done

```



## peak calling for rep123 merged

### d049_DMSO_BPK25_CTCF_ChIP_Seq

```bash
# data preperation
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir macs2_peak_calling_2_for_CTCF_rep123_together
cd macs2_peak_calling_2_for_CTCF_rep123_together
# 上传 bam 连接构造脚本并运行
bash 7_make_H9_d049_DMSO_BPK25_CTCF_ChIP_Seq_rep123_merged_bam_file_links.sh
ls -1

# d0_DMSO_IP_CTCF_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d0_DMSO_IP_CTCF_rep123.merged.sorted.bam -c hESC_H9_d0_DMSO_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_DMSO_IP_CTCF_rep123_TvsC -n d0_DMSO_IP_CTCF_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d0_DMSO_IP_CTCF_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d0_DMSO_IP_CTCF_rep123_TvsC.runlog
# d0_BPK25_IP_CTCF_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d0_BPK25_IP_CTCF_rep123.merged.sorted.bam -c hESC_H9_d0_BPK25_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_BPK25_IP_CTCF_rep123_TvsC -n d0_BPK25_IP_CTCF_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d0_BPK25_IP_CTCF_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d0_BPK25_IP_CTCF_rep123_TvsC.runlog
# d4_DMSO_IP_CTCF_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d4_DMSO_IP_CTCF_rep123.merged.sorted.bam -c hESC_H9_d4_DMSO_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_DMSO_IP_CTCF_rep123_TvsC -n d4_DMSO_IP_CTCF_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d4_DMSO_IP_CTCF_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d4_DMSO_IP_CTCF_rep123_TvsC.runlog
# d4_BPK25_IP_CTCF_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d4_BPK25_IP_CTCF_rep123.merged.sorted.bam -c hESC_H9_d4_BPK25_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_BPK25_IP_CTCF_rep123_TvsC -n d4_BPK25_IP_CTCF_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d4_BPK25_IP_CTCF_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d4_BPK25_IP_CTCF_rep123_TvsC.runlog
# d9_DMSO_IP_CTCF_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d9_DMSO_IP_CTCF_rep123.merged.sorted.bam -c hESC_H9_d9_DMSO_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_DMSO_IP_CTCF_rep123_TvsC -n d9_DMSO_IP_CTCF_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d9_DMSO_IP_CTCF_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d9_DMSO_IP_CTCF_rep123_TvsC.runlog
# d9_BPK25_IP_CTCF_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d9_BPK25_IP_CTCF_rep123.merged.sorted.bam -c hESC_H9_d9_BPK25_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_BPK25_IP_CTCF_rep123_TvsC -n d9_BPK25_IP_CTCF_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d9_BPK25_IP_CTCF_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d9_BPK25_IP_CTCF_rep123_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_CTCF_rep123_together
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./d0_BPK25_IP_CTCF_rep123_TvsC/d0_BPK25_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# 1296
# ./d0_DMSO_IP_CTCF_rep123_TvsC/d0_DMSO_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# 1843
# ./d4_BPK25_IP_CTCF_rep123_TvsC/d4_BPK25_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# 14500
# ./d4_DMSO_IP_CTCF_rep123_TvsC/d4_DMSO_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# 20675
# ./d9_BPK25_IP_CTCF_rep123_TvsC/d9_BPK25_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# 3795
# ./d9_DMSO_IP_CTCF_rep123_TvsC/d9_DMSO_IP_CTCF_rep123_TvsC_peaks.narrowPeak
# 8469

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./d0_BPK25_IP_CTCF_rep123_TvsC/d0_BPK25_IP_CTCF_rep123_TvsC_summits.bed
# 1551
# ./d0_DMSO_IP_CTCF_rep123_TvsC/d0_DMSO_IP_CTCF_rep123_TvsC_summits.bed
# 2329
# ./d4_BPK25_IP_CTCF_rep123_TvsC/d4_BPK25_IP_CTCF_rep123_TvsC_summits.bed
# 17812
# ./d4_DMSO_IP_CTCF_rep123_TvsC/d4_DMSO_IP_CTCF_rep123_TvsC_summits.bed
# 23330
# ./d9_BPK25_IP_CTCF_rep123_TvsC/d9_BPK25_IP_CTCF_rep123_TvsC_summits.bed
# 4595
# ./d9_DMSO_IP_CTCF_rep123_TvsC/d9_DMSO_IP_CTCF_rep123_TvsC_summits.bed
# 8842

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_d049_BPK25_DMSO_CTCF_ChIP_Seq_macs2_results_summary `find . -name *peaks.xls`


cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_CTCF_rep123_together
# 将bdg转为bigwig用于可视化
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done

```



#### d0_DMSO_BPK25 peak intersect

```bash
# d0_DMSO_BPK25 两组 peak 交集
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_CTCF_rep123_together
for i in `find . -name *d0*narrowPeak |sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
bedtools intersect -a ./d0_DMSO_IP_CTCF_rep123_TvsC/d0_DMSO_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -b ./d0_BPK25_IP_CTCF_rep123_TvsC/d0_BPK25_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 385
bedtools intersect -a ./d0_BPK25_IP_CTCF_rep123_TvsC/d0_BPK25_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -b ./d0_DMSO_IP_CTCF_rep123_TvsC/d0_DMSO_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 398

# 两组 summit extB50 交集数量
bedtools slop -b 50 -i ./d0_DMSO_IP_CTCF_rep123_TvsC/d0_DMSO_IP_CTCF_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d0_DMSO_IP_CTCF_rep123_TvsC/d0_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
bedtools slop -b 50 -i ./d0_BPK25_IP_CTCF_rep123_TvsC/d0_BPK25_IP_CTCF_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d0_BPK25_IP_CTCF_rep123_TvsC/d0_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
# 求交集
bedtools intersect -a ./d0_DMSO_IP_CTCF_rep123_TvsC/d0_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_CTCF_rep123_TvsC/d0_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 452
bedtools intersect -a ./d0_BPK25_IP_CTCF_rep123_TvsC/d0_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_DMSO_IP_CTCF_rep123_TvsC/d0_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 452

# 构造 三类 summit
# DMSO only summit
bedtools intersect -a ./d0_DMSO_IP_CTCF_rep123_TvsC/d0_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_CTCF_rep123_TvsC/d0_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d0_DMSO_IP_CTCF_rep123_TvsC/d0_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_CTCF_rep123_TvsC/d0_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d0_IP_CTCF.DMSO_only.summit.bed
# common summit
bedtools intersect -a ./d0_DMSO_IP_CTCF_rep123_TvsC/d0_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_CTCF_rep123_TvsC/d0_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a ./d0_DMSO_IP_CTCF_rep123_TvsC/d0_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_CTCF_rep123_TvsC/d0_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >d0_IP_CTCF.DB_common.summit.bed
# BPK25 only summit
bedtools intersect -a ./d0_BPK25_IP_CTCF_rep123_TvsC/d0_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_DMSO_IP_CTCF_rep123_TvsC/d0_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d0_BPK25_IP_CTCF_rep123_TvsC/d0_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_DMSO_IP_CTCF_rep123_TvsC/d0_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d0_IP_CTCF.BPK25_only.summit.bed

wc -l d0_IP_CTCF*summit.bed
# 1099 d0_IP_CTCF.BPK25_only.summit.bed
#  455 d0_IP_CTCF.DB_common.summit.bed
# 1877 d0_IP_CTCF.DMSO_only.summit.bed

```

#### d4_DMSO_BPK25 peak intersect

```bash
# d4_DMSO_BPK25 两组 peak 交集
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_CTCF_rep123_together
for i in `find . -name *d4_*narrowPeak |sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
bedtools intersect -a ./d4_DMSO_IP_CTCF_rep123_TvsC/d4_DMSO_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -b ./d4_BPK25_IP_CTCF_rep123_TvsC/d4_BPK25_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 7932
bedtools intersect -a ./d4_BPK25_IP_CTCF_rep123_TvsC/d4_BPK25_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -b ./d4_DMSO_IP_CTCF_rep123_TvsC/d4_DMSO_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 7825

# 两组 summit extB50 交集数量
bedtools slop -b 50 -i ./d4_DMSO_IP_CTCF_rep123_TvsC/d4_DMSO_IP_CTCF_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d4_DMSO_IP_CTCF_rep123_TvsC/d4_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
bedtools slop -b 50 -i ./d4_BPK25_IP_CTCF_rep123_TvsC/d4_BPK25_IP_CTCF_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d4_BPK25_IP_CTCF_rep123_TvsC/d4_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
# 求交集
bedtools intersect -a ./d4_DMSO_IP_CTCF_rep123_TvsC/d4_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_CTCF_rep123_TvsC/d4_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 8575
bedtools intersect -a ./d4_BPK25_IP_CTCF_rep123_TvsC/d4_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_DMSO_IP_CTCF_rep123_TvsC/d4_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 8566

# 构造 三类 summit
# DMSO only summit
bedtools intersect -a ./d4_DMSO_IP_CTCF_rep123_TvsC/d4_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_CTCF_rep123_TvsC/d4_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d4_DMSO_IP_CTCF_rep123_TvsC/d4_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_CTCF_rep123_TvsC/d4_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d4_IP_CTCF.DMSO_only.summit.bed
# common summit
bedtools intersect -a ./d4_DMSO_IP_CTCF_rep123_TvsC/d4_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_CTCF_rep123_TvsC/d4_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a ./d4_DMSO_IP_CTCF_rep123_TvsC/d4_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_CTCF_rep123_TvsC/d4_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >d4_IP_CTCF.DB_common.summit.bed
# BPK25 only summit
bedtools intersect -a ./d4_BPK25_IP_CTCF_rep123_TvsC/d4_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_DMSO_IP_CTCF_rep123_TvsC/d4_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d4_BPK25_IP_CTCF_rep123_TvsC/d4_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_DMSO_IP_CTCF_rep123_TvsC/d4_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d4_IP_CTCF.BPK25_only.summit.bed

wc -l d4_IP_CTCF*summit.bed
#  9246 d4_IP_CTCF.BPK25_only.summit.bed
#  8618 d4_IP_CTCF.DB_common.summit.bed
# 14755 d4_IP_CTCF.DMSO_only.summit.bed

```

#### d9_DMSO_BPK25 peak intersect

```bash
# d9_DMSO_BPK25 两组 peak 交集
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_CTCF_rep123_together
for i in `find . -name *d9_*narrowPeak |sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
bedtools intersect -a ./d9_DMSO_IP_CTCF_rep123_TvsC/d9_DMSO_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -b ./d9_BPK25_IP_CTCF_rep123_TvsC/d9_BPK25_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 1025
bedtools intersect -a ./d9_BPK25_IP_CTCF_rep123_TvsC/d9_BPK25_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -b ./d9_DMSO_IP_CTCF_rep123_TvsC/d9_DMSO_IP_CTCF_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 1014

# 两组 summit extB50 交集数量
bedtools slop -b 50 -i ./d9_DMSO_IP_CTCF_rep123_TvsC/d9_DMSO_IP_CTCF_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d9_DMSO_IP_CTCF_rep123_TvsC/d9_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
bedtools slop -b 50 -i ./d9_BPK25_IP_CTCF_rep123_TvsC/d9_BPK25_IP_CTCF_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d9_BPK25_IP_CTCF_rep123_TvsC/d9_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed
# 求交集
bedtools intersect -a ./d9_DMSO_IP_CTCF_rep123_TvsC/d9_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_CTCF_rep123_TvsC/d9_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 994
bedtools intersect -a ./d9_BPK25_IP_CTCF_rep123_TvsC/d9_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_DMSO_IP_CTCF_rep123_TvsC/d9_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 992

# 构造 三类 summit
# DMSO only summit
bedtools intersect -a ./d9_DMSO_IP_CTCF_rep123_TvsC/d9_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_CTCF_rep123_TvsC/d9_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d9_DMSO_IP_CTCF_rep123_TvsC/d9_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_CTCF_rep123_TvsC/d9_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d9_IP_CTCF.DMSO_only.summit.bed
# common summit
bedtools intersect -a ./d9_DMSO_IP_CTCF_rep123_TvsC/d9_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_CTCF_rep123_TvsC/d9_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a ./d9_DMSO_IP_CTCF_rep123_TvsC/d9_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_CTCF_rep123_TvsC/d9_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >d9_IP_CTCF.DB_common.summit.bed
# BPK25 only summit
bedtools intersect -a ./d9_BPK25_IP_CTCF_rep123_TvsC/d9_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_DMSO_IP_CTCF_rep123_TvsC/d9_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d9_BPK25_IP_CTCF_rep123_TvsC/d9_BPK25_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_DMSO_IP_CTCF_rep123_TvsC/d9_DMSO_IP_CTCF_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d9_IP_CTCF.BPK25_only.summit.bed

wc -l d9_IP_CTCF*summit.bed
# 3603 d9_IP_CTCF.BPK25_only.summit.bed
#  999 d9_IP_CTCF.DB_common.summit.bed
# 7848 d9_IP_CTCF.DMSO_only.summit.bed

```



### d049_DMSO_BPK25_MBD3_ChIP_Seq

```bash
# data preperation
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir macs2_peak_calling_2_for_MBD3_rep123_together
cd macs2_peak_calling_2_for_MBD3_rep123_together
# 上传 bam 连接构造脚本并运行
bash 7_make_H9_d049_DMSO_BPK25_MBD3_ChIP_Seq_rep123_merged_bam_file_links.sh
ls -1

# d0_DMSO_IP_MBD3_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d0_DMSO_IP_MBD3_rep123.merged.sorted.bam -c hESC_H9_d0_DMSO_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_DMSO_IP_MBD3_rep123_TvsC -n d0_DMSO_IP_MBD3_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d0_DMSO_IP_MBD3_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d0_DMSO_IP_MBD3_rep123_TvsC.runlog
# d0_BPK25_IP_MBD3_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d0_BPK25_IP_MBD3_rep123.merged.sorted.bam -c hESC_H9_d0_BPK25_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_BPK25_IP_MBD3_rep123_TvsC -n d0_BPK25_IP_MBD3_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d0_BPK25_IP_MBD3_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d0_BPK25_IP_MBD3_rep123_TvsC.runlog
# d4_DMSO_IP_MBD3_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d4_DMSO_IP_MBD3_rep123.merged.sorted.bam -c hESC_H9_d4_DMSO_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_DMSO_IP_MBD3_rep123_TvsC -n d4_DMSO_IP_MBD3_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d4_DMSO_IP_MBD3_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d4_DMSO_IP_MBD3_rep123_TvsC.runlog
# d4_BPK25_IP_MBD3_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d4_BPK25_IP_MBD3_rep123.merged.sorted.bam -c hESC_H9_d4_BPK25_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_BPK25_IP_MBD3_rep123_TvsC -n d4_BPK25_IP_MBD3_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d4_BPK25_IP_MBD3_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d4_BPK25_IP_MBD3_rep123_TvsC.runlog
# d9_DMSO_IP_MBD3_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d9_DMSO_IP_MBD3_rep123.merged.sorted.bam -c hESC_H9_d9_DMSO_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_DMSO_IP_MBD3_rep123_TvsC -n d9_DMSO_IP_MBD3_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d9_DMSO_IP_MBD3_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d9_DMSO_IP_MBD3_rep123_TvsC.runlog
# d9_BPK25_IP_MBD3_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d9_BPK25_IP_MBD3_rep123.merged.sorted.bam -c hESC_H9_d9_BPK25_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_BPK25_IP_MBD3_rep123_TvsC -n d9_BPK25_IP_MBD3_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d9_BPK25_IP_MBD3_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d9_BPK25_IP_MBD3_rep123_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_MBD3_rep123_together
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./d0_BPK25_IP_MBD3_rep123_TvsC/d0_BPK25_IP_MBD3_rep123_TvsC_peaks.narrowPeak
# 590
# ./d0_DMSO_IP_MBD3_rep123_TvsC/d0_DMSO_IP_MBD3_rep123_TvsC_peaks.narrowPeak
# 7511
# ./d4_BPK25_IP_MBD3_rep123_TvsC/d4_BPK25_IP_MBD3_rep123_TvsC_peaks.narrowPeak
# 15456
# ./d4_DMSO_IP_MBD3_rep123_TvsC/d4_DMSO_IP_MBD3_rep123_TvsC_peaks.narrowPeak
# 46132
# ./d9_BPK25_IP_MBD3_rep123_TvsC/d9_BPK25_IP_MBD3_rep123_TvsC_peaks.narrowPeak
# 44297
# ./d9_DMSO_IP_MBD3_rep123_TvsC/d9_DMSO_IP_MBD3_rep123_TvsC_peaks.narrowPeak
# 44052

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./d0_BPK25_IP_MBD3_rep123_TvsC/d0_BPK25_IP_MBD3_rep123_TvsC_summits.bed
# 738
# ./d0_DMSO_IP_MBD3_rep123_TvsC/d0_DMSO_IP_MBD3_rep123_TvsC_summits.bed
# 8539
# ./d4_BPK25_IP_MBD3_rep123_TvsC/d4_BPK25_IP_MBD3_rep123_TvsC_summits.bed
# 17954
# ./d4_DMSO_IP_MBD3_rep123_TvsC/d4_DMSO_IP_MBD3_rep123_TvsC_summits.bed
# 49283
# ./d9_BPK25_IP_MBD3_rep123_TvsC/d9_BPK25_IP_MBD3_rep123_TvsC_summits.bed
# 48184
# ./d9_DMSO_IP_MBD3_rep123_TvsC/d9_DMSO_IP_MBD3_rep123_TvsC_summits.bed
# 46422

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_d049_BPK25_DMSO_MBD3_ChIP_Seq_macs2_results_summary `find . -name *peaks.xls`


cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_MBD3_rep123_together
# 将bdg转为bigwig用于可视化
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done

```



#### d0_DMSO_BPK25 peak intersect

```bash
# d0_DMSO_BPK25 两组 peak 交集
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_MBD3_rep123_together
for i in `find . -name *d0*narrowPeak |sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
bedtools intersect -a ./d0_DMSO_IP_MBD3_rep123_TvsC/d0_DMSO_IP_MBD3_rep123_TvsC_peaks.narrowPeak.bed -b ./d0_BPK25_IP_MBD3_rep123_TvsC/d0_BPK25_IP_MBD3_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 443
bedtools intersect -a ./d0_BPK25_IP_MBD3_rep123_TvsC/d0_BPK25_IP_MBD3_rep123_TvsC_peaks.narrowPeak.bed -b ./d0_DMSO_IP_MBD3_rep123_TvsC/d0_DMSO_IP_MBD3_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 432

# 两组 summit extB50 交集数量
bedtools slop -b 50 -i ./d0_DMSO_IP_MBD3_rep123_TvsC/d0_DMSO_IP_MBD3_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d0_DMSO_IP_MBD3_rep123_TvsC/d0_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed
bedtools slop -b 50 -i ./d0_BPK25_IP_MBD3_rep123_TvsC/d0_BPK25_IP_MBD3_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d0_BPK25_IP_MBD3_rep123_TvsC/d0_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed
# 求交集
bedtools intersect -a ./d0_DMSO_IP_MBD3_rep123_TvsC/d0_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_MBD3_rep123_TvsC/d0_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 527
bedtools intersect -a ./d0_BPK25_IP_MBD3_rep123_TvsC/d0_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_DMSO_IP_MBD3_rep123_TvsC/d0_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 523

# 构造 三类 summit
# DMSO only summit
bedtools intersect -a ./d0_DMSO_IP_MBD3_rep123_TvsC/d0_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_MBD3_rep123_TvsC/d0_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d0_DMSO_IP_MBD3_rep123_TvsC/d0_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_MBD3_rep123_TvsC/d0_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d0_IP_MBD3.DMSO_only.summit.bed
# common summit
bedtools intersect -a ./d0_DMSO_IP_MBD3_rep123_TvsC/d0_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_MBD3_rep123_TvsC/d0_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a ./d0_DMSO_IP_MBD3_rep123_TvsC/d0_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_MBD3_rep123_TvsC/d0_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >d0_IP_MBD3.DB_common.summit.bed
# BPK25 only summit
bedtools intersect -a ./d0_BPK25_IP_MBD3_rep123_TvsC/d0_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_DMSO_IP_MBD3_rep123_TvsC/d0_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d0_BPK25_IP_MBD3_rep123_TvsC/d0_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_DMSO_IP_MBD3_rep123_TvsC/d0_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d0_IP_MBD3.BPK25_only.summit.bed

wc -l d0_IP_MBD3*summit.bed
#  215 d0_IP_MBD3.BPK25_only.summit.bed
#  527 d0_IP_MBD3.DB_common.summit.bed
# 8012 d0_IP_MBD3.DMSO_only.summit.bed

```

#### d4_DMSO_BPK25 peak intersect

```bash
# d4_DMSO_BPK25 两组 peak 交集
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_MBD3_rep123_together
for i in `find . -name *d4_*narrowPeak |sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
bedtools intersect -a ./d4_DMSO_IP_MBD3_rep123_TvsC/d4_DMSO_IP_MBD3_rep123_TvsC_peaks.narrowPeak.bed -b ./d4_BPK25_IP_MBD3_rep123_TvsC/d4_BPK25_IP_MBD3_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 5412
bedtools intersect -a ./d4_BPK25_IP_MBD3_rep123_TvsC/d4_BPK25_IP_MBD3_rep123_TvsC_peaks.narrowPeak.bed -b ./d4_DMSO_IP_MBD3_rep123_TvsC/d4_DMSO_IP_MBD3_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 5359

# 两组 summit extB50 交集数量
bedtools slop -b 50 -i ./d4_DMSO_IP_MBD3_rep123_TvsC/d4_DMSO_IP_MBD3_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d4_DMSO_IP_MBD3_rep123_TvsC/d4_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed
bedtools slop -b 50 -i ./d4_BPK25_IP_MBD3_rep123_TvsC/d4_BPK25_IP_MBD3_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d4_BPK25_IP_MBD3_rep123_TvsC/d4_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed
# 求交集
bedtools intersect -a ./d4_DMSO_IP_MBD3_rep123_TvsC/d4_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_MBD3_rep123_TvsC/d4_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 5374
bedtools intersect -a ./d4_BPK25_IP_MBD3_rep123_TvsC/d4_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_DMSO_IP_MBD3_rep123_TvsC/d4_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 5363

# 构造 三类 summit
# DMSO only summit
bedtools intersect -a ./d4_DMSO_IP_MBD3_rep123_TvsC/d4_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_MBD3_rep123_TvsC/d4_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d4_DMSO_IP_MBD3_rep123_TvsC/d4_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_MBD3_rep123_TvsC/d4_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d4_IP_MBD3.DMSO_only.summit.bed
# common summit
bedtools intersect -a ./d4_DMSO_IP_MBD3_rep123_TvsC/d4_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_MBD3_rep123_TvsC/d4_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a ./d4_DMSO_IP_MBD3_rep123_TvsC/d4_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_MBD3_rep123_TvsC/d4_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >d4_IP_MBD3.DB_common.summit.bed
# BPK25 only summit
bedtools intersect -a ./d4_BPK25_IP_MBD3_rep123_TvsC/d4_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_DMSO_IP_MBD3_rep123_TvsC/d4_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d4_BPK25_IP_MBD3_rep123_TvsC/d4_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_DMSO_IP_MBD3_rep123_TvsC/d4_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d4_IP_MBD3.BPK25_only.summit.bed

wc -l d4_IP_MBD3*summit.bed
# 12591 d4_IP_MBD3.BPK25_only.summit.bed
#  5401 d4_IP_MBD3.DB_common.summit.bed
# 43909 d4_IP_MBD3.DMSO_only.summit.bed

```

#### d9_DMSO_BPK25 peak intersect

```bash
# d9_DMSO_BPK25 两组 peak 交集
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_MBD3_rep123_together
for i in `find . -name *d9_*narrowPeak |sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
bedtools intersect -a ./d9_DMSO_IP_MBD3_rep123_TvsC/d9_DMSO_IP_MBD3_rep123_TvsC_peaks.narrowPeak.bed -b ./d9_BPK25_IP_MBD3_rep123_TvsC/d9_BPK25_IP_MBD3_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 4332
bedtools intersect -a ./d9_BPK25_IP_MBD3_rep123_TvsC/d9_BPK25_IP_MBD3_rep123_TvsC_peaks.narrowPeak.bed -b ./d9_DMSO_IP_MBD3_rep123_TvsC/d9_DMSO_IP_MBD3_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 4412

# 两组 summit extB50 交集数量
bedtools slop -b 50 -i ./d9_DMSO_IP_MBD3_rep123_TvsC/d9_DMSO_IP_MBD3_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d9_DMSO_IP_MBD3_rep123_TvsC/d9_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed
bedtools slop -b 50 -i ./d9_BPK25_IP_MBD3_rep123_TvsC/d9_BPK25_IP_MBD3_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d9_BPK25_IP_MBD3_rep123_TvsC/d9_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed
# 求交集
bedtools intersect -a ./d9_DMSO_IP_MBD3_rep123_TvsC/d9_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_MBD3_rep123_TvsC/d9_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 4168
bedtools intersect -a ./d9_BPK25_IP_MBD3_rep123_TvsC/d9_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_DMSO_IP_MBD3_rep123_TvsC/d9_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 4199

# 构造 三类 summit
# DMSO only summit
bedtools intersect -a ./d9_DMSO_IP_MBD3_rep123_TvsC/d9_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_MBD3_rep123_TvsC/d9_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d9_DMSO_IP_MBD3_rep123_TvsC/d9_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_MBD3_rep123_TvsC/d9_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d9_IP_MBD3.DMSO_only.summit.bed
# common summit
bedtools intersect -a ./d9_DMSO_IP_MBD3_rep123_TvsC/d9_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_MBD3_rep123_TvsC/d9_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a ./d9_DMSO_IP_MBD3_rep123_TvsC/d9_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_MBD3_rep123_TvsC/d9_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >d9_IP_MBD3.DB_common.summit.bed
# BPK25 only summit
bedtools intersect -a ./d9_BPK25_IP_MBD3_rep123_TvsC/d9_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_DMSO_IP_MBD3_rep123_TvsC/d9_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d9_BPK25_IP_MBD3_rep123_TvsC/d9_BPK25_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_DMSO_IP_MBD3_rep123_TvsC/d9_DMSO_IP_MBD3_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d9_IP_MBD3.BPK25_only.summit.bed

wc -l d9_IP_MBD3*summit.bed
# 43985 d9_IP_MBD3.BPK25_only.summit.bed
#  4208 d9_IP_MBD3.DB_common.summit.bed
# 42254 d9_IP_MBD3.DMSO_only.summit.bed

```



### d049_DMSO_BPK25_TET1_ChIP_Seq

```bash
# data preperation
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir macs2_peak_calling_2_for_TET1_rep123_together
cd macs2_peak_calling_2_for_TET1_rep123_together
# 上传 bam 连接构造脚本并运行
bash 7_make_H9_d049_DMSO_BPK25_TET1_ChIP_Seq_rep123_merged_bam_file_links.sh
ls -1

# d0_DMSO_IP_TET1_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d0_DMSO_IP_TET1_rep123.merged.sorted.bam -c hESC_H9_d0_DMSO_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_DMSO_IP_TET1_rep123_TvsC -n d0_DMSO_IP_TET1_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d0_DMSO_IP_TET1_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d0_DMSO_IP_TET1_rep123_TvsC.runlog
# d0_BPK25_IP_TET1_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d0_BPK25_IP_TET1_rep123.merged.sorted.bam -c hESC_H9_d0_BPK25_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_BPK25_IP_TET1_rep123_TvsC -n d0_BPK25_IP_TET1_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d0_BPK25_IP_TET1_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d0_BPK25_IP_TET1_rep123_TvsC.runlog
# d4_DMSO_IP_TET1_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d4_DMSO_IP_TET1_rep123.merged.sorted.bam -c hESC_H9_d4_DMSO_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_DMSO_IP_TET1_rep123_TvsC -n d4_DMSO_IP_TET1_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d4_DMSO_IP_TET1_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d4_DMSO_IP_TET1_rep123_TvsC.runlog
# d4_BPK25_IP_TET1_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d4_BPK25_IP_TET1_rep123.merged.sorted.bam -c hESC_H9_d4_BPK25_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_BPK25_IP_TET1_rep123_TvsC -n d4_BPK25_IP_TET1_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d4_BPK25_IP_TET1_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d4_BPK25_IP_TET1_rep123_TvsC.runlog
# d9_DMSO_IP_TET1_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d9_DMSO_IP_TET1_rep123.merged.sorted.bam -c hESC_H9_d9_DMSO_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_DMSO_IP_TET1_rep123_TvsC -n d9_DMSO_IP_TET1_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d9_DMSO_IP_TET1_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d9_DMSO_IP_TET1_rep123_TvsC.runlog
# d9_BPK25_IP_TET1_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d9_BPK25_IP_TET1_rep123.merged.sorted.bam -c hESC_H9_d9_BPK25_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_BPK25_IP_TET1_rep123_TvsC -n d9_BPK25_IP_TET1_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d9_BPK25_IP_TET1_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d9_BPK25_IP_TET1_rep123_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_TET1_rep123_together
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./d0_BPK25_IP_TET1_rep123_TvsC/d0_BPK25_IP_TET1_rep123_TvsC_peaks.narrowPeak
# 72148
# ./d0_DMSO_IP_TET1_rep123_TvsC/d0_DMSO_IP_TET1_rep123_TvsC_peaks.narrowPeak
# 8577
# ./d4_BPK25_IP_TET1_rep123_TvsC/d4_BPK25_IP_TET1_rep123_TvsC_peaks.narrowPeak
# 84203
# ./d4_DMSO_IP_TET1_rep123_TvsC/d4_DMSO_IP_TET1_rep123_TvsC_peaks.narrowPeak
# 45019
# ./d9_BPK25_IP_TET1_rep123_TvsC/d9_BPK25_IP_TET1_rep123_TvsC_peaks.narrowPeak
# 43694
# ./d9_DMSO_IP_TET1_rep123_TvsC/d9_DMSO_IP_TET1_rep123_TvsC_peaks.narrowPeak
# 24572

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./d0_BPK25_IP_TET1_rep123_TvsC/d0_BPK25_IP_TET1_rep123_TvsC_summits.bed
# 77784
# ./d0_DMSO_IP_TET1_rep123_TvsC/d0_DMSO_IP_TET1_rep123_TvsC_summits.bed
# 10009
# ./d4_BPK25_IP_TET1_rep123_TvsC/d4_BPK25_IP_TET1_rep123_TvsC_summits.bed
# 89856
# ./d4_DMSO_IP_TET1_rep123_TvsC/d4_DMSO_IP_TET1_rep123_TvsC_summits.bed
# 47511
# ./d9_BPK25_IP_TET1_rep123_TvsC/d9_BPK25_IP_TET1_rep123_TvsC_summits.bed
# 46382
# ./d9_DMSO_IP_TET1_rep123_TvsC/d9_DMSO_IP_TET1_rep123_TvsC_summits.bed
# 26760

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_d049_BPK25_DMSO_TET1_ChIP_Seq_macs2_results_summary `find . -name *peaks.xls`


cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_TET1_rep123_together
# 将bdg转为bigwig用于可视化
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done

```



#### d0_DMSO_BPK25 peak intersect

```bash
# d0_DMSO_BPK25 两组 peak 交集
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_TET1_rep123_together
for i in `find . -name *d0*narrowPeak |sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
bedtools intersect -a ./d0_DMSO_IP_TET1_rep123_TvsC/d0_DMSO_IP_TET1_rep123_TvsC_peaks.narrowPeak.bed -b ./d0_BPK25_IP_TET1_rep123_TvsC/d0_BPK25_IP_TET1_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 4865
bedtools intersect -a ./d0_BPK25_IP_TET1_rep123_TvsC/d0_BPK25_IP_TET1_rep123_TvsC_peaks.narrowPeak.bed -b ./d0_DMSO_IP_TET1_rep123_TvsC/d0_DMSO_IP_TET1_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 4637

# 两组 summit extB50 交集数量
bedtools slop -b 50 -i ./d0_DMSO_IP_TET1_rep123_TvsC/d0_DMSO_IP_TET1_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d0_DMSO_IP_TET1_rep123_TvsC/d0_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed
bedtools slop -b 50 -i ./d0_BPK25_IP_TET1_rep123_TvsC/d0_BPK25_IP_TET1_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d0_BPK25_IP_TET1_rep123_TvsC/d0_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed
# 求交集
bedtools intersect -a ./d0_DMSO_IP_TET1_rep123_TvsC/d0_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_TET1_rep123_TvsC/d0_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 4641
bedtools intersect -a ./d0_BPK25_IP_TET1_rep123_TvsC/d0_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_DMSO_IP_TET1_rep123_TvsC/d0_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 4626

# 构造 三类 summit
# DMSO only summit
bedtools intersect -a ./d0_DMSO_IP_TET1_rep123_TvsC/d0_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_TET1_rep123_TvsC/d0_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d0_DMSO_IP_TET1_rep123_TvsC/d0_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_TET1_rep123_TvsC/d0_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d0_IP_TET1.DMSO_only.summit.bed
# common summit
bedtools intersect -a ./d0_DMSO_IP_TET1_rep123_TvsC/d0_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_TET1_rep123_TvsC/d0_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a ./d0_DMSO_IP_TET1_rep123_TvsC/d0_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_TET1_rep123_TvsC/d0_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >d0_IP_TET1.DB_common.summit.bed
# BPK25 only summit
bedtools intersect -a ./d0_BPK25_IP_TET1_rep123_TvsC/d0_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_DMSO_IP_TET1_rep123_TvsC/d0_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d0_BPK25_IP_TET1_rep123_TvsC/d0_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_DMSO_IP_TET1_rep123_TvsC/d0_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d0_IP_TET1.BPK25_only.summit.bed

wc -l d0_IP_TET1*summit.bed
# 73158 d0_IP_TET1.BPK25_only.summit.bed
#  4677 d0_IP_TET1.DB_common.summit.bed
#  5368 d0_IP_TET1.DMSO_only.summit.bed

```

#### d4_DMSO_BPK25 peak intersect

```bash
# d4_DMSO_BPK25 两组 peak 交集
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_TET1_rep123_together
for i in `find . -name *d4_*narrowPeak |sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
bedtools intersect -a ./d4_DMSO_IP_TET1_rep123_TvsC/d4_DMSO_IP_TET1_rep123_TvsC_peaks.narrowPeak.bed -b ./d4_BPK25_IP_TET1_rep123_TvsC/d4_BPK25_IP_TET1_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 14784
bedtools intersect -a ./d4_BPK25_IP_TET1_rep123_TvsC/d4_BPK25_IP_TET1_rep123_TvsC_peaks.narrowPeak.bed -b ./d4_DMSO_IP_TET1_rep123_TvsC/d4_DMSO_IP_TET1_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 14873

# 两组 summit extB50 交集数量
bedtools slop -b 50 -i ./d4_DMSO_IP_TET1_rep123_TvsC/d4_DMSO_IP_TET1_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d4_DMSO_IP_TET1_rep123_TvsC/d4_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed
bedtools slop -b 50 -i ./d4_BPK25_IP_TET1_rep123_TvsC/d4_BPK25_IP_TET1_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d4_BPK25_IP_TET1_rep123_TvsC/d4_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed
# 求交集
bedtools intersect -a ./d4_DMSO_IP_TET1_rep123_TvsC/d4_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_TET1_rep123_TvsC/d4_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 11241
bedtools intersect -a ./d4_BPK25_IP_TET1_rep123_TvsC/d4_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_DMSO_IP_TET1_rep123_TvsC/d4_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 11267

# 构造 三类 summit
# DMSO only summit
bedtools intersect -a ./d4_DMSO_IP_TET1_rep123_TvsC/d4_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_TET1_rep123_TvsC/d4_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d4_DMSO_IP_TET1_rep123_TvsC/d4_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_TET1_rep123_TvsC/d4_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d4_IP_TET1.DMSO_only.summit.bed
# common summit
bedtools intersect -a ./d4_DMSO_IP_TET1_rep123_TvsC/d4_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_TET1_rep123_TvsC/d4_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a ./d4_DMSO_IP_TET1_rep123_TvsC/d4_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_TET1_rep123_TvsC/d4_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >d4_IP_TET1.DB_common.summit.bed
# BPK25 only summit
bedtools intersect -a ./d4_BPK25_IP_TET1_rep123_TvsC/d4_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_DMSO_IP_TET1_rep123_TvsC/d4_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d4_BPK25_IP_TET1_rep123_TvsC/d4_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_DMSO_IP_TET1_rep123_TvsC/d4_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d4_IP_TET1.BPK25_only.summit.bed

wc -l d4_IP_TET1*summit.bed
# 78589 d4_IP_TET1.BPK25_only.summit.bed
# 11308 d4_IP_TET1.DB_common.summit.bed
# 36270 d4_IP_TET1.DMSO_only.summit.bed

```

#### d9_DMSO_BPK25 peak intersect

```bash
# d9_DMSO_BPK25 两组 peak 交集
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_TET1_rep123_together
for i in `find . -name *d9_*narrowPeak |sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
bedtools intersect -a ./d9_DMSO_IP_TET1_rep123_TvsC/d9_DMSO_IP_TET1_rep123_TvsC_peaks.narrowPeak.bed -b ./d9_BPK25_IP_TET1_rep123_TvsC/d9_BPK25_IP_TET1_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 2669
bedtools intersect -a ./d9_BPK25_IP_TET1_rep123_TvsC/d9_BPK25_IP_TET1_rep123_TvsC_peaks.narrowPeak.bed -b ./d9_DMSO_IP_TET1_rep123_TvsC/d9_DMSO_IP_TET1_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 2667

# 两组 summit extB50 交集数量
bedtools slop -b 50 -i ./d9_DMSO_IP_TET1_rep123_TvsC/d9_DMSO_IP_TET1_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d9_DMSO_IP_TET1_rep123_TvsC/d9_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed
bedtools slop -b 50 -i ./d9_BPK25_IP_TET1_rep123_TvsC/d9_BPK25_IP_TET1_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d9_BPK25_IP_TET1_rep123_TvsC/d9_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed
# 求交集
bedtools intersect -a ./d9_DMSO_IP_TET1_rep123_TvsC/d9_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_TET1_rep123_TvsC/d9_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 1828
bedtools intersect -a ./d9_BPK25_IP_TET1_rep123_TvsC/d9_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_DMSO_IP_TET1_rep123_TvsC/d9_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 1815

# 构造 三类 summit
# DMSO only summit
bedtools intersect -a ./d9_DMSO_IP_TET1_rep123_TvsC/d9_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_TET1_rep123_TvsC/d9_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d9_DMSO_IP_TET1_rep123_TvsC/d9_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_TET1_rep123_TvsC/d9_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d9_IP_TET1.DMSO_only.summit.bed
# common summit
bedtools intersect -a ./d9_DMSO_IP_TET1_rep123_TvsC/d9_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_TET1_rep123_TvsC/d9_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a ./d9_DMSO_IP_TET1_rep123_TvsC/d9_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_TET1_rep123_TvsC/d9_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >d9_IP_TET1.DB_common.summit.bed
# BPK25 only summit
bedtools intersect -a ./d9_BPK25_IP_TET1_rep123_TvsC/d9_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_DMSO_IP_TET1_rep123_TvsC/d9_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d9_BPK25_IP_TET1_rep123_TvsC/d9_BPK25_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_DMSO_IP_TET1_rep123_TvsC/d9_DMSO_IP_TET1_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d9_IP_TET1.BPK25_only.summit.bed

wc -l d9_IP_TET1*summit.bed
  # 44567 d9_IP_TET1.BPK25_only.summit.bed
  #  1837 d9_IP_TET1.DB_common.summit.bed
  # 24932 d9_IP_TET1.DMSO_only.summit.bed

```



### d049_DMSO_BPK25_TET2_ChIP_Seq

```bash
# data preperation
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir macs2_peak_calling_2_for_TET2_rep123_together
cd macs2_peak_calling_2_for_TET2_rep123_together
# 上传 bam 连接构造脚本并运行
bash 7_make_H9_d049_DMSO_BPK25_TET2_ChIP_Seq_rep123_merged_bam_file_links.sh
ls -1

# d0_DMSO_IP_TET2_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d0_DMSO_IP_TET2_rep123.merged.sorted.bam -c hESC_H9_d0_DMSO_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_DMSO_IP_TET2_rep123_TvsC -n d0_DMSO_IP_TET2_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d0_DMSO_IP_TET2_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d0_DMSO_IP_TET2_rep123_TvsC.runlog
# d0_BPK25_IP_TET2_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d0_BPK25_IP_TET2_rep123.merged.sorted.bam -c hESC_H9_d0_BPK25_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_BPK25_IP_TET2_rep123_TvsC -n d0_BPK25_IP_TET2_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d0_BPK25_IP_TET2_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d0_BPK25_IP_TET2_rep123_TvsC.runlog
# d4_DMSO_IP_TET2_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d4_DMSO_IP_TET2_rep123.merged.sorted.bam -c hESC_H9_d4_DMSO_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_DMSO_IP_TET2_rep123_TvsC -n d4_DMSO_IP_TET2_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d4_DMSO_IP_TET2_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d4_DMSO_IP_TET2_rep123_TvsC.runlog
# d4_BPK25_IP_TET2_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d4_BPK25_IP_TET2_rep123.merged.sorted.bam -c hESC_H9_d4_BPK25_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_BPK25_IP_TET2_rep123_TvsC -n d4_BPK25_IP_TET2_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d4_BPK25_IP_TET2_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d4_BPK25_IP_TET2_rep123_TvsC.runlog
# d9_DMSO_IP_TET2_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d9_DMSO_IP_TET2_rep123.merged.sorted.bam -c hESC_H9_d9_DMSO_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_DMSO_IP_TET2_rep123_TvsC -n d9_DMSO_IP_TET2_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d9_DMSO_IP_TET2_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d9_DMSO_IP_TET2_rep123_TvsC.runlog
# d9_BPK25_IP_TET2_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d9_BPK25_IP_TET2_rep123.merged.sorted.bam -c hESC_H9_d9_BPK25_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_BPK25_IP_TET2_rep123_TvsC -n d9_BPK25_IP_TET2_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d9_BPK25_IP_TET2_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d9_BPK25_IP_TET2_rep123_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_TET2_rep123_together
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./d0_BPK25_IP_TET2_rep123_TvsC/d0_BPK25_IP_TET2_rep123_TvsC_peaks.narrowPeak
# 31522
# ./d0_DMSO_IP_TET2_rep123_TvsC/d0_DMSO_IP_TET2_rep123_TvsC_peaks.narrowPeak
# 3485
# ./d4_BPK25_IP_TET2_rep123_TvsC/d4_BPK25_IP_TET2_rep123_TvsC_peaks.narrowPeak
# 13324
# ./d4_DMSO_IP_TET2_rep123_TvsC/d4_DMSO_IP_TET2_rep123_TvsC_peaks.narrowPeak
# 36633
# ./d9_BPK25_IP_TET2_rep123_TvsC/d9_BPK25_IP_TET2_rep123_TvsC_peaks.narrowPeak
# 1928
# ./d9_DMSO_IP_TET2_rep123_TvsC/d9_DMSO_IP_TET2_rep123_TvsC_peaks.narrowPeak
# 15250

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./d0_BPK25_IP_TET2_rep123_TvsC/d0_BPK25_IP_TET2_rep123_TvsC_summits.bed
# 34450
# ./d0_DMSO_IP_TET2_rep123_TvsC/d0_DMSO_IP_TET2_rep123_TvsC_summits.bed
# 4006
# ./d4_BPK25_IP_TET2_rep123_TvsC/d4_BPK25_IP_TET2_rep123_TvsC_summits.bed
# 14687
# ./d4_DMSO_IP_TET2_rep123_TvsC/d4_DMSO_IP_TET2_rep123_TvsC_summits.bed
# 38606
# ./d9_BPK25_IP_TET2_rep123_TvsC/d9_BPK25_IP_TET2_rep123_TvsC_summits.bed
# 2030
# ./d9_DMSO_IP_TET2_rep123_TvsC/d9_DMSO_IP_TET2_rep123_TvsC_summits.bed
# 15734

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_d049_BPK25_DMSO_TET2_ChIP_Seq_macs2_results_summary `find . -name *peaks.xls`

cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_TET2_rep123_together
# 将bdg转为bigwig用于可视化
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done

```



#### d0_DMSO_BPK25 peak intersect

```bash
# d0_DMSO_BPK25 两组 peak 交集
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_TET2_rep123_together
for i in `find . -name *d0*narrowPeak |sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
bedtools intersect -a ./d0_DMSO_IP_TET2_rep123_TvsC/d0_DMSO_IP_TET2_rep123_TvsC_peaks.narrowPeak.bed -b ./d0_BPK25_IP_TET2_rep123_TvsC/d0_BPK25_IP_TET2_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 1817
bedtools intersect -a ./d0_BPK25_IP_TET2_rep123_TvsC/d0_BPK25_IP_TET2_rep123_TvsC_peaks.narrowPeak.bed -b ./d0_DMSO_IP_TET2_rep123_TvsC/d0_DMSO_IP_TET2_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 1770

# 两组 summit extB50 交集数量
bedtools slop -b 50 -i ./d0_DMSO_IP_TET2_rep123_TvsC/d0_DMSO_IP_TET2_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d0_DMSO_IP_TET2_rep123_TvsC/d0_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed
bedtools slop -b 50 -i ./d0_BPK25_IP_TET2_rep123_TvsC/d0_BPK25_IP_TET2_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d0_BPK25_IP_TET2_rep123_TvsC/d0_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed
# 求交集
bedtools intersect -a ./d0_DMSO_IP_TET2_rep123_TvsC/d0_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_TET2_rep123_TvsC/d0_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 1772
bedtools intersect -a ./d0_BPK25_IP_TET2_rep123_TvsC/d0_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_DMSO_IP_TET2_rep123_TvsC/d0_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 1776

# 构造 三类 summit
# DMSO only summit
bedtools intersect -a ./d0_DMSO_IP_TET2_rep123_TvsC/d0_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_TET2_rep123_TvsC/d0_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d0_DMSO_IP_TET2_rep123_TvsC/d0_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_TET2_rep123_TvsC/d0_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d0_IP_TET2.DMSO_only.summit.bed
# common summit
bedtools intersect -a ./d0_DMSO_IP_TET2_rep123_TvsC/d0_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_TET2_rep123_TvsC/d0_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a ./d0_DMSO_IP_TET2_rep123_TvsC/d0_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_TET2_rep123_TvsC/d0_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >d0_IP_TET2.DB_common.summit.bed
# BPK25 only summit
bedtools intersect -a ./d0_BPK25_IP_TET2_rep123_TvsC/d0_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_DMSO_IP_TET2_rep123_TvsC/d0_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d0_BPK25_IP_TET2_rep123_TvsC/d0_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_DMSO_IP_TET2_rep123_TvsC/d0_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d0_IP_TET2.BPK25_only.summit.bed

wc -l d0_IP_TET2*summit.bed
# 32674 d0_IP_TET2.BPK25_only.summit.bed
#  1793 d0_IP_TET2.DB_common.summit.bed
#  2234 d0_IP_TET2.DMSO_only.summit.bed

```

#### d4_DMSO_BPK25 peak intersect

```bash
# d4_DMSO_BPK25 两组 peak 交集
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_TET2_rep123_together
for i in `find . -name *d4_*narrowPeak |sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
bedtools intersect -a ./d4_DMSO_IP_TET2_rep123_TvsC/d4_DMSO_IP_TET2_rep123_TvsC_peaks.narrowPeak.bed -b ./d4_BPK25_IP_TET2_rep123_TvsC/d4_BPK25_IP_TET2_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 1901
bedtools intersect -a ./d4_BPK25_IP_TET2_rep123_TvsC/d4_BPK25_IP_TET2_rep123_TvsC_peaks.narrowPeak.bed -b ./d4_DMSO_IP_TET2_rep123_TvsC/d4_DMSO_IP_TET2_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 1910

# 两组 summit extB50 交集数量
bedtools slop -b 50 -i ./d4_DMSO_IP_TET2_rep123_TvsC/d4_DMSO_IP_TET2_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d4_DMSO_IP_TET2_rep123_TvsC/d4_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed
bedtools slop -b 50 -i ./d4_BPK25_IP_TET2_rep123_TvsC/d4_BPK25_IP_TET2_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d4_BPK25_IP_TET2_rep123_TvsC/d4_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed
# 求交集
bedtools intersect -a ./d4_DMSO_IP_TET2_rep123_TvsC/d4_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_TET2_rep123_TvsC/d4_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 1485
bedtools intersect -a ./d4_BPK25_IP_TET2_rep123_TvsC/d4_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_DMSO_IP_TET2_rep123_TvsC/d4_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 1479

# 构造 三类 summit
# DMSO only summit
bedtools intersect -a ./d4_DMSO_IP_TET2_rep123_TvsC/d4_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_TET2_rep123_TvsC/d4_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d4_DMSO_IP_TET2_rep123_TvsC/d4_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_TET2_rep123_TvsC/d4_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d4_IP_TET2.DMSO_only.summit.bed
# common summit
bedtools intersect -a ./d4_DMSO_IP_TET2_rep123_TvsC/d4_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_TET2_rep123_TvsC/d4_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a ./d4_DMSO_IP_TET2_rep123_TvsC/d4_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_TET2_rep123_TvsC/d4_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >d4_IP_TET2.DB_common.summit.bed
# BPK25 only summit
bedtools intersect -a ./d4_BPK25_IP_TET2_rep123_TvsC/d4_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_DMSO_IP_TET2_rep123_TvsC/d4_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d4_BPK25_IP_TET2_rep123_TvsC/d4_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_DMSO_IP_TET2_rep123_TvsC/d4_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d4_IP_TET2.BPK25_only.summit.bed

wc -l d4_IP_TET2*summit.bed
# 13208 d4_IP_TET2.BPK25_only.summit.bed
#  1487 d4_IP_TET2.DB_common.summit.bed
# 37121 d4_IP_TET2.DMSO_only.summit.bed

```

#### d9_DMSO_BPK25 peak intersect

```bash
# d9_DMSO_BPK25 两组 peak 交集
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_TET2_rep123_together
for i in `find . -name *d9_*narrowPeak |sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
bedtools intersect -a ./d9_DMSO_IP_TET2_rep123_TvsC/d9_DMSO_IP_TET2_rep123_TvsC_peaks.narrowPeak.bed -b ./d9_BPK25_IP_TET2_rep123_TvsC/d9_BPK25_IP_TET2_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 643
bedtools intersect -a ./d9_BPK25_IP_TET2_rep123_TvsC/d9_BPK25_IP_TET2_rep123_TvsC_peaks.narrowPeak.bed -b ./d9_DMSO_IP_TET2_rep123_TvsC/d9_DMSO_IP_TET2_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 646

# 两组 summit extB50 交集数量
bedtools slop -b 50 -i ./d9_DMSO_IP_TET2_rep123_TvsC/d9_DMSO_IP_TET2_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d9_DMSO_IP_TET2_rep123_TvsC/d9_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed
bedtools slop -b 50 -i ./d9_BPK25_IP_TET2_rep123_TvsC/d9_BPK25_IP_TET2_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d9_BPK25_IP_TET2_rep123_TvsC/d9_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed
# 求交集
bedtools intersect -a ./d9_DMSO_IP_TET2_rep123_TvsC/d9_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_TET2_rep123_TvsC/d9_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 491
bedtools intersect -a ./d9_BPK25_IP_TET2_rep123_TvsC/d9_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_DMSO_IP_TET2_rep123_TvsC/d9_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 490

# 构造 三类 summit
# DMSO only summit
bedtools intersect -a ./d9_DMSO_IP_TET2_rep123_TvsC/d9_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_TET2_rep123_TvsC/d9_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d9_DMSO_IP_TET2_rep123_TvsC/d9_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_TET2_rep123_TvsC/d9_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d9_IP_TET2.DMSO_only.summit.bed
# common summit
bedtools intersect -a ./d9_DMSO_IP_TET2_rep123_TvsC/d9_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_TET2_rep123_TvsC/d9_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a ./d9_DMSO_IP_TET2_rep123_TvsC/d9_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_TET2_rep123_TvsC/d9_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >d9_IP_TET2.DB_common.summit.bed
# BPK25 only summit
bedtools intersect -a ./d9_BPK25_IP_TET2_rep123_TvsC/d9_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_DMSO_IP_TET2_rep123_TvsC/d9_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d9_BPK25_IP_TET2_rep123_TvsC/d9_BPK25_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_DMSO_IP_TET2_rep123_TvsC/d9_DMSO_IP_TET2_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d9_IP_TET2.BPK25_only.summit.bed

wc -l d9_IP_TET2*summit.bed
#  1540 d9_IP_TET2.BPK25_only.summit.bed
#   493 d9_IP_TET2.DB_common.summit.bed
# 15243 d9_IP_TET2.DMSO_only.summit.bed

```



### d049_DMSO_BPK25_TET3_ChIP_Seq

```bash
# data preperation
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir macs2_peak_calling_2_for_TET3_rep123_together
cd macs2_peak_calling_2_for_TET3_rep123_together
# 上传 bam 连接构造脚本并运行
bash 7_make_H9_d049_DMSO_BPK25_TET3_ChIP_Seq_rep123_merged_bam_file_links.sh
ls -1

# d0_DMSO_IP_TET3_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d0_DMSO_IP_TET3_rep123.merged.sorted.bam -c hESC_H9_d0_DMSO_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_DMSO_IP_TET3_rep123_TvsC -n d0_DMSO_IP_TET3_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d0_DMSO_IP_TET3_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d0_DMSO_IP_TET3_rep123_TvsC.runlog
# d0_BPK25_IP_TET3_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d0_BPK25_IP_TET3_rep123.merged.sorted.bam -c hESC_H9_d0_BPK25_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d0_BPK25_IP_TET3_rep123_TvsC -n d0_BPK25_IP_TET3_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d0_BPK25_IP_TET3_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d0_BPK25_IP_TET3_rep123_TvsC.runlog
# d4_DMSO_IP_TET3_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d4_DMSO_IP_TET3_rep123.merged.sorted.bam -c hESC_H9_d4_DMSO_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_DMSO_IP_TET3_rep123_TvsC -n d4_DMSO_IP_TET3_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d4_DMSO_IP_TET3_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d4_DMSO_IP_TET3_rep123_TvsC.runlog
# d4_BPK25_IP_TET3_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d4_BPK25_IP_TET3_rep123.merged.sorted.bam -c hESC_H9_d4_BPK25_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d4_BPK25_IP_TET3_rep123_TvsC -n d4_BPK25_IP_TET3_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d4_BPK25_IP_TET3_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d4_BPK25_IP_TET3_rep123_TvsC.runlog
# d9_DMSO_IP_TET3_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d9_DMSO_IP_TET3_rep123.merged.sorted.bam -c hESC_H9_d9_DMSO_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_DMSO_IP_TET3_rep123_TvsC -n d9_DMSO_IP_TET3_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d9_DMSO_IP_TET3_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d9_DMSO_IP_TET3_rep123_TvsC.runlog
# d9_BPK25_IP_TET3_rep123_TvsC
nohup macs2 callpeak -t hESC_H9_d9_BPK25_IP_TET3_rep123.merged.sorted.bam -c hESC_H9_d9_BPK25_input_rep123.merged.sorted.bam --format BAMPE -g hs --keep-dup 1 --outdir d9_BPK25_IP_TET3_rep123_TvsC -n d9_BPK25_IP_TET3_rep123_TvsC -B --call-summits --cutoff-analysis >macs2_callpeak_for_d9_BPK25_IP_TET3_rep123_TvsC.runlog 2>&1 &
tail -f macs2_callpeak_for_d9_BPK25_IP_TET3_rep123_TvsC.runlog

# 统计peak数量
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_TET3_rep123_together
for i in `find . -name *narrowPeak|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./d0_BPK25_IP_TET3_rep123_TvsC/d0_BPK25_IP_TET3_rep123_TvsC_peaks.narrowPeak
# 15244
# ./d0_DMSO_IP_TET3_rep123_TvsC/d0_DMSO_IP_TET3_rep123_TvsC_peaks.narrowPeak
# 37154
# ./d4_BPK25_IP_TET3_rep123_TvsC/d4_BPK25_IP_TET3_rep123_TvsC_peaks.narrowPeak
# 12757
# ./d4_DMSO_IP_TET3_rep123_TvsC/d4_DMSO_IP_TET3_rep123_TvsC_peaks.narrowPeak
# 33121
# ./d9_BPK25_IP_TET3_rep123_TvsC/d9_BPK25_IP_TET3_rep123_TvsC_peaks.narrowPeak
# 321
# ./d9_DMSO_IP_TET3_rep123_TvsC/d9_DMSO_IP_TET3_rep123_TvsC_peaks.narrowPeak
# 25543

# 统计 peak summit 数量
for i in `find . -name **summits.bed|sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq | wc -l ;done
# ./d0_BPK25_IP_TET3_rep123_TvsC/d0_BPK25_IP_TET3_rep123_TvsC_summits.bed
# 16977
# ./d0_DMSO_IP_TET3_rep123_TvsC/d0_DMSO_IP_TET3_rep123_TvsC_summits.bed
# 40647
# ./d4_BPK25_IP_TET3_rep123_TvsC/d4_BPK25_IP_TET3_rep123_TvsC_summits.bed
# 13738
# ./d4_DMSO_IP_TET3_rep123_TvsC/d4_DMSO_IP_TET3_rep123_TvsC_summits.bed
# 36387
# ./d9_BPK25_IP_TET3_rep123_TvsC/d9_BPK25_IP_TET3_rep123_TvsC_summits.bed
# 355
# ./d9_DMSO_IP_TET3_rep123_TvsC/d9_DMSO_IP_TET3_rep123_TvsC_summits.bed
# 28493

# 也可使用 multiqc 对 macs2 的 peak结果进行汇总
multiqc --no-megaqc-upload -m macs2 -o all_d049_BPK25_DMSO_TET3_ChIP_Seq_macs2_results_summary `find . -name *peaks.xls`

cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_TET3_rep123_together
# 将bdg转为bigwig用于可视化
for i in `find . -name *.bdg|sort`;do echo $i; sort -k1,1 -k2,2n $i >${i%.bdg}.sorted.bdg; bedGraphToBigWig ${i%.bdg}.sorted.bdg /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes ${i%.bdg}.bigwig;done
# bigwig 重命名

# 删除 bdg 文件
for i in `find . -name *.bdg|sort`;do echo $i;done
for i in `find . -name *.bdg|sort`;do echo $i; rm $i;done
for i in `find . -name *.bdg|sort`;do echo $i;done

```



#### d0_DMSO_BPK25 peak intersect

```bash
# d0_DMSO_BPK25 两组 peak 交集
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_TET3_rep123_together
for i in `find . -name *d0*narrowPeak |sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
bedtools intersect -a ./d0_DMSO_IP_TET3_rep123_TvsC/d0_DMSO_IP_TET3_rep123_TvsC_peaks.narrowPeak.bed -b ./d0_BPK25_IP_TET3_rep123_TvsC/d0_BPK25_IP_TET3_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 6690
bedtools intersect -a ./d0_BPK25_IP_TET3_rep123_TvsC/d0_BPK25_IP_TET3_rep123_TvsC_peaks.narrowPeak.bed -b ./d0_DMSO_IP_TET3_rep123_TvsC/d0_DMSO_IP_TET3_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 6954

# 两组 summit extB50 交集数量
bedtools slop -b 50 -i ./d0_DMSO_IP_TET3_rep123_TvsC/d0_DMSO_IP_TET3_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d0_DMSO_IP_TET3_rep123_TvsC/d0_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed
bedtools slop -b 50 -i ./d0_BPK25_IP_TET3_rep123_TvsC/d0_BPK25_IP_TET3_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d0_BPK25_IP_TET3_rep123_TvsC/d0_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed
# 求交集
bedtools intersect -a ./d0_DMSO_IP_TET3_rep123_TvsC/d0_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_TET3_rep123_TvsC/d0_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 6600
bedtools intersect -a ./d0_BPK25_IP_TET3_rep123_TvsC/d0_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_DMSO_IP_TET3_rep123_TvsC/d0_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 6589

# 构造 三类 summit
# DMSO only summit
bedtools intersect -a ./d0_DMSO_IP_TET3_rep123_TvsC/d0_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_TET3_rep123_TvsC/d0_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d0_DMSO_IP_TET3_rep123_TvsC/d0_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_TET3_rep123_TvsC/d0_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d0_IP_TET3.DMSO_only.summit.bed
# common summit
bedtools intersect -a ./d0_DMSO_IP_TET3_rep123_TvsC/d0_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_TET3_rep123_TvsC/d0_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a ./d0_DMSO_IP_TET3_rep123_TvsC/d0_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_BPK25_IP_TET3_rep123_TvsC/d0_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >d0_IP_TET3.DB_common.summit.bed
# BPK25 only summit
bedtools intersect -a ./d0_BPK25_IP_TET3_rep123_TvsC/d0_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_DMSO_IP_TET3_rep123_TvsC/d0_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d0_BPK25_IP_TET3_rep123_TvsC/d0_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d0_DMSO_IP_TET3_rep123_TvsC/d0_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d0_IP_TET3.BPK25_only.summit.bed

wc -l d0_IP_TET3*summit.bed
# 10388 d0_IP_TET3.BPK25_only.summit.bed
#  6651 d0_IP_TET3.DB_common.summit.bed
# 34047 d0_IP_TET3.DMSO_only.summit.bed

```

#### d4_DMSO_BPK25 peak intersect

```bash
# d4_DMSO_BPK25 两组 peak 交集
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_TET3_rep123_together
for i in `find . -name *d4_*narrowPeak |sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
bedtools intersect -a ./d4_DMSO_IP_TET3_rep123_TvsC/d4_DMSO_IP_TET3_rep123_TvsC_peaks.narrowPeak.bed -b ./d4_BPK25_IP_TET3_rep123_TvsC/d4_BPK25_IP_TET3_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 6232
bedtools intersect -a ./d4_BPK25_IP_TET3_rep123_TvsC/d4_BPK25_IP_TET3_rep123_TvsC_peaks.narrowPeak.bed -b ./d4_DMSO_IP_TET3_rep123_TvsC/d4_DMSO_IP_TET3_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 6277

# 两组 summit extB50 交集数量
bedtools slop -b 50 -i ./d4_DMSO_IP_TET3_rep123_TvsC/d4_DMSO_IP_TET3_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d4_DMSO_IP_TET3_rep123_TvsC/d4_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed
bedtools slop -b 50 -i ./d4_BPK25_IP_TET3_rep123_TvsC/d4_BPK25_IP_TET3_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d4_BPK25_IP_TET3_rep123_TvsC/d4_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed
# 求交集
bedtools intersect -a ./d4_DMSO_IP_TET3_rep123_TvsC/d4_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_TET3_rep123_TvsC/d4_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 4345
bedtools intersect -a ./d4_BPK25_IP_TET3_rep123_TvsC/d4_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_DMSO_IP_TET3_rep123_TvsC/d4_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 4354

# 构造 三类 summit
# DMSO only summit
bedtools intersect -a ./d4_DMSO_IP_TET3_rep123_TvsC/d4_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_TET3_rep123_TvsC/d4_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d4_DMSO_IP_TET3_rep123_TvsC/d4_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_TET3_rep123_TvsC/d4_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d4_IP_TET3.DMSO_only.summit.bed
# common summit
bedtools intersect -a ./d4_DMSO_IP_TET3_rep123_TvsC/d4_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_TET3_rep123_TvsC/d4_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a ./d4_DMSO_IP_TET3_rep123_TvsC/d4_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_BPK25_IP_TET3_rep123_TvsC/d4_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >d4_IP_TET3.DB_common.summit.bed
# BPK25 only summit
bedtools intersect -a ./d4_BPK25_IP_TET3_rep123_TvsC/d4_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_DMSO_IP_TET3_rep123_TvsC/d4_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d4_BPK25_IP_TET3_rep123_TvsC/d4_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d4_DMSO_IP_TET3_rep123_TvsC/d4_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d4_IP_TET3.BPK25_only.summit.bed

wc -l d4_IP_TET3*summit.bed
#  9384 d4_IP_TET3.BPK25_only.summit.bed
#  4379 d4_IP_TET3.DB_common.summit.bed
# 32042 d4_IP_TET3.DMSO_only.summit.bed

```

#### d9_DMSO_BPK25 peak intersect

```bash
# d9_DMSO_BPK25 两组 peak 交集
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_TET3_rep123_together
for i in `find . -name *d9_*narrowPeak |sort`;do echo $i; cut -f 1,2,3 $i |sort -k1,1 -k2,2n |uniq >${i}.bed ;done
bedtools intersect -a ./d9_DMSO_IP_TET3_rep123_TvsC/d9_DMSO_IP_TET3_rep123_TvsC_peaks.narrowPeak.bed -b ./d9_BPK25_IP_TET3_rep123_TvsC/d9_BPK25_IP_TET3_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 175
bedtools intersect -a ./d9_BPK25_IP_TET3_rep123_TvsC/d9_BPK25_IP_TET3_rep123_TvsC_peaks.narrowPeak.bed -b ./d9_DMSO_IP_TET3_rep123_TvsC/d9_DMSO_IP_TET3_rep123_TvsC_peaks.narrowPeak.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 181

# 两组 summit extB50 交集数量
bedtools slop -b 50 -i ./d9_DMSO_IP_TET3_rep123_TvsC/d9_DMSO_IP_TET3_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d9_DMSO_IP_TET3_rep123_TvsC/d9_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed
bedtools slop -b 50 -i ./d9_BPK25_IP_TET3_rep123_TvsC/d9_BPK25_IP_TET3_rep123_TvsC_summits.bed -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes |cut -f 1,2,3 |sort -k1,1 -k2,2n |uniq >./d9_BPK25_IP_TET3_rep123_TvsC/d9_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed
# 求交集
bedtools intersect -a ./d9_DMSO_IP_TET3_rep123_TvsC/d9_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_TET3_rep123_TvsC/d9_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 159
bedtools intersect -a ./d9_BPK25_IP_TET3_rep123_TvsC/d9_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_DMSO_IP_TET3_rep123_TvsC/d9_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -wa |sort -k1,1 -k2,2n |uniq |wc -l
# 158

# 构造 三类 summit
# DMSO only summit
bedtools intersect -a ./d9_DMSO_IP_TET3_rep123_TvsC/d9_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_TET3_rep123_TvsC/d9_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d9_DMSO_IP_TET3_rep123_TvsC/d9_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_TET3_rep123_TvsC/d9_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d9_IP_TET3.DMSO_only.summit.bed
# common summit
bedtools intersect -a ./d9_DMSO_IP_TET3_rep123_TvsC/d9_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_TET3_rep123_TvsC/d9_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed |head
bedtools intersect -a ./d9_DMSO_IP_TET3_rep123_TvsC/d9_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_BPK25_IP_TET3_rep123_TvsC/d9_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed | awk 'BEGIN{FS="\t";OFS="\t"}{center = int(($2+$3)/2); print $1, center, center + 1}' |sort -k1,1 -k2,2n  |uniq >d9_IP_TET3.DB_common.summit.bed
# BPK25 only summit
bedtools intersect -a ./d9_BPK25_IP_TET3_rep123_TvsC/d9_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_DMSO_IP_TET3_rep123_TvsC/d9_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -wao |head
bedtools intersect -a ./d9_BPK25_IP_TET3_rep123_TvsC/d9_BPK25_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -b ./d9_DMSO_IP_TET3_rep123_TvsC/d9_DMSO_IP_TET3_rep123_TvsC_summits.extB50.sorted.bed -wao | awk 'BEGIN{FS="\t";OFS="\t"}{ if($7==0){center = int(($2+$3)/2); print $1, center, center + 1}}' |sort -k1,1 -k2,2n  |uniq >d9_IP_TET3.BPK25_only.summit.bed

wc -l d9_IP_TET3*summit.bed
 #   197 d9_IP_TET3.BPK25_only.summit.bed
 #   159 d9_IP_TET3.DB_common.summit.bed
 # 28334 d9_IP_TET3.DMSO_only.summit.bed

```



# reads enrichment heatmap use BPK25_vs_DMSO diff summits

##  d0_CTCF_ChIP_BPK25_vs_DMSO

### use individual bam of each repeat

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir plotHeatMap_use_CTCF_d0_BPK25_vs_DMSO_diff_summit
cd plotHeatMap_use_CTCF_d0_BPK25_vs_DMSO_diff_summit
# 创建 bam 和 bai 文件链接 
# 上传之前构造的脚本并运行
bash 6_make_H9_d0_DMSO_BPK25_CTCF_ChIP_Seq_sorted_bam_file_links.sh
ls -1
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_CTCF_rep123_together/d0_IP_CTCF.DMSO_only.summit.bed d0_IP_CTCF.DMSO_only.summit.bed
ln -s ../macs2_peak_calling_2_for_CTCF_rep123_together/d0_IP_CTCF.BPK25_only.summit.bed d0_IP_CTCF.BPK25_only.summit.bed
ln -s ../macs2_peak_calling_2_for_CTCF_rep123_together/d0_IP_CTCF.DB_common.summit.bed d0_IP_CTCF.DB_common.summit.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_CTCF_bdgdiff/d0_CTCF.ROI_DMSO.summit.sorted.bed d0_CTCF.ROI_DMSO.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_CTCF_bdgdiff/d0_CTCF.RIB_DB.summit.sorted.bed d0_CTCF.RIB_DB.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_CTCF_bdgdiff/d0_CTCF.ROI_BPK25.summit.sorted.bed d0_CTCF.ROI_BPK25.summit.sorted.bed
#
ls -1
wc -l *.bed
#  187 d0_CTCF.RIB_DB.summit.sorted.bed
#    3 d0_CTCF.ROI_BPK25.summit.sorted.bed
# 1265 d0_CTCF.ROI_DMSO.summit.sorted.bed
# 1099 d0_IP_CTCF.BPK25_only.summit.bed
#  455 d0_IP_CTCF.DB_common.summit.bed
# 1877 d0_IP_CTCF.DMSO_only.summit.bed

```

#### use bamCompare  to create log2 transformed bigwig files

```bash
# 构造 7_run_bam_comCompare_for_d0_DMSO_BPK25_CTCF_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# # DMSO 
# bamCompare -b1 hESC_H9_d0_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam -b2 hESC_H9_d0_DMSO_input_rep1.bowtie2.sorted.bam -o hESC_H9_d0_DMSO_IP_CTCF_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d0_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam -b2 hESC_H9_d0_DMSO_input_rep2.bowtie2.sorted.bam -o hESC_H9_d0_DMSO_IP_CTCF_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d0_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam -b2 hESC_H9_d0_DMSO_input_rep3.bowtie2.sorted.bam -o hESC_H9_d0_DMSO_IP_CTCF_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# # BPK25
# bamCompare -b1 hESC_H9_d0_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam -b2 hESC_H9_d0_BPK25_input_rep1.bowtie2.sorted.bam -o hESC_H9_d0_BPK25_IP_CTCF_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d0_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam -b2 hESC_H9_d0_BPK25_input_rep2.bowtie2.sorted.bam -o hESC_H9_d0_BPK25_IP_CTCF_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d0_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam -b2 hESC_H9_d0_BPK25_input_rep3.bowtie2.sorted.bam -o hESC_H9_d0_BPK25_IP_CTCF_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
nohup bash 7_run_bam_comCompare_for_d0_DMSO_BPK25_CTCF_each_repeat.sh >7_run_bam_comCompare_for_d0_DMSO_BPK25_CTCF_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bam_comCompare_for_d0_DMSO_BPK25_CTCF_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d0_DMSO_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_CTCF_rep3_TvsC.bigwig -R d0_IP_CTCF.DMSO_only.summit.bed d0_IP_CTCF.DB_common.summit.bed d0_IP_CTCF.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_CTCF_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d0_DMSO_BPK25_CTCF_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_CTCF_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d0_CTCF_DMSO_only(1877)" "d0_CTCF_DB_common(455)" "d0_CTCF_BPK25_only(1099)" --samplesLabel "d0_DMSO_CTCF_rep1_TvsC" "d0_DMSO_CTCF_rep2_TvsC" "d0_DMSO_CTCF_rep3_TvsC" "d0_BPK25_CTCF_rep1_TvsC" "d0_BPK25_CTCF_rep2_TvsC" "d0_BPK25_CTCF_rep3_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d0_DMSO_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_CTCF_rep3_TvsC.bigwig -R d0_CTCF.ROI_DMSO.summit.sorted.bed d0_CTCF.RIB_DB.summit.sorted.bed d0_CTCF.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_CTCF_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d0_DMSO_BPK25_CTCF_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_CTCF_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d0_CTCF_ROI_DMSO(1265)" "d0_CTCF_RIB_DB(187)" "d0_CTCF_ROI_BPK25(3)" --samplesLabel "d0_DMSO_CTCF_rep1_TvsC" "d0_DMSO_CTCF_rep2_TvsC" "d0_DMSO_CTCF_rep3_TvsC" "d0_BPK25_CTCF_rep1_TvsC" "d0_BPK25_CTCF_rep2_TvsC" "d0_BPK25_CTCF_rep3_TvsC"

```

### use rep123 bam together in one

#### make rep123 merge bam file links

```bash
# 上传之前构造好的脚本并运行
bash 7_make_H9_d049_DMSO_BPK25_CTCF_ChIP_Seq_rep123_merged_bam_file_links.sh
rm hESC_H9_d[49]*rep123*.bam
rm hESC_H9_d[49]*rep123*.bai

```

#### make log2 trans bigwig by bamCompare

```bash
ls -1 *_rep123.merged.sorted.bam
# hESC_H9_d0_BPK25_input_rep123.merged.sorted.bam
# hESC_H9_d0_BPK25_IP_CTCF_rep123.merged.sorted.bam
# hESC_H9_d0_DMSO_input_rep123.merged.sorted.bam
# hESC_H9_d0_DMSO_IP_CTCF_rep123.merged.sorted.bam
bamCompare -b1 hESC_H9_d0_BPK25_IP_CTCF_rep123.merged.sorted.bam -b2 hESC_H9_d0_BPK25_input_rep123.merged.sorted.bam -o H9_d0_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
bamCompare -b1 hESC_H9_d0_DMSO_IP_CTCF_rep123.merged.sorted.bam -b2 hESC_H9_d0_DMSO_input_rep123.merged.sorted.bam -o H9_d0_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S H9_d0_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig -R d0_IP_CTCF.DMSO_only.summit.bed d0_IP_CTCF.DB_common.summit.bed d0_IP_CTCF.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_CTCF_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d0_DMSO_BPK25_CTCF_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_CTCF_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d0_CTCF_DMSO_only(1877)" "d0_CTCF_DB_common(455)" "d0_CTCF_BPK25_only(1099)" --samplesLabel "d0_DMSO_CTCF_rep123_TvsC" "d0_BPK25_CTCF_rep123_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S H9_d0_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig -R d0_CTCF.ROI_DMSO.summit.sorted.bed d0_CTCF.RIB_DB.summit.sorted.bed d0_CTCF.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_CTCF_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d0_DMSO_BPK25_CTCF_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_CTCF_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d0_CTCF_ROI_DMSO(1265)" "d0_CTCF_RIB_DB(187)" "d0_CTCF_ROI_BPK25(3)" --samplesLabel "d0_DMSO_CTCF_rep123_TvsC" "d0_BPK25_CTCF_rep123_TvsC"

```



##  d4_CTCF_ChIP_BPK25_vs_DMSO

### use individual bam of each repeat

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir plotHeatMap_use_CTCF_d4_BPK25_vs_DMSO_diff_summit
cd plotHeatMap_use_CTCF_d4_BPK25_vs_DMSO_diff_summit
# 创建 bam 和 bai 文件链接 
# 上传之前构造的脚本并运行
bash 6_make_H9_d4_DMSO_BPK25_CTCF_ChIP_Seq_sorted_bam_file_links.sh
ls -1
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_CTCF_rep123_together/d4_IP_CTCF.DMSO_only.summit.bed d4_IP_CTCF.DMSO_only.summit.bed
ln -s ../macs2_peak_calling_2_for_CTCF_rep123_together/d4_IP_CTCF.BPK25_only.summit.bed d4_IP_CTCF.BPK25_only.summit.bed
ln -s ../macs2_peak_calling_2_for_CTCF_rep123_together/d4_IP_CTCF.DB_common.summit.bed d4_IP_CTCF.DB_common.summit.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_CTCF_bdgdiff/d4_CTCF.ROI_DMSO.summit.sorted.bed d4_CTCF.ROI_DMSO.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_CTCF_bdgdiff/d4_CTCF.RIB_DB.summit.sorted.bed d4_CTCF.RIB_DB.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_CTCF_bdgdiff/d4_CTCF.ROI_BPK25.summit.sorted.bed d4_CTCF.ROI_BPK25.summit.sorted.bed
#
ls -1
wc -l *.bed
#  3211 d4_CTCF.RIB_DB.summit.sorted.bed
#    37 d4_CTCF.ROI_BPK25.summit.sorted.bed
#  1456 d4_CTCF.ROI_DMSO.summit.sorted.bed
#  9246 d4_IP_CTCF.BPK25_only.summit.bed
#  8618 d4_IP_CTCF.DB_common.summit.bed
# 14755 d4_IP_CTCF.DMSO_only.summit.bed

```

#### use bamCompare  to create log2 transformed bigwig files

```bash
# 构造 7_run_bam_comCompare_for_d4_DMSO_BPK25_CTCF_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# # DMSO 
# bamCompare -b1 hESC_H9_d4_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam -b2 hESC_H9_d4_DMSO_input_rep1.bowtie2.sorted.bam -o hESC_H9_d4_DMSO_IP_CTCF_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d4_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam -b2 hESC_H9_d4_DMSO_input_rep2.bowtie2.sorted.bam -o hESC_H9_d4_DMSO_IP_CTCF_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d4_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam -b2 hESC_H9_d4_DMSO_input_rep3.bowtie2.sorted.bam -o hESC_H9_d4_DMSO_IP_CTCF_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# # BPK25
# bamCompare -b1 hESC_H9_d4_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam -b2 hESC_H9_d4_BPK25_input_rep1.bowtie2.sorted.bam -o hESC_H9_d4_BPK25_IP_CTCF_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d4_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam -b2 hESC_H9_d4_BPK25_input_rep2.bowtie2.sorted.bam -o hESC_H9_d4_BPK25_IP_CTCF_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d4_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam -b2 hESC_H9_d4_BPK25_input_rep3.bowtie2.sorted.bam -o hESC_H9_d4_BPK25_IP_CTCF_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
nohup bash 7_run_bam_comCompare_for_d4_DMSO_BPK25_CTCF_each_repeat.sh >7_run_bam_comCompare_for_d4_DMSO_BPK25_CTCF_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bam_comCompare_for_d4_DMSO_BPK25_CTCF_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d4_DMSO_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_CTCF_rep3_TvsC.bigwig -R d4_IP_CTCF.DMSO_only.summit.bed d4_IP_CTCF.DB_common.summit.bed d4_IP_CTCF.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_CTCF_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d4_DMSO_BPK25_CTCF_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_CTCF_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d4_CTCF_DMSO_only(14755)" "d4_CTCF_DB_common(8618)" "d4_CTCF_BPK25_only(9246)" --samplesLabel "d4_DMSO_CTCF_rep1_TvsC" "d4_DMSO_CTCF_rep2_TvsC" "d4_DMSO_CTCF_rep3_TvsC" "d4_BPK25_CTCF_rep1_TvsC" "d4_BPK25_CTCF_rep2_TvsC" "d4_BPK25_CTCF_rep3_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d4_DMSO_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_CTCF_rep3_TvsC.bigwig -R d4_CTCF.ROI_DMSO.summit.sorted.bed d4_CTCF.RIB_DB.summit.sorted.bed d4_CTCF.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_CTCF_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d4_DMSO_BPK25_CTCF_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_CTCF_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d4_CTCF_ROI_DMSO(1456)" "d4_CTCF_RIB_DB(3211)" "d4_CTCF_ROI_BPK25(37)" --samplesLabel "d4_DMSO_CTCF_rep1_TvsC" "d4_DMSO_CTCF_rep2_TvsC" "d4_DMSO_CTCF_rep3_TvsC" "d4_BPK25_CTCF_rep1_TvsC" "d4_BPK25_CTCF_rep2_TvsC" "d4_BPK25_CTCF_rep3_TvsC"

```

### use rep123 bam together in one

#### make rep123 merge bam file links

```bash
# 上传之前构造好的脚本并运行
bash 7_make_H9_d049_DMSO_BPK25_CTCF_ChIP_Seq_rep123_merged_bam_file_links.sh
rm hESC_H9_d[09]*rep123*.bam
rm hESC_H9_d[09]*rep123*.bai

```

#### make log2 trans bigwig by bamCompare

```bash
ls -1 *_rep123.merged.sorted.bam
# hESC_H9_d4_BPK25_input_rep123.merged.sorted.bam
# hESC_H9_d4_BPK25_IP_CTCF_rep123.merged.sorted.bam
# hESC_H9_d4_DMSO_input_rep123.merged.sorted.bam
# hESC_H9_d4_DMSO_IP_CTCF_rep123.merged.sorted.bam
bamCompare -b1 hESC_H9_d4_BPK25_IP_CTCF_rep123.merged.sorted.bam -b2 hESC_H9_d4_BPK25_input_rep123.merged.sorted.bam -o H9_d4_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
bamCompare -b1 hESC_H9_d4_DMSO_IP_CTCF_rep123.merged.sorted.bam -b2 hESC_H9_d4_DMSO_input_rep123.merged.sorted.bam -o H9_d4_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S H9_d4_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig -R d4_IP_CTCF.DMSO_only.summit.bed d4_IP_CTCF.DB_common.summit.bed d4_IP_CTCF.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_CTCF_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d4_DMSO_BPK25_CTCF_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_CTCF_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d4_CTCF_DMSO_only(14755)" "d4_CTCF_DB_common(8618)" "d4_CTCF_BPK25_only(9246)" --samplesLabel "d4_DMSO_CTCF_rep123_TvsC" "d4_BPK25_CTCF_rep123_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S H9_d4_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig -R d4_CTCF.ROI_DMSO.summit.sorted.bed d4_CTCF.RIB_DB.summit.sorted.bed d4_CTCF.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_CTCF_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d4_DMSO_BPK25_CTCF_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_CTCF_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d4_CTCF_ROI_DMSO(1456)" "d4_CTCF_RIB_DB(3211)" "d4_CTCF_ROI_BPK25(37)" --samplesLabel "d4_DMSO_CTCF_rep123_TvsC" "d4_BPK25_CTCF_rep123_TvsC"

```



##  d9_CTCF_ChIP_BPK25_vs_DMSO

### use individual bam of each repeat

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir plotHeatMap_use_CTCF_d9_BPK25_vs_DMSO_diff_summit
cd plotHeatMap_use_CTCF_d9_BPK25_vs_DMSO_diff_summit
# 创建 bam 和 bai 文件链接 
# 上传之前构造的脚本并运行
bash 6_make_H9_d9_DMSO_BPK25_CTCF_ChIP_Seq_sorted_bam_file_links.sh
ls -1
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_CTCF_rep123_together/d9_IP_CTCF.DMSO_only.summit.bed d9_IP_CTCF.DMSO_only.summit.bed
ln -s ../macs2_peak_calling_2_for_CTCF_rep123_together/d9_IP_CTCF.BPK25_only.summit.bed d9_IP_CTCF.BPK25_only.summit.bed
ln -s ../macs2_peak_calling_2_for_CTCF_rep123_together/d9_IP_CTCF.DB_common.summit.bed d9_IP_CTCF.DB_common.summit.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_CTCF_bdgdiff/d9_CTCF.ROI_DMSO.summit.sorted.bed d9_CTCF.ROI_DMSO.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_CTCF_bdgdiff/d9_CTCF.RIB_DB.summit.sorted.bed d9_CTCF.RIB_DB.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_CTCF_bdgdiff/d9_CTCF.ROI_BPK25.summit.sorted.bed d9_CTCF.ROI_BPK25.summit.sorted.bed
#
ls -1
wc -l *.bed
#   92 d9_CTCF.RIB_DB.summit.sorted.bed
# 2487 d9_CTCF.ROI_BPK25.summit.sorted.bed
#  298 d9_CTCF.ROI_DMSO.summit.sorted.bed
# 3603 d9_IP_CTCF.BPK25_only.summit.bed
#  999 d9_IP_CTCF.DB_common.summit.bed
# 7848 d9_IP_CTCF.DMSO_only.summit.bed

```

#### use bamCompare  to create log2 transformed bigwig files

```bash
# 构造 7_run_bam_comCompare_for_d9_DMSO_BPK25_CTCF_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# # DMSO 
# bamCompare -b1 hESC_H9_d9_DMSO_IP_CTCF_rep1.bowtie2.sorted.bam -b2 hESC_H9_d9_DMSO_input_rep1.bowtie2.sorted.bam -o hESC_H9_d9_DMSO_IP_CTCF_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d9_DMSO_IP_CTCF_rep2.bowtie2.sorted.bam -b2 hESC_H9_d9_DMSO_input_rep2.bowtie2.sorted.bam -o hESC_H9_d9_DMSO_IP_CTCF_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d9_DMSO_IP_CTCF_rep3.bowtie2.sorted.bam -b2 hESC_H9_d9_DMSO_input_rep3.bowtie2.sorted.bam -o hESC_H9_d9_DMSO_IP_CTCF_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# # BPK25
# bamCompare -b1 hESC_H9_d9_BPK25_IP_CTCF_rep1.bowtie2.sorted.bam -b2 hESC_H9_d9_BPK25_input_rep1.bowtie2.sorted.bam -o hESC_H9_d9_BPK25_IP_CTCF_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d9_BPK25_IP_CTCF_rep2.bowtie2.sorted.bam -b2 hESC_H9_d9_BPK25_input_rep2.bowtie2.sorted.bam -o hESC_H9_d9_BPK25_IP_CTCF_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d9_BPK25_IP_CTCF_rep3.bowtie2.sorted.bam -b2 hESC_H9_d9_BPK25_input_rep3.bowtie2.sorted.bam -o hESC_H9_d9_BPK25_IP_CTCF_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
nohup bash 7_run_bam_comCompare_for_d9_DMSO_BPK25_CTCF_each_repeat.sh >7_run_bam_comCompare_for_d9_DMSO_BPK25_CTCF_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bam_comCompare_for_d9_DMSO_BPK25_CTCF_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d9_DMSO_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_CTCF_rep3_TvsC.bigwig -R d9_IP_CTCF.DMSO_only.summit.bed d9_IP_CTCF.DB_common.summit.bed d9_IP_CTCF.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_CTCF_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d9_DMSO_BPK25_CTCF_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_CTCF_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d9_CTCF_DMSO_only(7848)" "d9_CTCF_DB_common(999)" "d9_CTCF_BPK25_only(3603)" --samplesLabel "d9_DMSO_CTCF_rep1_TvsC" "d9_DMSO_CTCF_rep2_TvsC" "d9_DMSO_CTCF_rep3_TvsC" "d9_BPK25_CTCF_rep1_TvsC" "d9_BPK25_CTCF_rep2_TvsC" "d9_BPK25_CTCF_rep3_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d9_DMSO_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_CTCF_rep3_TvsC.bigwig -R d9_CTCF.ROI_DMSO.summit.sorted.bed d9_CTCF.RIB_DB.summit.sorted.bed d9_CTCF.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_CTCF_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d9_DMSO_BPK25_CTCF_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_CTCF_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d9_CTCF_ROI_DMSO(298)" "d9_CTCF_RIB_DB(92)" "d9_CTCF_ROI_BPK25(2487)" --samplesLabel "d9_DMSO_CTCF_rep1_TvsC" "d9_DMSO_CTCF_rep2_TvsC" "d9_DMSO_CTCF_rep3_TvsC" "d9_BPK25_CTCF_rep1_TvsC" "d9_BPK25_CTCF_rep2_TvsC" "d9_BPK25_CTCF_rep3_TvsC"

```

### use rep123 bam together in one

#### make rep123 merge bam file links

```bash
# 上传之前构造好的脚本并运行
bash 7_make_H9_d049_DMSO_BPK25_CTCF_ChIP_Seq_rep123_merged_bam_file_links.sh
rm hESC_H9_d[04]*rep123*.bam
rm hESC_H9_d[04]*rep123*.bai

```

#### make log2 trans bigwig by bamCompare

```bash
ls -1 *_rep123.merged.sorted.bam
# hESC_H9_d9_BPK25_input_rep123.merged.sorted.bam
# hESC_H9_d9_BPK25_IP_CTCF_rep123.merged.sorted.bam
# hESC_H9_d9_DMSO_input_rep123.merged.sorted.bam
# hESC_H9_d9_DMSO_IP_CTCF_rep123.merged.sorted.bam
bamCompare -b1 hESC_H9_d9_BPK25_IP_CTCF_rep123.merged.sorted.bam -b2 hESC_H9_d9_BPK25_input_rep123.merged.sorted.bam -o H9_d9_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
bamCompare -b1 hESC_H9_d9_DMSO_IP_CTCF_rep123.merged.sorted.bam -b2 hESC_H9_d9_DMSO_input_rep123.merged.sorted.bam -o H9_d9_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S H9_d9_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig -R d9_IP_CTCF.DMSO_only.summit.bed d9_IP_CTCF.DB_common.summit.bed d9_IP_CTCF.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_CTCF_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d9_DMSO_BPK25_CTCF_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_CTCF_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d9_CTCF_DMSO_only(7848)" "d9_CTCF_DB_common(999)" "d9_CTCF_BPK25_only(3603)" --samplesLabel "d9_DMSO_CTCF_rep123_TvsC" "d9_BPK25_CTCF_rep123_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S H9_d9_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig -R d9_CTCF.ROI_DMSO.summit.sorted.bed d9_CTCF.RIB_DB.summit.sorted.bed d9_CTCF.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_CTCF_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d9_DMSO_BPK25_CTCF_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_CTCF_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d9_CTCF_ROI_DMSO(298)" "d9_CTCF_RIB_DB(92)" "d9_CTCF_ROI_BPK25(2487)" --samplesLabel "d9_DMSO_CTCF_rep123_TvsC" "d9_BPK25_CTCF_rep123_TvsC"

```



##  d0_MBD3_ChIP_BPK25_vs_DMSO

### use individual bam of each repeat

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir plotHeatMap_use_MBD3_d0_BPK25_vs_DMSO_diff_summit
cd plotHeatMap_use_MBD3_d0_BPK25_vs_DMSO_diff_summit
# 创建 bam 和 bai 文件链接 
# 上传之前构造的脚本并运行
bash 6_make_H9_d0_DMSO_BPK25_MBD3_ChIP_Seq_sorted_bam_file_links.sh
ls -1
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_MBD3_rep123_together/d0_IP_MBD3.DMSO_only.summit.bed d0_IP_MBD3.DMSO_only.summit.bed
ln -s ../macs2_peak_calling_2_for_MBD3_rep123_together/d0_IP_MBD3.BPK25_only.summit.bed d0_IP_MBD3.BPK25_only.summit.bed
ln -s ../macs2_peak_calling_2_for_MBD3_rep123_together/d0_IP_MBD3.DB_common.summit.bed d0_IP_MBD3.DB_common.summit.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_MBD3_bdgdiff/d0_MBD3.ROI_DMSO.summit.sorted.bed d0_MBD3.ROI_DMSO.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_MBD3_bdgdiff/d0_MBD3.RIB_DB.summit.sorted.bed d0_MBD3.RIB_DB.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_MBD3_bdgdiff/d0_MBD3.ROI_BPK25.summit.sorted.bed d0_MBD3.ROI_BPK25.summit.sorted.bed
#
ls -1
wc -l *.bed
#  215 d0_IP_MBD3.BPK25_only.summit.bed
#  527 d0_IP_MBD3.DB_common.summit.bed
# 8012 d0_IP_MBD3.DMSO_only.summit.bed
#  182 d0_MBD3.RIB_DB.summit.sorted.bed
#  136 d0_MBD3.ROI_BPK25.summit.sorted.bed
# 4354 d0_MBD3.ROI_DMSO.summit.sorted.bed

```

#### use bamCompare  to create log2 transformed bigwig files

```bash
# 构造 7_run_bam_comCompare_for_d0_DMSO_BPK25_MBD3_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# # DMSO 
# bamCompare -b1 hESC_H9_d0_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam -b2 hESC_H9_d0_DMSO_input_rep1.bowtie2.sorted.bam -o hESC_H9_d0_DMSO_IP_MBD3_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d0_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam -b2 hESC_H9_d0_DMSO_input_rep2.bowtie2.sorted.bam -o hESC_H9_d0_DMSO_IP_MBD3_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d0_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam -b2 hESC_H9_d0_DMSO_input_rep3.bowtie2.sorted.bam -o hESC_H9_d0_DMSO_IP_MBD3_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# # BPK25
# bamCompare -b1 hESC_H9_d0_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam -b2 hESC_H9_d0_BPK25_input_rep1.bowtie2.sorted.bam -o hESC_H9_d0_BPK25_IP_MBD3_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d0_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam -b2 hESC_H9_d0_BPK25_input_rep2.bowtie2.sorted.bam -o hESC_H9_d0_BPK25_IP_MBD3_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d0_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam -b2 hESC_H9_d0_BPK25_input_rep3.bowtie2.sorted.bam -o hESC_H9_d0_BPK25_IP_MBD3_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
nohup bash 7_run_bam_comCompare_for_d0_DMSO_BPK25_MBD3_each_repeat.sh >7_run_bam_comCompare_for_d0_DMSO_BPK25_MBD3_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bam_comCompare_for_d0_DMSO_BPK25_MBD3_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d0_DMSO_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_MBD3_rep3_TvsC.bigwig -R d0_IP_MBD3.DMSO_only.summit.bed d0_IP_MBD3.DB_common.summit.bed d0_IP_MBD3.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_MBD3_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d0_DMSO_BPK25_MBD3_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_MBD3_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d0_MBD3_DMSO_only(8012)" "d0_MBD3_DB_common(527)" "d0_MBD3_BPK25_only(215)" --samplesLabel "d0_DMSO_MBD3_rep1_TvsC" "d0_DMSO_MBD3_rep2_TvsC" "d0_DMSO_MBD3_rep3_TvsC" "d0_BPK25_MBD3_rep1_TvsC" "d0_BPK25_MBD3_rep2_TvsC" "d0_BPK25_MBD3_rep3_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d0_DMSO_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_MBD3_rep3_TvsC.bigwig -R d0_MBD3.ROI_DMSO.summit.sorted.bed d0_MBD3.RIB_DB.summit.sorted.bed d0_MBD3.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_MBD3_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d0_DMSO_BPK25_MBD3_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_MBD3_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d0_MBD3_ROI_DMSO(4354)" "d0_MBD3_RIB_DB(182)" "d0_MBD3_ROI_BPK25(136)" --samplesLabel "d0_DMSO_MBD3_rep1_TvsC" "d0_DMSO_MBD3_rep2_TvsC" "d0_DMSO_MBD3_rep3_TvsC" "d0_BPK25_MBD3_rep1_TvsC" "d0_BPK25_MBD3_rep2_TvsC" "d0_BPK25_MBD3_rep3_TvsC"

```

### use rep123 bam together in one

#### make rep123 merge bam file links

```bash
# 上传之前构造好的脚本并运行
bash 7_make_H9_d049_DMSO_BPK25_MBD3_ChIP_Seq_rep123_merged_bam_file_links.sh
rm hESC_H9_d[49]*rep123*.bam
rm hESC_H9_d[49]*rep123*.bai

```

#### make log2 trans bigwig by bamCompare

```bash
ls -1 *_rep123.merged.sorted.bam
# hESC_H9_d0_BPK25_input_rep123.merged.sorted.bam
# hESC_H9_d0_BPK25_IP_MBD3_rep123.merged.sorted.bam
# hESC_H9_d0_DMSO_input_rep123.merged.sorted.bam
# hESC_H9_d0_DMSO_IP_MBD3_rep123.merged.sorted.bam
bamCompare -b1 hESC_H9_d0_BPK25_IP_MBD3_rep123.merged.sorted.bam -b2 hESC_H9_d0_BPK25_input_rep123.merged.sorted.bam -o H9_d0_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
bamCompare -b1 hESC_H9_d0_DMSO_IP_MBD3_rep123.merged.sorted.bam -b2 hESC_H9_d0_DMSO_input_rep123.merged.sorted.bam -o H9_d0_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S H9_d0_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig -R d0_IP_MBD3.DMSO_only.summit.bed d0_IP_MBD3.DB_common.summit.bed d0_IP_MBD3.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_MBD3_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d0_DMSO_BPK25_MBD3_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_MBD3_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d0_MBD3_DMSO_only(8012)" "d0_MBD3_DB_common(527)" "d0_MBD3_BPK25_only(215)" --samplesLabel "d0_DMSO_MBD3_rep123_TvsC" "d0_BPK25_MBD3_rep123_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S H9_d0_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig -R d0_MBD3.ROI_DMSO.summit.sorted.bed d0_MBD3.RIB_DB.summit.sorted.bed d0_MBD3.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_MBD3_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d0_DMSO_BPK25_MBD3_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_MBD3_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d0_MBD3_ROI_DMSO(4354)" "d0_MBD3_RIB_DB(182)" "d0_MBD3_ROI_BPK25(136)" --samplesLabel "d0_DMSO_MBD3_rep123_TvsC" "d0_BPK25_MBD3_rep123_TvsC"

```



##  d4_MBD3_ChIP_BPK25_vs_DMSO

### use individual bam of each repeat

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir plotHeatMap_use_MBD3_d4_BPK25_vs_DMSO_diff_summit
cd plotHeatMap_use_MBD3_d4_BPK25_vs_DMSO_diff_summit
# 创建 bam 和 bai 文件链接 
# 上传之前构造的脚本并运行
bash 6_make_H9_d4_DMSO_BPK25_MBD3_ChIP_Seq_sorted_bam_file_links.sh
ls -1
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_MBD3_rep123_together/d4_IP_MBD3.DMSO_only.summit.bed d4_IP_MBD3.DMSO_only.summit.bed
ln -s ../macs2_peak_calling_2_for_MBD3_rep123_together/d4_IP_MBD3.BPK25_only.summit.bed d4_IP_MBD3.BPK25_only.summit.bed
ln -s ../macs2_peak_calling_2_for_MBD3_rep123_together/d4_IP_MBD3.DB_common.summit.bed d4_IP_MBD3.DB_common.summit.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_MBD3_bdgdiff/d4_MBD3.ROI_DMSO.summit.sorted.bed d4_MBD3.ROI_DMSO.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_MBD3_bdgdiff/d4_MBD3.RIB_DB.summit.sorted.bed d4_MBD3.RIB_DB.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_MBD3_bdgdiff/d4_MBD3.ROI_BPK25.summit.sorted.bed d4_MBD3.ROI_BPK25.summit.sorted.bed
#
ls -1
wc -l *.bed
# 12591 d4_IP_MBD3.BPK25_only.summit.bed
#  5401 d4_IP_MBD3.DB_common.summit.bed
# 43909 d4_IP_MBD3.DMSO_only.summit.bed
#  1405 d4_MBD3.RIB_DB.summit.sorted.bed
#   281 d4_MBD3.ROI_BPK25.summit.sorted.bed
#  5541 d4_MBD3.ROI_DMSO.summit.sorted.bed

```

#### use bamCompare  to create log2 transformed bigwig files

```bash
# 构造 7_run_bam_comCompare_for_d4_DMSO_BPK25_MBD3_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# # DMSO 
# bamCompare -b1 hESC_H9_d4_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam -b2 hESC_H9_d4_DMSO_input_rep1.bowtie2.sorted.bam -o hESC_H9_d4_DMSO_IP_MBD3_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d4_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam -b2 hESC_H9_d4_DMSO_input_rep2.bowtie2.sorted.bam -o hESC_H9_d4_DMSO_IP_MBD3_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d4_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam -b2 hESC_H9_d4_DMSO_input_rep3.bowtie2.sorted.bam -o hESC_H9_d4_DMSO_IP_MBD3_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# # BPK25
# bamCompare -b1 hESC_H9_d4_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam -b2 hESC_H9_d4_BPK25_input_rep1.bowtie2.sorted.bam -o hESC_H9_d4_BPK25_IP_MBD3_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d4_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam -b2 hESC_H9_d4_BPK25_input_rep2.bowtie2.sorted.bam -o hESC_H9_d4_BPK25_IP_MBD3_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d4_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam -b2 hESC_H9_d4_BPK25_input_rep3.bowtie2.sorted.bam -o hESC_H9_d4_BPK25_IP_MBD3_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
nohup bash 7_run_bam_comCompare_for_d4_DMSO_BPK25_MBD3_each_repeat.sh >7_run_bam_comCompare_for_d4_DMSO_BPK25_MBD3_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bam_comCompare_for_d4_DMSO_BPK25_MBD3_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d4_DMSO_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_MBD3_rep3_TvsC.bigwig -R d4_IP_MBD3.DMSO_only.summit.bed d4_IP_MBD3.DB_common.summit.bed d4_IP_MBD3.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_MBD3_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d4_DMSO_BPK25_MBD3_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_MBD3_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d4_MBD3_DMSO_only(43909)" "d4_MBD3_DB_common(5401)" "d4_MBD3_BPK25_only(12591)" --samplesLabel "d4_DMSO_MBD3_rep1_TvsC" "d4_DMSO_MBD3_rep2_TvsC" "d4_DMSO_MBD3_rep3_TvsC" "d4_BPK25_MBD3_rep1_TvsC" "d4_BPK25_MBD3_rep2_TvsC" "d4_BPK25_MBD3_rep3_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d4_DMSO_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_MBD3_rep3_TvsC.bigwig -R d4_MBD3.ROI_DMSO.summit.sorted.bed d4_MBD3.RIB_DB.summit.sorted.bed d4_MBD3.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_MBD3_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d4_DMSO_BPK25_MBD3_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_MBD3_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d4_MBD3_ROI_DMSO(5541)" "d4_MBD3_RIB_DB(1405)" "d4_MBD3_ROI_BPK25(281)" --samplesLabel "d4_DMSO_MBD3_rep1_TvsC" "d4_DMSO_MBD3_rep2_TvsC" "d4_DMSO_MBD3_rep3_TvsC" "d4_BPK25_MBD3_rep1_TvsC" "d4_BPK25_MBD3_rep2_TvsC" "d4_BPK25_MBD3_rep3_TvsC"

```

### use rep123 bam together in one

#### make rep123 merge bam file links

```bash
# 上传之前构造好的脚本并运行
bash 7_make_H9_d049_DMSO_BPK25_MBD3_ChIP_Seq_rep123_merged_bam_file_links.sh
rm hESC_H9_d[09]*rep123*.bam
rm hESC_H9_d[09]*rep123*.bai

```

#### make log2 trans bigwig by bamCompare

```bash
ls -1 *_rep123.merged.sorted.bam
# hESC_H9_d4_BPK25_input_rep123.merged.sorted.bam
# hESC_H9_d4_BPK25_IP_MBD3_rep123.merged.sorted.bam
# hESC_H9_d4_DMSO_input_rep123.merged.sorted.bam
# hESC_H9_d4_DMSO_IP_MBD3_rep123.merged.sorted.bam
bamCompare -b1 hESC_H9_d4_BPK25_IP_MBD3_rep123.merged.sorted.bam -b2 hESC_H9_d4_BPK25_input_rep123.merged.sorted.bam -o H9_d4_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
bamCompare -b1 hESC_H9_d4_DMSO_IP_MBD3_rep123.merged.sorted.bam -b2 hESC_H9_d4_DMSO_input_rep123.merged.sorted.bam -o H9_d4_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S H9_d4_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig -R d4_IP_MBD3.DMSO_only.summit.bed d4_IP_MBD3.DB_common.summit.bed d4_IP_MBD3.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_MBD3_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d4_DMSO_BPK25_MBD3_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_MBD3_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d4_MBD3_DMSO_only(43909)" "d4_MBD3_DB_common(5401)" "d4_MBD3_BPK25_only(12591)" --samplesLabel "d4_DMSO_MBD3_rep123_TvsC" "d4_BPK25_MBD3_rep123_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S H9_d4_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig -R d4_MBD3.ROI_DMSO.summit.sorted.bed d4_MBD3.RIB_DB.summit.sorted.bed d4_MBD3.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_MBD3_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d4_DMSO_BPK25_MBD3_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_MBD3_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d4_MBD3_ROI_DMSO(5541)" "d4_MBD3_RIB_DB(1405)" "d4_MBD3_ROI_BPK25(281)" --samplesLabel "d4_DMSO_MBD3_rep123_TvsC" "d4_BPK25_MBD3_rep123_TvsC"

```



##  d9_MBD3_ChIP_BPK25_vs_DMSO

### use individual bam of each repeat

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir plotHeatMap_use_MBD3_d9_BPK25_vs_DMSO_diff_summit
cd plotHeatMap_use_MBD3_d9_BPK25_vs_DMSO_diff_summit
# 创建 bam 和 bai 文件链接 
# 上传之前构造的脚本并运行
bash 6_make_H9_d9_DMSO_BPK25_MBD3_ChIP_Seq_sorted_bam_file_links.sh
ls -1
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_MBD3_rep123_together/d9_IP_MBD3.DMSO_only.summit.bed d9_IP_MBD3.DMSO_only.summit.bed
ln -s ../macs2_peak_calling_2_for_MBD3_rep123_together/d9_IP_MBD3.BPK25_only.summit.bed d9_IP_MBD3.BPK25_only.summit.bed
ln -s ../macs2_peak_calling_2_for_MBD3_rep123_together/d9_IP_MBD3.DB_common.summit.bed d9_IP_MBD3.DB_common.summit.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_MBD3_bdgdiff/d9_MBD3.ROI_DMSO.summit.sorted.bed d9_MBD3.ROI_DMSO.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_MBD3_bdgdiff/d9_MBD3.RIB_DB.summit.sorted.bed d9_MBD3.RIB_DB.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_MBD3_bdgdiff/d9_MBD3.ROI_BPK25.summit.sorted.bed d9_MBD3.ROI_BPK25.summit.sorted.bed
#
ls -1
wc -l *.bed
# 43985 d9_IP_MBD3.BPK25_only.summit.bed
#  4208 d9_IP_MBD3.DB_common.summit.bed
# 42254 d9_IP_MBD3.DMSO_only.summit.bed
#   906 d9_MBD3.RIB_DB.summit.sorted.bed
#   798 d9_MBD3.ROI_BPK25.summit.sorted.bed
#  3585 d9_MBD3.ROI_DMSO.summit.sorted.bed

```

#### use bamCompare  to create log2 transformed bigwig files

```bash
# 构造 7_run_bam_comCompare_for_d9_DMSO_BPK25_MBD3_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# # DMSO 
# bamCompare -b1 hESC_H9_d9_DMSO_IP_MBD3_rep1.bowtie2.sorted.bam -b2 hESC_H9_d9_DMSO_input_rep1.bowtie2.sorted.bam -o hESC_H9_d9_DMSO_IP_MBD3_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d9_DMSO_IP_MBD3_rep2.bowtie2.sorted.bam -b2 hESC_H9_d9_DMSO_input_rep2.bowtie2.sorted.bam -o hESC_H9_d9_DMSO_IP_MBD3_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d9_DMSO_IP_MBD3_rep3.bowtie2.sorted.bam -b2 hESC_H9_d9_DMSO_input_rep3.bowtie2.sorted.bam -o hESC_H9_d9_DMSO_IP_MBD3_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# # BPK25
# bamCompare -b1 hESC_H9_d9_BPK25_IP_MBD3_rep1.bowtie2.sorted.bam -b2 hESC_H9_d9_BPK25_input_rep1.bowtie2.sorted.bam -o hESC_H9_d9_BPK25_IP_MBD3_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d9_BPK25_IP_MBD3_rep2.bowtie2.sorted.bam -b2 hESC_H9_d9_BPK25_input_rep2.bowtie2.sorted.bam -o hESC_H9_d9_BPK25_IP_MBD3_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d9_BPK25_IP_MBD3_rep3.bowtie2.sorted.bam -b2 hESC_H9_d9_BPK25_input_rep3.bowtie2.sorted.bam -o hESC_H9_d9_BPK25_IP_MBD3_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
nohup bash 7_run_bam_comCompare_for_d9_DMSO_BPK25_MBD3_each_repeat.sh >7_run_bam_comCompare_for_d9_DMSO_BPK25_MBD3_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bam_comCompare_for_d9_DMSO_BPK25_MBD3_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d9_DMSO_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_MBD3_rep3_TvsC.bigwig -R d9_IP_MBD3.DMSO_only.summit.bed d9_IP_MBD3.DB_common.summit.bed d9_IP_MBD3.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_MBD3_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d9_DMSO_BPK25_MBD3_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_MBD3_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d9_MBD3_DMSO_only(42254)" "d9_MBD3_DB_common(4208)" "d9_MBD3_BPK25_only(43985)" --samplesLabel "d9_DMSO_MBD3_rep1_TvsC" "d9_DMSO_MBD3_rep2_TvsC" "d9_DMSO_MBD3_rep3_TvsC" "d9_BPK25_MBD3_rep1_TvsC" "d9_BPK25_MBD3_rep2_TvsC" "d9_BPK25_MBD3_rep3_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d9_DMSO_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_MBD3_rep3_TvsC.bigwig -R d9_MBD3.ROI_DMSO.summit.sorted.bed d9_MBD3.RIB_DB.summit.sorted.bed d9_MBD3.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_MBD3_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d9_DMSO_BPK25_MBD3_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_MBD3_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d9_MBD3_ROI_DMSO(3585)" "d9_MBD3_RIB_DB(798)" "d9_MBD3_ROI_BPK25(906)" --samplesLabel "d9_DMSO_MBD3_rep1_TvsC" "d9_DMSO_MBD3_rep2_TvsC" "d9_DMSO_MBD3_rep3_TvsC" "d9_BPK25_MBD3_rep1_TvsC" "d9_BPK25_MBD3_rep2_TvsC" "d9_BPK25_MBD3_rep3_TvsC"

```

### use rep123 bam together in one

#### make rep123 merge bam file links

```bash
# 上传之前构造好的脚本并运行
bash 7_make_H9_d049_DMSO_BPK25_MBD3_ChIP_Seq_rep123_merged_bam_file_links.sh
rm hESC_H9_d[04]*rep123*.bam
rm hESC_H9_d[04]*rep123*.bai

```

#### make log2 trans bigwig by bamCompare

```bash
ls -1 *_rep123.merged.sorted.bam
# hESC_H9_d9_BPK25_input_rep123.merged.sorted.bam
# hESC_H9_d9_BPK25_IP_MBD3_rep123.merged.sorted.bam
# hESC_H9_d9_DMSO_input_rep123.merged.sorted.bam
# hESC_H9_d9_DMSO_IP_MBD3_rep123.merged.sorted.bam
bamCompare -b1 hESC_H9_d9_BPK25_IP_MBD3_rep123.merged.sorted.bam -b2 hESC_H9_d9_BPK25_input_rep123.merged.sorted.bam -o H9_d9_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
bamCompare -b1 hESC_H9_d9_DMSO_IP_MBD3_rep123.merged.sorted.bam -b2 hESC_H9_d9_DMSO_input_rep123.merged.sorted.bam -o H9_d9_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S H9_d9_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig -R d9_IP_MBD3.DMSO_only.summit.bed d9_IP_MBD3.DB_common.summit.bed d9_IP_MBD3.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_MBD3_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d9_DMSO_BPK25_MBD3_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_MBD3_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d9_MBD3_DMSO_only(42254)" "d9_MBD3_DB_common(4208)" "d9_MBD3_BPK25_only(43985)" --samplesLabel "d9_DMSO_MBD3_rep123_TvsC" "d9_BPK25_MBD3_rep123_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S H9_d9_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig -R d9_MBD3.ROI_DMSO.summit.sorted.bed d9_MBD3.RIB_DB.summit.sorted.bed d9_MBD3.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_MBD3_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d9_DMSO_BPK25_MBD3_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_MBD3_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d9_MBD3_ROI_DMSO(3585)" "d9_MBD3_RIB_DB(798)" "d9_MBD3_ROI_BPK25(906)" --samplesLabel "d9_DMSO_MBD3_rep123_TvsC" "d9_BPK25_MBD3_rep123_TvsC"

```

##  d0_TET1_ChIP_BPK25_vs_DMSO

### use individual bam of each repeat

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir plotHeatMap_use_TET1_d0_BPK25_vs_DMSO_diff_summit
cd plotHeatMap_use_TET1_d0_BPK25_vs_DMSO_diff_summit
# 创建 bam 和 bai 文件链接 
# 上传之前构造的脚本并运行
bash 6_make_H9_d0_DMSO_BPK25_TET1_ChIP_Seq_sorted_bam_file_links.sh
ls -1
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_TET1_rep123_together/d0_IP_TET1.DMSO_only.summit.bed d0_IP_TET1.DMSO_only.summit.bed
ln -s ../macs2_peak_calling_2_for_TET1_rep123_together/d0_IP_TET1.BPK25_only.summit.bed d0_IP_TET1.BPK25_only.summit.bed
ln -s ../macs2_peak_calling_2_for_TET1_rep123_together/d0_IP_TET1.DB_common.summit.bed d0_IP_TET1.DB_common.summit.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_TET1_bdgdiff/d0_TET1.ROI_DMSO.summit.sorted.bed d0_TET1.ROI_DMSO.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_TET1_bdgdiff/d0_TET1.RIB_DB.summit.sorted.bed d0_TET1.RIB_DB.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_TET1_bdgdiff/d0_TET1.ROI_BPK25.summit.sorted.bed d0_TET1.ROI_BPK25.summit.sorted.bed
#
ls -1
wc -l *.bed
# 73158 d0_IP_TET1.BPK25_only.summit.bed
#  4677 d0_IP_TET1.DB_common.summit.bed
#  5368 d0_IP_TET1.DMSO_only.summit.bed
#  2336 d0_TET1.RIB_DB.summit.sorted.bed
#  3175 d0_TET1.ROI_BPK25.summit.sorted.bed
#   256 d0_TET1.ROI_DMSO.summit.sorted.bed

```

#### use bamCompare  to create log2 transformed bigwig files

```bash
# 构造 7_run_bam_comCompare_for_d0_DMSO_BPK25_TET1_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# # DMSO 
# bamCompare -b1 hESC_H9_d0_DMSO_IP_TET1_rep1.bowtie2.sorted.bam -b2 hESC_H9_d0_DMSO_input_rep1.bowtie2.sorted.bam -o hESC_H9_d0_DMSO_IP_TET1_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d0_DMSO_IP_TET1_rep2.bowtie2.sorted.bam -b2 hESC_H9_d0_DMSO_input_rep2.bowtie2.sorted.bam -o hESC_H9_d0_DMSO_IP_TET1_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d0_DMSO_IP_TET1_rep3.bowtie2.sorted.bam -b2 hESC_H9_d0_DMSO_input_rep3.bowtie2.sorted.bam -o hESC_H9_d0_DMSO_IP_TET1_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# # BPK25
# bamCompare -b1 hESC_H9_d0_BPK25_IP_TET1_rep1.bowtie2.sorted.bam -b2 hESC_H9_d0_BPK25_input_rep1.bowtie2.sorted.bam -o hESC_H9_d0_BPK25_IP_TET1_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d0_BPK25_IP_TET1_rep2.bowtie2.sorted.bam -b2 hESC_H9_d0_BPK25_input_rep2.bowtie2.sorted.bam -o hESC_H9_d0_BPK25_IP_TET1_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d0_BPK25_IP_TET1_rep3.bowtie2.sorted.bam -b2 hESC_H9_d0_BPK25_input_rep3.bowtie2.sorted.bam -o hESC_H9_d0_BPK25_IP_TET1_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
nohup bash 7_run_bam_comCompare_for_d0_DMSO_BPK25_TET1_each_repeat.sh >7_run_bam_comCompare_for_d0_DMSO_BPK25_TET1_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bam_comCompare_for_d0_DMSO_BPK25_TET1_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d0_DMSO_IP_TET1_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET1_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET1_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET1_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET1_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET1_rep3_TvsC.bigwig -R d0_IP_TET1.DMSO_only.summit.bed d0_IP_TET1.DB_common.summit.bed d0_IP_TET1.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET1_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET1_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET1_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d0_TET1_DMSO_only(5368)" "d0_TET1_DB_common(4677)" "d0_TET1_BPK25_only(73158)" --samplesLabel "d0_DMSO_TET1_rep1_TvsC" "d0_DMSO_TET1_rep2_TvsC" "d0_DMSO_TET1_rep3_TvsC" "d0_BPK25_TET1_rep1_TvsC" "d0_BPK25_TET1_rep2_TvsC" "d0_BPK25_TET1_rep3_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d0_DMSO_IP_TET1_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET1_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET1_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET1_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET1_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET1_rep3_TvsC.bigwig -R d0_TET1.ROI_DMSO.summit.sorted.bed d0_TET1.RIB_DB.summit.sorted.bed d0_TET1.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET1_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET1_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET1_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d0_TET1_ROI_DMSO(256)" "d0_TET1_RIB_DB(2336)" "d0_TET1_ROI_BPK25(3175)" --samplesLabel "d0_DMSO_TET1_rep1_TvsC" "d0_DMSO_TET1_rep2_TvsC" "d0_DMSO_TET1_rep3_TvsC" "d0_BPK25_TET1_rep1_TvsC" "d0_BPK25_TET1_rep2_TvsC" "d0_BPK25_TET1_rep3_TvsC"

```

### use rep123 bam together in one

#### make rep123 merge bam file links

```bash
# 上传之前构造好的脚本并运行
bash 7_make_H9_d049_DMSO_BPK25_TET1_ChIP_Seq_rep123_merged_bam_file_links.sh
rm hESC_H9_d[49]*rep123*.bam
rm hESC_H9_d[49]*rep123*.bai

```

#### make log2 trans bigwig by bamCompare

```bash
ls -1 *_rep123.merged.sorted.bam
# hESC_H9_d0_BPK25_input_rep123.merged.sorted.bam
# hESC_H9_d0_BPK25_IP_TET1_rep123.merged.sorted.bam
# hESC_H9_d0_DMSO_input_rep123.merged.sorted.bam
# hESC_H9_d0_DMSO_IP_TET1_rep123.merged.sorted.bam
bamCompare -b1 hESC_H9_d0_BPK25_IP_TET1_rep123.merged.sorted.bam -b2 hESC_H9_d0_BPK25_input_rep123.merged.sorted.bam -o H9_d0_BPK25_IP_TET1_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
bamCompare -b1 hESC_H9_d0_DMSO_IP_TET1_rep123.merged.sorted.bam -b2 hESC_H9_d0_DMSO_input_rep123.merged.sorted.bam -o H9_d0_DMSO_IP_TET1_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S H9_d0_DMSO_IP_TET1_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_TET1_rep123_merged_TvsC.bigwig -R d0_IP_TET1.DMSO_only.summit.bed d0_IP_TET1.DB_common.summit.bed d0_IP_TET1.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET1_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET1_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET1_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d0_TET1_DMSO_only(5368)" "d0_TET1_DB_common(4677)" "d0_TET1_BPK25_only(73158)" --samplesLabel "d0_DMSO_TET1_rep123_TvsC" "d0_BPK25_TET1_rep123_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S H9_d0_DMSO_IP_TET1_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_TET1_rep123_merged_TvsC.bigwig -R d0_TET1.ROI_DMSO.summit.sorted.bed d0_TET1.RIB_DB.summit.sorted.bed d0_TET1.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET1_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET1_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET1_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d0_TET1_ROI_DMSO(256)" "d0_TET1_RIB_DB(2336)" "d0_TET1_ROI_BPK25(3175)" --samplesLabel "d0_DMSO_TET1_rep123_TvsC" "d0_BPK25_TET1_rep123_TvsC"

```



##  d4_TET1_ChIP_BPK25_vs_DMSO

### use individual bam of each repeat

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir plotHeatMap_use_TET1_d4_BPK25_vs_DMSO_diff_summit
cd plotHeatMap_use_TET1_d4_BPK25_vs_DMSO_diff_summit
# 创建 bam 和 bai 文件链接 
# 上传之前构造的脚本并运行
bash 6_make_H9_d4_DMSO_BPK25_TET1_ChIP_Seq_sorted_bam_file_links.sh
ls -1
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_TET1_rep123_together/d4_IP_TET1.DMSO_only.summit.bed d4_IP_TET1.DMSO_only.summit.bed
ln -s ../macs2_peak_calling_2_for_TET1_rep123_together/d4_IP_TET1.BPK25_only.summit.bed d4_IP_TET1.BPK25_only.summit.bed
ln -s ../macs2_peak_calling_2_for_TET1_rep123_together/d4_IP_TET1.DB_common.summit.bed d4_IP_TET1.DB_common.summit.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_TET1_bdgdiff/d4_TET1.ROI_DMSO.summit.sorted.bed d4_TET1.ROI_DMSO.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_TET1_bdgdiff/d4_TET1.RIB_DB.summit.sorted.bed d4_TET1.RIB_DB.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_TET1_bdgdiff/d4_TET1.ROI_BPK25.summit.sorted.bed d4_TET1.ROI_BPK25.summit.sorted.bed
#
ls -1
wc -l *.bed
# 78589 d4_IP_TET1.BPK25_only.summit.bed
# 11308 d4_IP_TET1.DB_common.summit.bed
# 36270 d4_IP_TET1.DMSO_only.summit.bed
#  2929 d4_TET1.RIB_DB.summit.sorted.bed
# 10891 d4_TET1.ROI_BPK25.summit.sorted.bed
#   676 d4_TET1.ROI_DMSO.summit.sorted.bed

```

#### use bamCompare  to create log2 transformed bigwig files

```bash
# 构造 7_run_bam_comCompare_for_d4_DMSO_BPK25_TET1_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# # DMSO 
# bamCompare -b1 hESC_H9_d4_DMSO_IP_TET1_rep1.bowtie2.sorted.bam -b2 hESC_H9_d4_DMSO_input_rep1.bowtie2.sorted.bam -o hESC_H9_d4_DMSO_IP_TET1_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d4_DMSO_IP_TET1_rep2.bowtie2.sorted.bam -b2 hESC_H9_d4_DMSO_input_rep2.bowtie2.sorted.bam -o hESC_H9_d4_DMSO_IP_TET1_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d4_DMSO_IP_TET1_rep3.bowtie2.sorted.bam -b2 hESC_H9_d4_DMSO_input_rep3.bowtie2.sorted.bam -o hESC_H9_d4_DMSO_IP_TET1_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# # BPK25
# bamCompare -b1 hESC_H9_d4_BPK25_IP_TET1_rep1.bowtie2.sorted.bam -b2 hESC_H9_d4_BPK25_input_rep1.bowtie2.sorted.bam -o hESC_H9_d4_BPK25_IP_TET1_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d4_BPK25_IP_TET1_rep2.bowtie2.sorted.bam -b2 hESC_H9_d4_BPK25_input_rep2.bowtie2.sorted.bam -o hESC_H9_d4_BPK25_IP_TET1_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d4_BPK25_IP_TET1_rep3.bowtie2.sorted.bam -b2 hESC_H9_d4_BPK25_input_rep3.bowtie2.sorted.bam -o hESC_H9_d4_BPK25_IP_TET1_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
nohup bash 7_run_bam_comCompare_for_d4_DMSO_BPK25_TET1_each_repeat.sh >7_run_bam_comCompare_for_d4_DMSO_BPK25_TET1_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bam_comCompare_for_d4_DMSO_BPK25_TET1_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d4_DMSO_IP_TET1_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET1_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET1_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET1_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET1_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET1_rep3_TvsC.bigwig -R d4_IP_TET1.DMSO_only.summit.bed d4_IP_TET1.DB_common.summit.bed d4_IP_TET1.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET1_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET1_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET1_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d4_TET1_DMSO_only(36270)" "d4_TET1_DB_common(11308)" "d4_TET1_BPK25_only(78589)" --samplesLabel "d4_DMSO_TET1_rep1_TvsC" "d4_DMSO_TET1_rep2_TvsC" "d4_DMSO_TET1_rep3_TvsC" "d4_BPK25_TET1_rep1_TvsC" "d4_BPK25_TET1_rep2_TvsC" "d4_BPK25_TET1_rep3_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d4_DMSO_IP_TET1_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET1_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET1_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET1_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET1_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET1_rep3_TvsC.bigwig -R d4_TET1.ROI_DMSO.summit.sorted.bed d4_TET1.RIB_DB.summit.sorted.bed d4_TET1.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET1_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET1_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET1_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d4_TET1_ROI_DMSO(676)" "d4_TET1_RIB_DB(2929)" "d4_TET1_ROI_BPK25(10891)" --samplesLabel "d4_DMSO_TET1_rep1_TvsC" "d4_DMSO_TET1_rep2_TvsC" "d4_DMSO_TET1_rep3_TvsC" "d4_BPK25_TET1_rep1_TvsC" "d4_BPK25_TET1_rep2_TvsC" "d4_BPK25_TET1_rep3_TvsC"

```

### use rep123 bam together in one

#### make rep123 merge bam file links

```bash
# 上传之前构造好的脚本并运行
bash 7_make_H9_d049_DMSO_BPK25_TET1_ChIP_Seq_rep123_merged_bam_file_links.sh
rm hESC_H9_d[09]*rep123*.bam
rm hESC_H9_d[09]*rep123*.bai

```

#### make log2 trans bigwig by bamCompare

```bash
ls -1 *_rep123.merged.sorted.bam
# hESC_H9_d4_BPK25_input_rep123.merged.sorted.bam
# hESC_H9_d4_BPK25_IP_TET1_rep123.merged.sorted.bam
# hESC_H9_d4_DMSO_input_rep123.merged.sorted.bam
# hESC_H9_d4_DMSO_IP_TET1_rep123.merged.sorted.bam
bamCompare -b1 hESC_H9_d4_BPK25_IP_TET1_rep123.merged.sorted.bam -b2 hESC_H9_d4_BPK25_input_rep123.merged.sorted.bam -o H9_d4_BPK25_IP_TET1_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
bamCompare -b1 hESC_H9_d4_DMSO_IP_TET1_rep123.merged.sorted.bam -b2 hESC_H9_d4_DMSO_input_rep123.merged.sorted.bam -o H9_d4_DMSO_IP_TET1_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S H9_d4_DMSO_IP_TET1_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_TET1_rep123_merged_TvsC.bigwig -R d4_IP_TET1.DMSO_only.summit.bed d4_IP_TET1.DB_common.summit.bed d4_IP_TET1.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET1_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET1_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET1_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d4_TET1_DMSO_only(36270)" "d4_TET1_DB_common(11308)" "d4_TET1_BPK25_only(78589)" --samplesLabel "d4_DMSO_TET1_rep123_TvsC" "d4_BPK25_TET1_rep123_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S H9_d4_DMSO_IP_TET1_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_TET1_rep123_merged_TvsC.bigwig -R d4_TET1.ROI_DMSO.summit.sorted.bed d4_TET1.RIB_DB.summit.sorted.bed d4_TET1.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET1_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET1_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET1_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d4_TET1_ROI_DMSO(676)" "d4_TET1_RIB_DB(2929)" "d4_TET1_ROI_BPK25(10891)" --samplesLabel "d4_DMSO_TET1_rep123_TvsC" "d4_BPK25_TET1_rep123_TvsC"

```



##  d9_TET1_ChIP_BPK25_vs_DMSO

### use individual bam of each repeat

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir plotHeatMap_use_TET1_d9_BPK25_vs_DMSO_diff_summit
cd plotHeatMap_use_TET1_d9_BPK25_vs_DMSO_diff_summit
# 创建 bam 和 bai 文件链接 
# 上传之前构造的脚本并运行
bash 6_make_H9_d9_DMSO_BPK25_TET1_ChIP_Seq_sorted_bam_file_links.sh
ls -1
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_TET1_rep123_together/d9_IP_TET1.DMSO_only.summit.bed d9_IP_TET1.DMSO_only.summit.bed
ln -s ../macs2_peak_calling_2_for_TET1_rep123_together/d9_IP_TET1.BPK25_only.summit.bed d9_IP_TET1.BPK25_only.summit.bed
ln -s ../macs2_peak_calling_2_for_TET1_rep123_together/d9_IP_TET1.DB_common.summit.bed d9_IP_TET1.DB_common.summit.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_TET1_bdgdiff/d9_TET1.ROI_DMSO.summit.sorted.bed d9_TET1.ROI_DMSO.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_TET1_bdgdiff/d9_TET1.RIB_DB.summit.sorted.bed d9_TET1.RIB_DB.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_TET1_bdgdiff/d9_TET1.ROI_BPK25.summit.sorted.bed d9_TET1.ROI_BPK25.summit.sorted.bed
#
ls -1
wc -l *.bed
# 44567 d9_IP_TET1.BPK25_only.summit.bed
#  1837 d9_IP_TET1.DB_common.summit.bed
# 24932 d9_IP_TET1.DMSO_only.summit.bed
#   132 d9_TET1.RIB_DB.summit.sorted.bed
# 10859 d9_TET1.ROI_BPK25.summit.sorted.bed
#    47 d9_TET1.ROI_DMSO.summit.sorted.bed

```

#### use bamCompare  to create log2 transformed bigwig files

```bash
# 构造 7_run_bam_comCompare_for_d9_DMSO_BPK25_TET1_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# # DMSO 
# bamCompare -b1 hESC_H9_d9_DMSO_IP_TET1_rep1.bowtie2.sorted.bam -b2 hESC_H9_d9_DMSO_input_rep1.bowtie2.sorted.bam -o hESC_H9_d9_DMSO_IP_TET1_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d9_DMSO_IP_TET1_rep2.bowtie2.sorted.bam -b2 hESC_H9_d9_DMSO_input_rep2.bowtie2.sorted.bam -o hESC_H9_d9_DMSO_IP_TET1_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d9_DMSO_IP_TET1_rep3.bowtie2.sorted.bam -b2 hESC_H9_d9_DMSO_input_rep3.bowtie2.sorted.bam -o hESC_H9_d9_DMSO_IP_TET1_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# # BPK25
# bamCompare -b1 hESC_H9_d9_BPK25_IP_TET1_rep1.bowtie2.sorted.bam -b2 hESC_H9_d9_BPK25_input_rep1.bowtie2.sorted.bam -o hESC_H9_d9_BPK25_IP_TET1_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d9_BPK25_IP_TET1_rep2.bowtie2.sorted.bam -b2 hESC_H9_d9_BPK25_input_rep2.bowtie2.sorted.bam -o hESC_H9_d9_BPK25_IP_TET1_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d9_BPK25_IP_TET1_rep3.bowtie2.sorted.bam -b2 hESC_H9_d9_BPK25_input_rep3.bowtie2.sorted.bam -o hESC_H9_d9_BPK25_IP_TET1_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
nohup bash 7_run_bam_comCompare_for_d9_DMSO_BPK25_TET1_each_repeat.sh >7_run_bam_comCompare_for_d9_DMSO_BPK25_TET1_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bam_comCompare_for_d9_DMSO_BPK25_TET1_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d9_DMSO_IP_TET1_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET1_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET1_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET1_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET1_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET1_rep3_TvsC.bigwig -R d9_IP_TET1.DMSO_only.summit.bed d9_IP_TET1.DB_common.summit.bed d9_IP_TET1.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET1_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET1_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET1_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d9_TET1_DMSO_only(24932)" "d9_TET1_DB_common(1837)" "d9_TET1_BPK25_only(44567)" --samplesLabel "d9_DMSO_TET1_rep1_TvsC" "d9_DMSO_TET1_rep2_TvsC" "d9_DMSO_TET1_rep3_TvsC" "d9_BPK25_TET1_rep1_TvsC" "d9_BPK25_TET1_rep2_TvsC" "d9_BPK25_TET1_rep3_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d9_DMSO_IP_TET1_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET1_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET1_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET1_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET1_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET1_rep3_TvsC.bigwig -R d9_TET1.ROI_DMSO.summit.sorted.bed d9_TET1.RIB_DB.summit.sorted.bed d9_TET1.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET1_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET1_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET1_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d9_TET1_ROI_DMSO(47)" "d9_TET1_RIB_DB(132)" "d9_TET1_ROI_BPK25(10859)" --samplesLabel "d9_DMSO_TET1_rep1_TvsC" "d9_DMSO_TET1_rep2_TvsC" "d9_DMSO_TET1_rep3_TvsC" "d9_BPK25_TET1_rep1_TvsC" "d9_BPK25_TET1_rep2_TvsC" "d9_BPK25_TET1_rep3_TvsC"

```

### use rep123 bam together in one

#### make rep123 merge bam file links

```bash
# 上传之前构造好的脚本并运行
bash 7_make_H9_d049_DMSO_BPK25_TET1_ChIP_Seq_rep123_merged_bam_file_links.sh
rm hESC_H9_d[04]*rep123*.bam
rm hESC_H9_d[04]*rep123*.bai

```

#### make log2 trans bigwig by bamCompare

```bash
ls -1 *_rep123.merged.sorted.bam
# hESC_H9_d9_BPK25_input_rep123.merged.sorted.bam
# hESC_H9_d9_BPK25_IP_TET1_rep123.merged.sorted.bam
# hESC_H9_d9_DMSO_input_rep123.merged.sorted.bam
# hESC_H9_d9_DMSO_IP_TET1_rep123.merged.sorted.bam
bamCompare -b1 hESC_H9_d9_BPK25_IP_TET1_rep123.merged.sorted.bam -b2 hESC_H9_d9_BPK25_input_rep123.merged.sorted.bam -o H9_d9_BPK25_IP_TET1_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
bamCompare -b1 hESC_H9_d9_DMSO_IP_TET1_rep123.merged.sorted.bam -b2 hESC_H9_d9_DMSO_input_rep123.merged.sorted.bam -o H9_d9_DMSO_IP_TET1_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S H9_d9_DMSO_IP_TET1_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_TET1_rep123_merged_TvsC.bigwig -R d9_IP_TET1.DMSO_only.summit.bed d9_IP_TET1.DB_common.summit.bed d9_IP_TET1.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET1_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET1_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET1_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d9_TET1_DMSO_only(24932)" "d9_TET1_DB_common(1837)" "d9_TET1_BPK25_only(44567)" --samplesLabel "d9_DMSO_TET1_rep123_TvsC" "d9_BPK25_TET1_rep123_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S H9_d9_DMSO_IP_TET1_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_TET1_rep123_merged_TvsC.bigwig -R d9_TET1.ROI_DMSO.summit.sorted.bed d9_TET1.RIB_DB.summit.sorted.bed d9_TET1.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET1_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET1_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET1_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d9_TET1_ROI_DMSO(47)" "d9_TET1_RIB_DB(132)" "d9_TET1_ROI_BPK25(10859)" --samplesLabel "d9_DMSO_TET1_rep123_TvsC" "d9_BPK25_TET1_rep123_TvsC"

```

##  d0_TET2_ChIP_BPK25_vs_DMSO

### use individual bam of each repeat

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir plotHeatMap_use_TET2_d0_BPK25_vs_DMSO_diff_summit
cd plotHeatMap_use_TET2_d0_BPK25_vs_DMSO_diff_summit
# 创建 bam 和 bai 文件链接 
# 上传之前构造的脚本并运行
bash 6_make_H9_d0_DMSO_BPK25_TET2_ChIP_Seq_sorted_bam_file_links.sh
ls -1
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_TET2_rep123_together/d0_IP_TET2.DMSO_only.summit.bed d0_IP_TET2.DMSO_only.summit.bed
ln -s ../macs2_peak_calling_2_for_TET2_rep123_together/d0_IP_TET2.BPK25_only.summit.bed d0_IP_TET2.BPK25_only.summit.bed
ln -s ../macs2_peak_calling_2_for_TET2_rep123_together/d0_IP_TET2.DB_common.summit.bed d0_IP_TET2.DB_common.summit.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_TET2_bdgdiff/d0_TET2.ROI_DMSO.summit.sorted.bed d0_TET2.ROI_DMSO.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_TET2_bdgdiff/d0_TET2.RIB_DB.summit.sorted.bed d0_TET2.RIB_DB.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_TET2_bdgdiff/d0_TET2.ROI_BPK25.summit.sorted.bed d0_TET2.ROI_BPK25.summit.sorted.bed
#
ls -1
wc -l *.bed
# 32674 d0_IP_TET2.BPK25_only.summit.bed
#  1793 d0_IP_TET2.DB_common.summit.bed
#  2234 d0_IP_TET2.DMSO_only.summit.bed
#  1000 d0_TET2.RIB_DB.summit.sorted.bed
#  2027 d0_TET2.ROI_BPK25.summit.sorted.bed
#   304 d0_TET2.ROI_DMSO.summit.sorted.bed

```

#### use bamCompare  to create log2 transformed bigwig files

```bash
# 构造 7_run_bam_comCompare_for_d0_DMSO_BPK25_TET2_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# # DMSO 
# bamCompare -b1 hESC_H9_d0_DMSO_IP_TET2_rep1.bowtie2.sorted.bam -b2 hESC_H9_d0_DMSO_input_rep1.bowtie2.sorted.bam -o hESC_H9_d0_DMSO_IP_TET2_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d0_DMSO_IP_TET2_rep2.bowtie2.sorted.bam -b2 hESC_H9_d0_DMSO_input_rep2.bowtie2.sorted.bam -o hESC_H9_d0_DMSO_IP_TET2_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d0_DMSO_IP_TET2_rep3.bowtie2.sorted.bam -b2 hESC_H9_d0_DMSO_input_rep3.bowtie2.sorted.bam -o hESC_H9_d0_DMSO_IP_TET2_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# # BPK25
# bamCompare -b1 hESC_H9_d0_BPK25_IP_TET2_rep1.bowtie2.sorted.bam -b2 hESC_H9_d0_BPK25_input_rep1.bowtie2.sorted.bam -o hESC_H9_d0_BPK25_IP_TET2_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d0_BPK25_IP_TET2_rep2.bowtie2.sorted.bam -b2 hESC_H9_d0_BPK25_input_rep2.bowtie2.sorted.bam -o hESC_H9_d0_BPK25_IP_TET2_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d0_BPK25_IP_TET2_rep3.bowtie2.sorted.bam -b2 hESC_H9_d0_BPK25_input_rep3.bowtie2.sorted.bam -o hESC_H9_d0_BPK25_IP_TET2_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
nohup bash 7_run_bam_comCompare_for_d0_DMSO_BPK25_TET2_each_repeat.sh >7_run_bam_comCompare_for_d0_DMSO_BPK25_TET2_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bam_comCompare_for_d0_DMSO_BPK25_TET2_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d0_DMSO_IP_TET2_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET2_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET2_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET2_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET2_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET2_rep3_TvsC.bigwig -R d0_IP_TET2.DMSO_only.summit.bed d0_IP_TET2.DB_common.summit.bed d0_IP_TET2.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET2_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET2_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET2_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d0_TET2_DMSO_only(2234)" "d0_TET2_DB_common(1793)" "d0_TET2_BPK25_only(32674)" --samplesLabel "d0_DMSO_TET2_rep1_TvsC" "d0_DMSO_TET2_rep2_TvsC" "d0_DMSO_TET2_rep3_TvsC" "d0_BPK25_TET2_rep1_TvsC" "d0_BPK25_TET2_rep2_TvsC" "d0_BPK25_TET2_rep3_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d0_DMSO_IP_TET2_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET2_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET2_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET2_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET2_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET2_rep3_TvsC.bigwig -R d0_TET2.ROI_DMSO.summit.sorted.bed d0_TET2.RIB_DB.summit.sorted.bed d0_TET2.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET2_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET2_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET2_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d0_TET2_ROI_DMSO(304)" "d0_TET2_RIB_DB(1000)" "d0_TET2_ROI_BPK25(2027)" --samplesLabel "d0_DMSO_TET2_rep1_TvsC" "d0_DMSO_TET2_rep2_TvsC" "d0_DMSO_TET2_rep3_TvsC" "d0_BPK25_TET2_rep1_TvsC" "d0_BPK25_TET2_rep2_TvsC" "d0_BPK25_TET2_rep3_TvsC"

```

### use rep123 bam together in one

#### make rep123 merge bam file links

```bash
# 上传之前构造好的脚本并运行
bash 7_make_H9_d049_DMSO_BPK25_TET2_ChIP_Seq_rep123_merged_bam_file_links.sh
rm hESC_H9_d[49]*rep123*.bam
rm hESC_H9_d[49]*rep123*.bai

```

#### make log2 trans bigwig by bamCompare

```bash
ls -1 *_rep123.merged.sorted.bam
# hESC_H9_d0_BPK25_input_rep123.merged.sorted.bam
# hESC_H9_d0_BPK25_IP_TET2_rep123.merged.sorted.bam
# hESC_H9_d0_DMSO_input_rep123.merged.sorted.bam
# hESC_H9_d0_DMSO_IP_TET2_rep123.merged.sorted.bam
bamCompare -b1 hESC_H9_d0_BPK25_IP_TET2_rep123.merged.sorted.bam -b2 hESC_H9_d0_BPK25_input_rep123.merged.sorted.bam -o H9_d0_BPK25_IP_TET2_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
bamCompare -b1 hESC_H9_d0_DMSO_IP_TET2_rep123.merged.sorted.bam -b2 hESC_H9_d0_DMSO_input_rep123.merged.sorted.bam -o H9_d0_DMSO_IP_TET2_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S H9_d0_DMSO_IP_TET2_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_TET2_rep123_merged_TvsC.bigwig -R d0_IP_TET2.DMSO_only.summit.bed d0_IP_TET2.DB_common.summit.bed d0_IP_TET2.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET2_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET2_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET2_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d0_TET2_DMSO_only(2234)" "d0_TET2_DB_common(1793)" "d0_TET2_BPK25_only(32674)" --samplesLabel "d0_DMSO_TET2_rep123_TvsC" "d0_BPK25_TET2_rep123_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S H9_d0_DMSO_IP_TET2_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_TET2_rep123_merged_TvsC.bigwig -R d0_TET2.ROI_DMSO.summit.sorted.bed d0_TET2.RIB_DB.summit.sorted.bed d0_TET2.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET2_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET2_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET2_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d0_TET2_ROI_DMSO(304)" "d0_TET2_RIB_DB(1000)" "d0_TET2_ROI_BPK25(2027)" --samplesLabel "d0_DMSO_TET2_rep123_TvsC" "d0_BPK25_TET2_rep123_TvsC"

```



##  d4_TET2_ChIP_BPK25_vs_DMSO

### use individual bam of each repeat

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir plotHeatMap_use_TET2_d4_BPK25_vs_DMSO_diff_summit
cd plotHeatMap_use_TET2_d4_BPK25_vs_DMSO_diff_summit
# 创建 bam 和 bai 文件链接 
# 上传之前构造的脚本并运行
bash 6_make_H9_d4_DMSO_BPK25_TET2_ChIP_Seq_sorted_bam_file_links.sh
ls -1
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_TET2_rep123_together/d4_IP_TET2.DMSO_only.summit.bed d4_IP_TET2.DMSO_only.summit.bed
ln -s ../macs2_peak_calling_2_for_TET2_rep123_together/d4_IP_TET2.BPK25_only.summit.bed d4_IP_TET2.BPK25_only.summit.bed
ln -s ../macs2_peak_calling_2_for_TET2_rep123_together/d4_IP_TET2.DB_common.summit.bed d4_IP_TET2.DB_common.summit.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_TET2_bdgdiff/d4_TET2.ROI_DMSO.summit.sorted.bed d4_TET2.ROI_DMSO.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_TET2_bdgdiff/d4_TET2.RIB_DB.summit.sorted.bed d4_TET2.RIB_DB.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_TET2_bdgdiff/d4_TET2.ROI_BPK25.summit.sorted.bed d4_TET2.ROI_BPK25.summit.sorted.bed
#
ls -1
wc -l *.bed
# 13208 d4_IP_TET2.BPK25_only.summit.bed
#  1487 d4_IP_TET2.DB_common.summit.bed
# 37121 d4_IP_TET2.DMSO_only.summit.bed
#   466 d4_TET2.RIB_DB.summit.sorted.bed
#  4150 d4_TET2.ROI_BPK25.summit.sorted.bed
#  1870 d4_TET2.ROI_DMSO.summit.sorted.bed

```

#### use bamCompare  to create log2 transformed bigwig files

```bash
# 构造 7_run_bam_comCompare_for_d4_DMSO_BPK25_TET2_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# # DMSO 
# bamCompare -b1 hESC_H9_d4_DMSO_IP_TET2_rep1.bowtie2.sorted.bam -b2 hESC_H9_d4_DMSO_input_rep1.bowtie2.sorted.bam -o hESC_H9_d4_DMSO_IP_TET2_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d4_DMSO_IP_TET2_rep2.bowtie2.sorted.bam -b2 hESC_H9_d4_DMSO_input_rep2.bowtie2.sorted.bam -o hESC_H9_d4_DMSO_IP_TET2_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d4_DMSO_IP_TET2_rep3.bowtie2.sorted.bam -b2 hESC_H9_d4_DMSO_input_rep3.bowtie2.sorted.bam -o hESC_H9_d4_DMSO_IP_TET2_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# # BPK25
# bamCompare -b1 hESC_H9_d4_BPK25_IP_TET2_rep1.bowtie2.sorted.bam -b2 hESC_H9_d4_BPK25_input_rep1.bowtie2.sorted.bam -o hESC_H9_d4_BPK25_IP_TET2_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d4_BPK25_IP_TET2_rep2.bowtie2.sorted.bam -b2 hESC_H9_d4_BPK25_input_rep2.bowtie2.sorted.bam -o hESC_H9_d4_BPK25_IP_TET2_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d4_BPK25_IP_TET2_rep3.bowtie2.sorted.bam -b2 hESC_H9_d4_BPK25_input_rep3.bowtie2.sorted.bam -o hESC_H9_d4_BPK25_IP_TET2_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
nohup bash 7_run_bam_comCompare_for_d4_DMSO_BPK25_TET2_each_repeat.sh >7_run_bam_comCompare_for_d4_DMSO_BPK25_TET2_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bam_comCompare_for_d4_DMSO_BPK25_TET2_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d4_DMSO_IP_TET2_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET2_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET2_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET2_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET2_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET2_rep3_TvsC.bigwig -R d4_IP_TET2.DMSO_only.summit.bed d4_IP_TET2.DB_common.summit.bed d4_IP_TET2.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET2_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET2_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET2_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d4_TET2_DMSO_only(37121)" "d4_TET2_DB_common(1487)" "d4_TET2_BPK25_only(13208)" --samplesLabel "d4_DMSO_TET2_rep1_TvsC" "d4_DMSO_TET2_rep2_TvsC" "d4_DMSO_TET2_rep3_TvsC" "d4_BPK25_TET2_rep1_TvsC" "d4_BPK25_TET2_rep2_TvsC" "d4_BPK25_TET2_rep3_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d4_DMSO_IP_TET2_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET2_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET2_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET2_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET2_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET2_rep3_TvsC.bigwig -R d4_TET2.ROI_DMSO.summit.sorted.bed d4_TET2.RIB_DB.summit.sorted.bed d4_TET2.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET2_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET2_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET2_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d4_TET2_ROI_DMSO(1870)" "d4_TET2_RIB_DB(466)" "d4_TET2_ROI_BPK25(4150)" --samplesLabel "d4_DMSO_TET2_rep1_TvsC" "d4_DMSO_TET2_rep2_TvsC" "d4_DMSO_TET2_rep3_TvsC" "d4_BPK25_TET2_rep1_TvsC" "d4_BPK25_TET2_rep2_TvsC" "d4_BPK25_TET2_rep3_TvsC"

```

### use rep123 bam together in one

#### make rep123 merge bam file links

```bash
# 上传之前构造好的脚本并运行
bash 7_make_H9_d049_DMSO_BPK25_TET2_ChIP_Seq_rep123_merged_bam_file_links.sh
rm hESC_H9_d[09]*rep123*.bam
rm hESC_H9_d[09]*rep123*.bai

```

#### make log2 trans bigwig by bamCompare

```bash
ls -1 *_rep123.merged.sorted.bam
# hESC_H9_d4_BPK25_input_rep123.merged.sorted.bam
# hESC_H9_d4_BPK25_IP_TET2_rep123.merged.sorted.bam
# hESC_H9_d4_DMSO_input_rep123.merged.sorted.bam
# hESC_H9_d4_DMSO_IP_TET2_rep123.merged.sorted.bam
bamCompare -b1 hESC_H9_d4_BPK25_IP_TET2_rep123.merged.sorted.bam -b2 hESC_H9_d4_BPK25_input_rep123.merged.sorted.bam -o H9_d4_BPK25_IP_TET2_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
bamCompare -b1 hESC_H9_d4_DMSO_IP_TET2_rep123.merged.sorted.bam -b2 hESC_H9_d4_DMSO_input_rep123.merged.sorted.bam -o H9_d4_DMSO_IP_TET2_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S H9_d4_DMSO_IP_TET2_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_TET2_rep123_merged_TvsC.bigwig -R d4_IP_TET2.DMSO_only.summit.bed d4_IP_TET2.DB_common.summit.bed d4_IP_TET2.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET2_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET2_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET2_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d4_TET2_DMSO_only(37121)" "d4_TET2_DB_common(1487)" "d4_TET2_BPK25_only(13208)" --samplesLabel "d4_DMSO_TET2_rep123_TvsC" "d4_BPK25_TET2_rep123_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S H9_d4_DMSO_IP_TET2_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_TET2_rep123_merged_TvsC.bigwig -R d4_TET2.ROI_DMSO.summit.sorted.bed d4_TET2.RIB_DB.summit.sorted.bed d4_TET2.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET2_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET2_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET2_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d4_TET2_ROI_DMSO(1870)" "d4_TET2_RIB_DB(466)" "d4_TET2_ROI_BPK25(4150)" --samplesLabel "d4_DMSO_TET2_rep123_TvsC" "d4_BPK25_TET2_rep123_TvsC"

```



##  d9_TET2_ChIP_BPK25_vs_DMSO

### use individual bam of each repeat

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir plotHeatMap_use_TET2_d9_BPK25_vs_DMSO_diff_summit
cd plotHeatMap_use_TET2_d9_BPK25_vs_DMSO_diff_summit
# 创建 bam 和 bai 文件链接 
# 上传之前构造的脚本并运行
bash 6_make_H9_d9_DMSO_BPK25_TET2_ChIP_Seq_sorted_bam_file_links.sh
ls -1
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_TET2_rep123_together/d9_IP_TET2.DMSO_only.summit.bed d9_IP_TET2.DMSO_only.summit.bed
ln -s ../macs2_peak_calling_2_for_TET2_rep123_together/d9_IP_TET2.BPK25_only.summit.bed d9_IP_TET2.BPK25_only.summit.bed
ln -s ../macs2_peak_calling_2_for_TET2_rep123_together/d9_IP_TET2.DB_common.summit.bed d9_IP_TET2.DB_common.summit.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_TET2_bdgdiff/d9_TET2.ROI_DMSO.summit.sorted.bed d9_TET2.ROI_DMSO.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_TET2_bdgdiff/d9_TET2.RIB_DB.summit.sorted.bed d9_TET2.RIB_DB.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_TET2_bdgdiff/d9_TET2.ROI_BPK25.summit.sorted.bed d9_TET2.ROI_BPK25.summit.sorted.bed
#
ls -1
wc -l *.bed
#  1540 d9_IP_TET2.BPK25_only.summit.bed
#   493 d9_IP_TET2.DB_common.summit.bed
# 15243 d9_IP_TET2.DMSO_only.summit.bed
#    58 d9_TET2.RIB_DB.summit.sorted.bed
#    13 d9_TET2.ROI_BPK25.summit.sorted.bed
#   519 d9_TET2.ROI_DMSO.summit.sorted.bed

```

#### use bamCompare  to create log2 transformed bigwig files

```bash
# 构造 7_run_bam_comCompare_for_d9_DMSO_BPK25_TET2_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# # DMSO 
# bamCompare -b1 hESC_H9_d9_DMSO_IP_TET2_rep1.bowtie2.sorted.bam -b2 hESC_H9_d9_DMSO_input_rep1.bowtie2.sorted.bam -o hESC_H9_d9_DMSO_IP_TET2_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d9_DMSO_IP_TET2_rep2.bowtie2.sorted.bam -b2 hESC_H9_d9_DMSO_input_rep2.bowtie2.sorted.bam -o hESC_H9_d9_DMSO_IP_TET2_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d9_DMSO_IP_TET2_rep3.bowtie2.sorted.bam -b2 hESC_H9_d9_DMSO_input_rep3.bowtie2.sorted.bam -o hESC_H9_d9_DMSO_IP_TET2_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# # BPK25
# bamCompare -b1 hESC_H9_d9_BPK25_IP_TET2_rep1.bowtie2.sorted.bam -b2 hESC_H9_d9_BPK25_input_rep1.bowtie2.sorted.bam -o hESC_H9_d9_BPK25_IP_TET2_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d9_BPK25_IP_TET2_rep2.bowtie2.sorted.bam -b2 hESC_H9_d9_BPK25_input_rep2.bowtie2.sorted.bam -o hESC_H9_d9_BPK25_IP_TET2_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d9_BPK25_IP_TET2_rep3.bowtie2.sorted.bam -b2 hESC_H9_d9_BPK25_input_rep3.bowtie2.sorted.bam -o hESC_H9_d9_BPK25_IP_TET2_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
nohup bash 7_run_bam_comCompare_for_d9_DMSO_BPK25_TET2_each_repeat.sh >7_run_bam_comCompare_for_d9_DMSO_BPK25_TET2_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bam_comCompare_for_d9_DMSO_BPK25_TET2_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d9_DMSO_IP_TET2_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET2_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET2_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET2_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET2_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET2_rep3_TvsC.bigwig -R d9_IP_TET2.DMSO_only.summit.bed d9_IP_TET2.DB_common.summit.bed d9_IP_TET2.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET2_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET2_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET2_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d9_TET2_DMSO_only(15243)" "d9_TET2_DB_common(493)" "d9_TET2_BPK25_only(1540)" --samplesLabel "d9_DMSO_TET2_rep1_TvsC" "d9_DMSO_TET2_rep2_TvsC" "d9_DMSO_TET2_rep3_TvsC" "d9_BPK25_TET2_rep1_TvsC" "d9_BPK25_TET2_rep2_TvsC" "d9_BPK25_TET2_rep3_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d9_DMSO_IP_TET2_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET2_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET2_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET2_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET2_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET2_rep3_TvsC.bigwig -R d9_TET2.ROI_DMSO.summit.sorted.bed d9_TET2.RIB_DB.summit.sorted.bed d9_TET2.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET2_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET2_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET2_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d9_TET2_ROI_DMSO(519)" "d9_TET2_RIB_DB(58)" "d9_TET2_ROI_BPK25(13)" --samplesLabel "d9_DMSO_TET2_rep1_TvsC" "d9_DMSO_TET2_rep2_TvsC" "d9_DMSO_TET2_rep3_TvsC" "d9_BPK25_TET2_rep1_TvsC" "d9_BPK25_TET2_rep2_TvsC" "d9_BPK25_TET2_rep3_TvsC"

```

### use rep123 bam together in one

#### make rep123 merge bam file links

```bash
# 上传之前构造好的脚本并运行
bash 7_make_H9_d049_DMSO_BPK25_TET2_ChIP_Seq_rep123_merged_bam_file_links.sh
rm hESC_H9_d[04]*rep123*.bam
rm hESC_H9_d[04]*rep123*.bai

```

#### make log2 trans bigwig by bamCompare

```bash
ls -1 *_rep123.merged.sorted.bam
# hESC_H9_d9_BPK25_input_rep123.merged.sorted.bam
# hESC_H9_d9_BPK25_IP_TET2_rep123.merged.sorted.bam
# hESC_H9_d9_DMSO_input_rep123.merged.sorted.bam
# hESC_H9_d9_DMSO_IP_TET2_rep123.merged.sorted.bam
bamCompare -b1 hESC_H9_d9_BPK25_IP_TET2_rep123.merged.sorted.bam -b2 hESC_H9_d9_BPK25_input_rep123.merged.sorted.bam -o H9_d9_BPK25_IP_TET2_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
bamCompare -b1 hESC_H9_d9_DMSO_IP_TET2_rep123.merged.sorted.bam -b2 hESC_H9_d9_DMSO_input_rep123.merged.sorted.bam -o H9_d9_DMSO_IP_TET2_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S H9_d9_DMSO_IP_TET2_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_TET2_rep123_merged_TvsC.bigwig -R d9_IP_TET2.DMSO_only.summit.bed d9_IP_TET2.DB_common.summit.bed d9_IP_TET2.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET2_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET2_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET2_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d9_TET2_DMSO_only(15243)" "d9_TET2_DB_common(493)" "d9_TET2_BPK25_only(1540)" --samplesLabel "d9_DMSO_TET2_rep123_TvsC" "d9_BPK25_TET2_rep123_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S H9_d9_DMSO_IP_TET2_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_TET2_rep123_merged_TvsC.bigwig -R d9_TET2.ROI_DMSO.summit.sorted.bed d9_TET2.RIB_DB.summit.sorted.bed d9_TET2.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET2_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET2_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET2_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d9_TET2_ROI_DMSO(519)" "d9_TET2_RIB_DB(58)" "d9_TET2_ROI_BPK25(13)" --samplesLabel "d9_DMSO_TET2_rep123_TvsC" "d9_BPK25_TET2_rep123_TvsC"

```

##  d0_TET3_ChIP_BPK25_vs_DMSO

### use individual bam of each repeat

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir plotHeatMap_use_TET3_d0_BPK25_vs_DMSO_diff_summit
cd plotHeatMap_use_TET3_d0_BPK25_vs_DMSO_diff_summit
# 创建 bam 和 bai 文件链接 
# 上传之前构造的脚本并运行
bash 6_make_H9_d0_DMSO_BPK25_TET3_ChIP_Seq_sorted_bam_file_links.sh
ls -1
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_TET3_rep123_together/d0_IP_TET3.DMSO_only.summit.bed d0_IP_TET3.DMSO_only.summit.bed
ln -s ../macs2_peak_calling_2_for_TET3_rep123_together/d0_IP_TET3.BPK25_only.summit.bed d0_IP_TET3.BPK25_only.summit.bed
ln -s ../macs2_peak_calling_2_for_TET3_rep123_together/d0_IP_TET3.DB_common.summit.bed d0_IP_TET3.DB_common.summit.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_TET3_bdgdiff/d0_TET3.ROI_DMSO.summit.sorted.bed d0_TET3.ROI_DMSO.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_TET3_bdgdiff/d0_TET3.RIB_DB.summit.sorted.bed d0_TET3.RIB_DB.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_TET3_bdgdiff/d0_TET3.ROI_BPK25.summit.sorted.bed d0_TET3.ROI_BPK25.summit.sorted.bed
#
ls -1
wc -l *.bed
# 10388 d0_IP_TET3.BPK25_only.summit.bed
#  6651 d0_IP_TET3.DB_common.summit.bed
# 34047 d0_IP_TET3.DMSO_only.summit.bed
#   633 d0_TET3.RIB_DB.summit.sorted.bed
#    51 d0_TET3.ROI_BPK25.summit.sorted.bed
#  3234 d0_TET3.ROI_DMSO.summit.sorted.bed

```

#### use bamCompare  to create log2 transformed bigwig files

```bash
# 构造 7_run_bam_comCompare_for_d0_DMSO_BPK25_TET3_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# # DMSO 
# bamCompare -b1 hESC_H9_d0_DMSO_IP_TET3_rep1.bowtie2.sorted.bam -b2 hESC_H9_d0_DMSO_input_rep1.bowtie2.sorted.bam -o hESC_H9_d0_DMSO_IP_TET3_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d0_DMSO_IP_TET3_rep2.bowtie2.sorted.bam -b2 hESC_H9_d0_DMSO_input_rep2.bowtie2.sorted.bam -o hESC_H9_d0_DMSO_IP_TET3_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d0_DMSO_IP_TET3_rep3.bowtie2.sorted.bam -b2 hESC_H9_d0_DMSO_input_rep3.bowtie2.sorted.bam -o hESC_H9_d0_DMSO_IP_TET3_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# # BPK25
# bamCompare -b1 hESC_H9_d0_BPK25_IP_TET3_rep1.bowtie2.sorted.bam -b2 hESC_H9_d0_BPK25_input_rep1.bowtie2.sorted.bam -o hESC_H9_d0_BPK25_IP_TET3_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d0_BPK25_IP_TET3_rep2.bowtie2.sorted.bam -b2 hESC_H9_d0_BPK25_input_rep2.bowtie2.sorted.bam -o hESC_H9_d0_BPK25_IP_TET3_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d0_BPK25_IP_TET3_rep3.bowtie2.sorted.bam -b2 hESC_H9_d0_BPK25_input_rep3.bowtie2.sorted.bam -o hESC_H9_d0_BPK25_IP_TET3_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
nohup bash 7_run_bam_comCompare_for_d0_DMSO_BPK25_TET3_each_repeat.sh >7_run_bam_comCompare_for_d0_DMSO_BPK25_TET3_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bam_comCompare_for_d0_DMSO_BPK25_TET3_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d0_DMSO_IP_TET3_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET3_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET3_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET3_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET3_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET3_rep3_TvsC.bigwig -R d0_IP_TET3.DMSO_only.summit.bed d0_IP_TET3.DB_common.summit.bed d0_IP_TET3.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET3_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET3_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET3_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d0_TET3_DMSO_only(34047)" "d0_TET3_DB_common(6651)" "d0_TET3_BPK25_only(10388)" --samplesLabel "d0_DMSO_TET3_rep1_TvsC" "d0_DMSO_TET3_rep2_TvsC" "d0_DMSO_TET3_rep3_TvsC" "d0_BPK25_TET3_rep1_TvsC" "d0_BPK25_TET3_rep2_TvsC" "d0_BPK25_TET3_rep3_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d0_DMSO_IP_TET3_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET3_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET3_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET3_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET3_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET3_rep3_TvsC.bigwig -R d0_TET3.ROI_DMSO.summit.sorted.bed d0_TET3.RIB_DB.summit.sorted.bed d0_TET3.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET3_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET3_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET3_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d0_TET3_ROI_DMSO(3234)" "d0_TET3_RIB_DB(633)" "d0_TET3_ROI_BPK25(51)" --samplesLabel "d0_DMSO_TET3_rep1_TvsC" "d0_DMSO_TET3_rep2_TvsC" "d0_DMSO_TET3_rep3_TvsC" "d0_BPK25_TET3_rep1_TvsC" "d0_BPK25_TET3_rep2_TvsC" "d0_BPK25_TET3_rep3_TvsC"

```

### use rep123 bam together in one

#### make rep123 merge bam file links

```bash
# 上传之前构造好的脚本并运行
bash 7_make_H9_d049_DMSO_BPK25_TET3_ChIP_Seq_rep123_merged_bam_file_links.sh
rm hESC_H9_d[49]*rep123*.bam
rm hESC_H9_d[49]*rep123*.bai

```

#### make log2 trans bigwig by bamCompare

```bash
ls -1 *_rep123.merged.sorted.bam
# hESC_H9_d0_BPK25_input_rep123.merged.sorted.bam
# hESC_H9_d0_BPK25_IP_TET3_rep123.merged.sorted.bam
# hESC_H9_d0_DMSO_input_rep123.merged.sorted.bam
# hESC_H9_d0_DMSO_IP_TET3_rep123.merged.sorted.bam
bamCompare -b1 hESC_H9_d0_BPK25_IP_TET3_rep123.merged.sorted.bam -b2 hESC_H9_d0_BPK25_input_rep123.merged.sorted.bam -o H9_d0_BPK25_IP_TET3_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
bamCompare -b1 hESC_H9_d0_DMSO_IP_TET3_rep123.merged.sorted.bam -b2 hESC_H9_d0_DMSO_input_rep123.merged.sorted.bam -o H9_d0_DMSO_IP_TET3_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S H9_d0_DMSO_IP_TET3_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_TET3_rep123_merged_TvsC.bigwig -R d0_IP_TET3.DMSO_only.summit.bed d0_IP_TET3.DB_common.summit.bed d0_IP_TET3.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET3_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET3_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET3_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d0_TET3_DMSO_only(34047)" "d0_TET3_DB_common(6651)" "d0_TET3_BPK25_only(10388)" --samplesLabel "d0_DMSO_TET3_rep123_TvsC" "d0_BPK25_TET3_rep123_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S H9_d0_DMSO_IP_TET3_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_TET3_rep123_merged_TvsC.bigwig -R d0_TET3.ROI_DMSO.summit.sorted.bed d0_TET3.RIB_DB.summit.sorted.bed d0_TET3.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET3_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET3_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d0_DMSO_BPK25_TET3_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d0_TET3_ROI_DMSO(3234)" "d0_TET3_RIB_DB(633)" "d0_TET3_ROI_BPK25(51)" --samplesLabel "d0_DMSO_TET3_rep123_TvsC" "d0_BPK25_TET3_rep123_TvsC"

```



##  d4_TET3_ChIP_BPK25_vs_DMSO

### use individual bam of each repeat

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir plotHeatMap_use_TET3_d4_BPK25_vs_DMSO_diff_summit
cd plotHeatMap_use_TET3_d4_BPK25_vs_DMSO_diff_summit
# 创建 bam 和 bai 文件链接 
# 上传之前构造的脚本并运行
bash 6_make_H9_d4_DMSO_BPK25_TET3_ChIP_Seq_sorted_bam_file_links.sh
ls -1
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_TET3_rep123_together/d4_IP_TET3.DMSO_only.summit.bed d4_IP_TET3.DMSO_only.summit.bed
ln -s ../macs2_peak_calling_2_for_TET3_rep123_together/d4_IP_TET3.BPK25_only.summit.bed d4_IP_TET3.BPK25_only.summit.bed
ln -s ../macs2_peak_calling_2_for_TET3_rep123_together/d4_IP_TET3.DB_common.summit.bed d4_IP_TET3.DB_common.summit.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_TET3_bdgdiff/d4_TET3.ROI_DMSO.summit.sorted.bed d4_TET3.ROI_DMSO.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_TET3_bdgdiff/d4_TET3.RIB_DB.summit.sorted.bed d4_TET3.RIB_DB.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_TET3_bdgdiff/d4_TET3.ROI_BPK25.summit.sorted.bed d4_TET3.ROI_BPK25.summit.sorted.bed
#
ls -1
wc -l *.bed
#  9384 d4_IP_TET3.BPK25_only.summit.bed
#  4379 d4_IP_TET3.DB_common.summit.bed
# 32042 d4_IP_TET3.DMSO_only.summit.bed
#  1772 d4_TET3.RIB_DB.summit.sorted.bed
#   136 d4_TET3.ROI_BPK25.summit.sorted.bed
#   396 d4_TET3.ROI_DMSO.summit.sorted.bed

```

#### use bamCompare  to create log2 transformed bigwig files

```bash
# 构造 7_run_bam_comCompare_for_d4_DMSO_BPK25_TET3_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# # DMSO 
# bamCompare -b1 hESC_H9_d4_DMSO_IP_TET3_rep1.bowtie2.sorted.bam -b2 hESC_H9_d4_DMSO_input_rep1.bowtie2.sorted.bam -o hESC_H9_d4_DMSO_IP_TET3_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d4_DMSO_IP_TET3_rep2.bowtie2.sorted.bam -b2 hESC_H9_d4_DMSO_input_rep2.bowtie2.sorted.bam -o hESC_H9_d4_DMSO_IP_TET3_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d4_DMSO_IP_TET3_rep3.bowtie2.sorted.bam -b2 hESC_H9_d4_DMSO_input_rep3.bowtie2.sorted.bam -o hESC_H9_d4_DMSO_IP_TET3_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# # BPK25
# bamCompare -b1 hESC_H9_d4_BPK25_IP_TET3_rep1.bowtie2.sorted.bam -b2 hESC_H9_d4_BPK25_input_rep1.bowtie2.sorted.bam -o hESC_H9_d4_BPK25_IP_TET3_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d4_BPK25_IP_TET3_rep2.bowtie2.sorted.bam -b2 hESC_H9_d4_BPK25_input_rep2.bowtie2.sorted.bam -o hESC_H9_d4_BPK25_IP_TET3_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d4_BPK25_IP_TET3_rep3.bowtie2.sorted.bam -b2 hESC_H9_d4_BPK25_input_rep3.bowtie2.sorted.bam -o hESC_H9_d4_BPK25_IP_TET3_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
nohup bash 7_run_bam_comCompare_for_d4_DMSO_BPK25_TET3_each_repeat.sh >7_run_bam_comCompare_for_d4_DMSO_BPK25_TET3_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bam_comCompare_for_d4_DMSO_BPK25_TET3_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash

```

### use rep123 bam together in one

#### make rep123 merge bam file links

```bash
# 上传之前构造好的脚本并运行
bash 7_make_H9_d049_DMSO_BPK25_TET3_ChIP_Seq_rep123_merged_bam_file_links.sh
rm hESC_H9_d[09]*rep123*.bam
rm hESC_H9_d[09]*rep123*.bai

```

#### make log2 trans bigwig by bamCompare

```bash
ls -1 *_rep123.merged.sorted.bam
# hESC_H9_d4_BPK25_input_rep123.merged.sorted.bam
# hESC_H9_d4_BPK25_IP_TET3_rep123.merged.sorted.bam
# hESC_H9_d4_DMSO_input_rep123.merged.sorted.bam
# hESC_H9_d4_DMSO_IP_TET3_rep123.merged.sorted.bam
bamCompare -b1 hESC_H9_d4_BPK25_IP_TET3_rep123.merged.sorted.bam -b2 hESC_H9_d4_BPK25_input_rep123.merged.sorted.bam -o H9_d4_BPK25_IP_TET3_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
bamCompare -b1 hESC_H9_d4_DMSO_IP_TET3_rep123.merged.sorted.bam -b2 hESC_H9_d4_DMSO_input_rep123.merged.sorted.bam -o H9_d4_DMSO_IP_TET3_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S H9_d4_DMSO_IP_TET3_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_TET3_rep123_merged_TvsC.bigwig -R d4_IP_TET3.DMSO_only.summit.bed d4_IP_TET3.DB_common.summit.bed d4_IP_TET3.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET3_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET3_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET3_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d4_TET3_DMSO_only(32042)" "d4_TET3_DB_common(4379)" "d4_TET3_BPK25_only(9384)" --samplesLabel "d4_DMSO_TET3_rep123_TvsC" "d4_BPK25_TET3_rep123_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S H9_d4_DMSO_IP_TET3_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_TET3_rep123_merged_TvsC.bigwig -R d4_TET3.ROI_DMSO.summit.sorted.bed d4_TET3.RIB_DB.summit.sorted.bed d4_TET3.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET3_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET3_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d4_DMSO_BPK25_TET3_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d4_TET3_ROI_DMSO(396)" "d4_TET3_RIB_DB(1772)" "d4_TET3_ROI_BPK25(136)" --samplesLabel "d4_DMSO_TET3_rep123_TvsC" "d4_BPK25_TET3_rep123_TvsC"

```



##  d9_TET3_ChIP_BPK25_vs_DMSO

### use individual bam of each repeat

#### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir plotHeatMap_use_TET3_d9_BPK25_vs_DMSO_diff_summit
cd plotHeatMap_use_TET3_d9_BPK25_vs_DMSO_diff_summit
# 创建 bam 和 bai 文件链接 
# 上传之前构造的脚本并运行
bash 6_make_H9_d9_DMSO_BPK25_TET3_ChIP_Seq_sorted_bam_file_links.sh
ls -1
# 创建 intersect 三类 summit 文件连接
ln -s ../macs2_peak_calling_2_for_TET3_rep123_together/d9_IP_TET3.DMSO_only.summit.bed d9_IP_TET3.DMSO_only.summit.bed
ln -s ../macs2_peak_calling_2_for_TET3_rep123_together/d9_IP_TET3.BPK25_only.summit.bed d9_IP_TET3.BPK25_only.summit.bed
ln -s ../macs2_peak_calling_2_for_TET3_rep123_together/d9_IP_TET3.DB_common.summit.bed d9_IP_TET3.DB_common.summit.bed
# 创建 bdgdiff 三类 summit 文件连接
ln -s ../macs2_peak_calling_3_for_TET3_bdgdiff/d9_TET3.ROI_DMSO.summit.sorted.bed d9_TET3.ROI_DMSO.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_TET3_bdgdiff/d9_TET3.RIB_DB.summit.sorted.bed d9_TET3.RIB_DB.summit.sorted.bed
ln -s ../macs2_peak_calling_3_for_TET3_bdgdiff/d9_TET3.ROI_BPK25.summit.sorted.bed d9_TET3.ROI_BPK25.summit.sorted.bed
#
ls -1
wc -l *.bed
#   197 d9_IP_TET3.BPK25_only.summit.bed
#   159 d9_IP_TET3.DB_common.summit.bed
# 28334 d9_IP_TET3.DMSO_only.summit.bed
#   146 d9_TET3.RIB_DB.summit.sorted.bed
#   118 d9_TET3.ROI_BPK25.summit.sorted.bed
#   248 d9_TET3.ROI_DMSO.summit.sorted.bed

```

#### use bamCompare  to create log2 transformed bigwig files

```bash
# 构造 7_run_bam_comCompare_for_d9_DMSO_BPK25_TET3_each_repeat.sh 内容如下
# # Compare for bam of each repeat
# # DMSO 
# bamCompare -b1 hESC_H9_d9_DMSO_IP_TET3_rep1.bowtie2.sorted.bam -b2 hESC_H9_d9_DMSO_input_rep1.bowtie2.sorted.bam -o hESC_H9_d9_DMSO_IP_TET3_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d9_DMSO_IP_TET3_rep2.bowtie2.sorted.bam -b2 hESC_H9_d9_DMSO_input_rep2.bowtie2.sorted.bam -o hESC_H9_d9_DMSO_IP_TET3_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d9_DMSO_IP_TET3_rep3.bowtie2.sorted.bam -b2 hESC_H9_d9_DMSO_input_rep3.bowtie2.sorted.bam -o hESC_H9_d9_DMSO_IP_TET3_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# # BPK25
# bamCompare -b1 hESC_H9_d9_BPK25_IP_TET3_rep1.bowtie2.sorted.bam -b2 hESC_H9_d9_BPK25_input_rep1.bowtie2.sorted.bam -o hESC_H9_d9_BPK25_IP_TET3_rep1_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d9_BPK25_IP_TET3_rep2.bowtie2.sorted.bam -b2 hESC_H9_d9_BPK25_input_rep2.bowtie2.sorted.bam -o hESC_H9_d9_BPK25_IP_TET3_rep2_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
# bamCompare -b1 hESC_H9_d9_BPK25_IP_TET3_rep3.bowtie2.sorted.bam -b2 hESC_H9_d9_BPK25_input_rep3.bowtie2.sorted.bam -o hESC_H9_d9_BPK25_IP_TET3_rep3_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
nohup bash 7_run_bam_comCompare_for_d9_DMSO_BPK25_TET3_each_repeat.sh >7_run_bam_comCompare_for_d9_DMSO_BPK25_TET3_each_repeat.sh.runlog 2>&1 &
tail -f 7_run_bam_comCompare_for_d9_DMSO_BPK25_TET3_each_repeat.sh.runlog

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d9_DMSO_IP_TET3_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET3_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET3_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET3_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET3_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET3_rep3_TvsC.bigwig -R d9_IP_TET3.DMSO_only.summit.bed d9_IP_TET3.DB_common.summit.bed d9_IP_TET3.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET3_samples_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET3_samples_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET3_samples_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d9_TET3_DMSO_only(28334)" "d9_TET3_DB_common(159)" "d9_TET3_BPK25_only(197)" --samplesLabel "d9_DMSO_TET3_rep1_TvsC" "d9_DMSO_TET3_rep2_TvsC" "d9_DMSO_TET3_rep3_TvsC" "d9_BPK25_TET3_rep1_TvsC" "d9_BPK25_TET3_rep2_TvsC" "d9_BPK25_TET3_rep3_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S hESC_H9_d9_DMSO_IP_TET3_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET3_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET3_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET3_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET3_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET3_rep3_TvsC.bigwig -R d9_TET3.ROI_DMSO.summit.sorted.bed d9_TET3.RIB_DB.summit.sorted.bed d9_TET3.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET3_samples_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET3_samples_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET3_samples_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d9_TET3_ROI_DMSO(248)" "d9_TET3_RIB_DB(146)" "d9_TET3_ROI_BPK25(118)" --samplesLabel "d9_DMSO_TET3_rep1_TvsC" "d9_DMSO_TET3_rep2_TvsC" "d9_DMSO_TET3_rep3_TvsC" "d9_BPK25_TET3_rep1_TvsC" "d9_BPK25_TET3_rep2_TvsC" "d9_BPK25_TET3_rep3_TvsC"

```

### use rep123 bam together in one

#### make rep123 merge bam file links

```bash
# 上传之前构造好的脚本并运行
bash 7_make_H9_d049_DMSO_BPK25_TET3_ChIP_Seq_rep123_merged_bam_file_links.sh
rm hESC_H9_d[04]*rep123*.bam
rm hESC_H9_d[04]*rep123*.bai

```

#### make log2 trans bigwig by bamCompare

```bash
ls -1 *_rep123.merged.sorted.bam
# hESC_H9_d9_BPK25_input_rep123.merged.sorted.bam
# hESC_H9_d9_BPK25_IP_TET3_rep123.merged.sorted.bam
# hESC_H9_d9_DMSO_input_rep123.merged.sorted.bam
# hESC_H9_d9_DMSO_IP_TET3_rep123.merged.sorted.bam
bamCompare -b1 hESC_H9_d9_BPK25_IP_TET3_rep123.merged.sorted.bam -b2 hESC_H9_d9_BPK25_input_rep123.merged.sorted.bam -o H9_d9_BPK25_IP_TET3_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads
bamCompare -b1 hESC_H9_d9_DMSO_IP_TET3_rep123.merged.sorted.bam -b2 hESC_H9_d9_DMSO_input_rep123.merged.sorted.bam -o H9_d9_DMSO_IP_TET3_rep123_merged_TvsC.bigwig -of bigwig --skipZeroOverZero -p "max" --exactScaling --skipNAs --extendReads --ignoreDuplicates --centerReads

```

#### computeMatrix and plotHeatmap

```bash
# use summits based on intersect
computeMatrix reference-point -p 32 --referencePoint center -S H9_d9_DMSO_IP_TET3_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_TET3_rep123_merged_TvsC.bigwig -R d9_IP_TET3.DMSO_only.summit.bed d9_IP_TET3.DB_common.summit.bed d9_IP_TET3.BPK25_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET3_rep123_merged_for_plotHeatMap.intersect.gz
plotHeatmap -m bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET3_rep123_merged_for_plotHeatMap.intersect.gz -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET3_rep123_merged_for_plotHeatMap.intersect.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d9_TET3_DMSO_only(28334)" "d9_TET3_DB_common(159)" "d9_TET3_BPK25_only(197)" --samplesLabel "d9_DMSO_TET3_rep123_TvsC" "d9_BPK25_TET3_rep123_TvsC"

# use summits generated by bdgdiff
computeMatrix reference-point -p 32 --referencePoint center -S H9_d9_DMSO_IP_TET3_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_TET3_rep123_merged_TvsC.bigwig -R d9_TET3.ROI_DMSO.summit.sorted.bed d9_TET3.RIB_DB.summit.sorted.bed d9_TET3.ROI_BPK25.summit.sorted.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET3_rep123_merged_for_plotHeatMap.bdgdiff.gz
plotHeatmap -m bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET3_rep123_merged_for_plotHeatMap.bdgdiff.gz -o bamCompare_matrix_of_H9_d9_DMSO_BPK25_TET3_rep123_merged_for_plotHeatMap.bdgdiff.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "d9_TET3_ROI_DMSO(248)" "d9_TET3_RIB_DB(146)" "d9_TET3_ROI_BPK25(118)" --samplesLabel "d9_DMSO_TET3_rep123_TvsC" "d9_BPK25_TET3_rep123_TvsC"

```





# heatmap_of_reads_enrichment_around_H9_WT_d49_CTCF_diff_summits_by_intersect

## d0_DMSO_BPK25_CTCF

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d0_DB_CTCF_use_WT_d49_diff_CTCF_summits
cd reads_heatmap_of_d0_DB_CTCF_use_WT_d49_diff_CTCF_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d4_only.summit.bed d49_IP_CTCF.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d49_common.summit.bed d49_IP_CTCF.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d9_only.summit.bed d49_IP_CTCF.d9_only.summit.bed
#
ls -1
wc -l *.bed
# 19095 d49_IP_CTCF.d4_only.summit.bed
#  1547 d49_IP_CTCF.d49_common.summit.bed
#  7511 d49_IP_CTCF.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d0_DB_CTCF_use_WT_d49_diff_CTCF_summits
# each repeat
ln -s ../plotHeatMap_use_CTCF_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d0_DMSO_IP_CTCF_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_CTCF_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_CTCF_d0_BPK25_vs_DMSO_diff_summit/H9_d0_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig H9_d0_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d0_BPK25_vs_DMSO_diff_summit/H9_d0_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 CTCF
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d0_DMSO_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_CTCF_rep3_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d0_DB_CTCF_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d0_DB_CTCF_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d0_DB_CTCF_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d0_DMSO_CTCF_rep1" "d0_DMSO_CTCF_rep2" "d0_DMSO_CTCF_rep3" "d0_BPK25_CTCF_rep1" "d0_BPK25_CTCF_rep2" "d0_BPK25_CTCF_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d0_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d0_DB_CTCF_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d0_DB_CTCF_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d0_DB_CTCF_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d0_DMSO_CTCF_rep123" "d0_BPK25_CTCF_rep123"

```

## d4_DMSO_BPK25_CTCF

### data preperation
```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d4_DB_CTCF_use_WT_d49_diff_CTCF_summits
cd reads_heatmap_of_d4_DB_CTCF_use_WT_d49_diff_CTCF_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d4_only.summit.bed d49_IP_CTCF.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d49_common.summit.bed d49_IP_CTCF.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d9_only.summit.bed d49_IP_CTCF.d9_only.summit.bed
#
ls -1
wc -l *.bed
# 19095 d49_IP_CTCF.d4_only.summit.bed
#  1547 d49_IP_CTCF.d49_common.summit.bed
#  7511 d49_IP_CTCF.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d4_DB_CTCF_use_WT_d49_diff_CTCF_summits
# each repeat
ln -s ../plotHeatMap_use_CTCF_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d4_DMSO_IP_CTCF_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_CTCF_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_CTCF_d4_BPK25_vs_DMSO_diff_summit/H9_d4_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig H9_d4_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d4_BPK25_vs_DMSO_diff_summit/H9_d4_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 CTCF
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d4_DMSO_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_CTCF_rep3_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d4_DB_CTCF_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d4_DB_CTCF_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d4_DB_CTCF_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d4_DMSO_CTCF_rep1" "d4_DMSO_CTCF_rep2" "d4_DMSO_CTCF_rep3" "d4_BPK25_CTCF_rep1" "d4_BPK25_CTCF_rep2" "d4_BPK25_CTCF_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d4_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d4_DB_CTCF_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d4_DB_CTCF_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d4_DB_CTCF_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d4_DMSO_CTCF_rep123" "d4_BPK25_CTCF_rep123"

```

## d9_DMSO_BPK25_CTCF

### data preperation
```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d9_DB_CTCF_use_WT_d49_diff_CTCF_summits
cd reads_heatmap_of_d9_DB_CTCF_use_WT_d49_diff_CTCF_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d4_only.summit.bed d49_IP_CTCF.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d49_common.summit.bed d49_IP_CTCF.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d9_only.summit.bed d49_IP_CTCF.d9_only.summit.bed
#
ls -1
wc -l *.bed
# 19095 d49_IP_CTCF.d4_only.summit.bed
#  1547 d49_IP_CTCF.d49_common.summit.bed
#  7511 d49_IP_CTCF.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d9_DB_CTCF_use_WT_d49_diff_CTCF_summits
# each repeat
ln -s ../plotHeatMap_use_CTCF_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d9_DMSO_IP_CTCF_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_CTCF_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_CTCF_d9_BPK25_vs_DMSO_diff_summit/H9_d9_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig H9_d9_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d9_BPK25_vs_DMSO_diff_summit/H9_d9_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 CTCF
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d9_DMSO_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_CTCF_rep3_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d9_DB_CTCF_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d9_DB_CTCF_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d9_DB_CTCF_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d9_DMSO_CTCF_rep1" "d9_DMSO_CTCF_rep2" "d9_DMSO_CTCF_rep3" "d9_BPK25_CTCF_rep1" "d9_BPK25_CTCF_rep2" "d9_BPK25_CTCF_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d9_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d9_DB_CTCF_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d9_DB_CTCF_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d9_DB_CTCF_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d9_DMSO_CTCF_rep123" "d9_BPK25_CTCF_rep123"

```

## d0_DMSO_BPK25_MBD3

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d0_DB_MBD3_use_WT_d49_diff_CTCF_summits
cd reads_heatmap_of_d0_DB_MBD3_use_WT_d49_diff_CTCF_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d4_only.summit.bed d49_IP_CTCF.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d49_common.summit.bed d49_IP_CTCF.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d9_only.summit.bed d49_IP_CTCF.d9_only.summit.bed
#
ls -1
wc -l *.bed
# 19095 d49_IP_CTCF.d4_only.summit.bed
#  1547 d49_IP_CTCF.d49_common.summit.bed
#  7511 d49_IP_CTCF.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d0_DB_MBD3_use_WT_d49_diff_CTCF_summits
# each repeat
ln -s ../plotHeatMap_use_MBD3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_MBD3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_MBD3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d0_DMSO_IP_MBD3_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_MBD3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_MBD3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_MBD3_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_MBD3_d0_BPK25_vs_DMSO_diff_summit/H9_d0_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig H9_d0_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d0_BPK25_vs_DMSO_diff_summit/H9_d0_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 CTCF
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d0_DMSO_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_MBD3_rep3_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d0_DB_MBD3_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d0_DB_MBD3_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d0_DB_MBD3_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d0_DMSO_MBD3_rep1" "d0_DMSO_MBD3_rep2" "d0_DMSO_MBD3_rep3" "d0_BPK25_MBD3_rep1" "d0_BPK25_MBD3_rep2" "d0_BPK25_MBD3_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d0_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d0_DB_MBD3_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d0_DB_MBD3_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d0_DB_MBD3_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d0_DMSO_MBD3_rep123" "d0_BPK25_MBD3_rep123"

```

## d4_DMSO_BPK25_MBD3

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d4_DB_MBD3_use_WT_d49_diff_CTCF_summits
cd reads_heatmap_of_d4_DB_MBD3_use_WT_d49_diff_CTCF_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d4_only.summit.bed d49_IP_CTCF.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d49_common.summit.bed d49_IP_CTCF.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d9_only.summit.bed d49_IP_CTCF.d9_only.summit.bed
#
ls -1
wc -l *.bed
# 19095 d49_IP_CTCF.d4_only.summit.bed
#  1547 d49_IP_CTCF.d49_common.summit.bed
#  7511 d49_IP_CTCF.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d4_DB_MBD3_use_WT_d49_diff_CTCF_summits
# each repeat
ln -s ../plotHeatMap_use_MBD3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_MBD3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_MBD3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d4_DMSO_IP_MBD3_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_MBD3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_MBD3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_MBD3_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_MBD3_d4_BPK25_vs_DMSO_diff_summit/H9_d4_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig H9_d4_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d4_BPK25_vs_DMSO_diff_summit/H9_d4_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 CTCF
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d4_DMSO_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_MBD3_rep3_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d4_DB_MBD3_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d4_DB_MBD3_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d4_DB_MBD3_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d4_DMSO_MBD3_rep1" "d4_DMSO_MBD3_rep2" "d4_DMSO_MBD3_rep3" "d4_BPK25_MBD3_rep1" "d4_BPK25_MBD3_rep2" "d4_BPK25_MBD3_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d4_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d4_DB_MBD3_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d4_DB_MBD3_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d4_DB_MBD3_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d4_DMSO_MBD3_rep123" "d4_BPK25_MBD3_rep123"

```

## d9_DMSO_BPK25_MBD3

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d9_DB_MBD3_use_WT_d49_diff_CTCF_summits
cd reads_heatmap_of_d9_DB_MBD3_use_WT_d49_diff_CTCF_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d4_only.summit.bed d49_IP_CTCF.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d49_common.summit.bed d49_IP_CTCF.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d9_only.summit.bed d49_IP_CTCF.d9_only.summit.bed
#
ls -1
wc -l *.bed
# 19095 d49_IP_CTCF.d4_only.summit.bed
#  1547 d49_IP_CTCF.d49_common.summit.bed
#  7511 d49_IP_CTCF.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d9_DB_MBD3_use_WT_d49_diff_CTCF_summits
# each repeat
ln -s ../plotHeatMap_use_MBD3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_MBD3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_MBD3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d9_DMSO_IP_MBD3_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_MBD3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_MBD3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_MBD3_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_MBD3_d9_BPK25_vs_DMSO_diff_summit/H9_d9_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig H9_d9_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d9_BPK25_vs_DMSO_diff_summit/H9_d9_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 CTCF
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d9_DMSO_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_MBD3_rep3_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d9_DB_MBD3_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d9_DB_MBD3_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d9_DB_MBD3_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d9_DMSO_MBD3_rep1" "d9_DMSO_MBD3_rep2" "d9_DMSO_MBD3_rep3" "d9_BPK25_MBD3_rep1" "d9_BPK25_MBD3_rep2" "d9_BPK25_MBD3_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d9_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d9_DB_MBD3_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d9_DB_MBD3_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d9_DB_MBD3_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d9_DMSO_MBD3_rep123" "d9_BPK25_MBD3_rep123"

```

## d0_DMSO_BPK25_TET1

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d0_DB_TET1_use_WT_d49_diff_CTCF_summits
cd reads_heatmap_of_d0_DB_TET1_use_WT_d49_diff_CTCF_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d4_only.summit.bed d49_IP_CTCF.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d49_common.summit.bed d49_IP_CTCF.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d9_only.summit.bed d49_IP_CTCF.d9_only.summit.bed
#
ls -1
wc -l *.bed
# 19095 d49_IP_CTCF.d4_only.summit.bed
#  1547 d49_IP_CTCF.d49_common.summit.bed
#  7511 d49_IP_CTCF.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d0_DB_TET1_use_WT_d49_diff_CTCF_summits
# each repeat
ln -s ../plotHeatMap_use_TET1_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_TET1_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET1_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_TET1_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET1_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_TET1_rep3_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET1_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_TET1_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET1_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_TET1_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET1_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_TET1_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET1_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_TET1_d0_BPK25_vs_DMSO_diff_summit/H9_d0_DMSO_IP_TET1_rep123_merged_TvsC.bigwig H9_d0_DMSO_IP_TET1_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d0_BPK25_vs_DMSO_diff_summit/H9_d0_BPK25_IP_TET1_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_TET1_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 CTCF
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d0_DMSO_IP_TET1_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET1_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET1_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET1_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET1_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET1_rep3_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d0_DB_TET1_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d0_DB_TET1_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d0_DB_TET1_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d0_DMSO_TET1_rep1" "d0_DMSO_TET1_rep2" "d0_DMSO_TET1_rep3" "d0_BPK25_TET1_rep1" "d0_BPK25_TET1_rep2" "d0_BPK25_TET1_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d0_DMSO_IP_TET1_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_TET1_rep123_merged_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d0_DB_TET1_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d0_DB_TET1_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d0_DB_TET1_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d0_DMSO_TET1_rep123" "d0_BPK25_TET1_rep123"

```

## d4_DMSO_BPK25_TET1

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d4_DB_TET1_use_WT_d49_diff_CTCF_summits
cd reads_heatmap_of_d4_DB_TET1_use_WT_d49_diff_CTCF_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d4_only.summit.bed d49_IP_CTCF.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d49_common.summit.bed d49_IP_CTCF.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d9_only.summit.bed d49_IP_CTCF.d9_only.summit.bed
#
ls -1
wc -l *.bed
# 19095 d49_IP_CTCF.d4_only.summit.bed
#  1547 d49_IP_CTCF.d49_common.summit.bed
#  7511 d49_IP_CTCF.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d4_DB_TET1_use_WT_d49_diff_CTCF_summits
# each repeat
ln -s ../plotHeatMap_use_TET1_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_TET1_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET1_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_TET1_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET1_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_TET1_rep3_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET1_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_TET1_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET1_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_TET1_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET1_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_TET1_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET1_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_TET1_d4_BPK25_vs_DMSO_diff_summit/H9_d4_DMSO_IP_TET1_rep123_merged_TvsC.bigwig H9_d4_DMSO_IP_TET1_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d4_BPK25_vs_DMSO_diff_summit/H9_d4_BPK25_IP_TET1_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_TET1_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 CTCF
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d4_DMSO_IP_TET1_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET1_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET1_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET1_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET1_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET1_rep3_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d4_DB_TET1_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d4_DB_TET1_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d4_DB_TET1_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d4_DMSO_TET1_rep1" "d4_DMSO_TET1_rep2" "d4_DMSO_TET1_rep3" "d4_BPK25_TET1_rep1" "d4_BPK25_TET1_rep2" "d4_BPK25_TET1_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d4_DMSO_IP_TET1_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_TET1_rep123_merged_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d4_DB_TET1_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d4_DB_TET1_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d4_DB_TET1_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d4_DMSO_TET1_rep123" "d4_BPK25_TET1_rep123"

```

## d9_DMSO_BPK25_TET1

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d9_DB_TET1_use_WT_d49_diff_CTCF_summits
cd reads_heatmap_of_d9_DB_TET1_use_WT_d49_diff_CTCF_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d4_only.summit.bed d49_IP_CTCF.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d49_common.summit.bed d49_IP_CTCF.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d9_only.summit.bed d49_IP_CTCF.d9_only.summit.bed
#
ls -1
wc -l *.bed
# 19095 d49_IP_CTCF.d4_only.summit.bed
#  1547 d49_IP_CTCF.d49_common.summit.bed
#  7511 d49_IP_CTCF.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d9_DB_TET1_use_WT_d49_diff_CTCF_summits
# each repeat
ln -s ../plotHeatMap_use_TET1_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_TET1_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET1_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_TET1_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET1_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_TET1_rep3_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET1_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_TET1_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET1_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_TET1_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET1_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_TET1_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET1_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_TET1_d9_BPK25_vs_DMSO_diff_summit/H9_d9_DMSO_IP_TET1_rep123_merged_TvsC.bigwig H9_d9_DMSO_IP_TET1_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d9_BPK25_vs_DMSO_diff_summit/H9_d9_BPK25_IP_TET1_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_TET1_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 CTCF
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d9_DMSO_IP_TET1_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET1_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET1_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET1_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET1_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET1_rep3_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d9_DB_TET1_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d9_DB_TET1_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d9_DB_TET1_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d9_DMSO_TET1_rep1" "d9_DMSO_TET1_rep2" "d9_DMSO_TET1_rep3" "d9_BPK25_TET1_rep1" "d9_BPK25_TET1_rep2" "d9_BPK25_TET1_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d9_DMSO_IP_TET1_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_TET1_rep123_merged_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d9_DB_TET1_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d9_DB_TET1_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d9_DB_TET1_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d9_DMSO_TET1_rep123" "d9_BPK25_TET1_rep123"

```

## d0_DMSO_BPK25_TET2

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d0_DB_TET2_use_WT_d49_diff_CTCF_summits
cd reads_heatmap_of_d0_DB_TET2_use_WT_d49_diff_CTCF_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d4_only.summit.bed d49_IP_CTCF.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d49_common.summit.bed d49_IP_CTCF.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d9_only.summit.bed d49_IP_CTCF.d9_only.summit.bed
#
ls -1
wc -l *.bed
# 19095 d49_IP_CTCF.d4_only.summit.bed
#  1547 d49_IP_CTCF.d49_common.summit.bed
#  7511 d49_IP_CTCF.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d0_DB_TET2_use_WT_d49_diff_CTCF_summits
# each repeat
ln -s ../plotHeatMap_use_TET2_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_TET2_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET2_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_TET2_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET2_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_TET2_rep3_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET2_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_TET2_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET2_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_TET2_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET2_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_TET2_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET2_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_TET2_d0_BPK25_vs_DMSO_diff_summit/H9_d0_DMSO_IP_TET2_rep123_merged_TvsC.bigwig H9_d0_DMSO_IP_TET2_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d0_BPK25_vs_DMSO_diff_summit/H9_d0_BPK25_IP_TET2_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_TET2_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 CTCF
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d0_DMSO_IP_TET2_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET2_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET2_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET2_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET2_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET2_rep3_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d0_DB_TET2_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d0_DB_TET2_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d0_DB_TET2_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d0_DMSO_TET2_rep1" "d0_DMSO_TET2_rep2" "d0_DMSO_TET2_rep3" "d0_BPK25_TET2_rep1" "d0_BPK25_TET2_rep2" "d0_BPK25_TET2_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d0_DMSO_IP_TET2_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_TET2_rep123_merged_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d0_DB_TET2_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d0_DB_TET2_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d0_DB_TET2_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d0_DMSO_TET2_rep123" "d0_BPK25_TET2_rep123"

```

## d4_DMSO_BPK25_TET2

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d4_DB_TET2_use_WT_d49_diff_CTCF_summits
cd reads_heatmap_of_d4_DB_TET2_use_WT_d49_diff_CTCF_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d4_only.summit.bed d49_IP_CTCF.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d49_common.summit.bed d49_IP_CTCF.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d9_only.summit.bed d49_IP_CTCF.d9_only.summit.bed
#
ls -1
wc -l *.bed
# 19095 d49_IP_CTCF.d4_only.summit.bed
#  1547 d49_IP_CTCF.d49_common.summit.bed
#  7511 d49_IP_CTCF.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d4_DB_TET2_use_WT_d49_diff_CTCF_summits
# each repeat
ln -s ../plotHeatMap_use_TET2_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_TET2_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET2_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_TET2_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET2_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_TET2_rep3_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET2_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_TET2_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET2_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_TET2_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET2_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_TET2_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET2_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_TET2_d4_BPK25_vs_DMSO_diff_summit/H9_d4_DMSO_IP_TET2_rep123_merged_TvsC.bigwig H9_d4_DMSO_IP_TET2_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d4_BPK25_vs_DMSO_diff_summit/H9_d4_BPK25_IP_TET2_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_TET2_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 CTCF
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d4_DMSO_IP_TET2_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET2_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET2_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET2_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET2_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET2_rep3_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d4_DB_TET2_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d4_DB_TET2_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d4_DB_TET2_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d4_DMSO_TET2_rep1" "d4_DMSO_TET2_rep2" "d4_DMSO_TET2_rep3" "d4_BPK25_TET2_rep1" "d4_BPK25_TET2_rep2" "d4_BPK25_TET2_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d4_DMSO_IP_TET2_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_TET2_rep123_merged_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d4_DB_TET2_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d4_DB_TET2_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d4_DB_TET2_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d4_DMSO_TET2_rep123" "d4_BPK25_TET2_rep123"

```

## d9_DMSO_BPK25_TET2

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d9_DB_TET2_use_WT_d49_diff_CTCF_summits
cd reads_heatmap_of_d9_DB_TET2_use_WT_d49_diff_CTCF_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d4_only.summit.bed d49_IP_CTCF.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d49_common.summit.bed d49_IP_CTCF.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d9_only.summit.bed d49_IP_CTCF.d9_only.summit.bed
#
ls -1
wc -l *.bed
# 19095 d49_IP_CTCF.d4_only.summit.bed
#  1547 d49_IP_CTCF.d49_common.summit.bed
#  7511 d49_IP_CTCF.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d9_DB_TET2_use_WT_d49_diff_CTCF_summits
# each repeat
ln -s ../plotHeatMap_use_TET2_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_TET2_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET2_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_TET2_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET2_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_TET2_rep3_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET2_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_TET2_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET2_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_TET2_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET2_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_TET2_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET2_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_TET2_d9_BPK25_vs_DMSO_diff_summit/H9_d9_DMSO_IP_TET2_rep123_merged_TvsC.bigwig H9_d9_DMSO_IP_TET2_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d9_BPK25_vs_DMSO_diff_summit/H9_d9_BPK25_IP_TET2_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_TET2_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 CTCF
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d9_DMSO_IP_TET2_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET2_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET2_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET2_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET2_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET2_rep3_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d9_DB_TET2_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d9_DB_TET2_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d9_DB_TET2_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d9_DMSO_TET2_rep1" "d9_DMSO_TET2_rep2" "d9_DMSO_TET2_rep3" "d9_BPK25_TET2_rep1" "d9_BPK25_TET2_rep2" "d9_BPK25_TET2_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d9_DMSO_IP_TET2_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_TET2_rep123_merged_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d9_DB_TET2_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d9_DB_TET2_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d9_DB_TET2_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d9_DMSO_TET2_rep123" "d9_BPK25_TET2_rep123"

```

## d0_DMSO_BPK25_TET3

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d0_DB_TET3_use_WT_d49_diff_CTCF_summits
cd reads_heatmap_of_d0_DB_TET3_use_WT_d49_diff_CTCF_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d4_only.summit.bed d49_IP_CTCF.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d49_common.summit.bed d49_IP_CTCF.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d9_only.summit.bed d49_IP_CTCF.d9_only.summit.bed
#
ls -1
wc -l *.bed
# 19095 d49_IP_CTCF.d4_only.summit.bed
#  1547 d49_IP_CTCF.d49_common.summit.bed
#  7511 d49_IP_CTCF.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d0_DB_TET3_use_WT_d49_diff_CTCF_summits
# each repeat
ln -s ../plotHeatMap_use_TET3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_TET3_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_TET3_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_TET3_rep3_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET3_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_TET3_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_TET3_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_TET3_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET3_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_TET3_d0_BPK25_vs_DMSO_diff_summit/H9_d0_DMSO_IP_TET3_rep123_merged_TvsC.bigwig H9_d0_DMSO_IP_TET3_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d0_BPK25_vs_DMSO_diff_summit/H9_d0_BPK25_IP_TET3_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_TET3_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 CTCF
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d0_DMSO_IP_TET3_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET3_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET3_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET3_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET3_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET3_rep3_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d0_DB_TET3_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d0_DB_TET3_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d0_DB_TET3_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d0_DMSO_TET3_rep1" "d0_DMSO_TET3_rep2" "d0_DMSO_TET3_rep3" "d0_BPK25_TET3_rep1" "d0_BPK25_TET3_rep2" "d0_BPK25_TET3_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d0_DMSO_IP_TET3_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_TET3_rep123_merged_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d0_DB_TET3_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d0_DB_TET3_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d0_DB_TET3_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d0_DMSO_TET3_rep123" "d0_BPK25_TET3_rep123"

```

## d4_DMSO_BPK25_TET3

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d4_DB_TET3_use_WT_d49_diff_CTCF_summits
cd reads_heatmap_of_d4_DB_TET3_use_WT_d49_diff_CTCF_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d4_only.summit.bed d49_IP_CTCF.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d49_common.summit.bed d49_IP_CTCF.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d9_only.summit.bed d49_IP_CTCF.d9_only.summit.bed
#
ls -1
wc -l *.bed
# 19095 d49_IP_CTCF.d4_only.summit.bed
#  1547 d49_IP_CTCF.d49_common.summit.bed
#  7511 d49_IP_CTCF.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d4_DB_TET3_use_WT_d49_diff_CTCF_summits
# each repeat
ln -s ../plotHeatMap_use_TET3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_TET3_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_TET3_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_TET3_rep3_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET3_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_TET3_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_TET3_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_TET3_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET3_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_TET3_d4_BPK25_vs_DMSO_diff_summit/H9_d4_DMSO_IP_TET3_rep123_merged_TvsC.bigwig H9_d4_DMSO_IP_TET3_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d4_BPK25_vs_DMSO_diff_summit/H9_d4_BPK25_IP_TET3_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_TET3_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 CTCF
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d4_DMSO_IP_TET3_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET3_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET3_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET3_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET3_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET3_rep3_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d4_DB_TET3_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d4_DB_TET3_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d4_DB_TET3_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d4_DMSO_TET3_rep1" "d4_DMSO_TET3_rep2" "d4_DMSO_TET3_rep3" "d4_BPK25_TET3_rep1" "d4_BPK25_TET3_rep2" "d4_BPK25_TET3_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d4_DMSO_IP_TET3_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_TET3_rep123_merged_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d4_DB_TET3_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d4_DB_TET3_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d4_DB_TET3_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d4_DMSO_TET3_rep123" "d4_BPK25_TET3_rep123"

```

## d9_DMSO_BPK25_TET3

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d9_DB_TET3_use_WT_d49_diff_CTCF_summits
cd reads_heatmap_of_d9_DB_TET3_use_WT_d49_diff_CTCF_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d4_only.summit.bed d49_IP_CTCF.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d49_common.summit.bed d49_IP_CTCF.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_CTCF_rep123_together/d49_IP_CTCF_rep123_TvsC.d9_only.summit.bed d49_IP_CTCF.d9_only.summit.bed
#
ls -1
wc -l *.bed
# 19095 d49_IP_CTCF.d4_only.summit.bed
#  1547 d49_IP_CTCF.d49_common.summit.bed
#  7511 d49_IP_CTCF.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d9_DB_TET3_use_WT_d49_diff_CTCF_summits
# each repeat
ln -s ../plotHeatMap_use_TET3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_TET3_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_TET3_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_TET3_rep3_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET3_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_TET3_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_TET3_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_TET3_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET3_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_TET3_d9_BPK25_vs_DMSO_diff_summit/H9_d9_DMSO_IP_TET3_rep123_merged_TvsC.bigwig H9_d9_DMSO_IP_TET3_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d9_BPK25_vs_DMSO_diff_summit/H9_d9_BPK25_IP_TET3_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_TET3_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 CTCF
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d9_DMSO_IP_TET3_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET3_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET3_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET3_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET3_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET3_rep3_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d9_DB_TET3_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d9_DB_TET3_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d9_DB_TET3_all_samples_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d9_DMSO_TET3_rep1" "d9_DMSO_TET3_rep2" "d9_DMSO_TET3_rep3" "d9_BPK25_TET3_rep1" "d9_BPK25_TET3_rep2" "d9_BPK25_TET3_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d9_DMSO_IP_TET3_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_TET3_rep123_merged_TvsC.bigwig -R d49_IP_CTCF.d4_only.summit.bed d49_IP_CTCF.d49_common.summit.bed d49_IP_CTCF.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d9_DB_TET3_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d9_DB_TET3_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.gz -o bamCompare_matrix_of_d9_DB_TET3_rep123_merged_for_plotHeatMap.WT_d49_CTCF_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "CTCF_WT_d4_only(19095)" "CTCF_WT_d49_common(1547)" "CTCF_WT_d9_only(7511)" --samplesLabel "d9_DMSO_TET3_rep123" "d9_BPK25_TET3_rep123"

```



# heatmap_of_reads_enrichment_around_H9_WT_d49_MBD3_diff_summits_by_intersect

## d0_DMSO_BPK25_CTCF

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d0_DB_CTCF_use_WT_d49_diff_MBD3_summits
cd reads_heatmap_of_d0_DB_CTCF_use_WT_d49_diff_MBD3_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d4_only.summit.bed d49_IP_MBD3.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d49_common.summit.bed d49_IP_MBD3.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d9_only.summit.bed d49_IP_MBD3.d9_only.summit.bed
#
ls -1
wc -l *.bed
  # 53190 d49_IP_MBD3.d4_only.summit.bed
  #  3179 d49_IP_MBD3.d49_common.summit.bed
  # 27986 d49_IP_MBD3.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d0_DB_CTCF_use_WT_d49_diff_MBD3_summits
# each repeat
ln -s ../plotHeatMap_use_CTCF_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d0_DMSO_IP_CTCF_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_CTCF_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_CTCF_d0_BPK25_vs_DMSO_diff_summit/H9_d0_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig H9_d0_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d0_BPK25_vs_DMSO_diff_summit/H9_d0_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 CTCF
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d0_DMSO_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_CTCF_rep3_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d0_DB_CTCF_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d0_DB_CTCF_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d0_DB_CTCF_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d0_DMSO_CTCF_rep1" "d0_DMSO_CTCF_rep2" "d0_DMSO_CTCF_rep3" "d0_BPK25_CTCF_rep1" "d0_BPK25_CTCF_rep2" "d0_BPK25_CTCF_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d0_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d0_DB_CTCF_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d0_DB_CTCF_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d0_DB_CTCF_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d0_DMSO_CTCF_rep123" "d0_BPK25_CTCF_rep123"

```

## d4_DMSO_BPK25_CTCF

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d4_DB_CTCF_use_WT_d49_diff_MBD3_summits
cd reads_heatmap_of_d4_DB_CTCF_use_WT_d49_diff_MBD3_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d4_only.summit.bed d49_IP_MBD3.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d49_common.summit.bed d49_IP_MBD3.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d9_only.summit.bed d49_IP_MBD3.d9_only.summit.bed
#
ls -1
wc -l *.bed
  # 53190 d49_IP_MBD3.d4_only.summit.bed
  #  3179 d49_IP_MBD3.d49_common.summit.bed
  # 27986 d49_IP_MBD3.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d4_DB_CTCF_use_WT_d49_diff_MBD3_summits
# each repeat
ln -s ../plotHeatMap_use_CTCF_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d4_DMSO_IP_CTCF_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_CTCF_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_CTCF_d4_BPK25_vs_DMSO_diff_summit/H9_d4_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig H9_d4_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d4_BPK25_vs_DMSO_diff_summit/H9_d4_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 CTCF
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d4_DMSO_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_CTCF_rep3_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d4_DB_CTCF_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d4_DB_CTCF_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d4_DB_CTCF_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d4_DMSO_CTCF_rep1" "d4_DMSO_CTCF_rep2" "d4_DMSO_CTCF_rep3" "d4_BPK25_CTCF_rep1" "d4_BPK25_CTCF_rep2" "d4_BPK25_CTCF_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d4_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d4_DB_CTCF_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d4_DB_CTCF_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d4_DB_CTCF_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d4_DMSO_CTCF_rep123" "d4_BPK25_CTCF_rep123"

```

## d9_DMSO_BPK25_CTCF

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d9_DB_CTCF_use_WT_d49_diff_MBD3_summits
cd reads_heatmap_of_d9_DB_CTCF_use_WT_d49_diff_MBD3_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d4_only.summit.bed d49_IP_MBD3.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d49_common.summit.bed d49_IP_MBD3.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d9_only.summit.bed d49_IP_MBD3.d9_only.summit.bed
#
ls -1
wc -l *.bed
  # 53190 d49_IP_MBD3.d4_only.summit.bed
  #  3179 d49_IP_MBD3.d49_common.summit.bed
  # 27986 d49_IP_MBD3.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d9_DB_CTCF_use_WT_d49_diff_MBD3_summits
# each repeat
ln -s ../plotHeatMap_use_CTCF_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d9_DMSO_IP_CTCF_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_CTCF_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_CTCF_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_CTCF_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_CTCF_d9_BPK25_vs_DMSO_diff_summit/H9_d9_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig H9_d9_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_CTCF_d9_BPK25_vs_DMSO_diff_summit/H9_d9_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 CTCF
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d9_DMSO_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_CTCF_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_CTCF_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_CTCF_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_CTCF_rep3_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d9_DB_CTCF_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d9_DB_CTCF_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d9_DB_CTCF_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d9_DMSO_CTCF_rep1" "d9_DMSO_CTCF_rep2" "d9_DMSO_CTCF_rep3" "d9_BPK25_CTCF_rep1" "d9_BPK25_CTCF_rep2" "d9_BPK25_CTCF_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d9_DMSO_IP_CTCF_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_CTCF_rep123_merged_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d9_DB_CTCF_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d9_DB_CTCF_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d9_DB_CTCF_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d9_DMSO_CTCF_rep123" "d9_BPK25_CTCF_rep123"

```

## d0_DMSO_BPK25_MBD3

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d0_DB_MBD3_use_WT_d49_diff_MBD3_summits
cd reads_heatmap_of_d0_DB_MBD3_use_WT_d49_diff_MBD3_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d4_only.summit.bed d49_IP_MBD3.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d49_common.summit.bed d49_IP_MBD3.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d9_only.summit.bed d49_IP_MBD3.d9_only.summit.bed
#
ls -1
wc -l *.bed
  # 53190 d49_IP_MBD3.d4_only.summit.bed
  #  3179 d49_IP_MBD3.d49_common.summit.bed
  # 27986 d49_IP_MBD3.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d0_DB_MBD3_use_WT_d49_diff_MBD3_summits
# each repeat
ln -s ../plotHeatMap_use_MBD3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_MBD3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_MBD3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d0_DMSO_IP_MBD3_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_MBD3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_MBD3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_MBD3_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_MBD3_d0_BPK25_vs_DMSO_diff_summit/H9_d0_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig H9_d0_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d0_BPK25_vs_DMSO_diff_summit/H9_d0_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 MBD3
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d0_DMSO_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_MBD3_rep3_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d0_DB_MBD3_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d0_DB_MBD3_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d0_DB_MBD3_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d0_DMSO_MBD3_rep1" "d0_DMSO_MBD3_rep2" "d0_DMSO_MBD3_rep3" "d0_BPK25_MBD3_rep1" "d0_BPK25_MBD3_rep2" "d0_BPK25_MBD3_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d0_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d0_DB_MBD3_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d0_DB_MBD3_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d0_DB_MBD3_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d0_DMSO_MBD3_rep123" "d0_BPK25_MBD3_rep123"

```

## d4_DMSO_BPK25_MBD3

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d4_DB_MBD3_use_WT_d49_diff_MBD3_summits
cd reads_heatmap_of_d4_DB_MBD3_use_WT_d49_diff_MBD3_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d4_only.summit.bed d49_IP_MBD3.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d49_common.summit.bed d49_IP_MBD3.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d9_only.summit.bed d49_IP_MBD3.d9_only.summit.bed
#
ls -1
wc -l *.bed
  # 53190 d49_IP_MBD3.d4_only.summit.bed
  #  3179 d49_IP_MBD3.d49_common.summit.bed
  # 27986 d49_IP_MBD3.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d4_DB_MBD3_use_WT_d49_diff_MBD3_summits
# each repeat
ln -s ../plotHeatMap_use_MBD3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_MBD3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_MBD3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d4_DMSO_IP_MBD3_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_MBD3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_MBD3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_MBD3_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_MBD3_d4_BPK25_vs_DMSO_diff_summit/H9_d4_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig H9_d4_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d4_BPK25_vs_DMSO_diff_summit/H9_d4_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 MBD3
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d4_DMSO_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_MBD3_rep3_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d4_DB_MBD3_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d4_DB_MBD3_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d4_DB_MBD3_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d4_DMSO_MBD3_rep1" "d4_DMSO_MBD3_rep2" "d4_DMSO_MBD3_rep3" "d4_BPK25_MBD3_rep1" "d4_BPK25_MBD3_rep2" "d4_BPK25_MBD3_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d4_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d4_DB_MBD3_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d4_DB_MBD3_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d4_DB_MBD3_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d4_DMSO_MBD3_rep123" "d4_BPK25_MBD3_rep123"

```

## d9_DMSO_BPK25_MBD3

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d9_DB_MBD3_use_WT_d49_diff_MBD3_summits
cd reads_heatmap_of_d9_DB_MBD3_use_WT_d49_diff_MBD3_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d4_only.summit.bed d49_IP_MBD3.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d49_common.summit.bed d49_IP_MBD3.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d9_only.summit.bed d49_IP_MBD3.d9_only.summit.bed
#
ls -1
wc -l *.bed
  # 53190 d49_IP_MBD3.d4_only.summit.bed
  #  3179 d49_IP_MBD3.d49_common.summit.bed
  # 27986 d49_IP_MBD3.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d9_DB_MBD3_use_WT_d49_diff_MBD3_summits
# each repeat
ln -s ../plotHeatMap_use_MBD3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_MBD3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_MBD3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d9_DMSO_IP_MBD3_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_MBD3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_MBD3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_MBD3_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_MBD3_d9_BPK25_vs_DMSO_diff_summit/H9_d9_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig H9_d9_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_MBD3_d9_BPK25_vs_DMSO_diff_summit/H9_d9_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 MBD3
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d9_DMSO_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_MBD3_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_MBD3_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_MBD3_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_MBD3_rep3_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d9_DB_MBD3_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d9_DB_MBD3_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d9_DB_MBD3_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d9_DMSO_MBD3_rep1" "d9_DMSO_MBD3_rep2" "d9_DMSO_MBD3_rep3" "d9_BPK25_MBD3_rep1" "d9_BPK25_MBD3_rep2" "d9_BPK25_MBD3_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d9_DMSO_IP_MBD3_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_MBD3_rep123_merged_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d9_DB_MBD3_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d9_DB_MBD3_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d9_DB_MBD3_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d9_DMSO_MBD3_rep123" "d9_BPK25_MBD3_rep123"

```

## d0_DMSO_BPK25_TET1

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d0_DB_TET1_use_WT_d49_diff_MBD3_summits
cd reads_heatmap_of_d0_DB_TET1_use_WT_d49_diff_MBD3_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d4_only.summit.bed d49_IP_MBD3.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d49_common.summit.bed d49_IP_MBD3.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d9_only.summit.bed d49_IP_MBD3.d9_only.summit.bed
#
ls -1
wc -l *.bed
  # 53190 d49_IP_MBD3.d4_only.summit.bed
  #  3179 d49_IP_MBD3.d49_common.summit.bed
  # 27986 d49_IP_MBD3.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d0_DB_TET1_use_WT_d49_diff_MBD3_summits
# each repeat
ln -s ../plotHeatMap_use_TET1_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_TET1_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET1_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_TET1_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET1_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_TET1_rep3_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET1_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_TET1_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET1_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_TET1_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET1_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_TET1_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET1_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_TET1_d0_BPK25_vs_DMSO_diff_summit/H9_d0_DMSO_IP_TET1_rep123_merged_TvsC.bigwig H9_d0_DMSO_IP_TET1_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d0_BPK25_vs_DMSO_diff_summit/H9_d0_BPK25_IP_TET1_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_TET1_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 TET1
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d0_DMSO_IP_TET1_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET1_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET1_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET1_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET1_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET1_rep3_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d0_DB_TET1_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d0_DB_TET1_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d0_DB_TET1_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d0_DMSO_TET1_rep1" "d0_DMSO_TET1_rep2" "d0_DMSO_TET1_rep3" "d0_BPK25_TET1_rep1" "d0_BPK25_TET1_rep2" "d0_BPK25_TET1_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d0_DMSO_IP_TET1_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_TET1_rep123_merged_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d0_DB_TET1_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d0_DB_TET1_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d0_DB_TET1_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d0_DMSO_TET1_rep123" "d0_BPK25_TET1_rep123"

```

## d4_DMSO_BPK25_TET1

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d4_DB_TET1_use_WT_d49_diff_MBD3_summits
cd reads_heatmap_of_d4_DB_TET1_use_WT_d49_diff_MBD3_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d4_only.summit.bed d49_IP_MBD3.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d49_common.summit.bed d49_IP_MBD3.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d9_only.summit.bed d49_IP_MBD3.d9_only.summit.bed
#
ls -1
wc -l *.bed
  # 53190 d49_IP_MBD3.d4_only.summit.bed
  #  3179 d49_IP_MBD3.d49_common.summit.bed
  # 27986 d49_IP_MBD3.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d4_DB_TET1_use_WT_d49_diff_MBD3_summits
# each repeat
ln -s ../plotHeatMap_use_TET1_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_TET1_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET1_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_TET1_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET1_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_TET1_rep3_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET1_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_TET1_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET1_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_TET1_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET1_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_TET1_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET1_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_TET1_d4_BPK25_vs_DMSO_diff_summit/H9_d4_DMSO_IP_TET1_rep123_merged_TvsC.bigwig H9_d4_DMSO_IP_TET1_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d4_BPK25_vs_DMSO_diff_summit/H9_d4_BPK25_IP_TET1_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_TET1_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 TET1
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d4_DMSO_IP_TET1_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET1_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET1_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET1_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET1_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET1_rep3_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d4_DB_TET1_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d4_DB_TET1_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d4_DB_TET1_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d4_DMSO_TET1_rep1" "d4_DMSO_TET1_rep2" "d4_DMSO_TET1_rep3" "d4_BPK25_TET1_rep1" "d4_BPK25_TET1_rep2" "d4_BPK25_TET1_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d4_DMSO_IP_TET1_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_TET1_rep123_merged_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d4_DB_TET1_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d4_DB_TET1_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d4_DB_TET1_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d4_DMSO_TET1_rep123" "d4_BPK25_TET1_rep123"

```

## d9_DMSO_BPK25_TET1

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d9_DB_TET1_use_WT_d49_diff_MBD3_summits
cd reads_heatmap_of_d9_DB_TET1_use_WT_d49_diff_MBD3_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d4_only.summit.bed d49_IP_MBD3.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d49_common.summit.bed d49_IP_MBD3.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d9_only.summit.bed d49_IP_MBD3.d9_only.summit.bed
#
ls -1
wc -l *.bed
  # 53190 d49_IP_MBD3.d4_only.summit.bed
  #  3179 d49_IP_MBD3.d49_common.summit.bed
  # 27986 d49_IP_MBD3.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d9_DB_TET1_use_WT_d49_diff_MBD3_summits
# each repeat
ln -s ../plotHeatMap_use_TET1_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_TET1_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET1_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_TET1_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET1_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_TET1_rep3_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET1_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_TET1_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET1_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_TET1_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET1_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_TET1_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET1_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_TET1_d9_BPK25_vs_DMSO_diff_summit/H9_d9_DMSO_IP_TET1_rep123_merged_TvsC.bigwig H9_d9_DMSO_IP_TET1_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_TET1_d9_BPK25_vs_DMSO_diff_summit/H9_d9_BPK25_IP_TET1_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_TET1_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 TET1
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d9_DMSO_IP_TET1_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET1_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET1_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET1_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET1_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET1_rep3_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d9_DB_TET1_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d9_DB_TET1_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d9_DB_TET1_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d9_DMSO_TET1_rep1" "d9_DMSO_TET1_rep2" "d9_DMSO_TET1_rep3" "d9_BPK25_TET1_rep1" "d9_BPK25_TET1_rep2" "d9_BPK25_TET1_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d9_DMSO_IP_TET1_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_TET1_rep123_merged_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d9_DB_TET1_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d9_DB_TET1_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d9_DB_TET1_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d9_DMSO_TET1_rep123" "d9_BPK25_TET1_rep123"

```

## d0_DMSO_BPK25_TET2

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d0_DB_TET2_use_WT_d49_diff_MBD3_summits
cd reads_heatmap_of_d0_DB_TET2_use_WT_d49_diff_MBD3_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d4_only.summit.bed d49_IP_MBD3.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d49_common.summit.bed d49_IP_MBD3.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d9_only.summit.bed d49_IP_MBD3.d9_only.summit.bed
#
ls -1
wc -l *.bed
  # 53190 d49_IP_MBD3.d4_only.summit.bed
  #  3179 d49_IP_MBD3.d49_common.summit.bed
  # 27986 d49_IP_MBD3.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d0_DB_TET2_use_WT_d49_diff_MBD3_summits
# each repeat
ln -s ../plotHeatMap_use_TET2_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_TET2_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET2_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_TET2_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET2_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_TET2_rep3_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET2_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_TET2_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET2_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_TET2_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET2_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_TET2_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET2_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_TET2_d0_BPK25_vs_DMSO_diff_summit/H9_d0_DMSO_IP_TET2_rep123_merged_TvsC.bigwig H9_d0_DMSO_IP_TET2_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d0_BPK25_vs_DMSO_diff_summit/H9_d0_BPK25_IP_TET2_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_TET2_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 TET2
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d0_DMSO_IP_TET2_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET2_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET2_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET2_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET2_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET2_rep3_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d0_DB_TET2_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d0_DB_TET2_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d0_DB_TET2_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d0_DMSO_TET2_rep1" "d0_DMSO_TET2_rep2" "d0_DMSO_TET2_rep3" "d0_BPK25_TET2_rep1" "d0_BPK25_TET2_rep2" "d0_BPK25_TET2_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d0_DMSO_IP_TET2_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_TET2_rep123_merged_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d0_DB_TET2_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d0_DB_TET2_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d0_DB_TET2_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d0_DMSO_TET2_rep123" "d0_BPK25_TET2_rep123"

```

## d4_DMSO_BPK25_TET2

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d4_DB_TET2_use_WT_d49_diff_MBD3_summits
cd reads_heatmap_of_d4_DB_TET2_use_WT_d49_diff_MBD3_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d4_only.summit.bed d49_IP_MBD3.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d49_common.summit.bed d49_IP_MBD3.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d9_only.summit.bed d49_IP_MBD3.d9_only.summit.bed
#
ls -1
wc -l *.bed
  # 53190 d49_IP_MBD3.d4_only.summit.bed
  #  3179 d49_IP_MBD3.d49_common.summit.bed
  # 27986 d49_IP_MBD3.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d4_DB_TET2_use_WT_d49_diff_MBD3_summits
# each repeat
ln -s ../plotHeatMap_use_TET2_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_TET2_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET2_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_TET2_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET2_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_TET2_rep3_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET2_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_TET2_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET2_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_TET2_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET2_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_TET2_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET2_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_TET2_d4_BPK25_vs_DMSO_diff_summit/H9_d4_DMSO_IP_TET2_rep123_merged_TvsC.bigwig H9_d4_DMSO_IP_TET2_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d4_BPK25_vs_DMSO_diff_summit/H9_d4_BPK25_IP_TET2_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_TET2_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 TET2
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d4_DMSO_IP_TET2_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET2_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET2_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET2_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET2_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET2_rep3_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d4_DB_TET2_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d4_DB_TET2_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d4_DB_TET2_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d4_DMSO_TET2_rep1" "d4_DMSO_TET2_rep2" "d4_DMSO_TET2_rep3" "d4_BPK25_TET2_rep1" "d4_BPK25_TET2_rep2" "d4_BPK25_TET2_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d4_DMSO_IP_TET2_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_TET2_rep123_merged_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d4_DB_TET2_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d4_DB_TET2_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d4_DB_TET2_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d4_DMSO_TET2_rep123" "d4_BPK25_TET2_rep123"

```

## d9_DMSO_BPK25_TET2

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d9_DB_TET2_use_WT_d49_diff_MBD3_summits
cd reads_heatmap_of_d9_DB_TET2_use_WT_d49_diff_MBD3_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d4_only.summit.bed d49_IP_MBD3.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d49_common.summit.bed d49_IP_MBD3.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d9_only.summit.bed d49_IP_MBD3.d9_only.summit.bed
#
ls -1
wc -l *.bed
  # 53190 d49_IP_MBD3.d4_only.summit.bed
  #  3179 d49_IP_MBD3.d49_common.summit.bed
  # 27986 d49_IP_MBD3.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d9_DB_TET2_use_WT_d49_diff_MBD3_summits
# each repeat
ln -s ../plotHeatMap_use_TET2_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_TET2_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET2_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_TET2_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET2_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_TET2_rep3_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET2_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_TET2_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET2_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_TET2_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET2_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_TET2_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET2_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_TET2_d9_BPK25_vs_DMSO_diff_summit/H9_d9_DMSO_IP_TET2_rep123_merged_TvsC.bigwig H9_d9_DMSO_IP_TET2_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_TET2_d9_BPK25_vs_DMSO_diff_summit/H9_d9_BPK25_IP_TET2_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_TET2_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 TET2
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d9_DMSO_IP_TET2_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET2_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET2_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET2_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET2_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET2_rep3_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d9_DB_TET2_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d9_DB_TET2_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d9_DB_TET2_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d9_DMSO_TET2_rep1" "d9_DMSO_TET2_rep2" "d9_DMSO_TET2_rep3" "d9_BPK25_TET2_rep1" "d9_BPK25_TET2_rep2" "d9_BPK25_TET2_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d9_DMSO_IP_TET2_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_TET2_rep123_merged_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d9_DB_TET2_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d9_DB_TET2_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d9_DB_TET2_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d9_DMSO_TET2_rep123" "d9_BPK25_TET2_rep123"

```

## d0_DMSO_BPK25_TET3

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d0_DB_TET3_use_WT_d49_diff_MBD3_summits
cd reads_heatmap_of_d0_DB_TET3_use_WT_d49_diff_MBD3_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d4_only.summit.bed d49_IP_MBD3.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d49_common.summit.bed d49_IP_MBD3.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d9_only.summit.bed d49_IP_MBD3.d9_only.summit.bed
#
ls -1
wc -l *.bed
  # 53190 d49_IP_MBD3.d4_only.summit.bed
  #  3179 d49_IP_MBD3.d49_common.summit.bed
  # 27986 d49_IP_MBD3.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d0_DB_TET3_use_WT_d49_diff_MBD3_summits
# each repeat
ln -s ../plotHeatMap_use_TET3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_TET3_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_TET3_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_DMSO_IP_TET3_rep3_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET3_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_TET3_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_TET3_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d0_BPK25_vs_DMSO_diff_summit/hESC_H9_d0_BPK25_IP_TET3_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET3_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_TET3_d0_BPK25_vs_DMSO_diff_summit/H9_d0_DMSO_IP_TET3_rep123_merged_TvsC.bigwig H9_d0_DMSO_IP_TET3_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d0_BPK25_vs_DMSO_diff_summit/H9_d0_BPK25_IP_TET3_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_TET3_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 TET3
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d0_DMSO_IP_TET3_rep1_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET3_rep2_TvsC.bigwig hESC_H9_d0_DMSO_IP_TET3_rep3_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET3_rep1_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET3_rep2_TvsC.bigwig hESC_H9_d0_BPK25_IP_TET3_rep3_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d0_DB_TET3_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d0_DB_TET3_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d0_DB_TET3_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d0_DMSO_TET3_rep1" "d0_DMSO_TET3_rep2" "d0_DMSO_TET3_rep3" "d0_BPK25_TET3_rep1" "d0_BPK25_TET3_rep2" "d0_BPK25_TET3_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d0_DMSO_IP_TET3_rep123_merged_TvsC.bigwig H9_d0_BPK25_IP_TET3_rep123_merged_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d0_DB_TET3_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d0_DB_TET3_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d0_DB_TET3_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d0_DMSO_TET3_rep123" "d0_BPK25_TET3_rep123"

```

## d4_DMSO_BPK25_TET3

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d4_DB_TET3_use_WT_d49_diff_MBD3_summits
cd reads_heatmap_of_d4_DB_TET3_use_WT_d49_diff_MBD3_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d4_only.summit.bed d49_IP_MBD3.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d49_common.summit.bed d49_IP_MBD3.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d9_only.summit.bed d49_IP_MBD3.d9_only.summit.bed
#
ls -1
wc -l *.bed
  # 53190 d49_IP_MBD3.d4_only.summit.bed
  #  3179 d49_IP_MBD3.d49_common.summit.bed
  # 27986 d49_IP_MBD3.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d4_DB_TET3_use_WT_d49_diff_MBD3_summits
# each repeat
ln -s ../plotHeatMap_use_TET3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_TET3_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_TET3_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_DMSO_IP_TET3_rep3_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET3_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_TET3_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_TET3_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d4_BPK25_vs_DMSO_diff_summit/hESC_H9_d4_BPK25_IP_TET3_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET3_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_TET3_d4_BPK25_vs_DMSO_diff_summit/H9_d4_DMSO_IP_TET3_rep123_merged_TvsC.bigwig H9_d4_DMSO_IP_TET3_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d4_BPK25_vs_DMSO_diff_summit/H9_d4_BPK25_IP_TET3_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_TET3_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 TET3
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d4_DMSO_IP_TET3_rep1_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET3_rep2_TvsC.bigwig hESC_H9_d4_DMSO_IP_TET3_rep3_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET3_rep1_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET3_rep2_TvsC.bigwig hESC_H9_d4_BPK25_IP_TET3_rep3_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d4_DB_TET3_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d4_DB_TET3_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d4_DB_TET3_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d4_DMSO_TET3_rep1" "d4_DMSO_TET3_rep2" "d4_DMSO_TET3_rep3" "d4_BPK25_TET3_rep1" "d4_BPK25_TET3_rep2" "d4_BPK25_TET3_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d4_DMSO_IP_TET3_rep123_merged_TvsC.bigwig H9_d4_BPK25_IP_TET3_rep123_merged_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d4_DB_TET3_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d4_DB_TET3_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d4_DB_TET3_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d4_DMSO_TET3_rep123" "d4_BPK25_TET3_rep123"

```

## d9_DMSO_BPK25_TET3

### data preperation

```bash
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input
mkdir reads_heatmap_of_d9_DB_TET3_use_WT_d49_diff_MBD3_summits
cd reads_heatmap_of_d9_DB_TET3_use_WT_d49_diff_MBD3_summits
# 创建 d49 diff summit 文件连接 (由于 bdgdiff的差异summit基本都在intersect的结果中，后续以 intersect 的结果为准 )
# 创建 intersect 三类 summit 文件连接
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d4_only.summit.bed d49_IP_MBD3.d4_only.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d49_common.summit.bed d49_IP_MBD3.d49_common.summit.bed
ln -s /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/02_WN_H9_d049_WT_ChIP_Seq_with_added_input/macs2_peak_calling_2_for_H9_d049_WT_MBD3_rep123_together/d49_IP_MBD3_rep123_TvsC.d9_only.summit.bed d49_IP_MBD3.d9_only.summit.bed
#
ls -1
wc -l *.bed
  # 53190 d49_IP_MBD3.d4_only.summit.bed
  #  3179 d49_IP_MBD3.d49_common.summit.bed
  # 27986 d49_IP_MBD3.d9_only.summit.bed

```

### use bamCompare  to create log2 transformed bigwig files
```bash
# use bamCompare  to create log2 transformed bigwig files
# 在 7.1 中均已计算了，直接构造软连接
cd /home/sunwenju/projects_on_940xa/2022081701_WN_ChIP_Seq_data/03_WN_H9_d049_BPK25_ChIP_Seq_with_added_input/reads_heatmap_of_d9_DB_TET3_use_WT_d49_diff_MBD3_summits
# each repeat
ln -s ../plotHeatMap_use_TET3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_TET3_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_TET3_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_DMSO_IP_TET3_rep3_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET3_rep3_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_TET3_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET3_rep1_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_TET3_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET3_rep2_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d9_BPK25_vs_DMSO_diff_summit/hESC_H9_d9_BPK25_IP_TET3_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET3_rep3_TvsC.bigwig
# rep123 merged
ln -s ../plotHeatMap_use_TET3_d9_BPK25_vs_DMSO_diff_summit/H9_d9_DMSO_IP_TET3_rep123_merged_TvsC.bigwig H9_d9_DMSO_IP_TET3_rep123_merged_TvsC.bigwig
ln -s ../plotHeatMap_use_TET3_d9_BPK25_vs_DMSO_diff_summit/H9_d9_BPK25_IP_TET3_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_TET3_rep123_merged_TvsC.bigwig
ls -1

```

### computeMatrix and plotHeatmap
```bash
# plot for d0 DMSO and BPK25 TET3
# use individual bam of each repeat
computeMatrix reference-point -p max --referencePoint center -S hESC_H9_d9_DMSO_IP_TET3_rep1_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET3_rep2_TvsC.bigwig hESC_H9_d9_DMSO_IP_TET3_rep3_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET3_rep1_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET3_rep2_TvsC.bigwig hESC_H9_d9_BPK25_IP_TET3_rep3_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d9_DB_TET3_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d9_DB_TET3_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d9_DB_TET3_all_samples_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d9_DMSO_TET3_rep1" "d9_DMSO_TET3_rep2" "d9_DMSO_TET3_rep3" "d9_BPK25_TET3_rep1" "d9_BPK25_TET3_rep2" "d9_BPK25_TET3_rep3"

# use rep123 bam together in one
computeMatrix reference-point -p max --referencePoint center -S H9_d9_DMSO_IP_TET3_rep123_merged_TvsC.bigwig H9_d9_BPK25_IP_TET3_rep123_merged_TvsC.bigwig -R d49_IP_MBD3.d4_only.summit.bed d49_IP_MBD3.d49_common.summit.bed d49_IP_MBD3.d9_only.summit.bed -b 1000 -a 1000 --sortRegions descend --missingDataAsZero --skipZeros --blackListFileName /home/sunwenju/3_genomeData/human/hg19/encode_hg19_black_list.sorted.bed -o bamCompare_matrix_of_d9_DB_TET3_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz
plotHeatmap -m bamCompare_matrix_of_d9_DB_TET3_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.gz -o bamCompare_matrix_of_d9_DB_TET3_rep123_merged_for_plotHeatMap.WT_d49_MBD3_diff_summits.plot.pdf --heatmapHeight 18 --heatmapWidth 5 --whatToShow "plot, heatmap and colorbar" --colorMap 'coolwarm' --refPointLabel "0" --xAxisLabel "distance to summits" --regionsLabel "MBD3_WT_d4_only(53190)" "MBD3_WT_d49_common(3179)" "MBD3_WT_d9_only(27986)" --samplesLabel "d9_DMSO_TET3_rep123" "d9_BPK25_TET3_rep123"

```















