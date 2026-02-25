

bash 0_make_and_rename_WN_hESC_H9_RNA_Seq_fq_file_links.sh
ls -1
ls -1 *.fq.gz |wc -l

mkdir 1_all_WN_H9_RNA_Seq_raw_data_fastqc_result
nohup fastqc -outdir 1_all_WN_H9_RNA_Seq_raw_data_fastqc_result --threads 72 *.fq.gz >1_all_WN_H9_RNA_Seq_raw_data_fastqc.runlog 2>&1 &
tail -f 1_all_WN_H9_RNA_Seq_raw_data_fastqc.runlog

# 
cd 1_all_WN_H9_RNA_Seq_raw_data_fastqc_result
multiqc ./*fastqc.zip --ignore *.html --no-megaqc-upload
cd ../


for i in `ls -1 *R1.fq.gz`;do echo 'hisat2 -q -t --threads 72 --un-gz '${i%.R1.fq.gz}.notAlign.unpairedReads.fq.gz' --un-conc-gz '${i%.R1.fq.gz}.notAlign.pairedReads.fq.gz' --new-summary --summary-file '${i%.R1.fq.gz}.hisat2.mapping.summary.txt' -x /home/sunwenju/3_genomeData/human/hg19/hisat2Index/hg19hisat2idx -1 '$i' -2 '${i%.R1.fq.gz}.R2.fq.gz' -S '${i%.R1.fq.gz}.hisat2.sam;done >2_all_WN_hESC_H9_RNA_Seq_data_do_hisat2_map2.sh
nohup bash 2_all_WN_hESC_H9_RNA_Seq_data_do_hisat2_map2.sh >2_all_WN_hESC_H9_RNA_Seq_data_do_hisat2_map2.sh.runlog 2>&1 &
tail -f 2_all_WN_hESC_H9_RNA_Seq_data_do_hisat2_map2.sh.runlog


nohup featureCounts -F GTF -t exon -g gene_id --extraAttributes gene_name,gene_type -O -s 1 -p -C -T 36 -a /home/sunwenju/3_genomeData/human/hg19/annotations/gencode.v19.annotation.gtf -o 0_WN_hESC_H9_RNA_Seq_samples.featureCounts.s1.txt hESC_H9_0day_RNA_Seq_rep1.hisat2.sorted.bam hESC_H9_0day_RNA_Seq_rep2.hisat2.sorted.bam hESC_H9_0day_RNA_Seq_rep3.hisat2.sorted.bam hESC_H9_4day_RNA_Seq_rep1.hisat2.sorted.bam hESC_H9_4day_RNA_Seq_rep2.hisat2.sorted.bam hESC_H9_4day_RNA_Seq_rep3.hisat2.sorted.bam hESC_H9_9day_RNA_Seq_rep1.hisat2.sorted.bam hESC_H9_9day_RNA_Seq_rep2.hisat2.sorted.bam hESC_H9_9day_RNA_Seq_rep3.hisat2.sorted.bam hESC_H9_14day_RNA_Seq_rep1.hisat2.sorted.bam hESC_H9_14day_RNA_Seq_rep2.hisat2.sorted.bam hESC_H9_14day_RNA_Seq_rep3.hisat2.sorted.bam >0_WN_hESC_H9_RNA_Seq_samples.featureCounts.s1.runlog 2>&1 &
tail -f 0_WN_hESC_H9_RNA_Seq_samples.featureCounts.s1.runlog

cut -f 1,7-20 0_WN_hESC_H9_RNA_Seq_samples.featureCounts.s1.txt >0_WN_hESC_H9_RNA_Seq_samples.featureCounts.s1.cut.txt
# Geneid	gene_name	gene_type	H9_0day_rep1	H9_0day_rep2	H9_0day_rep3	H9_4day_rep1	H9_4day_rep2	H9_4day_rep3	H9_9day_rep1	H9_9day_rep2	H9_9day_rep3	H9_14day_rep1	H9_14day_rep2	H9_14day_rep3

# R
# DESeq2 DEGS / GSEA / GO + KEGG 
# 1_get_ge_and_degs_used_featureCounts_and_DESeq2.Comparison1_d4_VS_d0.r
# 1_get_ge_and_degs_used_featureCounts_and_DESeq2.Comparison2_d9_VS_d0.r
# 1_get_ge_and_degs_used_featureCounts_and_DESeq2.Comparison3_d14_VS_d0.r
# 1_get_ge_and_degs_used_featureCounts_and_DESeq2.Comparison4_d9_VS_d4.r
# 1_get_ge_and_degs_used_featureCounts_and_DESeq2.Comparison5_d14_VS_d4.r
# 1_get_ge_and_degs_used_featureCounts_and_DESeq2.Comparison6_d14_VS_d9.r

