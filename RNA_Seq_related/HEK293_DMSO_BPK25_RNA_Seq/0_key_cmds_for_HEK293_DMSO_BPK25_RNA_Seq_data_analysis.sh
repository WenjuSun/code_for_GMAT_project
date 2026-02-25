

bash 0_make_and_rename_HEK293_BPK25_RNA_Seq_raw_fq_gz_file_links.sh
ls -1
ls -1 *.fq.gz |wc -l

# 
mkdir 1.1_all_HEK293_RNA_Seq_raw_data_fastqc_result
nohup fastqc -outdir 1.1_all_HEK293_RNA_Seq_raw_data_fastqc_result --threads 128 HEK293*R[12].fq.gz >1.1_all_HEK293_RNA_Seq_raw_data_fastqc.runlog 2>&1 &
tail -f 1.1_all_HEK293_RNA_Seq_raw_data_fastqc.runlog

#
cd 1.1_all_HEK293_RNA_Seq_raw_data_fastqc_result
multiqc ./*fastqc.zip --ignore *.html --no-megaqc-upload --interactive
cd ../


for i in `ls -1 HEK293*_RNA_Seq_rep*.R[12].fq.gz|sort`;do echo $i; zcat $i|grep -E '^[ACGTN]{15,}$' |wc -l; zcat $i |grep -E '^[ACGTN]{15,}$' |grep "AGATCGGAAGAGC" |wc -l;done
for i in `ls -1 HEK293*_RNA_Seq_rep*.R1.fq.gz|sort`;do echo 'cutadapt -j 128 -u 13 -U -13 -q 20,20 --poly-a --trim-n -m 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -a "G{20};min_overlap=12" -A "G{20};min_overlap=12" -o '${i%.R1.fq.gz}.R1.cutadapt.fq.gz' -p '${i%.R1.fq.gz}.R2.cutadapt.fq.gz' '$i' '${i%.R1.fq.gz}.R2.fq.gz;done >1.2_all_HEK293_RNA_Seq_raw_data_run_cutadapt.sh
head 1.2_all_HEK293_RNA_Seq_raw_data_run_cutadapt.sh
nohup bash 1.2_all_HEK293_RNA_Seq_raw_data_run_cutadapt.sh >1.2_all_HEK293_RNA_Seq_raw_data_run_cutadapt.sh.runlog 2>&1 &
tail -f 1.2_all_HEK293_RNA_Seq_raw_data_run_cutadapt.sh.runlog


mkdir 1.3_all_HEK293_RNA_Seq_cutadapt_fastqc_result
nohup fastqc -outdir 1.3_all_HEK293_RNA_Seq_cutadapt_fastqc_result --threads 128 HEK293_*.cutadapt.fq.gz >1.3_all_HEK293_RNA_Seq_cutadapt_fastqc.runlog 2>&1 &
tail -f 1.3_all_HEK293_RNA_Seq_cutadapt_fastqc.runlog

#
cd 1.3_all_HEK293_RNA_Seq_cutadapt_fastqc_result
multiqc ./*fastqc.zip --ignore *.html --no-megaqc-upload --interactive
cd ../


for i in `ls -1 HEK293*.R1.cutadapt.fq.gz|sort`;do echo 'hisat2 -q -t --threads 64 --un-gz '${i%.R1.cutadapt.fq.gz}.notAlign.unpairedReads.fq.gz' --un-conc-gz '${i%.R1.cutadapt.fq.gz}.notAlign.pairedReads.fq.gz' --new-summary --summary-file '${i%.R1.cutadapt.fq.gz}.hisat2.mapping.summary.txt' -x /home/sunwenju/3_genomeData/human/hg19/hisat2Index/hg19hisat2idx -1 '$i' -2 '${i%.R1.cutadapt.fq.gz}.R2.cutadapt.fq.gz' |samtools view -S -b -h -F 4 -@ 32 |samtools sort -@ 32 -o '${i%.R1.cutadapt.fq.gz}.hisat2.sorted.bam; echo 'samtools index -@ 32 '${i%.R1.cutadapt.fq.gz}.hisat2.sorted.bam;done >2_all_HEK293_RNA_Seq_data_do_hisat2_map2.sh
head 2_all_HEK293_RNA_Seq_data_do_hisat2_map2.sh
nohup bash 2_all_HEK293_RNA_Seq_data_do_hisat2_map2.sh >2_all_HEK293_RNA_Seq_data_do_hisat2_map2.sh.runlog 2>&1 &
tail -f 2_all_HEK293_RNA_Seq_data_do_hisat2_map2.sh.runlog

nohup featureCounts -F GTF -t exon -g gene_id --extraAttributes gene_name,gene_type -O -s 0 -p -C -T 36 -a /home/sunwenju/3_genomeData/human/hg19/annotations/gencode.v19.annotation.gtf -o 0_all_HEK293_BPK25_RNA_Seq_samples.featureCounts.hg19.exonO.s0.txt HEK293_WT_DMSO_RNA_Seq_rep1.hisat2.sorted.bam HEK293_WT_DMSO_RNA_Seq_rep2.hisat2.sorted.bam HEK293_WT_DMSO_RNA_Seq_rep3.hisat2.sorted.bam HEK293_WT_BPK25_6h_RNA_Seq_rep1.hisat2.sorted.bam HEK293_WT_BPK25_6h_RNA_Seq_rep2.hisat2.sorted.bam HEK293_WT_BPK25_6h_RNA_Seq_rep3.hisat2.sorted.bam HEK293_WT_BPK25_6h_RNA_Seq_rep4.hisat2.sorted.bam HEK293_WT_BPK25_6h_RNA_Seq_rep5.hisat2.sorted.bam HEK293_WT_BPK25_6h_RNA_Seq_rep6.hisat2.sorted.bam HEK293_WT_BPK25_12h_RNA_Seq_rep1.hisat2.sorted.bam HEK293_WT_BPK25_12h_RNA_Seq_rep2.hisat2.sorted.bam HEK293_WT_BPK25_12h_RNA_Seq_rep3.hisat2.sorted.bam HEK293_WT_BPK25_12h_RNA_Seq_rep4.hisat2.sorted.bam HEK293_WT_BPK25_12h_RNA_Seq_rep5.hisat2.sorted.bam HEK293_WT_BPK25_12h_RNA_Seq_rep6.hisat2.sorted.bam >0_all_HEK293_BPK25_RNA_Seq_samples.featureCounts.hg19.exonO.s0.runlog 2>&1 &
tail -f 0_all_HEK293_BPK25_RNA_Seq_samples.featureCounts.hg19.exonO.s0.runlog


# R
# DESeq2 DEGS / GSEA / GO + KEGG 
# 0_get_ge_and_degs_used_fCt_and_DESeq2.c1_BPK25_6h_vs_DMSO.r
# 0_get_ge_and_degs_used_fCt_and_DESeq2.c2_BPK25_12h_vs_DMSO.r
# 0_get_ge_and_degs_used_fCt_and_DESeq2.c3_BPK25_12h_vs_6h.r
