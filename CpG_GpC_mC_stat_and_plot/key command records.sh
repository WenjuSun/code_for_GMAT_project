
#
nohup bismark_methylation_extractor 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.sortedn.merged.bam -o 293T_ERT_CTCF_24h_pd_rep12_dedup_merged_bam_bme_gw_CX_zb --no_overlap --no_header --comprehensive --cytosine_report --CX --zero_based --genome_folder /home/sunwenju/projects/20210000_CTCF_CpG_methylation/genome_data_for_Bismark/ >293T_ERT_CTCF_24h_pd_rep12_dedup_merged_bam_bme_gw_CX_zb.runlog 2>&1 &
nohup bismark_methylation_extractor 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.sortedn.merged.bam -o 293T_ERT_CTCF_48h_pd_rep12_dedup_merged_bam_bme_gw_CX_zb --no_overlap --no_header --comprehensive --cytosine_report --CX --zero_based --genome_folder /home/sunwenju/projects/20210000_CTCF_CpG_methylation/genome_data_for_Bismark/ >293T_ERT_CTCF_48h_pd_rep12_dedup_merged_bam_bme_gw_CX_zb.runlog 2>&1 &
# tail -f 293T_ERT_CTCF_24h_pd_rep12_dedup_merged_bam_bme_gw_CX_zb.runlog
# tail -f 293T_ERT_CTCF_48h_pd_rep12_dedup_merged_bam_bme_gw_CX_zb.runlog

# get genome wide CpG and GpC methylation information for each sample:
nohup python3 ../make_CpG_GpC_mC_site_info_from_bme_results.py -i ./293T_ERT_CTCF_24h_pd_rep12_dedup_merged_bam_bme_gw_CX_zb/293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.sortedn.merged.CX_report.txt -m 1 -o ./293T_ERT_CTCF_24h_pd_rep12_dedup_merged_bam_bme_gw_CX_zb/293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged >./293T_ERT_CTCF_24h_pd_rep12_dedup_merged_bam_bme_gw_CX_zb/293T_ERT_CTCF_24h_pd_rep12.bismark.dedup.merged.bme.makeCpG.GpC.runlog 2>&1 &
nohup python3 ../make_CpG_GpC_mC_site_info_from_bme_results.py -i ./293T_ERT_CTCF_48h_pd_rep12_dedup_merged_bam_bme_gw_CX_zb/293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.sortedn.merged.CX_report.txt -m 1 -o ./293T_ERT_CTCF_48h_pd_rep12_dedup_merged_bam_bme_gw_CX_zb/293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged >./293T_ERT_CTCF_48h_pd_rep12_dedup_merged_bam_bme_gw_CX_zb/293T_ERT_CTCF_48h_pd_rep12.bismark.dedup.merged.bme.makeCpG.GpC.runlog 2>&1 &
# tail -f ./293T_ERT_CTCF_24h_pd_rep12_dedup_merged_bam_bme_gw_CX_zb/293T_ERT_CTCF_24h_pd_rep12.bismark.dedup.merged.bme.makeCpG.GpC.runlog
# tail -f ./293T_ERT_CTCF_48h_pd_rep12_dedup_merged_bam_bme_gw_CX_zb/293T_ERT_CTCF_48h_pd_rep12.bismark.dedup.merged.bme.makeCpG.GpC.runlog


# 3.2_make_all_reads_mC_pos_info_and_stat_result.sh：
# for i in `ls -1 *dedup.sortedn.merged.bam`; do echo $i; time samtools view $i -o ${i}.convert.sam; time python3 /home/sunwenju/projects/20210000_CTCF_CpG_methylation/get_CpG_GpC_mC_pos_info_and_stat_for_every_read_from_bam_mapping_results.py -i ${i}.convert.sam -o ${i}.convert.sam.reads_mC.txt;done
nohup bash 3.2_make_all_reads_mC_pos_info_and_stat_result.sh >3.2_make_all_reads_mC_pos_info_and_stat_result.sh.runlog 2>&1 &
tail -f 3.2_make_all_reads_mC_pos_info_and_stat_result.sh.runlog
# chrom, g_start, g_end, new_name, mapq, strand, read_seq, read_str, read_mc, expanded_mdz, xr_tag, xg_tag, md_tag, read_len, cpg_mc_pos_info, gpc_mc_pos_info, all_mC_count, cpg_count, gpc_count, cpg_mc_count, gpc_mc_count, cpg_mc_ratio, gpc_mc_ratio
# cut：chrom  star    end strand  read_len    all_mC_count    CpG_count   GpC_count   CpG_mC_count    GpC_mC_count    CpG_mC_ratio    GpC_mC_ratio
# 4.2_make_all_reads_mC_pos_info_filter2_for_150nt_6bins_stat.sh：
# cut -f 1-3,6,14,17-23 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.sortedn.merged.bam.convert.sam.reads_mC.txt  | awk 'BEGIN{FS="\t";OFS="\t";print "chrom", "star", "end", "strand", "read_len", "all_mC_count", "CpG_count", "GpC_count", "CpG_mC_count", "GpC_mC_count", "CpG_mC_ratio", "GpC_mC_ratio"}{if(($7!=0 && $8!=0)&&($9!=0 || $10!=0)){print $0}}' >293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.sortedn.merged.bam.convert.sam.reads_mC.filter2.txt
# cut -f 1-3,6,14,17-23 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.sortedn.merged.bam.convert.sam.reads_mC.txt  | awk 'BEGIN{FS="\t";OFS="\t";print "chrom", "star", "end", "strand", "read_len", "all_mC_count", "CpG_count", "GpC_count", "CpG_mC_count", "GpC_mC_count", "CpG_mC_ratio", "GpC_mC_ratio"}{if(($7!=0 && $8!=0)&&($9!=0 || $10!=0)){print $0}}' >293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.sortedn.merged.bam.convert.sam.reads_mC.filter2.txt
nohup bash 4.2_make_all_reads_mC_pos_info_filter2_for_150nt_6bins_stat.sh >4.2_make_all_reads_mC_pos_info_filter2_for_150nt_6bins_stat.sh.runlog 2>&1 &
tail -f 4.2_make_all_reads_mC_pos_info_filter2_for_150nt_6bins_stat.sh.runlog

cd /home/sunwenju/projects/20210000_CTCF_CpG_methylation/analysis_based_on_rep_merged_dedup_bam
mkdir 3_all_samples_reads_CpG_GpC_filtered2_sites_for_6bins_stat
mv *.merged.bam.convert.sam.reads_mC.filter2.txt ./3_all_samples_reads_CpG_GpC_filtered2_sites_for_6bins_stat/

# 4.3_make_reads_CpG_GpC_data_for_eB150nt_intersect.sh：
# for i in `ls -1 *bam.convert.sam.reads_mC.txt`;do echo $i; cut -f 1-3,6,7,8,9,11,12,18,19,20,21 $i | awk 'BEGIN{FS="\t";OFS="\t";print "chrom", "star", "end", "strand", "gSeq","read_str", "read_mC_str", "xr_tag", "xg_tag", "CpG_count", "GpC_count", "CpG_mC_count", "GpC_mC_count"}{if(($10!=0 && $11!=0)&&($12!=0 || $13!=0)){print $0}}' >${i%.txt}.filtered2.forIntersect.txt;done
nohup bash 4.3_make_reads_CpG_GpC_data_for_eB150nt_intersect.sh >4.3_make_reads_CpG_GpC_data_for_eB150nt_intersect.sh.runlog 2>&1 &
tail -f 4.3_make_reads_CpG_GpC_data_for_eB150nt_intersect.sh.runlog

cd /home/sunwenju/projects/20210000_CTCF_CpG_methylation/analysis_based_on_rep_merged_dedup_bam
mkdir 4_all_samples_reads_CpG_GpC_filtered2_sites_for_eB150nt_intersect
mv *.merged.bam.convert.sam.reads_mC.filtered2.forIntersect.txt ./4_all_samples_reads_CpG_GpC_filtered2_sites_for_eB150nt_intersect/



# 
nohup python3 make_meta_profile_for_point_frank.py -a CTCF.summit.all.sorted.bed -b 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.txt -c 1 -d 5000 -o 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.flank5000.stat.txt >293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.flank5000.stat.runlog 2>&1 &
nohup python3 make_meta_profile_for_point_frank.py -a CTCF.summit.all.sorted.bed -b 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.txt -c 1 -d 5000 -o 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.flank5000.stat.txt >293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.flank5000.stat.runlog 2>&1 &
nohup python3 make_meta_profile_for_point_frank.py -a CTCF.summit.all.sorted.bed -b 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.txt -c 1 -d 5000 -o 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.flank5000.stat.txt >293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.flank5000.stat.runlog 2>&1 &
nohup python3 make_meta_profile_for_point_frank.py -a CTCF.summit.all.sorted.bed -b 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.txt -c 1 -d 5000 -o 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.flank5000.stat.txt >293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.flank5000.stat.runlog 2>&1 &

cut -f 1,2,3,4,5,6 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.flank5000.stat.txt >293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.flank5000.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.flank5000.stat.txt >293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.flank5000.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.flank5000.stat.txt >293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.flank5000.stat.cut1_6.txt
cut -f 1,2,3,4,5,6 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.flank5000.stat.txt >293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.flank5000.stat.cut1_6.txt


nohup python3 make_meta_profile_for_point_frank.py -a CTCF.summit.peak_height_quartile_1.sorted.bed -b 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.txt -c 1 -d 5000 -o 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq1.flank5k.stat.txt >293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq1.flank5k.stat.runlog 2>&1 &
nohup python3 make_meta_profile_for_point_frank.py -a CTCF.summit.peak_height_quartile_2.sorted.bed -b 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.txt -c 1 -d 5000 -o 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq2.flank5k.stat.txt >293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq2.flank5k.stat.runlog 2>&1 &
nohup python3 make_meta_profile_for_point_frank.py -a CTCF.summit.peak_height_quartile_3.sorted.bed -b 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.txt -c 1 -d 5000 -o 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq3.flank5k.stat.txt >293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq3.flank5k.stat.runlog 2>&1 &
nohup python3 make_meta_profile_for_point_frank.py -a CTCF.summit.peak_height_quartile_4.sorted.bed -b 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.txt -c 1 -d 5000 -o 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq4.flank5k.stat.txt >293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq4.flank5k.stat.runlog 2>&1 &
nohup python3 make_meta_profile_for_point_frank.py -a CTCF.summit.peak_height_quartile_1.sorted.bed -b 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.txt -c 1 -d 5000 -o 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq1.flank5k.stat.txt >293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq1.flank5k.stat.runlog 2>&1 &
nohup python3 make_meta_profile_for_point_frank.py -a CTCF.summit.peak_height_quartile_2.sorted.bed -b 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.txt -c 1 -d 5000 -o 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq2.flank5k.stat.txt >293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq2.flank5k.stat.runlog 2>&1 &
nohup python3 make_meta_profile_for_point_frank.py -a CTCF.summit.peak_height_quartile_3.sorted.bed -b 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.txt -c 1 -d 5000 -o 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq3.flank5k.stat.txt >293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq3.flank5k.stat.runlog 2>&1 &
nohup python3 make_meta_profile_for_point_frank.py -a CTCF.summit.peak_height_quartile_4.sorted.bed -b 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.txt -c 1 -d 5000 -o 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq4.flank5k.stat.txt >293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq4.flank5k.stat.runlog 2>&1 &
nohup python3 make_meta_profile_for_point_frank.py -a CTCF.summit.peak_height_quartile_1.sorted.bed -b 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.txt -c 1 -d 5000 -o 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq1.flank5k.stat.txt >293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq1.flank5k.stat.runlog 2>&1 &
nohup python3 make_meta_profile_for_point_frank.py -a CTCF.summit.peak_height_quartile_2.sorted.bed -b 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.txt -c 1 -d 5000 -o 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq2.flank5k.stat.txt >293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq2.flank5k.stat.runlog 2>&1 &
nohup python3 make_meta_profile_for_point_frank.py -a CTCF.summit.peak_height_quartile_3.sorted.bed -b 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.txt -c 1 -d 5000 -o 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq3.flank5k.stat.txt >293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq3.flank5k.stat.runlog 2>&1 &
nohup python3 make_meta_profile_for_point_frank.py -a CTCF.summit.peak_height_quartile_4.sorted.bed -b 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.txt -c 1 -d 5000 -o 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq4.flank5k.stat.txt >293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq4.flank5k.stat.runlog 2>&1 &
nohup python3 make_meta_profile_for_point_frank.py -a CTCF.summit.peak_height_quartile_1.sorted.bed -b 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.txt -c 1 -d 5000 -o 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq1.flank5k.stat.txt >293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq1.flank5k.stat.runlog 2>&1 &
nohup python3 make_meta_profile_for_point_frank.py -a CTCF.summit.peak_height_quartile_2.sorted.bed -b 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.txt -c 1 -d 5000 -o 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq2.flank5k.stat.txt >293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq2.flank5k.stat.runlog 2>&1 &
nohup python3 make_meta_profile_for_point_frank.py -a CTCF.summit.peak_height_quartile_3.sorted.bed -b 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.txt -c 1 -d 5000 -o 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq3.flank5k.stat.txt >293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq3.flank5k.stat.runlog 2>&1 &
nohup python3 make_meta_profile_for_point_frank.py -a CTCF.summit.peak_height_quartile_4.sorted.bed -b 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.txt -c 1 -d 5000 -o 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq4.flank5k.stat.txt >293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq4.flank5k.stat.runlog 2>&1 &


# cut1_6_for_CTCF_4_class_summit_meta_plot.sh
# cut -f 1,2,3,4,5,6 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq1.flank5k.stat.txt >293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq1.flank5k.stat.cut1_6.txt
# cut -f 1,2,3,4,5,6 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq2.flank5k.stat.txt >293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq2.flank5k.stat.cut1_6.txt
# cut -f 1,2,3,4,5,6 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq3.flank5k.stat.txt >293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq3.flank5k.stat.cut1_6.txt
# cut -f 1,2,3,4,5,6 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq4.flank5k.stat.txt >293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq4.flank5k.stat.cut1_6.txt
# cut -f 1,2,3,4,5,6 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq1.flank5k.stat.txt >293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq1.flank5k.stat.cut1_6.txt
# cut -f 1,2,3,4,5,6 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq2.flank5k.stat.txt >293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq2.flank5k.stat.cut1_6.txt
# cut -f 1,2,3,4,5,6 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq3.flank5k.stat.txt >293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq3.flank5k.stat.cut1_6.txt
# cut -f 1,2,3,4,5,6 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq4.flank5k.stat.txt >293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq4.flank5k.stat.cut1_6.txt
# cut -f 1,2,3,4,5,6 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq1.flank5k.stat.txt >293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq1.flank5k.stat.cut1_6.txt
# cut -f 1,2,3,4,5,6 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq2.flank5k.stat.txt >293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq2.flank5k.stat.cut1_6.txt
# cut -f 1,2,3,4,5,6 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq3.flank5k.stat.txt >293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq3.flank5k.stat.cut1_6.txt
# cut -f 1,2,3,4,5,6 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq4.flank5k.stat.txt >293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_CpG_result.c1.CTCF.phq4.flank5k.stat.cut1_6.txt
# cut -f 1,2,3,4,5,6 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq1.flank5k.stat.txt >293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq1.flank5k.stat.cut1_6.txt
# cut -f 1,2,3,4,5,6 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq2.flank5k.stat.txt >293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq2.flank5k.stat.cut1_6.txt
# cut -f 1,2,3,4,5,6 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq3.flank5k.stat.txt >293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq3.flank5k.stat.cut1_6.txt
# cut -f 1,2,3,4,5,6 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq4.flank5k.stat.txt >293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.c1.CTCF.phq4.flank5k.stat.cut1_6.txt
cd /home/sunwenju/projects/20210000_CTCF_CpG_methylation/analysis_based_on_rep_merged_dedup_bam/peak_summit_4_class_flank_stat_and_meta_plot
bash cut1_6_for_CTCF_4_class_summit_meta_plot.sh


# make_all_protein_6bins_reads_CpG_and_GpC_stat.sh ：
# for i in `ls -1 *CTCF*.convert.sam.reads_mC.filter2.txt`; do echo $i; sed '1d' $i >${i%.txt}.nohead.txt; bedtools intersect -a ${i%.txt}.nohead.txt -b CTCF.binding_region.bin1.sorted.bed -f 0.5 -wao >${i%.txt}.intersect.CTCF.bin1.txt; rm ${i%.txt}.nohead.txt; awk 'BEGIN{FS="\t";OFS="\t";print "chrom", "star", "end", "strand", "read_len", "all_mC_count", "CpG_count", "GpC_count", "CpG_mC_count", "GpC_mC_count", "CpG_mC_ratio", "GpC_mC_ratio"}{if($16!=0){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' ${i%.txt}.intersect.CTCF.bin1.txt >${i%.txt}.inCTCF.bin1.txt; rm ${i%.txt}.intersect.CTCF.bin1.txt;done
# for i in `ls -1 *CTCF*.convert.sam.reads_mC.filter2.txt`; do echo $i; sed '1d' $i >${i%.txt}.nohead.txt; bedtools intersect -a ${i%.txt}.nohead.txt -b CTCF.binding_region.bin2.sorted.bed -f 0.5 -wao >${i%.txt}.intersect.CTCF.bin2.txt; rm ${i%.txt}.nohead.txt; awk 'BEGIN{FS="\t";OFS="\t";print "chrom", "star", "end", "strand", "read_len", "all_mC_count", "CpG_count", "GpC_count", "CpG_mC_count", "GpC_mC_count", "CpG_mC_ratio", "GpC_mC_ratio"}{if($16!=0){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' ${i%.txt}.intersect.CTCF.bin2.txt >${i%.txt}.inCTCF.bin2.txt; rm ${i%.txt}.intersect.CTCF.bin2.txt;done
# for i in `ls -1 *CTCF*.convert.sam.reads_mC.filter2.txt`; do echo $i; sed '1d' $i >${i%.txt}.nohead.txt; bedtools intersect -a ${i%.txt}.nohead.txt -b CTCF.binding_region.bin3.sorted.bed -f 0.5 -wao >${i%.txt}.intersect.CTCF.bin3.txt; rm ${i%.txt}.nohead.txt; awk 'BEGIN{FS="\t";OFS="\t";print "chrom", "star", "end", "strand", "read_len", "all_mC_count", "CpG_count", "GpC_count", "CpG_mC_count", "GpC_mC_count", "CpG_mC_ratio", "GpC_mC_ratio"}{if($16!=0){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' ${i%.txt}.intersect.CTCF.bin3.txt >${i%.txt}.inCTCF.bin3.txt; rm ${i%.txt}.intersect.CTCF.bin3.txt;done
# for i in `ls -1 *CTCF*.convert.sam.reads_mC.filter2.txt`; do echo $i; sed '1d' $i >${i%.txt}.nohead.txt; bedtools intersect -a ${i%.txt}.nohead.txt -b CTCF.binding_region.bin4.sorted.bed -f 0.5 -wao >${i%.txt}.intersect.CTCF.bin4.txt; rm ${i%.txt}.nohead.txt; awk 'BEGIN{FS="\t";OFS="\t";print "chrom", "star", "end", "strand", "read_len", "all_mC_count", "CpG_count", "GpC_count", "CpG_mC_count", "GpC_mC_count", "CpG_mC_ratio", "GpC_mC_ratio"}{if($16!=0){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' ${i%.txt}.intersect.CTCF.bin4.txt >${i%.txt}.inCTCF.bin4.txt; rm ${i%.txt}.intersect.CTCF.bin4.txt;done
# for i in `ls -1 *CTCF*.convert.sam.reads_mC.filter2.txt`; do echo $i; sed '1d' $i >${i%.txt}.nohead.txt; bedtools intersect -a ${i%.txt}.nohead.txt -b CTCF.binding_region.bin5.sorted.bed -f 0.5 -wao >${i%.txt}.intersect.CTCF.bin5.txt; rm ${i%.txt}.nohead.txt; awk 'BEGIN{FS="\t";OFS="\t";print "chrom", "star", "end", "strand", "read_len", "all_mC_count", "CpG_count", "GpC_count", "CpG_mC_count", "GpC_mC_count", "CpG_mC_ratio", "GpC_mC_ratio"}{if($16!=0){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' ${i%.txt}.intersect.CTCF.bin5.txt >${i%.txt}.inCTCF.bin5.txt; rm ${i%.txt}.intersect.CTCF.bin5.txt;done
# for i in `ls -1 *CTCF*.convert.sam.reads_mC.filter2.txt`; do echo $i; sed '1d' $i >${i%.txt}.nohead.txt; bedtools intersect -a ${i%.txt}.nohead.txt -b CTCF.binding_region.bin6.sorted.bed -f 0.5 -wao >${i%.txt}.intersect.CTCF.bin6.txt; rm ${i%.txt}.nohead.txt; awk 'BEGIN{FS="\t";OFS="\t";print "chrom", "star", "end", "strand", "read_len", "all_mC_count", "CpG_count", "GpC_count", "CpG_mC_count", "GpC_mC_count", "CpG_mC_ratio", "GpC_mC_ratio"}{if($16!=0){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' ${i%.txt}.intersect.CTCF.bin6.txt >${i%.txt}.inCTCF.bin6.txt; rm ${i%.txt}.intersect.CTCF.bin6.txt;done
nohup bash make_all_protein_6bins_reads_CpG_and_GpC_stat.sh >make_all_protein_6bins_reads_CpG_and_GpC_stat.sh.runlog 2>&1 &
tail -f make_all_protein_6bins_reads_CpG_and_GpC_stat.sh.runlog


# intersection of the interval formed by CTCF summits ± 150
# 1_make_intersect_data_for_sorted_plot.sh：
# for i in `ls -1 *CTCF*.convert.sam.reads_mC.filtered2.forIntersect.sorted.txt`; do echo $i; sed '1d' $i >${i%.txt}.nohead.txt; bedtools intersect -a ${i%.txt}.nohead.txt -b CTCF.summit.all.sorted.extend.b150.bed -f 0.5 -wao >${i%.forIntersect.sorted.txt}.intersect.CTCF.eb150.txt; rm ${i%.txt}.nohead.txt; awk 'BEGIN{FS="\t";OFS="\t";print "chrom", "star", "end", "strand", "gSeq", "read_str", "read_mC_str", "xr_tag", "xg_tag", "CpG_count", "GpC_count", "CpG_mC_count", "GpC_mC_count", "ctcf_bed_chr", "ctcf_bed_start", "ctcf_bed_end", "ctcf_bed_name", "overlap_len"}{if($15!=-1){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}}' ${i%.forIntersect.sorted.txt}.intersect.CTCF.eb150.txt >${i%.forIntersect.sorted.txt}.inCTCF.eb150.txt; rm ${i%.forIntersect.sorted.txt}.intersect.CTCF.eb150.txt;done
nohup bash 1_make_intersect_data_for_sorted_plot.sh >1_make_intersect_data_for_sorted_plot.sh.runlog 2>&1 &
tail -f 1_make_intersect_data_for_sorted_plot.sh.runlog


# 2_make_stat_data_for_sorted_plot.sh：
# for i in `ls -1 *inCTCF.eb150.txt`;do echo $i; python3 make_data_for_TF_region_plot_from_reads_mC.filtered.2.inTF.eb150.py -a CTCF.summit.all.sorted.extend.b150.bed.fa -b $i -o ${i%.txt}.stat_and_plot.txt;done
nohup bash 2_make_stat_data_for_sorted_plot.sh >2_make_stat_data_for_sorted_plot.sh.runlog 2>&1 &
tail -f 2_make_stat_data_for_sorted_plot.sh.runlog

# region_name
# region_sequence
# reads_count
# read_pos_list
# read_seq_list
# read_ov_mc_str_list
# read_strand_list
# cpg_pos_list
# gpc_pos_list,
# cpg_count_list
# gpc_count_list
# cpg_mc_count_list
# gpc_mc_count_list
# cpg_mc_ratio_list
# gpc_mc_ratio_list,
# cor_cpg_gpc_count
# pv_cpg_gpc_count
# cor_cpg_gpc_mc_count
# pv_cpg_gpc_mc_count
# cor_cpg_gpc_mc_ratio
# pv_cpg_gpc_mc_ratio

#
for i in `ls -1 *plot_top*.txt`;do echo $i; python3 plot_tf_region_reads.v2.py -a $i -o ${i%.txt}.pdf;done


###############################################################################
# GpC mC peak calling and stat
###############################################################################
# GpmC_island_identification_new
cd /home/sunwenju/projects/20210000_CTCF_CpG_methylation/analysis_based_on_rep_merged_dedup_bam/GpmC_island_identification_new/
# 
ln -s /home/sunwenju/projects/20210000_CTCF_CpG_methylation/analysis_based_on_rep_merged_dedup_bam/2_all_samples_genome_wide_CpG_GpC_sites_info/293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.txt ./293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.txt
ln -s /home/sunwenju/projects/20210000_CTCF_CpG_methylation/analysis_based_on_rep_merged_dedup_bam/2_all_samples_genome_wide_CpG_GpC_sites_info/293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.txt ./293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.txt

# *bismark.all_GpC_result.txt ：bde6 + 4
# <chrom>   <start> <end>   <new_name>  <m_percentage>  <strand>    <gpc_type>  <coverage>  <mc_site_count> <unmc_site_count>
# chr5  11832   11933   GpC_site_1  0.0000  -   GCA 2   0   2
awk 'BEGIN{FS="\t";OFS="\t"}{if($5>0){print $0}}' 293T_ERT_CTCF_24h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.txt | sort -k1,1 -k2,2n >293T_ERT_CTCF_24h_pd_rep12.bismark.dedup.merged.all_GpmC.sorted.txt
awk 'BEGIN{FS="\t";OFS="\t"}{if($5>0){print $0}}' 293T_ERT_CTCF_48h_pd_rep12.bismark_bt2_pe.dedup.merged.all_GpC_result.txt | sort -k1,1 -k2,2n >293T_ERT_CTCF_48h_pd_rep12.bismark.dedup.merged.all_GpmC.sorted.txt

# 
bedtools slop -i 293T_ERT_CTCF_24h_pd_rep12.bismark.dedup.merged.all_GpmC.sorted.txt -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes -b 50 |sort -k1,1 -k2,2n >293T_ERT_CTCF_24h_pd_rep12.bismark.dedup.merged.all_GpmC.extB50.sorted.txt
bedtools slop -i 293T_ERT_CTCF_48h_pd_rep12.bismark.dedup.merged.all_GpmC.sorted.txt -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes -b 50 |sort -k1,1 -k2,2n >293T_ERT_CTCF_48h_pd_rep12.bismark.dedup.merged.all_GpmC.extB50.sorted.txt


# *bismark.all_GpC_result.txt ：bde6 + 4
# <chrom>   <start> <end>   <new_name>  <m_percentage>  <strand>    <gpc_type>  <coverage>  <mc_site_count> <unmc_site_count>
# chr5  11832   11933   GpC_site_1  0.0000  -   GCA 2   0   2
bedtools merge -c 4,6,7,8,9,10,5,4 -o collapse,collapse,collapse,mean,mean,mean,mean,count -i 293T_ERT_CTCF_24h_pd_rep12.bismark.dedup.merged.all_GpmC.extB50.sorted.txt |sort -k1,1 -k2,2n | awk 'BEGIN{FS=" ";OFS="\t"}{print $1,$2,$3,"cp_"NR,$3-$2,".",$4,$5,$6,$7,$8,$9,$10,$11}' >293T_ERT_CTCF_24h_pd_rep12.bismark.dedup.merged.all_GpmC.extB50.merged.sorted.txt
bedtools merge -c 4,6,7,8,9,10,5,4 -o collapse,collapse,collapse,mean,mean,mean,mean,count -i 293T_ERT_CTCF_48h_pd_rep12.bismark.dedup.merged.all_GpmC.extB50.sorted.txt |sort -k1,1 -k2,2n | awk 'BEGIN{FS=" ";OFS="\t"}{print $1,$2,$3,"cp_"NR,$3-$2,".",$4,$5,$6,$7,$8,$9,$10,$11}' >293T_ERT_CTCF_48h_pd_rep12.bismark.dedup.merged.all_GpmC.extB50.merged.sorted.txt

# <chrom>   <peak_start>    <peak_end>  <candidate_peak_name>   <peak_len>  <peak_strand.> <site_name_list> <strand_list>   <type_list> <mean_cov>  <mean_mc_rc>    <mean_unmc_rc>  <mean_mc_ratio> <mc_count>
#
cut -f 1,2,3,4 293T_ERT_CTCF_24h_pd_rep12.bismark.dedup.merged.all_GpmC.extB50.merged.sorted.txt >293T_ERT_CTCF_24h_pd_rep12.bismark.dedup.merged.all_GpmC.extB50.merged.sorted.cut1_4.bed
cut -f 1,2,3,4 293T_ERT_CTCF_48h_pd_rep12.bismark.dedup.merged.all_GpmC.extB50.merged.sorted.txt >293T_ERT_CTCF_48h_pd_rep12.bismark.dedup.merged.all_GpmC.extB50.merged.sorted.cut1_4.bed

# GpC peak calling
nohup time python3 make_GpC_peaks_stat_info_new2.py -g /home/sunwenju/3_genomeData/human/hg19/hg19.fa -r 293T_ERT_CTCF_24h_pd_rep12.bismark.dedup.merged.all_GpmC.sorted.txt -m 293T_ERT_CTCF_24h_pd_rep12.bismark.dedup.merged.all_GpmC.extB50.merged.sorted.txt -b 293T_ERT_CTCF_24h_pd_rep12.bismark.dedup.merged.all_GpmC.extB50.merged.sorted.cut1_4.bed -o 293T_ERT_CTCF_24h_pd_rep12.bismark.dedup.merged.all_GpmC.extB50.merged.sorted.peaks_with_stat_result.txt >293T_ERT_CTCF_24h_pd_rep12.bismark.dedup.merged.all_GpmC.extB50.merged.sorted.peaks_stat.runlog 2>&1 &
nohup time python3 make_GpC_peaks_stat_info_new2.py -g /home/sunwenju/3_genomeData/human/hg19/hg19.fa -r 293T_ERT_CTCF_48h_pd_rep12.bismark.dedup.merged.all_GpmC.sorted.txt -m 293T_ERT_CTCF_48h_pd_rep12.bismark.dedup.merged.all_GpmC.extB50.merged.sorted.txt -b 293T_ERT_CTCF_48h_pd_rep12.bismark.dedup.merged.all_GpmC.extB50.merged.sorted.cut1_4.bed -o 293T_ERT_CTCF_48h_pd_rep12.bismark.dedup.merged.all_GpmC.extB50.merged.sorted.peaks_with_stat_result.txt >293T_ERT_CTCF_48h_pd_rep12.bismark.dedup.merged.all_GpmC.extB50.merged.sorted.peaks_stat.runlog 2>&1 &
# tail -f 293T_ERT_CTCF_24h_pd_rep12.bismark.dedup.merged.all_GpmC.extB50.merged.sorted.peaks_stat.runlog
# tail -f 293T_ERT_CTCF_48h_pd_rep12.bismark.dedup.merged.all_GpmC.extB50.merged.sorted.peaks_stat.runlog
# wc -l *peaks_with_stat_result.withFDR.txt

# 不同过滤条件下 GpC peak 数量统计
nohup bash 293T_ERT_CTCF_24h_pd_rep12.peaks_with_stat.withFDR.diff_filter.stats.sh >293T_ERT_CTCF_24h_pd_rep12.peaks_with_stat.withFDR.diff_filter.stats.sh.runlog 2>&1 &
nohup bash 293T_ERT_CTCF_48h_pd_rep12.peaks_with_stat.withFDR.diff_filter.stats.sh >293T_ERT_CTCF_48h_pd_rep12.peaks_with_stat.withFDR.diff_filter.stats.sh.runlog 2>&1 &
# tail -f 293T_ERT_CTCF_24h_pd_rep12.peaks_with_stat.withFDR.diff_filter.stats.sh.runlog
# tail -f 293T_ERT_CTCF_48h_pd_rep12.peaks_with_stat.withFDR.diff_filter.stats.sh.runlog

# make 24h 48h GpC peaks bed
sed '1d' 293T_ERT_CTCF_24h_pd_rep12.bismark.dedup.merged.all_GpmC.extB50.merged.sorted.peaks_with_stat_result.withFDR.txt | sort -k1,1 -k2,2n | awk 'BEGIN{FS="\t";OFS="\t";n=0}{if($54<0.05){n++; print $1,$2,$3,"gpc_"n,$5,$6,$54,$10,$13}}' >293T_ERT_CTCF_24h_pd_rep12.all_GpmC.peaks.fdr0.05.filtered.sorted.bedPlus.txt
sed '1d' 293T_ERT_CTCF_48h_pd_rep12.bismark.dedup.merged.all_GpmC.extB50.merged.sorted.peaks_with_stat_result.withFDR.txt | sort -k1,1 -k2,2n | awk 'BEGIN{FS="\t";OFS="\t";n=0}{if($54<0.05){n++; print $1,$2,$3,"gpc_"n,$5,$6,$54,$10,$13}}' >293T_ERT_CTCF_48h_pd_rep12.all_GpmC.peaks.fdr0.05.filtered.sorted.bedPlus.txt
# data structure
# "chrom"
# "peak_start"
# "peak_end"
# "np_name"
# "peak_len"
# "peak_strand"
# "lmd_max_fdr"
# "mean_cov"
# "mean_mc_ratio"
wc -l 293T_ERT_CTCF_24h_pd_rep12.all_GpmC.peaks.fdr0.05.filtered.sorted.bedPlus.txt
# 92518
wc -l 293T_ERT_CTCF_48h_pd_rep12.all_GpmC.peaks.fdr0.05.filtered.sorted.bedPlus.txt
# 39771

cut -f 1-6 293T_ERT_CTCF_24h_pd_rep12.all_GpmC.peaks.fdr0.05.filtered.sorted.bedPlus.txt >293T_ERT_CTCF_24h_pd_rep12.all_GpmC.peaks.fdr0.05.filtered.sorted.bed
cut -f 1-6 293T_ERT_CTCF_48h_pd_rep12.all_GpmC.peaks.fdr0.05.filtered.sorted.bedPlus.txt >293T_ERT_CTCF_48h_pd_rep12.all_GpmC.peaks.fdr0.05.filtered.sorted.bed
wc -l 293T_ERT_CTCF_24h_pd_rep12.all_GpmC.peaks.fdr0.05.filtered.sorted.bed
# 92518
wc -l 293T_ERT_CTCF_48h_pd_rep12.all_GpmC.peaks.fdr0.05.filtered.sorted.bed
# 39771

