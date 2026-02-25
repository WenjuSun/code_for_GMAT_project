#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re
import subprocess
import logging
from scipy import stats


def create_custom_logger(name):
    formatter = logging.Formatter(fmt='[%(asctime)s][%(levelname)s][%(module)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    return logger


logger = create_custom_logger('root')


def get_rev_com_seq(upper_seq_str):
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
    rev_seq_list = list(reversed(upper_seq_str))
    rev_com_seq_list = [complement_dict[k] for k in rev_seq_list]
    rev_com_seq = ''.join(rev_com_seq_list)
    return rev_com_seq


def get_read_ov_mc_info(ge_seq_str, read_str, read_mc, xr_tag, xg_tag, strand):
    # input : ge_seq_str, read_str, read_mc, xr_tag, xg_tag, strand
    # output : read_len, cpg_mc_pos_list, gpc_mc_pos_list, all_mC_count, CpG_count, GpC_count, CpG_mC_count, GpC_mC_count, CpG_mC_ratio, GpC_mC_ratio
    ge_seq_base_list = list(ge_seq_str)
    read_base_list = list(read_str)
    read_mc_list = list(read_mc)
    read_convert = xr_tag.split(':')[2]
    genome_convert = xg_tag.split(':')[2]
    read_len = len(read_str)
    all_mC_count = 0
    for c in read_mc_list:
        if c.isupper():
            all_mC_count += 1
    # CpG = WCG = ACG + TCG = [AT]CG
    # GpC = GCH = GCA + GCC + GCT = GC[ACT]
    cpg_match = re.finditer('[AT]CG', ge_seq_str)
    cpg_match_groups = []
    cpg_c_pos_list = []
    for cpg_m in cpg_match:
        cpg_match_groups.append(cpg_m.group())
        cpg_c_pos_list.append(cpg_m.start())
    cpg_count = len(cpg_match_groups)

    gpc_match = re.finditer('GC[ACT]', ge_seq_str)
    gpc_match_groups = []
    gpc_c_post_list = []
    for gpc_m in gpc_match:
        gpc_match_groups.append(gpc_m.group())
        gpc_c_post_list.append(gpc_m.start())
    gpc_count = len(gpc_match_groups)

    cpg_mc_count = 0
    gpc_mc_count = 0
    cpg_mc_pos_list = []
    gpc_mc_pos_list = []
    for i in range(len(read_base_list)):
        rb = read_base_list[i]
        gb = ge_seq_base_list[i+1]
        if read_convert == 'CT' and genome_convert == 'CT' and strand == '+' and rb == 'C' and gb == 'C':
            if not read_mc_list[i].isupper():
                print(i, read_mc)
            if i in cpg_c_pos_list:
                cpg_mc_count += 1
                cpg_mc_pos_list.append(str(i+1))
            if i in gpc_c_post_list:
                gpc_mc_count += 1
                gpc_mc_pos_list.append(str(i+1))
        if read_convert == 'CT' and genome_convert == 'GA' and strand == '-' and rb == 'C' and gb == 'C':
            if not read_mc_list[i].isupper():
                print(i, read_mc)
            if i in cpg_c_pos_list:
                cpg_mc_count += 1
                cpg_mc_pos_list.append(str(i+1))
            if i in gpc_c_post_list:
                gpc_mc_count += 1
                gpc_mc_pos_list.append(str(i+1))
        if read_convert == 'GA' and genome_convert == 'CT' and strand == '-' and rb == 'G' and gb == 'G':
            if not read_mc_list[i].isupper():
                print(i, read_mc)
            if i in cpg_c_pos_list:
                cpg_mc_count += 1
                cpg_mc_pos_list.append(str(i+1))
            if i in gpc_c_post_list:
                gpc_mc_count += 1
                gpc_mc_pos_list.append(str(i+1))
        if read_convert == 'GA' and genome_convert == 'GA' and strand == '+' and rb == 'G' and gb == 'G':
            if not read_mc_list[i].isupper():
                print(i, read_mc)
            if i in cpg_c_pos_list:
                cpg_mc_count += 1
                cpg_mc_pos_list.append(str(i+1))
            if i in gpc_c_post_list:
                gpc_mc_count += 1
                gpc_mc_pos_list.append(str(i+1))
    if cpg_mc_count != len(cpg_mc_pos_list):
        print('cpg_mc_count is not equal to the length of cpg_mc_pos_list')
    if gpc_mc_count != len(gpc_mc_pos_list):
        print('gpc_mc_count is not equal to the length of gpc_mc_pos_list')
    cpg_mc_pos_info = ','.join(cpg_mc_pos_list)
    gpc_mc_pos_info = ','.join(gpc_mc_pos_list)
    cpg_mc_ratio = 0.00
    if cpg_count != 0:
        cpg_mc_ratio = format(cpg_mc_count / cpg_count * 100, ".2f")
    gpc_mc_ratio = 0.00
    if gpc_count != 0:
        gpc_mc_ratio = format(gpc_mc_count / gpc_count * 100, '.2f')
    read_mc_stat = "\t".join([str(read_len), cpg_mc_pos_info, gpc_mc_pos_info, str(all_mC_count), str(cpg_count), str(gpc_count), str(cpg_mc_count), str(gpc_mc_count), str(cpg_mc_ratio), str(gpc_mc_ratio)])
    return read_mc_stat


def my_parse_args():
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-a', '--sequence', type=str, required=True, help='ctcf eb150 tab fasta sequence')
    my_parser.add_argument('-b', '--reads', type=str, required=True, help='reads mC filtered 2 inCTCF eb150 result file')
    my_parser.add_argument('-o', '--output', type=str, required=True, help='output file name')
    return my_parser.parse_args()


def main():
    my_args = my_parse_args()
    seq_file = my_args.sequence
    reads_file = my_args.reads
    output_file = my_args.output

    logger.info('command line: python3 make_data_for_ctcf_region_plot_from_reads_mC.filtered.2.inCTCF.eb150.py -a {0} -b {1} -o {2}'.format(seq_file, reads_file, output_file))

    ctcf_seq_dcit = {}
    with open(seq_file, 'r') as seqin:
        for line in seqin:
            name, seq = line.rstrip().split("\t")
            ctcf_seq_dcit[name] = seq.upper()

    out_put_dict = {}
    with open(reads_file, 'r') as readin:
        title = readin.readline()
        for line in readin:
            items = line.rstrip().split("\t")
            chrom = items[0]
            start = int(items[1])
            end = int(items[2])
            strand = items[3]
            gSeq = items[4]
            read_str = items[5]
            read_mC_str = items[6]
            xr_tag = items[7]
            xg_tag = items[8]
            CpG_count = items[9]
            GpC_count = items[10]
            CpG_mC_count = items[11]
            GpC_mC_count = items[12]
            ctcf_bed_chr = items[13]
            ctcf_bed_start = int(items[14])
            ctcf_bed_end = int(items[15])
            ctcf_bed_name = items[16]
            overlap_len = int(items[17])

            if strand == '-':
                gSeq = get_rev_com_seq(items[4])
                read_str = get_rev_com_seq(items[5])
                read_mC_str = ''.join(list(reversed(items[6])))

            read_pos = 0
            gseq_ov = ''
            read_ov_seq = ''
            read_ov_mc_str = ''
            if start < ctcf_bed_start:
                read_pos = 0
                gseq_ov = gSeq[-int(overlap_len)-2:]
                read_ov_seq = read_str[-int(overlap_len):]
                read_ov_mc_str = read_mC_str[-int(overlap_len):]
            elif end > ctcf_bed_end:
                read_pos = 301 - int(overlap_len)
                gseq_ov = gSeq[:int(overlap_len)+2]
                read_ov_seq = read_str[:int(overlap_len)]
                read_ov_mc_str = read_mC_str[:int(overlap_len)]
            else:
                read_pos = int(start) - int(ctcf_bed_start)
                gseq_ov = gSeq
                read_ov_seq = read_str
                read_ov_mc_str = read_mC_str

            if strand == '-':
                gseq_ov = get_rev_com_seq(gseq_ov)
                read_ov_seq = get_rev_com_seq(read_ov_seq)
                read_ov_mc_str = ''.join(list(reversed(read_ov_mc_str)))

            # if len(gseq_ov) != len(read_ov_seq)+2 or len(gseq_ov) != len(read_ov_mc_str)+2:
            #     print(gseq_ov, read_ov_seq, read_ov_mc_str, xr_tag, xg_tag, strand, sep="\n")
            # print(gseq_ov, ' '+read_ov_seq, ' '+read_ov_mc_str, xr_tag, xg_tag, strand, sep="\n")
            read_ov_mc_info = get_read_ov_mc_info(gseq_ov, read_ov_seq, read_ov_mc_str, xr_tag, xg_tag, strand)
            # print(gseq_ov, ' '+read_ov_seq, ' '+read_ov_mc_str, xr_tag, xg_tag, strand, read_ov_mc_info, sep="\n")
            read_ov_mc_info_list = read_ov_mc_info.split("\t")
            all_mC_count = read_ov_mc_info_list[3]
            cpg_count = read_ov_mc_info_list[4]
            gpc_count = read_ov_mc_info_list[5]
            cpg_mc_count = read_ov_mc_info_list[6]
            gpc_mc_count = read_ov_mc_info_list[7]
            cpg_mc_ratio = read_ov_mc_info_list[8]
            gpc_mc_ratio = read_ov_mc_info_list[9]

            cpg_mc_p_list = []
            gpc_mc_p_list = []
            a_list = read_ov_mc_info_list[1].split(',')
            b_list = read_ov_mc_info_list[2].split(',')
            for x in a_list:
                if x != '':
                    cpg_mc_p_list.append(x)
            for y in b_list:
                if y != '':
                    gpc_mc_p_list.append(y)

            if strand == '-':
                gseq_ov = get_rev_com_seq(gseq_ov)
                read_ov_seq = get_rev_com_seq(read_ov_seq)
                read_ov_mc_str = ''.join(list(reversed(read_ov_mc_str)))
                # print(cpg_mc_p_list, gpc_mc_p_list)
                cpg_mc_p_list = []
                gpc_mc_p_list = []
                for x in a_list:
                    if x != '':
                        cpg_mc_p_list.append(str(overlap_len-int(x)+1))
                for y in b_list:
                    if y != '':
                        gpc_mc_p_list.append(str(overlap_len-int(y)+1))
            # print(ctcf_bed_name, gseq_ov, read_ov_seq, read_ov_mc_str, strand, cpg_mc_p_list, gpc_mc_p_list,
            #       all_mC_count, cpg_count, gpc_count, cpg_mc_count, gpc_mc_count, cpg_mc_ratio, gpc_mc_ratio, sep="\n")

            if ctcf_bed_name in ctcf_seq_dcit:
                if ctcf_bed_name not in out_put_dict:
                    out_put_dict[ctcf_bed_name] = {}
                    out_put_dict[ctcf_bed_name]['region_seq'] = ctcf_seq_dcit[ctcf_bed_name]
                    out_put_dict[ctcf_bed_name]['read_pos_list'] = []
                    out_put_dict[ctcf_bed_name]['read_pos_list'].append(str(read_pos))
                    out_put_dict[ctcf_bed_name]['read_seq_list'] = []
                    out_put_dict[ctcf_bed_name]['read_seq_list'].append(read_ov_seq)
                    out_put_dict[ctcf_bed_name]['read_ov_mc_str'] = []
                    out_put_dict[ctcf_bed_name]['read_ov_mc_str'].append(read_ov_mc_str)
                    out_put_dict[ctcf_bed_name]['read_strand_list'] = []
                    out_put_dict[ctcf_bed_name]['read_strand_list'].append(strand)
                    out_put_dict[ctcf_bed_name]['cpg_mc_pos_list'] = []
                    out_put_dict[ctcf_bed_name]['cpg_mc_pos_list'].append(','.join(cpg_mc_p_list))
                    out_put_dict[ctcf_bed_name]['gpc_mc_pos_list'] = []
                    out_put_dict[ctcf_bed_name]['gpc_mc_pos_list'].append(','.join(gpc_mc_p_list))
                    out_put_dict[ctcf_bed_name]['all_mC_count_list'] = []
                    out_put_dict[ctcf_bed_name]['all_mC_count_list'].append(all_mC_count)
                    out_put_dict[ctcf_bed_name]['cpg_count_list'] = []
                    out_put_dict[ctcf_bed_name]['cpg_count_list'].append(cpg_count)
                    out_put_dict[ctcf_bed_name]['gpc_count_list'] = []
                    out_put_dict[ctcf_bed_name]['gpc_count_list'].append(gpc_count)
                    out_put_dict[ctcf_bed_name]['cpg_mc_count_list'] = []
                    out_put_dict[ctcf_bed_name]['cpg_mc_count_list'].append(cpg_mc_count)
                    out_put_dict[ctcf_bed_name]['gpc_mc_count_list'] = []
                    out_put_dict[ctcf_bed_name]['gpc_mc_count_list'].append(gpc_mc_count)
                    out_put_dict[ctcf_bed_name]['cpg_mc_ratio_list'] = []
                    out_put_dict[ctcf_bed_name]['cpg_mc_ratio_list'].append(cpg_mc_ratio)
                    out_put_dict[ctcf_bed_name]['gpc_mc_ratio_list'] = []
                    out_put_dict[ctcf_bed_name]['gpc_mc_ratio_list'].append(gpc_mc_ratio)
                else:
                    out_put_dict[ctcf_bed_name]['read_pos_list'].append(str(read_pos))
                    out_put_dict[ctcf_bed_name]['read_seq_list'].append(read_ov_seq)
                    out_put_dict[ctcf_bed_name]['read_ov_mc_str'].append(read_ov_mc_str)
                    out_put_dict[ctcf_bed_name]['read_strand_list'].append(strand)
                    out_put_dict[ctcf_bed_name]['cpg_mc_pos_list'].append(','.join(cpg_mc_p_list))
                    out_put_dict[ctcf_bed_name]['gpc_mc_pos_list'].append(','.join(gpc_mc_p_list))
                    out_put_dict[ctcf_bed_name]['cpg_count_list'].append(cpg_count)
                    out_put_dict[ctcf_bed_name]['gpc_count_list'].append(gpc_count)
                    out_put_dict[ctcf_bed_name]['cpg_mc_count_list'].append(cpg_mc_count)
                    out_put_dict[ctcf_bed_name]['gpc_mc_count_list'].append(gpc_mc_count)
                    out_put_dict[ctcf_bed_name]['cpg_mc_ratio_list'].append(cpg_mc_ratio)
                    out_put_dict[ctcf_bed_name]['gpc_mc_ratio_list'].append(gpc_mc_ratio)

    for region_name in out_put_dict:
        region_sequence = out_put_dict[region_name]['region_seq']
        reads_count = len(out_put_dict[region_name]['read_pos_list'])
        read_pos_list = ';'.join(out_put_dict[region_name]['read_pos_list'])
        read_seq_list = ';'.join(out_put_dict[region_name]['read_seq_list'])
        read_ov_mc_str_list = ';'.join(out_put_dict[region_name]['read_ov_mc_str'])
        read_strand_list = ';'.join(out_put_dict[region_name]['read_strand_list'])
        cpg_pos_list = ';'.join(out_put_dict[region_name]['cpg_mc_pos_list'])
        gpc_pos_list = ';'.join(out_put_dict[region_name]['gpc_mc_pos_list'])
        cpg_count_list = ';'.join(out_put_dict[region_name]['cpg_count_list'])
        gpc_count_list = ';'.join(out_put_dict[region_name]['gpc_count_list'])
        cpg_mc_count_list = ';'.join(out_put_dict[region_name]['cpg_mc_count_list'])
        gpc_mc_count_list = ';'.join(out_put_dict[region_name]['gpc_mc_count_list'])
        cpg_mc_ratio_list = ';'.join(out_put_dict[region_name]['cpg_mc_ratio_list'])
        gpc_mc_ratio_list = ';'.join(out_put_dict[region_name]['gpc_mc_ratio_list'])

        # spearman correlation stat for CpG_count and GpC_count
        cor_cpg_gpc_count, pv_cpg_gpc_count = stats.spearmanr(out_put_dict[region_name]['cpg_count_list'], out_put_dict[region_name]['gpc_count_list'])
        # spearman correlation stat for CpG_mC_count and GpC_mC_count
        cor_cpg_gpc_mc_count, pv_cpg_gpc_mc_count = stats.spearmanr(out_put_dict[region_name]['cpg_mc_count_list'], out_put_dict[region_name]['gpc_mc_count_list'])
        # spearman correlation stat for CpG_mC_ratio and GpC_mC_ratio
        cor_cpg_gpc_mc_ratio, pv_cpg_gpc_mc_ratio = stats.spearmanr(out_put_dict[region_name]['cpg_mc_ratio_list'], out_put_dict[region_name]['gpc_mc_ratio_list'])

        print(region_name, region_sequence, reads_count, read_pos_list, read_seq_list, read_ov_mc_str_list, read_strand_list, cpg_pos_list, gpc_pos_list,
              cpg_count_list, gpc_count_list,cpg_mc_count_list, gpc_mc_count_list, cpg_mc_ratio_list, gpc_mc_ratio_list,
              cor_cpg_gpc_count, pv_cpg_gpc_count, cor_cpg_gpc_mc_count, pv_cpg_gpc_mc_count, cor_cpg_gpc_mc_ratio, pv_cpg_gpc_mc_ratio, sep="\t", end="\n")

    logger.info('All things have done!')


if __name__ == '__main__':
    main()
