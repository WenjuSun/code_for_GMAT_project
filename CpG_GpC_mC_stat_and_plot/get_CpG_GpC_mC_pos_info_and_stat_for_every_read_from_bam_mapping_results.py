#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re
import subprocess
import logging


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


def make_temp_bed(input_file, temp_bed, chrom_szie_dcit):
    logger.info('making temp bed file...')
    with open(input_file, 'r') as sam_in, open(temp_bed, 'w') as temp:
        for line in sam_in:
            items = line.rstrip().split("\t")
            # <QNAME>	<FLAG>	<RNAME>	<POS>	<MAPQ>	<CIGAR>	<RNEXT>	<PNEXT>	<TLEN>	<SEQ>	<QUAL>	<NM-TAG>	<MD-TAG>	<XM-TAG>	<XR-TAG>	<XG-TAG>
            read_name = items[0]
            flag = items[1]
            chrom = items[2]
            g_start = int(items[3]) - 1
            mapq = items[4]
            cigar = items[5]

            strand = '+'
            if int(flag) & 16 == 16:
                strand = '-'

            read_n = ''
            if int(flag) & 64 == 64:
                read_n = '1'
            elif int(flag) & 128 == 128:
                read_n = '2'
            new_name = read_name + '/'+read_n
            cigar_c = [x for x in re.split("\d+", cigar) if x]
            cigar_n = [x for x in re.split("[A-Z]", cigar) if x]
            # print(cigar_c, cigar_n)
            cigar_dict = {}
            for i in range(len(cigar_c)):
                if cigar_c[i] not in cigar_dict:
                    cigar_dict[cigar_c[i]] = int(cigar_n[i])
                else:
                    cigar_dict[cigar_c[i]] += int(cigar_n[i])

            map_g_len = cigar_dict['M']
            if 'D' in cigar_dict:
                map_g_len = cigar_dict['M'] + cigar_dict['D']
            g_end = g_start + map_g_len
            # bam to bed out put
            # print(chrom, g_start, g_end, new_name, mapq, strand, sep="\t", end="\n")
            if len(cigar_dict.keys()) == 1 and 'M' in cigar_dict.keys():
                start = g_start - 1
                end = g_end + 1
                if start < 0:
                    start = 0
                if end > int(chrom_szie_dcit[chrom]):
                    end = int(chrom_szie_dcit[chrom])
                print(chrom, start, end, new_name, mapq, strand, sep="\t", end="\n", file=temp)


def get_fasta(bed6_file, temp_fasta):
    logger.info('getting fasta for temp bed file...')
    cmd = 'bedtools getfasta -fi /home/sunwenju/3_genomeData/human/hg19/hg19.fa -bed ' + bed6_file + ' -fo ' + temp_fasta + ' -name -tab -s'
    # logger.info(cmd)
    os.system(cmd)


def make_expanded_mdz_str(md_tag, strand):
    md_n = [x for x in re.split("[A-Z]+", md_tag) if x]
    md_c = [x for x in re.split("\d+", md_tag) if x]
    md_str = ''
    # print(md_n, md_c)
    for i in range(len(md_c)):
        md_str += ('.'*int(md_n[i])+md_c[i])
    md_str += '.'*int(md_n[-1])
    if strand == '-':
        md_str = ''.join(list(reversed(md_str)))
    # print(md_str)
    return md_str


def get_read_mc_stat(ge_seq_str, read_str, read_mc, xr_tag, xg_tag, strand):
    # input : ge_seq_str, read_str, read_mc, xr_tag, xg_tag, strand
    # output : read_len, cpg_mc_pos_info, gpc_mc_pos_info, all_mC_count, CpG_count, GpC_count, CpG_mC_count, GpC_mC_count, CpG_mC_ratio, GpC_mC_ratio
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
                cpg_mc_pos_list.append(str(i))
            if i in gpc_c_post_list:
                gpc_mc_count += 1
                gpc_mc_pos_list.append(str(i))
        if read_convert == 'CT' and genome_convert == 'GA' and strand == '-' and rb == 'C' and gb == 'C':
            if not read_mc_list[i].isupper():
                print(i, read_mc)
            if i in cpg_c_pos_list:
                cpg_mc_count += 1
                cpg_mc_pos_list.append(str(i))
            if i in gpc_c_post_list:
                gpc_mc_count += 1
                gpc_mc_pos_list.append(str(i))
        if read_convert == 'GA' and genome_convert == 'CT' and strand == '-' and rb == 'G' and gb == 'G':
            if not read_mc_list[i].isupper():
                print(i, read_mc)
            if i in cpg_c_pos_list:
                cpg_mc_count += 1
                cpg_mc_pos_list.append(str(i))
            if i in gpc_c_post_list:
                gpc_mc_count += 1
                gpc_mc_pos_list.append(str(i))
        if read_convert == 'GA' and genome_convert == 'GA' and strand == '+' and rb == 'G' and gb == 'G':
            if not read_mc_list[i].isupper():
                print(i, read_mc)
            if i in cpg_c_pos_list:
                cpg_mc_count += 1
                cpg_mc_pos_list.append(str(i))
            if i in gpc_c_post_list:
                gpc_mc_count += 1
                gpc_mc_pos_list.append(str(i))
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
    my_parser.add_argument('-i', '--input', type=str, required=True, help='input sam format mapping result')
    my_parser.add_argument('-o', '--output', type=str, required=True, help='output file name')
    return my_parser.parse_args()


def main():
    my_args = my_parse_args()
    input_file = my_args.input
    output_file = my_args.output

    logger.info('command line: python3 get_CpG_GpC_mC_pos_info_and_stat_for_every_read_from_bam_mapping_results.py -i {0} -o {1}'.format(input_file, output_file))
    cs_file = '/home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes'
    cs_dict = {}
    with open(cs_file, 'r') as  csin:
        for line in csin:
            chrom, size = line.rstrip().split("\t")
            cs_dict[chrom] = size
    temp_bed = input_file+".temp.bed"
    temp_fasta = input_file+".temp.fa"
    make_temp_bed(input_file, temp_bed, cs_dict)
    get_fasta(temp_bed, temp_fasta)
    os.system("rm "+temp_bed)
    seq_dict = {}
    with open(temp_fasta, 'r') as fin:
        for line in fin:
            seq_name, seq_str = line.rstrip().split("\t")
            seq_str = seq_str.upper()
            seq_dict[seq_name] = seq_str
    os.system("rm "+temp_fasta)

    logger.info('getting mC stat for every read...')
    with open(input_file, 'r') as sam_in, open(output_file, 'w') as fout:
        for line in sam_in:
            items = line.rstrip().split("\t")
            # <QNAME>   <FLAG>  <RNAME> <POS>   <MAPQ>  <CIGAR> <RNEXT> <PNEXT> <TLEN>  <SEQ>   <QUAL>  <NM-TAG>    <MD-TAG>    <XM-TAG>    <XR-TAG>    <XG-TAG>
            read_name = items[0]
            flag = items[1]
            chrom = items[2]
            g_start = int(items[3]) - 1
            mapq = items[4]
            cigar = items[5]
            sequence = items[9]
            md_tag = items[12].split(':')[2]
            xm_tag = items[13].split(':')[2]
            xr_tag = items[14]
            xg_tag = items[15]

            strand = '+'
            if int(flag) & 16 == 16:
                strand = '-'

            read_n = ''
            if int(flag) & 64 == 64:
                read_n = '1'
            elif int(flag) & 128 == 128:
                read_n = '2'
            new_name = read_name + '/' + read_n
            cigar_c = [x for x in re.split("\d+", cigar) if x]
            cigar_n = [x for x in re.split("[A-Z]", cigar) if x]
            # print(cigar_c, cigar_n)
            cigar_dict = {}
            for i in range(len(cigar_c)):
                if cigar_c[i] not in cigar_dict:
                    cigar_dict[cigar_c[i]] = int(cigar_n[i])
                else:
                    cigar_dict[cigar_c[i]] += int(cigar_n[i])

            map_g_len = cigar_dict['M']
            if 'D' in cigar_dict:
                map_g_len = cigar_dict['M'] + cigar_dict['D']
            g_end = g_start + map_g_len
            # bam to bed out put
            # print(chrom, g_start, g_end, new_name, mapq, strand, sep="\t", end="\n")
            read_str = sequence
            read_mc = xm_tag
            if strand == '-':
                read_str = get_rev_com_seq(sequence)
                read_mc = ''.join(list(reversed(xm_tag)))

            # only perfect mapping, without insert or deletion
            if len(cigar_dict.keys()) == 1 and 'M' in cigar_dict.keys() and g_start >= 1 and g_end <= int(cs_dict[chrom]):
                expanded_mdz = make_expanded_mdz_str(md_tag, strand)
                read_mC_stat_str = get_read_mc_stat(seq_dict[new_name], read_str, read_mc, xr_tag, xg_tag, strand)
                print(chrom, g_start, g_end, new_name, mapq, strand, seq_dict[new_name], read_str, read_mc, expanded_mdz, xr_tag, xg_tag, md_tag, read_mC_stat_str, sep="\t", end="\n", file=fout)
                # if xr_tag == 'XR:Z:CT' and xg_tag == 'XG:Z:CT':
                # if xr_tag == 'XR:Z:CT' and xg_tag == 'XG:Z:GA':
                # if xr_tag == 'XR:Z:GA' and xg_tag == 'XG:Z:CT':
                # if xr_tag == 'XR:Z:GA' and xg_tag == 'XG:Z:GA':
                #     print(chrom, g_start, g_end, new_name, mapq, strand, sep="\t", end="\n")
                #     print(seq_dict[new_name], end="\n")
                #     print(' ' + read_str + ' ', end="\n")
                #     print(' ' + read_mc + ' ', end="\n")
                #     print(' ' + expanded_mdz + ' ', end="\n")
                #     print(' ' + read_mC_stat_str + ' ', end="\n")
                #     print(xr_tag, xg_tag, md_tag, sep="\t", end="\n\n")

    logger.info('All things have done!')


if __name__ == '__main__':
    main()
