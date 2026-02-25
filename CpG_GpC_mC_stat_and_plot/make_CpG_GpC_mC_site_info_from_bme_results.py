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


def bme_to_filterd_bed(bme_result, filterd_bed):
    logger.info('do bme_to_filterd_bed for {0} ...'.format(bme_result))
    # 
    convert_cmd = '''awk 'BEGIN{FS="\\t";OFS="\\t";n=0}{if($4!=0 || $5!=0){n++;print $1, $2, $2+1, "mC_site_"n, $4, $3, $5}}' ''' + bme_result + ' >' + filterd_bed
    logger.info(convert_cmd)
    os.system(convert_cmd)


def get_fasta_for_c_site(all_mc_site_info, all_mc_site_info_fa):
    logger.info('get fasta sequence for {0} ...'.format(all_mc_site_info))
    get_cmd = 'bedtools getfasta -fi /home/sunwenju/3_genomeData/human/hg19/hg19.fa -bed '+all_mc_site_info+' -fo '+all_mc_site_info_fa+' -name -tab -s'
    logger.info(get_cmd)
    os.system(get_cmd)


def make_c_site_extend_bed(all_mC_site_info, all_mC_site_info_extend):
    logger.info('make extend bed for {0} ...'.format(all_mC_site_info))
    make_cmd = '''awk 'BEGIN{FS="\\t";OFS="\\t";n=0}{print $1, $2-1, $3+1, $4, $5, $6, $7}' ''' + all_mC_site_info + ' >' + all_mC_site_info_extend
    logger.info(make_cmd)
    os.system(make_cmd)


def get_fasta_for_c_extend(all_mC_site_info_extend, all_mC_site_info_extend_fa):
    logger.info('get c site extend sequence for {0} ...'.format(all_mC_site_info_extend))
    get_extend_seq_cmd = 'bedtools getfasta -fi /home/sunwenju/3_genomeData/human/hg19/hg19.fa -bed '+all_mC_site_info_extend+' -fo '+all_mC_site_info_extend_fa+' -name -tab -s'
    logger.info(get_extend_seq_cmd)
    os.system(get_extend_seq_cmd)


def make_mc_class_result(c_site_extend_fa, c_site_info, all_CpG_result, all_GpC_result, mini_coverage):
    logger.info('make CpG and GpC site info output result for {0} ...'.format(c_site_info))
    cpg_site_dict = {}
    gpc_site_dict = {}
    with open(c_site_extend_fa, 'r') as cfin:
        for line in cfin:
            site_name, sequence = line.rstrip().split("\t")
            sequence = sequence.upper()
            # CpG = WCG = ACG + TCG = [AT]CG
            # GpC = GCH = GCA + GCC + GCT = GC[ACT]
            if sequence == 'ACG' or sequence == 'TCG':
                cpg_site_dict[site_name] = sequence
            elif sequence == 'GCA' or sequence == 'GCC' or sequence == 'GCT':
                gpc_site_dict[site_name] = sequence
    with open(c_site_info, 'r') as csin, open(all_CpG_result, 'w') as cpg_out, open(all_GpC_result, 'w') as gpc_out:
        m = 0
        n = 0
        for line in csin:
            items = line.rstrip().split("\t")
            chrom = items[0]
            start = items[1]
            end = items[2]
            site_name = items[3]
            mc_site_count = int(items[4])
            strand = items[5]
            unmc_site_count = int(items[6])
            coverage = mc_site_count + unmc_site_count
            if coverage >= mini_coverage:
                m_percentage = "{:.4f}".format(mc_site_count/coverage)
                if site_name in cpg_site_dict:
                    m += 1
                    new_name = 'CpG_site_'+str(m)
                    print(chrom, start, end, new_name, m_percentage, strand, cpg_site_dict[site_name], coverage, mc_site_count, unmc_site_count, sep="\t", end="\n", file=cpg_out)
                if site_name in gpc_site_dict:
                    n += 1
                    new_name = 'GpC_site_' + str(n)
                    print(chrom, start, end, new_name, m_percentage, strand, gpc_site_dict[site_name], coverage,
                          mc_site_count, unmc_site_count, sep="\t", end="\n", file=gpc_out)


def my_parse_args():
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-i', '--input', type=str, required=True, help='bismark methylation extractor genome wide all C site stat result')
    my_parser.add_argument('-m', '--mincov', type=int, default=6, help='mini reads coverage for a C site')
    my_parser.add_argument('-o', '--output', type=str, required=True, help='prefix str for output file, base name')
    return my_parser.parse_args()


def main():
    my_args = my_parse_args()
    input_file = my_args.input
    mini_coverage = my_args.mincov
    base_name = my_args.output

    logger.info('command line: python3 make_CpG_GpC_mC_site_info_from_bme_results.py -i {0} -m {1} -o {2}'.format(input_file, mini_coverage, base_name))

    # filter 0 mC site and convert bismark methylation extractor genome wide all C site stat result to bed6+
    # <chromosome>	<start>	<end>	<name>	<count methylated>	<strand>	<count unmethylated>
    all_mC_site_info = base_name + '.all_mC_site_info.txt'
    bme_to_filterd_bed(input_file, all_mC_site_info)

    # get fasta for all_mC_site_info
    all_mC_site_info_fa = base_name + '.all_mC_site_info.txt.fa'
    get_fasta_for_c_site(all_mC_site_info, all_mC_site_info_fa)

    # make all mC site info extend +- 1nt
    all_mC_site_info_extend = base_name + '.all_mC_site_info.entend.txt'
    make_c_site_extend_bed(all_mC_site_info, all_mC_site_info_extend)

    # get fasta for all_mC_site_info_extend
    all_mC_site_info_extend_fa = base_name + '.all_mC_site_info.entend.txt.fa'
    get_fasta_for_c_extend(all_mC_site_info_extend, all_mC_site_info_extend_fa)

    # make GpC and CpG output result
    # CpG = WCG = ACG + TCG = [AT]CG
    # GpC = GCH = GCA + GCC + GCT = GC[ACT]
    all_CpG_result = base_name + '.all_CpG_result.txt'
    all_GpC_result = base_name + '.all_GpC_result.txt'
    make_mc_class_result(all_mC_site_info_extend_fa, all_mC_site_info, all_CpG_result, all_GpC_result, mini_coverage)
    logger.info('Run completed!')


if __name__ == '__main__':
    main()
