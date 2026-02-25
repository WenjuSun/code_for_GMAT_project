#!/home/sunwenju/1_softwares/Python-3.7.0/bin/python3
# -*- encoding: utf-8 -*-
'''
@File    :   0_GpmC_peak_calling.v3.py
@Time    :   2020/11/13 16:53:11
@Author  :   Wenju Sun
@Version :   1.0
@Contact :   wenju.sun@qq.com
@License :   Copyright (C) - All Rights Reserved Sun Wenju
@Desc    :   None
'''

import os
import re
import logging
import argparse
import subprocess
import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import Counter
from multiprocessing import Pool
from scipy.stats import poisson
from statsmodels.stats.multitest import fdrcorrection


def create_custom_logger(name):
    formatter = logging.Formatter(fmt='[%(asctime)s][%(levelname)s][%(module)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    return logger

logger = create_custom_logger('root')


def make_ext_bed(raw_bed, extB_len, ext_bed):
    """
    docstring
    """
    logger.info('making extend bed for {0}, pid:{1}, ppid:{2} ...'.format(ext_bed, os.getpid(), os.getppid()))
    cmd = 'bedtools slop -i {0} -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes -b {1} |sort -k1,1 -k2,2n >{2}'.format(raw_bed, extB_len, ext_bed)
    subprocess.call(cmd, shell=True)


def get_fasta_from_bed(ext_bed, ext_fasta):
    """
    docstring
    """
    logger.info('getting fasta for {0}, pid:{1}, ppid:{2} ...'.format(ext_bed, os.getpid(), os.getppid()))
    cmd = 'bedtools getfasta -fi /home/sunwenju/3_genomeData/human/hg19/hg19.fa -bed {0} -fo {1} -name -tab'.format(ext_bed, ext_fasta)
    subprocess.call(cmd, shell=True)


def making_intersect_results(extBn_bed, all_raw_gpmc_bed, intersect_results):
    """
    docstring
    """
    logger.info('making intersect results for {0}, pid:{1}, ppid:{2} ...'.format(extBn_bed, os.getpid(), os.getppid()))
    cmd = 'bedtools intersect -a {0} -b {1} -wao >{2}'.format(extBn_bed, all_raw_gpmc_bed, intersect_results)
    subprocess.call(cmd, shell=True)


def get_genome_gpc_count(genome_fasta):
    """
    input  : genome sequence
    return : genome GpC site count (with reverse complement strand sequence stat)
    """
    logger.info('getting genome all gpc count ...')
    genome_gpc_count = 0
    for seq_record in SeqIO.parse(genome_fasta, "fasta"):
        sequence = seq_record.seq.upper()
        rev_comp_seq = sequence.reverse_complement()
        genome_gpc_count += sequence.count('GC')
        genome_gpc_count -= sequence.count('GCG')
        genome_gpc_count += rev_comp_seq.count('GC')
        genome_gpc_count -= rev_comp_seq.count('GCG')
    return genome_gpc_count


def get_all_gpmc_count(raw_gpmc_file):
    """
    input : all GpmC sites info file
    return: genome GpmC sites count
    """
    logger.info('getting all gpmc count ...')
    gpmc_count = 0
    with open(raw_gpmc_file, 'r') as raw_in:
        gpmc_count = len(raw_in.readlines()) # if file bigger than 10G, try other methods
    return gpmc_count


def get_gpc_count_from_fasta(tab_fasta_file, gpc_dict):
    """
    docstring
    """
    logger.info('making gpc dict for {0}, pid:{1}, ppid:{2} ...'.format(tab_fasta_file, os.getpid(), os.getppid()))
    gpc_dict = {}
    with open(tab_fasta_file, 'r') as tfin:
        for line in tfin:
            seq_name, seq_str = line.rstrip().split("\t")
            sequence = seq_str.upper()
            rev_comp_seq = get_reverse_complement(sequence)
            gpc_dict[seq_name] = 0
            gpc_dict[seq_name] += sequence.count('GC')
            gpc_dict[seq_name] -= sequence.count('GCG')
            gpc_dict[seq_name] += rev_comp_seq.count('GC')
            gpc_dict[seq_name] -= rev_comp_seq.count('GCG')
    return gpc_dict


def get_reverse_complement(upper_seq_str):
    """
    docstring
    """
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
    rev_seq_list = list(reversed(upper_seq_str))
    rev_com_seq_list = [complement_dict[k] for k in rev_seq_list]
    rev_com_seq = ''.join(rev_com_seq_list)
    return rev_com_seq


def get_gpmc_count_from_intersect(intersect_result, gpmc_dict):
    """
    docstring
    """
    logger.info('making gpmc dict for {0}, pid:{1}, ppid:{2} ...'.format(intersect_result, os.getpid(), os.getppid()))
    # file format, nohead, 
    # <ext_chrom> <ext_start> <ext_end> <ext_name> <mc_chrom> <mc_start> <mc_end> <mc_site_name> <mc_ratio> <mc_strand> <mc_type> <cov> <mc_rc> <unmc_rc> <overlap>
    gpmc_dict = {}
    my_data = pd.read_csv(intersect_result, sep="\t", header=None, usecols=[3])
    gpmc_dict = dict(Counter(my_data[3].values))
    return gpmc_dict


def get_gpmc_pv_dict(gpc_dict, gpmc_dict, data_set):
    """
    docstring
    """
    logger.info('making gpmc p value dict for {0}, pid:{1}, ppid:{2}'.format(data_set, os.getpid(), os.getppid()))
    gpmc_pv_dict = {}
    for name in gpc_dict:
        gpmc_pv = format(float(gpmc_dict[name])/float(gpc_dict[name]), '.8f')
        gpmc_pv_dict[name] = gpmc_pv
    return gpmc_pv_dict


def my_parse_args():
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-g', '--genome', type=str, default="/home/sunwenju/3_genomeData/human/hg19/hg19.fa", help='input genome fasta file')
    my_parser.add_argument('-r', '--raw_gpmc', type=str, required=True, help='raw GpmC file with all GpmC sites info')
    my_parser.add_argument('-m', '--merged_gpmc', type=str, required=True, help='merged GpmC file with candidate peaks')
    my_parser.add_argument('-b', '--merged_peaks_bed', type=str, required=True, help='merged GpmC peaks bed4')
    my_parser.add_argument('-o', '--output', type=str, required=True, help='output file name')
    return my_parser.parse_args()


def main():
    my_args = my_parse_args()
    genome_fasta = my_args.genome
    raw_gpmc_file = my_args.raw_gpmc
    merged_gpmc_file = my_args.merged_gpmc
    mp_bed = my_args.merged_peaks_bed
    output_file = my_args.output
    out_with_fdr = output_file.replace(".txt", "") + '.withFDR.txt'
    logger.info('command line: python3 make_GpC_peaks_stat_info_new2.py -g {0} -r {1} -m {2} -f {3} -o {4}'.format(genome_fasta, raw_gpmc_file, merged_gpmc_file, mp_bed, output_file))

    base_name = mp_bed.replace(".bed", "")
    mp_fasta = mp_bed + ".fa"

    l3_bed = base_name + ".exto3L.bed"
    l3_fasta = base_name + ".exto3L.bed.fa"
    l3_intersect = base_name + ".exto3L.bed.intersect.txt"
    logger.info(l3_bed+"\t"+l3_fasta+"\t"+l3_intersect)

    k5_bed = base_name + ".extB5k.bed"
    k5_fasta = base_name + ".extB5k.bed.fa"
    k5_intersect = base_name + ".extB5k.bed.intersect.txt"
    logger.info(k5_bed+"\t"+k5_fasta+"\t"+k5_intersect)

    k10_bed = base_name + ".extB10k.bed"
    k10_fasta = base_name + ".extB10k.bed.fa"
    k10_intersect = base_name + ".extB10k.bed.intersect.txt"
    logger.info(k10_bed+"\t"+k10_bed+"\t"+k10_intersect)

    k50_bed = base_name + ".extB50k.bed"
    k50_fasta = base_name + ".extB50k.bed.fa"
    k50_intersect = base_name + ".extB50k.bed.intersect.txt"
    logger.info(k50_bed+"\t"+k50_fasta+"\t"+k50_intersect)

    # make extend bed
    logger.info('making extend bed for exto3L, pid:{0}, ppid:{1} ...'.format(os.getpid(), os.getppid()))
    cmd = 'bedtools slop -i {0} -g /home/sunwenju/3_genomeData/human/hg19/hg19.chrom.sizes -b 1.0 -pct |sort -k1,1 -k2,2n >{1}'.format(mp_bed, l3_bed)
    subprocess.call(cmd, shell=True)
    ## extB 5k, 10k, 50k
    with Pool(3) as pool:
        mp_extend_result = pool.starmap_async(make_ext_bed, [(mp_bed, 5000, k5_bed), (mp_bed, 10000, k10_bed), (mp_bed, 50000, k50_bed)])
        mp_extend_result.wait()
    try:
        bed_files = mp_extend_result.get()
    except TimeoutError:
        logger.error(TimeoutError.__name__)

    # get fasta
    with Pool(5) as pool:
        ext_fa_result = pool.starmap_async(get_fasta_from_bed, [(mp_bed, mp_fasta), (l3_bed, l3_fasta), (k5_bed, k5_fasta), (k10_bed, k10_fasta), (k50_bed, k50_fasta)])
        ext_fa_result.wait()
    try:
        fasta_files = ext_fa_result.get()
    except TimeoutError:
        logger.error(TimeoutError.__name__)

    # intesect
    with Pool(4) as pool:
        intersect_result = pool.starmap_async(making_intersect_results, [(l3_bed, raw_gpmc_file, l3_intersect), (k5_bed, raw_gpmc_file, k5_intersect), (k10_bed, raw_gpmc_file, k10_intersect), (k50_bed, raw_gpmc_file, k50_intersect)])
        intersect_result.wait()
    try:
        intersect_files = intersect_result.get()
    except TimeoutError:
        logger.error(TimeoutError.__name__)

    # get gpc count
    genome_gpc_count = get_genome_gpc_count(genome_fasta)
    # genome_gpc_count = 233643193
    mp_fa_gpc_dict = {}
    l3_fa_gpc_dict = {}
    k5_fa_gpc_dict = {}
    k10_fa_gpc_dict = {}
    k50_fa_gpc_dict = {}
    with Pool(5) as pool:
        get_gpc_result = pool.starmap_async(get_gpc_count_from_fasta, [(mp_fasta, mp_fa_gpc_dict), (l3_fasta, l3_fa_gpc_dict), (k5_fasta, k5_fa_gpc_dict), (k10_fasta, k10_fa_gpc_dict), (k50_fasta, k50_fa_gpc_dict)])
        get_gpc_result.wait()
    try:
        mp_fa_gpc_dict, l3_fa_gpc_dict, k5_fa_gpc_dict, k10_fa_gpc_dict, k50_fa_gpc_dict = get_gpc_result.get()
    except TimeoutError:
        logger.error(TimeoutError.__name__)

    # get gpmc count
    all_gpcm_count = get_all_gpmc_count(raw_gpmc_file)
    # all_gpcm_count = 506580
    l3_intersect_gpmc_dict = {}
    k5_intersect_gpmc_dict = {}
    k10_intersect_gpmc_dict = {}
    k50_intersect_gpmc_dict = {}
    with Pool(4) as pool:
        get_gpmc_result = pool.starmap_async(get_gpmc_count_from_intersect, [(l3_intersect, l3_intersect_gpmc_dict), (k5_intersect, k5_intersect_gpmc_dict), (k10_intersect, k10_intersect_gpmc_dict), (k50_intersect, k50_intersect_gpmc_dict)])
        get_gpmc_result.wait()
    try:
        l3_intersect_gpmc_dict, k5_intersect_gpmc_dict, k10_intersect_gpmc_dict, k50_intersect_gpmc_dict = get_gpmc_result.get()
    except TimeoutError:
        logger.error(TimeoutError.__name__)

    # get probability
    genome_gpmc_prob = float(format(all_gpcm_count/genome_gpc_count, '.8f'))
    l3_gpmc_pv_dict = {}
    k5_gpmc_pv_dict = {}
    k10_gpmc_pv_dict = {}
    k50_gpmc_pv_dict = {}
    with Pool(4) as pool:
        get_gpmc_pv_result = pool.starmap_async(get_gpmc_pv_dict, [(l3_fa_gpc_dict, l3_intersect_gpmc_dict, 'l3'), (k5_fa_gpc_dict, k5_intersect_gpmc_dict, 'k5'), (k10_fa_gpc_dict, k10_intersect_gpmc_dict, 'k10'), (k50_fa_gpc_dict, k50_intersect_gpmc_dict, 'k50')])
        get_gpmc_pv_result.wait()
    try:
        l3_gpmc_pv_dict, k5_gpmc_pv_dict, k10_gpmc_pv_dict, k50_gpmc_pv_dict = get_gpmc_pv_result.get()
    except TimeoutError:
        logger.error(TimeoutError.__name__)

    # prepare for output
    logger.info('processing and making stat info for output ...')
    with open(merged_gpmc_file, 'r') as merge_in, open(output_file, 'w') as fout:
        # print output title
        print('chrom', 'peak_start', 'peak_end', 'peak_name', 'peak_len', 'peak_strand', 'site_name_list', 'strand_list', 
            'type_list', 'mean_cov', 'mean_mc_rc', 'mean_unmc_rc', 'mean_mc_ratio', 'peak_gpc_count', 'peak_gpmc_count', 
            'genome_gpc_count', 'all_gpcm_count', 'genome_gpmc_prob', 'lambda_genome', 'peak_pv_based_lambda_genome', 
            'l3_fa_gpc', 'l3_intersect_gpmc', 'l3_gpmc_prob', 'lambda_l3', 'peak_pv_based_lambda_l3', 
            'k5_fa_gpc', 'k5_intersect_gpmc', 'k5_gpmc_prob', 'lambda_k5', 'peak_pv_based_lambda_k5', 
            'k10_fa_gpc', 'k10_intersect_gpmc', 'k10_gpmc_prob', 'lambda_k10', 'peak_pv_based_lambda_k10', 
            'k50_fa_gpc', 'k50_intersect_gpmc', 'k50_gpmc_prob', 'lambda_k50', 'peak_pv_based_lambda_k50', 
            'lambda_max', 'max_lambda_from', 'peak_pv_based_lambda_max', sep="\t", end="\n", file=fout)
        for line in merge_in:
            items = line.rstrip().split("\t")
            # <chrom> <peak_start> <peak_end> <peak_name> <peak_len> <peak_strand> <site_name_list> <strand_list> <type_list> <mean_cov> <mean_mc_rc> <mean_unmc_rc> <mean_mc_ratio> <mc_count>
            name = items[3]
            mc_count = int(items[13])
            peak_gpc_count = int(mp_fa_gpc_dict[name])
            l3_gpmc_prob = float(l3_gpmc_pv_dict[name])
            k5_gpmc_prob = float(k5_gpmc_pv_dict[name])
            k10_gpmc_prob = float(k10_gpmc_pv_dict[name])
            k50_gpmc_prob = float(k50_gpmc_pv_dict[name])
            # get lambda
            lambda_genome = float(peak_gpc_count * genome_gpmc_prob)
            lambda_l3 = float(peak_gpc_count * l3_gpmc_prob)
            lambda_k5 = float(peak_gpc_count * k5_gpmc_prob)
            lambda_k10 = float(peak_gpc_count * k10_gpmc_prob)
            lambda_k50 = float(peak_gpc_count * k50_gpmc_prob)
            lambda_info = ['g', 'k5', 'k10', 'k50']
            lambda_list = [lambda_genome, lambda_k5, lambda_k10, lambda_k50]
            lambda_max = max(lambda_list)
            max_lambda_from = lambda_info[lambda_list.index(lambda_max)]
            #
            peak_pv_based_lambda_genome = poisson.sf(mc_count-1, lambda_genome)
            peak_pv_based_lambda_l3 = poisson.sf(mc_count-1, lambda_l3)
            peak_pv_based_lambda_k5 = poisson.sf(mc_count-1, lambda_k5)
            peak_pv_based_lambda_k10 = poisson.sf(mc_count-1, lambda_k10)
            peak_pv_based_lambda_k50 = poisson.sf(mc_count-1, lambda_k50)
            peak_pv_based_lambda_max = poisson.sf(mc_count-1, lambda_max)
            # output
            basic_info = '\t'.join(items[0:13])
            print(basic_info, peak_gpc_count, mc_count, 
                genome_gpc_count, all_gpcm_count, genome_gpmc_prob, lambda_genome, peak_pv_based_lambda_genome, 
                l3_fa_gpc_dict[name], l3_intersect_gpmc_dict[name], l3_gpmc_prob, lambda_l3, peak_pv_based_lambda_l3,
                k5_fa_gpc_dict[name], k5_intersect_gpmc_dict[name], k5_gpmc_prob, lambda_k5, peak_pv_based_lambda_k5,
                k10_fa_gpc_dict[name], k10_intersect_gpmc_dict[name], k10_gpmc_prob, lambda_k10, peak_pv_based_lambda_k10,
                k50_fa_gpc_dict[name], k50_intersect_gpmc_dict[name], k50_gpmc_prob, lambda_k50, peak_pv_based_lambda_k50,
                lambda_max, max_lambda_from, peak_pv_based_lambda_max,
                sep="\t", end="\n", file=fout)

    # add fdr output
    logger.info('add fdr info for results output')
    df = pd.read_csv(output_file, sep='\t')
    # ref：https://www.statsmodels.org/devel/generated/statsmodels.stats.multitest.fdrcorrection.html#statsmodels.stats.multitest.fdrcorrection
    genome_rejected, genome_fdr = fdrcorrection(np.array(df.peak_pv_based_lambda_genome))
    l3_rejected, l3_fdr = fdrcorrection(np.array(df.peak_pv_based_lambda_l3))
    k5_rejected, k5_fdr = fdrcorrection(np.array(df.peak_pv_based_lambda_k5))
    k10_rejected, k10_fdr = fdrcorrection(np.array(df.peak_pv_based_lambda_k10))
    k50_rejected, k50_fdr = fdrcorrection(np.array(df.peak_pv_based_lambda_k50))
    lmd_max_rejected, lmd_max_fdr = fdrcorrection(np.array(df.peak_pv_based_lambda_max))
    # 
    col_name=df.columns.tolist()
    col_name.insert(20, 'genome_fdr')
    col_name.insert(21, 'genome_rejected')
    col_name.insert(27, 'l3_fdr')
    col_name.insert(28, 'l3_rejected')
    col_name.insert(34, 'k5_fdr')
    col_name.insert(35, 'k5_rejected')
    col_name.insert(41, 'k10_fdr')
    col_name.insert(42, 'k10_rejected')
    col_name.insert(48, 'k50_fdr')
    col_name.insert(49, 'k50_rejected')
    col_name.insert(53, 'lmd_max_fdr')
    col_name.insert(54, 'lmd_max_rejected')
    # 
    df=df.reindex(columns=col_name)
    df['genome_fdr']=genome_fdr
    df['genome_rejected']=genome_rejected
    df['l3_fdr']=l3_fdr
    df['l3_rejected']=l3_rejected
    df['k5_fdr']=k5_fdr
    df['k5_rejected']=k5_rejected
    df['k10_fdr']=k10_fdr
    df['k10_rejected']=k10_rejected
    df['k50_fdr']=k50_fdr
    df['k50_rejected']=k50_rejected
    df['lmd_max_fdr']=lmd_max_fdr
    df['lmd_max_rejected']=lmd_max_rejected
    df.to_csv(out_with_fdr, sep='\t', header=True, index=False)

    logger.info('Run completed!')


if __name__ == '__main__':
    main()
