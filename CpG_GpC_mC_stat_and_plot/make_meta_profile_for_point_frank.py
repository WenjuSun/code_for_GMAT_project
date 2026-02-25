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


def my_parse_args():
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-a', '--anno', type=str, required=True, help='point annotation bed file')
    my_parser.add_argument('-b', '--bme', type=str, required=True, help='CpG or GpC methylation site info bed plus file')
    my_parser.add_argument('-c', '--minicover', type=int, default=1, help='mini reads coverage for a C site')
    my_parser.add_argument('-d', '--distance', type=int, default=5000, help='max distance to the point annotation')
    my_parser.add_argument('-o', '--output', type=str, required=True, help='statistical result output file name')
    return my_parser.parse_args()


def main():
    my_args = my_parse_args()
    point_anno = my_args.anno
    bme_mc = my_args.bme
    mini_cover = my_args.minicover
    max_distance = my_args.distance
    output_file = my_args.output

    logger.info('command line: python3 make_meta_profile_for_point_frank.py -a {0} -b {1} -c {2} -d {3} -o {4}'.format(point_anno, bme_mc, mini_cover, max_distance, output_file))

    anno_dict = {}
    with open(point_anno, 'r') as afin:
        for line in afin:
            chrom, start, end = line.rstrip().split("\t")
            if chrom not in anno_dict:
                anno_dict[chrom] = []
                anno_dict[chrom].append(int(end))
            else:
                anno_dict[chrom].append(int(end))

    d_stat_dict = {}
    for distance in range(-max_distance, max_distance+1):
        d_stat_dict[distance] = {}
        d_stat_dict[distance]['m_percentage'] = []
        d_stat_dict[distance]['coverage'] = []
        d_stat_dict[distance]['mc_site_count'] = []
        d_stat_dict[distance]['unmc_site_count'] = []

    with open(bme_mc, 'r') as bfin:
        for line in bfin:
            # chrom, start, end, name, m_percentage, strand, context, coverage, mc_site_count, unmc_site_count
            chrom, start, end, name, m_percentage, strand, context, coverage, mc_site_count, unmc_site_count = line.rstrip().split("\t")
            if int(coverage) < mini_cover:
                continue
            if chrom in anno_dict:
                for summit in anno_dict[chrom]:
                    dist_to_summit = int(end) - summit
                    if - max_distance <= dist_to_summit <= max_distance:
                            d_stat_dict[dist_to_summit]['m_percentage'].append(float(m_percentage))
                            d_stat_dict[dist_to_summit]['coverage'].append(int(coverage))
                            d_stat_dict[dist_to_summit]['mc_site_count'].append(int(mc_site_count))
                            d_stat_dict[dist_to_summit]['unmc_site_count'].append(int(unmc_site_count))

    with open(output_file, 'w') as out:
        print('distance', 'stack_n', 'avg_mC_ratio', 'avg_coverage', 'avg_mC_count', 'avg_unmC_count', 'mC_ratio_list', 'coverage_list', 'mC_count_list', 'unmC_count_list', sep="\t", end="\n", file=out)
        for distance in range(-max_distance, max_distance+1):
            n = len(d_stat_dict[distance]['m_percentage'])
            avg_m_ratio = 0
            avg_cover = 0
            avg_mc_sc = 0
            avg_unmc_sc = 0
            if n != 0:
                avg_m_ratio = format(sum(d_stat_dict[distance]['m_percentage'])/n, '.4f')
                avg_cover = format(sum(d_stat_dict[distance]['coverage']) / n, '.4f')
                avg_mc_sc = format(sum(d_stat_dict[distance]['mc_site_count']) / n, '.4f')
                avg_unmc_sc = format(sum(d_stat_dict[distance]['unmc_site_count']) / n, '.4f')

            print(distance, n, avg_m_ratio, avg_cover, avg_mc_sc, avg_unmc_sc,
                      d_stat_dict[distance]['m_percentage'], d_stat_dict[distance]['coverage'],
                      d_stat_dict[distance]['mc_site_count'], d_stat_dict[distance]['unmc_site_count'], sep="\t", end="\n", file=out)

    logger.info('Run completed!')


if __name__ == '__main__':
    main()
