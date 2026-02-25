

import argparse
from reportlab.lib import colors
from reportlab.pdfgen import canvas
from reportlab.lib.units import mm
from reportlab.graphics.shapes import *
from reportlab.graphics import renderPDF
from reportlab.lib.colors import *


def my_parse_args():
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-a', '--plot', type=str, required=True, help='ctcf region info for plot')
    my_parser.add_argument('-o', '--output', type=str, required=True, help='output pdf plot file name')
    return my_parser.parse_args()


def main():
    my_args = my_parse_args()
    info_file = my_args.plot
    pdf_file = my_args.output

    g_seq = ''
    read_pos_list = []
    read_seq_list = []
    read_strand_list = []
    cpg_mc_pos_list = []
    gpc_mc_pos_list = []


    with open(info_file, 'r') as fin:
        for line in fin:
            items = line.rstrip().split("\t")
            g_seq = items[1]
            read_pos_list = items[3].split(';')
            read_seq_list = items[4].split(';')
            read_strand_list = items[6].split(';')
            cpg_mc_pos_list = items[7].split(';')
            gpc_mc_pos_list = items[8].split(';')
            # print("cpg_mc_pos_list:", cpg_mc_pos_list)
            # print("gpc_mc_pos_list:", gpc_mc_pos_list)

    grey20transparent = Color( 0.5, 0.5, 0.5, alpha=0.2)
    red20transparent = Color(1, 0, 0, alpha=0.2)
    red70transparent = Color(1, 0, 0, alpha=0.7)
    blue20transparent = Color( 0, 0, 1, alpha=0.2)
    blue70transparent = Color( 0, 0, 1, alpha=0.7)
    #
    cLen = 4.8
    bBind = 9
    tBind = 10
    page_len = 1444.8
    page_height = 10*(1+len(read_seq_list))+1

    cv = canvas.Canvas(pdf_file, pagesize=(page_len, page_height))
    # plot all sequence
    cv.setFont("Courier", 8)
    cv.drawString(0, page_height-8, g_seq)
    for i in range(len(read_seq_list)):
        cv.drawString(float(read_pos_list[i])*cLen, page_height-8-(i+1)*tBind, read_seq_list[i])
    # plot all background fill color
    cv.setFillColor(grey20transparent)
    cv.rect(0, page_height - 10, len(g_seq) * cLen, bBind, stroke=False, fill=True)
    for i in range(len(read_seq_list)):
        if read_strand_list[i] == '+':
            cv.setFillColor(red20transparent)
            cv.rect(float(read_pos_list[i])*cLen, page_height-10-(i+1)*tBind, len(read_seq_list[i]) * cLen, bBind, stroke=False, fill=True)
        else:
            cv.setFillColor(blue20transparent)
            cv.rect(float(read_pos_list[i])*cLen, page_height-10-(i+1)*tBind, len(read_seq_list[i]) * cLen, bBind, stroke=False, fill=True)

    # plot all cpg mc sites
    for i in range(len(read_seq_list)):
        read_cpg_mc_pos_list = cpg_mc_pos_list[i].split(',')
        for cpg_mc_pos in read_cpg_mc_pos_list:
            if cpg_mc_pos != '':
                # print(i+1, 'cpg_mc_pos:', cpg_mc_pos, sep="\t", end="\n")
                cv.setFillColor(red70transparent)
                cv.rect((int(read_pos_list[i])+int(cpg_mc_pos)-1)*cLen, page_height-10-(i+1)*tBind, 1*cLen, bBind, stroke=False, fill=True)
                # print((int(read_pos_list[i])+int(cpg_mc_pos)-1)*cLen, page_height-10-(i+1)*tBind, 1*cLen, bBind,)

    # plot all gpc mc sites
    for i in range(len(read_seq_list)):
        read_gpc_mc_pos_list = gpc_mc_pos_list[i].split(',')
        for gpc_mc_pos in read_gpc_mc_pos_list:
            if gpc_mc_pos != '':
                # print(i+1, 'gpc_mc_pos:', gpc_mc_pos, sep="\t", end="\n")
                cv.setFillColor(blue70transparent)
                cv.rect((int(read_pos_list[i])+int(gpc_mc_pos)-1)*cLen, page_height-10-(i+1)*tBind, 1*cLen, bBind, stroke=False, fill=True)
                # print((int(read_pos_list[i])+int(gpc_mc_pos)-1)*cLen, page_height-10-(i+1)*tBind, 1*cLen, bBind,)

    cv.save()


if __name__ == '__main__':
    main()

