

import argparse
from reportlab.lib import colors
from reportlab.pdfgen import canvas
from reportlab.lib.units import mm
from reportlab.graphics.shapes import *
from reportlab.graphics import renderPDF
from reportlab.lib.colors import *


def get_rev_com_seq(upper_seq_str):
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
    rev_seq_list = list(reversed(upper_seq_str))
    rev_com_seq_list = [complement_dict[k] for k in rev_seq_list]
    rev_com_seq = ''.join(rev_com_seq_list)
    return rev_com_seq


def get_com_seq(upper_seq_str):
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
    com_seq_list = [complement_dict[k] for k in list(upper_seq_str)]
    com_seq = ''.join(com_seq_list)
    return com_seq


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

    balck = Color( 0, 0, 0, alpha=1)
    grey20transparent = Color( 0.5, 0.5, 0.5, alpha=0.2)
    red20transparent = Color(1, 0, 0, alpha=0.2)
    red70transparent = Color(1, 0, 0, alpha=0.7)
    green70transparent = Color(0, 1, 0, alpha=0.7)
    blue20transparent = Color( 0, 0, 1, alpha=0.2)
    blue70transparent = Color( 0, 0, 1, alpha=0.7)
    #
    cLen = 4.8
    bBind = 9
    tBind = 10
    page_len = 1444.8
    page_height = 10*(2+len(read_seq_list))+1

    cv = canvas.Canvas(pdf_file, pagesize=(page_len, page_height))
    # plot all sequence
    cv.setFont("Courier", 8)
    cv.setFillColor(balck)
    cv.drawString(0, page_height-8, g_seq)
    cv.drawString(0, page_height-8-tBind, get_com_seq(g_seq))
    for i in range(len(read_seq_list)):
        if read_strand_list[i] == '+':
            cv.drawString(float(read_pos_list[i])*cLen, page_height-8-(i+2)*tBind, read_seq_list[i])
        elif read_strand_list[i] == '-':
            cv.drawString(float(read_pos_list[i])*cLen, page_height-8-(i+2)*tBind, get_com_seq(read_seq_list[i]))
    
    # plot CTCF summit
    cv.setFillColor(green70transparent)
    cv.rect(float(150)*cLen, page_height - 10, cLen, bBind, stroke=False, fill=True)
    cv.rect(float(150)*cLen, page_height - 20, cLen, bBind, stroke=False, fill=True)

    # make background fill '<>'
    cv.setFont("Courier", 8)
    cv.setFillColor(grey20transparent)
    for i in range(len(read_seq_list)):
        start = int(read_pos_list[i])+len(read_seq_list[i])
        str_len = 301 - start
        if read_strand_list[i] == '+':
            cv.drawString(0, page_height-8-(i+2)*tBind, '>'*int(read_pos_list[i]))
            cv.drawString(float(start)*cLen, page_height-8-(i+2)*tBind, '>'*int(str_len))
        elif read_strand_list[i] == '-':
            cv.drawString(0, page_height-8-(i+2)*tBind, '<'*int(read_pos_list[i]))
            cv.drawString(float(start)*cLen, page_height-8-(i+2)*tBind, '<'*int(str_len))
    
    # plot all background fill color
    cv.setFillColor(red20transparent)
    cv.rect(0, page_height - 10, len(g_seq) * cLen, bBind, stroke=False, fill=True)
    cv.setFillColor(blue20transparent)
    cv.rect(0, page_height - 20, len(g_seq) * cLen, bBind, stroke=False, fill=True)
    for i in range(len(read_seq_list)):
        if read_strand_list[i] == '+':
            cv.setFillColor(red20transparent)
            cv.rect(float(read_pos_list[i])*cLen, page_height-10-(i+2)*tBind, len(read_seq_list[i]) * cLen, bBind, stroke=False, fill=True)
        else:
            cv.setFillColor(blue20transparent)
            cv.rect(float(read_pos_list[i])*cLen, page_height-10-(i+2)*tBind, len(read_seq_list[i]) * cLen, bBind, stroke=False, fill=True)
    
    # plot all cpg mc sites
    for i in range(len(read_seq_list)):
        read_cpg_mc_pos_list = cpg_mc_pos_list[i].split(',')
        for cpg_mc_pos in read_cpg_mc_pos_list:
            if cpg_mc_pos != '':
                # print(i+2, 'cpg_mc_pos:', cpg_mc_pos, sep="\t", end="\n")
                cv.setFillColor(red70transparent)
                cv.rect((int(read_pos_list[i])+int(cpg_mc_pos)-1)*cLen, page_height-10-(i+2)*tBind, 1*cLen, bBind, stroke=False, fill=True)
                # print((int(read_pos_list[i])+int(cpg_mc_pos)-1)*cLen, page_height-10-(i+2)*tBind, 1*cLen, bBind,)

    # plot all gpc mc sites
    for i in range(len(read_seq_list)):
        read_gpc_mc_pos_list = gpc_mc_pos_list[i].split(',')
        for gpc_mc_pos in read_gpc_mc_pos_list:
            if gpc_mc_pos != '':
                # print(i+2, 'gpc_mc_pos:', gpc_mc_pos, sep="\t", end="\n")
                cv.setFillColor(blue70transparent)
                cv.rect((int(read_pos_list[i])+int(gpc_mc_pos)-1)*cLen, page_height-10-(i+2)*tBind, 1*cLen, bBind, stroke=False, fill=True)
                # print((int(read_pos_list[i])+int(gpc_mc_pos)-1)*cLen, page_height-10-(i+2)*tBind, 1*cLen, bBind,)

    cv.save()

if __name__ == '__main__':
    main()

