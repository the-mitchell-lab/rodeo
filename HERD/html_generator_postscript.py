#==============================================================================
# Copyright (C) 2017 Bryce L. Kille
# University of Illinois
# Department of Chemistry
#
# Copyright (C) 2017 Christopher J. Schwalen
# University of Illinois
# Department of Chemistry
#
# Copyright (C) 2026 Shravan R. Dommaraju
# Vanderbilt University
# Department of Biochemistry
#
# Copyright (C) 2026 Douglas A. Mitchell
# Vanderbilt University
# Department of Biochemistry
#
# License: GNU Affero General Public License v3 or later
# Complete license availabel in the accompanying LICENSE.txt.
# or <http://www.gnu.org/licenses/>.
#
# This file is part of RODEO.
#
# RODEO is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RODEO is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#==============================================================================

import argparse, sys, os, re, csv, glob, time
import xml.parsers.expat
import numpy as np
import pandas as pd
import math
from Bio import Seq, SeqIO
from decimal import Decimal

def write_record(main_html, record, filter_list, print_precursors, curr_query="", peptide_type=""):
    #RESULTS FOR xxxx
    #DRAW ORF
    #DRAW ORF scale
    #PUT LINK TO nuc SEQUENCE
    #TABLE of CDS
    #TABLE of ORFs
    if record.locus is not -1:
        main_html.write('<h2 id="%s"> Results for Locus %s in Record %s [%s]\n' % (record.cluster_accession, record.locus, record.cluster_accession, record.cluster_genus_species))
    else:
        main_html.write('<h2 id="%s"> Results for %s [%s]\n' % (record.query_accession_id, record.query_accession_id, record.cluster_genus_species))
    main_html.write('<a href="#header"><small><small>back to top</small></small></a></h2>') #TODO keep for single?
    main_html.write('<p></p>') 
    draw_orf_diagram(main_html, record, filter_list)
    main_html.write('<p></p>') 
    main_html.write('<a href="https://www.ncbi.nlm.nih.gov/nuccore/%s">Link to nucleotide sequence</a>' % (record.cluster_accession))
    draw_cds_table(main_html, record)
    if print_precursors:
        draw_orf_table(main_html, record, peptide_type)

def compress_sequence(sequence, threshold=100):
    if len(sequence) > threshold:
        return sequence[:(threshold//2)] +  "..." + sequence[-(threshold//2):]
    else:
        return sequence

def write_header(html_file):
    html_file.write("""
    <html>
    <head>
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css">
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap-theme.min.css">
    </head>
    <style media="screen" type="text/css">
     
    .square {
      width: 54px;
      height: 14px;
      background-color: white;
      outline: #ffffff solid 1px;
      text-align: center;
      line-height: 14px;
      font-size: 12px;
    }
     
    table {
       font-size: 11px;
    }
     
    </style>
     
    <script src='https://img.jgi.doe.gov//js/overlib.js'></script>
    <div class="container">
    <h1 align="center" id="header">RODEO</h1>
    """)

def get_fill_color(cds, filter_list):
    for i in range(len(cds.pfam_descr_list)):
        if cds.pfam_descr_list[i][3].upper() in filter_list:
            return "red"
        if cds.pfam_descr_list[i][0].upper() in filter_list:
            return "red"
    return "white"

def draw_CDS_arrow(main_html, cds, filter_list, sub_by, scale_factor):
    fill_color = get_fill_color(cds, filter_list)
    start = cds.start
    end = cds.end
    #HMM info
    if len(cds.pfam_descr_list) == 0:
        pfamID = "No Pfam match"
        pfam_desc = ""
    else:
        pfamID = cds.pfam_descr_list[0][0].split('.')[0] #No need for version?
        pfam_desc = cds.pfam_descr_list[0][1]
    main_html.write('<polygon points=\"')
    arrow_wid = int((start - sub_by) * scale_factor)
    arrow_wid3 = int((end - sub_by) * scale_factor)
    if abs(arrow_wid3 - arrow_wid) < 40:
        arrow_wid2 = (arrow_wid + arrow_wid3) / 2 #middle?
    else:
        if start < end:
            arrow_wid2 = arrow_wid3 - 20
        else:
            arrow_wid2 = arrow_wid3 + 20
        
    str_arrow_wid = str(arrow_wid)
    str_arrow_wid2 = str(arrow_wid2)
    str_arrow_wid3 = str(arrow_wid3)
    main_html.write(str_arrow_wid + ",10 " + str_arrow_wid + ",40 "\
                        + str_arrow_wid2 + ",40 " + str_arrow_wid2 + ",50 "+ str_arrow_wid3 \
                        + ",25 " + str_arrow_wid2 + ",0 " + str_arrow_wid2 + ",10 " + str_arrow_wid + ",10")
    main_html.write('" style="fill:' + fill_color)
    main_html.write(';stroke:black;stroke-width:.5" onMouseOver="return overlib(')
    main_html.write("'" + cds.accession_id + " - " + pfamID + " : " + pfam_desc + "'")
    main_html.write(')" onMouseOut="return nd()"/>')

def draw_orf_diagram(main_html, record, filter_list):
    main_html.write('<h3>Architecture</h3>\n')
    main_html.write('<svg width="1060" height="53">')
    bsc_start = min(record.CDSs[0].start, record.CDSs[0].end)
    bsc_end = max(record.CDSs[-1].end, record.CDSs[-1].start)
    sub_by = max(bsc_start - 500, 0)
    scale_factor = (660./(bsc_end - bsc_start))
    for cds in record.CDSs:
        draw_CDS_arrow(main_html, cds, filter_list, sub_by, scale_factor)
    main_html.write('</svg>')
    bar_length = scale_factor * 1000
    bar_legx = bar_length + 5
    main_html.write('<svg width="500" height="23">')
    main_html.write('<polygon points="')
    main_html.write("0,10 %f, 10" % (bar_length))
    main_html.write('" style="fill:white;stroke:black;stroke-width:.5" />')
    main_html.write('<text x="%f"y="12"' % (bar_legx))
    main_html.write("""font-famil="sans-serif"
                    font-size="10px"
                    text_anchor="right"
                    fill="black">1000 nucleotides</text>""")
    main_html.write('<polygon points="')
    main_html.write("0,5 0,15")
    main_html.write('" style="fill:white;stroke:black;stroke-width:.5" />')
    main_html.write('<polygon points="%f,5 %f,15' % (bar_length, bar_length))
    main_html.write('" style="fill:white;stroke:black;stroke-width:.5" />')
    main_html.write('</svg>')
    
def draw_cds_table(main_html, record):
    main_html.write("""<br><br><table class="table table-condensed">
  <tbody>
    <tr>
      <th scope="col">Accession</th>
      <th scope="col">start</th>
      <th scope="col">end</th>
      <th scope="col">direction</th>
      <th scope="col">length (aa)</th>
      <th scope="col">Pfam/HMM</th>
      <th scope="col">name</th>
      <th scope="col">description</th>
      <th scope="col">E-value</th>    
    </tr>""")
    for cds in record.CDSs:
        main_html.write("<tr>\n")
        if cds.meta:
            prot_or_nucc = "nuccore"
            acc = record.cluster_accession
            main_html.write("""\t<td><a href='https://www.ncbi.nlm.nih.gov/%s/%s?from=%s&to=%s'>%s</a></td>
            <td>%s</td> 
            <td>%d</td>
            <td>%s</td>
            <td>%d</td>""" % (prot_or_nucc, acc, min(cds.start, cds.end), max(cds.start, cds.end), cds.accession_id, cds.start, cds.end, cds.direction, round(abs(cds.end-cds.start)/3, 0)))
        else:
            prot_or_nucc = "protein"
            acc = cds.accession_id
            main_html.write("""\t<td><a href='https://www.ncbi.nlm.nih.gov/%s/%s'>%s</a></td>
            <td>%s</td> 
            <td>%d</td>
            <td>%s</td>
            <td>%d</td>""" % (prot_or_nucc, acc, cds.accession_id, cds.start, cds.end, cds.direction, round(abs(cds.end-cds.start)/3, 0)))
        if len(cds.pfam_descr_list) == 0:
            main_html.write("<td>NO PFAM MATCH</td>")
            main_html.write("<td>-</td>")
            main_html.write("<td>-</td>")
            main_html.write("<td>-</td>")
        else:
            if cds.pfam_descr_list[0][0][:2] == "PF":
                main_html.write("<td><a href='https://www.ebi.ac.uk/interpro/entry/pfam/%s'>%s</a>" % (cds.pfam_descr_list[0][0], cds.pfam_descr_list[0][0]))
            elif cds.pfam_descr_list[0][0][:4] == "TIGR":
                main_html.write("<td><a href='https://www.ebi.ac.uk/interpro/entry/ncbifam/%s'>%s</a>" % (cds.pfam_descr_list[0][0], cds.pfam_descr_list[0][0]))
            else:
                main_html.write("<td>%s" % (cds.pfam_descr_list[0][0]))
            n = 5
            for pfamid, _, _, _, in cds.pfam_descr_list[1:n]:
                if pfamid[:2] == "PF":
                    main_html.write("<br><a href='https://www.ebi.ac.uk/interpro/entry/pfam/%s'>%s</a>" % (pfamid, pfamid))
                elif pfamid[:4] == "TIGR":
                    main_html.write("<br><a href='https://www.ebi.ac.uk/interpro/entry/ncbifam/%s'>%s</a>" % (pfamid, pfamid))
                else:
                    main_html.write("<br>%s" % (pfamid))
                    
            main_html.write("</td><td>%s"%(cds.pfam_descr_list[0][3]))
            for _, _, _, name in cds.pfam_descr_list[1:n]:
                main_html.write("<br>%s"%(name))
                
            descr = cds.pfam_descr_list[0][1]
            main_html.write("</td><td>%s" % (descr))
            for _, descr, _, _, in cds.pfam_descr_list[1:n]:
                main_html.write("<br>%s" % (descr))
                
            e_val = cds.pfam_descr_list[0][2]
            main_html.write("</td><td>%.2E" % Decimal(e_val))
            for _, _, e_val, _, in cds.pfam_descr_list[1:n]:
                main_html.write("<br>%.2E" % Decimal(e_val))
            
            main_html.write("</td>")
        main_html.write("</tr>") 
    main_html.write("</tbody></table><p></p>")


def draw_orf_table(main_html, record, peptide_type):
    main_html.write("""<table class="table table-bordered">
              <tbody>
                <tr>
                <th scope="col">index</th>""")
    main_html.write("""
          <th scope="col">peptide</th>""")
    main_html.write("""
      <th scope="col">start</th>
      <th scope="col">end</th>
      <th scope="col">dir</th>""")
    if peptide_type != '':
        main_html.write("""<th scope="col">score</th>""")
    main_html.write("""\n</tr>""")
    
    index = 1
    prev_end = 0
    rowcolor = 1
    for orf in record.orfs:
        if prev_end != orf.end:
            rowcolor = rowcolor * -1
        if rowcolor == 1:
            main_html.write("<tr>\n")
        else:
            main_html.write('<tr style="background-color:#E8E8E8">\n')
        main_html.write("<td>%d</td>" % (index))
        main_html.write('<td style="text-align:right">%s</td>' % (compress_sequence(orf.sequence)))
        main_html.write("<td>%d</td>" % (orf.start))
        main_html.write("<td>%d</td>" % (orf.end))
        main_html.write("<td>%s</td>" % (orf.direction))
        if peptide_type != '':
            main_html.write("<td>%s</td>" % orf.score)
        main_html.write("</tr>\n")
        prev_end = orf.end
        index += 1
    main_html.write("</tbody></table><p></p>")

class Sub_Seq(object):
    """Useful for storing subsequences and their coordinates"""
    def __init__(self, start, end, accession_id=None, meta=False, seq="", score=0):
        self.start = start
        self.end = end
        if int(self.end) > int(self.start):
            self.direction = "+"
        else:
            self.direction = "-"
        self.accession_id = accession_id
        self.pfam_descr_list = []
        self.meta = meta
        self.sequence = seq
        self.score = score


class My_Record(object):
    """ """
    #TODO get genus and species frecom gb file
    def __init__(self, query_accession_id):
        self.query_accession_id = query_accession_id
        self.cluster_accession = ""
        self.cluster_genus_species = ""
        self.CDSs = []
        self.orfs = []
        self.locus = -1

def __main__():

# Document starts with argument parsing
    parser = argparse.ArgumentParser("RODEO HTML Postscript")
    parser.add_argument('folder_input', type=str)
    parser.add_argument('OFO', type=str)
    parser.add_argument('-s', '--single', type=int, default=0)
    parser.add_argument('-print', '--print_precursors', action='store_true', default=False)
    parser.add_argument('-pt', '--peptide_type', type=str, default="")
    parser.add_argument('-f', '--locus_filter', nargs='*', default=[])
    parser.add_argument('-e', '--exclude', nargs='*', default=[])
    args = parser.parse_args()

    st = time.time()
    folders=[]
    if args.single:
        folders.append(args.folder_input)
    else:
        folders = glob.glob(args.folder_input)

    pre_filter_list = {}
    pre_exclude_list = {}
    filter_list = []


    if len(args.locus_filter) > 0:
        for folder in folders:
            if not os.path.isdir(folder):
                continue
            curr_file = folder.strip()+"/main_co_occur.csv"
            reader = csv.reader(open(curr_file,'r'))
            lines = list(reader)
            i = 1
            meta = False
            if "Locus" in lines[0]:
                meta = True
            while i < len(lines):
                attr = lines[i]
                if meta:
                    s = 8
                    temp_name = "".join((attr[0] + "_" + attr[1] + "_" + attr[3]).split())
                else:
                    s = 7
                    temp_name = "".join((attr[0] + "_---_" + attr[2]).split())
                for keyword in args.locus_filter:
                    for t in range(s, len(attr)-2, 4):
                    #for t in range(s, min(s+6, len(attr))-2, 4):
                        if keyword in str(attr[t]) and Decimal(attr[t+3]) < 1E-5:
                            if (temp_name) not in pre_filter_list.keys():
                                pre_filter_list[temp_name] = []
                            pre_filter_list[temp_name].append(keyword)
                for keyword in args.exclude:
                    for t in range(s, len(attr)-2, 4):
                        if keyword in str(attr[t]) and Decimal(attr[t+3]) < 1E-5:
                            if (temp_name) not in pre_exclude_list.keys():
                                pre_exclude_list[temp_name] = []
                            pre_exclude_list[temp_name].append(keyword)
                i = i+1

        for hit in set(pre_filter_list.keys()).difference(set(pre_exclude_list.keys())):
            #if all(x in pre_filter_list[hit] for x in args.locus_filter):
            if any(x in pre_filter_list[hit] for x in args.locus_filter):
                filter_list.append(hit)
    records = {}
    #metabait = ["PF02624__YcaO", "PF05147__LANC_like", "PF13471__Transglut_core3", "PF05114__DUF692", "PF04738__Lant_dehydr_N", "TIGR03793__TIGR03793",  "Graspetide_synthetase__Graspetide_synthetase", "PF00590__TP_methylase", "TIGR03798__leader_Nif11", "PF05402__PqqD"] 

    for folder in folders:
        if not os.path.isdir(folder):
            continue
        curr_file = folder.strip()+"/main_co_occur.csv"
        reader = csv.reader(open(curr_file,'r'))
        lines = list(reader)
        i = 1
        meta = False
        if "Locus" in lines[0]:
            meta = True
        while i < len(lines):
            attr = lines[i]
            if meta:
                curr_query = "".join((attr[0] + "_" + attr[1] + "_" + attr[3]).split())
            else:
                curr_query = "".join((attr[0] + "_---_" + attr[2]).split())
            if len(args.locus_filter) > 0:
                if curr_query not in filter_list:
                    i = i+1
                    continue
            if curr_query not in records.keys():
                records[curr_query] = My_Record(curr_query)
                if meta:
                    records[curr_query].cluster_accession = attr[3]
                    records[curr_query].cluster_genus_species = attr[2]
                    records[curr_query].locus = attr[1]
                else:
                    records[curr_query].cluster_accession = attr[2]
                    records[curr_query].cluster_genus_species = attr[1]
            if meta:
                curr_protein = attr[4] + "_" + attr[1] + "_" + attr[3]      # change to 3 if not meta
                curr_cds = Sub_Seq(int(attr[5]), int(attr[6]), accession_id=attr[4], meta=meta)
            else:
                curr_protein = attr[3] + "_---_" + attr[2]
                curr_cds = Sub_Seq(int(attr[4]), int(attr[5]), accession_id=attr[3])
            j = lines[0].index('PfamID1')
            while j < len(attr)-2:
                pfam_descr = [attr[j], attr[j+2], attr[j+3], attr[j+1]]
                curr_cds.pfam_descr_list.append(pfam_descr)
                j = j+4
            records[curr_query].CDSs.append(curr_cds)
            i = i+1

    records_to_print = []
    main_csv = open(args.OFO+"_results.csv", 'w')
    if args.print_precursors:
        i = -1
        for folder in folders:
            if args.peptide_type:
                curr_file = folder.strip()+args.peptide_type+"/"+args.peptide_type+"_features.csv"
            else:
                curr_file = folder.strip()+"/main_results.csv"
            try:
                reader = csv.reader(open(curr_file,'r'))
            except:
                continue
            lines = list(reader)
            if i == -1:
                main_csv.write(",".join(lines[0]))
            i = 1
            if "Locus" in lines[0]:
                meta = True
            else:
                meta = False
            while i < len(lines):
                attr = lines[i]
                if meta:
                    curr_query = "".join((attr[0][:-11] + "_" + attr[1] + "_" + attr[3]).split())
                    j = 4
                else:
                    curr_query = "".join((attr[0][:-11] + "_---_" + attr[2]).split())
                    j = 3
                if curr_query in filter_list:
                    main_csv.write("\n"+",".join(attr))
                    if args.peptide_type == "boro":
                        curr_orf = Sub_Seq(int(attr[j+4]), int(attr[j+5]), meta=meta, seq=attr[j], score=attr[j+10])
                        if attr[j+11] == "Y":
                            main_csv.write("\n"+",".join(attr))
                            records_to_print.append(curr_query)
                            (records[curr_query].orfs).append(curr_orf)
                    elif args.peptide_type == "grasp":
                        curr_orf = Sub_Seq(int(attr[j+3]), int(attr[j+4]), meta=meta, seq=attr[j], score=attr[j+8])
                        if attr[j+9] == "Y":
                            main_csv.write("\n"+",".join(attr))
                            records_to_print.append(curr_query)
                            (records[curr_query].orfs).append(curr_orf)
                    elif args.peptide_type:
                        curr_orf = Sub_Seq(int(attr[j+2]), int(attr[j+3]), meta=meta, seq=attr[j]+attr[j+1], score=attr[j+4])
                        if attr[j+5] == "Y":
                            main_csv.write("\n"+",".join(attr))
                            records_to_print.append(curr_query)
                            (records[curr_query].orfs).append(curr_orf)
                    else:
                        curr_orf = Sub_Seq(int(attr[j]), int(attr[j+1]), meta=meta, seq=attr[j+3])
                        if re.search(r'Y..P', curr_orf.sequence):
                            main_csv.write("\n"+",".join(attr))
                            records_to_print.append(curr_query)
                        (records[curr_query].orfs).append(curr_orf)
                i = i+1

    else:
        records_to_print = records.keys()

    main_html = open(args.OFO+".html", 'w')
    write_header(main_html)
    for record in set(records_to_print):
        write_record(main_html, records[record], args.locus_filter, args.print_precursors, curr_query=record, peptide_type=args.peptide_type)

    end = time.time()
    print(str(end-st))

if __name__=="__main__":
    __main__()

