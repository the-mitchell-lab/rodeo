#==============================================================================
# Copyright (C) 2017 Bryce L. Kille
# University of Illinois
# Department of Chemistry
#
# Copyright (C) 2017 Christopher J. Schwalen
# University of Illinois
# Department of Chemistry
#
# Copyright (C) 2017 Douglas A. Mitchell
# University of Illinois
# Department of Chemistry
#
# License: GNU Affero General Public License v3 or later
# Complete license availabel in the accompanying LICENSE.txt.
# or <http://www.gnu.org/licenses/>.
#
# This file is part of RODEO2.
#
# RODEO2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RODEO2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#==============================================================================
import datetime
from decimal import Decimal
from rodeo_main import VERSION

index = 0

def compress_sequence(sequence, threshold=100):
    if len(sequence) > threshold:
        return sequence[:(threshold//2)] +  "..." + sequence[-(threshold//2):]
    else:
        return sequence

def write_header(html_file, master_conf, peptide_type):
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
    <h1 align="center" id="header">RODEO2</h1>
    <div class="row">
         <div class="col-md-5">
            <h3>Parameters</h3>
                 <table class="table table-condensed" style="width:100%;">
                        <tr><th scope="row">Run Time</th><td>""" + datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y") + """</td></tr>
                        """)
    html_file.write('<tr><th scope="row">Version</th><td>' + VERSION +  "</td></tr>\n")
    html_file.write('<tr><th scope="row">Gene Window</th><td>+/-' + str(master_conf['general']['variables']['fetch_n'])+ "  " + master_conf['general']['variables']['fetch_type'] + "</td></tr>\n")
    html_file.write("""
                        <tr><th scope="row">Peptide Range</th><td>""" + str(master_conf[peptide_type]['variables']['precursor_min']) + "-" + str(master_conf[peptide_type]['variables']['precursor_max']) + """ aa</td></tr>
                        <tr><th scope="row">Fetch Distance</th><td>""" + str(master_conf['general']['variables']['overlap']) + """bp</td></tr>
                        <tr><th scope="row">Peptide Type</th><td> %s </td></tr>
                  </table>
         </div>
         <div class="col-md-6">
                    <div class="panel panel-default">""" % (peptide_type))
    write_legend(html_file, master_conf, peptide_type)

def write_legend(main_html, conf, peptide_type):
    main_html.write("""<div class="panel-heading">Annotation Legend</div>
                <table class="table table-striped">
                        <tr>
                            <th>Appearance</th>
                            <th>Accession/Name</th>
                        </tr>""")
    for pfam in conf[peptide_type]['pfam_colors'].keys():
        main_html.write("<td><div class='square' style='outline-color:black;background-color:{}; color:grey;'>{}</div></td>\n".format(conf[peptide_type]['pfam_colors'][pfam], ""))
        main_html.write("<td>{}</td>\n".format(pfam))
        main_html.write('</tr>\n')
    if peptide_type != 'general':
        main_html.write("<td><div class='square' style='outline-color:black;background-color:{}; color:grey;'>{}</div></td>\n".format("black", ""))
        main_html.write("<td>{}</td>\n".format("predicted " + peptide_type +" peptide"))
        main_html.write('</tr>\n')
    main_html.write("""</table>
                        </div></div>
                        </div>""")

def get_fill_color(cds, peptide_conf):
    if "RRE" in peptide_conf["pfam_colors"].keys() and cds.isRRE:
        return peptide_conf["pfam_colors"]["RRE"]
    for i in range(len(cds.pfam_descr_list)):
        if cds.pfam_descr_list[i][3].upper() in peptide_conf['pfam_colors'].keys():
            return peptide_conf['pfam_colors'][cds.pfam_descr_list[i][3].upper()]
        if cds.pfam_descr_list[i][0].upper() in peptide_conf['pfam_colors'].keys():
            return peptide_conf['pfam_colors'][cds.pfam_descr_list[i][0].upper()]
    return "white"

def draw_CDS_arrow(main_html, cds, peptide_conf, sub_by, scale_factor):
    fill_color = get_fill_color(cds, peptide_conf)
    start = cds.start
    end = cds.end
    #HMM info
    if len(cds.pfam_descr_list) == 0:
        if cds.inferred:
            pfamID = "INFERRED GENE"
        else:
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

def draw_orf_arrow(main_html, orf, sub_by, scale_factor, index):
    fill_color = "rgba(0, 0, 0, {}) ".format((min(orf.confidence**2,1)))
    start = orf.start
    end = orf.end
    arrow_wid = int((start - sub_by) * scale_factor)
    arrow_wid3 = int((end - sub_by) * scale_factor)
    if abs(arrow_wid3 - arrow_wid) < 40:
        arrow_wid2 = (arrow_wid + arrow_wid3) / 2 #middle?
    else:
        if start < end:
            arrow_wid2 = arrow_wid3 - 20
        else:
            arrow_wid2 = arrow_wid3 + 20
    letter_x = (arrow_wid + arrow_wid3) / 2
    main_html.write('<polygon points="')
    main_html.write("%d,10 %d,40 %d,40 %d,50 %d,25 %d,0 %d,10 %d,10" % (arrow_wid, arrow_wid, arrow_wid2, arrow_wid2, arrow_wid3, arrow_wid2, arrow_wid2, arrow_wid))
    main_html.write('" style="fill:' + fill_color + ';stroke:black;stroke-width:.5" />' )
    main_html.write('<text x="%d"y="32" font-family="sans-serif" font-size="12px" text-anchor="middle" fill="grey">' % (letter_x))
    main_html.write(str(index))
    main_html.write("</text>")

def draw_orf_diagram(main_html, peptide_conf, record, peptide_type):
    main_html.write('<h3>Architecture</h3>\n')
    main_html.write('<svg width="1060" height="53">')
    bsc_start = min(record.CDSs[0].start, record.CDSs[0].end)
    bsc_end = max(record.CDSs[-1].end, record.CDSs[-1].start)
    sub_by = bsc_start - 500
    scale_factor = (660./(bsc_end - bsc_start))
    for cds in record.CDSs:
        draw_CDS_arrow(main_html, cds, peptide_conf, sub_by, scale_factor)
    main_html.write('</svg>')
    main_html.write('<svg width="1060" height="53">')
    index = 0
    if peptide_type != 'general':
        for ripp in record.ripps[peptide_type]:
            if ripp.score <= ripp.CUTOFF // 2:
                continue
            index += 1
            if peptide_conf['variables']['precursor_min'] <= len(ripp.sequence) <= peptide_conf['variables']['precursor_max'] or \
                            ("M" in ripp.sequence[-peptide_conf['variables']['precursor_max']:]):
                            draw_orf_arrow(main_html, ripp, sub_by, scale_factor, index)
    
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
        if cds.inferred:
            prot_or_nucc = "nuccore"
            acc = record.cluster_accession
        else:
            prot_or_nucc = "protein"
            acc = cds.accession_id
        main_html.write("""\t<td><a href='https://www.ncbi.nlm.nih.gov/%s/%s'>%s</a></td>
            <td>%s</td> 
            <td>%d</td>
            <td>%s</td>
            <td>%d</td>""" % (prot_or_nucc, acc, cds.accession_id, cds.start, cds.end, cds.direction, len(cds.sequence)))
        if len(cds.pfam_descr_list) == 0:
            if cds.inferred:
                main_html.write("<td>INFERRED GENE</td>")
            else:
                main_html.write("<td>NO PFAM MATCH</td>")
            main_html.write("<td>-</td>")
            main_html.write("<td>-</td>")
            main_html.write("<td>-</td>")
        else:
            if cds.pfam_descr_list[0][0][:2] == "PF":
                main_html.write("<td><a href='https://www.ebi.ac.uk/interpro/entry/pfam/%s'>%s</a>" % (cds.pfam_descr_list[0][0], cds.pfam_descr_list[0][0]))
            elif cds.pfam_descr_list[0][0][:4] == "TIGR":
                main_html.write("<td><a href='https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/%s'>%s</a>" % (cds.pfam_descr_list[0][0], cds.pfam_descr_list[0][0]))
            else:
                main_html.write("<td>%s" % (cds.pfam_descr_list[0][0]))
#            main_html.write("-%s"%(cds.pfam_descr_list[0][3]))
            n = 5
            for pfamid, _, _, _, in cds.pfam_descr_list[1:n]:
                if pfamid[:2] == "PF":
                    main_html.write("<br><a href='https://www.ebi.ac.uk/interpro/entry/pfam/%s'>%s</a>" % (pfamid, pfamid))
                elif pfamid[:4] == "TIGR":
                    main_html.write("<br><a href='https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/%s'>%s</a>" % (pfamid, pfamid))
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
       

def draw_orf_table(main_html, record, peptide_type, master_conf):
    print_precursors = master_conf[peptide_type]['variables']['print_precursors']
    if not print_precursors:
        return
    main_html.write("""<table class="table table-bordered">
              <tbody>
                <tr>
                <th scope="col">index</th>""")
    if peptide_type in ["lasso", "lanthi", "lanthi1", "lanthi2", "lanthi3", "lanthi4", "sacti", "thio", "linar"]:
        main_html.write("""
              <th scope="col">leader</th>
              <th scope="col">core</th>""")
    elif peptide_type == "grasp":
        main_html.write("""
              <th scope="col">first half</th>
              <th scope="col">second half</th>""")
    else:
        main_html.write("""
          <th scope="col">peptide</th>""")
    main_html.write("""
      <th scope="col">start</th>
      <th scope="col">end</th>
      <th scope="col">dir</th>""")
    if peptide_type != 'general':
      main_html.write("""<th scope="col">score</th>""")
    main_html.write("""\n</tr>""")
    
    index = 1
    if peptide_type == 'general':
        prev_end = 0
        rowcolor = 1
        for orf in record.intergenic_orfs:
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
            main_html.write("</tr>\n")
            prev_end = orf.end
            index += 1
            
    elif peptide_type in ["lasso", "lanthi", "lanthi1", "lanthi2", "lanthi3", "lanthi4", "sacti", "thio", "grasp", "linar", "boro"]:
        prev_end = 0
        rowcolor = 1

        for ripp in record.ripps[peptide_type]:
            if ripp.score <= ripp.CUTOFF // 2:
                continue
            if not master_conf[peptide_type]['variables']['precursor_min'] <= len(ripp.sequence) <= master_conf[peptide_type]['variables']['precursor_max'] and \
                            not ("M" in ripp.sequence[-master_conf[peptide_type]['variables']['precursor_max']:]):
                            continue
            if prev_end != ripp.end:
                rowcolor = rowcolor * -1
            if rowcolor == 1:
                main_html.write("<tr>\n")
            else:
                main_html.write('<tr style="background-color:#E8E8E8">\n')
            main_html.write("<td>%d</td>" % (index))
            if print_precursors:

                if peptide_type in ["lasso", "lanthi", "lanthi1", "lanthi2", "lanthi3", "lanthi4", "sacti", "thio", "grasp", "linar"]:
                    main_html.write('<td style="text-align:right">%s</td>' % (compress_sequence(ripp.leader,50)))
                    main_html.write('<td style="text-align:left">%s</td>' % (compress_sequence(ripp.core, 50)))
                else:
                    main_html.write('<td style="text-align:right">%s</td>' % (compress_sequence(ripp.sequence)))
            main_html.write("<td>%d</td>" % (ripp.start))
            main_html.write("<td>%d</td>" % (ripp.end))
            main_html.write("<td>%s</td>" % (ripp.direction))
            main_html.write("<td>%d</td>" % (ripp.score))
            main_html.write("</tr>\n")
            prev_end = ripp.end
            index += 1
    main_html.write("</tbody></table>")
    
        

def write_table_of_contents(main_html, queries, meta=False):
    global index
    index = 0
    main_html.write("<h3> Input Queries (click to navigate)</h3>")
    main_html.write('<ul style="list-style-type:none">')
    if meta == False:
        for query in queries:
            main_html.write('<li><a href="#%s">%s</a></li>' % (query, query))
    main_html.write("</ul>\n")

def write_failed_query(main_html, query, message):
    global index
    main_html.write('<h2 id="%s"> Results for %s\n' % (query, query))
    main_html.write('<a href="#header"><small><small>back to top</small></small></a></h2>') #TODO keep for single?
    main_html.write('<h2 id="num{}"><a href="#num{}"><small><small>previous</small></small></a><small><small>\t\t\t\t\t-\t\t\t\t\t</small></small><a href="#num{}"><small><small>next</small></small></a></h2>'.format(index, index-1, index+1))
    main_html.write('<p></p>') # TODO why
    main_html.write(message)
    main_html.write('<p></p>')
    index += 1
    
def write_record(main_html, master_conf, record, peptide_type):
    #RESULTS FOR xxxx
    #DRAW ORF
    #DRAW ORF scale
    #PUT LINK TO nuc SEQUENCE
    #TABLE of CDS
    #TABLE of ORFs
    global index
    if record.bait_iteration != -1:
        main_html.write('<h2 id="%s"> Results for %s - locus %s for gene at %s in %s [%s]' % 
                                   (record.query_accession_id, record.query_accession_id, str(record.bait_iteration + 1), 
                                    record.cds_start_list[record.bait_iteration], record.cluster_accession, record.cluster_genus_species))
    else:
        main_html.write('<h2 id="%s"> Results for %s [%s]\n' % (record.query_accession_id, record.query_accession_id, record.cluster_genus_species))
    main_html.write('<a href="#header"><small><small>back to top</small></small></a></h2>') #TODO keep for single?
    main_html.write('<h2 id="num{}"><a href="#num{}"><small><small>previous</small></small></a><small><small>\t\t\t\t\t-\t\t\t\t\t</small></small><a href="#num{}"><small><small>next</small></small></a></h2>'.format(index, index-1, index+1))
    main_html.write('<p></p>') # TODO why
    if any(record.CDSs):
        draw_orf_diagram(main_html, master_conf[peptide_type], record, peptide_type)
    main_html.write('<p></p>') 
    main_html.write('<a href="https://www.ncbi.nlm.nih.gov/nuccore/%s">Link to nucleotide sequence</a>' % (record.cluster_accession))
    draw_cds_table(main_html, record)
    draw_orf_table(main_html, record, peptide_type, master_conf)
    index += 1


