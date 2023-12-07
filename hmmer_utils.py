# -*- coding: utf-8 -*-
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

import os
import csv
import subprocess
import logging 
import socket
from rodeo_main import VERBOSITY

WEB_TOOL = False
if socket.gethostname() == "rodeo.scs.illinois.edu":
    WEB_TOOL = True
    

logger = logging.getLogger(__name__)
logger.setLevel(VERBOSITY)
# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(VERBOSITY)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)
pid_prefix = None    

def _generate_fasta_from_string(query_string):
    handle = open(pid_prefix+"fasta_file.tmp.fasta", "w+")
    handle.write(">temp_string\n" + query_string)

def _generate_hmmer(hmm_file):
    hmmer_process = subprocess.call(["hmmscan", "-o", pid_prefix+"hmm_out.tmp.tab", "--cpu=1", "--noali",
                                      "--domtblout", pid_prefix+"pFamInfo.tmp.tab", hmm_file,
                                      pid_prefix+"fasta_file.tmp.fasta"])    

#TODO try/except blocks for parsing to double check format
def _parse(n, e_cutoff):
    ret_list = []
    NUM_COLUMNS = 23
    pfam_accession = ''
    with open(pid_prefix+'pFamInfo.tmp.tab') as handle:
        for line in csv.reader(handle, delimiter=' '):
            if line[0][0] == '#':
                continue
            tmp_list = []
            for item in line:
                if item != '':
                    tmp_list.append(item)
            if tmp_list[1] != pfam_accession: #Check to make sure its not the same as the last one
                pfam_accession = tmp_list[1].split('.')[0]
#                if pfam_accession == '-':
#                    pfam_accession = "N/A"
                pfam_name = tmp_list[0]
                description = ' '.join(tmp_list[NUM_COLUMNS-1:])
                e_val = tmp_list[6]
                if float(e_val) < e_cutoff and len(ret_list) < n: #Will go through one extra line but oh well #TODO
                    found = False
                    for i in range(len(ret_list)):
                        (pfam_accession2, description2, e_val2, name2) = ret_list[i]
                        if pfam_accession == pfam_accession2 and pfam_name == name2:
                            found = True
                            if float(e_val) < float(e_val2):
                                ret_list[i] = (pfam_accession, description2, float(e_val2), pfam_name)
                            break
                    if not found:
                        ret_list.append((pfam_accession, description, float(e_val), pfam_name))
                else:
                    break
    ret_list_final = []
    for acc, desc, e_val, name in ret_list:
        in_list = False
        for final_acc, _, _, final_name in ret_list_final:
            if acc == final_acc and name == final_name:
                in_list = True
        if not in_list:
            ret_list_final.append((acc, desc, e_val, name))
    return ret_list_final

def get_hmmer_info(query, primary_hmm, cust_hmm, output_dir, n=5, e_cutoff=.001, query_is_accession=False): #TODO handle lists of accessions
    """Returns top n hmmscan hits with e_values lower than e_cutoff"""
    try:
        global pid_prefix
        pid_prefix = output_dir + "/tmp/" + str(os.getpid())
        _generate_fasta_from_string(query)
        pfam_desc_list = []
        try:
            if primary_hmm:
                _generate_hmmer(primary_hmm)
                pfam_desc_list = _parse(n, e_cutoff) 
        except OSError:
            logger.error("Couldn't find %s" % (primary_hmm)) 
        
        for hmm in cust_hmm:
            if WEB_TOOL:
                try:
                    for f in os.listdir(hmm):
                        if f[-4:].lower() != '.hmm' and f[-4:].lower() != '.lib':
                            continue
                        if f == 'Pfam-A.hmm':
                            continue
                        _generate_hmmer(hmm + f)
                        pfam_desc_list += _parse(n, e_cutoff)
                except OSError:
                    logger.error("Couldn't find %s" % (hmm))
            else:
                try:
                    _generate_hmmer(hmm)
                    add_list = _parse(n, e_cutoff)
                    pfam_desc_list += add_list
                except OSError:
                    logger.error("Couldn't find %s" % (hmm))
        
        if len(pfam_desc_list) > 0:
            pfam_desc_list = sorted(pfam_desc_list, key=lambda entry: entry[2])
            pfam_desc_list = sorted(pfam_desc_list, key=lambda entry: entry[0])
            pfam_temp_list = [pfam_desc_list[0]]
            (curr_pfam, curr_desc, curr_eval, curr_name) = pfam_desc_list[0]
            for i in range(len(pfam_desc_list)):
                if pfam_desc_list[i][0] == pfam_temp_list[-1][0] and pfam_desc_list[i][3] == pfam_temp_list[-1][3]:
                    continue
                pfam_temp_list.append(pfam_desc_list[i])
                
            pfam_desc_list = sorted(pfam_temp_list, key=lambda entry: entry[2])
        # os.remove(pid_prefix+'pFamInfo.tmp.tab')
        # os.remove(pid_prefix+"fasta_file.tmp.fasta")
        # os.remove(pid_prefix+"hmm_out.tmp.tab")
        return pfam_desc_list
    except KeyboardInterrupt:
        try:
            os.remove(pid_prefix+'pFamInfo.tmp.tab')
        except:
            pass
        try:
            os.remove(pid_prefix+"fasta_file.tmp.fasta")
        except:
            pass
        try:
            os.remove(pid_prefix+"hmm_out.tmp.tab")
        except:
            pass
        logger.critical("SIGINT recieved during HMMScan")
        raise KeyboardInterrupt
