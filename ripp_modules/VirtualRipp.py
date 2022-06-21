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

import csv
import os
import subprocess
import importlib
import logging
import random
import string
import socket
import ripp_modules.SvmClassify as svmc
from rodeo_main import VERBOSITY
import pathlib
FILE_DIR = pathlib.Path(__file__).parent.absolute()

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
index = 0

def execute(commands, input=None):
    "Execute commands in a system-independent manner"

    if input is not None:
        stdin_redir = subprocess.PIPE
    else:
        stdin_redir = None

    try:
        proc = subprocess.Popen(commands, stdin=stdin_redir,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, 
                                shell=True)
        out, err = proc.communicate(input=input)
        retcode = proc.returncode
        return out, err, retcode
    except OSError as e:
        logger.error(e)
        raise e


def ripp_write_rows(output_dir, peptide_type, accession_id, genus_species, list_of_rows, feature_count=5):
    dir_prefix = output_dir + '/{}/'.format(peptide_type)
    global index
    features_csv_file = open(dir_prefix + "temp_features.csv", 'a')
    svm_csv_file = open(os.path.join(FILE_DIR, "{}/svm/fitting_set.csv".format(peptide_type)), 'a')
    features_writer = csv.writer(features_csv_file)
    svm_writer = csv.writer(svm_csv_file)
    for row in list_of_rows:
        features_writer.writerow([accession_id, genus_species] + row[0:feature_count] + ["valid_precursor_placeholder", index, ''] + row[feature_count:])
        svm_writer.writerow([index, ''] + row[feature_count:]) #Don't include accession_id, leader, core sequence, start, end, or score
        index += 1
        
def run_svm(output_dir, peptide_type, cutoff, feature_count=5):
    runner = svmc.SVMRunner(peptide_type)
    runner.run_svm()
    svm_output_reader = csv.reader(open(os.path.join(FILE_DIR, "{}/svm/fitting_results.csv".format(peptide_type))))
    final_output_writer = csv.writer(open(output_dir + "/{}/{}_features.csv".format(peptide_type, peptide_type), 'w'))
    features_reader = csv.reader(open(output_dir + "/{}/temp_features.csv".format(peptide_type)))
    header_row = next(features_reader) #skip header
    final_output_writer.writerow(header_row)
    for row, svm_output_line in zip(features_reader, svm_output_reader):
        svm_output = svm_output_line[1]
        row[feature_count+4] = svm_output
        if int(svm_output) == 1:
            row[feature_count+1] = int(row[feature_count+1]) + 10
        if int(row[feature_count+1]) >= cutoff: #CUTOFF
            row[feature_count+2] = 'Y'
        else:
            row[feature_count+2] = 'N'
        final_output_writer.writerow(row)        

def get_match_score(a, b):
    if a == b:
        return 1
    elif (a in ['D', 'E'] and b in ['D', 'E']) \
    or (a in ['S', 'T'] and b in ['S', 'T']):
        return .5
    else:
        return 0

def get_repeat_score(repeats):
    score = 0
    for i in range(len(repeats[0])):
        char_to_match = repeats[0][i]
        if char_to_match not in ['D', 'E', 'T', 'S', 'K']:
            continue
        count_same = 0
        for r in repeats:
            count_same += get_match_score(char_to_match, r[i])
        if count_same == len(repeats):
            score += 1
    return score
            
    
def parse_radar_output(radar_output):
    score = 0
    i = 0
    repeat_count = 0
    while i < len(radar_output):
        line = radar_output[i]
        if len(line.split('|')) > 1:
            if line.split('|')[0] == "No. of Repeats":
                i += 1
                line = radar_output[i]
                repeat_count = int(line.split('|')[0])
                i += 2
                repeats = []
                for repeat_idx in range(repeat_count):
                    line = radar_output[i]
                    repeats.append(line.split('\t')[-1])
                    i += 1
                if len(repeats) > 0:
                    score = max(score, get_repeat_score(repeats))
        i += 1
    return score
                    
def get_radar_score(sequence):
    #TODO change to temp file
        pid = str(os.getpid())
        try:
            with open("/tmp/" + pid + "RADAR.fasta", 'w+') as tfile:
                tfile.write(">query\n%s" % (sequence))
                command = ["radar.py -a /tmp/" + pid + "RADAR.fasta"]
            try:
                out, err, retcode = execute(command)
            except OSError:
                logger.error("Could not run RADAR")
                try:
                    os.remove(pid+"RADAR.fasta")
                except OSError:
                    pass
                return 0
            if retcode != 0:
                logger.error('RADAR returned %d: %r', retcode,
                                err)
                logger.error(sequence)
                return 0
            try:
                os.remove("/tmp/" + pid + "RADAR.fasta")
            except OSError:
                    pass
        except KeyboardInterrupt:
            try:
                os.remove("/tmp/" + pid + "RADAR.fasta")
                return
            except OSError:
                pass
        # RADAR output is inconsistent. Most of the time it will print
        # two lines of "no results", but other times just 1...
        lines = out.decode("utf-8", errors="replace").split('\n')
        if (len(lines) < 4): 
            return 0
        results = lines[3].split('|')
        if len(results) > 1:
            return parse_radar_output(lines[1:])
        else:
            return 0

        
class VirtualRipp(object):
    def __init__(self,
                 start,
                 end,
                 sequence,
                 upstream_sequence,
                 pfam_2_coords,
                 has_rre):
        self.start = start
        self.end = end
        self.sequence = sequence
        self.upstream_sequence = upstream_sequence
        self.has_rre = has_rre
        #This line is here to ensure that the type is always set.
        #No peptide type should ever be left virtual. 
        self.peptide_type = 'virtual' 
        self.score  = 0
        self.pfam_2_coords = pfam_2_coords
        if start < end:
            self.direction = '+'
        else:
            self.direction = '-'
        self.valid_split = True #set to false if no valid split

#        self.set_leader_core()
#        self.set_monoisotopic_mass()
#        self.csv_columns = [self.leader, self.core, self.start, self.end]

    def run_svm(self, output_dir):
        try:
            rmod = importlib.import_module(self.peptide_type)
            from rmod.svm import svm_classify as svm
            from rmod import CUTOFF as cutoff
        except:
            logger.error("{} not a valid peptide type".format(self.peptide_type))
            
        svm.run_svm()
        svm_output_reader = csv.reader(open(os.path.join(FILE_DIR, self.peptide_type, "svm/fitting_results.csv")))
        final_output_writer = csv.writer(open(output_dir + "/" + self.peptide_type + '/'\
                                              + self.peptide_type + "_features.csv", 'w'))
        features_reader = csv.reader(open(output_dir + "/" + self.peptide_type + '/'\
                                              + self.peptide_type + "/temp_features.csv"))
        header_row = features_reader.next() #skip header
        final_output_writer.writerow(header_row)
        
        #need to make ubiquitous 
        for row in features_reader:
            svm_output = svm_output_reader.next()[1]
            row.append(svm_output)
            if int(svm_output) == 1:
                row[-2] = int(row[-2]) + 10
            if int(row[-2]) > cutoff: #CUTOFF
                row[7] = 'Y'
            else:
                row[7] = 'N'
            final_output_writer.writerow(row)
        
    def run_fimo_simple(self, query_motif_file=None):
    #TODO change to temp file
        if not query_motif_file:
            query_motif_file = os.path.join(FILE_DIR, self.peptide_type, self.peptide_type + "_fimo.txt")
        pid = str(os.getpid())
        try:
            with open("/tmp/" + pid + "FIMO.seq", 'w+') as tfile:
                tfile.write(">query\n%s" % (self.sequence))
            if WEB_TOOL:
                command = ["/home/ubuntu/meme/bin/fimo --text --verbosity 1 " + query_motif_file + ' ' + "/tmp/" + pid + "FIMO.seq"]
            else:
                command = ["fimo --text --verbosity 1 " + query_motif_file + ' ' + "/tmp/" + pid + "FIMO.seq"]
            try:
                out, err, retcode = execute(command)
            except OSError:
                logger.error("Could not run FIMO on %s" % (self.peptide_type))
                try:
                    os.remove(pid+"TEMP.seq")
                except OSError:
                    pass
                return ""
            if retcode != 0:
                logger.error('FIMO returned %d: %r while searching %r', retcode,
                                err, query_motif_file)
                return ""
            try:
                os.remove("/tmp/" + pid + "FIMO.seq")
            except OSError:
                    pass
        except KeyboardInterrupt:
            try:
                os.remove("/tmp/" + pid + "FIMO.seq")
                return
            except OSError:
                pass
        return out.decode("utf-8")

        
    def get_min_dist(self, coords_list):
        if coords_list == []:
            return 66666
        min_dist = abs(self.start-coords_list[0][0])
        for coord in coords_list:
            min_dist = min(abs(self.start-coord[0]), abs(self.end-coord[0]),
                           abs(self.start-coord[1]), abs(self.end-coord[1]),
                           min_dist)
        return min_dist
