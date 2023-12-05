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

import csv, subprocess, os, statistics
import logging
import hmmer_utils
from rodeo_main import VERBOSITY

logger = logging.getLogger(__name__)
logger.setLevel(VERBOSITY)
# create console handler and set level
ch = logging.StreamHandler()
ch.setLevel(VERBOSITY)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)

def run_prodigal(record, whole_contig=False):
    try:
        #print(record.query_short)
        if whole_contig:
            prod_prefix = "/tmp/" + record.query_short + "_" + record.random_tag
            prod_os_file = open(prod_prefix+"os.txt", 'w')
            prod_use_file = open(prod_prefix+"prod.fasta", 'w')
            prod_use_file.write(">%s %s %s\n%s" % (record.query_accession_id + "_" + record.random_tag, record.genus, record.species, 
                            record.cluster_sequence))
        else:
            prod_prefix = "/tmp/" + record.query_short + "_" + record.random_tag
            prod_os_file = open(prod_prefix+"os.txt", 'w')
            prod_use_file = open(prod_prefix+"prod.fasta", 'w')
            prod_use_file.write(">%s %s %s\n%s" % (record.query_accession_id  + "_" + record.random_tag, record.genus, record.species, 
        					record.cluster_sequence[record.window_start:record.window_end]))
        prod_use_file.close()
        if record.prod_window_end-record.prod_window_start < 100000 or whole_contig:
            process = subprocess.Popen(["prodigal", "-i", prod_prefix+"prod.fasta", "-o", prod_prefix+"output.txt", 
            							"-p", "meta", "-s", prod_prefix+"orfs.tsv", "-q"], stdout=prod_os_file, stderr=prod_os_file)
            process.wait()
        else:
            prod_train_file = open(prod_prefix+"train.fasta", 'w')
            prod_train_file.write(">%s %s %s\n%s" % (record.query_accession_id  + "_" + record.random_tag, record.genus, record.species, 
            						record.cluster_sequence[record.prod_window_start:record.prod_window_end]))
            prod_train_file.close()
            process = subprocess.Popen(["prodigal", "-i", prod_prefix+"train.fasta", "-o", prod_prefix+"output.txt", 
            							"-t", prod_prefix+"train.txt", "-q"], stdout=prod_os_file, stderr=prod_os_file)
            process.wait()
            process = subprocess.Popen(["prodigal", "-i", prod_prefix+"prod.fasta", "-o", prod_prefix+"output.txt", 
            							"-t", prod_prefix+"train.txt", "-s", prod_prefix+"orfs.tsv", "-q"], stdout=prod_os_file, stderr=prod_os_file)
            process.wait()
             
        try:
            os.remove(prod_prefix+"train.fasta")
        except:
            pass
        try:
            os.remove(prod_prefix+"output.txt")
        except:
            pass
        try:
            os.remove(prod_prefix+"prod.fasta")
        except:
            pass
        try:
            os.remove(prod_prefix+"train.txt")
        except:
            pass
        try:
        	os.remove(prod_prefix+"os.txt")
        except:
        	pass
    except KeyboardInterrupt:
        try:
            os.remove(prod_prefix+"train.fasta")
        except:
            pass
        try:
            os.remove(prod_prefix+"output.txt")
        except:
            pass
        try:
            os.remove(prod_prefix+"prod.fasta")
        except:
            pass
        try:
            os.remove(prod_prefix+"train.txt")
        except:
            pass
        try:
            os.remove(prod_prefix+"orfs.tsv")
        except:
            pass
        try:
        	os.remove(prod_prefix+"os.txt")
        except:
        	pass
        logger.critical("SIGINT recieved during Prodigal")
        raise KeyboardInterrupt
