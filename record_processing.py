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

from entrez_utils import get_gb_handles, get_record_from_gb_handle
import logging
from rodeo_main import VERBOSITY, QUEUE_CAP
import traceback
import sys

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

class ErrorReport(object):
    
    def __init__(self, query, error_message):
        self.query = query
        self.query_accession_id = query
        self.error_message = error_message
        
def process_record_worker(unprocessed_records_q, processed_records_q, args, master_conf, ripp_modules, index):
    try:
#        my_id = str(os.getpid())
        my_id = str(index + 1)
        logger.debug("Worker process %s started" % (my_id))
        record = unprocessed_records_q.get()
        while record != QUEUE_CAP:
            if type(record) == ErrorReport:
                processed_records_q.put(record)
                record = unprocessed_records_q.get()
                continue
            try:
                logger.info("Worker process %s is processing %s" % (my_id, record.query_accession_id))
                if master_conf['general']['variables']['fetch_type'].lower() == 'cds':
                    record.trim_to_n_orfs(master_conf['general']['variables']['fetch_n'], master_conf['general']['variables']['fetch_distance'])
                elif master_conf['general']['variables']['fetch_type'].lower() == 'nucs':
                    record.trim_to_n_nucleotides(master_conf['general']['variables']['fetch_n'])
                if "grasp" in args.peptide_types:
                    record.run_radar()
                record.annotate_w_RREFinder()
                record.annotate_w_hmmer(master_conf['general']['variables']['pfam_dir'], args.custom_hmm, 
                                        min_length=master_conf['general']['variables']['precursor_min'], 
                                        max_length=master_conf['general']['variables']['precursor_max'])
                record.set_intergenic_seqs(min_length=master_conf['general']['variables']['precursor_min'], 
                                           max_length=master_conf['general']['variables']['precursor_max'])
                record.set_intergenic_orfs(min_aa_seq_length=master_conf['general']['variables']['precursor_min'], 
                                           max_aa_seq_length=master_conf['general']['variables']['precursor_max'],
                                           overlap=master_conf['general']['variables']['overlap']) 
                for peptide_type in args.peptide_types:
                    module = ripp_modules[peptide_type]
                    record.filter_RREs_and_HMMs(hmm_list=list(master_conf[peptide_type]["pfam_colors"].keys()))
                    record.set_ripps(module, master_conf)
                    record.score_ripps(module, master_conf['general']['variables']['pfam_dir'], args.custom_hmm)
                    record.color_ripps(module)
                logger.debug("Worker process %s finished processing %s" % (my_id, record.query_accession_id))
                processed_records_q.put(record)
            except KeyboardInterrupt:
                raise KeyboardInterrupt
            except Exception as e:
                logger.error("ERROR FOR %s" % (record.query_accession_id))
                logger.error(e)
                traceback.print_exc(file=sys.stdout)
                processed_records_q.put(ErrorReport(record.query_accession_id, str(e)))
                logger.error("Worker process %s is moving on" % (my_id))

            record = unprocessed_records_q.get()
        
        logger.debug("Worker process %s pulled queue cap" % (my_id))
        unprocessed_records_q.put(QUEUE_CAP) #Replace cap for other threads
        processed_records_q.put(QUEUE_CAP)
        return
    except KeyboardInterrupt:
#        pid = str(os.getpid())
#        for f in glob.glob("tmp_files/" + pid + "*"):
#            os.remove(f)
        logger.critical("KeyboardInterrupt recieved during record processing")
        return
    
    
def fill_request_queue(queries, processed_records_q, unprocessed_records_q, args, master_conf, ripp_modules):
    try:
        for query in queries:
            logger.debug("Fetching %s data" % query)
            if '.gbk' != query[-4:] and '.gb' != query[-3:]: #accession_id
                gb_handles = get_gb_handles(query)
                nuccore_accession = query
                if type(gb_handles) is int:
                    if gb_handles == -1:
                        error_message = "No results in protein db for Esearch on %s" % (query)
                    elif gb_handles == -2:
                        error_message = "No results in nuccore db for value obtained from protein db"
                    elif gb_handles == -3:
                        error_message = "Any response failure from Entrez database (error on database side)"
                    else:
                        error_message = "Unknown Entrez error."
                    unprocessed_records_q.put(ErrorReport(query, error_message))
                    continue
            else:#gbk file
                nuccore_accession = query.split('\t')[0]
                try:
                    gb_handles = [open(query.split('\t')[1])]
                except OSError as e:
                    error_message = "Error opening %s" % (query)
                    logger.error(e)
                    unprocessed_records_q.put(ErrorReport(query, error_message))
                    continue
            for handle in gb_handles:
                record = get_record_from_gb_handle(handle, nuccore_accession)
                if type(record) is int:
                    if record == -1:
                        error_message = "Couldn't process %s Genbank filestream. May be corrupt."\
                          % (query)
                    else:
                        error_message = "Unknown error"
                    unprocessed_records_q.put(ErrorReport(query, error_message))
                    if not master_conf['general']['variables']['evaluate_all']:
                        break
                    else:
                        continue
                logger.debug("Putting %s on the queue" % (record.query_accession_id))
                unprocessed_records_q.put(record)
                if not master_conf['general']['variables']['evaluate_all']:
                    break
        unprocessed_records_q.put(QUEUE_CAP)
    except KeyboardInterrupt:
        logger.critical("KeyboardInterrupt recieved during record fetching")
        return
    except EOFError:
        logger.critical("EOFError recieved during record fetching")
        return

