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

from entrez_utils import get_gb_handles, get_record_from_gb_handle
from My_Record import Sub_Seq
import logging
from rodeo_main import VERBOSITY, QUEUE_CAP
import traceback
import sys
import prodigal_processing
import os

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
                if record.bait_iteration > -1:
                    logger.info("Worker process %s is processing locus %d in %s" % (my_id, record.bait_iteration+1, record.cluster_accession))
                else:
                    logger.info("Worker process %s is processing %s" % (my_id, record.query_accession_id))
                if args.prodigal:
                    record.trim_for_prodigal()
                if args.meta:
                    if len(record.cds_start_list) > record.bait_iteration:
                        record.trim_to_n_nucleotides_nuc(10000, record.bait_iteration)
                elif master_conf['general']['variables']['fetch_type'].lower() == 'cds':
                    record.trim_to_n_orfs(master_conf['general']['variables']['fetch_n'], master_conf['general']['variables']['fetch_distance'])
                elif master_conf['general']['variables']['fetch_type'].lower() == 'nucs':
                    record.trim_to_n_nucleotides(master_conf['general']['variables']['fetch_n'])
                if ("grasp" in record.peptide_types and args.meta==False) or ("grasp" in args.peptide_types and args.meta==False):
                    record.run_radar(args.output_dir)
                record.pfam_2_evalue = []
                if ("boro" in args.peptide_types and args.meta==False) or ("grasp" in args.peptide_types and args.meta==False) or ("cyclo" in args.peptide_types and args.meta==False) or "boro" in record.peptide_types or "grasp" in record.peptide_types or "cyclo" in record.peptide_types:
                    record.get_evalue(master_conf['general']['variables']['pfam_dir'], args.custom_hmm, args.output_dir)
                record.annotate_w_hmmer(master_conf['general']['variables']['pfam_dir'], args.custom_hmm, args.output_dir,
                                        min_length=master_conf['general']['variables']['precursor_min'], 
                                        max_length=master_conf['general']['variables']['precursor_max'])
                if args.meta:
                    k = 0
                    while k < len(record.CDSs):
                        if record.CDSs[k].accession_id.isdigit() and ((record.CDSs[k].pfam_descr_list == [] and record.CDSs[k].score < 50) or    #double-check this number before broad use through some data analysis
                            any((x.end == record.CDSs[k].end or x.start == record.CDSs[k].start) for x in record.CDSs[:k])):
                            del record.CDSs[k]
                        else:
                            k += 1
                #if args.megarun == False:
                #    record.annotate_w_RREFinder()
                record.set_intergenic_seqs(min_length=master_conf['general']['variables']['precursor_min'], 
                                           max_length=master_conf['general']['variables']['precursor_max'])
                record.set_intergenic_orfs(min_aa_seq_length=master_conf['general']['variables']['precursor_min'], 
                                           max_aa_seq_length=master_conf['general']['variables']['precursor_max'],
                                           overlap=master_conf['general']['variables']['overlap'], output_dir=args.output_dir)
                if args.prodigal:
                    prodigal_processing.run_prodigal(record, args.output_dir)
                #temp_peptide_types = args.peptide_types
                #print(args.peptide_types)
                #if args.megarun:
                #    temp_peptide_types = record.peptide_types
                for peptide_type in args.peptide_types:
                    if args.meta and peptide_type not in record.peptide_types:
                        continue
                    module = ripp_modules[peptide_type]
                    if peptide_type == "lasso":
                        record.filter_RREs_and_HMMs(hmm_list=list(master_conf[peptide_type]["pfam_colors"].keys()))
                    elif peptide_type == "boro" and master_conf[peptide_type]['variables']['skip_mt']:
                        record.filter_RREs_and_HMMs(hmm_list=["PF00590", "NMT_1", "NMT_2", "BoroMT"])
                    #Filter out enzymes in grasp BGCs from scoring
                    elif peptide_type == "grasp":
                        record.filter_RREs_and_HMMs(hmm_list=["Graspetide_synthetase", "TIGR04187", "TIGR04192", "TIGR04193", "TIGR04188", "TIGR04364", "PF01135", "TIGR00080", "PF00583", "PF13673", "PF13302", "PF13523"])
                    record.set_ripps(module, master_conf, args.output_dir)
                    record.score_ripps(module, master_conf['general']['variables']['pfam_dir'], args.custom_hmm, args.output_dir)
                    record.color_ripps(module)
                if record.bait_iteration > -1:
                    logger.debug("Worker process %s finished processing locus %d in %s" % (my_id, record.bait_iteration+1, record.query_accession_id))
                else:
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
            if '.gbk' != query[-4:] and '.gb' != query[-3:] and '.fna' != query[-4:] and '.fa' != query[-3:] and '.fasta' != query[-6:] and '.seq' != query[-4:] and '.gbff' != query[-5:]: #accession_id
                gb_handles = get_gb_handles(query, master_conf, meta=args.meta)
                nuccore_accession = query
                if type(gb_handles) is int:
                    if gb_handles == -1:
                        if args.meta:
                            error_message = "No results in nuccore db for Esearch on %s" % (query)
                        else:
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
                if args.meta:
                    if '.fna' == query[-4:] or '.fa' == query[-3:] or '.fasta' == query[-6:]: 
                        records = get_record_from_gb_handle(handle, query, master_conf, meta=True, fasta=True, self_annotate=args.self_annotate, megarun=args.megarun)
                    else:
                        records = get_record_from_gb_handle(handle, query, master_conf, meta=True, self_annotate=args.self_annotate, megarun=args.megarun)
                    try:
                        record = records[0]
                    except:
                        record = -2
                else:
                    record = get_record_from_gb_handle(handle, nuccore_accession, master_conf)
                if type(record) is int:
                    if record == -1:
                        error_message = "Couldn't process %s Genbank filestream. May be corrupt."\
                          % (query)
                    elif record == -2:
                        error_message = "No records of sufficient length detected. Couldn't process %s filestream." % (query)
                    else:
                        error_message = "Unknown error"
                    unprocessed_records_q.put(ErrorReport(query, error_message))
                    if not master_conf['general']['variables']['evaluate_all']:
                        break
                    else:
                        continue
                logger.debug("Putting %s on the queue" % (record.query_accession_id))
                if args.meta:
                    for record in records:
                        temp_peptide_types = []
                        if record.CDSs == []:
                            record.set_hypothetical_cds(50, 20000, 1000, args.output_dir)
                        else:
                            record.hypothetical_cds = record.CDSs
                        if args.bait_list:
                            record.annotate_w_hmmer_nuc(args.bait_list, args.output_dir, args.peptide_types)
                        if len(record.cluster_sequence) < 20000000:
                            i = 1
                            if record.CDSs == []:
                                for entry in record.hypothetical_cds:
                                    cds = Sub_Seq(seq_type='CDS', seq=entry.sequence, start=entry.start, end=entry.end, direction=entry.direction, accession_id=str(i), score=entry.score)
                                    cds.inferred = False
                                    record.CDSs.append(cds)
                                    i += 1
                        try:
                            os.remove("%s/tmp/%s_%sorfs.tsv" % (args.output_dir, record.query_short, record.random_tag))
                        except:
                            pass
                        if record.cds_start_list == [] and not args.bait_list:
                            unprocessed_records_q.put(record)
                        else:
                            for i in range(0,len(record.cds_start_list)):
                                record.bait_iteration = i
                                record.peptide_types = list(set(record.peptide_type_list[record.bait_iteration]))
                                unprocessed_records_q.put(record)
                else:
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
