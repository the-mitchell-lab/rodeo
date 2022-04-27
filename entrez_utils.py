# -*- coding: utf-8 -*-
#==============================================================================
# NCBI utilities policies are available here:
# https://www.ncbi.nlm.nih.gov/home/about/policies/
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

import traceback
from Bio import Entrez, SeqIO
from My_Record import My_Record, Sub_Seq
import logging
from rodeo_main import VERBOSITY
import time
from timeout_decorator import timeout, TimeoutError

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

Entrez.email ='kille2@illinois.edu' 

@timeout(300)
def get_gb_handles(prot_accession_id):
    """Returns a list of .gb/.gbk filestreams from protein db accession.
    
    ERROR CODES:
        -1 = No results in protein db for Esearch on prot_accession_id
        -2 = No results in nuccore db for value obtained from protein db
        -3 = Any response failure from Entrez database (error on database side)
    """
    tries = 3 #Max number of times to try the database before it gives up
    logger.info("Fetching %s data from GenBank." % prot_accession_id)
    for i in range(tries):
        try:
            record = Entrez.read(Entrez.esearch("protein",term=prot_accession_id))
            total_count = record["Count"]
            if int(total_count) < 1:
                return -1
                
            IdList = record["IdList"]
            time.sleep(0.5)
            link_records = Entrez.read(Entrez.elink(dbfrom="protein",db="nuccore",id=IdList))
            nuccore_ids=[]
            if len(link_records[0]['LinkSetDb']) == 0:
                return -2
            for record in link_records[0]['LinkSetDb'][0]['Link']:
                nuccore_ids.append(record['Id']) 
                
            search_handle     = Entrez.epost(db="nuccore", id=",".join(nuccore_ids))
            search_results    = Entrez.read(search_handle)
            webenv,query_key  = search_results["WebEnv"], search_results["QueryKey"] 
            
            batchSize = 1
            
            handles = []
            for start in range(len(nuccore_ids)):
                time.sleep(0.5)
                orig_handle = Entrez.efetch(db="nuccore", dbfrom="protein", rettype="gbwithparts", 
                                               retmode="text", retstart=start, retmax=batchSize, 
                                               webenv=webenv, query_key=query_key)
                handles.append(orig_handle)
            return handles
        except KeyboardInterrupt:
            raise KeyboardInterrupt
        except TimeoutError:
            logger.error("Timeout while reaching genbank for %s." % (prot_accession_id))
            return -3
        except Exception as e:
            time.sleep(0.5)
            logger.error("Failed to fetch record for %s." % (prot_accession_id))
            logger.error(e)
            pass
        time.sleep(2)
    return -3

#gb_handle should only be a handle to ONE query 
@timeout(300)
def get_record_from_gb_handle(gb_handle, nuccore_accession_id):
    """Takes an input gb_filestream and query accession_id.
    Returns a record containing basic information about the query.
    i.e. CDSs, genus/species, sequence.
    
    ERROR CODES:
        -1 = Couldn't process .gb filestream for some reason...
    """
    logger.info("Parsing %s handle." % nuccore_accession_id)
    try:
        try:
            gb_record = SeqIO.parse(gb_handle, "genbank")
        except KeyboardInterrupt:
            raise KeyboardInterrupt
        except Exception as e:
            logger.error("Error parsing GenBank handle for %s" % (nuccore_accession_id))
            logger.error(e)
        any_record = False
        accession_found = False
        for record in gb_record: #Should only be one record in gb_record
#            if not first_record:
#                logger.info("Multiple records in gb_record for %s." % (nuccore_accession_id))
            any_record = True
            ret_record = My_Record(nuccore_accession_id)
            ret_record.cluster_accession = record.id
            ret_record.cluster_sequence = record.seq
            ret_record.cluster_genus_species = record.annotations['organism'] #TODO verify
            for feature in record.features:
                if feature.type == 'CDS':
                    inferred = False
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    seq = "X" * int((abs(start-end)/3))
                    direction = 1
                    if end < start:
                        direction = -1
                    locus_tag = ""
                    accession_id = ""
                    if 'locus_tag' in feature.qualifiers.keys():
                        locus_tag = feature.qualifiers['locus_tag'][0]
                    if 'protein_id' in feature.qualifiers.keys():
                        accession_id = feature.qualifiers['protein_id'][0]
                        direction = feature.strand
                        if direction == -1:
                            tmp = start
                            start = end
                            end = tmp
                    if accession_id == "":
                        if locus_tag == "":
                            continue
                        else:
                            inferred = True
                            if 'inference' in feature.qualifiers.keys():
                                inference = feature.qualifiers['inference'][0]
                                splits = inference.split(':')
                                if splits[-2] == "RefSeq": 
                                    accession_id = splits[-1]
                                else:
                                    accession_id = locus_tag
                            else:
                                accession_id = locus_tag
                    if nuccore_accession_id.split('.')[0].upper() == accession_id.split('.')[0].upper() or \
                       nuccore_accession_id.split('.')[0].upper() == locus_tag.split('.')[0].upper():
                        ret_record.query_index = len(ret_record.CDSs)
                        accession_found = True
                    if 'translation' in feature.qualifiers.keys():
                        seq = feature.qualifiers['translation'][0]
                    
                    cds = Sub_Seq(seq_type='CDS', seq=seq, start=start, end=end, direction=direction, accession_id=accession_id)
                    if inferred:
                        cds.inferred = True
                    else:
                        cds.inferred = False
                    ret_record.CDSs.append(cds)
            if accession_found and any_record and len(ret_record.cluster_sequence) > 0 and len(ret_record.CDSs) > 0:
                logger.debug("Record made for %s" % (ret_record.cluster_accession))
                
                return ret_record
        if not accession_found:
            logger.error("Accession %s not found." % (nuccore_accession_id))
            return -1
        if any_record and len(ret_record.cluster_sequence) > 0 and len(ret_record.CDSs) > 0:
            logger.debug("Record made for %s" % (ret_record.cluster_accession))
            return ret_record
        else:
            logger.error("Corrupted GenBank record for %s" % (nuccore_accession_id))
            traceback.print_exc()
            return -1
    except KeyboardInterrupt:
        raise KeyboardInterrupt
    except TimeoutError:
        logger.error("Timeout while reaching genbank for %s" % (nuccore_accession_id))
        return -1
    except Exception as e:
        logger.error("Corrupted GenBank record for %s" % (nuccore_accession_id))
        logger.error(e)
        traceback.print_exc()
        return -1
    return -1
