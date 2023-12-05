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

import hmmer_utils, time
import csv, subprocess, os, re
import logging 
import random
import string
import pathlib
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ripp_modules.VirtualRipp import execute
import prodigal_processing
from rodeo_main import VERBOSITY
from ripp_modules.VirtualRipp import get_radar_score

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

# PFAMs that should never be added as potential ripps
no_ripp_pfams = []

class Sub_Seq(object):
        """Useful for storing subsequences and their coordinates"""
        def __init__(self, seq_type, seq, start, end, direction, accession_id=None, score=None):
            self.start = start
            self.end = end
            if direction == 1 or direction == "+":
                self.direction = "+"
            else:
                self.direction = "-"
            self.sequence = seq
            self.accession_id = accession_id
            self.radar_score = -1
            self.upstream_sequence ="xxxxx"
            self.type = seq_type ##aa, nt etc.
            self.isRRE = False
            self.pfam_descr_list = []
            self.score = score

class My_Record(object):
    """ """
    #TODO get genus and species frecom gb file
    def __init__(self, query_accession_id):
        self.query_accession_id = query_accession_id
        self.query_short = query_accession_id.replace(".", "_").replace("\t", "_").replace("/", "_")[:20] #split(".")[0].split("\t")[0]
        self.random_tag = ''.join(random.choices(string.ascii_uppercase + string.digits, k=10))
        self.peptide_types = []
        self.cluster_accession = ""
        self.cluster_sequence = ""
        self.cluster_length = ""
        self.query_index = -1
        self.genus = ""
        self.species = ""
        self.CDSs = []
        self.intergenic_seqs = []
        self.intergenic_orfs = []
        self.hypothetical_seqs = []
        self.hypothetical_cds = []
        self.cds_start_list = []
        self.cds_end_list = []
        self.bait_iteration = -1
        self.prod_window_start = 0
        self.prod_window_end = 0
        self.window_start = 0
        self.window_end = 0
        self.start_codons = ['ATG','GTG', 'TTG']
        self.stop_codons = ['TAA','TAG','TGA']
        self.ripps = {}
        
    
    def _get_query_index(self):
        """Get the index of the query CDS in the list of CDSs"""
        query_index = 0
        for cds in self.CDSs:
            if cds.accession_id == self.query_accession_id:
                return query_index
            else:
                query_index += 1
                
    def _clean_CDSs(self):
        """Called after trimming to get rid of CDSs not in the window"""
        i = 0
        while i < len(self.CDSs):
            #TODO entire CDS outside of window
            if self.CDSs[i].end < self.window_start or self.CDSs[i].start < self.window_start \
                or self.CDSs[i].start > self.window_end or self.CDSs[i].end > self.window_end:
                del self.CDSs[i]
            else:
                i += 1
    
    def trim_for_prodigal(self, n=50000):
        """Trim the window down to -n nucleotides of the start of the 
        query CDS and +n nucleotides of the end of the CDS"""
        query_index = self.query_index
        if query_index == -1:
            return
        self.prod_window_start = max(0, self.CDSs[query_index].start - n)
        self.prod_window_end = min(len(self.cluster_sequence), 
                              self.CDSs[query_index].end + n)

    #TODO cutoff or keep if in middle of gene?
    def trim_to_n_nucleotides(self, n):
        """Trim the window down to -n nucleotides of the start of the 
        query CDS and +n nucleotides of the end of the CDS"""
        query_index = self.query_index
        self.fetch_n = n
        if query_index == -1:
            return
        self.window_start = max(0, min(self.CDSs[query_index].start, self.CDSs[query_index].end) - n)
        self.window_end = min(len(self.cluster_sequence), 
                              max(self.CDSs[query_index].start, self.CDSs[query_index].end) + n)
        self._clean_CDSs()

    def trim_to_n_nucleotides_nuc(self, n, iteration):
        """Trim the window down to -n nucleotides of the start of the 
        query CDS and +n nucleotides of the end of the CDS"""
        self.fetch_n = n
        try:
            if self.bait_iteration == -1:
                self.window_start = 0
                self.window_end = min(2*n, len(self.cluster_sequence))
            else:
                self.window_start = max(0, min(self.cds_start_list[self.bait_iteration]-n, self.cds_end_list[self.bait_iteration]-n))
                self.window_end = min(len(self.cluster_sequence), 
                                      max(self.cds_start_list[self.bait_iteration]+n, self.cds_end_list[self.bait_iteration]+n))
        except OSError:
            logger.error("Error finding bait hmm hit")

        self._clean_CDSs()
        
    def trim_to_n_orfs(self, n, fetch_distance):
        """Trims the window to +- n CDSs from the query CDS, adding on "fetch_distance"
        nucleotides to each end of the window past the end CDSs"""
        query_index = self.query_index
        self.fetch_n = n
        if query_index == -1:
            return
        first_cds = max(0, query_index - n)
        last_cds = min(len(self.CDSs)-1, query_index + n)
        self.window_start = min(self.CDSs[first_cds].start, self.CDSs[first_cds].end,
                                self.CDSs[last_cds].start, self.CDSs[last_cds].end)
        self.window_end = max(self.CDSs[first_cds].start,  self.CDSs[first_cds].end,
                              self.CDSs[last_cds].start, self.CDSs[last_cds].end)
        self.window_start = max(0, self.window_start - fetch_distance)
        self.window_end = min(len(self.cluster_sequence), 
                              self.window_end + fetch_distance)
        self._clean_CDSs()
    
    def run_radar(self):
        for CDS in self.CDSs:
            CDS.radar_score = get_radar_score(CDS.sequence)

    def get_evalue(self, primary_hmm, cust_hmm):
        self.pfam_2_evalue = []
        for prots in self.CDSs:
            try:
                evalue_temp = hmmer_utils.get_hmmer_info(prots.sequence, primary_hmm, cust_hmm)
            except:
                logger.error("Unable to obtain results from HMMER. Please check provided Pfam path")
                evalue_temp = []
            for pfam_acc, desc, e_val, name in evalue_temp:
                self.pfam_2_evalue.append((prots.accession_id, pfam_acc, e_val, prots.start, prots.end, int(abs(prots.start - prots.end)/3)))

    def annotate_w_hmmer(self, primary_hmm, cust_hmm, min_length, max_length):
        self.pfam_2_coords = {}
        for CDS in self.CDSs:
            try:
                CDS.pfam_descr_list = hmmer_utils.get_hmmer_info(CDS.sequence, primary_hmm, cust_hmm) #Possible input for n and e_cutoff here
            except:
                logger.error("Unable to obtain results from HMMER. Please check provided Pfam path")
                CDS.pfam_descr_list = []
            if min_length <= len(CDS.sequence) <= max_length: # len(CDS.pfam_descr_list) == 0 and 
                if not any(any(fam in annot[0] for fam in no_ripp_pfams) for annot in CDS.pfam_descr_list) and not CDS.isRRE:
                    self.intergenic_orfs.append(CDS)
            for annot in CDS.pfam_descr_list:
                if any(fam in annot[0] for fam in ["PF14404", "PF14406", "PF14407", "PF14408", "PF14409", "PF12559" ,"TIGR04186"]):
                    self.intergenic_orfs.append(CDS)
                    continue
                if annot[0] not in self.pfam_2_coords.keys(): #annot[0] is the PF* key
                    self.pfam_2_coords[annot[0]] = []
                self.pfam_2_coords[annot[0]].append((CDS.start, CDS.end))

    def run_RREFinder(self, sequence, name, working_dir):
    #TODO change to temp file
        try:
            with open(working_dir + "RRE.fasta", 'w+') as tfile:
                tfile.write(">%s\n%s\n" % (name, sequence))
            command = ["RRE.py rre -o " + working_dir + " -i " + working_dir +  "RRE.fasta -t fasta -m precision"]
            try:
                out, err, retcode = execute(command)
            except OSError:
                logger.error("Could not run RREFinder")
                try:
                    os.remove(working_dir +"RRE.fasta")
                except OSError:
                    pass
                return 
            if retcode != 0:
                logger.error('RREFinder returned %d: %r', retcode,
                                err)
                return 
            # try:
                # os.remove(working_dir + "RRE.fasta")
            # except OSError:
                    # pass
        except KeyboardInterrupt:
            try:
                os.remove(working_dir + "RRE.fasta")
                return
            except OSError:
                pass
        ret = []
        try: 
            with open(working_dir + "/rre/rre_rrefam_results.txt", 'r') as rre_results:
                results = csv.DictReader(rre_results, delimiter='\t')
                for line in results:
                    ret.append(dict(line))
        except:
            ret = None
        return ret
           
    def has_RRE(self, sequence, name, evalue_thresh=1):
        #inactivate RRE finder requirement on local install by uncommenting following line:
        #return False
        
        random_tag = ''.join(random.choices(string.ascii_uppercase + string.digits, k=5))
        working_dir = "/tmp/RRE_" + random_tag + "_" + name +  "/"
        pathlib.Path(working_dir).mkdir(parents=True, exist_ok=True)
        results = self.run_RREFinder(sequence, name, working_dir)
        if results:
            for line in results:
                if float(line["E-value"]) < evalue_thresh:
                    return True
        return False

    def annotate_w_RREFinder(self):
        self.rre_present = False
        for CDS in self.CDSs:
            CDS.isRRE = self.has_RRE(CDS.sequence, CDS.accession_id) #Possible input for n and e_cutoff here
            self.rre_present = self.rre_present or CDS.isRRE
            # except:
                # logger.error("Unable to obtain results from RREFinder. Please check provided Pfam path")

    def annotate_w_hmmer_nuc(self, bait_list): #min_length, max_length):
        self.pfam_2_coords = {}
        peptide_type_list = []
        for CDS in self.hypothetical_cds:
            CDS.pfam_descr_list = hmmer_utils.get_hmmer_info(CDS.sequence, [], bait_list) #Possible input for n and e_cutoff here
            #if min_length <= len(CDS.sequence) <= max_length: # len(CDS.pfam_descr_list) == 0 and 
            #    self.intergenic_orfs.append(CDS)
            for annot in CDS.pfam_descr_list:
                if any(annot[0]):
                    if annot[0] == "PF05402":
                        CDS.isRRE = True
                    if float(annot[2]) < 1.0E-5:
                        self.cds_start_list.append(CDS.start)
                        self.cds_end_list.append(CDS.end)
                        if annot[0] in ["PF13471", "PF00733"]:
                            peptide_type_list.append(["lasso"])
                        elif annot[0] in ["TIGR04184"]:
                            peptide_type_list.append(["grasp"])
                        elif annot[0] in ["PF05147"]:
                            peptide_type_list.append(["lanthi1", "lanthi2", "lanthi3", "lanthi4"])
                        elif annot[0] in ["PF14028"]
                            peptide_type_list.append(["thio", "lanthi1"])
                        elif annot[0] in ["BorosinMT"]
                            peptide_type_list.append(["boro"])
                        else:
                            peptide_type_list.append([""])
        return peptide_type_list

    def set_intergenic_seqs(self, min_length, max_length):
        """Sets the sequences between called CDSs"""
        #First need to check if we have trimmed our sequence yet
        MIN_CUTOFF = 75 #Minimum number of intergenic nucs to be considered for ORF scanning        
        if self.window_end == 0:
            self.window_end = len(self.cluster_sequence)
        start = self.window_start
        for cds in self.CDSs:
            end = min(cds.start, cds.end)
            if end-start >= MIN_CUTOFF:
                #end == start could happen if the first cds starts at 0
                nt_seq = self.cluster_sequence[start:end]
                intergenic_sequence = Sub_Seq(seq_type='IGS',
                                                   seq=nt_seq, 
                                                   start=start,
                                                   end=end,
                                                   direction=0) #direction doesnt matter
                self.intergenic_seqs.append(intergenic_sequence)
            start = max(cds.end, cds.start)
        nt_seq = self.cluster_sequence[start:]
        end = self.window_end
        if end > start and abs(end-start) >= MIN_CUTOFF:
            intergenic_sequence = Sub_Seq(seq_type='IGS',
                                                   seq=nt_seq, 
                                                   start=start,
                                                   end=end,
                                                   direction=0) #direction doesnt matter
            self.intergenic_seqs.append(intergenic_sequence)
        return
        
    def set_intergenic_orfs(self, min_aa_seq_length, max_aa_seq_length, overlap):
        """Examines intergenic sequences to determine whether or not there are ORFs
        that code for valid aa sequences"""
        identifier = ''.join(random.choices(string.ascii_uppercase + string.digits, k=5))
        SeqIO.write([SeqRecord(
            #Seq(self.cluster_sequence[self.window_start:self.window_end]),
            self.cluster_sequence[self.window_start:self.window_end],
            id=self.cluster_accession,
            name="prodigal")], f"/tmp/{self.cluster_accession}_{identifier}.fasta", "fasta")
        os.system(f"prodigal -q -i /tmp/{self.cluster_accession}_{identifier}.fasta -p meta -a /tmp/{self.cluster_accession}_{identifier}_ig.fasta -f sco -o /dev/null")
        all_orfs = SeqIO.parse(f"/tmp/{self.cluster_accession}_{identifier}_ig.fasta", "fasta")
        for orf_record in all_orfs:
            if not (min_aa_seq_length < len(orf_record.seq) - 1 < max_aa_seq_length):
                continue
            _, start, stop, strand, _ = orf_record.description.split('#' )
            start = self.window_start + int(start)
            stop =  self.window_start + int(stop)
            strand = int(strand)
            if strand == -1:
                start, stop = stop, start
                stop = stop - 1
            else:
                start = start - 1
            assert abs(start - stop) % 3 == 0
            found_overlap = False
            for cds in self.CDSs:
                cds_s = min(cds.start, cds.end)
                cds_e = max(cds.start, cds.end)
                orf_s = min(start, stop)
                orf_e = max(start, stop)
                if cds_s == orf_s or cds_e == orf_e:
                    found_overlap = True
                    break
                elif cds_s < orf_s < cds_e and cds_e - orf_s > overlap:
                    found_overlap = True
                    break
                elif cds_s < orf_e < cds_e and orf_e - cds_s > overlap:
                    found_overlap = True
                    break
                    
            if not found_overlap:
                self.intergenic_orfs.append(Sub_Seq("ORF", 
                    str(orf_record.seq)[:-1], 
                    start,
                    stop,
                    strand))
                        
        self.intergenic_orfs.sort(key=lambda seq: seq.start)
        #Get rid of duplicates. Duplicate ORFs will appear when the overlap is
        #set such that two intergenic sequences are expanded to a point where 
        #they share nucleotides with eachother
        i = 1
        while i < len(self.intergenic_orfs):
            #print(self.intergenic_orfs[i].sequence)
            if self.intergenic_orfs[i].start == self.intergenic_orfs[i-1].start:
                del self.intergenic_orfs[i]
            i += 1


    def set_hypothetical_seqs(self):
        """Sets the sequences between called CDSs"""
        #First need to check if we have trimmed our sequence yet
        MIN_CUTOFF = 75 #Minimum number of intergenic nucs to be considered for ORF scanning        
        if self.window_end == 0:
            self.window_end = len(self.cluster_sequence)
        start = self.window_start
        for cds in self.CDSs:
            end = min(cds.start, cds.end)
            if end-start >= MIN_CUTOFF:
                #end == start could happen if the first cds starts at 0
                nt_seq = self.cluster_sequence[start:end]
                intergenic_sequence = Sub_Seq(seq_type='IGS',
                                                   seq=nt_seq, 
                                                   start=start,
                                                   end=end,
                                                   direction=0) #direction doesnt matter
                self.intergenic_seqs.append(intergenic_sequence)
            start = max(cds.end, cds.start)
        nt_seq = self.cluster_sequence[start:]
        end = self.window_end
        if end > start and abs(end-start) >= MIN_CUTOFF:
            intergenic_sequence = Sub_Seq(seq_type='IGS',
                                                   seq=nt_seq, 
                                                   start=start,
                                                   end=end,
                                                   direction=0) #direction doesnt matter
            self.hypothetical_seqs.append(intergenic_sequence)
        return

    def set_hypothetical_cds(self, min_aa_seq_length, max_aa_seq_length, overlap):
        """Examines intergenic sequences to determine whether or not there are ORFs
        that code for valid aa sequences"""
        nt_seq, nt_seq_rev = self.cluster_sequence, self.cluster_sequence.reverse_complement()
        prodigal_processing.run_prodigal(record=self, whole_contig=True)
        try: 
            prod_file = open("/tmp/%s_%sorfs.tsv" % (self.query_short, self.random_tag), 'r')
        except:
            time.sleep(2)
            prod_file = open("/tmp/%s_%sorfs.tsv" % (self.query_short, self.random_tag), 'r')
            pass
        prod_results = prod_file.readlines()
        prod_file.close()
        dup_removed_rows = {}

        if len(prod_results) > 5:
            for line in prod_results[5:]:
                tmp_line = line.split("\t")
                if len(tmp_line) != 13:
                    continue
                start, end = self.find_orf_coordinates(int(tmp_line[0]), int(tmp_line[1]), tmp_line[2])
                row = [start, end, tmp_line[2], tmp_line[3]]
                if float(tmp_line[3]) < -10:
                    continue
                if tmp_line[2]+str(end) in dup_removed_rows:
                    if float(dup_removed_rows[tmp_line[2]+str(end)][3]) < float(row[3]):
                        dup_removed_rows[tmp_line[2]+str(end)] = row
                else:
                    dup_removed_rows[tmp_line[2]+str(end)] = row

        for key in dup_removed_rows:
            if dup_removed_rows[key][2] == "+":
                nt_subsequence = nt_seq[dup_removed_rows[key][0]:dup_removed_rows[key][1]]
                if nt_subsequence[-3:] in ["TAG", "TAA", "TGA"]:
                    nt_subsequence = nt_subsequence[:-3]
                aa_sequence = nt_subsequence.translate(11)
                self.hypothetical_cds.append(Sub_Seq('ORF', aa_sequence, dup_removed_rows[key][0], dup_removed_rows[key][1], direction=1, score=float(dup_removed_rows[key][3])))
            else:
                nt_subsequence = nt_seq_rev[len(nt_seq_rev)-dup_removed_rows[key][0]: len(nt_seq_rev)-dup_removed_rows[key][1]]
                if nt_subsequence[-3:] in ["TAG", "TAA", "TGA"]:
                    nt_subsequence = nt_subsequence[:-3]
                aa_sequence = nt_subsequence.translate(11)
                self.hypothetical_cds.append(Sub_Seq('ORF', aa_sequence, dup_removed_rows[key][0], dup_removed_rows[key][1], direction=-1, score=float(dup_removed_rows[key][3])))
        
        self.hypothetical_cds.sort(key=lambda seq: abs(seq.end-seq.start), reverse=True)
        time.sleep(0.1)               
        self.hypothetical_cds.sort(key=lambda seq: seq.end)
        #Get rid of duplicates. Duplicate ORFs will appear when the overlap is
        #set such that two intergenic sequences are expanded to a point where 
        #they share nucleotides with eachother
        i = 1
        prev_len = len(self.hypothetical_cds)
        while i < len(self.hypothetical_cds):
            #print(self.intergenic_orfs[i].sequence)
            j = 0
            #print(i, j, len(self.hypothetical_cds))
            while i+j < len(self.hypothetical_cds):
                if self.hypothetical_cds[i+j].end == self.hypothetical_cds[i-1].end:
                    #print (self.hypothetical_cds[i+j].__dict__, self.hypothetical_cds[i-1].__dict__, "deleted")
                    del self.hypothetical_cds[i+j]
                    j += 1
                else:
                    break
            i += 1+j
            if i >= len(self.hypothetical_cds) and len(self.hypothetical_cds) != prev_len:
                i = 1
                prev_len = len(self.hypothetical_cds)
        return
    
    def filter_RREs_and_HMMs(self, hmm_list):
        ret_list = []
        for cds in self.intergenic_orfs:
            if not any(any(fam in annot[0] for fam in hmm_list) for annot in cds.pfam_descr_list) and not cds.isRRE:
                ret_list.append(cds)
        self.intergenic_orfs = ret_list

    def set_ripps(self, module, master_conf):
        logger.debug("Setting %s ripps for %s" % (module.peptide_type, self.query_accession_id))
        self.ripps[module.peptide_type] = []
        for orf in self.intergenic_orfs:
            if module.peptide_type == "grasp":
                orf.radar_score = get_radar_score(orf.sequence)
            if master_conf[module.peptide_type]['variables']['precursor_min'] <= len(orf.sequence) <=  master_conf[module.peptide_type]['variables']['precursor_max'] \
                    or ("M" in orf.sequence[-master_conf[module.peptide_type]['variables']['precursor_max']:]) \
                    or (module.peptide_type == "grasp" and len(orf.sequence) < 400)\
                    or (module.peptide_type == "grasp" and orf.radar_score > 0 and len(orf.sequence) > 400):
                if module.peptide_type == "sacti":
                    ripp = module.Ripp(orf.start, orf.end, str(orf.sequence), orf.upstream_sequence, self.pfam_2_coords, self.rre_present)
                
                elif module.peptide_type == "boro" or "grasp":
                    ripp = module.Ripp(orf.start, orf.end, str(orf.sequence), orf.upstream_sequence, self.pfam_2_coords, self.pfam_2_evalue) 
                
                else:
                    ripp = module.Ripp(orf.start, orf.end, str(orf.sequence), orf.upstream_sequence, self.pfam_2_coords)
                if module.peptide_type == "grasp":
                    ripp.radar_score = orf.radar_score
                ripp.radar_score  = orf.radar_score
                if ripp.valid_split or master_conf[module.peptide_type]['variables']['exhaustive']:
                    self.ripps[module.peptide_type].append(ripp)
                
    def score_ripps(self, module, pfam_hmm, cust_hmm):
        logger.debug("Scoring %s ripps for %s" % (module.peptide_type, self.query_accession_id))
        for ripp in self.ripps[module.peptide_type]:
            ripp.set_score(pfam_hmm, cust_hmm)
    
    def color_ripps(self, module):
        logger.debug("Setting confidence for %s ripps for %s" % (module.peptide_type, self.query_accession_id))
        for ripp in self.ripps[module.peptide_type]:
            ripp.confidence = float(ripp.score)/(ripp.CUTOFF)
                
    def print_info(self):
        print("="*50)
        counter = 0
        print("CDSs")
        for sub_seq in self.CDSs:
            print(counter)
            counter+=1
            if sub_seq.accession_id == self.query_accession_id:
                print("QUERY Accession:  " + sub_seq.accession_id)
            else:
                print("Accession:  " + sub_seq.accession_id)
            print("Coords:  " + str(sub_seq.start) + " to " +  str(sub_seq.end))
        print("="*50)
        print("IGSs")
        for sub_seq in self.intergenic_seqs:
            print(counter)
            counter+=1
            print("Coords:  " + str(sub_seq.start) + " to " +  str(sub_seq.end))
        print("="*50)
        print("Intergenic ORFs")
        for sub_seq in self.intergenic_orfs:
            print(counter)
            counter+=1
            if sub_seq.end < sub_seq.start:
                strand = -1
            else:
                strand = 1
            print("Potential ORF of length " + str(len(sub_seq.sequence)-1) + 
                  " found at " + str(sub_seq.start) + ":" + str(sub_seq.end) +
                  " on strand " + str(strand))
            print(sub_seq.sequence + '\n')
        print("="*50)


    def find_prod_coordinates(self, beg, end):
        if(beg<end):
            return(beg-self.window_start+1, end-self.window_start)
        else:
            return(beg-self.window_start+2, end-self.window_start-1)

    def find_orf_coordinates(self, beg, end, direction):
        if direction == "+":
            return(beg-1, end)
        else:
            return(end, beg-1)

        
    
def update_score_w_svm(output_dir, records):
        """Order should be preserved. Goes through file and updates scores"""
        for peptide_type in records[0].ripps.keys():
            score_reader = csv.reader(open(output_dir + '/' + peptide_type + '/' +\
                                           peptide_type + '_features.csv')) 
            header = next(score_reader)
            score_col = 6
            try:
                score_col = header.index("Total Score")
            except ValueError:
                logger.error("Temporary CSV format invalid. No column named \"Total Score\". Score results are most likely invalid.")
            
            score_reader_done = False
            total_ripps = 0
            for record in records:
                total_ripps += len(record.ripps[peptide_type])
                for ripp in record.ripps[peptide_type]:
                    if not score_reader_done:
                        try:
                            line = next(score_reader)
                        except KeyboardInterrupt:
                            raise KeyboardInterrupt
                        except Exception as e:
                            import traceback as tb
                            tb.print_exc()
                            print(e)
                            score_reader_done = True
                            logger.warning("Mismatch in RiPP count and length of CSV. Score results are most likely invalid")
                            print(total_ripps)
                            return
                    ripp.score = int(line[score_col])
                    ripp.confidence = float(ripp.score)/(ripp.CUTOFF)

