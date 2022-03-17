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

import hmmer_utils
import csv, subprocess, os
import logging 
from rodeo_main import VERBOSITY
from ripp_modules.VirtualRipp import get_radar_score
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random, string

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

class Sub_Seq(object):
        """Useful for storing subsequences and their coordinates"""
        def __init__(self, seq_type, seq, start, end, direction, accession_id=None):
            self.start = start
            self.end = end
            if direction == 1:
                self.direction = "+"
            else:
                self.direction = "-"
            self.sequence = seq
            self.accession_id = accession_id
            self.radar_score = -1
            self.upstream_sequence ="xxxxx"
            self.type = seq_type ##aa, nt etc.

class My_Record(object):
    """ """
    #TODO get genus and species frecom gb file
    def __init__(self, query_accession_id):
        self.query_accession_id = query_accession_id
        self.query_short = query_accession_id.split(".")[0]
        self.cluster_accession = ""
        self.cluster_sequence = ""
        self.cluster_length = ""
        self.query_index = -1
        self.genus = ""
        self.species = ""
        self.CDSs = []
        self.intergenic_seqs = []
        self.intergenic_orfs = []
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
            
    def annotate_w_hmmer(self, primary_hmm, cust_hmm, min_length, max_length):
        self.pfam_2_coords = {}
        for CDS in self.CDSs:
            try:
                CDS.pfam_descr_list = hmmer_utils.get_hmmer_info(CDS.sequence, primary_hmm, cust_hmm) #Possible input for n and e_cutoff here
            except:
                logger.error("Unable to obtain results from HMMER. Please check provided Pfam path")
                CDS.pfam_descr_list = []
            if min_length <= len(CDS.sequence) <= max_length: # len(CDS.pfam_descr_list) == 0 and 
                self.intergenic_orfs.append(CDS)
            for annot in CDS.pfam_descr_list:
                if any(fam in annot[0] for fam in ["PF14404", "PF14406", "PF14407", "PF14408", "PF14409", "PF12559" ,"TIGR04186"]):
                    self.intergenic_orfs.append(CDS)
                    continue
                if annot[0] not in self.pfam_2_coords.keys(): #annot[0] is the PF* key
                    self.pfam_2_coords[annot[0]] = []
                self.pfam_2_coords[annot[0]].append((CDS.start, CDS.end))
                

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
            Seq(self.cluster_sequence[self.window_start:self.window_end]),
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
                if cds_s == orf_s and cds_e == orf_e:
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
        # with open(f"{self.cluster_accession}_ig.txt") as orf_f:
            # for line in orf_f:
                # if line[0] != ">":
                    # continue
                # idx, start, end, strand = line.rstrip().split('_')
                # if strand = '-':
                    # start, end = self.window_start + int(end), self.window_start + int(start)
        # print(prodigal_record)
            
        # for strand, sequence in [(1, self.cluster_sequence),
                                         # (-1, self.cluster_sequence.reverse_complement())]:
            # for intergenic_seq in self.intergenic_seqs:
                # #Do start codons iteratively to ensure you don't skip any
                # #In theory, should only be a 3x slowdown in runtime MAX.
                # for start_codon in self.start_codons:
                    # if strand == -1:
                        # next_start = max(len(sequence) - intergenic_seq.end - overlap, 
                                         # len(sequence)-self.window_end)
                        # intergenic_seq_end = min(len(sequence) - intergenic_seq.start + overlap,
                                                 # len(sequence) - self.window_start)
                    # else:
                        # next_start = max(intergenic_seq.start - overlap, self.window_start)
                        # intergenic_seq_end = min(self.window_end, intergenic_seq.end + overlap)
                    # start = 0
                    # #Stay in the loop until we can't find a stop codon
                    # while start < intergenic_seq_end:
                        # start = sequence.find(start_codon, next_start)
                        # if start == -1:
                            # #Can't find any more of this codon
                            # break
                        # next_start = start + 1
                        # end = start 
                        # found_stop = False
                        # while end < len(sequence) and end < intergenic_seq_end:
                            # codon = sequence[end:end+3]
                            # if str(codon) in self.stop_codons:
                                # found_stop = True
                                # break
                            # else:
                                # end = end + 3
                        # #If we didn't find a stop codon, do what?
                        # if not found_stop:
                            # continue
                        # if not (min_aa_seq_length < (end-start)/3 < max_aa_seq_length):
                            # continue
                        # nt_subsequence = sequence[start:end] #Add 3 if you want to inclue stop codon
                        # aa_sequence = nt_subsequence.translate(11)
                        # #Get nucleotide coords for original strand
                        # if strand == -1:
                            # old_end = end
                            # end = len(sequence) - start
                            # old_start = len(sequence) - old_end - 3
                            # potential_orf = Sub_Seq('ORF',
                                                         # aa_sequence,
                                                         # end,
                                                         # old_start,
                                                         # -1)
                        # else:
                            # end = end + 3
                            # potential_orf = Sub_Seq('ORF',
                                                         # aa_sequence,
                                                         # start,
                                                         # end,
                                                         # 1)
                        # if potential_orf.start > potential_orf.end:
                            # upstream_sequence = str(self.cluster_sequence[potential_orf.start+4:potential_orf.start+13].reverse_complement())
                        # else:
                            # upstream_sequence = str(self.cluster_sequence[potential_orf.start-13:potential_orf.start+4])
                        # potential_orf.upstream_sequence = upstream_sequence
                        # self.intergenic_orfs.append(potential_orf)
                        
        # self.intergenic_orfs.sort(key=lambda seq: seq.start)
        # #Get rid of duplicates. Duplicate ORFs will appear when the overlap is
        # #set such that two intergenic sequences are expanded to a point where 
        # #they share nucleotides with eachother
        # i = 1
        # while i < len(self.intergenic_orfs):
            # #print(self.intergenic_orfs[i].sequence)
            # if self.intergenic_orfs[i].start == self.intergenic_orfs[i-1].start:
                # del self.intergenic_orfs[i]
            # i += 1
        
       
        return
    
    def set_ripps(self, module, master_conf):
        logger.debug("Setting %s ripps for %s" % (module.peptide_type, self.query_accession_id))
        self.ripps[module.peptide_type] = []
        for orf in self.intergenic_orfs:
            if module.peptide_type == "grasp":
                orf.radar_score = get_radar_score(orf.sequence)
            if master_conf[module.peptide_type]['variables']['precursor_min'] <= len(orf.sequence) <=  master_conf[module.peptide_type]['variables']['precursor_max'] \
                    or ("M" in orf.sequence[-master_conf[module.peptide_type]['variables']['precursor_max']:]) \
                    or (module.peptide_type == "grasp" and orf.radar_score > 0 and len(orf.sequence) < 400):
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

