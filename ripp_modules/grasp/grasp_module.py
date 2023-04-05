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
#/
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
#
#==============================================================================
# Special thanks goes to AntiSmash team, whose antiSMASH-rodeo repository at
# https://bitbucket.org/mmedema/antismash-rodeo/ provided the backbone code for 
# a great deal of the heuristic calculations.
#==============================================================================

import csv
import os
import re
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from ripp_modules.VirtualRipp import VirtualRipp
import hmmer_utils
from collections import defaultdict
import pathlib
FILE_DIR = pathlib.Path(__file__).parent.absolute()

peptide_type = "grasp"
CUTOFF = 12
index = 0


def write_csv_headers(output_dir):
    dir_prefix = output_dir + '/grasp/'
    if not os.path.exists(dir_prefix):
        os.makedirs(dir_prefix)
    svm_headers = 'Precursor Index,classification, Radar Score, <1100 from ATP-grasp ligase, <300 from ATP-grasp ligase, >2000 from ATP-grasp ligase, 2nd half contains >8% Asp residue (acceptor site), 2n half of precursor contains >7% Thr residues (donor site), Full length precursor contains >7% Lys residues (donor site), Second half of precursor contains >7% Pro residues, Precursor contains >0 Cys residues, Second half of precursor contains <3 acceptor residues (Asp + Glu), Second half of precursor contains <3 donor residues (Ser + Thr + Lys), Second half of precursor contains more donor residues than acceptor residues, % acceptor residues in second half of precursor > % acceptor residues in first half of precursor,  % acceptor residues (Asp + Glu) in second half of precursor >13%, % donor residues (Ser + Thr + Lys) in second half of precursor >18%, Precursor ends with Leu, Precursor ends with Asp Val or Ser, First half of the precursor contains a “PFxL” motif, Precursor and the ATP-grasp protein (MdnC homolog) are encoded on same strand, Gene cluster contains one of the following: PF0005 PF06472 PF00664 PF03412 TIGR03796 TIGR01846 TIGR03797 TIGR00954 TIGR02203 TIGR02204 TIGR03375 (ABC transporters that co-occur frequently), Gene cluster contains one of the following: PF13302 PF00583 PF13523 (acetyltransferases that co-occur frequently), A local gene product hits TIGR04188 (methyltransferase), Precursor hits PF12559 (serine endopeptidase inhibitor) or TIGR04186 (GRASP_targ), Precursor hits PF14404 (strep_pep) PF14406 (bacteroid_pep) PF14407 (frankia_pep) PF14408 (actino_pep) or PF14409 (herpeto_pep), Acceptor compressionC index > Donor compressionC index*, Calculated charge at pH7 of second half of precursor < charge of first half of precursor, Precursor peptide contains sequence motif #1 “PFFAxxL”, Count of sequence motif #2 “TxKxPSDx(E/D)(E/D)”, Precursor peptide contains sequence motif #3 “LFIx(D/E)L”, Count of sequence motif #4 “xTxxxxEx” in precursor, Precursor contains of sequence motif #5 “KPYxxxYxE” in precursor, Count of sequence motif “TxxTxxxExxDxD” in precursor, Count of sequence motif “TxxTxTxExxDxxE” in precursor, Count of sequence motif “TKTxxxxExDD” in precursor,Count of sequence motif “TxTxTxxxxExxDxD” in precursor, Minimum distance from MdnB C D TgnB and PsnB homologs  (nt), Precursor length, Estimated first half charge, Estimated second half charge, Estimated precursor charge,Absolute value of second half charge,Absolute value of first half charge,Absolute value of precursor charge,FIRST HALF A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,Aromatics,Neg charged,Pos charged,Charged,Aliphatic,Hydroxyl,SECOND HALF A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,Aromatics,Neg charged,Pos charged,Charged,Aliphatic,Hydroxyl,PRECURSOR A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V, LAST RESIDUE A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,Aromatics,Neg charged,Pos charged,Charged,Aliphatic,Hydroxyl, acceptor compression index, donor compression index, total motifs hit,meme1,meme2,meme3,meme4,meme5,meme6,meme7,meme8,meme9,meme10,meme11,meme12,meme13,meme14,meme15,meme16,meme17,meme18,meme19,meme20,meme21,meme22, sum of meme scores, no motifs present '
    svm_headers = svm_headers.split(',')
    features_headers = ['Accession_id', 'Genus/Species', 'First half', 'Second Half', 'Start', 'End', "Total Score", "Valid Precursor"] + svm_headers 
    features_csv_file = open(dir_prefix + "temp_features.csv", 'w')
    svm_csv_file = open("{}/svm/fitting_set.csv".format(FILE_DIR), 'w')
    features_writer = csv.writer(features_csv_file)
    svm_writer = csv.writer(svm_csv_file)
    features_writer.writerow(features_headers)
    svm_writer.writerow(svm_headers)#Don't include accession_id, genus/species,
                                        #leader, core sequence, score, or svm classification
    
class Ripp(VirtualRipp):
    def __init__(self, 
                 start, 
                 end, 
                 sequence,
                 upstream_sequence,
                 pfam_2_coords):
        super(Ripp, self).__init__(start, 
                                     end, 
                                     sequence,
                                     upstream_sequence,
                                     pfam_2_coords)
        self.peptide_type = 'grasp'
        self.set_split()
        self.set_monoisotopic_mass()
        self.csv_columns = [self.leader, self.core, self.start, self.end]
        self.CUTOFF = CUTOFF
        
    def set_split(self):
        #TODO add regexes
        self.split_index = int(.5*len(self.sequence))
        
        self.leader = self.sequence[0:self.split_index]
        self.core = self.sequence[self.split_index:]
                
    def get_fimo_score(self, fimo_file):
        fimo_output = self.run_fimo_simple(fimo_file)
        fimo_motifs = []
        fimo_motifs = [line.partition("\t")[0] for line in fimo_output.split("\n")[1:] if "\t" in line]
        fimo_scores = defaultdict(int) # = {int(line.split("\t")[0]): float(line.split("\t")[6]) for line in fimo_output.split("\n") if "\t" in line and line.partition("\t")[0].isdigit()}
        for line in fimo_output.split("\n"):
            if not ("\t" in line and "motif_id" not in line):
                continue
            fimo_scores[line.split("\t")[0]] = float(line.split("\t")[6])
        #Calculate score
        motif_score = 0
        return fimo_motifs, fimo_scores
    
    def set_monoisotopic_mass(self):
#        print(self.core)
        if "B" in self.core:
            print(self.core)
            print("AA sequence contains Asx. Currently reviewing how to assign weights to such translations.")
        monoisotopic_mass = ProteinAnalysis(self.core.replace('X', '').replace('B', ''), monoisotopic=True).molecular_weight()
        self._monoisotopic_weight = monoisotopic_mass
    
    def compression_index(self, aa_list):
        positions = []
        for i, aa in enumerate(self.sequence):
            if aa in aa_list:
                positions.append(i)
        if len(positions) == 0:
            return 0
        return np.mean([np.mean(positions), np.std(positions), max(positions),\
                        min(positions), max(positions)-min(positions)])
            
        
    def set_score(self, pfam_dir, cust_hmm):
        scoring_csv_columns = [self.radar_score]
        self.score = 0
#        mdn_bcd = ["TIGR04185", "TIGR04184", "PF00583", "PF13673", "PF13302", "PF13523", "TIGR04187", "TgnB", "PsnB", "MvdD", "MvdC"]
        atp_grasp = ["TIGR04184", "TIGR04185", "TIGR04187"]
        atp_grasp_coords = []
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in atp_grasp):
                atp_grasp_coords += self.pfam_2_coords[pfam]
        min_distance_atp_grasp = self.get_min_dist(atp_grasp_coords)
        if min_distance_atp_grasp is None:
            min_distance_atp_grasp = 66666
        within_300 = False
        within_1100 = False
        within_2000 = False
        
        if min_distance_atp_grasp < 2000:
            within_2000 = True
            if min_distance_atp_grasp < 1100:
                within_1100 = True
                if min_distance_atp_grasp < 300:
                    within_300 = True
                    self.score += 2
                else:
                    self.score += 1
        else:
            self.score -= 1
                
        if within_1100:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if within_300:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if not within_2000:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        per_asp = self.core.count("D")/float(len(self.core))
        if per_asp > 0.08:
            scoring_csv_columns.append(1)
            self.score += 1
        else:
            scoring_csv_columns.append(0)
            
        per_thr = self.core.count("T")/float(len(self.core))
        if per_thr > 0.07:
            scoring_csv_columns.append(1)
            self.score += 1
        else:
            scoring_csv_columns.append(0)
            
        per_lys = self.sequence.count("K")/float(len(self.sequence))
        if per_lys > 0.07:
            scoring_csv_columns.append(1)
            self.score += 1
        else:
            scoring_csv_columns.append(0)
            
        #
        per_pro = self.core.count("P")/float(len(self.core))
        if per_pro > 0.07:
            scoring_csv_columns.append(1)
            self.score += 1
        else:
            scoring_csv_columns.append(0)
            
        if self.sequence.count("C") > 0:
            scoring_csv_columns.append(1)
            self.score += -2
        else:
            scoring_csv_columns.append(0)
            
        acceptor_core_count = self.core.count('E') + self.core.count('D')
        if acceptor_core_count < 3:
            self.score += -1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

        donor_core_count = self.core.count('S') + self.core.count('T') + self.core.count('L')
        if donor_core_count < 3:
            self.score += -1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
           
        if acceptor_core_count < donor_core_count:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        acceptor_leader_count = self.leader.count('E') + self.leader.count('D')
        if acceptor_core_count/float(len(self.core)) > acceptor_leader_count/float(len(self.leader)):
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        if acceptor_core_count/float(len(self.core)) > 0.13:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        if donor_core_count/float(len(self.core)) > 0.18:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        if self.sequence[-1] in "L":
            self.score += 2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        
        if self.sequence[-1] in "DVS":
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        
        match = re.search('(PF.L)', self.leader)
        if match is None:
            scoring_csv_columns.append(0)
        else:
            self.score += 1
            scoring_csv_columns.append(1)
            
        if self.start < self.end:
            direction = 1
        else:
            direction = -1
        same_dir = False
        for tag in atp_grasp:
            if tag in self.pfam_2_coords.keys():
                for coord in self.pfam_2_coords[tag]:
                    if (coord[0] < coord[1] and direction == 1) or coord[0] > coord[1] and direction == -1:
                        same_dir = True
        if same_dir:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
                
            
        ABC_trans = False
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in ["PF0005", "PF06472", "PF00664", "PF03412", 
                                           "TIGR03796", "TIGR01846", "TIGR03797", 
                                           "TIGR00954", "TIGR02203", "TIGR02204", 
                                           "TIGR03375"]):
                ABC_trans = True
        if ABC_trans:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        acetyl_t = False    
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in ["PF13302", "PF00583", "PF13523"]):
                acetyl_t = True
        if acetyl_t:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        methyl_t = False    
        for pfam in self.pfam_2_coords.keys():
            if "TIGR04188" in pfam:
                methyl_t = True
        if methyl_t:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        precursor_hmm_info = hmmer_utils.get_hmmer_info(self.sequence, pfam_dir, cust_hmm)
        pfams = []
        for pfam_dot, _, _, _, in precursor_hmm_info:
            pfams.append(pfam_dot.split('.')[0])
            
        if "PF12559" in pfams or "TIGR04186" in pfams:
            scoring_csv_columns.append(1)
            self.score += 5
        else:
            scoring_csv_columns.append(0)
        
        precursor_hit = False
        targets = ["PF14404", "PF14406", "PF14407", "PF14408", "PF14409"]
        for tar in targets:
            if tar in pfams:
                precursor_hit = True
        if precursor_hit in pfams:
            scoring_csv_columns.append(1)
            self.score += 5
        else:
            scoring_csv_columns.append(0) 
        donor_c_index = self.compression_index(['S', 'T', 'K'])
        acceptor_c_index = self.compression_index(['D','G'])
        if acceptor_c_index > donor_c_index:
            scoring_csv_columns.append(1)
            self.score += 1
        else:
            scoring_csv_columns.append(0)  
        
        
        charge_dict = {"E": -1, "D": -1, "K": 1, "H": 1, "R": 1}
        leader_charge = sum([charge_dict[aa] for aa in self.leader if aa in charge_dict])
        core_charge = sum([charge_dict[aa] for aa in self.core if aa in charge_dict])
        if leader_charge > core_charge:
            scoring_csv_columns.append(1)
            self.score += 1
        else:
            scoring_csv_columns.append(0) 
        
        #MOTIFS
        all_fimo_motifs, all_fimo_score = [], {}
        # switch the # on the two lines below if your computer username has a space in it. # RAL
        #fimo_motifs, fimo_scores = self.get_fimo_score("{}/MEME/Graspetide_1/meme.txt".format("./ripp_modules/grasp"))
        fimo_motifs, fimo_scores = self.get_fimo_score("{}/MEME/Graspetide_1/meme.txt".format(FILE_DIR))
        all_motifs = [
            'KNKNSKKKPFFASFLEKQVKDPETVKGG',
             'VTMKYPSDGDE',
             'DITIPERDVITKPIIDTVTSP',
             'KAVPFFARFLEEQKEEESSEE',
             'VTTLAIGEEGG',
             'MKNHPLTQNQPANAKPVTLAESQLFEVSGGLRARLSESVPV',
             'MHPFSLEKAQLENVVGGVQEV',
             'LPDILS',
             'STTITETIENSDEDE',
             'MKPYIMNYSETIEJS',
             'SLAPSNCERKYEISNEKQSAM',
             'NLVYDHITQLSYJDKEKTVKAINFGPR',
             'YGZTLITKTREGIDRSERS',
             'NQKSPLLKYRIECDCNETAELPNKSYDFVSGLWVDESGEP',
             'MKQKNPILYFR',
             'ANHLLELVTKTKVQQERDDDS',
             'ASPFLARFATPRTGEDNJPGYYSYEQDMWVVETPQGEVPII',
             'LELLTKTYAKVESDD',
             'NDPKPFILNYATERVHSNIPTLGGHYCDDEHMWVDDESNTPIINLVDKNC',
             'LLMGTETFTKTYNEGSDSDNB',
             'MKPLILQFTEKPTEENIDYSQIEYSKKLN',
             'SVDKKTGLPAI']
        motif_1a = "KNKNSKKKPFFASFLEKQVKDPETVKGG"
        motif_1b = "KAVPFFARFLEEQKEEESSEE"
        motif_2 = "VTMKYPSDGDE"
        motif_3 = "?"
        motif_4 = "VTTLAIGEEGG"
        motif_5 = "MKPYIMNYSETIEJS"
        motif_6 = "STTITETIENSDEDE"
        motif_7 = "YGZTLITKTREGIDRSERS"
        motif_8 = "ANHLLELVTKTKVQQERDDDS"
        motif_9 = "LLMGTETFTKTYNEGSDSDNB"
        if motif_1a in fimo_motifs or motif_1b in fimo_motifs:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if motif_2 in fimo_motifs:
            self.score += fimo_motifs.count(motif_2)
            scoring_csv_columns.append(fimo_motifs.count(motif_2))
        else:
            scoring_csv_columns.append(0)
        # TODO MISSING MEME 3
        scoring_csv_columns.append(0)
        all_fimo_motifs.extend(fimo_motifs)
        all_fimo_score.update(fimo_scores)
        # switch the # on the two lines below if your computer username has a space in it. # RAL
        #fimo_motifs, fimo_scores = self.get_fimo_score("{}/MEME/Graspetide_2/meme.txt".format("./ripp_modules/grasp"))
        fimo_motifs, fimo_scores = self.get_fimo_score("{}/MEME/Graspetide_2/meme.txt".format(FILE_DIR))
        if motif_4 in fimo_motifs:
            self.score += fimo_motifs.count(motif_4)
            scoring_csv_columns.append(fimo_motifs.count(motif_4))
        else:
            scoring_csv_columns.append(0)
        all_fimo_motifs.extend(fimo_motifs)
        all_fimo_score.update(fimo_scores)
        # switch the # on the two lines below if your computer username has a space in it. # RAL
        #fimo_motifs, fimo_scores = self.get_fimo_score("{}/MEME/Graspetide_3/meme.txt".format("./ripp_modules/grasp"))
        fimo_motifs, fimo_scores = self.get_fimo_score("{}/MEME/Graspetide_3/meme.txt".format(FILE_DIR))
        if motif_5 in fimo_motifs:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if motif_6 in fimo_motifs:
            self.score += fimo_motifs.count(motif_6)
            scoring_csv_columns.append(fimo_motifs.count(motif_6))
        else:
            scoring_csv_columns.append(0)
        all_fimo_motifs.extend(fimo_motifs)
        all_fimo_score.update(fimo_scores)
         # switch the # on the two lines below if your computer username has a space in it. # RAL
        #fimo_motifs, fimo_scores = self.get_fimo_score("{}/MEME/Graspetide_4/meme.txt".format("./ripp_modules/grasp"))   
        fimo_motifs, fimo_scores = self.get_fimo_score("{}/MEME/Graspetide_4/meme.txt".format(FILE_DIR))
        if motif_7 in fimo_motifs:
            self.score += fimo_motifs.count(motif_7)
            scoring_csv_columns.append(fimo_motifs.count(motif_7))
        else:
            scoring_csv_columns.append(0)
        all_fimo_motifs.extend(fimo_motifs)
        all_fimo_score.update(fimo_scores)
        # switch the # on the two lines below if your computer username has a space in it. # RAL
        #fimo_motifs, fimo_scores = self.get_fimo_score("{}/MEME/Graspetide_5/meme.txt".format("./ripp_modules/grasp"))   
        fimo_motifs, fimo_scores = self.get_fimo_score("{}/MEME/Graspetide_5/meme.txt".format(FILE_DIR))
        if motif_8 in fimo_motifs:
            self.score += fimo_motifs.count(motif_8)
            scoring_csv_columns.append(fimo_motifs.count(motif_8))
        else:
            scoring_csv_columns.append(0)
        all_fimo_motifs.extend(fimo_motifs)
        all_fimo_score.update(fimo_scores)
        # switch the # on the two lines below if your computer username has a space in it. # RAL
        #fimo_motifs, fimo_scores = self.get_fimo_score("{}/MEME/Graspetide_6/meme.txt".format("./ripp_modules/grasp"))   
        fimo_motifs, fimo_scores = self.get_fimo_score("{}/MEME/Graspetide_6/meme.txt".format(FILE_DIR))
        if motif_9 in fimo_motifs:
            self.score += fimo_motifs.count(motif_9)
            scoring_csv_columns.append(fimo_motifs.count(motif_9))
        else:
            scoring_csv_columns.append(0)
        all_fimo_motifs.extend(fimo_motifs)
        all_fimo_score.update(fimo_scores)
        
        ##SVM SCORING SECTION
        scoring_csv_columns.append(min_distance_atp_grasp)
        scoring_csv_columns.append(len(self.sequence))
        scoring_csv_columns.append(leader_charge)
        scoring_csv_columns.append(core_charge)
        scoring_csv_columns.append(leader_charge + core_charge)
        scoring_csv_columns.append(abs(core_charge))
        scoring_csv_columns.append(abs(leader_charge))
        scoring_csv_columns.append(abs(leader_charge + core_charge))
        #Counts of AAs in leader
        scoring_csv_columns += [self.leader.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
        #Aromatics in leader
        scoring_csv_columns.append(sum([self.leader.count(aa) for aa in "FWY"]))
        #Neg charged in leader
        scoring_csv_columns.append(sum([self.leader.count(aa) for aa in "DE"]))
        #Pos charged in leader
        scoring_csv_columns.append(sum([self.leader.count(aa) for aa in "RK"]))
        #Charged in leader
        scoring_csv_columns.append(sum([self.leader.count(aa) for aa in "RKDE"]))
        #Aliphatic in leader
        scoring_csv_columns.append(sum([self.leader.count(aa) for aa in "GAVLMI"]))
        #Hydroxyl in leader
        scoring_csv_columns.append(sum([self.leader.count(aa) for aa in "ST"]))
        #Counts of AAs in core
        scoring_csv_columns += [self.core.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
        #Aromatics in core
        scoring_csv_columns.append(sum([self.core.count(aa) for aa in "FWY"]))
        #Neg charged in core
        scoring_csv_columns.append(sum([self.core.count(aa) for aa in "DE"]))
        #Pos charged in core
        scoring_csv_columns.append(sum([self.core.count(aa) for aa in "RK"]))
        #Charged in core
        scoring_csv_columns.append(sum([self.core.count(aa) for aa in "RKDE"]))
        #Aliphatic in core
        scoring_csv_columns.append(sum([self.core.count(aa) for aa in "GAVLMI"]))
        #Hydroxyl in core
        scoring_csv_columns.append(sum([self.core.count(aa) for aa in "ST"]))
        
        #Counts of AAs in leader+core
        scoring_csv_columns += [self.sequence.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"] #Temp to work with current training CSV
        #Aromatics in precursor
        scoring_csv_columns.append(sum([self.sequence.count(aa) for aa in "FWY"]))
        #Neg charged in precursor
        scoring_csv_columns.append(sum([self.sequence.count(aa) for aa in "DE"]))
        #Pos charged in precursor
        scoring_csv_columns.append(sum([self.sequence.count(aa) for aa in "RK"]))
        #Charged in precursor
        scoring_csv_columns.append(sum([self.sequence.count(aa) for aa in "RKDE"]))
        #Aliphatic in precursor
        scoring_csv_columns.append(sum([self.sequence.count(aa) for aa in "GAVLMI"]))
        #Hydroxyl in precursor
        scoring_csv_columns.append(sum([self.sequence.count(aa) for aa in "ST"]))
        #Counts (0 or 1) of amino acids within last AA position of sequence
        scoring_csv_columns += [self.sequence[-1].count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
        
        #acceptor compression index TODO
        scoring_csv_columns.append(acceptor_c_index)
        #donor compression index TODO
        scoring_csv_columns.append(donor_c_index)
        
        self.fimo_motifs = all_fimo_motifs
        self.fimo_scores = all_fimo_score
        #Total motifs hit
        scoring_csv_columns.append(len(all_fimo_motifs))
        #Motif scores
        scoring_csv_columns += [all_fimo_score[motif] if motif in all_fimo_score else 0 for motif in all_motifs]
        
        #Sum of MEME scores
        scoring_csv_columns.append(sum([all_fimo_score[motif] if motif in all_fimo_score else 0 for motif in all_motifs]))
        #No Motifs?
        if len(all_fimo_motifs) == 0:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if self.leader[0] != 'M':
            self.score += -1
        self.csv_columns += [self.score] +  scoring_csv_columns

