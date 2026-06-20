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
CUTOFF = 20
index = 0


def write_csv_headers(output_dir, meta=False):
    dir_prefix = output_dir + '/grasp/'
    if not os.path.exists(dir_prefix):
        os.makedirs(dir_prefix)
    svm_headers = 'Precursor Index,classification, Radar Score, <1100 from ATP-grasp ligase, <300 from ATP-grasp ligase, >2000 from ATP-grasp ligase, 2nd half contains >8% Asp residue (acceptor site), 2n half of precursor contains >8% Thr residues (donor site), Second half of precursor contains >7% Pro residues, Second half of precursor contains <3 acceptor residues (Asp + Glu), Second half of precursor contains <3 donor residues (Ser + Thr + Lys), Second half of precursor contains more donor residues than acceptor residues, % acceptor residues in second half of precursor > % acceptor residues in first half of precursor,  % acceptor residues (Asp + Glu) in second half of precursor >13%, % donor residues (Ser + Thr + Lys) in second half of precursor >18%, Precursor ends with Asp, First half of the precursor contains a “PFxL” motif, Precursor and the ATP-grasp protein (MdnC homolog) are encoded on same strand, Gene cluster contains one of the following: PF0005 PF06472 PF00664 PF03412 TIGR03796 TIGR01846 TIGR03797 TIGR00954 TIGR02203 TIGR02204 TIGR03375 (ABC transporters that co-occur frequently), Gene cluster contains one of the following: PF13302 PF00583 PF13523 (acetyltransferases that co-occur frequently), A local gene product hits one of the following: TIGR04188 TIGR04364 PF01135 TIGR00080 (methyltransferase), Precursor hits PF12559 (serine endopeptidase inhibitor) or TIGR04186 (GRASP_targ), Precursor hits PF14404 (strep_pep) PF14406 (bacteroid_pep) PF14407 (frankia_pep) PF14408 (actino_pep) or PF14409 (herpeto_pep), Acceptor compressionC index > Donor compressionC index*, Calculated charge at pH7 of second half of precursor < charge of first half of precursor, Acceptor regex match, Donor regex match, Both acceptor & donor regex match, Precursor peptide contains sequence motif #1 “Group1 leader”, motif #2 "Group1 core",motif #3 "Group2 leader",motif #4 "Group2 core",motif #5 "Group3 leader",motif #6 "Group3 core",motif #7 "Group4 core",motif #8 "Group5 core",motif #9 "Group6 core",motif #10 "Group7 core",motif #11 "Group11 core",motif #12 "Group13 core",motif #13 "Group16 core",motif #14 "Group21 core",Count of sequence motif #1 “Group1 leader”, motif #2 "Group1 core",motif #3 "Group2 leader",motif #4 "Group2 core",motif #5 "Group3 leader",motif #6 "Group3 core",motif #7 "Group4 core",motif #8 "Group5 core",motif #9 "Group6 core",motif #10 "Group7 core",motif #11 "Group11 core",motif #12 "Group13 core",motif #13 "Group16 core",motif #14 "Group21 core",Total motifs in sequence, No motifs,Minimum distance from MdnB C D TgnB and PsnB homologs  (nt), Precursor length, Estimated first half charge, Estimated second half charge, Estimated precursor charge,Absolute value of second half charge,Absolute value of first half charge,Absolute value of precursor charge,FIRST HALF A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,Aromatics,Neg charged,Pos charged,Charged,Aliphatic,Hydroxyl,SECOND HALF A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,Aromatics,Neg charged,Pos charged,Charged,Aliphatic,Hydroxyl,PRECURSOR A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,Aromatics,Neg charged,Pos charged,Charged,Aliphatic,Hydroxyl, acceptor compression index, donor compression index'
    svm_headers = svm_headers.split(',')
    if meta: 
        features_headers = ['Accession_id', 'Locus', 'Genus/Species', 'Nucleotide_Acc', 'Precursor peptide','First half', 'Second Half','Start', 'End', 'Best graspetide synthetase', 'Best graspetide synthetase evalue', 'Multiple graspetide synthetases', "Total Score", "Valid Precursor"] + svm_headers 
    else:
        features_headers = ['Accession_id', 'Genus/Species', 'Precursor peptide','First half', 'Second Half','Start', 'End', 'Best graspetide synthetase', 'Best graspetide synthetase evalue', 'Multiple graspetide synthetases', "Total Score", "Valid Precursor"] + svm_headers 
    features_csv_file = open(dir_prefix + "temp_features.csv", 'w')
    svm_csv_file = open("{}fitting_set.csv".format(dir_prefix), 'w')
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
                 pfam_2_coords,
                 output_dir,
                 pfam_2_evalue):
        super(Ripp, self).__init__(start, 
                                     end, 
                                     sequence,
                                     upstream_sequence,
                                     pfam_2_coords,
                                     output_dir,
                                     pfam_2_evalue)
        self.peptide_type = 'grasp'
        self.set_split()
        self.set_monoisotopic_mass()
        self.best_grasp_acc()
        self.csv_columns = [self.sequence, self.leader, self.core, self.start, self.end, self.best_grasp_syn, self.best_grasp_syn_evalue, self.gra_syn_list]
        self.CUTOFF = CUTOFF

        #pfam_2_evalue has a format of: [protein accession, HMM/PFAM hit, e value, start, end, AA length]
        
    def set_split(self):
        self.split_index = int(.5*len(self.sequence))
        
        self.leader = self.sequence[0:self.split_index]
        self.core = self.sequence[self.split_index:]
    
    def best_grasp_acc(self):
# Best scoring graspetide synthetase
        best_grasp_syn = 0
        best_grasp_syn_evalue = 1
        for gra in self.pfam_2_evalue:
            if gra[1] == "Graspetide_synthetase" and gra[2] < best_grasp_syn_evalue:
                best_grasp_syn = gra[0]
                best_grasp_syn_evalue = gra[2]
        self.best_grasp_syn = best_grasp_syn
        self.best_grasp_syn_evalue = best_grasp_syn_evalue
        
# Multiple graspetide synthetases in BGC
        multi_gra_syn = []
        gra_syn_list = ""
        i = 0
        for gra in self.pfam_2_evalue:
            if gra[1] == "Graspetide_synthetase" and gra[0] not in multi_gra_syn:
                multi_gra_syn.append(gra[0])
                i += 1
        if i >= 2:
            gra_syn_list = str(multi_gra_syn)
            gra_syn_list = gra_syn_list.replace(",", " ")
        self.gra_syn_list = gra_syn_list 
        
        
        
    def get_fimo_score(self): 
        fimo_output = self.run_fimo_simple()
        fimo_motifs = []
        fimo_motifs = [line.partition("\t")[0] for line in fimo_output.split("\n")[1:] if "\t" in line]
        fimo_scores = defaultdict(int) # = {int(line.split("\t")[0]): float(line.split("\t")[6]) for line in fimo_output.split("\n") if "\t" in line and line.partition("\t")[0].isdigit()}
        for line in fimo_output.split("\n"):
            if not ("\t" in line and "motif_id" not in line):
                continue
            fimo_scores[line.split("\t")[0]] = float(line.split("\t")[6])
        return fimo_motifs, fimo_scores
    
    def set_monoisotopic_mass(self):
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
        scoring_csv_columns = []
        scoring_csv_columns = [self.radar_score]
        self.score = 0
#        mdn_bcd = ["TIGR04185", "TIGR04184", "PF00583", "PF13673", "PF13302", "PF13523", "TIGR04187", "TgnB", "PsnB", "MvdD", "MvdC"]
        atp_grasp = ["Graspetide_synthetase", "TIGR04187", "TIGR04192"]
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
        if per_thr > 0.08:
            scoring_csv_columns.append(1)
            self.score += 1
        else:
            scoring_csv_columns.append(0)
        
        per_pro = self.core.count("P")/float(len(self.core))
        if per_pro > 0.07:
            scoring_csv_columns.append(1)
            self.score += 1
        else:
            scoring_csv_columns.append(0)
            
        acceptor_core_count = self.core.count('E') + self.core.count('D')
        if acceptor_core_count < 3:
            self.score += -1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

        donor_core_count = self.core.count('S') + self.core.count('T') + self.core.count('K') 
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
        
        if self.sequence[-1] in "D":
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
            self.score += 2
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
            if any(fam in pfam for fam in ["TIGR04188", "TIGR04364", "PF01135", "TIGR00080"]):
                methyl_t = True
        if methyl_t:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        precursor_hmm_info = hmmer_utils.get_hmmer_info(self.sequence, pfam_dir, cust_hmm, self.output_dir)
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
        if precursor_hit:
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
        
        
        #Acceptor Regex
        acc_match = re.search('([DE].{0,3}[DE].{0,3}[DE])', self.sequence)
        if acc_match is None:
            scoring_csv_columns.append(0)
        else:
            scoring_csv_columns.append(1)
        #Donor Regex
        don_match = re.search('([TSK].{0,3}[TSK].{0,3}[TSK])', self.sequence)
        if don_match is None:
            scoring_csv_columns.append(0)
        else:
            scoring_csv_columns.append(1)
        #both regex
        if acc_match and don_match:
            scoring_csv_columns.append(1)
            self.score += 3
        else:
            scoring_csv_columns.append(0)
            
        #FIMO motifs
        fimo_motifs, fimo_scores = self.get_fimo_score()
        self.fimo_motifs = sorted(fimo_motifs)
        self.fimo_scores = sorted(fimo_scores)
                   
        # Is the motif present in the sequence 
        j = 1        
        for j in range(1, 15):
            mot_scor1 = 0
            for motif1 in fimo_motifs:
                motif_split1 = motif1.split("-")
                if int(motif_split1[0]) == j:
                    mot_scor1 = 1 
            scoring_csv_columns.append(mot_scor1)

        # Number of each motif in sequence              
        k = 1        
        for k in range(1, 15):
            mot_scor2 = 0
            for motif2 in fimo_motifs:
                motif_split2 = motif2.split("-")
                if int(motif_split2[0]) == k:
                    mot_scor2 += 1
            # add 1 point if a leader motif found
            if (k == 1 or 3 or 5) and (mot_scor2 > 0):
                self.score += 1
            # add a point for each core motif found
            else:
                self.score += mot_scor2
            scoring_csv_columns.append(mot_scor2) 
        
        # How many total motifs are present in the sequence
        scoring_csv_columns.append(len(fimo_motifs))
        
        # No Motifs
        if len(fimo_motifs) == 0:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

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
        
        #acceptor compression index TODO
        scoring_csv_columns.append(acceptor_c_index)
        #donor compression index TODO
        scoring_csv_columns.append(donor_c_index)
        
        self.csv_columns += [self.score] +  scoring_csv_columns

