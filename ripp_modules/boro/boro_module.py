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
#
#==============================================================================
# Special thanks goes to AntiSmash team, whose antiSMASH-rodeo repository at
# https://bitbucket.org/mmedema/antismash-rodeo/ provided the backbone code for 
# a great deal of the heuristic calculations.
#==============================================================================

import csv
import os
from pickle import FALSE
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from ripp_modules.VirtualRipp import VirtualRipp
import hmmer_utils
from collections import defaultdict
import pathlib
FILE_DIR = pathlib.Path(__file__).parent.absolute()
peptide_type = "boro"
CUTOFF = 14 #TODO change cutoff to match svm
index = 0


#sets up the csv file with correct headers

def write_csv_headers(output_dir):
    dir_prefix = output_dir + '/boro/'
    if not os.path.exists(dir_prefix):
        os.makedirs(dir_prefix)
    svm_headers = 'Precursor Index, Classification, Precursor is a Methyltransferase, <1000 nt from mtase, <200 nt from mtase, <100 nt from mtase, >3000 nt from mtase, mtase (PF00590) present, mtase hits BoroMT hmm, mtase hits NMT_1(DCLFAD) hmm, mtase hits NMT_2(YGHP) hmm, LigA-like domain (PF07746) present, nearby domain hits a BBD hmm, precursor & mtase same direction, mtase & BBD both present, GGDEF present, acetyltransferase present, peptidase present, precursor has mtase domain (PF00590), precursor hits a borosin mtase hmm, precursor has a LigA-like domain (PF07746), precursor hits BBD_A hmm, precursor hits BBD_B hmm, precursor hits a legionellales hmm, precursor hits type VI/VIII hmm, precursor hits mtase & BBD hmm, precursor hits a legionellales & BBD hmm, precursor hits type VI/VIII & BBD hmm, precursor hits a legionellales & mtase hmm, precursor hits type VI/VII & mtase hmm, precursor is > 900 AA & hits a BBD hmm, precursor is > 700 AA & hits a BBD hmm, precursor is < 200 AA & hits a BBD hmm, precursor is < 400 AA & hits a BBD hmm, FIMO motif 1 present, motif 2, motif 3, motif 4, motif 5, motif 6, motif 7, motif 8, motif 9, motif 10, motif 11, motif 12, motif 13, motif 14, motif 15, motif 16, motif 17, motif 18, motif 19, motif 20, motif 21, motif 22, motif 23, motif 24, motif 25, motif 26, motif 27, Number of motifs present, Occurances of motif 1, motif 2, motif 3, motif 4, motif 5, motif 6, motif 7, motif 8, motif 9, motif 10, motif 11, motif 12, motif 13, motif 14, motif 15, motif 16, motif 17, motif 18, motif 19, motif 20, motif 21, motif 22, motif 23, motif 24, motif 25, motif 26, motif 27, No Motifs, mtase distance, precursor length, core charge, precursor charge, abs(core charge), abs(precursor charge), CORE COUNT A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V, PERCENT A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V, aromatics, neg, pos, charged, aliphatic, hydroxyl, isoelectric point, PRECURSOR count A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V, PERCENT A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V, aromatics, neg, pos, charged, aliphatic, hydroxyl, isoelectric point'
    svm_headers = svm_headers.split(',')
    features_headers = ['Accession_id', 'Genus/Species', 'Sequence', 'Region1', 'Region2', 'Region3', 'Start', 'End', 'Total Score',"Valid Precursor",] + svm_headers
#TODO close all these write headers at the end of the document
    features_csv_file = open(dir_prefix + "temp_features.csv", 'w')
    svm_csv_file = open("{}/svm/fitting_set.csv".format(FILE_DIR), 'w')
    features_writer = csv.writer(features_csv_file)
    svm_writer = csv.writer(svm_csv_file)
    features_writer.writerow(features_headers)
    svm_writer.writerow(svm_headers)#Don't include accession_id, genus/species,leader, core sequence, score, or svm classification
                                        
class Ripp(VirtualRipp):
# sets up borosins as a class and separates out the final 120 residues or middle 120 residues depending on peptide length

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
        
        self.peptide_type = 'boro'
        self.set_split()
        self.csv_columns = [self.sequence, self.region1, self.region2, self.region3, self.start, self.end]
        self.CUTOFF = CUTOFF

    
    def set_split(self):
# splits the precursor into parts based on length. Region 2 is evaluated as the potential "core", although most of the score evaluates the entire precursor.    
        if int(len(self.sequence)) < 120:
            self.region1 = "N/A"
            self.region2 = self.sequence
            self.region3 = "N/A"
        elif int(len(self.sequence)) <= 400:
            self.split_index = int(len(self.sequence)-120) 
        
            self.region1 = self.sequence[0:self.split_index]
            self.region2 = self.sequence[self.split_index:]
            self.region3 = "N/A"
        else:
            self.split_index1 = int((.5*len(self.sequence))-60)
            self.split_index2 = int((.5*len(self.sequence))+60)  
            self.region1 = self.sequence[0:self.split_index1]
            self.region2 = self.sequence[self.split_index1:self.split_index2]
            self.region3 = self.sequence[self.split_index2:]
        self.core = self.region2

# gets the fimo score from the motifs in the fimo file (fimo file is a STREME output)
    def get_fimo_score(self):   
        fimo_output = self.run_fimo_simple()
        fimo_motifs = [line.partition("\t")[0] for line in fimo_output.split("\n")[1:] if "\t" in line]
        fimo_scores = defaultdict(int) # = {int(line.split("\t")[0]): float(line.split("\t")[6]) for line in fimo_output.split("\n") if "\t" in line and line.partition("\t")[0].isdigit()}
        for line in fimo_output.split("\n")[1:]:
            if not ("\t" in line and "motif_id" not in line):
                continue
            fimo_scores[line.split("\t")[0]] = float(line.split("\t")[6])
 # if the precursor has a 
        fimo_motif_score = 0
        fimo_motif_count = 0
        for hits in fimo_motifs:
            fimo_motif_count += 1 
        if fimo_motif_count == 0:
            fimo_motif_score = -3
        elif fimo_motif_count == 1:
            fimo_motif_score = 1
        elif fimo_motif_count == 2:
            fimo_motif_score = 2
        elif fimo_motif_count >= 4:
            fimo_motif_score = 3
        return fimo_motifs, fimo_motif_score, fimo_scores

# NEIGHBORHOOD SCORING SECTION
    def set_score(self, pfam_dir, cust_hmm): 
        scoring_csv_columns = []
        self.score = 0

# stores pfams that precursor hits 
        precursor_hmm_info = hmmer_utils.get_hmmer_info(self.sequence, pfam_dir, cust_hmm)      
        pfams = []
        for pfam_dot, _, _, _, in precursor_hmm_info:
            pfams.append(pfam_dot.split('.')[0])

# Is the precursor a methyltransferase
        methyltransferase = False
        targets = ["PF00590", "NMT_1", "NMT_2", "BoroMT"]
        for tar in targets:
            if tar in pfams:
                methyltransferase = True
        if methyltransferase:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)  

            
# adjusts score based on precursor distance from mtase
        mtase = ["PF00590"] 
        mtase_coords = []
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in mtase):
                mtase_coords += self.pfam_2_coords[pfam]
        min_distance_mtase = self.get_min_dist(mtase_coords)
        if min_distance_mtase is None:
            min_distance_mtase = 66666 
        within_100 = False
        within_200 = False
        within_1000 = False
        within_3000 = False
# 3 points if within 100 residues, 2 points if within 200, 1 point if within 1000, -4 points if NOT within 3000
        if min_distance_mtase < 3000:
            within_3000 = True
            if min_distance_mtase < 1000:
                within_1000 = True
                if min_distance_mtase < 200:
                    within_200 = True
                    if min_distance_mtase < 100:
                        within_100 = True
                        self.score += 3
                    else:
                        self.score += 2
                else: 
                    self.score += 1
        else:
            self.score += -4
                
        if within_1000:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if within_200:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if within_100:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if not within_3000:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# if anything nearby has an mtase
        mtPfam = False
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in ["PF00590"]):
                mtPfam = True
        if mtPfam:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)   

# if the nearby mtase matches one or more of the hmms (additive)
        mtHmm = False
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in ["BoroMT"]):
                mtHmm = True
        if mtHmm:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

        dclfad = False
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in ["NMT_1"]):
                dclfad = True
        if dclfad:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

        yghp = False
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in ["NMT_2"]):
                yghp = True
        if yghp:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# if anything nearby has PF07746 nearby
        bbd = False
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in ["PF07746"]):
                bbd = True
        if bbd:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# if the nearby BBD actually matches one or both of the hmms (not additive)        
        BBDhmm = False
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in ["BBD_A", "BBD_B"]):
                BBDhmm = True
        if BBDhmm:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
       
# if the precursor is encoded on the same strand as the mtase
        if self.start < self.end:
            direction = 1
        else:
            direction = -1        
        same_dir = False
        for tag in mtase:
            if tag in self.pfam_2_coords.keys():
                for coord in self.pfam_2_coords[tag]:
                    if (coord[0] < coord[1] and direction == 1) or coord[0] > coord[1] and direction == -1:
                        same_dir = True
        if same_dir:
            self.score += 2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
# if BBD and mtPFAM both occur in BGC
        BBD_PFAM = False
        if BBDhmm and mtPfam == True:
            BBD_PFAM = True
        if BBD_PFAM:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
# if ggdef co-occurs       
        ggdef = False    
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in ["PF00990"]):
                ggdef = True
        if ggdef:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

#if an acetyltransferase co-occurs        
        acetyl_t = False    
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in ["PF13302", "PF00583"]):
                acetyl_t = True
        if acetyl_t:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# if a peptidase of interest co-occurs
        peptidase = False
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in ["PF06167", "PF13650", "PF04389", "PF00246", "PF02225", "PF00768", "PF18027", "PF13529", "PF00026", "PF05649", "PF01431"]):
                peptidase = True
        if peptidase:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# PRECURSOR SCORING SECTION  
# if the precursor has a mt domain    
        if "PF00590" in pfams:
            scoring_csv_columns.append(1)
            self.score += -2
        else:
            scoring_csv_columns.append(0)

# if the precursor hits any of the borosin mt hmms
        precursor_mt_hit = False
        targets = ["NMT_1", "NMT_2", "BoroMT"]
        for tar in targets:
            if tar in pfams:
                precursor_mt_hit = True
        if precursor_mt_hit:
            scoring_csv_columns.append(1)
            self.score += -2
        else:
            scoring_csv_columns.append(0)  

# if the precursor has a LigA-like domain or hits a BBD hmm
        precursor_bbd_hit = False
        targets = ["PF07746", "BBD_A", "BBD_B"]
        for tar in targets:
            if tar in pfams:
                precursor_bbd_hit = True
                self.score += 1

# if the precursor has a LigA-like domain
        if "PF07746" in pfams:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# if the precursor hits the first BBD hmm
        if "BBD_A" in pfams:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# if the precursor hits the second BBD hmm
        if "BBD_B" in pfams:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)   


# if the precursor hits either of the legionellales hmms
        precursor_leg_hit = False
        targets = ["BoroLegA", "BoroLegB"]
        for tar in targets:
            if tar in pfams:
                precursor_leg_hit = True
        if precursor_leg_hit == True: 
            scoring_csv_columns.append(1)
            self.score += 3
        else:
            scoring_csv_columns.append(0) 

# if the precursor hits the type VI/VIII hmm
        precursor_6_8_hit = False
        if "Boro_6_8" in pfams:
            scoring_csv_columns.append(1)
            self.score += 3
            precursor_6_8_hit = True
        else:
            scoring_csv_columns.append(0)

# if precursor hits mtase & BBD hmm
        if precursor_mt_hit and precursor_bbd_hit:
             scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)                                   
                                 
# if precursor hits legionellales & BBD hmm
        if precursor_leg_hit and precursor_bbd_hit:
             scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)                                 
                                 
# if precursor hits type VI/VIII & BBD hmm
        if precursor_6_8_hit and precursor_bbd_hit:
             scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)                                 
#                                   
# if precursor hits legionellales & mtase hmm
        if precursor_leg_hit and precursor_mt_hit:
             scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)                                 
                                 
# if precursor hits type VI/VIII & mtase hmm
        if precursor_6_8_hit and precursor_mt_hit:
             scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)         

# if precursor hits is >900AA & hits BBD hmm
        if int(len(self.sequence)) >= 900 and precursor_bbd_hit:
             scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0) 

# if precursor hits is >700AA & hits BBD hmm
        if int(len(self.sequence)) >= 700 and precursor_bbd_hit:
             scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0) 
            
# if precursor hits is <400AA & hits BBD hmm
        if int(len(self.sequence)) <= 400 and precursor_bbd_hit:
             scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# if precursor hits is <200AA & hits BBD hmm
        if int(len(self.sequence)) <= 200 and precursor_bbd_hit:
             scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
                    
# AA % of "core" (region 2) for DEIVT
        per_asp = self.core.count("D")/float(len(self.core))
        if per_asp > 0.10: 
            self.score += 1
        per_glu =  self.core.count("E")/float(len(self.core))   
        if per_glu > 0.10:
            self.score += 1             
        per_ile =  self.core.count("I")/float(len(self.core))   
        if per_ile > 0.10:
            self.score += 1         
        per_val =  self.core.count("V")/float(len(self.core))   
        if per_val > 0.10:
            self.score += 1 
        per_thr = self.core.count("T")/float(len(self.core))
        if per_thr > 0.10:
            self.score += 1

#fimo section
        fimo_motifs, motif_score, fimo_scores = self.get_fimo_score()
        self.fimo_motifs = sorted(fimo_motifs)
        self.fimo_scores = sorted(fimo_scores)
        self.score += motif_score
                   
# Is the motif present in the sequence 
        j = 1        
        for j in range(1, 28): 
            mot_scor1 = 0
            for motif1 in fimo_motifs:
                motif_split1 = motif1.split("-")
                if int(motif_split1[0]) == j:
                    mot_scor1 = 1 
            scoring_csv_columns.append(mot_scor1)

# How many total motifs are present in the sequence
        scoring_csv_columns.append(len(fimo_motifs))  

# Number of each motif in sequence              
        k = 1        
        for k in range(1, 28):
            mot_scor2 = 0
            for motif2 in fimo_motifs:
                motif_split2 = motif2.split("-")
                if int(motif_split2[0]) == k:
                    mot_scor2 += 1
            scoring_csv_columns.append(mot_scor2) 

# No Motifs
        if len(fimo_motifs) == 0:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)


# SVM SCORING SECTION
        scoring_csv_columns.append(min_distance_mtase)
        scoring_csv_columns.append(len(self.sequence))

        charge_dict = {"E": -1, "D": -1, "K": 1, "H": 1, "R": 1}
        precursor_charge = sum([charge_dict[aa] for aa in self.sequence if aa in charge_dict]) 
        core_charge = sum([charge_dict[aa] for aa in self.core if aa in charge_dict])
        scoring_csv_columns.append(core_charge)
        scoring_csv_columns.append(precursor_charge)
        scoring_csv_columns.append(abs(core_charge))
        scoring_csv_columns.append(abs(precursor_charge))

        #Counts of AAs in core
        scoring_csv_columns += [self.core.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
        #Percentages of AAs in core
        core_analysis = ProteinAnalysis(self.core, monoisotopic=True) 
        for aa in "ARDNCQEGHILKMFPSTWYV":
            scoring_csv_columns.append(core_analysis.get_amino_acids_percent()[aa])
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
        #isoelectric point in core
        scoring_csv_columns.append(core_analysis.isoelectric_point())
        
        #Counts of AAs in precursor
        scoring_csv_columns += [self.sequence.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
        #Percentages of AAs in precursor
        seq_analysis = ProteinAnalysis(self.sequence, monoisotopic=True)
        for aa in "ARDNCQEGHILKMFPSTWYV":
            scoring_csv_columns.append(seq_analysis.get_amino_acids_percent()[aa])
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
        #isoelectric point in precursor
        scoring_csv_columns.append(seq_analysis.isoelectric_point())

        self.csv_columns += [self.score] +  scoring_csv_columns
 

        
     


