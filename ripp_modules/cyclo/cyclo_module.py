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
import math
FILE_DIR = pathlib.Path(__file__).parent.absolute()

peptide_type = "cyclo"
CUTOFF = 16
index = 0


def write_csv_headers(output_dir, meta=False):
    dir_prefix = output_dir + '/cyclo/'
    if not os.path.exists(dir_prefix):
        os.makedirs(dir_prefix)
    svm_headers = "Precursor index,classification, rSAM present,rSAM distance <2000,rSAM distance <1000,rSAM distance >5000,ansme evalue,ansme_score,ansme score >40,ansme score >100,sulfatase present,HexxH domain present,precuror retrieved by Wpr RRE,precursor length <70,precursor length <90,precursor length >110,precursor length >200,1 cys,2 cys,3+ cys,>2 cyclophane motifs,Number of cyclophane motifs,Percent Aromatic,Triceptide precursor hmms,precursor no hmm match,rgg/shp,gene synteny,common transporter,rsam hmm,Inferred gene,fimo motif 1, motif 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,number of motifs present, occurances of 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,rSAM distance,precursor length,core charge,precursor charge,last_third charge,abs(core charge),abs(precursor charge),abs(last third charge),CORE COUNT A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,PERCENT A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,aromatics,neg,pos,charged,aliphatic,hydroxyl,isoelectric point,PRECURSOR count A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,PERCENT A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,aromatics,neg,pos,charged,aliphatic,hydroxyl,isoelectric point,LAST THIRD COUNT A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,PERCENT A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,aromatics,neg,pos,charged,aliphatic,hydroxyl,isoelectric point"
    svm_headers = svm_headers.split(',')
    if meta: 
        features_headers = ['Accession_id', 'Locus', 'Genus/Species', 'Nucleotide_Acc', 'Precursor sequence','Leader', 'Core','Start', 'End', 'Total Score', 'Valid Precursor'] + svm_headers 
    else:
        features_headers = ['Accession_id', 'Genus/Species', 'Precursor Sequence','Leader', 'Core','Start', 'End', 'Nearby rSAM Enzyme', 'Multiple rSAM Enzymes', 'Total Score', 'Valid Precursor'] + svm_headers 
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
        self.peptide_type = 'cyclo'
        self.set_split()
        self.find_nearby_rsam()
        self.csv_columns = [self.sequence, self.leader, self.core, self.start, self.end, self.nearby_rsam, self.multi_rsam]
        self.CUTOFF = CUTOFF

        #pfam_2_evalue has a format of: [protein accession, HMM/PFAM hit, e value, start, end, AA length]
        
    def set_split(self):
        self.split_index = int(.5*len(self.sequence))
        self.third_index = int(.67*len(self.sequence))
        
        self.leader = self.sequence[0:self.split_index]
        self.core = self.sequence[self.split_index:]
        self.last_third = self.sequence[self.third_index:]

# Find closest rSAM and return the accessions and rsam_distance
    def find_nearby_rsam(self): 
        rsam_present = False
        general_rsam = 'PF04055'
        ansme_rsam = 'TIGR03942'
        min_dist_dict = {}
        for rsam in self.pfam_2_evalue:
            if rsam[1] == general_rsam or rsam[1] == ansme_rsam:
                min_dist = abs(self.start-0)
                min_dist = min(abs(self.start-rsam[3]), abs(self.end-rsam[3]),
                        abs(self.start-rsam[4]), abs(self.end-rsam[4]),
                        min_dist)
                min_dist_dict[min_dist] = rsam[0]
                rsam_present = True

        if rsam_present == True:
            min_rsam_distance = min(min_dist_dict.keys())
            closest_rsam = min_dist_dict[min_rsam_distance]
            self.nearby_rsam = closest_rsam
        else:
            closest_rsam = 'NA'
            self.nearby_rsam = closest_rsam

# Multiple rSAM enzymes in BGC
        multiple_rsam = False
        multiple_rsam_list = []  
        multi_rsam = 0
        i = 0
        for pfam in self.pfam_2_evalue:
            if pfam[1] == "PF04055" and pfam[0] not in multiple_rsam_list:
                multiple_rsam_list.append(pfam[0])
                i += 1
        if i >= 2:
            multiple_rsam = True
            multi_rsam = str(multiple_rsam_list)
            multi_rsam = multi_rsam.replace(",", " ")
        self.multi_rsam = multi_rsam

# gets the fimo score from the motifs in the fimo file (fimo file is a STREME output)
    def get_fimo_score(self):   
        fimo_output = self.run_fimo_simple()
        fimo_motifs = [line.partition("\t")[0] for line in fimo_output.split("\n")[1:] if "\t" in line]
        fimo_scores = defaultdict(int) # = {int(line.split("\t")[0]): float(line.split("\t")[6]) for line in fimo_output.split("\n") if "\t" in line and line.partition("\t")[0].isdigit()}
        for line in fimo_output.split("\n")[1:]:
            if not ("\t" in line and "motif_id" not in line):
                continue
            fimo_scores[line.split("\t")[0]] = float(line.split("\t")[6])

# if the precursor has a motif
        fimo_motif_score = 0
        fimo_motif_count = 0
        for hits in fimo_motifs:
            fimo_motif_count += 1
        if fimo_motif_count == 0:
            fimo_motif_score = -3
        elif fimo_motif_count >= 1:
            fimo_motif_score += 1 # score
        elif fimo_motif_count >= 2:
            fimo_motif_score +=1 # score
        elif fimo_motif_count >= 4:
            fimo_motif_score +=1 # score
        elif fimo_motif_count >= 12:
            fimo_motif_score +=1 # score   
        return fimo_motifs, fimo_motif_score, fimo_scores
             
#scoring/svm section
    def set_score(self, pfam_dir, cust_hmm):
        scoring_csv_columns = []
        self.score = 0

# Stores precursor PFAM & HMM hits 
        precursor_hmm_info = hmmer_utils.get_hmmer_info(self.sequence, pfam_dir, cust_hmm, self.output_dir)      
        prec_pfams = []
        for pfam_dot, _, pf_evalue, _, in precursor_hmm_info:
            prec_pfams.append(pfam_dot.split('.')[0])

#rSAM nearby (PF04055 is rSAM domain pfam)
        general_rsam = 'PF04055'
        if general_rsam in self.pfam_2_coords.keys():
            self.score += 3
            dist = self.get_min_dist(self.pfam_2_coords[general_rsam])
            scoring_csv_columns.append(1)
        else: 
            dist = 999999
            scoring_csv_columns.append(0)

#rSAM within 2000 nt?
        if dist < 2000:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

#rSAM within 1000 nt?
        if dist < 1000:
            self.score += 2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        
#rSAM more than 5000 nt away
        if dist > 5000:
            self.score -= 3
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

#anSME rsam pfam nearby - most triceptide rSAM enzymes are retrived by the ansme tigrfam due to the presence of a SPASM domain
        ansme_rsam_present = False
        ansme_rsam = "TIGR03942"
        ansme_dict = {}
        for prot in self.pfam_2_evalue:
            if prot[1] == ansme_rsam:
                min_dist = abs(self.start-0)
                #for coord in nearby_rsam_dict[rsam[0]]:
                min_dist = min(abs(self.start-prot[3]), abs(self.end-prot[3]),
                        abs(self.start-prot[4]), abs(self.end-prot[4]),
                        min_dist)
                ansme_dict[min_dist] = prot[2]
                ansme_rsam_present = True

        if ansme_rsam_present == True:
            min_ansme_distance = min(ansme_dict.keys())
            ansme_evalue = ansme_dict[min_ansme_distance]   
        else:
            ansme_evalue = 1
        
#Append the anSME score onto the output CSV file
        scoring_csv_columns.append(ansme_evalue)

        ansme_score = -(math.log10(ansme_evalue))
        scoring_csv_columns.append(ansme_score)
        if ansme_score >= 40 and ansme_score < 100:
            self.score += 2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

        if ansme_score >= 100:
            self.score -= 2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

#Check if a sulfatase is nearby (the truw substrate of anSMEs, filter these from cyclophane/triceptide rSAMs which don't commonly have sulfatases)
        general_sulfatase = 'PF00884'
        if general_sulfatase in self.pfam_2_coords.keys():
            self.score -= 3
            scoring_csv_columns.append(1)
        else: 
            scoring_csv_columns.append(0)

#Check if any of the proteins have a HexxH domain (common on FxS and other triceptides)
        HExxH_domain = 'TIGR04267'
        if HExxH_domain in self.pfam_2_coords.keys():
            self.score += 1
            scoring_csv_columns.append(1)
        else: 
            scoring_csv_columns.append(0)
        
#Penalize for being retrieved by the Wpr RRE
        if 'Wpr_RRE' in prec_pfams:
            self.score -= 3
            scoring_csv_columns.append(1)
        else: 
            scoring_csv_columns.append(0)

#Calculate the precursor's Length
        precursor_length = len(self.sequence)
        
        if precursor_length < 70:
            self.score += 2
            scoring_csv_columns.append(1)
        else: 
            scoring_csv_columns.append(0)

        if precursor_length < 90:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

        if precursor_length > 110:
            self.score -= 2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

        if precursor_length >200:
            self.score -= 4
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

#Count the number of Cys in the potential precursor (One Cys)
        if self.sequence.count("C") == 1:
            self.score -= 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
#Two Cys
        if self.sequence.count("C") == 2:
            self.score -= 3
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
#Three or more Cys
        if self.sequence.count("C") >= 3:
            self.score -= 2*self.sequence.count("C")
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        
#number of potential three-membered cyclophane motifs
        motif = ('[WHFY].[ARNDCEQHILKMFPSTWYV]', 2)
        motif_count = len(re.findall(motif[0], self.core))
        #if re.search(motif[0], self.core) != None:
        if motif_count >= 2:
            self.score += motif[1]
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)  
        scoring_csv_columns.append(motif_count)

#Precent of W/H/Y/F in last third of precursor core 
        per_aromatic = (self.last_third.count("H") + self.last_third.count("F") + self.last_third.count("Y") + self.last_third.count("W"))/float(len(self.last_third))
        if per_aromatic >= 0.13:
            scoring_csv_columns.append(1)
            self.score += 3
        else:
            scoring_csv_columns.append(0)
            
# Precursor pfam/hmm hits (only a handful of hmms were included because having too many affected the SVM training)
        prec_cyclo_hit = False
        prec_cyclo_list = ["TIGR04268", "TIGR04260", "TIGR04495"]
        for hmm in prec_cyclo_list:
            if hmm in prec_pfams:
                prec_cyclo_hit = True
        if prec_cyclo_hit:
            self.score += 3 # score
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0) 

# Precursor does not match any pfams/hmm
        if len(prec_pfams) == 0:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0) 
        
#rgg/shp
        if 'TIGR01716' in self.pfam_2_coords.keys():
            rgg_dist = self.get_min_dist(self.pfam_2_coords['TIGR01716'])
        else:
            rgg_dist = 99999

        if rgg_dist < 500:
            self.score += 2
            scoring_csv_columns.append(1)
        else: 
            scoring_csv_columns.append(0)
        
#Gene synteny (rSAM and precursor in same direction)
        if "PF04055" not in self.pfam_2_coords.keys():
            scoring_csv_columns.append(0)
        elif np.sign(self.start - self.end) == np.sign(self.pfam_2_coords["PF04055"][0][0]-self.pfam_2_coords["PF04055"][0][1]):
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

#common transporters
        trans_hmm = False
        trans_list = ['TIGR03796', 'TIGR01843', 'TIGR03797', 'TIGR01193', 'TIGR03794', 'TIGR03434', 'PF13437']
        for trans in trans_list:
            if trans in self.pfam_2_coords.keys():
                trans_hmm = True
                
        if trans_hmm == True:
            self.score += 2
            scoring_csv_columns.append(1)
        else: 
            scoring_csv_columns.append(0)        
        
# ticeptide rsam hmms
        rsam_hmm = False
        cyclophane_rsam_hmms = ["NF040899", "NF041707", "NF041718", "TIGR04080", "TIGR4083", "TIGR04261", "TIGR04269", "TIGR04403", "TIGR04496", "Darobactin_rSAM", "Wpr_rSAM"]
        for hmm in cyclophane_rsam_hmms:
            if hmm in self.pfam_2_coords.keys():
                tri_rsam_dist = self.get_min_dist(self.pfam_2_coords[hmm])
            else:
                tri_rsam_dist = 99999
            if tri_rsam_dist < 2000:
                    rsam_hmm = True
        if rsam_hmm == True: 
            self.score += 3 # score
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)   

#Potential precursor is an NCBI inferred gene
        if 'X' in self.sequence:
            self.score -= 10
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)      

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

    
# SVM SCORING SECTION
        scoring_csv_columns.append(dist)
        scoring_csv_columns.append(precursor_length)

        charge_dict = {"E": -1, "D": -1, "K": 1, "H": 1, "R": 1}
        precursor_charge = sum([charge_dict[aa] for aa in self.sequence if aa in charge_dict]) 
        core_charge = sum([charge_dict[aa] for aa in self.core if aa in charge_dict])
        last_third_charge = sum([charge_dict[aa] for aa in self.core if aa in charge_dict])
        scoring_csv_columns.append(core_charge)
        scoring_csv_columns.append(precursor_charge)
        scoring_csv_columns.append(last_third_charge)
        scoring_csv_columns.append(abs(core_charge))
        scoring_csv_columns.append(abs(precursor_charge))
        scoring_csv_columns.append(abs(last_third_charge))

        #Counts of AAs in core
        scoring_csv_columns += [self.core.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
        #Percentages of AAs in core
        core_analysis = ProteinAnalysis(self.core, monoisotopic=True) 
        for aa in "ARDNCQEGHILKMFPSTWYV":
            scoring_csv_columns.append(core_analysis.amino_acids_percent[aa])
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
            scoring_csv_columns.append(seq_analysis.amino_acids_percent[aa])
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

        #Counts of AAs in core
        scoring_csv_columns += [self.last_third.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
        #Percentages of AAs in core
        last_third_analysis = ProteinAnalysis(self.last_third, monoisotopic=True) 
        for aa in "ARDNCQEGHILKMFPSTWYV":
            scoring_csv_columns.append(last_third_analysis.amino_acids_percent[aa])
        #Aromatics in core
        scoring_csv_columns.append(sum([self.last_third.count(aa) for aa in "FWY"]))
        #Neg charged in core
        scoring_csv_columns.append(sum([self.last_third.count(aa) for aa in "DE"]))
        #Pos charged in core
        scoring_csv_columns.append(sum([self.last_third.count(aa) for aa in "RK"]))
        #Charged in core
        scoring_csv_columns.append(sum([self.last_third.count(aa) for aa in "RKDE"]))
        #Aliphatic in core
        scoring_csv_columns.append(sum([self.last_third.count(aa) for aa in "GAVLMI"]))
        #Hydroxyl in core
        scoring_csv_columns.append(sum([self.last_third.count(aa) for aa in "ST"]))
        #isoelectric point in core
        scoring_csv_columns.append(last_third_analysis.isoelectric_point())

        self.csv_columns += [self.score] +  scoring_csv_columns
