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

import My_Record
import entrez_utils
import record_processing
import math

import hmmer_utils
from collections import defaultdict
import pathlib
FILE_DIR = pathlib.Path(__file__).parent.absolute()
peptide_type = "boro"
CUTOFF = 16
 
index = 0


#sets up the csv file with correct headers

def write_csv_headers(output_dir):
    dir_prefix = output_dir + '/boro/'
    if not os.path.exists(dir_prefix):
        os.makedirs(dir_prefix)
    svm_headers = "Precursor Index,Classification,Precursor has a methyltransferase domain,Precursor has a BBD domain,Precursor's Length,Best Boro MT length,Precursor is best Boro MT,Best BBD length,Precursor is best BBD,MT & BBD fused,MT & GGDEF fused,Nearby PF00590 domain,Nearby BorosinMT hmm hit,Nearby NMT_1 (DCLFAD) hmm hit,Nearby NMT_2 (YGHP) hmm hit,Nearby PF03819 (MazG) or PF12643 (MazG-like) domain,Best BorosinMT score [-log(e-value)],BorosinMT score >= 40,BorosinMT score >= 30,BorosinMT score <= 25,BorosinMT score <=15,Best BorosinMT hit's PF00590 (TP methyltransferase) score [-log(e-value)],BorosinMT hit's PF00590 score >= 35,BorosinMT hit's PF00590 score >= 25,Best BorosinMT hit's PF03819 (MazG) score [-log(e-value)],BorosinMT hit's PF03819 (MazG) score >= 35,BorosinMT hit's PF03819 (MazG) score >= 25,BorosinMT hit's PF03819 (MazG) score >= 15,BorosinMT hit's PF03819 (MazG) score >= 5,Nearby PF07746 (LigA),Nearby BBD hmm hit,PF07746 (LigA) score [-log(e-value)],BBD_A score [-log(e-value)],BBD_B score [-log(e-value)],BBD_C score [-log(e-value)],BBD_old score [-log(e-value)],Precursor within 100 nucleotides of a MT,Precursor within 200 nucleotides of a MT,Precursor within 1000 nucleotides of a MT,Precursor is farther than 3000 nucleotides from a MT,Precursor & MT on same strand,Both BBD and PF00590 domains nearby,Nearby GGDEF,Nearby acetyltransferase,Nearby peptidase,Nearby legionellales hmm hit,Nearby type VI/VII/VIII hmm hit,precursor has a PF00590 (TP methyltransferase) domain,precursors hits the BorosinMT hmm,precursor hits the NMT_1 hmm,precursor hits the NMT_2 hmm,precursor has a LigA-like domain (PF07746),precursor hits BBD_A hmm,precursor hits BBD_B hmm,precursor hits BBD_C hmm,precursor hits BBD_old hmm,precursor hits a legionellales hmm,precursor hits type VI/VII/VIII hmm,precursor hits MT & BBD hmm,precursor hits MT hmm and nearby type 6-8 hmm hit,precursor hits a legionellales & BBD hmm,precursor hits type VI/VII/VIII & BBD hmm,precursor hits a legionellales & MT hmm,precursor hits type VI/VII/VIII & MT hmm,precursor is > 900 AA & hits a BBD hmm,precursor is > 700 AA & hits a BBD hmm,precursor is < 400 AA & hits a BBD hmm,precursor is < 250 AA & hits a BBD hmm,precursor is < 100 AA and hits a BBD hmm,precursor is < 400 AA and hits the type VI/VII/VIII hmm,precursor is < 100 AA and hits a legionellales hmm,FIMO motif 1 present,motif 2,motif 3,motif 4,motif 5,motif 6,motif 7,motif 8,motif 9,motif 10,motif 11,motif 12,motif 13,motif 14,motif 15,motif 16,motif 17,motif 18,motif 19,motif 20,motif 21,motif 22,motif 23,motif 24,motif 25,motif 26,motif 27,Number of motifs present,Occurances of motif 1,motif 2,motif 3,motif 4,motif 5,motif 6,motif 7,motif 8,motif 9,motif 10,motif 11,motif 12,motif 13,motif 14,motif 15,motif 16,motif 17,motif 18,motif 19,motif 20,motif 21,motif 22,motif 23,motif 24,motif 25,motif 26,motif 27,No Motifs,MT distance,core charge,precursor charge,abs(core charge),abs(precursor charge),CORE COUNT A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,PERCENT A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,aromatics,neg,pos,charged,aliphatic,hydroxyl,isoelectric point,PRECURSOR count A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,PERCENT A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,aromatics,neg,pos,charged,aliphatic,hydroxyl,isoelectric point"
    svm_headers = svm_headers.split(',')
    features_headers = ['Accession_id', 'Genus/Species', 'Sequence', 'Region1', 'Region2', 'Region3', 'Start', 'End', 'Best Borosin MT accession', 'Best BBD accession', 'Multiple Borosin MT', 'Multiple BBDs', 'Total Score',"Valid Precursor",] + svm_headers
#TODO close all these write headers at the end of the document
    features_csv_file = open(dir_prefix + "temp_features.csv", 'w')
    svm_csv_file = open("{}fitting_set.csv".format(dir_prefix), 'w')
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
                 pfam_2_coords,
                 pfam_2_evalue):
        
        super(Ripp, self).__init__(start, 
                                     end, 
                                     sequence,
                                     upstream_sequence,
                                     pfam_2_coords,
                                     pfam_2_evalue)
                                     
        self.peptide_type = 'boro'
        self.set_split()
        self.best_boro_acc()
        self.csv_columns = [self.sequence, self.region1, self.region2, self.region3, self.start, self.end, self.best_boro_mt, self.best_boro_bbd, self.p_mul_mt, self.p_mul_bbd]
        self.CUTOFF = CUTOFF
        #pfam_2_evalue has a format of: [protein accession, HMM/PFAM hit, e value, start, end, AA length]

        
        
    
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
    
    def best_boro_acc(self):
# Best scoring Borosin methyltransferase
        best_boro_mt1 = 0
        best_boroMT_evalue1 = 1
        for mts in self.pfam_2_evalue:
            if mts[1] == "BorosinMT" and mts[2] < best_boroMT_evalue1:
                best_boro_mt1 = mts[0]
                best_boroMT_evalue1 = mts[2]
        self.best_boro_mt = best_boro_mt1
        
# Best scoring BBD
        bbd_list = ["PF07746", "BBD_A", "BBD_B", "BBD_C", "BBD_old"]
        best_boro_bbd1 = 0
        best_boro_bbd_evalue1 = 1
        for bbds in self.pfam_2_evalue:
            if bbds[1] in bbd_list and bbds[2] < best_boro_bbd_evalue1:
                best_boro_bbd1 = bbds[0]
                best_boro_bbd_evalue1 = bbds[2]
        self.best_boro_bbd = best_boro_bbd1
                
# Multiple borosin methyltransferases in BGC
        multi_boro_mt = False
        multi_boro_mt_list = []
        p_mul_mt = 0
        i = 0
        for mts in self.pfam_2_evalue:
            if mts[1] == "BorosinMT" and mts[0] not in multi_boro_mt_list:
                multi_boro_mt_list.append(mts[0])
                i += 1
        if i >= 2:
            multi_boro_mt = True
            p_mul_mt = str(multi_boro_mt_list)
            p_mul_mt = p_mul_mt.replace(",", " ")
        self.p_mul_mt = p_mul_mt 
        
# Multiple BBDs in BGC
        multi_boro_bbd = False
        bbd_list = ["BBD_A", "BBD_B", "BBD_C", "BBD_old", "PF07746"]
        multi_boro_bbd_list = []
        p_mul_bbd = 0
        i = 0
        for bbds in self.pfam_2_evalue:
            if bbds[1] in bbd_list and bbds[0] not in multi_boro_bbd_list:
                multi_boro_bbd_list.append(bbds[0])
                i += 1
        if i >= 2:
            multi_boro_bbd = True
            p_mul_bbd = str(multi_boro_bbd_list)
            p_mul_bbd = p_mul_bbd.replace(",", " ")
        self.p_mul_bbd = p_mul_bbd
        
        
        
# gets the fimo score from the motifs in the fimo file (fimo file is a STREME output)
    def get_fimo_score(self):   
        fimo_output = self.run_fimo_simple()
        fimo_motifs = [line.partition("\t")[0] for line in fimo_output.split("\n")[1:] if "\t" in line]
        fimo_scores = defaultdict(int) # = {int(line.split("\t")[0]): float(line.split("\t")[6]) for line in fimo_output.split("\n") if "\t" in line and line.partition("\t")[0].isdigit()}
        for line in fimo_output.split("\n")[1:]:
            if not ("\t" in line and "motif_id" not in line):
                continue
            fimo_scores[line.split("\t")[0]] = float(line.split("\t")[6])
 # if the precursor has a motif # TODO
        fimo_motif_score = 0
        fimo_motif_count = 0
        for hits in fimo_motifs:
            fimo_motif_count += 1
        if fimo_motif_count == 0:
            fimo_motif_score = -3
        elif fimo_motif_count >= 3:
            fimo_motif_score += 1 # score
        elif fimo_motif_count >= 5:
            fimo_motif_score +=1 # score
        elif fimo_motif_count >= 7:
            fimo_motif_score +=1 # score
        return fimo_motifs, fimo_motif_score, fimo_scores

# SCORING & SVM SECTION
    def set_score(self, pfam_dir, cust_hmm): 
        scoring_csv_columns = []
        self.score = 0

# Stores precursor PFAM & HMM hits 

        precursor_hmm_info = hmmer_utils.get_hmmer_info(self.sequence, pfam_dir, cust_hmm)      
        prec_pfams = []
        for pfam_dot, _, pf_evalue, _, in precursor_hmm_info:
            prec_pfams.append(pfam_dot.split('.')[0])
            
# Precursor is a methyltransferase
        prec_mt_hit = False
        mt_list = ["BorosinMT", "NMT_1", "NMT_2", "PF00590"]
        for tar in mt_list:
            if tar in prec_pfams:
                prec_mt_hit = True
        if prec_mt_hit:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0) 
        
# Precursor has a BBD domain
        prec_bbd_hit = False
        bbd_list = ["BBD_A", "BBD_B", "BBD_C", "BBD_old", "PF07746"]
        for tar in bbd_list:
            if tar in prec_pfams:
                prec_bbd_hit = True
        if prec_bbd_hit:
            self.score += 3 # score
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
# Precursor's Length
        scoring_csv_columns.append(len(self.sequence))
        
# Best scoring Borosin methyltransferase
        boromt_hmm = False
        best_boro_mt = ""
        best_boroMT_evalue = 1
        best_boro_mt_len = 0
        best_boro_mt_start = 0
        best_boro_mt_end = 0
        
        for mts in self.pfam_2_evalue:
            if mts[1] == "BorosinMT" and mts[2] < best_boroMT_evalue:
                best_boro_mt = mts[0]
                best_boroMT_evalue = mts[2]
                best_boro_mt_start = mts[3]
                best_boro_mt_end = mts[4]
                best_boro_mt_len = mts[5]
                boromt_hmm = True
            
        if not boromt_hmm:
            self.score -+ 5 # score

# Best scoring Borosin methyltransferase AA length
        if boromt_hmm:
            scoring_csv_columns.append(best_boro_mt_len)
        else:
            scoring_csv_columns.append(0)

# Precursor is best scoring Borosin methyltransferase
        if boromt_hmm and (best_boro_mt_start == self.start) and (best_boro_mt_end == self.end):
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# Best scoring BBD 
        borobbd_hmm = False
        bbd_list = ["PF07746", "BBD_A", "BBD_B", "BBD_C", "BBD_old"]
        best_boro_bbd = ""
        best_boro_bbd_evalue = 1
        best_boro_bbd_len = 0
        best_boro_bbd_start = 0
        best_boro_bbd_end = 0

        for bbds in self.pfam_2_evalue:
            if bbds[1] in bbd_list and bbds[2] < best_boro_bbd_evalue:
                best_boro_bbd = bbds[0]
                best_boro_bbd_evalue = bbds[2]
                best_boro_bbd_start = bbds[3]
                best_boro_bbd_end = bbds[4]
                best_boro_bbd_len = bbds[5]
                borobbd_hmm = True
        if not borobbd_hmm:
            self.score -= 2

# Best scoring BBD AA length
        if borobbd_hmm:
            scoring_csv_columns.append(best_boro_bbd_len)
        else:
            scoring_csv_columns.append(0)

# Precursor is best scoring BBD
        if borobbd_hmm and (best_boro_bbd_start == self.start) and (best_boro_bbd_end == self.end):
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# MT & BBD fused in BGC
        fused_mt_bbd = False
        mt_list = ["BorosinMT", "NMT_1", "NMT_2", "PF00590"]
        bbd_list = ["BBD_A", "BBD_B", "BBD_C", "BBD_old", "PF07746"]
        mt_temp = []
        bbd_temp = []
        for prots in self.pfam_2_evalue:
            if prots[1] in mt_list and prots[0] not in mt_temp:
                mt_temp.append(prots[0])
            if prots[1] in bbd_list and prots[0] not in bbd_temp:
                bbd_temp.append(prots[0])
        for mts in mt_temp:
            if mts in bbd_temp:
                fused_mt_bbd = True
        if fused_mt_bbd:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
# MT & GGDEF fused in BGC
        fused_mt_ggdef = False
        mt_list = ["BorosinMT", "NMT_1", "NMT_2"]
        ggdef_list = ["PF00990", "PF17853"]
        mt_temp = []
        ggdef_temp = []
        for prots in self.pfam_2_evalue:
            if prots[1] in mt_list and prots[0] not in mt_temp:
                mt_temp.append(prots[0])
            if prots[1] in ggdef_list and prots[0] not in ggdef_temp:
                ggdef_temp.append(prots[0])
        for mts in mt_temp:
            if mts in ggdef_temp:
                fused_mt_ggdef = True
                
        if fused_mt_ggdef:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# if anything nearby has a PF00590 (TP methyltransferase)
        tpPfam = False
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in ["PF00590"]):
                tpPfam = True
        if tpPfam:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)   

# if anything nearby hits the BorosinMT hmm
        mtHmm = False
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in ["BorosinMT"]):
                mtHmm = True
        if mtHmm:
            self.score += 1 # score
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
# if anything nearby hits the NMT_1 (DCLFAD) hmm
        dclfad = False
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in ["NMT_1"]):
                dclfad = True
        if dclfad:
            self.score += 1 # score
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
# if anything nearby hits the NMT_2 (YGHP) hmm
        yghp = False
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in ["NMT_2"]):
                yghp = True
        if yghp:
            self.score += 1 # score
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
# if anything nearby has a PF03819 (MazG) or PF12643 (MazG-like) domain
        mazgPfam = False
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in ["PF03819", "PF12643"]):
                mazgPfam = True
        if mazgPfam:
            self.score -= 4 # score
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)   
            
# Methyltransferase BorosinMT hmm -log(evalue)
        boro_mt_score = -(math.log10(best_boroMT_evalue))
        scoring_csv_columns.append(boro_mt_score)
        
# Methyltransferase hits BorosinMT hmm with a -log(evalue) >= 40
        boro_mt_40 = False
        if boro_mt_score >= 40:
            scoring_csv_columns.append(1)
            self.score += 1 # score
            boro_mt_40 = True
        else:
            scoring_csv_columns.append(0)
        
# Methyltransferase hits BorosinMT hmm with a -log(evalue) >= 30
        boro_mt_30 = False
        if boro_mt_score >= 30:
            scoring_csv_columns.append(1)
            self.score += 1 # score
            boro_mt_30 = True
        else:
            scoring_csv_columns.append(0)
        
# Methyltransferase hits BorosinMT hmm with a -log(evalue) =< 25
        boro_mt_25 = False
        if boro_mt_score <= 25:
            scoring_csv_columns.append(1)
            self.score -= 2 # score
            boro_mt_25 = True
        else:
            scoring_csv_columns.append(0)
            
# Methyltransferase hits BorosinMT hmm with a -log(evalue) =< 15
        boro_mt_15 = False
        if boro_mt_score <= 15:
            scoring_csv_columns.append(1)
            self.score -= 3 # score
            boro_mt_15 = True
        else:
            scoring_csv_columns.append(0)
        
#Best BorosinMT other mt pfam hits evalues
        tp_evalue = False
        mazg_evalue = False
        best_boro_tp_evalue = 1
        best_boro_mazg_evalue = 1
        neglog_tp = 1
        neglog_mazg = 1
        
        if boromt_hmm == True:
            for mts in self.pfam_2_evalue:
                if mts[0] == best_boro_mt and mts[1] == "PF00590":
                    best_boro_tp_evalue = mts[2]
                    tp_evalue = True
                
                if mts[0] == best_boro_mt and mts[1] == "PF03819":
                    best_boro_mazg_evalue = mts[2]
                    mazg_evalue = True
                    
#Best BorosinMT hit's PF00590 -log(evalue)
        record_pf00590 = False
        # if the BorosinMT hmm is hit and that protein has a PF00590 domain
        if boromt_hmm == True and tp_evalue == True:
            neglog_tp = -(math.log10(best_boro_tp_evalue))
            scoring_csv_columns.append(neglog_tp)
            record_pf00590 = True
        # if the BorosinMT hmm is hit and that protein doesn't have a PF00590 domain
        elif boromt_hmm == True and tp_evalue == False:
            neglog_tp = 0
            scoring_csv_columns.append(0)
            record_pf00590 = True
        # if the BorosinMT hmm isn't hit
        else:
            neglog_tp = -1
            scoring_csv_columns.append(-1)

#Best BorosinMT hit's PF00590 -log(evalue) >= 35
        if record_pf00590:
            if neglog_tp >= 35:
                self.score -= 3 # score
                scoring_csv_columns.append(1)
            else:
                scoring_csv_columns.append(0)
        else:
            scoring_csv_columns.append(-1)
        

#Best BorosinMT hit's PF00590 -log(evalue) >= 25
        if record_pf00590:
            if neglog_tp >= 25:
                self.score -= 2 # score
                scoring_csv_columns.append(1)
            else:
                scoring_csv_columns.append(0)
        else:
            scoring_csv_columns.append(-1)
            
#Best BorosinMT hit's PF03819 -log(evalue)
        record_pf03819 = False
        # if the BorosinMT hmm is hit and that protein has a PF03819 domain
        if boromt_hmm == True and mazg_evalue == True:
            neglog_mazg = -(math.log10(best_boro_mazg_evalue))
            scoring_csv_columns.append(neglog_mazg)
            record_pf03819 = True
        # if the BorosinMT hmm is hit and that protein doesn't have a PF00590 domain
        elif boromt_hmm == True and mazg_evalue == False:
            neglog_mazg = 0
            scoring_csv_columns.append(0)
            record_pf03819 = True
        # if the BorosinMT hmm isn't hit
        else:
            neglog_tp = -1
            scoring_csv_columns.append(-1)
        

#Best BorosinMT hit's PF03819 -log(evalue) >= 35
        if record_pf03819:
            if neglog_mazg >= 35:
                self.score -= 5 # score
                scoring_csv_columns.append(1)
            else:
                scoring_csv_columns.append(0)
        else:
            scoring_csv_columns.append(-1)

#Best BorosinMT hit's PF03819 -log(evalue) >= 25
        if record_pf03819:
            if neglog_mazg >= 25:
                self.score -= 3 # score
                scoring_csv_columns.append(1)
            else:
                scoring_csv_columns.append(0)
        else:
            scoring_csv_columns.append(-1)

#Best BorosinMT hit's PF03819 -log(evalue) >= 15
        if record_pf03819:
            if neglog_mazg >= 15:
                self.score -= 1 # score
                scoring_csv_columns.append(1)
            else:
                scoring_csv_columns.append(0)
        else:
            scoring_csv_columns.append(-1)

#Best BorosinMT hit's PF03819 -log(evalue) >= 5
        if record_pf03819:
            if neglog_mazg >= 5:
                scoring_csv_columns.append(1)
            else:
                scoring_csv_columns.append(0)
        else:
            scoring_csv_columns.append(-1)
        

# if anything nearby has a PF07746 domain (LigA)
        bbd = False
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in ["PF07746"]):
                bbd = True
        if bbd:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# if the nearby BBD actually matches one or both of the hmms (not additive)        
        BBDhmm = False
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in ["BBD_A", "BBD_B", "BBD_C", "BBD_old"]):
                BBDhmm = True
        if BBDhmm:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
# BBD or ligA score
        if BBDhmm or bbd:
            self.score += 2 # score
            
# LigA & BBD best hit's e-values
        ligA_evalue = False
        bbdA_evalue = False
        bbdB_evalue = False
        bbdC_evalue = False
        bbdOld_evalue = False
        best_LigA_evalue = 1
        best_bbdA_evalue = 1
        best_bbdB_evalue = 1
        best_bbdC_evalue = 1
        best_bbdOld_evalue = 1
        neglog_LigA_evalue = 1
        neglog_bbdA_evalue = 1
        neglog_bbdB_evalue = 1
        neglog_bbdC_evalue = 1
        neglog_bbdOld_evalue = 1
        
        if borobbd_hmm == True:
            for bbds in self.pfam_2_evalue:
                if bbds[0] == best_boro_bbd and bbds[1] == "PF07746":
                    best_LigA_evalue = bbds[2]
                    ligA_evalue = True
                if bbds[0] == best_boro_bbd and bbds[1] == "BBD_A":
                    best_bbdA_evalue = bbds[2]
                    neglog_bbdA_evalue = -(math.log10(best_bbdA_evalue))
                    bbdA_evalue = True
                if bbds[0] == best_boro_bbd and bbds[1] == "BBD_B":
                    best_bbdB_evalue = bbds[2]
                    neglog_bbdB_evalue = -(math.log10(best_bbdB_evalue))
                    bbdB_evalue = True
                if bbds[0] == best_boro_bbd and bbds[1] == "BBD_C":
                    best_bbdC_evalue = bbds[2]
                    neglog_bbdC_evalue = -(math.log10(best_bbdC_evalue))
                    bbdC_evalue = True
                if bbds[0] == best_boro_bbd and bbds[1] == "BBD_old":
                    best_bbdOld_evalue = bbds[2]
                    neglog_bbdOld_evalue = -(math.log10(best_bbdOld_evalue))
                    bbdOld_evalue = True

# PF07746 (LigA) -log(evalue)
        # if a LigA domain or BBD hmm is present and that protein has a PF07746 (LigA) domain
        record_pf07746 = False
        if borobbd_hmm == True and ligA_evalue == True:
            neglog_LigA_evalue = -(math.log10(best_LigA_evalue))
            scoring_csv_columns.append(neglog_LigA_evalue)
            record_pf07746 = True
            
        #if a LigA domain or BBD hmm is present and that protein doesn't have a PF07746 (LigA) domain
        elif borobbd_hmm == True and ligA_evalue == False:
            neglog_LigA_evalue = 0
            scoring_csv_columns.append(0)
            record_pf07746 = True
            
        # if the LigA domain isn't hit
        else:
            neglog_LigA_evalue = -1
            scoring_csv_columns.append(-1)
            
# BBD_A hmm -log(evalue)
        # if a LigA domain or BBD hmm is present and that protein hits the BBD_A hmm
        record_bbdA = False
        if borobbd_hmm == True and bbdA_evalue == True:
            neglog_bbdA_evalue = -(math.log10(best_bbdA_evalue))
            scoring_csv_columns.append(neglog_bbdA_evalue)
            record_bbdA = True
            
        #if a LigA domain or BBD hmm is present and that protein doesn't hit the BBD_A hmm
        elif borobbd_hmm == True and bbdA_evalue == False:
            neglog_bbdA_evalue = 0
            scoring_csv_columns.append(0)
            record_bbdA = True
            
        # if the BBD_A hmm isn't hit
        else:
            neglog_bbdA_evalue = -1
            scoring_csv_columns.append(-1)
        
# BBD_B hmm -log(evalue)
        # if a LigA domain or BBD hmm is present and that protein hits the BBD_B hmm
        record_bbdB = False
        if borobbd_hmm == True and bbdB_evalue == True:
            neglog_bbdB_evalue = -(math.log10(best_bbdB_evalue))
            scoring_csv_columns.append(neglog_bbdB_evalue)
            record_bbdB = True
            
        #if a LigA domain or BBD hmm is present and that protein doesn't hit the BBD_B hmm
        elif borobbd_hmm == True and bbdB_evalue == False:
            neglog_bbdB_evalue = 0
            scoring_csv_columns.append(0)
            record_bbdB = True
            
        # if the BBD_B hmm isn't hit
        else:
            neglog_bbdB_evalue = -1
            scoring_csv_columns.append(-1)
        
# BBD_C hmm -log(evalue)
        # if a LigA domain or BBD hmm is present and that protein hits the BBD_C hmm
        record_bbdC = False
        if borobbd_hmm == True and bbdC_evalue == True:
            neglog_bbdC_evalue = -(math.log10(best_bbdC_evalue))
            scoring_csv_columns.append(neglog_bbdC_evalue)
            record_bbdC = True
            
        #if a LigA domain or BBD hmm is present and that protein doesn't hit the BBD_C hmm
        elif borobbd_hmm == True and bbdC_evalue == False:
            neglog_bbdC_evalue = 0
            scoring_csv_columns.append(0)
            record_bbdC = True
            
        # if the BBD_C hmm isn't hit
        else:
            neglog_bbdC_evalue = -1
            scoring_csv_columns.append(-1)
            
# BBD_old hmm -log(evalue)
        # if a LigA domain or BBD hmm is present and that protein hits the BBD_old hmm
        record_bbdOld = False
        if borobbd_hmm == True and bbdOld_evalue == True:
            neglog_bbdOld_evalue = -(math.log10(best_bbdOld_evalue))
            scoring_csv_columns.append(neglog_bbdOld_evalue)
            record_bbdOld = True
            
        #if a LigA domain or BBD hmm is present and that protein doesn't hit the BBD_old hmm
        elif borobbd_hmm == True and bbdOld_evalue == False:
            neglog_bbdOld_evalue = 0
            scoring_csv_columns.append(0)
            record_bbdOld = True
            
        # if the BBD_old hmm isn't hit
        else:
            neglog_bbdOld_evalue = -1
            scoring_csv_columns.append(-1)

# # adjusts score based on precursor distance from MT
        MT = ["PF00590", "BorosinMT", "NMT_1", "NMT_2" ] 
        MT_coords = []
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in MT):
                MT_coords += self.pfam_2_coords[pfam]
        min_distance_MT = self.get_min_dist(MT_coords)
        if min_distance_MT is None:
            min_distance_MT = 66666 
        within_100 = False
        within_200 = False
        within_1000 = False
        within_3000 = False
        if min_distance_MT < 3000:
            within_3000 = True
            if min_distance_MT < 1000:
                within_1000 = True
                if min_distance_MT < 200:
                    within_200 = True
                    if min_distance_MT < 100:
                        within_100 = True
        if within_100:
            self.score += 1 # score
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if within_200:
            self.score += 1 # score
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if within_1000:
            self.score += 1 # score
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if not within_3000:
            self.score -= 2 # score
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
# if the precursor is encoded on the same strand as the MT
        if self.start < self.end:
            direction = 1
        else:
            direction = -1        
        same_dir = False
        for tag in MT:
            if tag in self.pfam_2_coords.keys():
                for coord in self.pfam_2_coords[tag]:
                    if (coord[0] < coord[1] and direction == 1) or coord[0] > coord[1] and direction == -1:
                        same_dir = True
        if same_dir:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# if BBD and PF00590 both occur in BGC
        BBD_PFAM = False
        if (BBDhmm or bbd) and tpPfam == True:
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
            
# if a legionellales hmm hit co-occurs
        leg_cooccur = False
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in ["BoroLegA", "BoroLegB"]):
                leg_cooccur = True
        if leg_cooccur:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# if a Boro_6_8 hmm hit co-occurs
        Boro_6_8_cooccur = False
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in ["Boro_6_8"]):
                Boro_6_8_cooccur = True
        if Boro_6_8_cooccur:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# if the precursor has a PF00590 domain
        if "PF00590" in prec_pfams:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# if the precursor hits the BorosinMT hmm
        if "BorosinMT" in prec_pfams:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# if the precursor hits the NMT_1 (DCLFAD) hmm
        if "NMT_1" in prec_pfams:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# if the precursor hits the NMT_2 (YGHP)hmm
        if "NMT_2" in prec_pfams:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# if the precursor has a LigA-like domain
        if "PF07746" in prec_pfams:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# if the precursor hits the first BBD hmm
        if "BBD_A" in prec_pfams:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# if the precursor hits the second BBD hmm
        if "BBD_B" in prec_pfams:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
# if the precursor hits the third BBD hmm
        if "BBD_C" in prec_pfams:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# if the precursor hits the old BBD hmm
        if "BBD_old" in prec_pfams:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
# if the precursor hits either of the legionellales hmms
        precursor_leg_hit = False
        targets = ["BoroLegA", "BoroLegB"]
        for tar in targets:
            if tar in prec_pfams:
                precursor_leg_hit = True
        if precursor_leg_hit == True: 
            self.score += 4 # score
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0) 

# if the precursor hits the type VI/VIII hmm
        precursor_6_8_hit = False
        if "Boro_6_8" in prec_pfams:
            self.score += 4 # score
            scoring_csv_columns.append(1)
            precursor_6_8_hit = True
        else:
            scoring_csv_columns.append(0)
            
# if precursor hits MT & BBD hmm
        if prec_mt_hit and prec_bbd_hit:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
# if precursor hits MT hmm and nearby type 6-8 hmm hit
        if Boro_6_8_cooccur and prec_mt_hit:
            self.score -= 10 # score
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
                                 
# if precursor hits legionellales & BBD hmm
        if precursor_leg_hit and prec_bbd_hit:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)                                 
                                 
# if precursor hits type VI/VIII & BBD hmm
        if precursor_6_8_hit and prec_bbd_hit:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)                                 
#                                   
# if precursor hits legionellales & MT hmm
        if precursor_leg_hit and prec_mt_hit:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)                                 
                                 
# if precursor hits type VI/VIII & MT hmm
        if precursor_6_8_hit and prec_mt_hit:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)         

# if precursor hits is >900AA & hits BBD hmm
        if int(len(self.sequence)) > 900 and prec_bbd_hit:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0) 

# if precursor hits is >700AA & hits BBD hmm
        if int(len(self.sequence)) > 700 and prec_bbd_hit:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0) 
            
# if precursor hits is <400AA & hits BBD hmm
        if int(len(self.sequence)) < 400 and prec_bbd_hit:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

# if precursor hits is <250AA & hits BBD hmm
        if int(len(self.sequence)) < 250 and prec_bbd_hit:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
# if precursor hits is <100AA & hits BBD hmm
        if int(len(self.sequence)) < 100 and prec_bbd_hit:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
# if precursor hits is <400AA & hits the Boro_6_8 hmm
        if int(len(self.sequence)) < 400 and precursor_6_8_hit:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
# if precursor hits is <100AA & hits a legionellales hmm
        if int(len(self.sequence)) < 100 and precursor_leg_hit:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
                    
# AA % of "core" (region 2) for DEIVT
        per_asp = self.core.count("D")/float(len(self.core))
        if per_asp > 0.10: 
            self.score += 1 # score
        per_glu =  self.core.count("E")/float(len(self.core))   
        if per_glu > 0.10:
            self.score += 1 # score           

        per_ile =  self.core.count("I")/float(len(self.core))   
        if per_ile > 0.10:
            self.score += 1 # score        

        per_val =  self.core.count("V")/float(len(self.core))   
        if per_val > 0.10:
            self.score += 1 # score

        per_thr = self.core.count("T")/float(len(self.core))
        if per_thr > 0.10:
            self.score += 1 # score


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
        scoring_csv_columns.append(min_distance_MT)

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
 
     


