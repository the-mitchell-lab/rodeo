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
import re
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from ripp_modules.VirtualRipp import VirtualRipp
import pathlib
FILE_DIR = pathlib.Path(__file__).parent.absolute()

peptide_type = "thio"
CUTOFF = 20
index = 0

def write_csv_headers(output_dir, meta=False):
    dir_prefix = output_dir + '/thio/'
    if not os.path.exists(dir_prefix):
        os.makedirs(dir_prefix)
    svm_headers = 'PK,Classification,Contains TOMM YcaO PF02624,Contains LanB Nterm PF04738,Contains LanB Cterm PF14028,Contains TOMM dehy PF00881,Contains rSAM MTase PF04055,Contains P450 PF00067,Contains ABC trans1 PF00005,Contains ABC trans2 PF01061,Contains ABC trans3 PF12698,Contains abhydrolase1 PF12697,Contains abhydrolase2 PF00561,CSS motif,CTT motif,SS motif,SSS motif,SSSS motif,CC motif,CCC motif,CCCC motif,TT motif,TTT motif,TTTT motif,No Cys core residues,No Ser core residues,No Thr core residues,Core mass < 2100,Sum of repeating Cys/Ser/Thr > 4,Avg heterocycle block length > 3,Leader net charge < 5,Leader net charge > 0,Leader contains a Cys?,Peptide terminates Cys/Ser/Thr,Core contains >= 2 positive residues,Heterocycle ratio > 0.4,Number of core repeating blocks,Number of core repeating Cys,Number of core repeating Ser,Number of core repeating Thr,Number of core heterocycle blocks,avg core heterocycle block length,Precursor peptide mass (unmodified),Leader peptide mass (unmodified),Core peptide mass (unmodified),Length of Precursor,Length of Leader,Length of Core,Leader / core ratio,Heterocycle residues/lenth of core,A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,Aromatics,Neg charged,Pos charged,Charged,Aliphatic,Hydroxyl,A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,Aromatics,Neg charged,Pos charged,Charged,Aliphatic,Hydroxyl,A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,Aromatics,Neg charged,Pos charged,Charged,Aliphatic,Hydroxyl'
    svm_headers = svm_headers.split(',')
    if meta:
        features_headers = ["Accession_id", "Locus", "Genus/Species/Code", 'Nucleotide_Acc', "Leader", "Core", "Start", "End", "Total Score", "Valid Precursor" ] + svm_headers
    else:
        features_headers = ["Accession_id", "Genus/Species/Code", "Leader", "Core", "Start", "End", "Total Score", "Valid Precursor" ] + svm_headers
    features_csv_file = open(dir_prefix + "temp_features.csv", 'w')
    svm_csv_file = open("{}fitting_set.csv".format(dir_prefix), 'w')
    features_writer = csv.writer(features_csv_file)
    svm_writer = csv.writer(svm_csv_file)
    features_writer.writerow(features_headers)
    svm_writer.writerow(svm_headers)

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
        self.peptide_type = 'thio'
        self.set_split()
#        self.set_monoisotopic_mass()
        self.csv_columns = [self.leader, self.core, self.start, self.end]
        self.CUTOFF = CUTOFF
        
    def set_split(self):
        """Try to identify cleavage site using regular expressions"""
        #Regular expressions; try 1 first, then 2, etc.
#        rex1 = re.compile('([I|V]AS)')
#        rex2 = re.compile('([G|A|S]AS)')
#        rex3 = re.compile('([S|A|T][Q|K|T][V|A|I|M]M[A|S]A)')
#        rex4 = re.compile('([L|M|I|V][P|S|T|V|M][E|D][T|M|V|G|N|I|L|A|S][G|A|T|S]A)')
#        rex5 = re.compile('(DL[T|S][V|E]T[M|L])')
#        end = 0
#        #For each regular expression, check if there is a match that is <10 AA from the end
#
#        if re.search(rex3,self.sequence) and len(re.split(rex3,self.sequence)[-1]) > 10:
#            start, end = [m.span() for m in rex3.finditer(self.sequence)][-1]
#        elif re.search(rex4,self.sequence) and len(re.split(rex4,self.sequence)[-1]) > 10:
#            start, end = [m.span() for m in rex4.finditer(self.sequence)][-1]
#        elif re.search(rex5,self.sequence) and len(re.split(rex5,self.sequence)[-1]) > 10:
#            start, end = [m.span() for m in rex5.finditer(self.sequence)][-1]
#        elif re.search(rex1,self.sequence) and len(re.split(rex1,self.sequence)[-1]) > 10:
#            start, end = [m.span() for m in rex1.finditer(self.sequence)][-1]
#            end -= 5
#        elif re.search(rex2,self.sequence) and len(re.split(rex2,self.sequence)[-1]) > 10:
#            start, end = [m.span() for m in rex2.finditer(self.sequence)][-1]
#            end -= 5 
#        self.split_index = end
#        self.leader = self.sequence[:end]
#        self.core = self.sequence[end:]
#        if len(self.leader) < 5 or len(self.core) < 5: #TAG CHRIS
#            self.valid_split = False
#            self.split = int(.25*len(self.sequence))
#            self.leader = self.sequence[:self.split]
#            self.core = self.sequence[self.split:]
        scores = [(1,int(.25*len(self.sequence)))]*3
        
        # switch the # on the two lines below if your computer username has a space in it. # RAL
        #fimo_output = self.run_fimo_simple("{}/berninamycin_fimo.txt".format("./ripp_modules/thio"))
        fimo_output = self.run_fimo_simple("{}/berninamycin_fimo.txt".format(FILE_DIR))
        fimo_output = fimo_output.split('\n')
        valid_split = False
        if len(fimo_output) > 1:
            for line in fimo_output[1:]:
                line = line.split('\t')
                if len(line) <= 1:
                    continue
                if float(line[7]) < scores[0][0]:
                    scores[0] = (float(line[7]), int(line[4]))
            valid_split = True  
        # switch the # on the two lines below if your computer username has a space in it.     # RAL
        #fimo_output = self.run_fimo_simple("{}/thio_fimo.txt".format("./ripp_modules/thio")).split('\n')
        fimo_output = self.run_fimo_simple("{}/thio_fimo.txt".format(FILE_DIR)).split('\n')
        if len(fimo_output) > 1:
            for line in fimo_output[1:]:
                line = line.split('\t')
                if len(line) <= 1:
                    continue
                if float(line[7]) < scores[1][0]:
                    scores[1] = (float(line[7]), int(line[4]))
            valid_split = True    
        # switch the # on the two lines below if your computer username has a space in it. # RAL
        #fimo_output = self.run_fimo_simple("{}/dhpip_fimo.txt".format("./ripp_modules/thio")).split('\n')
        fimo_output = self.run_fimo_simple("{}/dhpip_fimo.txt".format(FILE_DIR)).split('\n')
        if len(fimo_output) > 1:
            for line in fimo_output[1:]:
                line = line.split('\t')
                if len(line) <= 1:
                    continue
                if float(line[7]) < scores[2][0]:
                    scores[2] = (float(line[7]), int(line[4]))
            valid_split = True            
        scores = sorted(scores)
        self.valid_split = valid_split
        self.split = scores[0][1]
        if self.split < 4 or self.split > len(self.sequence) - 4:
            self.split = int(.25*len(self.sequence))
            self.valid_split = False
        self.leader = self.sequence[:self.split]
        self.core = self.sequence[self.split:]
            
        
        
    def set_score(self, pfam_hmm, cust_hmm):
        tabs = []
        score = 0
        pfams = []
        for pfam_dot in self.pfam_2_coords.keys():
            pfams.append(pfam_dot.split('.')[0])
        
        #Contains TOMM YcaO (PF02624)
        if "PF02624" in pfams:
            score += 2
            tabs.append(1)
        else:
            tabs.append(0)
        #TODO ask Chris why not just use pfam id... TAG Chris
        #Contains LanB N-terminal domain (PF04738)
        if "Lant_dehyd_N" in pfams or "PF04738" in pfams or "tsrC" in pfams:
            score += 2
            tabs.append(1)
        else:
            tabs.append(0)
        #Contains LanB C-terminal domain (PF14028)
        if "PF14028" in pfams:
            score += 2
            tabs.append(1)
        else:
            tabs.append(0)
        #Contains TOMM dehydrogenase (PF00881)
        if "PF00881" in pfams:
            score += 2
            tabs.append(1)
        else:
            tabs.append(0)
        #Contains rSAM methyltransferase (PF04055)
        if "PF04055" in pfams:
            score += 2
            tabs.append(1)
        else:
            tabs.append(0)
        #Contains P450 (PF00067)
        if "PF00067" in pfams:
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)
        #Contains ABC transporter or abhydrolase
        abc_transp_abhydrolases = ["PF00005", "PF01061", "PF12698", "PF12697", "PF00561"]
        for dom in abc_transp_abhydrolases:
            if dom in pfams:
                score += 1
                tabs.append(1)
            else:
                tabs.append(0)
        #CSS/CTT, SS/SSS/SSS, CC/CCC/CCCC, TT/TT/TTTT motifs
        motifs = (('[C][S]{2,}', 1), ('[C][T]{2,}', 1), ('[S]{2,}', 1), ('[S]{3,}', 1), \
           ('[S]{4,}', 2), ('[C]{2,}', 1), ('[C]{3,}', 1), ('[C]{4,}', 2), \
           ('[T]{2,}', 1), ('[T]{3,}', 1), ('[T]{4,}', 2))
        for motif in motifs:
            if re.search(motif[0], self.core) != None:
                score += motif[1]
                tabs.append(1)
            else:
                tabs.append(0)
        #No Cys/Ser/Thr core residues
        for aa in "CST":
            if aa not in self.core:
                score -= 2
                tabs.append(1)
            else:
                tabs.append(0)    
        #Mass of core peptide (unmodified) < 2100
        if "X" in self.core:
            Xs = self.core.count("X")
            noXcore = self.core.replace("X", "") #Remove Xs, as they do not work with molecular weight calculation
            core_analysis = ProteinAnalysis(noXcore, monoisotopic=True)
        else:
            core_analysis = ProteinAnalysis(self.core, monoisotopic=True)
        if core_analysis.molecular_weight() < 2100:
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)
        #Sum of repeating Cys/Ser/Thr > 4
        number_of_repeating_CST, number_of_repeat_blocks, avg_heteroblock_length, number_of_heteroblocks = thioscout(self.core)
        if sum([int(nr) for nr in number_of_repeating_CST.split(", ")]) > 4:
            score += 2
            tabs.append(1)
        else:
            tabs.append(0)
        #Avg heterocycle block length > 3
        if type(avg_heteroblock_length) == float and avg_heteroblock_length > 3:
            score += 2
            tabs.append(1)
        else:
            tabs.append(0)
        #Leader net charge < 5
        charge_dict = {"E": -1, "D": -1, "K": 1, "R": 1}
        leader_charge = sum([charge_dict[aa] for aa in self.leader if aa in charge_dict])
        if leader_charge < 5:
            score += 2
            tabs.append(1)
        else:
            tabs.append(0)
        #Leader net charge > 0
        if leader_charge > 0:
            score -= 2
            tabs.append(1)
        else:
            tabs.append(0)
        #Leader contains a Cys
        if "C" in self.leader:
            score -= 1
            tabs.append(1)
        else:
            tabs.append(0)
        #Peptide terminates Cys/Ser/Thr
        try:
            if self.core[-1] in ["C", "S", "T"]:
                score += 1
                tabs.append(1)
            else:
                tabs.append(0)
        except:
            print(self.core)
            print(self.leader)
            print(self.sequence)
        #Core contains >= 2 positive residues
        if sum([self.core.count(aa) for aa in "RK"]) >= 2:
            score -= 1
            tabs.append(1)
        else:
            tabs.append(0)
        #Number of heterocyclizable residues to core ratio > 0.4
        if float(sum([self.core.count(aa) for aa in "CST"])) / len(self.core) >= 0.4:
            score += 2
            tabs.append(1)
        else:
            tabs.append(0)
    
        #Now svm heuristics
        columns = []
        #append classification and index
        columns += tabs
        #Number repeating blocks of heterocyclizable residues in core
        number_of_repeating_CST, number_of_repeat_blocks, avg_heteroblock_length, number_of_heteroblocks = thioscout(self.core)
        columns.append(number_of_repeat_blocks)
        #Number of core repeating Cys
        columns.append(int(number_of_repeating_CST.split(", ")[0]))
        #Number of core repeating Ser
        columns.append(int(number_of_repeating_CST.split(", ")[1]))
        #Number of core repeating Thr
        columns.append(int(number_of_repeating_CST.split(", ")[2]))
        #Number of blocks of heterocyclizable residues in core
        columns.append(number_of_heteroblocks)
        #Average core heterocycle block length
        if avg_heteroblock_length == "nan":
            columns.append(0)
        else:
            columns.append(avg_heteroblock_length)
        #Precursor peptide mass (unmodified)
        if "X" in self.sequence:
            Xs = self.sequence.count("X")
            noXprecursor = self.sequence.replace("X", "") #Remove Xs, as they do not work with molecular weight calculation
            precursor_analysis = ProteinAnalysis(noXprecursor, monoisotopic=True)
            columns.append(float(precursor_analysis.molecular_weight()) + 110 * Xs)
        else:
            precursor_analysis = ProteinAnalysis(self.sequence, monoisotopic=True)
            columns.append(float(precursor_analysis.molecular_weight()))
        #Unmodified leader peptide mass
        if "X" in self.leader:
            Xs = self.leader.count("X")
            noXleader = self.leader.replace("X", "") #Remove Xs, as they do not work with molecular weight calculation
            leader_analysis = ProteinAnalysis(noXleader, monoisotopic=True)
            columns.append(float(leader_analysis.molecular_weight()) + 110 * Xs)
        else:
            leader_analysis = ProteinAnalysis(self.leader, monoisotopic=True)
            columns.append(float(leader_analysis.molecular_weight()))
        #Unmodified core peptide mass
        if "X" in self.core:
            Xs = self.core.count("X")
            noXcore = self.core.replace("X", "") #Remove Xs, as they do not work with molecular weight calculation
            core_analysis = ProteinAnalysis(noXcore, monoisotopic=True)
            columns.append(float(core_analysis.molecular_weight()) + 110 * Xs)
        else:
            core_analysis = ProteinAnalysis(self.core, monoisotopic=True)
            columns.append(float(core_analysis.molecular_weight()))
        #Length of Precursor
        columns.append(len(self.sequence))
        #Length of Leader
        columns.append(len(self.leader))
        #Length of Core
        columns.append(len(self.core))
        #Ratio of length of leader / length of core
        columns.append(float(len(self.core)) / float(len(self.leader)))
 
        #Ratio of heterocyclizable  residues / length of core
        columns.append(float(sum([self.core.count(aa) for aa in "CST"])) / len(self.core))
        #Estimated core charge at neutral pH
        #charge_dict = {"E": -1, "D": -1, "K": 1, "H": 1, "R": 1}
        #columns.append(sum([charge_dict[aa] for aa in core if aa in charge_dict]))
        #Estimated leader charge at neutral pH
        #columns.append(sum([charge_dict[aa] for aa in leader if aa in charge_dict]))
        #Estimated precursor charge at neutral pH
        #columns.append(sum([charge_dict[aa] for aa in leader+core if aa in charge_dict]))
        #Absolute value of core charge at neutral pH
        #columns.append(abs(sum([charge_dict[aa] for aa in core if aa in charge_dict])))
        #Absolute value of leader charge at neutral pH
        #columns.append(abs(sum([charge_dict[aa] for aa in leader if aa in charge_dict])))
        #Absolute value of precursor charge at neutral pH
        #columns.append(abs(sum([charge_dict[aa] for aa in leader+core if aa in charge_dict])))
        #Number in leader of each amino acid
        columns += [self.leader.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
        #Aromatics in leader
        columns.append(sum([self.leader.count(aa) for aa in "FWY"]))
        #Neg charged in leader
        columns.append(sum([self.leader.count(aa) for aa in "DE"]))
        #Pos charged in leader
        columns.append(sum([self.leader.count(aa) for aa in "RK"]))
        #Charged in leader
        columns.append(sum([self.leader.count(aa) for aa in "RKDE"]))
        #Aliphatic in leader
        columns.append(sum([self.leader.count(aa) for aa in "GAVLMI"]))
        #Hydroxyl in leader
        columns.append(sum([self.leader.count(aa) for aa in "ST"]))
        #Counts of AAs in core
        columns += [self.core.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
        #Aromatics in core
        columns.append(sum([self.core.count(aa) for aa in "FWY"]))
        #Neg charged in core
        columns.append(sum([self.core.count(aa) for aa in "DE"]))
        #Pos charged in core
        columns.append(sum([self.core.count(aa) for aa in "RK"]))
        #Charged in core
        columns.append(sum([self.core.count(aa) for aa in "RKDE"]))
        #Aliphatic in core
        columns.append(sum([self.core.count(aa) for aa in "GAVLMI"]))
        #Hydroxyl in core
        columns.append(sum([self.core.count(aa) for aa in "ST"]))
        #Counts of AAs in entire precursor (leader+core)
        columns += [self.sequence.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
        #Aromatics in precursor
        columns.append(sum([self.sequence.count(aa) for aa in "FWY"]))
        #Neg charged in precursor
        columns.append(sum([self.sequence.count(aa) for aa in "DE"]))
        #Pos charged in precursor
        columns.append(sum([self.sequence.count(aa) for aa in "RK"]))
        #Charged in precursor
        columns.append(sum([self.sequence.count(aa) for aa in "RKDE"]))
        #Aliphatic in precursor
        columns.append(sum([self.sequence.count(aa) for aa in "GAVLMI"]))
        #Hydroxyl in precursor
        columns.append(sum([self.sequence.count(aa) for aa in "ST"]))
        self.score = score
        self.csv_columns += [score] + columns
        return
        
        
def thioscout(core):
    """ThioScout function from Chris Schwalen to count repeat blocks"""
    #rex1 repeating Cys Ser Thr residues
    rex1 = re.compile('[C]{2,}|[S]{2,}|[T]{2,}')
    #rex2 contiguous cyclizable residues
    rex2 = re.compile('[C|S|T]{2,}')
    rexout1 = re.findall(rex1,core)
    number_of_repeat_blocks = len(rexout1)
    temp = "".join(rexout1)
    number_of_repeating_CST = str([temp.count("C"), temp.count("S"), temp.count("T")]).strip("[]")
    rexout2 = re.findall(rex2,core)
    number_of_heteroblocks = len(rexout2)
    if len(rexout2) == 0:
        rexout2 = np.nan #TODO ask chris what to do here TAG Chris
        avg_heteroblock_length = 0
    else:
        rexout2 = np.mean([len(x) for x in rexout2])
        avg_heteroblock_length = str(rexout2).strip("[]")
    return number_of_repeating_CST, number_of_repeat_blocks, avg_heteroblock_length, number_of_heteroblocks
    
