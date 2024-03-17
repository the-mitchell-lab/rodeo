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
import hmmer_utils
import pathlib
FILE_DIR = pathlib.Path(__file__).parent.absolute()

peptide_type = "linar"
CUTOFF = 12
index = 0

def write_csv_headers(output_dir, meta=False):
    dir_prefix = output_dir + '/linar/'
    if not os.path.exists(dir_prefix):
        os.makedirs(dir_prefix)
    svm_headers = 'PK,Classification,Contains ABC Transporter,Contains n-methyl transferase,Contains flavin decarboxylase,Has CXXC motif and flavin decarboxylase in BGC,Has GST motif and flavin decarboxylase,Gene cluster and precursor in same direction,Charge of leader <1 at pH 7,Core contains >2 Cysteine,Leader contains 0 Cysteine,Core % Aliphatic and Dhb residues,Core begins with XTP motif,Leader contains GXG motif,Leader contains LXD motif,Leader contains FAN motif,Distance from LinG homolog,Length of precursor,Leader Percent A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,% Aromatics,% Negative,% Positive,% Charged,% Aliphatic,% Hydroxyl,Core percent A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,% Aromatics,% Negative,% Positive,% Charged,% Aliphatic,% Hydroxyl,Relative charge of core,Relative charge of leader,Relative charge of precursor,Absolute charge of core,Absolute charge of Leader,Absolute charge of precursor' #,Motif 1,2,3,4,5,6,7,8,Number of motifs hit,Motif score 1,2,3,4,5,6,7,8,Total motif score,No motifs hit'
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
        self.peptide_type = 'linar'
        self.set_split()
#        self.set_monoisotopic_mass()
        self.csv_columns = [self.leader, self.core, self.start, self.end]
        self.CUTOFF = CUTOFF
        
    def set_split(self):
        
        # Linaridin cleavage site set by scanning FIMO identified motifs from linar_cutsites.
        # If no motif is identified, then regular expressions are used to predict the cleavage site.
        # If still no cleavage site is identified, a default cleavage site halfway through the sequence is selected.
        score = [1, int(.50*len(self.sequence))]
        # switch the # on the two lines below if your computer username has a space in it.     # RAL
        #fimo_output = self.run_fimo_simple("{}/linar_cutsites.txt".format("./ripp_modules/linar"))
        fimo_output = self.run_fimo_simple("{}/linar_cutsites.txt".format(FILE_DIR))
        fimo_output = fimo_output.split('\n')
        valid_split = False
        if len(fimo_output) > 1:
            for line in fimo_output[1:]:
                line = line.split('\t')
                if len(line) <= 1:
                    continue
                if float(line[7]) < score[0]:
                        score = [float(line[7]), int(line[4])-4]
            valid_split = True

        if valid_split == False:
            core_starts = [".TP", ".T[A|V|L|I|F|T]", ".[A|V|L|I|P|C]P"]     # These motifs were derived from analysis of predicted linaridin sequences from other   
            leader_ends = ["P..", "[A|V|L|I|F].."]                          # previous datasets and refined as further linaridins were identified with this tool

            end = 0
            valid_split = False
            for leader_end in leader_ends:
                for core_start in core_starts:

                    # Checks for all iterations of cleavage site at least 10 aa from either terminus
                    if re.search(leader_end + core_start, self.sequence) and len(re.split(leader_end + core_start, self.sequence)[-1]) > 10:
                        rex = re.compile(leader_end + core_start)
                        for m in rex.finditer(self.sequence):
                            _, score[1] = m.span()
                        score[1] -= 3
                        valid_split = True
                        break
                if valid_split == True:
                    break

        self.split = score[1]
        self.leader = self.sequence[:self.split]
        self.core = self.sequence[self.split:]
        if len(self.leader) < 5 or len(self.core) < 5: 
            valid_split = False

        self.valid_split = valid_split
        self.split = score[1] 
        if self.split < 10 or self.split > len(self.sequence) - 10:
            self.split = int(.5*len(self.sequence))
            self.valid_split = False
        self.leader = self.sequence[:self.split]
        self.core = self.sequence[self.split:]
            
    def get_fimo_score(self):
        fimo_output = str(self.run_fimo_simple())
        fimo_motifs = []
        fimo_motifs = [int(line.split("\t")[1].partition("-")[2]) for line in fimo_output.split("\n") if "\t" in line and "motif" not in line]
        fimo_scores = {int(line.split("\t")[1].partition("-")[2]): float(line.split("\t")[6]) for line in fimo_output.split("\n") if "\t" in line and "motif" not in line}

        #Calculate score
        motif_score = 0
        if any in fimo_motifs:
            motif_score += 2
        else:
            motif_score -= 1

        return fimo_motifs, motif_score, fimo_scores    
        
    def set_score(self, pfam_hmm, cust_hmm):
        tabs = []
        score = 0
        pfams = []
        for pfam_dot in self.pfam_2_coords.keys():
            pfams.append(pfam_dot.split('.')[0])

        #Contains abc transporters 
        abc_transporters = ["PF00005", "PF02163", "PF00664"]
        if set(abc_transporters).intersection(pfams):
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)

        #Gene cluster contains an n-methyl transferase
        n_methyl_transf = ["PF08241", "PF08242", "PF13489", "PF13649", "PF13847"]
        if set(n_methyl_transf).intersection(pfams):
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)

        #Gene cluster contains a Flavin decarboxylase (save flavin decarboylase for other scoring)
        flavin_decarboxylase = ["TIGR00521", "TIGR02113", "PF02441"]
        flavin_decarboxylase_present = 0
        if set(flavin_decarboxylase).intersection(pfams):
            score += 1
            flavin_decarboxylase_present += 1
            tabs.append(1)
        else:
            tabs.append(0)

        #core has CXXC motif and a flavin decarboxylase is in the BGC
        CXXC = re.compile("C[A-Z]{2}C")
        if CXXC.search(self.core) != None and flavin_decarboxylase_present == 1:
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)

        #precursor has GST motif and a flavin decarboxylase is in the BGC
        GST = re.compile("GST")
        if GST.search(self.sequence) != None and flavin_decarboxylase_present == 1:
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)

        #Gene cluster and precursor in same direction
        if "LinE" not in self.pfam_2_coords.keys():
            tabs.append(0)
        elif np.sign(self.start - self.end) == np.sign(self.pfam_2_coords["LinE"][0][0]-self.pfam_2_coords["LinE"][0][1]):
            score += 2
            tabs.append(1)
        else:
            tabs.append(0)

        #Positive hit for LinA.hmm file at e-value < 0.001
        # switch the # on the two lines below if your computer username has a space in it.     # RAL
        #precursor_pfam = hmmer_utils.get_hmmer_info(self.sequence, "{}/hmms/LinA.hmm".format("./ripp_modules/linar"), "", n=1)
        precursor_pfam = hmmer_utils.get_hmmer_info(self.sequence, "{}/hmms/LinA.hmm".format(FILE_DIR), "", self.output_dir, n=1)
        if precursor_pfam and len(self.sequence) < 120:
            score += 10

        #Adding distance heuristic scoring
        dist = 1000000
        for pfam in self.pfam_2_coords.keys():
            for fam in ["LinE", "LinG", "LinL"]:
                if fam in pfam[0]:
                    dist = min(self.get_min_dist(self.pfam_2_coords[pfam[0]]), dist)
        if dist < 300:
            score += 1
        if dist < 600:
            score += 1
        if dist > 2200:
            score -= 1

        #Length of Precursor is between 50 and 70 aa
        if 50 <= len(self.sequence) <= 70:
            score += 2
        if len(self.sequence) > 120:
            score -= 5
        if len(self.sequence) > 150:
            score -= 5

        #Length of precursor is between 71 and 100 aa
        if 71 <= len(self.sequence) <= 100:
            score += 1

        #Charge of leader <1 at pH 7
        charge_dict = {"E": -1, "D": -1, "K": 1, "R": 1}
        leader_charge = sum([charge_dict[aa] for aa in self.leader if aa in charge_dict])
        if leader_charge < 1:
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)

        #Charge of core has value greater than 6 or less then -6 at pH 7
        core_charge = sum([charge_dict[aa] for aa in self.core if aa in charge_dict])
        if abs(core_charge) > 6:
            score -= 1

        #Core is >15% alanine
        if float(self.core.count("A"))/len(self.core) > .15:
            score += 1

        #Core is >20% alanine
        if float(self.core.count("A"))/len(self.core) > .2:
            score += 1

        #Core is >10% Threonine
        if float(self.core.count("T"))/len(self.core) > .1:
            score += 1

        #Core is >15% Threonine
        if float(self.core.count("T"))/len(self.core) > .15:
            score += 1

        #Core contains 0 threonines
        if not "T" in self.core:
            score -= 2

        #Core contains >7% Valine
        if float(self.core.count("V"))/len(self.core) > .07:
            score += 1

        #Core contains >12% Valine
        if float(self.core.count("V"))/len(self.core) > .12:
            score += 1

        #Core contains >2 Cysteine
        if self.core.count("C") > 2:
            score -= 1
            tabs.append(1)
        else:
            tabs.append(0)

        #Leader contains 0 Cysteine
        if "C" in self.leader:
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)

        #Core Aliphatic + Dhb < 50%
        tabs.append(float(sum([self.core.count(aa) for aa in "GAVLMIT"]))/len(self.core))
        if sum([self.core.count(aa) for aa in "GAVLMIT"])/len(self.core) < 0.5:
            score -= 1

        #Core Aliphatic + Dhb > 60%
        if sum([self.core.count(aa) for aa in "GAVLMIT"])/len(self.core) > 0.6:
            score += 1

        #Core Aliphatic + Dhb > 75%
        if sum([self.core.count(aa) for aa in "GAVLMIT"])/len(self.core) > 0.75:
            score += 1

        #Core begins "XTP"
        XTP = re.compile(".TP")
        if XTP.match(self.core):
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)

        #Leader contains "GXG"
        GXG = re.compile("G.G")
        if GXG.search(self.leader):
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)

        #Leader contains "LXD"
        LXD = re.compile("L.D")
        if LXD.search(self.leader):
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)

        #Leader contains "FAN"
        FAN = re.compile("FAN")
        if FAN.search(self.leader):
            score += 1
            tabs.append(1)
        else:
            tabs.append(0)


        #Now svm heuristics
        columns = []
        #append classification and index
        columns += tabs
        #Minimum distance from LinG homolog
        dist = 10000
        for pfam in self.pfam_2_coords.keys():
            if "LinG" in pfam:
                dist = self.get_min_dist(self.pfam_2_coords[pfam])
        columns.append(dist)
        #Length of Precursor
        columns.append(len(self.sequence))

        # % in leader of each amino acid
        leader_length = len(self.leader)
        columns += [float(self.leader.count(aa))/leader_length for aa in "ARDNCQEGHILKMFPSTWYV"]
        # % Aromatics in leader
        columns.append(sum([float(self.leader.count(aa))/leader_length for aa in "FWY"]))
        # % Neg charged in leader
        columns.append(sum([float(self.leader.count(aa))/leader_length for aa in "DE"]))
        # % Pos charged in leader
        columns.append(sum([float(self.leader.count(aa))/leader_length for aa in "RK"]))
        # % Charged in leader
        columns.append(sum([float(self.leader.count(aa))/leader_length for aa in "RKDE"]))
        # % Aliphatic in leader
        columns.append(sum([float(self.leader.count(aa))/leader_length for aa in "GAVLMI"]))
        # % Hydroxyl in leader
        columns.append(sum([float(self.leader.count(aa))/leader_length for aa in "ST"]))

        # % of AAs in core
        core_length = len(self.core)
        columns += [float(self.core.count(aa))/core_length for aa in "ARDNCQEGHILKMFPSTWYV"]
        # % Aromatics in core
        columns.append(sum([float(self.core.count(aa))/core_length for aa in "FWY"]))
        # % Neg charged in core
        columns.append(sum([float(self.core.count(aa))/core_length for aa in "DE"]))
        # % Pos charged in core
        columns.append(sum([float(self.core.count(aa))/core_length for aa in "RK"]))
        # % Charged in core
        columns.append(sum([float(self.core.count(aa))/core_length for aa in "RKDE"]))
        # % Aliphatic in core
        columns.append(sum([float(self.core.count(aa))/core_length for aa in "GAVLMI"]))
        # % Hydroxyl in core
        columns.append(sum([float(self.core.count(aa))/core_length for aa in "ST"]))

        columns.append(sum([charge_dict[aa] for aa in self.core if aa in charge_dict]))
        columns.append(sum([charge_dict[aa] for aa in self.leader if aa in charge_dict]))
        columns.append(sum([charge_dict[aa] for aa in self.sequence if aa in charge_dict]))
        columns.append(abs(sum([charge_dict[aa] for aa in self.core if aa in charge_dict])))
        columns.append(abs(sum([charge_dict[aa] for aa in self.leader if aa in charge_dict])))
        columns.append(abs(sum([charge_dict[aa] for aa in self.sequence if aa in charge_dict])))

        fimo_motifs, motif_score, fimo_scores = self.get_fimo_score()
        self.fimo_motifs = fimo_motifs
        self.fimo_scores = fimo_scores
        self.score += motif_score
        self.score = score
        self.csv_columns += [score] + columns
        return
