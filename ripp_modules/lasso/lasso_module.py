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
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from ripp_modules.VirtualRipp import VirtualRipp
import pathlib
FILE_DIR = pathlib.Path(__file__).parent.absolute()

peptide_type = "lasso"
CUTOFF = 13
index = 0


def write_csv_headers(output_dir):
    dir_prefix = output_dir + '/lasso/'
    if not os.path.exists(dir_prefix):
        os.makedirs(dir_prefix)
    svm_headers = 'Precursor Index,classification,Calcd. Lasso Mass (Da), Distance,Within 500 nt?,Within 150 nt?,Further than 1000 nt?,Core has 2 or 4 Cys?,Leader longer than core?,Lasso ring length (-1 if not plausible),Leader has GxxxxxT motif?,Core starts with G?,Ratio leader/core < 2 and > 0.5,Core starts with Cys and even number of Cys?,No Gly in core?,Core has at least 1 aromatic aa?,Core has at least 2 aromatic aa?,Core has odd number of Cys?,Leader has Trp?,Leader has Lys?,Leader has Cys?,Cluster has PF00733?,Cluster has PF05402?,Cluster has PF13471?,Leader has LxxxxxT motif?,Core has adjacent identical aas (doubles)?,Core length (aa),Leader length (aa),Precursor length (aa),Leader/core ratio,Number of Gly in first 9 aa of core?,Number of Pro in first 9 aa of core?,CORE starts with A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,Estimated core charge,Estimated leader charge,Estimated precursor charge,Absolute value of core charge,Absolute value of leader charge,Absolute value of precursor charge,LEADER A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,Aromatics,Neg charged,Pos charged,Charged,Aliphatic,Hydroxyl,CORE A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,Aromatics,Aromatics (last 10aa),Neg charged,Pos charged,Charged,Aliphatic,Hydroxyl,FIRST CORE RESIDUE A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V, PRECURSOR A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,Aromatics,Neg charged,Pos charged,Charged,Aliphatic,Hydroxyl,Motif1?,Motif2?,Motif3?,Motif4?,Motif5?,Motif6?,Motif7?,Motif8?,Motif9?,Motif10?,Motif11?,Motif12?,Motif13?,Motif14?,Motif15?,Motif16?,Total motifs hit,"Score Motif1","Score Motif2","Score Motif3","Score Motif4","Score Motif5","Score Motif6","Score Motif7","Score Motif8","Score Motif9","Score Motif10","Score Motif11","Score Motif12","Score Motif13","Score Motif14","Score Motif15","Score Motif16",Sum of MEME scores,No Motifs?,Alternate Start Codon?'
    svm_headers = svm_headers.split(',')
    features_headers = ['Accession_id', 'Genus/Species', 'Leader', 'Core', 'Start', 'End' , "Total Score", "Valid Precursor"] + svm_headers
    features_csv_file = open(dir_prefix + "temp_features.csv", 'w')
    svm_csv_file = open(os.path.join(FILE_DIR, "svm/fitting_set.csv"), 'w')
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
        self.peptide_type = 'lasso'
        self.set_split()
        self.set_monoisotopic_mass()
        self.csv_columns = [self.leader, self.core, self.start, self.end]
        self.CUTOFF = CUTOFF


    def set_split(self):
        scores = [(1,int(.25*len(self.sequence)))]*3
        fimo_output = self.run_fimo_simple("{}/wxxp_fimo.txt".format(FILE_DIR))
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
        fimo_output = self.run_fimo_simple("{}/yxxp_fimo.txt".format(FILE_DIR)).split('\n')
        if len(fimo_output) > 1:
            for line in fimo_output[1:]:
                line = line.split('\t')
                if len(line) <= 1:
                    continue
                if float(line[7]) < scores[1][0]:
                    scores[1] = (float(line[7]), int(line[4]) + (1 if line[1] in ("MEME-2", "MEME-4") else 0))
            valid_split = True
        fimo_output = self.run_fimo_simple("{}/zoops20.txt".format(FILE_DIR)).split('\n')
        if len(fimo_output) > 1:
            for line in fimo_output[1:]:
                line = line.split('\t')
                if len(line) <= 1 or line[1] not in ("MEME-1", "MEME-2"):
                    continue
                if float(line[7]) < scores[2][0]:
                    scores[2] = (float(line[7]), int(line[4]))
            valid_split = True
        if valid_split and scores[2][0] == min(score for score, pos in scores) and 4 <= scores[2][1] <= len(self.sequence) - 4:
            pass
            # print("Split fused!")
        scores = sorted(scores)
        self.valid_split = valid_split
        self.split = scores[0][1]
        if self.split < 4 or self.split > len(self.sequence) - 4:
            self.split = int(.25*len(self.sequence))
            self.valid_split = False
        self.leader = self.sequence[:self.split]
        self.core = self.sequence[self.split:]
        if not self.valid_split:
            self.set_split2()



    def set_split2(self):
        #TODO add more regexes
        #change initial check to include L/V/I rather than only L
        #added YxxP instances without L/G but remaining in frame
        #remove ExxP and added WxxP in the three possible registers - could make it two seperate lines if prioritizing W in -15
        #added the new recognition site - most common LxxLG xxxxx T x but L can also be I/V others too but keep it simple for now
        #change the loop-tail size to prioritze size 5-15 the most common sizes amoung known lassos then 15-25 25 being the largest known then 26-30
        match = re.search('(Y..P.(L|V|I)...G.....T)',self.sequence)
        motif = 1
        if match is None:
            match = re.search('(Y..P...........T)', self.sequence)
            motif = 1
        if match is None:
            match = re.search('(W..P.[A-Z]{1-3}.......T)', self.sequence)
            motif = 1
        if match is None:
            match = re.search('((L|V|I)..(L|V|I)G.....T)', self.sequence)
            motif = 1
        if match is None:
            match = re.search('(Y..P.(L|V|I)...G.....(S|V|A))', self.sequence)
            motif = 1
        if match is None:
            motif = 2
            match = re.search('(T[A-Z]{7,10}(D|E)[A-Z]{5,15}\*)', self.sequence + '*')
            if match is None:
                match = re.search('(T[A-Z]{7,10}(D|E)[A-Z]{15,25}\*)', self.sequence + '*')
            if match is None:
                match = re.search('(T[A-Z]{7,10}(D|E)[A-Z]{26,30}\*)', self.sequence + '*')
            if match is None:
                match = re.search('(T[A-Z]{10,18}(D|E)[A-Z]{5,15}\*)', self.sequence + '*')
            if match is None:
                match = re.search('(T[A-Z]{10,18}(D|E)[A-Z]{15,25}\*)', self.sequence + '*')
            if match is None:
                match = re.search('(T[A-Z]{10,18}(D|E)[A-Z]{26,30}\*)', self.sequence + '*')
        if match is not None:
            if motif == 1:
                self.split_index = match.end() + 1
            elif motif == 2:
                self.split_index = match.start() + 2
        else:
            self.split_index = -1
        if self.split_index == -1 or abs(len(self.sequence)-self.split_index) < 5:
            self.valid_split = False
            self.split_index = int(.5*len(self.sequence))

        self.leader = self.sequence[0:self.split_index]
        self.core = self.sequence[self.split_index:]


    def get_fimo_score(self):
        fimo_output = str(self.run_fimo_simple())
        fimo_motifs = []
        fimo_motifs = [int(line.partition("\t")[0]) for line in fimo_output.split("\n") if "\t" in line and line.partition("\t")[0].isdigit()]
        fimo_scores = {int(line.split("\t")[0]): float(line.split("\t")[6]) for line in fimo_output.split("\n") if "\t" in line and line.partition("\t")[0].isdigit()}
        #Calculate score
        motif_score = 0
        if 2 in fimo_motifs:
            motif_score += 4
        elif len(fimo_motifs) > 0:
            motif_score += 2
        else:
            motif_score += -1
        return fimo_motifs, motif_score, fimo_scores

    def set_monoisotopic_mass(self):
        self._set_number_bridges()
        CC_mass = 2*self._num_bridges
        # dehydration indicative of cyclization
        bond = 18.02
#        print(self.core)
        if "B" in self.core:
            print(self.core)
            print("AA sequence contains Asx. Currently reviewing how to assign weights to such translations.")
        if "J" in self.core:
            print(self.core)
            print("AA sequence contains 'J'. Currently reviewing how to assign weights to such translations.")
        if "Z" in self.core:
            print(self.core)
            print("AA sequence contains 'Z'. Currently reviewing how to assign weights to such translations.")
        try:
            monoisotopic_mass = ProteinAnalysis(self.core.replace('X', '').replace('B', '').replace('J', '').replace('Z', ''), monoisotopic=True).molecular_weight()
        except ValueError:
            print("ERROR assigning molecular mass to\n {} \n Assigning a mass of 0.".format(self.sequence))

        self._monoisotopic_weight = monoisotopic_mass + CC_mass - bond

    def _set_number_bridges(self):
        '''
        Predict the lassopeptide number of disulfide bridges
        '''
        self._num_bridges = 0
        if self.core.count("C") == 2:
            self._num_bridges = 1
        if self.core.count("C") >= 4:
            self._num_bridges = 2
        return self._num_bridges

    def set_score(self, pfam_dir, cust_hmm):
        scoring_csv_columns = []

        if not (25 <= len(self.sequence) <= 80):
            self.score -= 2
        if len(self.sequence) > 100:
            self.score -= 5
        if len(self.sequence) > 150:
            self.score -= 10

        if self.leader[-2] == "T":
            self.score += 4

        if any(c == self.leader[-2] for c in "VIMS"):
            self.score += 1

        if "P" == self.leader[-1]:
            self.score -= 4
        if "P" == self.core[0]:
            self.score -= 4



        if self.core[0] in "KRDE":
            self.score -= 1

        if sum(map(self.core[:10].count, ['G','P'])) > sum(map(self.core[-10:].count, ["G","P"])):
            self.score += 0

        matches = list(re.finditer('[YW]..P', self.leader))
        if matches:
            self.score += 2
            for m in matches:
                if len(self.leader) - m.start() in [-17, -16, -15]:
                    self.score += 3
                    break


        cPFam = "PF00733"
        ePFam = "PF05402"
        bPFam = "PF13471"
        self.score = 0
        scoring_csv_columns.append(self._monoisotopic_weight)

        bsp_coords = []
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in [cPFam, ePFam, bPFam]):
                bsp_coords += self.pfam_2_coords[pfam]
        min_distance = self.get_min_dist(bsp_coords)
        if min_distance is None:
            scoring_csv_columns.append(999999)
        else:
            scoring_csv_columns.append(min_distance)

        has_bPFam = False
        has_cPFam = False
        has_ePFam = False
        within_500 = False
        within_150 = False
        within_1000 = False
        cyclase_same_strand= False
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in [cPFam, ePFam, bPFam]):
                if bPFam in pfam:
                    has_bPFam = True
                if ePFam in pfam:
                    has_ePFam = True
                if cPFam in pfam:
                    if self.start < self.end and self.pfam_2_coords[pfam][0] < self.pfam_2_coords[pfam][0] or\
                       self.start > self.end and self.pfam_2_coords[pfam][0] > self.pfam_2_coords[pfam][0]:
                       cyclase_same_strand = True

                    has_cPFam = True
                dist = self.get_min_dist(self.pfam_2_coords[pfam])
                if dist < 1000:
                    within_1000 = True
                    if dist < 500:
                        within_500 = True
                        if dist < 150:
                            within_150 = True
        if within_500:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if within_150:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

        if not within_1000:
            self.score += -2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)


        if self.core.count("C") == 2 or self.core.count("C") == 4:
            self.score += 2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

        if len(self.leader) > len(self.core):
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

        # Ring length
        ts = []
        range_priorities = [range(8, 10), range(7, 8), range(10, 11), range(6, 7), range(11, 12)]
        for r in range_priorities:
            for c in "DE":
                r = list(r)
                t = self.core.find(c, r[0], r[-1] + 1)
                if t != -1:
                    break
            if t != -1:
                break


        if t in [8, 9]:
            self.score += 6
        elif t == 7:
            self.score += 4
        elif t == 10:
            self.score += 2
        elif t in [6, 11]:
            pass
        else:
            self.score += -6
        scoring_csv_columns.append(t)


        #Check for GxxxxxT motif
        match = re.search('G.{5}T', self.leader)
        if match is not None:
            self.score += 2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

        if self.core[0] == "G":
            self.score += 2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

        if "G" not in self.core[:12]:
            self.score -= 4

        #check if peptide and lasso cyclase are on same strand +1
        # if cyclase_same_strand:
            # self.score += 1
            # scoring_csv_columns.append(1)
        # else:
            # scoring_csv_columns.append(0)

        if 0.5 < len(self.leader)/float(len(self.core)) < 2:
            self.score += 2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        
        if self.core.count("C") % 2 == 1:
            self.score -= 2

        if self.core[0] == 'C' and self.core.count('C') % 2 == 0:
            self.score += 0
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

        if "G" not in self.core[:2]:
            self.score += -4
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Check for aromatic residues
        if any(aa in self.core for aa in ["F", "Y", "W"]):
            self.score += 2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

        if sum(aa in self.core for aa in ["F", "Y", "W"]) >= 2:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

        if self.core.count("C") % 2 == 1:
            # self.score += -2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

        if "W" in self.leader:
            # TODO "check for tryptophan"
            self.score += -1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if "K" in self.leader:
            self.score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if "C" in self.leader:
            # TODO "leader has cystine"
            self.score += -2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

        if has_cPFam:
            scoring_csv_columns.append(1)
        else:
            self.score += -5
            scoring_csv_columns.append(0)
        if has_ePFam:
            scoring_csv_columns.append(1)
        else:
            self.score += -5
            scoring_csv_columns.append(0)
        if has_bPFam:
            scoring_csv_columns.append(1)
        else:
            self.score += -5
            scoring_csv_columns.append(0)

        match = re.search('L.{5}T', self.leader)
        if match is not None:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

        #'Core has adjacent identical aas'
        prev_aa = ''
        adjacent_aas= False
        for aa in self.core:
            if aa == prev_aa:
                adjacent_aas = True
                break
            else:
                prev_aa = aa
        if adjacent_aas:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)

        scoring_csv_columns.append(len(self.core))
        scoring_csv_columns.append(len(self.leader))
        scoring_csv_columns.append(len(self.sequence))
        scoring_csv_columns.append(float(len(self.leader))/len(self.core))
        scoring_csv_columns.append(self.core[:9].count('G'))
        scoring_csv_columns.append(self.core[:9].count('P'))

        #Counts of AAs in leader
        scoring_csv_columns += [int(self.core[0] == aa) for aa in "ARDNCQEGHILKMFPSTWYV"]

        charge_dict = {"E": -1, "D": -1, "K": 1, "H": 1, "R": 1}
        scoring_csv_columns.append(sum([charge_dict[aa] for aa in self.core if aa in charge_dict]))
        #Estimated leader charge
        scoring_csv_columns.append(sum([charge_dict[aa] for aa in self.leader if aa in charge_dict]))
        #Estimated precursor charge
        scoring_csv_columns.append(sum([charge_dict[aa] for aa in self.sequence if aa in charge_dict]))
        #Absolute value of core charge
        scoring_csv_columns.append(abs(sum([charge_dict[aa] for aa in self.core if aa in charge_dict])))
        #Absolute value of leader charge
        scoring_csv_columns.append(abs(sum([charge_dict[aa] for aa in self.leader if aa in charge_dict])))
        #Absolute value of precursor charge
        scoring_csv_columns.append(abs(sum([charge_dict[aa] for aa in self.sequence if aa in charge_dict])))
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
        #Aromatics in last 10 of core
        scoring_csv_columns.append(sum([self.core[-10:].count(aa) for aa in "FWY"]))
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
        #Counts (0 or 1) of amino acids within first AA position of core sequence
        scoring_csv_columns += [self.core[0].count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
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


        fimo_motifs, motif_score, fimo_scores = self.get_fimo_score()
        self.fimo_motifs = fimo_motifs
        self.fimo_scores = fimo_scores
        self.score += motif_score
        #Motifs
        scoring_csv_columns += [1 if motif in fimo_motifs else 0 for motif in range(1, 17)]
        #Total motifs hit
        scoring_csv_columns.append(len(fimo_motifs))
        #Motif scores
        scoring_csv_columns += [fimo_scores[motif] if motif in fimo_motifs else 0 for motif in range(1, 17)]
        #Sum of MEME scores
        scoring_csv_columns.append(sum([fimo_scores[motif] if motif in fimo_motifs else 0 for motif in range(1, 17)]))
        #No Motifs?
        if len(fimo_motifs) == 0:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if self.leader[0] != 'M':
            scoring_csv_columns.append(1)
            self.score += -1
        else:
            scoring_csv_columns.append(0)
        self.csv_columns += [self.score] +  scoring_csv_columns

