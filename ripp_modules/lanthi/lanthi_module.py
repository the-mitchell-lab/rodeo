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
import hmmer_utils
import pathlib
FILE_DIR = pathlib.Path(__file__).parent.absolute()

peptide_type = "lanthi"
CUTOFF = 14
index = 0

def write_csv_headers(output_dir):
    dir_prefix = output_dir + '/lanthi/'
    if not os.path.exists(dir_prefix):
        os.makedirs(dir_prefix)
    svm_headers = 'PK,Classification,F?LD,S????C,T????C,S?????C,T?????C,Within 500 nt?,Cluster contains PF04738,Cluster contains PF05147,Cluster LACKS PF04738,Cluster LACKS PF05147,Cluster contains PF14028,Cluster contains PF00082,Cluster contains PF03412,Cluster contains PF00005,Cluster contains PF02624,Cluster contains PF00899,Cluster contains PF02052,Cluster contains PF08130,Precursor mass < 4000,Core mass < 2000,Peptide hits cl03420 (Gallidermin),Peptide hits TIGR03731 (lantibio_gallid),Peptide hits cl22812 (lanti_SCO0268),Peptide hits TIGR04363 (LD_lanti_pre),Peptide hits cl06940 (Antimicrobial18),Peptide hits PF02052 (Gallidermin),Peptide hits PF08130 (Antimicrobial18),Precursor peptide mass (unmodified),Leader peptide mass (unmodified),Core peptide mass (unmodified),Length of Leader,Length of Core,Length of precursor,Leader / core ratio,Core >= 35,Has repeating C motifs (not in last 3 residues),Leader > 4 neg charge motifs,Leader net neg charge,Leader FxLD,C-terminal CC,core DGCGxTC motif,core SFNS motif,core SxxLC motif,core CTxGC motif,core TPGC motif,core SFNS?C,Core Cys <3,Core Cys <2,No Core Cys residues,No Core Ser residues,No Core Thr residues,LS max >4,LS max <3,LS 4-membered ring >2,LS 5-membered ring >2,LS 6-membered ring,LS 7-membered ring,LS 8-membered ring,MEME/FIMO,LS max ring number,LS lan4,LS lan5,LS lan6,LS lan7,LS lan8,Ratio of Cys to sum of Ser/Thr,Ratio of Cys/Ser/Thr to len of core,Log10 MEME/FIMO score,log10 MEME motif 1,MEME motif 2,MEME motif 3,MEME motif 4,MEME motif 5,A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,Aromatics,Neg charged,Pos charged,Charged,Aliphatic,Hydroxyl,A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,Aromatics,Neg charged,Pos charged,Charged,Aliphatic,Hydroxyl,A,R,D,N,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,Aromatics,Neg charged,Pos charged,Charged,Aliphatic,Hydroxyl'  
    svm_headers = svm_headers.split(',')
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
                 pfam_2_evalue):
        super(Ripp, self).__init__(start, 
                                     end, 
                                     sequence,
                                     upstream_sequence,
                                     pfam_2_coords,
                                     pfam_2_evalue)
        self.peptide_type = 'lanthi'
        self.set_split()
#        self.set_monoisotopic_mass()
        self.csv_columns = [self.leader, self.core, self.start, self.end]
        self.CUTOFF = CUTOFF
        
    def set_split(self):
        """Try to identify cleavage site using regular expressions"""
        #Regular expressions; try 1 first, then 2, etc.
        rex1 = re.compile('F?LD')
        rex2 = re.compile('[LF]?LQ')
        end = 0
        #For regular expression, check if there is a match that is >10 AA from the end
        if re.search(rex1, self.sequence) and len(re.split(rex1, self.sequence)[-1]) > 10:
            start, end = [m.span() for m in rex1.finditer(self.sequence)][-1]
            end += 16 
        elif re.search(rex2, self.sequence) and len(re.split(rex2,self.sequence)[-1]) > 10:
            start, end = [m.span() for m in rex2.finditer(self.sequence)][-1]
            end += 15
        self.split_index = end
        self.leader = self.sequence[:end]
        self.core = self.sequence[end:]
        if len(self.leader) < 5 or len(self.core) < 5: #TAG CHRIS
            self.valid_split = False
            self.split = int(.25*len(self.sequence))
            self.leader = self.sequence[:self.split]
            self.core = self.sequence[self.split:]
        
    
    def get_fimo_score(self):
        fimo_output = self.run_fimo_simple()
        fimo_motifs = [int(line.partition("\t")[0]) for line in fimo_output.split("\n") if "\t" in line and line.partition("\t")[0].isdigit()]
        fimo_scores = {int(line.split("\t")[0]): float(line.split("\t")[6]) for line in fimo_output.split("\n") if "\t" in line and line.partition("\t")[0].isdigit()}
        return fimo_motifs, fimo_scores
    
    def set_score(self, pfam_hmm, cust_hmm):
        scoring_csv_columns = []
        score = 0
        #Location of F?LD motif
        if re.search('F[ARNDBCEQZGHILKMFPSTWYV]LD', self.leader) != None:
            score += 2
            scoring_csv_columns.append(re.search('F[ARNDBCEQZGHILKMFPSTWYV]LD', self.leader).span()[0])
        else:
            scoring_csv_columns.append(0)
        #Core residue position of Sx4C motif
        if re.search('S[ARNDBCEQZGHILKMFPSTWYV]{4}C', self.core) != None:
            scoring_csv_columns.append(re.search('S[ARNDBCEQZGHILKMFPSTWYV]{4}C', self.core).span()[0])
        else:
            scoring_csv_columns.append(0)
        #Core residue position of Tx4C motif
        if re.search('T[ARNDBCEQZGHILKMFPSTWYV]{4}C', self.core) != None:
            scoring_csv_columns.append(re.search('T[ARNDBCEQZGHILKMFPSTWYV]{4}C', self.core).span()[0])
        else:
            scoring_csv_columns.append(0)
        #Core residue position of Sx5C motif
        if re.search('S[ARNDBCEQZGHILKMFPSTWYV]{5}C', self.core) != None:
            scoring_csv_columns.append(re.search('S[ARNDBCEQZGHILKMFPSTWYV]{5}C', self.core).span()[0])
        else:
            scoring_csv_columns.append(0)
        #Core residue position of Tx5C motif
        if re.search('T[ARNDBCEQZGHILKMFPSTWYV]{5}C', self.core) != None:
            scoring_csv_columns.append(re.search('T[ARNDBCEQZGHILKMFPSTWYV]{5}C', self.core).span()[0])
        else:
            scoring_csv_columns.append(0)
        
        LANC_like = "PF05147"
        Lant_dehyd_C = "PF04738"
        pfams = []
        bsp_coords = []
        for pfam in self.pfam_2_coords.keys():
            if any(fam in pfam for fam in [LANC_like, Lant_dehyd_C]):
                bsp_coords += self.pfam_2_coords[pfam]
        min_distance = self.get_min_dist(bsp_coords)
        if min_distance < 500:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        
        for pfam_dot in self.pfam_2_coords.keys():
            pfams.append(pfam_dot.split('.')[0])
            
        #Cluster contains LanB dehydratase domain (Lant_dehyd_C)
        if "PF04738" in pfams:
            score += 2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Cluster contains Lan C cyclase domain (LANC_like)
        if "PF05147" in pfams:
            score += 2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Cluster LACKS LanB dehydratase domain (Lant_dehyd_C) 
        if "PF04738" not in pfams:
            score -= 2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Cluster LACKS Lan C cyclase domain (LANC_like)
        if "PF05147" not in pfams:
            score -= 2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Cluster contains LanB dehydratase elimination C-terminal domain (PF14028)
        if "PF14028" in pfams:
            score += 2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Cluster contains S8 peptidase subtilase (Peptidase_S8)
        if "PF00082" in pfams:
            score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Cluster contains C39 peptidase (Peptidase_C39)
        if "PF03412" in pfams:
            score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Cluster contains ABC transporter (PF00005)
        if "PF00005" in pfams:
            score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Cluster contains YcaO-like protein (YcaO)
        if "PF02624" in pfams:
            score -= 4
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Cluster contains ThiF-like protein (ThiF)
        if "PF00899" in pfams:
            score -= 4
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Cluster contains PF02052 (Gallidermin)
        if "PF02052" in pfams:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Cluster contains Antimicr18
        if "PF8130" in pfams:
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Precursor peptide mass < 4000 Da
        precursor_analysis = ProteinAnalysis(self.sequence.replace('X', ''), monoisotopic=True)
        if precursor_analysis.molecular_weight() < 4000:
            score -= 3
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Core peptide mass < 2000 Da
        core_analysis = ProteinAnalysis(self.core.replace('X', ''), monoisotopic=True)
        if core_analysis.molecular_weight() < 2000:
            score -= 3
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        
        precursor_hmm_info = hmmer_utils.get_hmmer_info(self.sequence, pfam_hmm, cust_hmm)
        
        pfams = []
        for pfam_dot, _, _, _, in precursor_hmm_info:
            pfams.append(pfam_dot.split('.')[0])
            
        precursor_hit = False
        #Precursor peptide hits gallidermin superfamily (cl03420) HMM
        if "cl03420" in pfams or "TIGR03731" in pfams or "PF02052" in pfams:
            precursor_hit = True
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Precursor peptide hits lantibio_gallid (TIGR03731) HMM
        if "TIGR03731" in pfams:
            precursor_hit = True
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Precursor peptide hits lanti_SCO0268 superfamily (cl22812) HMM
        if  "cl22812" in pfams:
            precursor_hit = True
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Precursor peptide hits LD_lanti_pre (TIGR04363) HMM
        if "TIGR04363" in pfams:
            precursor_hit = True
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Precursor peptide hits Antimicrobial18 (cl06940) HMM
        if  "cl06940" in pfams:
            precursor_hit = True
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Precursor peptide hits gallidermin (PF02052) HMM
        if  "PF02052" in pfams:
            precursor_hit = True
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #; precursor peptide hits Antimicrobial18 (PF08130) HMM
        if  "PF08130" in pfams:
            precursor_hit = True
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        if precursor_hit:
            score += 3
        
        #Precursor peptide mass (unmodified)
        if "X" in self.sequence:
            Xs = self.sequence.count("X")
            noXprecursor = self.sequence.replace("X", "") #Remove Xs, as they do not work with molecular weight calculation
            precursor_analysis = ProteinAnalysis(noXprecursor, monoisotopic=True)
            scoring_csv_columns.append(float(precursor_analysis.molecular_weight()) + 110 * Xs)
        else:
            precursor_analysis = ProteinAnalysis(self.sequence, monoisotopic=True)
            scoring_csv_columns.append(float(precursor_analysis.molecular_weight()))
        #Unmodified leader peptide mass
        if "X" in self.leader:
            Xs = self.leader.count("X")
            noXleader = self.leader.replace("X", "") #Remove Xs, as they do not work with molecular weight calculation
            leader_analysis = ProteinAnalysis(noXleader, monoisotopic=True)
            scoring_csv_columns.append(float(leader_analysis.molecular_weight()) + 110 * Xs)
        else:
            leader_analysis = ProteinAnalysis(self.leader, monoisotopic=True)
            scoring_csv_columns.append(float(leader_analysis.molecular_weight()))
        #Unmodified core peptide mass
        if "X" in self.core:
            Xs = self.core.count("X")
            noXcore = self.core.replace("X", "") #Remove Xs, as they do not work with molecular weight calculation
            core_analysis = ProteinAnalysis(noXcore, monoisotopic=True)
            scoring_csv_columns.append(float(core_analysis.molecular_weight()) + 110 * Xs)
        else:
            core_analysis = ProteinAnalysis(self.core, monoisotopic=True)
            scoring_csv_columns.append(float(core_analysis.molecular_weight()))
        
        #Length of leader peptide
        scoring_csv_columns.append(len(self.leader))
        #Length of core peptide
        scoring_csv_columns.append(len(self.core))
        #Length of precursor peptide
        scoring_csv_columns.append(len(self.sequence))
        #Ratio of length of leader peptide / length of core peptide
        scoring_csv_columns.append(float(len(self.leader) / float(len(self.core))))
        #Core peptide ≥ 35 residues
        if len(self.core) >= 35:
            score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Core peptide contains CC motif (not in last 3 residues)
        if re.search('CC', self.core[:-3]) != None:
            score -= 3
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Leader peptide has > 4 negatively charge motifs
        if sum([self.leader.count(aa) for aa in "DE"]) > 4:
            score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Leader peptide has net negative charge
        charge_dict = {"E": -1, "D": -1, "K": 1, "R": 1}
        if sum([charge_dict[aa] for aa in self.leader if aa in charge_dict]) < 0:
            score += 1
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Leader residue position of FxLD motif
        if re.search('F[ARNDBCEQZGHILKMFPSTWYV]LD', self.leader) != None:
            scoring_csv_columns.append(re.search('F[ARNDBCEQZGHILKMFPSTWYV]LD', self.leader).span()[0])
        else:
            scoring_csv_columns.append(0)
        #Core peptide contains C-terminal CC (within last 3 residues)
        if re.search('CC', self.core[-3:]) != None:
            score += 2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
            
        #Core peptide contains DGCGxTC / SFNS / SxxLC / CTxGC / TPGC / SFNSxC motifs
        motifs = (('DGCG[ARNDBCEQZGHILKMFPSTWYV]TC', 2), ('SFNS', 2), \
           ('S[ARNDBCEQZGHILKMFPSTWYV]{2}LC', 2),('CT[ARNDBCEQZGHILKMFPSTWYV]{1}GC', 1), \
           ('TPGC', 1), ('SFNS[ARNDBCEQZGHILKMFPSTWYV]C', 1))
        for motif in motifs:
            if re.search(motif[0], self.core) != None:
                score += motif[1]
                scoring_csv_columns.append(1)
            else:
                scoring_csv_columns.append(0)
        #Core peptide contains < 2 or < 3 Cys
        if self.core.count("C") < 2:
            score -= 6
            scoring_csv_columns += [1,1]
        elif self.core.count("C") < 3:
            score -= 3
            scoring_csv_columns += [1,0]
        else:
            scoring_csv_columns += [0,0]
        #No Cys/Ser/Thr in core peptide
        for aa, penalty in [("C", -10), ("S", -4), ("T", -4)]:
            if aa not in self.core:
                score += penalty
                scoring_csv_columns.append(1)
            else:
                scoring_csv_columns.append(0)
                
        numringlist, strings, profile = self.lanscout([[self.core]])
        if numringlist[0] > 4:
            score += 2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Lanthionine regex maximum ring number < 3
        if numringlist[0] < 3:
            score -= 2
            scoring_csv_columns.append(1)
        else:
            scoring_csv_columns.append(0)
        #Lanthionine regex 4-membered ring/5-membered ring/6-membered ring/7-membered ring/8-membered ring
        scores = [2, 2, 2, 2, 1]
        scorepos = 0
        for ringsize in profile[0].split(", ")[:2]:
            if ringsize != "0" and ringsize != "1" and ringsize != "2":
                score += scores[scorepos]
                scoring_csv_columns.append(1)
            else:
                scoring_csv_columns.append(0)
            scorepos += 1
        for ringsize in profile[0].split(", ")[2:]:
            if ringsize != "0":
                score += scores[scorepos]
                scoring_csv_columns.append(1)
            else:
                scoring_csv_columns.append(0)
            scorepos += 1
        self.scoring_csv_columns = scoring_csv_columns
        self.score = score
        
        self.svm_columns = []
        self.svm_columns += self.scoring_csv_columns
        
        fimo_motifs, fimo_scores  = self.get_fimo_score()
        if len(fimo_motifs) > 0:
            self.svm_columns.append(1)
        else:
            self.svm_columns.append(0)
        self.svm_columns.append(numringlist[0])
        #Lanthionine regex 4-membered ring count
        self.svm_columns.append(int(profile[0].split(", ")[0]))
        #Lanthionine regex 5-membered ring count
        self.svm_columns.append(int(profile[0].split(", ")[1]))
        #Lanthionine regex 6-membered ring count
        self.svm_columns.append(int(profile[0].split(", ")[2]))
        #Lanthionine regex 7-membered ring count
        self.svm_columns.append(int(profile[0].split(", ")[3]))
        #Lanthionine regex 8-membered ring count
        self.svm_columns.append(int(profile[0].split(", ")[4]))
        #Ratio of number of Cys in core peptide to sum of Ser/Thr in core peptide
        if "S" in self.core or "T" in self.core:
            self.svm_columns.append(self.core.count("C") / float(self.core.count("S") + self.core.count("T")))
        else:
            self.svm_columns.append(1.0)
        #Ratio of number of Cys/Ser/Thr to length of core peptide
        self.svm_columns.append(float(self.core.count("S") + self.core.count("T") + self.core.count("C")) / len(self.core))
        
#        ##log10 of scores
        self.svm_columns.append(sum(fimo_scores))
        
        if 1 in fimo_motifs:
            self.svm_columns.append(fimo_scores[1])
        else:
            self.svm_columns.append(0)
        #log10 p-value MEME motif 2
        if 2 in fimo_motifs:
            self.svm_columns.append(fimo_scores[2])
        else:
            self.svm_columns.append(0)
        #log10 p-value MEME motif 3
        if 3 in fimo_motifs:
            self.svm_columns.append(fimo_scores[3])
        else:
            self.svm_columns.append(0)
        #log10 p-value MEME motif 4
        if 4 in fimo_motifs:
            self.svm_columns.append(fimo_scores[4])
        else:
            self.svm_columns.append(0)
        #log10 p-value MEME motif 5
        if 5 in fimo_motifs:
            self.svm_columns.append(fimo_scores[5])
        else:
            self.svm_columns.append(0)
        #Number in leader of each amino acid
        self.svm_columns += [self.leader.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
        #Number in leader of each amino acid type (aromatic, aliphatic, hydroxyl, basic, acidic)
        self.svm_columns.append(sum([self.leader.count(aa) for aa in "FWY"]))
        self.svm_columns.append(sum([self.leader.count(aa) for aa in "DE"]))
        self.svm_columns.append(sum([self.leader.count(aa) for aa in "RK"]))
        self.svm_columns.append(sum([self.leader.count(aa) for aa in "RKDE"]))
        self.svm_columns.append(sum([self.leader.count(aa) for aa in "GAVLMI"]))
        self.svm_columns.append(sum([self.leader.count(aa) for aa in "ST"]))
        #Number in core of each amino acid
        self.svm_columns += [self.core.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
        #Number in core of each amino acid type (aromatic, aliphatic, hydroxyl, basic, acidic)
        self.svm_columns.append(sum([self.core.count(aa) for aa in "FWY"]))
        self.svm_columns.append(sum([self.core.count(aa) for aa in "DE"]))
        self.svm_columns.append(sum([self.core.count(aa) for aa in "RK"]))
        self.svm_columns.append(sum([self.core.count(aa) for aa in "RKDE"]))
        self.svm_columns.append(sum([self.core.count(aa) for aa in "GAVLMI"]))
        self.svm_columns.append(sum([self.core.count(aa) for aa in "ST"]))
        #Number in entire precursor of each amino acid
        self.svm_columns += [self.sequence.count(aa) for aa in "ARDNCQEGHILKMFPSTWYV"]
        #Number in entire precursor of each amino acid type (aromatic, aliphatic, hydroxyl, basic, acidic)
        self.svm_columns.append(sum([self.sequence.count(aa) for aa in "FWY"]))
        self.svm_columns.append(sum([self.sequence.count(aa) for aa in "DE"]))
        self.svm_columns.append(sum([self.sequence.count(aa) for aa in "RK"]))
        self.svm_columns.append(sum([self.sequence.count(aa) for aa in "RKDE"]))
        self.svm_columns.append(sum([self.sequence.count(aa) for aa in "GAVLMI"]))
        self.svm_columns.append(sum([self.sequence.count(aa) for aa in "ST"]))
        self.svm_columns = [self.score] + self.svm_columns
        self.csv_columns += self.svm_columns
        
                
    def lanscout(self, seq):
        #nis = 'ITSISLCTPGCKTGALMGCNMKTATCHCSIHVSK'
        #define lanthionine ring with a c-terminal Cys
        lanlower = 2
        lanupper = 6
        lan = '(?=([T|S].{%d,%d}C))' % (lanlower, lanupper)
        db = []
        strings = []
        numringlist = []
        seq = [num for elem in seq for num in elem]
        for n in range(0, len(seq)):
            core = str(seq[n][:])
            strings.append(core)
            totrings = re.compile(lan, re.I).findall(core)
            size = []
            for i in range(0, len(totrings)):
                size.append(len(totrings[i]))
            db.append(size)
            numrings = len(totrings)
            numringlist.append(numrings)
    
        profile = []
        for i in range(0,len(db)):
            #print db[i]
            temp = []
            for j in range(lanlower+2,lanupper+3):
                temp.append(db[i].count(j))
            profile.append(temp)
    
        for i in range(0,len(profile)):
            #profile[i] = str(profile[i])
            profile[i]=str(profile[i]).strip('[]')
        return numringlist, strings, profile
