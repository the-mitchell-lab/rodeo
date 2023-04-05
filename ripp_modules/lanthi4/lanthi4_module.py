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
import numpy as np
from ripp_modules.VirtualRipp import VirtualRipp
import hmmer_utils
import pathlib
FILE_DIR = pathlib.Path(__file__).parent.absolute()

peptide_type = "lanthi4"
CUTOFF = 20
index = 0

def write_csv_headers(output_dir):
    dir_prefix = output_dir + '/{}/'.format(peptide_type)
    if not os.path.exists(dir_prefix):
        os.makedirs(dir_prefix)
    svm_headers = 'PK,Classification,1/5_frac_c, 1/5_frac_s_t, 1/5_frac_s_t_c, 2/5_frac_c, 2/5_frac_s_t, 2/5_frac_s_t_c, 3/5_frac_c, 3/5_frac_s_t, 3/5_frac_s_t_c, 4/5_frac_c, 4/5_frac_s_t, 4/5_frac_s_t_c, 5/5_frac_c, 5/5_frac_s_t, 5/5_frac_s_t_c, S+T, aliphatic, aromatic, charged, core_A, core_AA, core_AD, core_AE, core_AG, core_AH, core_AI, core_AK, core_AM, core_AP, core_AQ, core_AR, core_AS, core_AT, core_AV, core_AW, core_AxA, core_AxG, core_AxH, core_AxK, core_AxL, core_AxM, core_AxP, core_AxQ, core_AxR, core_AxS, core_AxV, core_AxW, core_AxxA, core_AxxD, core_AxxE, core_AxxG, core_AxxH, core_AxxI, core_AxxL, core_AxxM, core_AxxP, core_AxxR, core_AxxS, core_AxxV, core_AxxxA, core_AxxxE, core_AxxxG, core_AxxxH, core_AxxxM, core_AxxxP, core_AxxxQ, core_AxxxR, core_AxxxS, core_AxxxT, core_AxxxW, core_AxxxxA, core_AxxxxD, core_AxxxxE, core_AxxxxF, core_AxxxxG, core_AxxxxH, core_AxxxxK, core_AxxxxM, core_AxxxxP, core_AxxxxQ, core_AxxxxR, core_AxxxxS, core_AxxxxW, core_AxxxxxA, core_AxxxxxC, core_AxxxxxD, core_AxxxxxG, core_AxxxxxH, core_AxxxxxL, core_AxxxxxP, core_AxxxxxQ, core_AxxxxxR, core_AxxxxxS, core_AxxxxxT, core_C, core_CA, core_CD, core_CE, core_CG, core_CI, core_CK, core_CL, core_CP, core_CR, core_CS, core_CT, core_CV, core_CxA, core_CxC, core_CxF, core_CxG, core_CxI, core_CxP, core_CxR, core_CxW, core_CxY, core_CxxA, core_CxxG, core_CxxI, core_CxxL, core_CxxN, core_CxxP, core_CxxR, core_CxxS, core_CxxT, core_CxxV, core_CxxxA, core_CxxxC, core_CxxxG, core_CxxxL, core_CxxxP, core_CxxxR, core_CxxxS, core_CxxxT, core_CxxxV, core_CxxxxA, core_CxxxxC, core_CxxxxD, core_CxxxxG, core_CxxxxI, core_CxxxxL, core_CxxxxP, core_CxxxxR, core_CxxxxT, core_CxxxxW, core_CxxxxxA, core_CxxxxxD, core_CxxxxxE, core_CxxxxxG, core_CxxxxxI, core_CxxxxxL, core_CxxxxxP, core_CxxxxxR, core_CxxxxxW, core_D, core_DA, core_DE, core_DG, core_DR, core_DS, core_DV, core_DxG, core_DxL, core_DxN, core_DxP, core_DxR, core_DxS, core_DxW, core_DxxG, core_DxxL, core_DxxR, core_DxxxA, core_DxxxP, core_DxxxR, core_DxxxS, core_DxxxT, core_DxxxV, core_DxxxxG, core_DxxxxP, core_DxxxxR, core_DxxxxxA, core_DxxxxxC, core_DxxxxxR, core_DxxxxxS, core_DxxxxxV, core_EC, core_EG, core_EP, core_EQ, core_ER, core_ES, core_ET, core_ExA, core_ExL, core_ExQ, core_ExR, core_ExV, core_ExxG, core_ExxL, core_ExxR, core_ExxS, core_ExxV, core_ExxxA, core_ExxxC, core_ExxxG, core_ExxxL, core_ExxxP, core_ExxxR, core_ExxxxC, core_ExxxxL, core_ExxxxR, core_ExxxxS, core_ExxxxT, core_ExxxxxC, core_ExxxxxL, core_ExxxxxP, core_ExxxxxR, core_ExxxxxT, core_FA, core_FC, core_FG, core_FL, core_FP, core_FR, core_FT, core_FV, core_FxA, core_FxC, core_FxL, core_FxP, core_FxR, core_FxxA, core_FxxD, core_FxxL, core_FxxP, core_FxxR, core_FxxT, core_FxxV, core_FxxxA, core_FxxxC, core_FxxxI, core_FxxxP, core_FxxxR, core_FxxxS, core_FxxxxE, core_FxxxxG, core_FxxxxN, core_FxxxxP, core_FxxxxR, core_FxxxxS, core_FxxxxT, core_FxxxxxC, core_FxxxxxG, core_FxxxxxL, core_FxxxxxP, core_FxxxxxR, core_G, core_GA, core_GD, core_GE, core_GG, core_GI, core_GL, core_GM, core_GP, core_GR, core_GS, core_GV, core_GW, core_GxA, core_GxE, core_GxF, core_GxG, core_GxH, core_GxI, core_GxL, core_GxP, core_GxR, core_GxV, core_GxW, core_GxxA, core_GxxC, core_GxxD, core_GxxG, core_GxxH, core_GxxK, core_GxxL, core_GxxP, core_GxxR, core_GxxS, core_GxxT, core_GxxV, core_GxxW, core_GxxxD, core_GxxxE, core_GxxxG, core_GxxxH, core_GxxxL, core_GxxxP, core_GxxxR, core_GxxxS, core_GxxxV, core_GxxxW, core_GxxxY, core_GxxxxA, core_GxxxxC, core_GxxxxD, core_GxxxxE, core_GxxxxG, core_GxxxxH, core_GxxxxL, core_GxxxxP, core_GxxxxQ, core_GxxxxR, core_GxxxxS, core_GxxxxT, core_GxxxxV, core_GxxxxW, core_GxxxxxA, core_GxxxxxC, core_GxxxxxD, core_GxxxxxE, core_GxxxxxF, core_GxxxxxG, core_GxxxxxH, core_GxxxxxN, core_GxxxxxP, core_GxxxxxR, core_GxxxxxS, core_GxxxxxV, core_GxxxxxW, core_H, core_HA, core_HH, core_HL, core_HP, core_HR, core_HS, core_HV, core_HxG, core_HxP, core_HxR, core_HxS, core_HxT, core_HxV, core_HxxA, core_HxxD, core_HxxL, core_HxxP, core_HxxR, core_HxxS, core_HxxV, core_HxxxA, core_HxxxG, core_HxxxP, core_HxxxR, core_HxxxS, core_HxxxV, core_HxxxxA, core_HxxxxG, core_HxxxxP, core_HxxxxR, core_HxxxxS, core_HxxxxV, core_HxxxxxA, core_HxxxxxG, core_HxxxxxP, core_HxxxxxR, core_HxxxxxS, core_HxxxxxV, core_HxxxxxY, core_I, core_IE, core_IG, core_IP, core_IR, core_IS, core_IT, core_IxA, core_IxC, core_IxG, core_IxI, core_IxP, core_IxQ, core_IxR, core_IxT, core_IxV, core_IxxA, core_IxxC, core_IxxE, core_IxxI, core_IxxP, core_IxxQ, core_IxxS, core_IxxxC, core_IxxxE, core_IxxxP, core_IxxxxA, core_IxxxxQ, core_IxxxxR, core_IxxxxxP, core_IxxxxxQ, core_IxxxxxR, core_IxxxxxS, core_K, core_KG, core_KP, core_KR, core_KS, core_KY, core_KxA, core_KxC, core_KxG, core_KxL, core_KxR, core_KxxA, core_KxxR, core_KxxS, core_KxxxF, core_KxxxP, core_KxxxR, core_KxxxxA, core_KxxxxP, core_KxxxxR, core_KxxxxS, core_KxxxxT, core_KxxxxxA, core_KxxxxxC, core_KxxxxxR, core_KxxxxxS, core_LF, core_LG, core_LH, core_LI, core_LK, core_LL, core_LN, core_LR, core_LS, core_LT, core_LV, core_LxA, core_LxH, core_LxK, core_LxL, core_LxN, core_LxP, core_LxR, core_LxS, core_LxT, core_LxxA, core_LxxC, core_LxxF, core_LxxK, core_LxxL, core_LxxP, core_LxxR, core_LxxS, core_LxxT, core_LxxV, core_LxxY, core_LxxxA, core_LxxxC, core_LxxxE, core_LxxxF, core_LxxxG, core_LxxxL, core_LxxxN, core_LxxxP, core_LxxxR, core_LxxxS, core_LxxxV, core_LxxxxA, core_LxxxxC, core_LxxxxG, core_LxxxxH, core_LxxxxI, core_LxxxxN, core_LxxxxP, core_LxxxxR, core_LxxxxS, core_LxxxxT, core_LxxxxxA, core_LxxxxxD, core_LxxxxxF, core_LxxxxxG, core_LxxxxxI, core_LxxxxxL, core_LxxxxxP, core_LxxxxxR, core_LxxxxxS, core_LxxxxxT, core_LxxxxxV, core_M, core_MP, core_MR, core_MS, core_MV, core_MW, core_MxR, core_MxxA, core_MxxP, core_MxxR, core_MxxS, core_MxxxR, core_MxxxS, core_MxxxxA, core_MxxxxR, core_MxxxxxR, core_NR, core_NT, core_NxA, core_NxP, core_NxR, core_NxV, core_NxxC, core_NxxP, core_NxxR, core_NxxxI, core_NxxxR, core_NxxxS, core_NxxxV, core_NxxxxG, core_NxxxxP, core_NxxxxR, core_NxxxxxC, core_NxxxxxI, core_NxxxxxP, core_NxxxxxR, core_NxxxxxS, core_NxxxxxV, core_P, core_PA, core_PC, core_PD, core_PE, core_PF, core_PG, core_PH, core_PK, core_PL, core_PP, core_PQ, core_PR, core_PS, core_PT, core_PV, core_PxA, core_PxC, core_PxD, core_PxE, core_PxG, core_PxH, core_PxK, core_PxL, core_PxN, core_PxP, core_PxQ, core_PxR, core_PxS, core_PxT, core_PxV, core_PxW, core_PxY, core_PxxA, core_PxxD, core_PxxG, core_PxxH, core_PxxI, core_PxxL, core_PxxP, core_PxxQ, core_PxxR, core_PxxS, core_PxxV, core_PxxW, core_PxxxA, core_PxxxC, core_PxxxD, core_PxxxE, core_PxxxF, core_PxxxG, core_PxxxH, core_PxxxL, core_PxxxP, core_PxxxQ, core_PxxxR, core_PxxxS, core_PxxxV, core_PxxxxA, core_PxxxxD, core_PxxxxF, core_PxxxxG, core_PxxxxH, core_PxxxxI, core_PxxxxL, core_PxxxxP, core_PxxxxQ, core_PxxxxR, core_PxxxxS, core_PxxxxT, core_PxxxxV, core_PxxxxxA, core_PxxxxxC, core_PxxxxxG, core_PxxxxxH, core_PxxxxxL, core_PxxxxxP, core_PxxxxxQ, core_PxxxxxR, core_PxxxxxS, core_PxxxxxV, core_PxxxxxW, core_Q, core_QA, core_QG, core_QP, core_QQ, core_QR, core_QV, core_QY, core_QxA, core_QxF, core_QxG, core_QxP, core_QxR, core_QxT, core_QxV, core_QxxC, core_QxxL, core_QxxR, core_QxxS, core_QxxxA, core_QxxxG, core_QxxxP, core_QxxxR, core_QxxxS, core_QxxxxA, core_QxxxxE, core_QxxxxI, core_QxxxxP, core_QxxxxR, core_QxxxxS, core_QxxxxxA, core_QxxxxxG, core_QxxxxxP, core_QxxxxxR, core_R, core_RA, core_RC, core_RD, core_RE, core_RG, core_RH, core_RI, core_RK, core_RL, core_RM, core_RN, core_RP, core_RQ, core_RR, core_RS, core_RT, core_RV, core_RW, core_RxA, core_RxC, core_RxD, core_RxE, core_RxF, core_RxG, core_RxH, core_RxI, core_RxK, core_RxL, core_RxM, core_RxN, core_RxP, core_RxQ, core_RxR, core_RxS, core_RxT, core_RxV, core_RxW, core_RxxA, core_RxxD, core_RxxE, core_RxxF, core_RxxG, core_RxxH, core_RxxI, core_RxxK, core_RxxL, core_RxxM, core_RxxP, core_RxxQ, core_RxxR, core_RxxS, core_RxxT, core_RxxV, core_RxxW, core_RxxxA, core_RxxxC, core_RxxxE, core_RxxxF, core_RxxxG, core_RxxxH, core_RxxxI, core_RxxxL, core_RxxxM, core_RxxxN, core_RxxxP, core_RxxxQ, core_RxxxR, core_RxxxS, core_RxxxT, core_RxxxV, core_RxxxY, core_RxxxxA, core_RxxxxC, core_RxxxxD, core_RxxxxE, core_RxxxxG, core_RxxxxH, core_RxxxxI, core_RxxxxL, core_RxxxxM, core_RxxxxN, core_RxxxxP, core_RxxxxQ, core_RxxxxR, core_RxxxxS, core_RxxxxT, core_RxxxxV, core_RxxxxW, core_RxxxxxA, core_RxxxxxC, core_RxxxxxD, core_RxxxxxE, core_RxxxxxF, core_RxxxxxG, core_RxxxxxH, core_RxxxxxI, core_RxxxxxK, core_RxxxxxL, core_RxxxxxM, core_RxxxxxN, core_RxxxxxP, core_RxxxxxQ, core_RxxxxxR, core_RxxxxxS, core_RxxxxxT, core_RxxxxxV, core_RxxxxxW, core_S, core_SA, core_SC, core_SF, core_SG, core_SH, core_SI, core_SK, core_SL, core_SM, core_SN, core_SP, core_SR, core_SS, core_ST, core_SV, core_SxD, core_SxE, core_SxG, core_SxH, core_SxI, core_SxK, core_SxL, core_SxP, core_SxR, core_SxS, core_SxT, core_SxV, core_SxW, core_SxY, core_SxxA, core_SxxD, core_SxxE, core_SxxF, core_SxxG, core_SxxH, core_SxxP, core_SxxR, core_SxxS, core_SxxT, core_SxxV, core_SxxxA, core_SxxxC, core_SxxxD, core_SxxxF, core_SxxxG, core_SxxxM, core_SxxxP, core_SxxxR, core_SxxxS, core_SxxxV, core_SxxxW, core_SxxxxA, core_SxxxxF, core_SxxxxG, core_SxxxxH, core_SxxxxI, core_SxxxxL, core_SxxxxP, core_SxxxxR, core_SxxxxS, core_SxxxxT, core_SxxxxW, core_SxxxxxA, core_SxxxxxD, core_SxxxxxF, core_SxxxxxG, core_SxxxxxH, core_SxxxxxI, core_SxxxxxK, core_SxxxxxL, core_SxxxxxM, core_SxxxxxP, core_SxxxxxR, core_SxxxxxS, core_T, core_TC, core_TH, core_TI, core_TL, core_TN, core_TP, core_TV, core_TW, core_TxC, core_TxE, core_TxI, core_TxL, core_TxP, core_TxT, core_TxxA, core_TxxC, core_TxxF, core_TxxI, core_TxxN, core_TxxR, core_TxxV, core_TxxW, core_TxxxC, core_TxxxE, core_TxxxF, core_TxxxG, core_TxxxK, core_TxxxP, core_TxxxR, core_TxxxT, core_TxxxV, core_TxxxxA, core_TxxxxC, core_TxxxxD, core_TxxxxF, core_TxxxxG, core_TxxxxI, core_TxxxxP, core_TxxxxQ, core_TxxxxR, core_TxxxxT, core_TxxxxV, core_TxxxxW, core_TxxxxY, core_TxxxxxC, core_TxxxxxD, core_TxxxxxG, core_TxxxxxI, core_TxxxxxL, core_TxxxxxP, core_TxxxxxR, core_TxxxxxS, core_TxxxxxT, core_V, core_VC, core_VD, core_VF, core_VG, core_VH, core_VI, core_VK, core_VL, core_VP, core_VQ, core_VR, core_VT, core_VxA, core_VxC, core_VxD, core_VxH, core_VxI, core_VxL, core_VxP, core_VxR, core_VxS, core_VxT, core_VxxC, core_VxxF, core_VxxH, core_VxxI, core_VxxL, core_VxxP, core_VxxR, core_VxxT, core_VxxV, core_VxxW, core_VxxxA, core_VxxxC, core_VxxxG, core_VxxxH, core_VxxxI, core_VxxxR, core_VxxxS, core_VxxxT, core_VxxxxC, core_VxxxxD, core_VxxxxG, core_VxxxxH, core_VxxxxK, core_VxxxxL, core_VxxxxP, core_VxxxxQ, core_VxxxxR, core_VxxxxS, core_VxxxxV, core_VxxxxW, core_VxxxxxA, core_VxxxxxC, core_VxxxxxF, core_VxxxxxG, core_VxxxxxH, core_VxxxxxL, core_VxxxxxP, core_VxxxxxQ, core_VxxxxxR, core_VxxxxxV, core_VxxxxxY, core_W, core_WA, core_WC, core_WP, core_WR, core_WS, core_WT, core_WxA, core_WxC, core_WxG, core_WxP, core_WxR, core_WxV, core_WxxR, core_WxxS, core_WxxV, core_WxxxA, core_WxxxC, core_WxxxG, core_WxxxI, core_WxxxK, core_WxxxP, core_WxxxQ, core_WxxxR, core_WxxxV, core_WxxxxG, core_WxxxxP, core_WxxxxR, core_WxxxxS, core_WxxxxT, core_WxxxxV, core_WxxxxxG, core_WxxxxxI, core_WxxxxxP, core_YC, core_YR, core_YT, core_YxC, core_YxR, core_YxS, core_YxxF, core_YxxxT, core_YxxxxC, core_YxxxxxC, core_YxxxxxL, core_YxxxxxR, core_length, core_pI, frac_S+T+C, frac_aliphatic, frac_aromatic, frac_c, frac_charge, frac_polar, frac_pos, neg, net_charge, peptide_length, pos'

    svm_headers = svm_headers.split(',')
    features_headers = ["Accession_id", "Genus/Species/Code", "Leader", "Core", "Start", "End", "Total Score", "Valid Precursor" ] + svm_headers
    features_csv_file = open(dir_prefix + "temp_features.csv", 'w')
    svm_csv_file = open("{}/svm/fitting_set.csv".format(FILE_DIR), 'w')
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
                 pfam_2_coords):
        super(Ripp, self).__init__(start, 
                                     end, 
                                     sequence,
                                     upstream_sequence,
                                     pfam_2_coords)
        self.peptide_type = peptide_type
        self.score = 0
        self.set_split()
#        self.set_monoisotopic_mass()
        self.csv_columns = [self.leader, self.core, self.start, self.end]
        self.CUTOFF = CUTOFF
    
    def get_fimo_score(self, fimo_file):
        fimo_output = self.run_fimo_simple(fimo_file)
        fimo_motifs = [int(line.partition("\t")[0]) for line in fimo_output.split("\n") if "\t" in line and line.partition("\t")[0].isdigit()]
        fimo_stops = []
        for line in fimo_output.split("\n"):
            line = line.split("\t")
            if len(line) > 3 and line[3].isdigit():
                fimo_stops.append(int(line[3]))
        fimo_scores = {int(line.split("\t")[0]): float(line.split("\t")[6]) for line in fimo_output.split("\n") if "\t" in line and line.partition("\t")[0].isdigit()}
        return fimo_motifs, fimo_scores, fimo_stops
    
    def set_split(self):
        """Try to identify cleavage site using regular expressions"""
        # switch the # on the two lines below if your computer username has a space in it. # RAL
        #f_scores = self.get_fimo_score("{}_leader_fimo.txt".format("./ripp_modules/lanthi4/lanthi4"))
        f_scores = self.get_fimo_score("{}/{}_leader_fimo.txt".format(FILE_DIR, peptide_type))  
        l_stop = sorted([0] + f_scores[2])[-1]
        if len(f_scores[0]) > 0:
            self.score += 10
        regex1 = re.search('(.[S|T].{2,7}C)', self.sequence[l_stop:])
        regex2 = re.search('(G[G|A|S])', self.sequence[l_stop:])
    
        if regex1 and regex2:
            if regex2.end() < regex1.start() and regex2.end()+l_stop > 10:
                self.split_index = regex2.end()+l_stop
            elif regex1.start() + l_stop > 10:
                self.split_index = regex1.start()+l_stop
            else:
                self.split_index =  (int)((len(self.sequence)-l_stop)/2)+l_stop
        elif regex1 and not regex2:
            if regex1.start() + l_stop > 10:
                self.split_index = regex1.start()+l_stop
            else:
                self.split_index = (int)((len(self.sequence)-l_stop)/2)+l_stop
        elif regex2 and not regex1:
            if regex2.end() + l_stop > 10:
                self.split_index = regex2.end()+l_stop
            else:
                self.split_index = (int)((len(self.sequence)-l_stop)/2)+l_stop
        else:
            self.split_index = (int)((len(self.sequence)-l_stop)/2)+l_stop
        self.leader = self.sequence[:self.split_index]
        self.core = self.sequence[self.split_index:]
        if len(self.core) < 10:
            self.split_index = len(self.sequence) // 2
            self.leader = self.sequence[:self.split_index]
            self.core = self.sequence[self.split_index:]
    
    def set_score(self, pfam_hmm, cust_hmm):
        scoring_csv_columns = []
        scoring_column_names = []
        c_count = self.core.count("C")
        s_count = self.core.count("S")
        t_count = self.core.count("T")
        fifth_size = len(self.sequence) // 5
        for i in range(5):
            scoring_csv_columns.append(self.sequence[i*fifth_size:(i+1)*fifth_size].count('C') / len(self.sequence))
            scoring_column_names.append("{}/5_frac_c".format(i+1))
                    
        for i in range(5):
            count = 0
            for c in "ST":
                count += self.sequence[i*fifth_size:(i+1)*fifth_size].count(c)
            scoring_csv_columns.append(count / len(self.sequence))
            scoring_column_names.append("{}/5_frac_s_t".format(i+1))
            
        for i in range(5):
            count = 0
            for c in "STC":
                count += self.sequence[i*fifth_size:(i+1)*fifth_size].count(c)
            scoring_csv_columns.append(count / len(self.sequence))
            scoring_column_names.append("{}/5_frac_s_t_c".format(i+1))

           
        scoring_csv_columns.append((c_count) / len(self.sequence))
        scoring_column_names.append("frac_c")
        
        scoring_csv_columns.append((s_count + t_count))
        scoring_column_names.append("S+T")
        
        scoring_csv_columns.append((s_count + t_count + c_count) / len(self.core))
        scoring_column_names.append("frac_S+T+C")
#        
        aliphatic_core_count = self.core.count("A") + self.core.count("I") + self.core.count("M") + self.core.count("L") + self.core.count("V")
        scoring_csv_columns.append(aliphatic_core_count)
        scoring_column_names.append("aliphatic")
        scoring_csv_columns.append(aliphatic_core_count / len(self.core))
        scoring_column_names.append("frac_aliphatic")
        
        aromatic_core_count = self.core.count("F") + self.core.count("W") + self.core.count("Y")
        scoring_csv_columns.append(aromatic_core_count)
        scoring_column_names.append("aromatic")
        scoring_csv_columns.append(aromatic_core_count / len(self.core))
        scoring_column_names.append("frac_aromatic")
        
        charged_count = self.core.count("D") + self.core.count("E") + self.core.count("K") + self.core.count("R")
        scoring_csv_columns.append(charged_count)
        scoring_column_names.append("charged")
        scoring_csv_columns.append(charged_count  / len(self.core))
        scoring_column_names.append("frac_charge")
        
        polar_count = c_count + s_count + t_count + self.core.count("H") + self.core.count("N") + self.core.count("Q")
#        scoring_csv_columns.append(polar_count)
#        scoring_column_names.append("polar")
        scoring_csv_columns.append(polar_count  / len(self.core))
        scoring_column_names.append("frac_polar")
        
        pos = self.core.count("R") + self.core.count("K")
        scoring_csv_columns.append(pos)
        scoring_column_names.append("pos")
        scoring_csv_columns.append(pos / len(self.core))
        scoring_column_names.append("frac_pos")
        
        
        neg = self.core.count("D") + self.core.count("E")
        scoring_csv_columns.append(neg)
        scoring_column_names.append("neg")
#        scoring_csv_columns.append(neg / len(self.core))
#        scoring_column_names.append("frac_neg")
        
        scoring_csv_columns.append(neg + pos)
        scoring_column_names.append("net_charge")
        
#        scoring_csv_columns.append(len(self.leader))
#        scoring_column_names.append("leader_length")
        
        scoring_csv_columns.append(len(self.sequence))
        scoring_column_names.append("peptide_length")
        
        scoring_csv_columns.append(len(self.core))
        scoring_column_names.append("core_length")
        
        core_seqs = ['core_IxV', 'core_SxxxxI', 'core_DxxxR', 'core_RxxW', 'core_WxP', 'core_IxQ', 'core_WxxxV', 'core_RxxxxxD', 'core_GxR', 'core_LxxP', 'core_FP', 'core_CxF', 'core_PxxxxxR', 'core_DxxxP', 'core_RD', 'core_HxxxxxS', 'core_LxxS', 'core_HxxxxxR', 'core_GxxxxR', 'core_FxxxxxP', 'core_FxxT', 'core_GxxL', 'core_NxxP', 'core_PxxxxxG', 'core_PxxxxxW', 'core_NxxC', 'core_SxxxxxS', 'core_YxxxxC', 'core_SxK', 'core_LxxxS', 'core_AxxxH', 'core_PxL', 'core_FxC', 'core_IxxxC', 'core_LxxxxxG', 'core_GxxC', 'core_AxxG', 'core_GxxxxxA', 'core_IxC', 'core_GxxxxxC', 'core_TxxxxxS', 'core_SK', 'core_RH', 'core_QA', 'core_LxxxxxR', 'core_LxT', 'core_PxxxxG', 'core_VK', 'core_CxxI', 'core_RxxP', 'core_SH', 'core_LxxxxxF', 'core_CG', 'core_KxC', 'core_GxxxxA', 'core_WxxxxT', 'core_A', 'core_GD', 'core_VxxxxxC', 'core_AxxxxxP', 'core_RxxxxS', 'core_SxY', 'core_LI', 'core_FxxxxxG', 'core_RxxxC', 'core_AxxxxxD', 'core_PxW', 'core_NxxxxG', 'core_GxxxxS', 'core_PG', 'core_AxW', 'core_GxxG', 'core_CS', 'core_GxxxxxN', 'core_WxxV', 'core_RE', 'core_D', 'core_DxL', 'core_YxxxxxC', 'core_NxxxxP', 'core_SxxD', 'core_SxxxxG', 'core_FG', 'core_AE', 'core_AxxxxxG', 'core_VxI', 'core_SxxR', 'core_AxxxxA', 'core_NxxxxxP', 'core_WR', 'core_GxxxxV', 'core_SxH', 'core_LT', 'core_PxxxxP', 'core_SxW', 'core_YxxxT', 'core_SN', 'core_RxxxxM', 'core_KxxxxxR', 'core_GS', 'core_AxxR', 'core_PxxQ', 'core_PxxxxL', 'core_R', 'core_RxxxxV', 'core_CP', 'core_GxA', 'core_IxxI', 'core_FxxV', 'core_LxH', 'core_VxxxxC', 'core_IxP', 'core_CxxxL', 'core_PxxxxS', 'core_SxxxxxF', 'core_NxxxxxS', 'core_TP', 'core_VP', 'core_HxV', 'core_PxH', 'core_VxP', 'core_FxxxxR', 'core_VC', 'core_RN', 'core_QP', 'core_CxxxxC', 'core_GP', 'core_CxxxxR', 'core_SV', 'core_VI', 'core_MxxP', 'core_PS', 'core_SxxxD', 'core_TxxxxV', 'core_VxxxxR', 'core_SxxxxxL', 'core_AxQ', 'core_GxxD', 'core_LxxT', 'core_TxxxG', 'core_KY', 'core_VxxxxH', 'core_DxxxxG', 'core_SxxxS', 'core_TxxA', 'core_TxxxxI', 'core_KxR', 'core_TxxxxxT', 'core_TxxxxxI', 'core_SxxxR', 'core_SF', 'core_NxxxR', 'core_LxxxxxD', 'core_IxxxxxS', 'core_AxxxG', 'core_VxC', 'core_HxxxxxY', 'core_VD', 'core_HA', 'core_FxR', 'core_RxxxxxP', 'core_KxxxxxS', 'core_TxxxxT', 'core_GxxxxxP', 'core_PxxxxH', 'core_GxxxxW', 'core_AxxxxG', 'core_TxP', 'core_VxxxxS', 'core_RxxV', 'core_NxA', 'core_HxxxP', 'core_TxxxxD', 'core_LxxxxA', 'core_FxxxxP', 'core_SxI', 'core_EC', 'core_RT', 'core_GxxxxxR', 'core_AxK', 'core_CxxxR', 'core_LN', 'core_QxxxxxP', 'core_AxxxxR', 'core_IxxxxxQ', 'core_WA', 'core_WxR', 'core_AP', 'core_LxN', 'core_PxxV', 'core_SxxxV', 'core_CxxxxxG', 'core_EG', 'core_AW', 'core_PxxxQ', 'core_PD', 'core_FT', 'core_RxxxxxS', 'core_FxxxxS', 'core_LxxxxI', 'core_RxC', 'core_LH', 'core_SxxxxF', 'core_MR', 'core_FxxxxxR', 'core_KxxA', 'core_GxG', 'core_PxxxR', 'core_SxxxxR', 'core_IG', 'core_PxxxxxS', 'core_T', 'core_WxG', 'core_AxxxxM', 'core_AxxA', 'core_PxS', 'core_QxT', 'core_GxxxD', 'core_RxxxxL', 'core_GI', 'core_P', 'core_SxS', 'core_CxxxxP', 'core_QxxxR', 'core_SxxxA', 'core_PF', 'core_ES', 'core_IP', 'core_LxxxN', 'core_VxxxxP', 'core_NxP', 'core_SxxxxxH', 'core_HxR', 'core_LxxxxH', 'core_VxxxxxQ', 'core_CxxxxD', 'core_AH', 'core_AxxxxW', 'core_ExxxxxC', 'core_AxxE', 'core_LxxxV', 'core_DxxxxxC', 'core_GxW', 'core_SG', 'core_RxA', 'core_TxxxxxL', 'core_ExxxR', 'core_CK', 'core_PxxxxT', 'core_CxxxG', 'core_WP', 'core_PxK', 'core_CD', 'core_LxxL', 'core_GxxV', 'core_DxS', 'core_AA', 'core_AxL', 'core_DS', 'core_HxxR', 'core_ExxG', 'core_GxxK', 'core_TxxxxY', 'core_GxxxP', 'core_C', 'core_GxxxV', 'core_PxxxxxV', 'core_WT', 'core_PxxxxxQ', 'core_LV', 'core_GA', 'core_VxL', 'core_KxA', 'core_RxxxxxK', 'core_GR', 'core_QxxxxR', 'core_RxxK', 'core_VxxxxV', 'core_TxxC', 'core_SC', 'core_AxxxxxQ', 'core_RxxD', 'core_WxxxxV', 'core_VxA', 'core_NxxxI', 'core_RxxxxH', 'core_TN', 'core_CxxxS', 'core_KxxR', 'core_GxxxW', 'core_TC', 'core_VxxC', 'core_RxxxxxM', 'core_RxxxxxH', 'core_RxxxxD', 'core_VL', 'core_CxxxxL', 'core_HH', 'core_KxxxxA', 'core_CxxxxxW', 'core_AD', 'core_TxxxxxD', 'core_YxS', 'core_SxxxxW', 'core_GW', 'core_ExxxC', 'core_WxxxxS', 'core_AxV', 'core_PxxxxF', 'core_TxxR', 'core_AxxxxxH', 'core_MxxS', 'core_SxxxxxR', 'core_HxxxxxP', 'core_SxxA', 'core_AxxM', 'core_SR', 'core_VxD', 'core_FxxP', 'core_DxG', 'core_TxxV', 'core_AxP', 'core_HxxS', 'core_FxxxC', 'core_IxxA', 'core_RW', 'core_SxxxxP', 'core_HP', 'core_AS', 'core_HS', 'core_GxxxL', 'core_KS', 'core_QxxS', 'core_RxxxxxF', 'core_HxxxR', 'core_HxxxV', 'core_SxxxxxA', 'core_CI', 'core_AxxxxxT', 'core_CxxT', 'core_WxxxxP', 'core_MxxxS', 'core_ExxxP', 'core_KxxxxxC', 'core_FxxxI', 'core_PxV', 'core_GxxxxxD', 'core_FxxxP', 'core_MV', 'core_RxR', 'core_GxxxxL', 'core_CV', 'core_LG', 'core_LK', 'core_S', 'core_SxxxxT', 'core_VQ', 'core_PxxL', 'core_VxxxG', 'core_PxxxF', 'core_LxL', 'core_RxxxxxL', 'core_CxR', 'core_TV', 'core_AxxxW', 'core_NxxR', 'core_AxxxM', 'core_AxxxE', 'core_CxxxxxA', 'core_KxxxxT', 'core_TxC', 'core_PH', 'core_RxxG', 'core_NxxxxxI', 'core_AxxxxxR', 'core_HxxD', 'core_RxxxR', 'core_RG', 'core_HxxxxP', 'core_RxxxY', 'core_KxxxxS', 'core_GxxA', 'core_RxxxxxR', 'core_H', 'core_HxT', 'core_CxG', 'core_GxxxxD', 'core_QxxL', 'core_LxxxC', 'core_GxxR', 'core_VxxxxQ', 'core_GxI', 'core_RxxxN', 'core_PP', 'core_PV', 'core_HxxV', 'core_PxA', 'core_AK', 'core_LF', 'core_CR', 'core_IxA', 'core_KxxS', 'core_KxxxP', 'core_CxxP', 'core_KG', 'core_PxxxxD', 'core_QY', 'core_LxxF', 'core_QxxxG', 'core_RxxxS', 'core_CxxxxxL', 'core_GE', 'core_HxxxxV', 'core_RxxxxxG', 'core_VxxxxL', 'core_MxxR', 'core_VxR', 'core_VxT', 'core_NxxxxxR', 'core_RxxxxP', 'core_YR', 'core_PxP', 'core_PxC', 'core_PxN', 'core_DxxxxxS', 'core_VxS', 'core_LxA', 'core_VxxxxxG', 'core_NxxxxxV', 'core_VH', 'core_RxW', 'core_AxxD', 'core_YxxF', 'core_QR', 'core_RxxE', 'core_RxxxxG', 'core_RxxF', 'core_SxxF', 'core_GxxxxP', 'core_VxxH', 'core_DV', 'core_VxxxA', 'core_RxxxxxE', 'core_RxxxxT', 'core_G', 'core_VxxxxxA', 'core_MxxxxA', 'core_FxxxS', 'core_SxxxxH', 'core_QxxxP', 'core_AxxxA', 'core_HxS', 'core_IxxS', 'core_AxxS', 'core_PxxW', 'core_TxxxxA', 'core_W', 'core_IE', 'core_RxN', 'core_HxxxS', 'core_MxxxxxR', 'core_KxxxF', 'core_LxxxxxA', 'core_LxxxE', 'core_NxxxS', 'core_CA', 'core_TxxxxP', 'core_LxP', 'core_GxxxxC', 'core_RM', 'core_GV', 'core_IR', 'core_DxxG', 'core_GxxxxxF', 'core_IxxxxQ', 'core_GxxxxQ', 'core_HxxxxG', 'core_ExxxxL', 'core_ExxxxT', 'core_CxA', 'core_HxxA', 'core_KxL', 'core_WxxxI', 'core_GxF', 'core_TxxxxW', 'core_WxxxxxG', 'core_AxxH', 'core_FxxxxN', 'core_QxxR', 'core_SxP', 'core_TxxxxxG', 'core_SxxxG', 'core_VxxL', 'core_SxxH', 'core_RS', 'core_PxxxxI', 'core_CL', 'core_DxxxxxA', 'core_IxxxxA', 'core_WxxxQ', 'core_IxxE', 'core_IxxxxxP', 'core_KR', 'core_GxP', 'core_ExA', 'core_TxxxC', 'core_QxG', 'core_GxxxxxG', 'core_LxxxxR', 'core_KxxxR', 'core_SxxxxL', 'core_VR', 'core_AxM', 'core_HxxP', 'core_HxxxxxA', 'core_PxxxxR', 'core_HxxxG', 'core_TxxxxF', 'core_HR', 'core_CxxxxW', 'core_AxxxxxC', 'core_ExxxxxP', 'core_HxxxxR', 'core_PxxxxxH', 'core_RxxxT', 'core_VF', 'core_LxxxxxL', 'core_PxxxH', 'core_QxxxxE', 'core_VxxxxG', 'core_GM', 'core_GxxxxxV', 'core_PR', 'core_KxG', 'core_GxxH', 'core_ExxxxR', 'core_RxxxxxA', 'core_FL', 'core_QxxxxxR', 'core_IS', 'core_RxH', 'core_CxxxxG', 'core_RxE', 'core_DxxR', 'core_TxxxR', 'core_GG', 'core_PxxxxxA', 'core_AxxxQ', 'core_FxxA', 'core_LxxA', 'core_RxG', 'core_FxxxxT', 'core_CxxxxxI', 'core_VxxP', 'core_PxxxxA', 'core_QxP', 'core_PxxxxQ', 'core_WxxxK', 'core_GxxP', 'core_HxxxxA', 'core_AM', 'core_SxxxC', 'core_PxxP', 'core_GxxT', 'core_GxxxxT', 'core_MxxxR', 'core_IxxxxxR', 'core_PxxxG', 'core_RxxQ', 'core_VxxI', 'core_IxxP', 'core_IxxC', 'core_LxxxxxT', 'core_GxxS', 'core_GxxxxG', 'core_PxxxL', 'core_FxxL', 'core_QxxxxxG', 'core_RxxS', 'core_TxxF', 'core_PxxH', 'core_WxxR', 'core_PxxA', 'core_HxxxxS', 'core_GxE', 'core_NT', 'core_ST', 'core_WxxxP', 'core_RxxxxA', 'core_RxxxxxW', 'core_CE', 'core_RxF', 'core_VxH', 'core_TxxxT', 'core_CxP', 'core_LxxxP', 'core_ExR', 'core_KxxxxP', 'core_GxxW', 'core_PxxG', 'core_TxxxxxR', 'core_VxxxxK', 'core_ExxxxxR', 'core_CxxxP', 'core_TxE', 'core_ExxxxS', 'core_RV', 'core_AxxL', 'core_SxT', 'core_FA', 'core_AxxxxE', 'core_GxxxxxS', 'core_SxxxM', 'core_AxxV', 'core_WxxxC', 'core_TxxxV', 'core_LxxxxxV', 'core_RxxxI', 'core_VxxxC', 'core_LL', 'core_LxR', 'core_IxxxP', 'core_TxxxxR', 'core_VxxV', 'core_PxQ', 'core_SxxS', 'core_YT', 'core_GxxxxxH', 'core_SxxxxxD', 'core_IxR', 'core_TxxxF', 'core_QxxxS', 'core_CxxxT', 'core_DxxxA', 'core_PL', 'core_MP', 'core_VxxxxxV', 'core_YC', 'core_CxxxC', 'core_PxR', 'core_TxxxE', 'core_WC', 'core_VxxxR', 'core_CxC', 'core_LxxxxT', 'core_PxxxxxC', 'core_WxxxxxP', 'core_LxxxA', 'core_CxxV', 'core_IxxxxR', 'core_VxxxT', 'core_CxxxxA', 'core_LxS', 'core_DxR', 'core_VxxxxxP', 'core_RxD', 'core_WxA', 'core_TL', 'core_RxxxxW', 'core_AxxP', 'core_WxxxA', 'core_CxxG', 'core_AT', 'core_TxT', 'core_RR', 'core_HxG', 'core_RL', 'core_SxxE', 'core_AxxxxxA', 'core_IT', 'core_RxxH', 'core_AxxI', 'core_AxxxxS', 'core_FxL', 'core_SxxxF', 'core_TxxxxxP', 'core_PxxxD', 'core_RxS', 'core_SxxxxxP', 'core_YxC', 'core_KxxxxR', 'core_WxxxxG', 'core_TI', 'core_RxxxV', 'core_ER', 'core_AxxxxxL', 'core_GxxxE', 'core_SxD', 'core_LxxK', 'core_PxxxxxL', 'core_LxxxxN', 'core_AxxxxH', 'core_PA', 'core_PC', 'core_GxH', 'core_SxR', 'core_LxxxR', 'core_QQ', 'core_CxxxxxE', 'core_IxG', 'core_RxxA', 'core_ExQ', 'core_AxxxxD', 'core_GxV', 'core_KxxxxxA', 'core_PxE', 'core_WS', 'core_AxxxT', 'core_DxP', 'core_ExxxxC', 'core_GxxxY', 'core_QxxxxP', 'core_AxxxxQ', 'core_YxxxxxL', 'core_SxxxxxK', 'core_HxxxxxV', 'core_RQ', 'core_QxxC', 'core_DG', 'core_RxM', 'core_RxxxxxI', 'core_RxxxxQ', 'core_WxxxxxI', 'core_TxxxK', 'core_AxxxxP', 'core_MxR', 'core_GxxxxxW', 'core_RxK', 'core_SxxxP', 'core_SS', 'core_RxxxxxT', 'core_RxxxF', 'core_CxxxV', 'core_DxxxT', 'core_PxxS', 'core_NxxxxR', 'core_ExxR', 'core_RxT', 'core_K', 'core_IxxQ', 'core_LxxR', 'core_TH', 'core_ExV', 'core_VxxxxxF', 'core_CxxxxxR', 'core_VxxxxxH', 'core_CxxxxT', 'core_YxR', 'core_QxxxxxA', 'core_GxxxxxE', 'core_DxN', 'core_ExL', 'core_GxxxxE', 'core_LxxxxC', 'core_RxP', 'core_ExxxxxT', 'core_CxxxxxP', 'core_QV', 'core_RK', 'core_AQ', 'core_TxxxxG', 'core_PxxxxV', 'core_TxxN', 'core_RxxM', 'core_VxxR', 'core_HV', 'core_SxxxxxM', 'core_CxxL', 'core_WxxxG', 'core_TxxxP', 'core_DxxxxR', 'core_PxxxxxP', 'core_SxxxxxG', 'core_CxxxxI', 'core_HxP', 'core_GxxxR', 'core_VxxxxxL', 'core_PxxxA', 'core_LxxxxP', 'core_CxxxA', 'core_GxxxxH', 'core_LxxxF', 'core_QG', 'core_TxxW', 'core_WxxS', 'core_YxxxxxR', 'core_PxxxV', 'core_SxxxxS', 'core_LS', 'core_TxI', 'core_TxL', 'core_SxE', 'core_RxL', 'core_AV', 'core_I', 'core_RxxI', 'core_RxxL', 'core_RxQ', 'core_VxxxxW', 'core_ExxxG', 'core_ExxxxxL', 'core_TxxI', 'core_FxxxA', 'core_PxY', 'core_WxV', 'core_AG', 'core_LxxxxG', 'core_AxA', 'core_LxxV', 'core_AR', 'core_SxV', 'core_GxxxG', 'core_AxxxR', 'core_EP', 'core_DxxxxP', 'core_LR', 'core_NxR', 'core_DxxL', 'core_RxxxE', 'core_SxxxxxI', 'core_AxH', 'core_VxxW', 'core_PQ', 'core_PxxxP', 'core_VxxxxxY', 'core_RxxxxE', 'core_FxP', 'core_AxxxxxS', 'core_FR', 'core_SxxV', 'core_NxV', 'core_LxxxxxI', 'core_CxxS', 'core_AxxxP', 'core_ExxV', 'core_QxxxxI', 'core_TW', 'core_RxxxxxV', 'core_VxxxxD', 'core_RP', 'core_KP', 'core_CxI', 'core_CxY', 'core_LxxxxxS', 'core_SxxxxA', 'core_LxxxG', 'core_DA', 'core_FxA', 'core_DxxxxxV', 'core_RxxxH', 'core_RxxxP', 'core_SxxG', 'core_QxxxxS', 'core_DxxxV', 'core_AxxxxK', 'core_RxxxxC', 'core_FxxR', 'core_AxxxxF', 'core_DxxxxxR', 'core_FC', 'core_FxxxxE', 'core_EQ', 'core_CxW', 'core_ET', 'core_AI', 'core_TxxxxxC', 'core_PxT', 'core_PxxD', 'core_GL', 'core_PxxxS', 'core_MxxA', 'core_RxxxxI', 'core_AxxxS', 'core_RxxxxxQ', 'core_SM', 'core_RC', 'core_SxL', 'core_SA', 'core_DxW', 'core_PxxxE', 'core_VxxT', 'core_PxxxC', 'core_AxS', 'core_DxxxS', 'core_AxG', 'core_LxxY', 'core_CT', 'core_QxA', 'core_RxxT', 'core_SxxP', 'core_Q', 'core_ExxL', 'core_PE', 'core_LxxxL', 'core_PxD', 'core_NR', 'core_ExxxA', 'core_FxxxxG', 'core_WxC', 'core_PK', 'core_ExxxL', 'core_QxR', 'core_RxV', 'core_RxxxM', 'core_QxxxA', 'core_VG', 'core_FxxxR', 'core_LxxC', 'core_LxxxxS', 'core_TxxxxQ', 'core_VxxxI', 'core_PxxR', 'core_RxI', 'core_LxK', 'core_RA', 'core_IxT', 'core_RxxxQ', 'core_GxL', 'core_ExxS', 'core_CxxA', 'core_PT', 'core_FxxxxxC', 'core_RxxxxxC', 'core_VxxF', 'core_QxF', 'core_RxxxxN', 'core_DR', 'core_GxxxS', 'core_DE', 'core_MxxxxR', 'core_IxI', 'core_FxxD', 'core_FV', 'core_NxxxV', 'core_VxxxS', 'core_M', 'core_HxxL', 'core_NxxxxxC', 'core_VT', 'core_QxxxxA', 'core_V', 'core_WxxxR', 'core_SP', 'core_AxR', 'core_VxxxxxR', 'core_WxxxxR', 'core_SL', 'core_RxxxL', 'core_TxxxxC', 'core_HL', 'core_CxxN', 'core_SI', 'core_HxxxxxG', 'core_MS', 'core_PxG', 'core_CxxxxxD', 'core_HxxxA', 'core_MW', 'core_SxxxW', 'core_LxxxxxP', 'core_IxxxE', 'core_RxxR', 'core_RxxxG', 'core_RI', 'core_GxxxH', 'core_RxxxA', 'core_PxxI', 'core_CxxR', 'core_FxxxxxL', 'core_VxxxH', 'core_SxG', 'core_SxxT', 'core_RxxxxR', 'core_RxxxxxN', 'core_QxV']
#        
#        
        for seq in core_seqs:
            regex = re.compile(seq[5:].replace("x", "."))
            scoring_csv_columns.append(len(regex.findall(self.core)))
            scoring_column_names.append(seq)
        core_pi = ProteinAnalysis(self.core).isoelectric_point()
        scoring_csv_columns.append(core_pi)
        scoring_column_names.append("core_pI")
        
        if c_count >= 2:
            self.score += 4
            
        if self.sequence[-3:] == "GCD":
            self.score += 2
        if self.sequence[-3:] == "LGS":
            self.score += 2
        
        scoring_csv_columns = [x for _, x in sorted(zip(scoring_column_names, scoring_csv_columns))]
#        print(sorted(scoring_column_names))
        self.csv_columns += [self.score] + scoring_csv_columns
        
                
  
