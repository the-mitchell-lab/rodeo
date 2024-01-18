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

peptide_type = "lanthi1"
CUTOFF = 20
index = 0

def write_csv_headers(output_dir, meta=False):
    dir_prefix = output_dir + '/{}/'.format(peptide_type)
    if not os.path.exists(dir_prefix):
        os.makedirs(dir_prefix)
    svm_headers = 'PK,Classification, 1/5_frac_c, 1/5_frac_s_t, 1/5_frac_s_t_c, 2/5_frac_c, 2/5_frac_s_t, 2/5_frac_s_t_c, 3/5_frac_c, 3/5_frac_s_t, 3/5_frac_s_t_c, 4/5_frac_s_t, 4/5_frac_s_t_c, 5/5_frac_c, 5/5_frac_s_t, 5/5_frac_s_t_c, S+T+C, aliphatic, aromatic, charged, core_A, core_AC, core_AE, core_AF, core_AH, core_AK, core_AL, core_AM, core_AP, core_AR, core_AT, core_AY, core_AxC, core_AxD, core_AxE, core_AxH, core_AxL, core_AxM, core_AxN, core_AxP, core_AxQ, core_AxR, core_AxT, core_AxxD, core_AxxH, core_AxxI, core_AxxL, core_AxxP, core_AxxR, core_AxxT, core_AxxW, core_AxxxD, core_AxxxE, core_AxxxH, core_AxxxI, core_AxxxL, core_AxxxM, core_AxxxP, core_AxxxQ, core_AxxxR, core_AxxxT, core_AxxxW, core_AxxxxC, core_AxxxxE, core_AxxxxF, core_AxxxxH, core_AxxxxI, core_AxxxxL, core_AxxxxP, core_AxxxxR, core_AxxxxT, core_AxxxxV, core_AxxxxxC, core_AxxxxxD, core_AxxxxxE, core_AxxxxxF, core_AxxxxxH, core_AxxxxxK, core_AxxxxxL, core_AxxxxxP, core_AxxxxxQ, core_AxxxxxR, core_C, core_CC, core_CD, core_CE, core_CF, core_CG, core_CK, core_CL, core_CM, core_CN, core_CP, core_CQ, core_CR, core_CT, core_CV, core_CY, core_CxE, core_CxF, core_CxG, core_CxH, core_CxN, core_CxQ, core_CxR, core_CxS, core_CxT, core_CxW, core_CxxC, core_CxxE, core_CxxF, core_CxxG, core_CxxH, core_CxxL, core_CxxN, core_CxxP, core_CxxQ, core_CxxR, core_CxxS, core_CxxT, core_CxxW, core_CxxxC, core_CxxxD, core_CxxxE, core_CxxxF, core_CxxxK, core_CxxxL, core_CxxxP, core_CxxxR, core_CxxxT, core_CxxxV, core_CxxxxA, core_CxxxxC, core_CxxxxE, core_CxxxxG, core_CxxxxH, core_CxxxxL, core_CxxxxM, core_CxxxxN, core_CxxxxP, core_CxxxxQ, core_CxxxxR, core_CxxxxS, core_CxxxxT, core_CxxxxW, core_CxxxxY, core_CxxxxxC, core_CxxxxxD, core_CxxxxxE, core_CxxxxxF, core_CxxxxxG, core_CxxxxxH, core_CxxxxxK, core_CxxxxxL, core_CxxxxxM, core_CxxxxxN, core_CxxxxxP, core_CxxxxxQ, core_CxxxxxR, core_CxxxxxT, core_D, core_DA, core_DC, core_DE, core_DF, core_DH, core_DI, core_DK, core_DL, core_DP, core_DQ, core_DR, core_DS, core_DV, core_DxA, core_DxF, core_DxG, core_DxH, core_DxI, core_DxL, core_DxN, core_DxQ, core_DxR, core_DxS, core_DxV, core_DxW, core_DxxD, core_DxxE, core_DxxF, core_DxxH, core_DxxI, core_DxxL, core_DxxP, core_DxxQ, core_DxxR, core_DxxS, core_DxxV, core_DxxW, core_DxxxA, core_DxxxD, core_DxxxE, core_DxxxF, core_DxxxG, core_DxxxI, core_DxxxL, core_DxxxN, core_DxxxP, core_DxxxQ, core_DxxxR, core_DxxxS, core_DxxxV, core_DxxxW, core_DxxxxA, core_DxxxxC, core_DxxxxD, core_DxxxxE, core_DxxxxF, core_DxxxxG, core_DxxxxH, core_DxxxxL, core_DxxxxR, core_DxxxxT, core_DxxxxxA, core_DxxxxxD, core_DxxxxxE, core_DxxxxxH, core_DxxxxxK, core_DxxxxxL, core_DxxxxxN, core_DxxxxxQ, core_DxxxxxR, core_DxxxxxS, core_DxxxxxT, core_DxxxxxY, core_E, core_EA, core_EC, core_EE, core_EF, core_EG, core_EH, core_EI, core_EK, core_EL, core_EM, core_EN, core_EP, core_EQ, core_ER, core_ES, core_ET, core_EV, core_EY, core_ExA, core_ExC, core_ExD, core_ExE, core_ExF, core_ExG, core_ExI, core_ExK, core_ExL, core_ExM, core_ExP, core_ExQ, core_ExR, core_ExS, core_ExT, core_ExV, core_ExxC, core_ExxD, core_ExxE, core_ExxH, core_ExxI, core_ExxK, core_ExxL, core_ExxM, core_ExxP, core_ExxQ, core_ExxR, core_ExxxD, core_ExxxF, core_ExxxH, core_ExxxL, core_ExxxN, core_ExxxP, core_ExxxR, core_ExxxT, core_ExxxW, core_ExxxxA, core_ExxxxC, core_ExxxxE, core_ExxxxF, core_ExxxxG, core_ExxxxH, core_ExxxxI, core_ExxxxK, core_ExxxxL, core_ExxxxN, core_ExxxxQ, core_ExxxxR, core_ExxxxT, core_ExxxxV, core_ExxxxxD, core_ExxxxxE, core_ExxxxxG, core_ExxxxxH, core_ExxxxxI, core_ExxxxxN, core_ExxxxxR, core_ExxxxxS, core_ExxxxxT, core_ExxxxxW, core_F, core_FA, core_FD, core_FE, core_FF, core_FH, core_FI, core_FK, core_FL, core_FM, core_FN, core_FP, core_FQ, core_FR, core_FS, core_FT, core_FY, core_FxD, core_FxE, core_FxF, core_FxH, core_FxI, core_FxK, core_FxL, core_FxM, core_FxN, core_FxP, core_FxQ, core_FxR, core_FxS, core_FxV, core_FxY, core_FxxA, core_FxxD, core_FxxE, core_FxxF, core_FxxH, core_FxxI, core_FxxK, core_FxxL, core_FxxM, core_FxxN, core_FxxP, core_FxxQ, core_FxxR, core_FxxS, core_FxxV, core_FxxY, core_FxxxA, core_FxxxE, core_FxxxF, core_FxxxG, core_FxxxH, core_FxxxI, core_FxxxK, core_FxxxL, core_FxxxM, core_FxxxN, core_FxxxP, core_FxxxQ, core_FxxxR, core_FxxxS, core_FxxxV, core_FxxxY, core_FxxxxA, core_FxxxxD, core_FxxxxE, core_FxxxxF, core_FxxxxH, core_FxxxxI, core_FxxxxK, core_FxxxxL, core_FxxxxM, core_FxxxxN, core_FxxxxP, core_FxxxxQ, core_FxxxxR, core_FxxxxS, core_FxxxxT, core_FxxxxV, core_FxxxxY, core_FxxxxxA, core_FxxxxxD, core_FxxxxxF, core_FxxxxxG, core_FxxxxxH, core_FxxxxxI, core_FxxxxxK, core_FxxxxxL, core_FxxxxxM, core_FxxxxxN, core_FxxxxxP, core_FxxxxxQ, core_FxxxxxR, core_FxxxxxS, core_FxxxxxT, core_FxxxxxV, core_FxxxxxY, core_GC, core_GD, core_GE, core_GG, core_GH, core_GN, core_GP, core_GQ, core_GR, core_GV, core_GW, core_GY, core_GxA, core_GxC, core_GxD, core_GxE, core_GxG, core_GxH, core_GxI, core_GxK, core_GxL, core_GxN, core_GxP, core_GxQ, core_GxR, core_GxT, core_GxV, core_GxW, core_GxY, core_GxxA, core_GxxC, core_GxxG, core_GxxH, core_GxxI, core_GxxL, core_GxxN, core_GxxP, core_GxxR, core_GxxT, core_GxxV, core_GxxW, core_GxxY, core_GxxxC, core_GxxxD, core_GxxxE, core_GxxxG, core_GxxxH, core_GxxxL, core_GxxxN, core_GxxxQ, core_GxxxR, core_GxxxS, core_GxxxT, core_GxxxY, core_GxxxxD, core_GxxxxE, core_GxxxxH, core_GxxxxI, core_GxxxxK, core_GxxxxL, core_GxxxxN, core_GxxxxP, core_GxxxxQ, core_GxxxxR, core_GxxxxT, core_GxxxxV, core_GxxxxW, core_GxxxxY, core_GxxxxxA, core_GxxxxxC, core_GxxxxxD, core_GxxxxxE, core_GxxxxxF, core_GxxxxxH, core_GxxxxxK, core_GxxxxxL, core_GxxxxxP, core_GxxxxxQ, core_GxxxxxR, core_GxxxxxS, core_GxxxxxT, core_GxxxxxY, core_H, core_HA, core_HC, core_HD, core_HE, core_HF, core_HG, core_HH, core_HI, core_HK, core_HL, core_HN, core_HP, core_HQ, core_HR, core_HS, core_HT, core_HV, core_HxA, core_HxC, core_HxD, core_HxE, core_HxF, core_HxG, core_HxH, core_HxI, core_HxK, core_HxL, core_HxP, core_HxQ, core_HxR, core_HxS, core_HxV, core_HxW, core_HxY, core_HxxA, core_HxxC, core_HxxD, core_HxxE, core_HxxF, core_HxxG, core_HxxI, core_HxxK, core_HxxL, core_HxxP, core_HxxQ, core_HxxR, core_HxxS, core_HxxV, core_HxxY, core_HxxxA, core_HxxxD, core_HxxxE, core_HxxxF, core_HxxxG, core_HxxxH, core_HxxxI, core_HxxxK, core_HxxxL, core_HxxxN, core_HxxxP, core_HxxxQ, core_HxxxR, core_HxxxS, core_HxxxT, core_HxxxV, core_HxxxW, core_HxxxxA, core_HxxxxC, core_HxxxxD, core_HxxxxE, core_HxxxxF, core_HxxxxG, core_HxxxxH, core_HxxxxI, core_HxxxxK, core_HxxxxL, core_HxxxxN, core_HxxxxP, core_HxxxxQ, core_HxxxxR, core_HxxxxS, core_HxxxxT, core_HxxxxV, core_HxxxxY, core_HxxxxxA, core_HxxxxxD, core_HxxxxxE, core_HxxxxxG, core_HxxxxxH, core_HxxxxxI, core_HxxxxxL, core_HxxxxxM, core_HxxxxxN, core_HxxxxxP, core_HxxxxxQ, core_HxxxxxR, core_HxxxxxS, core_HxxxxxV, core_HxxxxxY, core_I, core_IC, core_ID, core_IE, core_IF, core_IG, core_IH, core_IL, core_IN, core_IP, core_IQ, core_IR, core_IT, core_IY, core_IxA, core_IxD, core_IxE, core_IxF, core_IxG, core_IxL, core_IxM, core_IxN, core_IxQ, core_IxR, core_IxT, core_IxY, core_IxxA, core_IxxC, core_IxxE, core_IxxF, core_IxxG, core_IxxH, core_IxxI, core_IxxK, core_IxxL, core_IxxM, core_IxxN, core_IxxP, core_IxxR, core_IxxT, core_IxxY, core_IxxxC, core_IxxxD, core_IxxxE, core_IxxxF, core_IxxxH, core_IxxxL, core_IxxxM, core_IxxxN, core_IxxxP, core_IxxxQ, core_IxxxR, core_IxxxT, core_IxxxY, core_IxxxxE, core_IxxxxF, core_IxxxxG, core_IxxxxI, core_IxxxxL, core_IxxxxM, core_IxxxxN, core_IxxxxP, core_IxxxxR, core_IxxxxT, core_IxxxxY, core_IxxxxxC, core_IxxxxxF, core_IxxxxxG, core_IxxxxxH, core_IxxxxxK, core_IxxxxxL, core_IxxxxxP, core_IxxxxxQ, core_IxxxxxR, core_IxxxxxS, core_IxxxxxW, core_K, core_KA, core_KC, core_KD, core_KF, core_KG, core_KH, core_KI, core_KL, core_KN, core_KP, core_KQ, core_KR, core_KS, core_KY, core_KxA, core_KxC, core_KxD, core_KxE, core_KxF, core_KxH, core_KxI, core_KxL, core_KxM, core_KxN, core_KxP, core_KxQ, core_KxR, core_KxS, core_KxT, core_KxV, core_KxY, core_KxxA, core_KxxD, core_KxxE, core_KxxF, core_KxxH, core_KxxI, core_KxxL, core_KxxN, core_KxxR, core_KxxxE, core_KxxxF, core_KxxxG, core_KxxxI, core_KxxxK, core_KxxxL, core_KxxxM, core_KxxxN, core_KxxxQ, core_KxxxR, core_KxxxS, core_KxxxT, core_KxxxV, core_KxxxW, core_KxxxY, core_KxxxxA, core_KxxxxD, core_KxxxxE, core_KxxxxF, core_KxxxxH, core_KxxxxI, core_KxxxxL, core_KxxxxP, core_KxxxxQ, core_KxxxxR, core_KxxxxS, core_KxxxxV, core_KxxxxY, core_KxxxxxC, core_KxxxxxD, core_KxxxxxE, core_KxxxxxF, core_KxxxxxH, core_KxxxxxI, core_KxxxxxK, core_KxxxxxL, core_KxxxxxN, core_KxxxxxP, core_KxxxxxQ, core_KxxxxxR, core_KxxxxxS, core_KxxxxxY, core_L, core_LA, core_LD, core_LE, core_LF, core_LG, core_LH, core_LI, core_LK, core_LL, core_LM, core_LN, core_LP, core_LQ, core_LR, core_LS, core_LT, core_LV, core_LW, core_LY, core_LxA, core_LxC, core_LxE, core_LxF, core_LxG, core_LxH, core_LxI, core_LxK, core_LxL, core_LxM, core_LxP, core_LxQ, core_LxR, core_LxS, core_LxV, core_LxW, core_LxY, core_LxxA, core_LxxF, core_LxxG, core_LxxH, core_LxxI, core_LxxK, core_LxxL, core_LxxM, core_LxxN, core_LxxP, core_LxxQ, core_LxxR, core_LxxS, core_LxxT, core_LxxV, core_LxxxA, core_LxxxC, core_LxxxD, core_LxxxE, core_LxxxF, core_LxxxH, core_LxxxI, core_LxxxK, core_LxxxL, core_LxxxM, core_LxxxN, core_LxxxP, core_LxxxQ, core_LxxxR, core_LxxxS, core_LxxxV, core_LxxxY, core_LxxxxA, core_LxxxxD, core_LxxxxE, core_LxxxxF, core_LxxxxG, core_LxxxxH, core_LxxxxI, core_LxxxxK, core_LxxxxL, core_LxxxxN, core_LxxxxP, core_LxxxxQ, core_LxxxxR, core_LxxxxS, core_LxxxxT, core_LxxxxV, core_LxxxxY, core_LxxxxxA, core_LxxxxxD, core_LxxxxxE, core_LxxxxxF, core_LxxxxxG, core_LxxxxxH, core_LxxxxxI, core_LxxxxxK, core_LxxxxxL, core_LxxxxxM, core_LxxxxxP, core_LxxxxxQ, core_LxxxxxR, core_LxxxxxS, core_LxxxxxT, core_LxxxxxV, core_LxxxxxY, core_M, core_MF, core_MH, core_MI, core_MK, core_ML, core_MQ, core_MR, core_MV, core_MW, core_MY, core_MxE, core_MxF, core_MxH, core_MxI, core_MxL, core_MxP, core_MxR, core_MxT, core_MxW, core_MxY, core_MxxA, core_MxxC, core_MxxF, core_MxxH, core_MxxI, core_MxxL, core_MxxM, core_MxxP, core_MxxQ, core_MxxR, core_MxxS, core_MxxV, core_MxxW, core_MxxxF, core_MxxxL, core_MxxxN, core_MxxxP, core_MxxxQ, core_MxxxR, core_MxxxS, core_MxxxV, core_MxxxW, core_MxxxY, core_MxxxxA, core_MxxxxC, core_MxxxxI, core_MxxxxK, core_MxxxxL, core_MxxxxP, core_MxxxxQ, core_MxxxxR, core_MxxxxS, core_MxxxxxA, core_MxxxxxD, core_MxxxxxF, core_MxxxxxL, core_MxxxxxM, core_MxxxxxP, core_MxxxxxQ, core_MxxxxxR, core_MxxxxxS, core_NA, core_NC, core_ND, core_NE, core_NF, core_NG, core_NH, core_NI, core_NK, core_NL, core_NQ, core_NR, core_NT, core_NV, core_NxA, core_NxC, core_NxE, core_NxG, core_NxH, core_NxI, core_NxK, core_NxL, core_NxP, core_NxR, core_NxS, core_NxW, core_NxY, core_NxxE, core_NxxG, core_NxxH, core_NxxI, core_NxxL, core_NxxP, core_NxxR, core_NxxT, core_NxxxD, core_NxxxE, core_NxxxI, core_NxxxN, core_NxxxP, core_NxxxQ, core_NxxxR, core_NxxxS, core_NxxxT, core_NxxxW, core_NxxxxA, core_NxxxxC, core_NxxxxD, core_NxxxxF, core_NxxxxH, core_NxxxxK, core_NxxxxL, core_NxxxxN, core_NxxxxQ, core_NxxxxR, core_NxxxxT, core_NxxxxY, core_NxxxxxC, core_NxxxxxD, core_NxxxxxE, core_NxxxxxH, core_NxxxxxI, core_NxxxxxK, core_NxxxxxL, core_NxxxxxM, core_NxxxxxN, core_NxxxxxP, core_NxxxxxR, core_NxxxxxS, core_NxxxxxT, core_NxxxxxY, core_P, core_PA, core_PD, core_PE, core_PF, core_PG, core_PH, core_PK, core_PL, core_PP, core_PQ, core_PR, core_PS, core_PT, core_PV, core_PxA, core_PxC, core_PxD, core_PxE, core_PxF, core_PxG, core_PxH, core_PxK, core_PxL, core_PxM, core_PxP, core_PxQ, core_PxR, core_PxS, core_PxT, core_PxV, core_PxY, core_PxxA, core_PxxC, core_PxxD, core_PxxE, core_PxxF, core_PxxG, core_PxxH, core_PxxL, core_PxxP, core_PxxQ, core_PxxR, core_PxxS, core_PxxT, core_PxxW, core_PxxxA, core_PxxxD, core_PxxxE, core_PxxxG, core_PxxxH, core_PxxxK, core_PxxxL, core_PxxxM, core_PxxxN, core_PxxxQ, core_PxxxR, core_PxxxS, core_PxxxxC, core_PxxxxD, core_PxxxxE, core_PxxxxF, core_PxxxxG, core_PxxxxH, core_PxxxxK, core_PxxxxL, core_PxxxxM, core_PxxxxN, core_PxxxxP, core_PxxxxQ, core_PxxxxR, core_PxxxxS, core_PxxxxT, core_PxxxxV, core_PxxxxY, core_PxxxxxA, core_PxxxxxC, core_PxxxxxD, core_PxxxxxE, core_PxxxxxF, core_PxxxxxG, core_PxxxxxH, core_PxxxxxI, core_PxxxxxK, core_PxxxxxL, core_PxxxxxM, core_PxxxxxN, core_PxxxxxP, core_PxxxxxQ, core_PxxxxxR, core_PxxxxxS, core_PxxxxxW, core_PxxxxxY, core_Q, core_QD, core_QE, core_QH, core_QI, core_QK, core_QL, core_QM, core_QQ, core_QR, core_QS, core_QT, core_QV, core_QY, core_QxA, core_QxD, core_QxE, core_QxF, core_QxG, core_QxH, core_QxK, core_QxL, core_QxP, core_QxQ, core_QxR, core_QxS, core_QxV, core_QxW, core_QxY, core_QxxA, core_QxxC, core_QxxD, core_QxxE, core_QxxG, core_QxxH, core_QxxI, core_QxxK, core_QxxL, core_QxxM, core_QxxN, core_QxxP, core_QxxQ, core_QxxR, core_QxxS, core_QxxT, core_QxxV, core_QxxY, core_QxxxA, core_QxxxC, core_QxxxD, core_QxxxE, core_QxxxF, core_QxxxG, core_QxxxH, core_QxxxI, core_QxxxL, core_QxxxM, core_QxxxP, core_QxxxQ, core_QxxxR, core_QxxxS, core_QxxxT, core_QxxxV, core_QxxxY, core_QxxxxA, core_QxxxxC, core_QxxxxD, core_QxxxxE, core_QxxxxF, core_QxxxxG, core_QxxxxH, core_QxxxxL, core_QxxxxM, core_QxxxxN, core_QxxxxP, core_QxxxxQ, core_QxxxxR, core_QxxxxS, core_QxxxxV, core_QxxxxY, core_QxxxxxA, core_QxxxxxC, core_QxxxxxE, core_QxxxxxF, core_QxxxxxG, core_QxxxxxH, core_QxxxxxI, core_QxxxxxK, core_QxxxxxL, core_QxxxxxN, core_QxxxxxP, core_QxxxxxQ, core_QxxxxxR, core_QxxxxxS, core_QxxxxxY, core_R, core_RA, core_RC, core_RD, core_RE, core_RF, core_RG, core_RH, core_RI, core_RK, core_RL, core_RM, core_RN, core_RP, core_RQ, core_RR, core_RS, core_RT, core_RV, core_RW, core_RY, core_RxA, core_RxC, core_RxD, core_RxE, core_RxF, core_RxG, core_RxH, core_RxI, core_RxK, core_RxL, core_RxM, core_RxN, core_RxP, core_RxQ, core_RxR, core_RxS, core_RxT, core_RxV, core_RxW, core_RxY, core_RxxA, core_RxxC, core_RxxD, core_RxxE, core_RxxF, core_RxxG, core_RxxH, core_RxxI, core_RxxK, core_RxxL, core_RxxM, core_RxxN, core_RxxP, core_RxxQ, core_RxxR, core_RxxS, core_RxxT, core_RxxV, core_RxxW, core_RxxY, core_RxxxA, core_RxxxC, core_RxxxD, core_RxxxE, core_RxxxF, core_RxxxG, core_RxxxH, core_RxxxI, core_RxxxK, core_RxxxL, core_RxxxM, core_RxxxN, core_RxxxP, core_RxxxQ, core_RxxxR, core_RxxxS, core_RxxxT, core_RxxxV, core_RxxxW, core_RxxxY, core_RxxxxA, core_RxxxxC, core_RxxxxD, core_RxxxxE, core_RxxxxF, core_RxxxxG, core_RxxxxH, core_RxxxxI, core_RxxxxK, core_RxxxxL, core_RxxxxM, core_RxxxxN, core_RxxxxP, core_RxxxxQ, core_RxxxxR, core_RxxxxS, core_RxxxxT, core_RxxxxV, core_RxxxxW, core_RxxxxY, core_RxxxxxA, core_RxxxxxC, core_RxxxxxD, core_RxxxxxE, core_RxxxxxF, core_RxxxxxG, core_RxxxxxH, core_RxxxxxI, core_RxxxxxK, core_RxxxxxL, core_RxxxxxM, core_RxxxxxN, core_RxxxxxP, core_RxxxxxQ, core_RxxxxxR, core_RxxxxxS, core_RxxxxxT, core_RxxxxxV, core_RxxxxxW, core_RxxxxxY, core_S, core_SE, core_SF, core_SL, core_SM, core_SN, core_SP, core_SQ, core_SR, core_SS, core_ST, core_SW, core_SY, core_SxD, core_SxE, core_SxF, core_SxH, core_SxI, core_SxL, core_SxM, core_SxP, core_SxQ, core_SxR, core_SxS, core_SxT, core_SxV, core_SxW, core_SxxA, core_SxxC, core_SxxE, core_SxxF, core_SxxH, core_SxxI, core_SxxL, core_SxxM, core_SxxP, core_SxxQ, core_SxxR, core_SxxT, core_SxxY, core_SxxxA, core_SxxxC, core_SxxxE, core_SxxxF, core_SxxxG, core_SxxxH, core_SxxxI, core_SxxxK, core_SxxxL, core_SxxxM, core_SxxxN, core_SxxxP, core_SxxxQ, core_SxxxR, core_SxxxV, core_SxxxW, core_SxxxxC, core_SxxxxD, core_SxxxxE, core_SxxxxF, core_SxxxxH, core_SxxxxI, core_SxxxxL, core_SxxxxN, core_SxxxxP, core_SxxxxQ, core_SxxxxR, core_SxxxxV, core_SxxxxY, core_SxxxxxA, core_SxxxxxD, core_SxxxxxE, core_SxxxxxF, core_SxxxxxH, core_SxxxxxI, core_SxxxxxL, core_SxxxxxM, core_SxxxxxN, core_SxxxxxP, core_SxxxxxQ, core_SxxxxxR, core_SxxxxxS, core_SxxxxxT, core_SxxxxxY, core_T, core_TC, core_TE, core_TI, core_TK, core_TL, core_TM, core_TN, core_TP, core_TR, core_TS, core_TT, core_TV, core_TW, core_TxA, core_TxC, core_TxD, core_TxE, core_TxF, core_TxG, core_TxH, core_TxI, core_TxK, core_TxL, core_TxP, core_TxR, core_TxS, core_TxT, core_TxV, core_TxW, core_TxxA, core_TxxC, core_TxxH, core_TxxI, core_TxxK, core_TxxL, core_TxxP, core_TxxQ, core_TxxR, core_TxxT, core_TxxV, core_TxxxA, core_TxxxC, core_TxxxD, core_TxxxE, core_TxxxH, core_TxxxI, core_TxxxP, core_TxxxQ, core_TxxxR, core_TxxxS, core_TxxxT, core_TxxxV, core_TxxxxA, core_TxxxxC, core_TxxxxE, core_TxxxxG, core_TxxxxH, core_TxxxxI, core_TxxxxL, core_TxxxxR, core_TxxxxS, core_TxxxxT, core_TxxxxV, core_TxxxxxA, core_TxxxxxC, core_TxxxxxD, core_TxxxxxE, core_TxxxxxG, core_TxxxxxI, core_TxxxxxK, core_TxxxxxM, core_TxxxxxN, core_TxxxxxQ, core_TxxxxxR, core_TxxxxxS, core_TxxxxxT, core_TxxxxxV, core_VC, core_VE, core_VF, core_VH, core_VK, core_VL, core_VM, core_VP, core_VQ, core_VR, core_VS, core_VT, core_VW, core_VY, core_VxA, core_VxC, core_VxD, core_VxF, core_VxH, core_VxL, core_VxM, core_VxN, core_VxP, core_VxQ, core_VxR, core_VxT, core_VxV, core_VxxD, core_VxxE, core_VxxF, core_VxxH, core_VxxK, core_VxxP, core_VxxR, core_VxxT, core_VxxxD, core_VxxxE, core_VxxxF, core_VxxxH, core_VxxxI, core_VxxxM, core_VxxxP, core_VxxxQ, core_VxxxR, core_VxxxT, core_VxxxV, core_VxxxY, core_VxxxxA, core_VxxxxD, core_VxxxxE, core_VxxxxF, core_VxxxxL, core_VxxxxP, core_VxxxxQ, core_VxxxxR, core_VxxxxS, core_VxxxxT, core_VxxxxY, core_VxxxxxC, core_VxxxxxD, core_VxxxxxF, core_VxxxxxP, core_VxxxxxQ, core_VxxxxxR, core_VxxxxxT, core_WC, core_WG, core_WL, core_WN, core_WQ, core_WR, core_WT, core_WxC, core_WxF, core_WxL, core_WxR, core_WxT, core_WxxI, core_WxxL, core_WxxR, core_WxxT, core_WxxV, core_WxxxC, core_WxxxF, core_WxxxG, core_WxxxH, core_WxxxI, core_WxxxP, core_WxxxR, core_WxxxT, core_WxxxxH, core_WxxxxK, core_WxxxxP, core_WxxxxR, core_WxxxxT, core_WxxxxxC, core_WxxxxxE, core_WxxxxxK, core_WxxxxxR, core_WxxxxxV, core_Y, core_YC, core_YD, core_YF, core_YH, core_YI, core_YK, core_YL, core_YP, core_YQ, core_YR, core_YxC, core_YxD, core_YxF, core_YxG, core_YxH, core_YxI, core_YxK, core_YxL, core_YxN, core_YxP, core_YxQ, core_YxR, core_YxV, core_YxY, core_YxxA, core_YxxF, core_YxxG, core_YxxH, core_YxxI, core_YxxL, core_YxxP, core_YxxQ, core_YxxR, core_YxxY, core_YxxxC, core_YxxxE, core_YxxxF, core_YxxxH, core_YxxxI, core_YxxxK, core_YxxxL, core_YxxxP, core_YxxxQ, core_YxxxR, core_YxxxxF, core_YxxxxG, core_YxxxxH, core_YxxxxI, core_YxxxxL, core_YxxxxM, core_YxxxxN, core_YxxxxP, core_YxxxxQ, core_YxxxxR, core_YxxxxxF, core_YxxxxxG, core_YxxxxxH, core_YxxxxxI, core_YxxxxxK, core_YxxxxxL, core_YxxxxxM, core_YxxxxxP, core_YxxxxxQ, core_YxxxxxR, core_YxxxxxS, core_YxxxxxT, core_YxxxxxV, core_length, core_pI, frac_aliphatic, frac_aromatic, frac_charge, frac_polar, frac_pos, leader_length, neg, peptide_length, polar, pos'

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
        #f_scores = self.get_fimo_score("{}_leader_fimo.txt".format("./ripp_modules/lanthi1/lanthi1"))
        f_scores = self.get_fimo_score("{}/{}_leader_fimo.txt".format(FILE_DIR, peptide_type))  
        l_stop = sorted([0] + f_scores[2])[-1]
        if len(f_scores[0]) > 0:
            self.score += 10
        regex1 = re.search('(.[S|T].{2,7}C)', self.sequence[l_stop:])
        regex2 = re.search('(G[G|A])', self.sequence[l_stop:])
    
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
        fifth_size = len(self.sequence) // 5
        for i in range(5):
            scoring_csv_columns.append(self.sequence[i*fifth_size:(i+1)*fifth_size].count('C') / len(self.sequence))
            scoring_column_names.append("{}/5_frac_c".format(i+1))
        for i in [0, 1, 3, 4]:
            count = 0
            for c in "STC":
                count += self.sequence[i*fifth_size:(i+1)*fifth_size].count(c)
            scoring_csv_columns.append(count / len(self.sequence))
            scoring_column_names.append("{}/5_frac_s_t_c".format(i+1))
            
        for i in [2, 3, 4]:
            count = 0
            for c in "ST":
                count += self.sequence[i*fifth_size:(i+1)*fifth_size].count(c)
            scoring_csv_columns.append(count / len(self.sequence))
            scoring_column_names.append("{}/5_frac_s_t".format(i+1))
           
        scoring_csv_columns.append((self.core.count("S") + self.core.count("T")) / len(self.sequence))
        scoring_column_names.append("S+T")
        
        scoring_csv_columns.append((self.core.count("S") + self.core.count("T") + self.core.count("C"))/ len(self.sequence))
        scoring_column_names.append("S+T+C")
        
        aliphatic_core_count = self.core.count("A") + self.core.count("I") + self.core.count("M") + self.core.count("L") + self.core.count("V")
        scoring_csv_columns.append(aliphatic_core_count)
        scoring_column_names.append("aliphatic")
        scoring_csv_columns.append(aliphatic_core_count / len(self.core))
        scoring_column_names.append("frac_aliphatic")
        
        aromatic_core_count = self.core.count("F") + self.core.count("W") + self.core.count("Y")
        scoring_csv_columns.append(aromatic_core_count)
        scoring_column_names.append("aromatic")
        scoring_csv_columns.append(aromatic_core_count  / len(self.core))
        scoring_column_names.append("frac_aromatic")
        
        charged_count = self.core.count("D") + self.core.count("E") + self.core.count("K") + self.core.count("R")
        scoring_csv_columns.append(charged_count)
        scoring_column_names.append("charged")
        scoring_csv_columns.append(charged_count  / len(self.core))
        scoring_column_names.append("frac_charge")
        
        polar_count = self.core.count("C") + self.core.count("S") + self.core.count("T") + self.core.count("H") + self.core.count("N") + self.core.count("Q")
        scoring_csv_columns.append(polar_count)
        scoring_column_names.append("polar")
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
        scoring_csv_columns.append(neg / len(self.core))
        scoring_column_names.append("frac_neg")
        
        scoring_csv_columns.append(len(self.leader))
        scoring_column_names.append("leader_length")
        
        scoring_csv_columns.append(len(self.sequence))
        scoring_column_names.append("peptide_length")
        
        scoring_csv_columns.append(len(self.core))
        scoring_column_names.append("core_length")
        
        core_seqs = ["core_CxxxC","core_C","core_TC","core_CxxxxC","core_TxxxxC","core_TxxxC","core_CxT","core_TxxC","core_CT","core_CxxT","core_SxC","core_TxxxxxC","core_TxxxxxT","core_DxC","core_SC","core_DxxxxxC","core_CxxxxxT","core_T","core_CxxxxxN","core_R","core_GxT","core_TxxxT","core_TxxT","core_SxxxxxC","core_TxD","core_CxxG","core_CxxxT","core_CK","core_CxxxxT","core_GxxC","core_CxS","core_DxxxxT","core_TxxxxT","core_CN","core_CxxxxN","core_TxG","core_NxxC","core_CxxxxxC","core_NxxxC","core_CxxC","core_CxC","core_CC","core_TxN","core_SxxxC","core_SxxC","core_CxxS","core_GxxN","core_DN","core_NT","core_CxxN","core_DxxC","core_TxxxN","core_CG","core_GC","core_CxxxxG","core_TxxG","core_L","core_GxxxT","core_NC","core_GxxxxC","core_TxT","core_DxxxxS","core_SxN","core_GxxxxxC","core_SxxxxT","core_SxxxxC","core_AR","core_CxN","core_RxR","core_KT","core_RxxR","core_RxxxR","core_RR","core_RS","core_TD","core_TxC","core_TxxN","core_TxxS","core_RxxA","core_TxV","core_RxxxxxR","core_AxxxxR","core_RxxxxR","core_IxxxxxC","core_AxR","core_DxxxC","core_TxxxxxG","core_GxxxC","core_RxxxA","core_CxxD","core_CxxxxxI","core_AxxR","core_SxxxR","core_AxxxxxR","core_RxxG","core_TT","core_RxG","core_SxxR","core_AxxxR","core_RxA","core_CxK","core_PxxxR","core_PxR","core_GxxR","core_PxxR","core_GxR","core_RP","core_P","core_KxxC","core_NxxxxxT","core_RxxxS","core_RxxxxxA","core_RA","core_RxxP","core_TN","core_RxP","core_N","core_RxxxxG","core_RxxS","core_GxC","core_CxxxxxS","core_RxxxxA","core_NxC","core_GR","core_RxxxxxP","core_CS","core_SxR","core_SxxxxR","core_SxxxT","core_TxxxxG","core_CxxxG","core_RxS","core_TxxxxN","core_GxxxR","core_PxxxxR","core_RxxxxS","core_SR","core_GxxxxR","core_RG","core_NxxxxC","core_RxxxP","core_NxxxxxG","core_RxxxG","core_RxxxxxS","core_CxxxD","core_KxC","core_PxxxxxR","core_RxxxxxG","core_SxxxxxR","core_GxxxxT","core_RxL","core_TI","core_RxxxxP","core_RL","core_LR","core_GxxxxxT","core_TG","core_GS","core_NxxT","core_KxxxxxC","core_TxxxxxN","core_LxxR","core_VR","core_PR","core_CxxxS","core_LxxxxxR","core_NxG","core_GxxxxxR","core_LxxxxR","core_TxxxxxS","core_GxxxxK","core_CxxxN","core_NxxxT","core_CxxxH","core_DxxxxxT","core_CxxxK","core_GxxxxxN","core_LxxxR","core_FxxxC","core_RxxV","core_NS","core_CR","core_MW","core_RxxxxL","core_SxxxxxT","core_GxN","core_TK","core_RxxL","core_W","core_AxP","core_LxR","core_DxxxxC","core_GxxxxN","core_KxxxxT","core_CxxxxS","core_KxT","core_DxxxS","core_TS","core_NxxxxxC","core_FxC","core_QxxxxC","core_RxV","core_RxxxxV","core_GxxT","core_RxxxL","core_VxxxxR","core_H","core_KxxxT","core_KxxT","core_CxD","core_VxR","core_GxxxS","core_LxxxxC","core_NxxxxT","core_PxxxxxC","core_VxxR","core_IxC","core_RxxxV","core_CxxxxxK","core_PA","core_TxK","core_TxE","core_GT","core_VxxxR","core_AC","core_A","core_DT","core_HR","core_CxG","core_FT","core_NxT","core_TxxxxxK","core_ExxxxT","core_DxT","core_LL","core_TxxR","core_NG","core_DxE","core_TxI","core_TxxxI","core_RxxT","core_RC","core_TxxxD","core_TxxxxK","core_IT","core_VT","core_YT","core_SP","core_RxxxxxL","core_PxxxxxP","core_RxxxxxV","core_SxxxxK","core_NxxxxP","core_RV","core_CxxK","core_SD","core_DC","core_PxxxP","core_LxxxL","core_IC","core_LP","core_GxxxN","core_AxL","core_VxxxxxC","core_CP","core_CxxR","core_PxxxxA","core_PxS","core_TxxxxxR","core_SxxxxN","core_AxxT","core_IxxxT","core_DxxT","core_TxA","core_PxL","core_CxxxR","core_HxR","core_VxxxxxR","core_AxxxxxN","core_DxxS","core_TxxxS","core_NxS","core_PxxA","core_RxxxxT","core_QR","core_PxxxxxA","core_QxxxxxC","core_LxxxxxL","core_AxxxxA","core_CI","core_PxxP","core_NxK","core_RxxxH","core_CxxxxR","core_DxxxxxD","core_D","core_AxxxA","core_TxxxR","core_PxC","core_TxxI","core_AxxxP","core_RxxH","core_SxxxD","core_KxxxxxG","core_PxxxxP","core_DD","core_AxxxxxP","core_LxxL","core_YxC","core_AxxA","core_AG","core_PS","core_LxxxxL","core_TxxxK","core_RxW","core_RxxxW","core_GxxxxxS","core_IxxT","core_LxxxxxA","core_ExxC","core_NxxxG","core_RT","core_TxxxxS","core_RxxxxxH","core_DF","core_HxxxxxK","core_WR","core_F","core_DxxxG","core_RxxxQ","core_SxxxL","core_RxxxT","core_RxxxxxQ","core_FxR","core_LxxxF","core_CxxxxxQ","core_VxxT","core_AxxxxxL","core_RxE","core_AxxxxP","core_RW","core_GP","core_WxxxxxR","core_QxxC","core_RxH","core_WxxxxR","core_YxxC","core_QxR","core_AP","core_CxxxxxR","core_CD","core_PxA","core_RxxxxH","core_PxxxxG","core_FxL","core_CA","core_HxxxxxR","core_PxxxxxV","core_FxxxR","core_V","core_VP","core_TxxxE","core_HxxxC","core_WxxR","core_AxxxxxA","core_RQ","core_RxxW","core_TxxxA","core_FxxL","core_VxC","core_AV","core_KxxxxxS","core_IxxxS","core_LxL","core_TxY","core_AxxxxL","core_DxxxT","core_VxK","core_DxxxK","core_CxxV","core_GxS","core_SxxxxxL","core_GxxxxP","core_CxxY","core_PxG","core_NxxxxxS","core_KxxD","core_CxxxxI","core_TxxxxR","core_QxxxxxR","core_ExxxxxR","core_ExxxxxC","core_IxxC","core_CxR","core_TxxV","core_WxR","core_SxxxN","core_M","core_LxxxxA","core_IxR","core_FxxxxxR","core_AL","core_GxxxxL","core_ExxT","core_DxxxxxR","core_TxxxxD","core_ED","core_LxG","core_CxE","core_IR","core_SxxxxL","core_RxxxxF","core_GxP","core_HxxxR","core_GxxxL","core_SxxP","core_DxxG","core_KxxxC","core_LxP","core_HxxxxR","core_ExxR","core_GxxI","core_LxS","core_FxxxxR","core_PxxxS","core_FL","core_AxxxxC","core_PxxL","core_ExR","core_VxA","core_QxxR","core_RxxxxxT","core_IxxxxC","core_SxP","core_GxA","core_PxxxxxL","core_PxxxxL","core_LxV","core_KxxxxD","core_DS","core_RxxxxxC","core_MxxT","core_GW","core_TQ","core_LS","core_AxxxL","core_LxxxxI","core_VL","core_HxxA","core_ER","core_DxxxR","core_LxxF","core_AxxxxxT","core_DP","core_VxxxxL","core_PxxxA","core_RxxxE","core_FP","core_KC","core_PH","core_TxxxG","core_LxxxxS","core_CE","core_FxxxxxL","core_VC","core_RxxxxD","core_VxG","core_HxxxxxP","core_HxxxP","core_IxxxxT","core_RxC","core_WxxxR","core_AxH","core_FF","core_SxxK","core_LxxxxP","core_RI","core_QxxxxR","core_LxH","core_LH","core_YC","core_FI","core_KxxxxG","core_GxV","core_RxxxxxE","core_QxxxxxL","core_RxxxxQ","core_RD","core_VxxxL","core_AxxxH","core_RxxxI","core_NxR","core_TE","core_LxxG","core_AxV","core_GxxxxS","core_FxxxxxT","core_PxP","core_NxxS","core_LxxxxxF","core_LV","core_RxxF","core_GD","core_GxxxP","core_SxxxxxN","core_ExxxxR","core_RxK","core_RH","core_WxA","core_ST","core_DR","core_TxxxxE","core_RxxxxxW","core_PxV","core_SxxxxP","core_WP","core_RxxxxxF","core_PxxxG","core_CxxxxxG","core_AxxL","core_RxxxxC","core_PT","core_NxxxP","core_TH","core_AH","core_QT","core_SxxxW","core_VxL","core_GxxxxxL","core_YR","core_FxxxxC","core_HxxxS","core_HxxL","core_VxxxxxV","core_GxxP","core_LxxxxxQ","core_CxP","core_DxR","core_KxxxE","core_HxxC","core_FxxxL","core_IF","core_QxxxxL","core_AxC","core_RxD","core_AxN","core_VxxxxxL","core_VxxG","core_FxxR","core_SxL","core_HxxxL","core_CxxxxxF","core_GH","core_WxxxxA","core_FxxxF","core_ExxxN","core_SxK","core_HxxxxA","core_RxxxxxD","core_GxxxxA","core_WA","core_TxR","core_PxxxL","core_KxxR","core_RxxxC","core_GE","core_DxxxxN","core_PxxxxW","core_WxxxxxA","core_KxxxxxL","core_RxQ","core_GL","core_LxxP","core_PD","core_TxxxxI","core_FxxC","core_TxxxQ","core_FxxS","core_KxxA","core_LxF","core_GxxL","core_RxT","core_AxxP","core_KxxxxC","core_AxxG","core_CxxxI","core_LxxV","core_RxF","core_GxxV","core_DxN","core_WxxxP","core_SxG","core_VxxxxP","core_CW","core_SxxxxA","core_MxxxxR","core_HxxS","core_GxxxxxP","core_ExxxD","core_HA","core_LxxxxH","core_FxxxxxF","core_NxxxxQ","core_SxxT","core_CxxxxA","core_FxxxxI","core_QxxxR","core_AxT","core_RxxxxW","core_KxxxxxR","core_SxW","core_SxxxA","core_YxxL","core_NxxV","core_KxS","core_IxxxxxL","core_RxxQ","core_ExT","core_PW","core_HxxxA","core_PL","core_KxR","core_GN","core_IxxxR","core_WxP","core_FxI","core_RxxxM","core_FR","core_PV","core_IxS","core_ExxxxxT","core_KxxxxR","core_RxxxxE","core_PF","core_HL","core_FxxxxL","core_RxxxF","core_SxS","core_WxxS","core_PxxxxxD","core_AxxxN","core_VxxxP","core_KxxxxE","core_CxV","core_MxR","core_HxxG","core_RK","core_RxY","core_LF","core_RM","core_VxxV","core_GxxxA","core_RN","core_FxxxI","core_RxxI","core_VxxL","core_FxF","core_AxxxxS","core_GxL","core_FxP","core_IxxxxR","core_IxT","core_DxxN","core_AxA","core_KxxxL","core_LxxxH","core_SxxxxxW","core_YL","core_RxxY","core_PP","core_KxG","core_RxxD","core_HxxxxL","core_PxxxT","core_TxxxF","core_AxxH","core_AW","core_IxxxxF","core_IxxxxxR","core_SxxL","core_LxxH","core_WxxxA","core_MxxR","core_YxR","core_QxxxL","core_RxxxxxN","core_RxxxD","core_SS","core_RxN","core_SxxD","core_LxxS","core_VxxxxxA","core_LxN","core_TxxD","core_LxxxP","core_DxxxxY","core_RE","core_ET","core_DxxxD","core_VxxxxV","core_VI","core_LxxxxxI","core_CxxQ","core_PN","core_MC","core_HC","core_SA","core_LxxxS","core_AxxxxxW","core_LxxxxV","core_NxxQ","core_VxxxG","core_VW","core_VxxxA","core_NxD","core_IxxxxxF","core_NxxxxA","core_LA","core_VxxxxxS","core_RxxM","core_YxxxxT","core_PxxxxV","core_VxxH","core_ExxxxxN","core_DxxxxxE","core_PxxxxxG","core_FxA","core_VxV","core_SxD","core_FxxxxxK","core_SxxxxW","core_PxT","core_SK","core_RF","core_RxxxxxM","core_AxxxxW","core_NxN","core_IxxF","core_TxS","core_PxW","core_TxxA","core_DH","core_TxxK","core_GxxxxxD","core_HxA","core_LxM","core_FxxxW","core_IxL","core_NxxG","core_KxL","core_SN","core_CxxxxD","core_AxxxxM","core_DG","core_MxxxR","core_FxxxxG","core_RxxxxxY","core_DxG","core_HxG","core_LxxxxG","core_LK","core_LxxxxxY","core_MV","core_LxY","core_TxxF","core_VxP","core_NxxxxN","core_GxxxK","core_IxxxxxT","core_PxxS","core_LxxxxxH","core_YxxxxxR","core_HP","core_AM","core_MxxP","core_HxxxxxL","core_NxxxxxD","core_VG","core_GxxxxH","core_SxA","core_LxxxA","core_NxxxxE","core_SxxxxxA","core_RxxxxM","core_NR","core_TxxxxxE","core_ExxxxxH","core_WxxP","core_GxxA","core_SxM","core_RxxxxxI","core_VF","core_SG","core_FQ","core_FxxN","core_DxxxN","core_CV","core_IxxA","core_QC","core_NxxR","core_TxxxxQ","core_FA","core_FxxxV","core_LxxxxxW","core_YxxxL","core_DxxxxL","core_RxxxxK","core_QxxF","core_DxxxxR","core_FxxxxxI","core_VxxxxxP","core_AxW","core_FxxxS","core_WxxxxxL","core_DxxxxxN","core_AT","core_LxQ","core_LxxxxxV","core_LxxxxxM","core_QxxxxxH","core_LxxK"]
#        
#        
        for seq in core_seqs:
            regex = re.compile(seq[5:].replace("x", "."))
            scoring_csv_columns.append(len(regex.findall(self.core)))
            scoring_column_names.append(seq)
        core_pi = ProteinAnalysis(self.core).isoelectric_point()
        scoring_csv_columns.append(core_pi)
        scoring_column_names.append("core_pI")
        
        if core_pi < 9:
            self.score += 4
        if self.core.count("C") >= 2:
            self.score += 4
            
        regex = re.compile("KL.L.K")
        if self.leader[-2:] == "GG" and len(regex.findall(self.leader)) > 0:
            self.score += 2
        
        scoring_csv_columns = [x for _, x in sorted(zip(scoring_column_names, scoring_csv_columns))]
        self.csv_columns += [self.score] + scoring_csv_columns
        
