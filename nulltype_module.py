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
#==============================================================================

import csv

def main_write_headers(output_folder, meta=False):
    if meta:
        headers = ['Query', 'Locus', 'Genus/Species', 'Nucleotide_acc', 'start', 'end', 'dir', 'AA_seq']
    else:
        headers = ['Query', 'Genus/Species', 'Nucleotide_acc', 'start', 'end', 'dir', 'AA_seq']
    features_writer = csv.writer(open(output_folder + "/main_results.csv", 'w'))
    features_writer.writerow(headers)
    
def main_write_row(output_folder, row):
    features_writer = csv.writer(open(output_folder + "/main_results.csv", 'a'))
    features_writer.writerow(row)
    
def co_occur_write_headers(output_folder, meta):
    if meta:
        headers = ['Query', 'Locus', 'Genus/Species', 'Nucleotide_acc', 'Protein_acc', 'start', 'end', 'dir',\
               'PfamID1', 'Name1', 'Description1', 'E-value1',
               'PfamID2', 'Name2', 'Description2', 'E-value2',
               'PfamID3', 'Name3', 'Description3', 'E-value3']
    else:
        headers = ['Query', 'Genus/Species', 'Nucleotide_acc', 'Protein_acc', 'start', 'end', 'dir',\
               'PfamID1', 'Name1', 'Description1', 'E-value1',
               'PfamID2', 'Name2', 'Description2', 'E-value2',
               'PfamID3', 'Name3', 'Description3', 'E-value3']
    features_writer = csv.writer(open(output_folder + "/main_co_occur.csv", 'w'))
    features_writer.writerow(headers)
    
def co_occur_write_row(output_folder, row):
    features_writer = csv.writer(open(output_folder +  "/main_co_occur.csv", 'a'))
    features_writer.writerow(row)

def prod_write_headers(output_folder):
    headers = ['ID', 'Query', 'Genus/Species', 'Nucleotide_acc', 'start', 'end', 'dir', 'Prod_Score', 'CodPot', 'StrtSc', 'Codon', 'RBSMot', 'Spacer', 'RBSScr', 'UpsScr', 'TypeScr', 'AA_seq']
    features_writer = csv.writer(open(output_folder + "/prodigal/prod_results.csv", 'w'))
    features_writer.writerow(headers)
    
def prod_write_row(output_folder, row):
    features_writer = csv.writer(open(output_folder + "/prodigal/prod_results.csv", 'a'))
    features_writer.writerow(row)