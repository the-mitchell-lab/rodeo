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

#DELIMITERS
BEGIN_RIPP = ">"
END_RIPP = ">>"
BEGIN_COLORS = "#BEGIN_COLORS"
END_COLORS = "#END_COLORS"
BEGIN_VARIABLES = "#BEGIN_VARIABLES"
END_VARIABLES = "#END_VARIABLES"



def parse_config_file(filename):
    conf_dict = {}
    with open(filename) as config_handle:
        ripp_type = ""
        in_color = False
        in_variable = False
        in_ripp = False
        for line in config_handle:
            line = line.strip()
            if line == "\n" or line == "":
                continue
            elif line[0] == BEGIN_RIPP and line != ">>":
                in_ripp = True
                ripp_type = line[1:].lower()
                conf_dict[ripp_type] = {}
                conf_dict[ripp_type]['pfam_colors'] = {}
                conf_dict[ripp_type]['variables'] = {}
            elif line == END_RIPP:
                in_ripp = False
            elif line == BEGIN_COLORS:
                in_color = True
            elif line == END_COLORS:
                in_color = False
            elif line == BEGIN_VARIABLES:
                in_variable = True
            elif line == END_VARIABLES:
                in_variable = False
            elif in_ripp:
                line = line.split(' ')
                key = line[0]
                val = line[1]
                if line[0] == 'int':
                    key = line[1]
                    val = int(line[2])
                elif line[0] == 'bool':
                    key = line[1]
                    if line[2].lower() == "yes" or line[2].lower() == "true":
                        val = True
                    else:
                        val = False
                
                if in_color:
                    conf_dict[ripp_type]['pfam_colors'][key.upper()] = val.lower()
                elif in_variable:
                    conf_dict[ripp_type]['variables'][key.lower()] = val
    return conf_dict

def merge_confs(conf_list):
    master_conf_dict = conf_list[0]
    for conf in conf_list[1:]:
        for ripp in conf.keys():
            current_dict = conf[ripp]
            for entry_type in current_dict.keys():
                current_dict_nested = current_dict[entry_type]
                for key in current_dict_nested:
                    master_conf_dict[ripp][entry_type][key] = current_dict_nested[key]
    return master_conf_dict

def merge_conf_and_arg(conf, args):    
    if args.precursor_max is not None:
        conf['general']['variables']['precursor_max'] = args.precursor_max
        for peptide_type in conf.keys():
                conf[peptide_type]['variables']['precursor_max'] = args.precursor_max
    
    if args.precursor_min is not None:
        conf['general']['variables']['precursor_min'] = args.precursor_min
        for peptide_type in conf.keys():
            conf[peptide_type]['variables']['precursor_min'] = args.precursor_min
    
    if args.overlap is not None:
        conf['general']['variables']['overlap'] = args.overlap
        
    if args.fetch_type is not None:
        conf['general']['variables']['fetch_type'] = args.fetch_type
    
    if args.fetch_n is not None:
        conf['general']['variables']['fetch_n'] = args.fetch_n
        
    if args.fetch_distance is not None:
        conf['general']['variables']['fetch_distance'] = args.fetch_distance
    
    if args.evaluate_all is not None:
        conf['general']['variables']['evaluate_all'] = args.evaluate_all
        
    if args.exhaustive is not None:
        for peptide_type in conf.keys():
            conf[peptide_type]['variables']['exhaustive'] = args.exhaustive
    
    if args.print_precursors is not None:
        for peptide_type in conf.keys():
            conf[peptide_type]['variables']['print_precursors'] = args.print_precursors
    return conf
        
