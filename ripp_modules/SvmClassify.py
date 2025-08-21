# -*- coding: utf-8 -*-
#==============================================================================
# Copyright (C) 2017 Bryce L. Kille
# University of Illinois
# Department of Chemistry
#
# Copyright (C) 2015 Jonathan I. Tietz
# University of Illinois
# Department of Chemistry
#
# Copyright (C) 2015 Christopher J. Schwalen
# University of Illinois
# Department of Chemistry
#
# Copyright (C) 2015 Parth S. Patel
# University of Illinois
# Department of Chemistry
#
# Copyright (C) 2015 Douglas A. Mitchell
# University of Illinois
# Department of Chemistry
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
#
'''
   SVM classification script
   
   Required input:
     a training file CSV
     a set of model data CSV (if fitting)
     
   Output:
     a CSV list of identifiers and classifications
     
   Note that all options must be hard-coded (the script does not take command-line arguments).
   
   This script assumes you have already optimized kernel parameters. Use the included optimization script if not.
   
   RECOMMENDATION:
     Input CSV should ideally have its primary key as column 0, its classification as column 1, and feature as columns [2,...,end]
     For the fitting CSV, there will be no classification; leave it blank if you want (it'll be ignored upon import)
'''

import csv
from sklearn import svm
from sklearn import preprocessing
#from scikit-learn import svm
#from scikit-learn import preprocessing
import pathlib
#from rodeo_main import output_dir
FILE_DIR = pathlib.Path(__file__).parent.absolute()


    #def output_dir_retrieve():
        

class SVMRunner(object):
    
    def __init__(self, peptide_type, output_dir):
        self.peptide_type = peptide_type
        self.output_dir = output_dir
        # CONFIGURATION OPTIONS
        ''' change these as desired '''
        prefix = '{}/{}/svm'.format(FILE_DIR, peptide_type)
        
        self.input_training_file = '{}/{}_training_set.csv'.format(prefix, peptide_type)     # the CSV containing the training set
        self.input_fitting_file = '{}/{}/fitting_set.csv'.format(output_dir, peptide_type)          # the CSV containing the data to be fitted
        self.output_filename = '{}/{}/fitting_results.csv'.format(output_dir, peptide_type)      # output filename; this will be a CSV with the first column being the primary key and the second being the classification
        
        self.primary_key_column = 0;            # the column of the CSV that contains the primary key (identifier) for each record
        self.classification_column = 1;         # the column of the CSV that contains the classification for each record
        self.csv_has_header = True;             # set to true if CSVs have header rows to ignore upon import
        
        self.kernel_option = 'rbf'
        self.class_weight_option = 'balanced'
        if peptide_type == 'grasp':
            self.C_option = 9.172484262400534
            self.gamma_option = 0.00022271232130782889
        elif peptide_type == 'lasso':
            self.C_option = 6.8986
            self.gamma_option = 1E-4
        elif peptide_type == 'thio':        
            self.C_option = 283117
            self.gamma_option = 1E-9
        elif peptide_type == 'lanthi':
            self.C_option = 48205.77
            self.gamma_option = 1E-9
        elif peptide_type == 'lanthi1':
            self.C_option = 283117
            self.gamma_option = 1E-7
        elif peptide_type == 'lanthi2':
            self.C_option = 8208
            self.gamma_option = 1E-6
        elif peptide_type == 'lanthi3':
            self.C_option = 283117
            self.gamma_option = 1E-7
        elif peptide_type == 'lanthi4':
            self.C_option = 283117
            self.gamma_option = 1E-7
        elif peptide_type == 'sacti':
            self.C_option = 17651
            self.gamma_option = 2.1E-7
        elif peptide_type == 'linar':
            self.C_option = 1.17462
            self.gamma_option = 0.001
        elif peptide_type == 'boro':
            self.C_option = 1.1746189430880192
            self.gamma_option = 0.0001
        elif peptide_type == 'cyclo':
            self.C_option = 77.18912785617864
            self.gamma_option = 0.00018171423466314444
        else:
            # default kernel values
            self.C_option = 2.8E10
            self.gamma_option = 1E-10

    def parse_CSV_to_dataset(self, csv_filename, dataset_type):
        '''Parse an input CSV into a data set
        
           Inputs:
                csv_filename            name of CSV file to be parsed
                dataset_type            either 'training' or 'fitting'
        '''
        dataset = []
        with open(csv_filename, 'r') as csvfile:
          csv_read = csv.reader(csvfile, delimiter=',', quotechar='"')
          if self.csv_has_header == True:
            next(csv_read,None)
            
          if self.classification_column < self.primary_key_column:
            pk_first = 0
          else:
            pk_first = 1
            
          for row in csv_read:
            temp_entry = []
            temp_entry.append(row.pop(self.primary_key_column))
            temp_entry.append(row.pop(self.classification_column - pk_first))
            for c in row:
              temp_entry.append(float(c))
            # remove all unclassified features if training
            if dataset_type == 'training':
              if int(temp_entry[1]) == 1 or int(temp_entry[1]) == 0:
                temp_entry[1] = int(temp_entry[1])
                dataset.append(temp_entry)
            if dataset_type == 'fitting':
              dataset.append(temp_entry)
        return dataset
        
    def write_to_csv(self, list_of_primary_keys, list_of_classifications, output_file):
        with open(self.output_filename, 'w') as csvfile:
          csv_write = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
          for i in range(len(list_of_primary_keys)):
            temp_row = [list_of_primary_keys[i],list_of_classifications[i]]
            csv_write.writerow(temp_row)
        return
        
    def run_svm(self):
    
        # parse data
        
        training_data_unrefined = self.parse_CSV_to_dataset(self.input_training_file, 'training')
        fitting_data_unrefined = self.parse_CSV_to_dataset(self.input_fitting_file, 'fitting')
        
        primary_key_list = []
        fitting_data_just_features = []
        training_data_just_features = []
        training_data_classifications = []
       
        for entry in training_data_unrefined:
            training_data_classifications.append(entry.pop(1))
            entry.pop(0)
            training_data_just_features.append(entry)  
        for entry in fitting_data_unrefined:
            primary_key_list.append(entry.pop(0))
            entry.pop(0)
            fitting_data_just_features.append(entry)
        
        # Scaling -- this ensures standardization of model and target data
        # training_data_refined = preprocessing.scale(training_data_just_features)
        scaler = preprocessing.StandardScaler().fit(training_data_just_features)
        training_data_refined = scaler.transform(training_data_just_features)
        fitting_data_refined = scaler.transform(fitting_data_just_features)
    
        # This creates the classifier and learns using the model data
        clf = svm.SVC(kernel=self.kernel_option,class_weight=self.class_weight_option,C=self.C_option,gamma=self.gamma_option)
        clf.fit(training_data_refined, training_data_classifications)
    
        # This classifies the input data as a list
        classification_list = clf.predict(fitting_data_refined)
        
        # Output results to file
        self.write_to_csv(primary_key_list, classification_list, self.output_filename)
        
        return

