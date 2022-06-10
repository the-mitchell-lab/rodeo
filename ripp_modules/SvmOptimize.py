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
   SVM optimization script
   
   Required input:
     a training file CSV
     
   Output:
     a CSV list of identifiers and classifications
     
   Note that all options must be hard-coded (the script does not take command-line arguments).
   
   RECOMMENDATION:
     Input CSV should ideally have its primary key as column 0, its classification as column 1, and feature as columns [2,...,end]
     For the fitting CSV, there will be no classification; leave it blank if you want (it'll be ignored upon import)
'''

import numpy as np
import csv
import sys

from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from sklearn import preprocessing
import scipy.stats as stats
from sklearn.utils.fixes import loguniform

from sklearn.model_selection import GridSearchCV, RandomizedSearchCV
from sklearn.model_selection import cross_val_score
from sklearn import metrics

from SvmClassify import SVMRunner

# CONFIGURATION OPTIONS
''' change these as desired '''


primary_key_column = 0;            # the column of the CSV that contains the primary key (identifier) for each record
classification_column = 1;         # the column of the CSV that contains the classification for each record
csv_has_header = True;             # set to true if CSVs have header rows to ignore upon import

# values for optimization
kernel_option = 'rbf'       # the options below are for rbf; feel free to modify for poly, linear, etc
C_max = 10
C_min = -1
C_base = 5
C_steps = 11
C_options = np.logspace(C_min,C_max,num=C_steps,base=C_base,dtype=float)
gamma_max = -2
gamma_min = -10
gamma_base = 10
gamma_steps = 9
gamma_options = np.logspace(gamma_min,gamma_max,num=gamma_steps,base=gamma_base,dtype=float)
class_weight_option = 'balanced'
folds_validation = [10]


def report(results, n_top=3):
    for i in range(1, n_top + 1):
        candidates = np.flatnonzero(results["rank_test_score"] == i)
        for candidate in candidates:
            print("Model with rank: {0}".format(i))
            print(
                "Mean validation score: {0:.3f} (std: {1:.3f})".format(
                    results["mean_test_score"][candidate],
                    results["std_test_score"][candidate],
                )
            )
            print("Parameters: {0}".format(results["params"][candidate]))
            print("")
    
    
def main():
    input_training_file = sys.argv[1]         # the CSV containing the training set
    # parse data
    
    print("Importing training and fitting data ...")
    runner = SVMRunner("default")
    training_data_unrefined = runner.parse_CSV_to_dataset(input_training_file, 'training')
    
    primary_key_list = []

    training_data_just_features = []
    training_data_classifications = []
   
    for entry in training_data_unrefined:
        training_data_classifications.append(entry.pop(1))
        entry.pop(0)
        training_data_just_features.append(entry)  
    print("(done)")
    print("Initiating learning and fitting")
    
    # Scaling -- this ensures standardization of model and target data
    # training_data_refined = preprocessing.scale(training_data_just_features)
    scaler = preprocessing.StandardScaler()#.fit(training_data_just_features)
    training_data_refined = scaler.fit_transform(training_data_just_features)
    print("Data has been scaled" )
    test_results = []
    
    # C_option = gamma_option = 1
    # for fold in folds_validation:
        # for C_option in C_options:
            # for gamma_option in gamma_options:
                # # print(C_option, gamma_option)
                # clf = svm.SVC(kernel=kernel_option,class_weight=class_weight_option,C=C_option,gamma=gamma_option)
                # # clf = RandomForestClassifier(n_estimators=250, min_samples_split=2)
                # clf.fit(training_data_refined, training_data_classifications)
                # prec = cross_val_score(clf, training_data_refined, training_data_classifications, cv=fold, scoring='precision')
                # recd = cross_val_score(clf, training_data_refined, training_data_classifications, cv=fold, scoring='recall')
                # f1we = cross_val_score(clf, training_data_refined, training_data_classifications, cv=fold, scoring='f1')
                # scor = prec.mean() * recd.mean()
                # print('Using %d fold, %.4E C, %.4E gamma: %0.4f precision, %0.4f recall, %0.4f f1, %0.4f score' % (fold, C_option, gamma_option, prec.mean(), recd.mean(), f1we.mean(), scor))
                # test_results.append([kernel_option,fold,C_option,gamma_option,class_weight_option,prec.mean(),recd.mean(),f1we.mean(),scor])
    param_distribution = {'C': loguniform(1e-2, 1e15),
         'gamma': loguniform(1e-10, .1),
         'kernel': ['rbf'],
         'class_weight':['balanced', None]}
    clf = svm.SVC(kernel=kernel_option,class_weight=class_weight_option)
    random_search = RandomizedSearchCV(clf, param_distributions=param_distribution, cv=10, n_jobs=40, n_iter=100, scoring="f1")
    random_search.fit(training_data_refined, training_data_classifications)
    report(random_search.cv_results_)

    
    # best = max(test_results, key=lambda x: x[-2])
    # print("C: {}, G: {}, f1: {}".format(best[2], best[3], best[-2]))
    
    return

if __name__ == '__main__':
  main()
