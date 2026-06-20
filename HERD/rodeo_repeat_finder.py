import argparse, sys, os, re, csv, glob, time
import numpy as np
import pandas as pd
import math
def __main__():

# Document starts with argument parsing
    parser = argparse.ArgumentParser("RODEO Co-occurrence calculator")
    parser.add_argument('folder_input', type=str)
    parser.add_argument('OFO', type=str)
    parser.add_argument('-n', '--number', type=int, default=3)
    parser.add_argument('-p', '--hmm_number', type=int, default=50)
    parser.add_argument('-s', '--single', type=int, default=0)
    parser.add_argument('-f', '--locus_filter', nargs='*', default=[])
    parser.add_argument('-print', '--print', action='store_true', default=False)
    args = parser.parse_args()

    st = time.time()

    if not os.path.exists(args.OFO):
        os.mkdir(args.OFO)
    co_occur_files = []
    results_files = []
    if args.single:
        co_occur_files.append(args.folder_input+"main_co_occur.csv")
        results_files.append(args.folder_input+"main_results.csv")
    else:
        co_occur_files = glob.glob(args.folder_input+"main_co_occur.csv")
        results_files = glob.glob(args.folder_input+"main_results.csv")

    pre_filter_list = {}
    filter_list = []

    if len(args.locus_filter) > 0:
        for file in co_occur_files:
            curr_file = open(file.strip(), 'r')
            reader = csv.reader(open(file,'r'))
            lines = list(reader)
            i = 1
            meta = False
            if "Locus" in lines[0]:
                meta = True
            while i < len(lines):
                attr = lines[i]
                for keyword in args.locus_filter:
                    if meta and keyword in ",".join(attr[:args.number*4+8]):
                        if (attr[0] + "_" + attr[1] + "_" + attr[3]) not in pre_filter_list.keys():
                            pre_filter_list[attr[0] + "_" + attr[1] + "_" + attr[3]] = []
                        pre_filter_list[attr[0] + "_" + attr[1] + "_" + attr[3]].append(keyword)
                    elif keyword in ",".join(attr[:args.number*4+7]):
                        if (attr[0] + "_---_" + attr[2]) not in pre_filter_list.keys():
                            pre_filter_list[attr[0] + "_---_" + attr[2]] = []
                        pre_filter_list[attr[0] + "_---_" + attr[2]].append(keyword)
                i = i+1

        for hit in pre_filter_list.keys():
            if all(x in args.locus_filter for x in pre_filter_list[hit]):
                filter_list.append(hit)

    log_output = open(args.OFO + '/log.txt', 'w')
    for arg in vars(args):
        log_output.write(arg + ":\t" + str(getattr(args, arg)) + "\n")

    file_output = open(args.OFO + '/repeats.txt', 'w')

    for file in results_files:
        seqs_to_output = {}
        curr_file = open(file.strip(), 'r')
        reader = csv.reader(open(file, 'r'))
        lines = list(reader)
        i = 1
        meta = False
        if "Locus" in lines[0]:
            meta = True
        while i < len(lines):
            attr = lines[i]
            derandom_tag = attr[0][:-11]
            if meta:
                curr_query = derandom_tag + "_" + attr[1] + "_" + attr[3]
            else:
                curr_query = attr[0] + "_---_" + attr[2]
            if (curr_query not in filter_list and len(args.locus_filter) > 0) or "X" in attr[-1]:
                i=i+1
                continue
            detected_repeat, repeat_occurrences = seq_detect_repeats(attr[-1])
            if len(repeat_occurrences) >= 4 and ("S" in detected_repeat or "T" in detected_repeat):# and re.search("[Y|F|W|V]..[P|L]", str(attr[-1])):
                if curr_query not in seqs_to_output.keys():
                    seqs_to_output[curr_query] = [detected_repeat, repeat_occurrences, curr_query, attr[-1]]
                elif len(repeat_occurrences) >= len(seqs_to_output[curr_query][1]) and len(attr[1]) < len(seqs_to_output[curr_query][3]):
                    seqs_to_output[curr_query] = [detected_repeat, repeat_occurrences, curr_query, attr[-1]]
            i = i+1
            #if repeat_occurrences and sum(1 for residue in detected_repeat if residue in {"T", "K", "S", "E", "D"})/len(detected_repeat) > 0.5:
        for seq in seqs_to_output.keys():
            print(seqs_to_output[seq][0], seqs_to_output[seq][1], seqs_to_output[seq][2], seqs_to_output[seq][3])
            file_output.write(">" + seqs_to_output[seq][2] + "_" + seqs_to_output[seq][0] + "_" + str(len(seqs_to_output[seq][1])) + "\n" + str(seqs_to_output[seq][3]) + "\n")
            

def seq_detect_repeats(seq, kmin=2, kmax=10, threshold=0.5):
    k = kmin
    detected_repeat = ""
    repeat_occurrences = []
    while k <= kmax:
        kmer_dict = index_seq(seq, k)
        repeat_word, repeat_indices = find_top_kmer(kmer_dict)
        if repeat_indices:
            total_repeat_len = len(repeat_word)*len(repeat_indices)
            coverage = total_repeat_len/float(max(repeat_indices)+k-1-min(repeat_indices))
            if total_repeat_len > max(12,len(repeat_occurrences)*len(detected_repeat)) and coverage >= threshold:
                detected_repeat = repeat_word
                repeat_occurrences = repeat_indices
        k = k+1
    return detected_repeat,repeat_occurrences
    

def index_seq(seq,k):
    kmer_dict = {}
    for i in range(len(seq)-k):
        word = str(seq[i:i+k])
        if word in kmer_dict:
            if i-k >= max(kmer_dict[word]):
                kmer_dict[word].append(i)
        else:
            kmer_dict[word] = [i]
    return kmer_dict

def find_top_kmer(kmer_dict):
    
    top_kmer_word = ""
    top_kmer_hits = []
    n_hits = 1
    for key in kmer_dict:
        if len(kmer_dict[key])>n_hits:
            top_kmer_word = key
            top_kmer_hits = kmer_dict[key]
            n_hits = len(kmer_dict[key]) * math.sqrt(len(key))

    return top_kmer_word,top_kmer_hits

if __name__=="__main__":
    __main__()