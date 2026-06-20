import argparse, sys, os, re, csv, glob, time
import networkx as nx
import xml.parsers.expat
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import hdbscan
import numpy as np
import pandas as pd
import math
from decimal import Decimal

def __main__():

# Document starts with argument parsing
    parser = argparse.ArgumentParser("RODEO Co-occurrence calculator")
    parser.add_argument('file_input', type=str)
    parser.add_argument('OFO', type=str)
    parser.add_argument('-n', '--number', type=int, default=3)
    parser.add_argument('-p', '--hmm_number', type=int, default=1000)
    parser.add_argument('-s', '--single', type=int, default=0)
    parser.add_argument('-f', '--locus_filter', nargs='*', default=[])
    parser.add_argument('-e', '--exclude', nargs='*', default=[])
    parser.add_argument('-print', '--print', action='store_true', default=False)
    args = parser.parse_args()

    st = time.time()

    if not os.path.exists(args.OFO):
        os.mkdir(args.OFO)
    files=[]
    if args.single:
        files.append(args.file_input)
    else:
        files = glob.glob(args.file_input)

    pre_filter_list = {}
    pre_exclude_list = {}
    filter_list = []

    if len(args.locus_filter) > 0:
        for file in files:
            curr_file = open(file.strip(), 'r')
            reader = csv.reader(open(file,'r'))
            lines = list(reader)
            i = 1
            meta = False
            if "Locus" in lines[0]:
                meta = True
            while i < len(lines):
                attr = lines[i]
                if meta:
                    s = 8
                    temp_name = attr[0] + "_" + attr[1] + "_" + attr[3]
                else:
                    s = 7
                    temp_name = attr[0] + "_---_" + attr[2]
                for keyword in args.locus_filter:
                    for t in range(s, len(attr)-2, 4):
                        if keyword in str(attr[t]) and Decimal(attr[t+3]) < 1E-30:
                            if (temp_name) not in pre_filter_list.keys():
                                pre_filter_list[temp_name] = []
                            pre_filter_list[temp_name].append(keyword)
                for keyword in args.exclude:
                    for t in range(s, len(attr)-2, 4):
                        if keyword in str(attr[t]) and Decimal(attr[t+3]) < 1E-30:
                            if (temp_name) not in pre_exclude_list.keys():
                                pre_exclude_list[temp_name] = []
                            pre_exclude_list[temp_name].append(keyword)
                i = i+1

        for hit in set(pre_filter_list.keys()).difference(set(pre_exclude_list.keys())):
            if all(x in pre_filter_list[hit] for x in args.locus_filter):
                filter_list.append(hit)

    log_output = open(args.OFO + '/log.txt', 'w')
    for arg in vars(args):
        log_output.write(arg + ":\t" + str(getattr(args, arg)) + "\n")

    allpfamhits = {}
    alldata = {}
    curr_BGC_combos = set()
    self_hit_matrix = {}
    bgc_otherhit_matrix = {}
    bgc_allhit_matrix = {}
    bgc_hit_counter = {}
    total_bgc_count = 0
    records = {}
    #metabait = ["PF02624__YcaO", "PF05147__LANC_like", "PF13471__Transglut_core3", "PF05114__DUF692", "PF04738__Lant_dehydr_N", "TIGR03793__TIGR03793",  "Graspetide_synthetase__Graspetide_synthetase", "PF00590__TP_methylase", "TIGR03798__leader_Nif11", "PF05402__PqqD"]
    
    for file in files:
        curr_file = open(file.strip(), 'r')
        reader = csv.reader(open(file,'r'))
        lines = list(reader)
        i = 1
        curr_query = ""
        curr_protein = ""
        curr_BGC_combos = set()
        curr_protein_combos = set()
        meta = False
        if "Locus" in lines[0]:
            meta = True
        while i < len(lines):
            attr = lines[i]
            if meta:
                curr_query = attr[0] + "_" + attr[1] + "_" + attr[3]
            else:
                curr_query = attr[0] + "_---_" + attr[2]
            if len(args.locus_filter) > 0:
                if curr_query not in filter_list:
                    i = i+1
                    continue
            if curr_query not in records.keys():
                records[curr_query] = My_Record(curr_query)
                if meta:
                    records[curr_query].cluster_accession = attr[3]
                    records[curr_query].cluster_genus_species = attr[2]
                    records[curr_query].locus = attr[1]
                else:
                    records[curr_query].cluster_accession = attr[2]
                    records[curr_query].cluster_genus_species = attr[1]
            if curr_query not in alldata:   # If the current query is not in the dataset, add it to the list and reset the last BGC, since a new BGC has started
                alldata[curr_query] = []
                total_bgc_count = total_bgc_count + 1
                for hit in list(curr_BGC_combos):
                    bgc_hit_counter[hit] = bgc_hit_counter[hit] + 1
                curr_BGC_combos = set()
            curr_protein_combos = set() # resets the temp protein data for the start of the next BGC
            if meta:
                curr_protein = attr[4] + "_" + attr[1] + "_" + attr[3]      # change to 3 if not meta
                curr_cds = Sub_Seq(int(attr[5]), int(attr[6]), accession_id=attr[4], meta=meta)
            else:
                curr_protein = attr[3] + "_---_" + attr[2]
                curr_cds = Sub_Seq(int(attr[4]), int(attr[5]), accession_id=attr[3])
            pfam_hits = []              # may be superfluous, but I didn't want to deal with extra hashing stuff for data handling later
            j = lines[0].index('PfamID1')
            old_BGC_combos = list(curr_BGC_combos)
            while j < len(attr)-2 and j < (4*args.number + 6):
                pfam_descr = [attr[j], attr[j+2], attr[j+3], attr[j+1]]
                curr_cds.pfam_descr_list.append(pfam_descr)
                j = j+4
                if j >= (4*args.number + 6):
                    continue
                duplicate_in_BGC = 0
                duplicate_in_protein = 0
                hit1 = (str(attr[j-4]) + "__" + str(attr[j-3])).strip()
                pfam_hits.append(hit1)         # add the new HMM hit to the list for the current protein (maintains flexibility, since array and not set)
                if hit1 in curr_protein_combos:
                    duplicate_in_protein = 1
                curr_protein_combos.add(hit1)
                if hit1 in curr_BGC_combos:
                    duplicate_in_BGC = 1
                curr_BGC_combos.add(hit1)
                if hit1 not in allpfamhits.keys():    # checking that all the 2D matrices are properly initialized basically
                    allpfamhits[hit1] = [attr[j-4], attr[j-3], attr[j-2], 0]
                if hit1 not in self_hit_matrix.keys():
                    self_hit_matrix[hit1] = {}
                    self_hit_matrix[hit1][hit1] = 0
                if hit1 not in bgc_hit_counter.keys():
                    bgc_hit_counter[hit1] = 0
                if hit1 not in bgc_otherhit_matrix.keys():
                    bgc_otherhit_matrix[hit1] = {}
                if hit1 not in bgc_allhit_matrix.keys():
                    bgc_allhit_matrix[hit1] = {}
                allpfamhits[hit1][3] = allpfamhits[hit1][3] + 1   # increment the counter for total occurrences
                       # reference the current protein HMM hit against the previous HMMs in the BGC
                for hit2 in list(old_BGC_combos):
                    if hit2 in bgc_otherhit_matrix[hit1].keys():
                        bgc_otherhit_matrix[hit1][hit2] = bgc_otherhit_matrix[hit1][hit2] + 1
                    else:
                        bgc_otherhit_matrix[hit1][hit2] = 1
                    if hit1 in bgc_otherhit_matrix[hit2].keys():
                        bgc_otherhit_matrix[hit2][hit1] = bgc_otherhit_matrix[hit2][hit1] + 1
                    else:
                        bgc_otherhit_matrix[hit2][hit1] = 1
                       # reference the current protein HMM list against itself for self_hit_matrix
                for hit2 in list(curr_protein_combos):
                    if hit2 in self_hit_matrix[hit1].keys():
                        self_hit_matrix[hit1][hit2] = self_hit_matrix[hit1][hit2] + 1
                    else:
                        self_hit_matrix[hit1][hit2] = 1
                    if hit1 is not hit2:
                        if hit1 in self_hit_matrix[hit2].keys():
                            self_hit_matrix[hit2][hit1] = self_hit_matrix[hit2][hit1] + 1
                        else:
                            self_hit_matrix[hit2][hit1] = 1
                       # reference the current protein HMM list against all HMMs in the BGC
                for hit2 in list(curr_BGC_combos):
                    if duplicate_in_BGC:
                        break
                    if hit2 in bgc_allhit_matrix[hit1].keys():
                        bgc_allhit_matrix[hit1][hit2] += 1
                    else:
                        bgc_allhit_matrix[hit1][hit2] = 1
                    if hit1 is not hit2:
                        if hit1 in bgc_allhit_matrix[hit2].keys():
                            bgc_allhit_matrix[hit2][hit1] += 1
                        else:
                            bgc_allhit_matrix[hit2][hit1] = 1

            curr_protein_data = {curr_protein:pfam_hits}
            alldata[curr_query].append(curr_protein_data)
            records[curr_query].CDSs.append(curr_cds)
            i = i+1

        for hit in list(curr_BGC_combos):
            bgc_hit_counter[hit] = bgc_hit_counter[hit] + 1

    if args.print:
        log_output.write("\n\nProtein queries")
        for query in alldata.keys():
            log_output.write("\n" + "\t".join(query.split("_")[-2:]))

    main_html = open(args.OFO+"/main_html.html", 'w')
    write_header(main_html)
    for record in records.keys():
        write_record(main_html, records[record], args.locus_filter)

    print_list = []
    for entry in allpfamhits.keys():                # convert data into printable format and sortable format
        print_list.append([allpfamhits[entry][0], allpfamhits[entry][1], allpfamhits[entry][2], allpfamhits[entry][3], entry])
    print_list.sort(key=lambda x: x[3], reverse=True)       # sort the HMM list by occurrences (greatest to lowest)
    OFO = open(args.OFO+"/BGC_co_occurrence.tsv", 'w')      
    OFO.write("Total Loci\t" + str(total_bgc_count)+ "\n\nHMM_ID\tHMM_name\tHMM_description\tHMM_Counts\tHMM_Percentage")
    k=0
    tophmms = []
    for line in print_list:
        OFO.write("\n"+"\t".join(str(k) for k in line[:4]) + "\t" + str(round(float(line[3])/total_bgc_count,4)))   #print data to to overall frequency table .tsv file
        if k < args.hmm_number:
            tophmms.append(line[4])
            k=k+1

    #metaset = list(set(metabait) & set(allpfamhits.keys()))

    OFO = open(args.OFO+"/proteinself_matrixdata.tsv", 'w')
    OFO2 = open(args.OFO+"/otherhit_matrixdata.tsv", 'w')
    OFO3 = open(args.OFO+"/allhit_matrixdata.tsv", 'w')
    OFO4 = open(args.OFO+"/phi_matrixdata.tsv", 'w')
    #OFO.write("\t"+"\t".join(metaset))
    #OFO2.write("\t"+"\t".join(metaset))
    #OFO3.write("\t"+"\t".join(metaset))
    OFO.write("\t"+"\t".join(tophmms))
    OFO2.write("\t"+"\t".join(tophmms))
    OFO3.write("\t"+"\t".join(tophmms))
    OFO4.write("\t"+"\t".join(tophmms))

    P = []
    Q = []
    R = []
    S = []
    for tophmm in tophmms:         
        curr_prot_row = []
        OFO_curr_prot_row = []
        curr_BGC_row = []
        OFO_curr_BGC_row = []
        total_BGC_row = []
        OFO_total_BGC_row = []
        phi_row = []
        OFO_phi_row = []
        #for tophmm2 in metaset:
        for tophmm2 in tophmms:
            #print(allpfamhits[tophmm], allpfamhits[tophmm2])
            if tophmm2 in self_hit_matrix[tophmm].keys():
                curr_prot_row.append(min(round(100*float(self_hit_matrix[tophmm][tophmm2])/allpfamhits[tophmm2][3], 2), 100.00))
                OFO_curr_prot_row.append(self_hit_matrix[tophmm][tophmm2])
            else:
                curr_prot_row.append(0.00)
                OFO_curr_prot_row.append(0)
            if tophmm2 in bgc_otherhit_matrix[tophmm].keys():
                #print(allpfamhits[tophmm], bgc_otherhit_matrix[tophmm2][tophmm], bgc_otherhit_matrix[tophmm][tophmm2], bgc_hit_counter[tophmm2])
                #print(100*float(bgc_otherhit_matrix[tophmm][tophmm2])/bgc_hit_counter[tophmm2])
                curr_BGC_row.append(min(round(100*float(bgc_otherhit_matrix[tophmm][tophmm2])/bgc_hit_counter[tophmm2], 2), 100.00))
                OFO_curr_BGC_row.append(bgc_otherhit_matrix[tophmm][tophmm2])
            else:
                curr_BGC_row.append(0.00)
                OFO_curr_BGC_row.append(0)
            if tophmm2 in bgc_allhit_matrix[tophmm].keys():
                total_BGC_row.append(min(round(100*float(bgc_allhit_matrix[tophmm][tophmm2])/(total_bgc_count), 2), 100.00))
                OFO_total_BGC_row.append(bgc_allhit_matrix[tophmm][tophmm2])
            else:
                total_BGC_row.append(0.00)
                OFO_total_BGC_row.append(0)
            try:
                phi = (total_bgc_count*bgc_allhit_matrix[tophmm][tophmm2]-bgc_hit_counter[tophmm2]*bgc_hit_counter[tophmm])/math.sqrt(bgc_hit_counter[tophmm2]*bgc_hit_counter[tophmm]*(total_bgc_count-bgc_hit_counter[tophmm2])*(total_bgc_count-bgc_hit_counter[tophmm]))
            except:
                phi = -10
            if phi > -9:
                phi_row.append(round(phi, 5))
                OFO_phi_row.append(round(phi, 5))
            else:
                phi_row.append(0.00)
                OFO_phi_row.append(0)


        OFO.write("\n"+tophmm+"\t"+"\t".join(str(x) for x in OFO_curr_prot_row))
        P.append(curr_prot_row)
        OFO2.write("\n"+tophmm+"\t"+"\t".join(str(x) for x in OFO_curr_BGC_row))
        Q.append(curr_BGC_row)
        OFO3.write("\n"+tophmm+"\t"+"\t".join(str(x) for x in OFO_total_BGC_row))
        R.append(total_BGC_row)
        OFO4.write("\n"+tophmm+"\t"+"\t".join(str(x) for x in OFO_phi_row))
        S.append(phi_row)

    similarity_matrix = cosine_similarity(np.stack(S))
    OFO5 = open(args.OFO + "/BiGBERD_similarity_matrix.tsv", 'w')
    OFO5.write("\t"+"\t".join(tophmms))
    for index in range(0, len(tophmms)):
        OFO5.write("\n"+tophmms[index] + "\t" + "\t".join(str(round(x, 4)) for x in similarity_matrix[index]))

    distarray = hdbscan.HDBSCAN(min_cluster_size=5, min_samples=2, cluster_selection_epsilon=0.8).fit(similarity_matrix)
    dendrogram = distarray.condensed_tree_.to_networkx()

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(similarity_matrix)

    pca = PCA(n_components=5)
    X_pca = pca.fit_transform(X_scaled)
    OFO5 = open(args.OFO + "/BiGBERD_pca_matrix.tsv", 'w')
    for index in range(0, len(tophmms)):
        OFO5.write("\n"+tophmms[index] + "\t" + "\t".join(str(round(x, 4)) for x in X_pca[index]))
    #X_train, X-test, y_train, y_test = train_test_split(X_pca, y, test_size=0.2, random_state=42)

    #model = LogisticRegression()
    #model.fit(X_train, y_train)

    #y_pred = model.predict(X_test)

    max_cluster_number = max(int(x) for x in distarray.labels_)
    cluster_counter = {}
    for cluster in range(0, max_cluster_number+1):
        cluster_counter[cluster] = (0, np.zeros(len(tophmms)))

    cluster_number = max_cluster_number
    for distindex in range(0, len(tophmms)):
        if int(distarray.labels_[distindex]) == -1:
            cluster_number = cluster_number + 1
            cluster_counter[cluster_number] = (1, np.array(S[distindex]))
        else:
            cluster_counter[distarray.labels_[distindex]] = (cluster_counter[distarray.labels_[distindex]][0] + 1, np.add(cluster_counter[distarray.labels_[distindex]][1], np.array(S[distindex])))
#            cluster_counter[distarray.labels_[distindex]][0] = cluster_counter[distarray.labels_[distindex]][0] + 1
#            cluster_counter[distarray.labels_[distindex]][1].append(tophmms[distindex])
#            cluster_counter[distarray.labels_[distindex]][2] = np.add(cluster_counter[distarray.labels_[distindex]][2], similarity_matrix[distindex])

    cluster_pos = []
    for cluster in range(0, cluster_number+1):
        cluster_pos.append(cluster_counter[cluster][1]/cluster_counter[cluster][0])

    similarity_matrix = cosine_similarity(np.stack(cluster_pos))
    OFO5 = open(args.OFO + "/BiGBERD_cluster_similarity_matrix.tsv", 'w')
    for index in range(0, cluster_number+1):
        OFO5.write("\n"+str(index) + "\t" + "\t".join(str(round(x, 4)) for x in similarity_matrix[index]))

    distarray_2 = hdbscan.HDBSCAN(min_cluster_size=5, min_samples=2, cluster_selection_epsilon=0.8).fit(cluster_pos)
    dendrogram_2 = distarray_2.condensed_tree_.to_networkx()

    OFO5 = open(args.OFO + "/BiGBERD_clusters.tsv", 'w')
    OFO5.write("HMM name\tCluster number\tMegacluster number")

    cluster_number = max_cluster_number
    #OFO5.write("\t"+"\t".join(tophmms))
    for distindex in range(0, len(tophmms)):
        if int(distarray.labels_[distindex]) == -1:
            cluster_number = cluster_number + 1
            distarray.labels_[distindex] = cluster_number
            OFO5.write("\n"+tophmms[distindex]+"\t"+ str(distarray.labels_[distindex]) + "\t" + str(distarray_2.labels_[cluster_number]))
            nx.set_node_attributes(dendrogram, {distindex: {"HMM index": int(distindex), "HMM": str(allpfamhits[tophmms[distindex]][0]), "HMM Name": str(allpfamhits[tophmms[distindex]][1]), "HMM Description": str(allpfamhits[tophmms[distindex]][2]), "Cluster Number": str(cluster_number), "Megacluster Number": str(distarray_2.labels_[cluster_number])}})
        else:
            OFO5.write("\n"+tophmms[distindex]+"\t"+ str(distarray.labels_[distindex]) + "\t" + str(distarray_2.labels_[distarray.labels_[distindex]]))
            nx.set_node_attributes(dendrogram, {distindex: {"HMM index": int(distindex), "HMM": str(allpfamhits[tophmms[distindex]][0]), "HMM Name": str(allpfamhits[tophmms[distindex]][1]), "HMM Description": str(allpfamhits[tophmms[distindex]][2]), "Cluster Number": str(distarray.labels_[distindex]), "Megacluster Number": str(distarray_2.labels_[distarray.labels_[distindex]])}})
    XGMMLWriter(open(args.OFO + "/BGCnetwork_test.xgmml", 'w'), dendrogram, "New Graph")

    #fig, ax = plt.subplots(figsize=(20,12))
    #annotated = set()
    #plt.scatter(X_pca[:, 0], X_pca[:, 1])
    
    for b in range(0,4):
        for c in range(b+1, 5):
            fig, ax = plt.subplots(figsize=(20,12))
            annotated = set()
            plt.scatter(X_pca[:, b], X_pca[:, c])
            for i in range(len(tophmms)):
                if distarray.labels_[i] not in annotated:
                    plt.annotate(distarray.labels_[i], (X_pca[i, b], X_pca[i, c]), fontsize=8)
                    annotated.add(distarray.labels_[i])
            plt.savefig(args.OFO+"/pca_matrixfig_PC"+str(b)+"_PC"+str(c)+".png", dpi=1200)
            plt.clf()
            
#figure generation
    if(args.hmm_number < 101):
        X = range(0,args.hmm_number)
        Y = range(0,args.hmm_number)

        fig, ax = plt.subplots(figsize=(16,14))
        plt.xlabel("Query HMM", fontsize=14)
        ax.xaxis.set_label_position('top')
        plt.ylabel("Target HMM", fontsize=14)
        plt.title("Percent of Proteins with Query HMM hit also containing Target HMM hit", pad=14, fontsize=20)
        ax.tick_params(direction='out', length=5, width=1.5, pad=1.5, top=True, labeltop=True, bottom=False, labelbottom=False)
        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(1)
            spine.set_color('black')
        ax.set_xticks(range(0,args.hmm_number))
        plotlabels = []
        for hmm in tophmms:
            if len(allpfamhits[hmm][0])<3:
                plotlabels.append(hmm)
            else:
                plotlabels.append(allpfamhits[hmm][0])
        ax.set_xticklabels(plotlabels)
        plt.xticks(rotation=45, ha='left', rotation_mode='anchor', fontsize=8)
        plt.yticks(fontsize=8)
        ax.set_yticks(range(0,args.hmm_number))
        ax.set_yticklabels(list(reversed(plotlabels)))

        
        plt.tight_layout()
        fig.set_tight_layout(True)

        temp_ax = ax.pcolormesh(X,Y,list(reversed(P)), shading='auto', cmap='Blues', vmin=1, vmax=100, lw=0.5, ec='white')
        cb = plt.colorbar(temp_ax, ax=ax, fraction = 0.03, aspect=8, label='Percent Frequency', pad=0.03)
        #cb.ax.text(0.5,1,"HMM total: %d\nHMM hits/protein: %d" % (args.hmm_number, args.number))
        fig.text(0.948, 0.97,"HMM hits/protein: %d\nHMM total: %d" % (args.number, args.hmm_number), ha='center', bbox=dict(facecolor='none', edgecolor='black', pad=5))
        plt.savefig(args.OFO+"/proteinself_matrixfig.png", dpi=600)

        plt.title("Percent of BGCs with Query HMM hit also containing discrete Target HMM hit", pad=14, fontsize=20)
        cb.remove()
        temp_ax = ax.pcolormesh(X,Y,list(reversed(Q)), shading='auto', cmap='Blues', vmin=1, vmax=100, lw=0.5, ec='white')
        cb = plt.colorbar(temp_ax, ax=ax, fraction = 0.03, aspect=8, label='Percent Frequency', pad=0.03)
        plt.savefig(args.OFO+"/BGCotherhit_matrixfig.png", dpi=600)

        plt.title("Percent of all BGCs with both Query HMM hit and Target HMM hit", pad=14, fontsize=20)
        cb.remove()
        temp_ax = ax.pcolormesh(X,Y,list(reversed(R)), shading='auto', cmap='Blues', vmin=0, vmax=100, lw=0.5, ec='white')
        cb = plt.colorbar(temp_ax, ax=ax, fraction = 0.03, aspect=8, label='Percent Frequency', pad = 0.03)
        plt.savefig(args.OFO+"/BGCallhit_matrixfig.png", dpi=600)

        plt.title("Phi coefficient for Query HMM and Target HMM correlation", pad=14, fontsize=20)
        cb.remove()
        temp_ax = ax.pcolormesh(X,Y,list(reversed(S)), shading='auto', cmap='bwr', vmin=-1, vmax=1, lw=0.5, ec='white')
        cb = plt.colorbar(temp_ax, ax=ax, fraction = 0.03, aspect=8, label='Phi value', pad = 0.03)
        plt.savefig(args.OFO+"/phi_matrixfig.png", dpi=600)

    end = time.time()
    print(str(end-st))

class XGMMLParserHelper(object):
    """
    """

    def __init__(self):
        """

        Arguments:
        - `graph`: Network X graph object
        """
        self._graph = nx.DiGraph()
        self._parser = xml.parsers.expat.ParserCreate()
        self._parser.StartElementHandler = self._start_element
        self._parser.EndElementHandler = self._end_element
        self._tagstack = list()

        self._network_att_el = dict()
        self._current_att_el = dict()
        self._current_list_att_el = list()
        self._current_obj = dict()

    def _start_element(self, tag, attr):
        """

        Arguments:
        - `self`:
        - `tag`:
        - `attr`:
        """

        self._tagstack.append(tag)

        if tag == 'graph':
            self._network_att_el = dict()

        if tag == 'node' or tag == 'edge':
            self._current_obj = dict(attr)

        if tag == 'att' and (self._tagstack[-2] == 'node' or
                             self._tagstack[-2] == 'edge'):
            if 'value' in attr:
                self._current_att_el = self._parse_att_el(self._current_att_el,
                                                          tag, attr)
            elif attr['type'] == 'list':
                self._current_list_name = attr['name']
                self._current_att_el[attr['name']] = list()

        if tag == 'att' and (self._tagstack[-2] == 'att'):
            self._current_list_att_el = dict(attr)
            if 'value' in attr:
                self._current_list_att_el = self._parse_att_el(
                    self._current_list_att_el, tag, attr)
                self._current_att_el[self._current_list_name].append(
                    self._current_list_att_el[attr['name']])

        if tag == 'att' and self._tagstack[-2] == 'graph':
            if 'value' in attr:
                self._network_att_el[attr['name']] = attr['value']

    def _parse_att_el(self, att_el, tag, attr):
        """

        Arguments:
        - `self`:
        - `att_el`: att element. Can be child of node, edge or another att.
        - `tag`:
        - `attr`:
        """

        if 'value' in attr:
            if attr['type'] == 'string':
                att_el[attr['name']] = attr['value']
            elif attr['type'] == 'real':
                att_el[attr['name']] = float(attr['value'])
            elif attr['type'] == 'integer':
                att_el[attr['name']] = int(attr['value'])
            elif attr['type'] == 'boolean':
                att_el[attr['name']] = bool(attr['value'])
            else:
                raise NotImplementedError(attr['type'])

            return att_el

    def _end_element(self, tag):
        """

        Arguments:
        - `self`:
        - `tag`:
        """

        if tag == 'node':
            if 'label' in self._current_obj:
                if 'label' in self._current_att_el:
                    self._current_att_el['@label'] = self._current_att_el['label']
                    del self._current_att_el['label']

                self._graph.add_node(self._current_obj['id'],
                                     label=self._current_obj['label'],
                                     **self._current_att_el)
            else:
                self._graph.add_node(self._current_obj['id'],
                                     **self._current_att_el)
        elif tag == 'edge':
            self._graph.add_edge(self._current_obj['source'],
                                 self._current_obj['target'],
                                 **self._current_att_el)

        self._tagstack.pop()

    def parseFile(self, file):
        """

        Arguments:
        - `self`:
        - `file`:
        """

        self._parser.ParseFile(file)

    def graph(self):
        """

        Arguments:
        - `self`:
        """

        return self._graph

    def graph_attributes(self):
        """

        Arguments:
        - `self`:
        """

        return self._network_att_el


def XGMMLReader(graph_file):
    """

    Arguments:
    - `file`:
    """

    parser = XGMMLParserHelper()
    parser.parseFile(graph_file)
    return parser.graph()


def XGMMLWriter(graph_file, graph, graph_name, directed=True):
    """

    Arguments:
    - `graph_file` output network file (file object)
    - `graph`: NetworkX Graph Object
    - `graph_name`: Name of the graph
    - `directed`: is directed or not
    """

    graph_file.write("""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<graph directed="{directed}"  xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns="http://www.cs.rpi.edu/XGMML">
 <att name="selected" value="1" type="boolean" />
 <att name="name" value="{0}" type="string"/>
 <att name="shared name" value="{0}" type="string"/>\n""".format(graph_name, directed=(1 if directed else 0)))

    def quote(text):
        """

        Arguments:
        - `text`:
        """

        return text.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

    def write_att_el(k, v, indent_count):
        indentation_string = ''
        for i in range(0, indent_count):
            indentation_string += '  '
        if isinstance(v, int):
            graph_file.write(
                indentation_string +
                '<att name="{}" value="{}" type="integer" />\n'.format(k, v))
        elif isinstance(v, np.int64):
            graph_file.write(
                indentation_string +
                '<att name="{}" value="{}" type="integer" />\n'.format(k, v))
        elif isinstance(v, bool):
            graph_file.write(
                indentation_string +
                '<att name="{}" value="{}" type="boolean" />\n'.format(k, v))
        elif isinstance(v, float):
            graph_file.write(
                indentation_string +
                '<att name="{}" value="{}" type="real" />\n'.format(k, v))
        #elif hasattr(v, '__iter__'):
        #    graph_file.write(
        #        indentation_string + '<att name="{}" type="list">\n'.format(k))
        #    for item in v:
        #        write_att_el(k, item, 3)
        #    graph_file.write(indentation_string + '</att>\n')
        else:
            graph_file.write(
                indentation_string +
                '<att name="{}" value="{}" type="string" />\n'.format(k,
                                                                    quote(v)))

    for onenode in graph.nodes(data=True):
        id = onenode[0]
        attr = dict(onenode[1])

        if 'label' in attr:
            label = attr['label']
            del attr['label']
        else:
            label = id

        graph_file.write(
            '  <node id="{id}" label="{label}">\n'.format(id=id, label=label))

        # Add color element
        if 'color' in attr:
            color = attr['color']
            del attr['color']
            graph_file.write(
                '  <graphics fill="{color}" />\n'.format(color=color))

        for k, v in iter(attr.items()):
            write_att_el(k, v, 2)

        graph_file.write('  </node>\n')

    for oneedge in graph.edges(data=True):
        #
        # The spec, http://cgi5.cs.rpi.edu/research/groups/pb/punin/public_html/XGMML/draft-xgmml.html#GlobalA,
        # requires an "id", even for edges. This id is supposed to be unique across the entire document, so it
        # can't be equal to one of the node ids. We're making the assumption that whoever created the graph
        # object knew about and respected the "uniqueness" requirement, and passed a suitable id as the attribute
        # "id" to the edge. If the creator of the graph *didn't* pass a unique id, the best I can come up with at
        # this moment is to just ignore the id requirement entirely.
        #
        if 'id' in oneedge[2]:
            edge_id = oneedge[2].pop("id", None)
            graph_file.write('  <edge id="{}" source="{}" target="{}">\n'.format(
                edge_id, oneedge[0], oneedge[1]))
        else:
            graph_file.write('  <edge source="{}" target="{}">\n'.format(
                oneedge[0], oneedge[1]))

        for k, v in iter(oneedge[2].items()):
            write_att_el(k, v, 2)
        graph_file.write('  </edge>\n')
    graph_file.write('</graph>\n')

def write_record(main_html, record, filter_list):
    #RESULTS FOR xxxx
    #DRAW ORF
    #DRAW ORF scale
    #PUT LINK TO nuc SEQUENCE
    #TABLE of CDS
    #TABLE of ORFs
    if record.locus is not -1:
        main_html.write('<h2 id="%s"> Results for Locus %s in Record %s [%s]\n' % (record.cluster_accession, record.locus, record.cluster_accession, record.cluster_genus_species))
    else:
        main_html.write('<h2 id="%s"> Results for %s [%s]\n' % (record.query_accession_id, record.query_accession_id, record.cluster_genus_species))
    main_html.write('<a href="#header"><small><small>back to top</small></small></a></h2>') #TODO keep for single?
    main_html.write('<p></p>') 
    draw_orf_diagram(main_html, record, filter_list)
    main_html.write('<p></p>') 
    main_html.write('<a href="https://www.ncbi.nlm.nih.gov/nuccore/%s">Link to nucleotide sequence</a>' % (record.cluster_accession))
    draw_cds_table(main_html, record)

def write_header(html_file):
    html_file.write("""
    <html>
    <head>
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css">
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap-theme.min.css">
    </head>
    <style media="screen" type="text/css">
     
    .square {
      width: 54px;
      height: 14px;
      background-color: white;
      outline: #ffffff solid 1px;
      text-align: center;
      line-height: 14px;
      font-size: 12px;
    }
     
    table {
       font-size: 11px;
    }
     
    </style>
     
    <script src='https://img.jgi.doe.gov//js/overlib.js'></script>
    <div class="container">
    <h1 align="center" id="header">RODEO2</h1>
    """)

def get_fill_color(cds, filter_list):
    for i in range(len(cds.pfam_descr_list)):
        if cds.pfam_descr_list[i][3].upper() in filter_list:
            return "red"
        if cds.pfam_descr_list[i][0].upper() in filter_list:
            return "red"
    return "white"

def draw_CDS_arrow(main_html, cds, filter_list, sub_by, scale_factor):
    fill_color = get_fill_color(cds, filter_list)
    start = cds.start
    end = cds.end
    #HMM info
    if len(cds.pfam_descr_list) == 0:
        pfamID = "No Pfam match"
        pfam_desc = ""
    else:
        pfamID = cds.pfam_descr_list[0][0].split('.')[0] #No need for version?
        pfam_desc = cds.pfam_descr_list[0][1]
    main_html.write('<polygon points=\"')
    arrow_wid = int((start - sub_by) * scale_factor)
    arrow_wid3 = int((end - sub_by) * scale_factor)
    if abs(arrow_wid3 - arrow_wid) < 40:
        arrow_wid2 = (arrow_wid + arrow_wid3) / 2 #middle?
    else:
        if start < end:
            arrow_wid2 = arrow_wid3 - 20
        else:
            arrow_wid2 = arrow_wid3 + 20
        
    str_arrow_wid = str(arrow_wid)
    str_arrow_wid2 = str(arrow_wid2)
    str_arrow_wid3 = str(arrow_wid3)
    main_html.write(str_arrow_wid + ",10 " + str_arrow_wid + ",40 "\
                        + str_arrow_wid2 + ",40 " + str_arrow_wid2 + ",50 "+ str_arrow_wid3 \
                        + ",25 " + str_arrow_wid2 + ",0 " + str_arrow_wid2 + ",10 " + str_arrow_wid + ",10")
    main_html.write('" style="fill:' + fill_color)
    main_html.write(';stroke:black;stroke-width:.5" onMouseOver="return overlib(')
    main_html.write("'" + cds.accession_id + " - " + pfamID + " : " + pfam_desc + "'")
    main_html.write(')" onMouseOut="return nd()"/>')

def draw_orf_diagram(main_html, record, filter_list):
    main_html.write('<h3>Architecture</h3>\n')
    main_html.write('<svg width="1060" height="53">')
    bsc_start = min(record.CDSs[0].start, record.CDSs[0].end)
    bsc_end = max(record.CDSs[-1].end, record.CDSs[-1].start)
    sub_by = max(bsc_start - 500, 0)
    scale_factor = (660./(bsc_end - bsc_start))
    for cds in record.CDSs:
        draw_CDS_arrow(main_html, cds, filter_list, sub_by, scale_factor)
    main_html.write('</svg>')
    bar_length = scale_factor * 1000
    bar_legx = bar_length + 5
    main_html.write('<svg width="500" height="23">')
    main_html.write('<polygon points="')
    main_html.write("0,10 %f, 10" % (bar_length))
    main_html.write('" style="fill:white;stroke:black;stroke-width:.5" />')
    main_html.write('<text x="%f"y="12"' % (bar_legx))
    main_html.write("""font-famil="sans-serif"
                    font-size="10px"
                    text_anchor="right"
                    fill="black">1000 nucleotides</text>""")
    main_html.write('<polygon points="')
    main_html.write("0,5 0,15")
    main_html.write('" style="fill:white;stroke:black;stroke-width:.5" />')
    main_html.write('<polygon points="%f,5 %f,15' % (bar_length, bar_length))
    main_html.write('" style="fill:white;stroke:black;stroke-width:.5" />')
    main_html.write('</svg>')
    
def draw_cds_table(main_html, record):
    main_html.write("""<br><br><table class="table table-condensed">
  <tbody>
    <tr>
      <th scope="col">Accession</th>
      <th scope="col">start</th>
      <th scope="col">end</th>
      <th scope="col">direction</th>
      <th scope="col">length (aa)</th>
      <th scope="col">Pfam/HMM</th>
      <th scope="col">name</th>
      <th scope="col">description</th>
      <th scope="col">E-value</th>    
    </tr>""")
    for cds in record.CDSs:
        main_html.write("<tr>\n")
        if cds.meta:
            prot_or_nucc = "nuccore"
            acc = record.cluster_accession
        else:
            prot_or_nucc = "protein"
            acc = cds.accession_id
        main_html.write("""\t<td><a href='https://www.ncbi.nlm.nih.gov/%s/%s'>%s</a></td>
            <td>%s</td> 
            <td>%d</td>
            <td>%s</td>
            <td>%d</td>""" % (prot_or_nucc, acc, cds.accession_id, cds.start, cds.end, cds.direction, round(abs(cds.end-cds.start)/3, 0)))
        if len(cds.pfam_descr_list) == 0:
            main_html.write("<td>NO PFAM MATCH</td>")
            main_html.write("<td>-</td>")
            main_html.write("<td>-</td>")
            main_html.write("<td>-</td>")
        else:
            if cds.pfam_descr_list[0][0][:2] == "PF":
                main_html.write("<td><a href='https://www.ebi.ac.uk/interpro/entry/pfam/%s'>%s</a>" % (cds.pfam_descr_list[0][0], cds.pfam_descr_list[0][0]))
            elif cds.pfam_descr_list[0][0][:4] == "TIGR":
                main_html.write("<td><a href='https://www.ebi.ac.uk/interpro/entry/ncbifam/%s'>%s</a>" % (cds.pfam_descr_list[0][0], cds.pfam_descr_list[0][0]))
            else:
                main_html.write("<td>%s" % (cds.pfam_descr_list[0][0]))
            n = 5
            for pfamid, _, _, _, in cds.pfam_descr_list[1:n]:
                if pfamid[:2] == "PF":
                    main_html.write("<br><a href='https://www.ebi.ac.uk/interpro/entry/pfam/%s'>%s</a>" % (pfamid, pfamid))
                elif pfamid[:4] == "TIGR":
                    main_html.write("<br><a href='https://www.ebi.ac.uk/interpro/entry/ncbifam/%s'>%s</a>" % (pfamid, pfamid))
                else:
                    main_html.write("<br>%s" % (pfamid))
                    
            main_html.write("</td><td>%s"%(cds.pfam_descr_list[0][3]))
            for _, _, _, name in cds.pfam_descr_list[1:n]:
                main_html.write("<br>%s"%(name))
                
            descr = cds.pfam_descr_list[0][1]
            main_html.write("</td><td>%s" % (descr))
            for _, descr, _, _, in cds.pfam_descr_list[1:n]:
                main_html.write("<br>%s" % (descr))
                
            e_val = cds.pfam_descr_list[0][2]
            main_html.write("</td><td>%.2E" % Decimal(e_val))
            for _, _, e_val, _, in cds.pfam_descr_list[1:n]:
                main_html.write("<br>%.2E" % Decimal(e_val))
            
            main_html.write("</td>")
        main_html.write("</tr>") 
    main_html.write("</tbody></table><p></p>")

class Sub_Seq(object):
    """Useful for storing subsequences and their coordinates"""
    def __init__(self, start, end, accession_id=None, meta=False):
        self.start = start
        self.end = end
        if int(self.end) > int(self.start):
            self.direction = "+"
        else:
            self.direction = "-"
        self.accession_id = accession_id
        self.pfam_descr_list = []
        self.meta = meta


class My_Record(object):
    """ """
    #TODO get genus and species frecom gb file
    def __init__(self, query_accession_id):
        self.query_accession_id = query_accession_id
        self.cluster_accession = ""
        self.cluster_genus_species = ""
        self.CDSs = []
        self.locus = -1


if __name__=="__main__":
    __main__()