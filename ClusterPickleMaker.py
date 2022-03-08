import os
import numpy as np
import scipy as sp
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import minimize
import math
from itertools import combinations
from scipy.sparse import csr_matrix
import pickle

dic_DNA={
    'A': 0,
    'C': 1,
    'G': 2,
    'T': 3,
    '-': 4,
    'N': 5
}

dic_DNA_binary={
    'A': 0,
    'C': 0,
    'G': 0,
    'T': 0,
    '-': 1,
    'N': 0
}


def isDNA(c):
    if c in {'A', 'C', 'G', 'T', 'N'}:
        return True
    else:
        return False

def jci(dis):
    return .75 * (1 - math.exp((-4.0* dis)/ 3 ))

def parseFasta(alignment_file="/Users/omersabary/Desktop/GenePlacement/3000/data/381/ref.fa"):
    list_alignment=[]
    species_list=[]
    dic_genes={}
    file_alignment=open(alignment_file)
    line_next=""
    for i, line in enumerate(file_alignment):
        if line.startswith('>'):
            species_list.append(line)
            if i!=0:
                list_alignment.append(line_next)
                dic_genes[line[1:-1]] = line_next.rstrip()
                line_next=""
            #continue
        else:
            line_next=line_next+line.rstrip()
            if(line_next.startswith('>')):
                print(line_next)
                exit(0)

    return list_alignment, species_list, dic_genes

def symbol_similarity(a, b):
    #print(a)
    #print(b)
    if isDNA(a) and isDNA(b):
        if a != b:
            return 1
        else:
            return 0
    else:
        return 0

def is_both_DNA(a, b):
    #print(a)
    #print(b)
    if isDNA(a) and isDNA(b):
        return 1
    else:
        return 0

def createLenArray(genome_dic, tree_dis_file, lens_name):
    df_tree = pd.read_csv(tree_dis_file)
    pairs = df_tree['pair'].to_list()
    #number_of_pairs= len(pairs)
    pairs_species=list(combinations(genome_dic.keys(),2))
    number_of_pairs=(len(pairs_species))
    number_of_pairs= len(pairs)
    #print(number_of_pairs)
    genome_len =len(list(genome_dic.values())[0])
    #hamming_dis_matrix=np.ndarray((number_of_pairs, genome_len))
    #hamming_dis_matrix=csr_matrix((number_of_pairs, genome_len), dtype=np.int8)
    #print(hamming_dis_matrix.toarray())
    lens=np.zeros(number_of_pairs)
    print(len(lens))
    print("sskj")
    #exit(0)
    for i, pair in enumerate(pairs):
        #first_spec=pair.__getitem__(0)
        first_spec = pair[2:12]
        #print(first_seq)
        #exit(0)
        first_seq=genome_dic[first_spec]
        #sec_spec=pair.__getitem__(1)
        sec_spec = pair[15:25]
        sec_seq=genome_dic[sec_spec]
        for j in range(genome_len):
            lens[i] = lens[i] + is_both_DNA(first_seq[j], sec_seq[j])
    np.save(lens_name, lens, allow_pickle=True)
    return lens
    #exit(0)
def createHammingDisMatrix(genome_dic, matrix_name, tree_dis_file):
    df_tree = pd.read_csv(tree_dis_file)
    pairs = df_tree['pair'].to_list()
    #number_of_pairs= len(pairs)
    pairs_species=list(combinations(genome_dic.keys(),2))
    number_of_pairs=(len(pairs_species))
    number_of_pairs= len(pairs)
    print(number_of_pairs)
    genome_len =len(list(genome_dic.values())[0])
    hamming_dis_matrix=np.ndarray((number_of_pairs, genome_len))
    #hamming_dis_matrix=csr_matrix((number_of_pairs, genome_len), dtype=np.int8)
    #print(hamming_dis_matrix.toarray())
    #exit(0)
    for i, pair in enumerate(pairs):
        #first_spec=pair.__getitem__(0)
        first_spec = pair[2:12]
        #print(first_seq)
        #exit(0)
        first_seq=genome_dic[first_spec]
        #sec_spec=pair.__getitem__(1)
        sec_spec = pair[15:25]
        sec_seq=genome_dic[sec_spec]
        for j in range(genome_len):
            hamming_dis_matrix[i][j]=symbol_similarity(first_seq[j], sec_seq[j])
    np.save(matrix_name, hamming_dis_matrix, allow_pickle=True)
    #exit(0)


def set_multipliers(sitemultipliers, multipliers, m, set_mult):
    sitemultipliers = sorted(sitemultipliers, reverse=True)
      # TODO: Taking this line out of the function
    for i in range(len(set_mult)):
        multipliers[multipliers == set_mult[i] ] = sitemultipliers[i]
    #print(multipliers[30:40])
    return multipliers

def target_function_L1(sitemultipliers, *args):
    print(sitemultipliers)
    list_args=list(args[0].values())
    multipliers=list_args[1]
    trees_distances=list_args[2]
    hm_matrix=list_args[3]
    lens = list_args[4]
    pairs = list_args[5]

    multipliers=set_multipliers(sitemultipliers, multipliers, len(sitemultipliers))
    print("max:")
    print(np.max(multipliers))
    multipliers = np.transpose(multipliers)
    print_mul = np.matmul(hm_matrix, multipliers)
    #hm_dis = (print_mul.sum())

    sum = 0
    for i, pair in enumerate(pairs):
        sum = sum + abs((trees_distances[i]) - (print_mul[i]/lens[i]) ) #TODO: ADD a Length
    print(sitemultipliers)
    print(np.median(sitemultipliers))

    return sum

def target_function_L2(sitemultipliers, *args):
    print(sitemultipliers)
    list_args=list(args[0].values())
    multipliers=list_args[1]
    trees_distances=list_args[2]
    hm_matrix=list_args[3]
    lens = list_args[4]
    pairs = list_args[5]
    set_mult = list_args[6]

    multipliers=set_multipliers(sitemultipliers, multipliers, len(sitemultipliers), set_mult)
    #print("max:")
    #print(np.max(multipliers))
    multipliers = np.transpose(multipliers)
    print_mul = np.matmul(hm_matrix, multipliers)
    hm_dis = (print_mul.sum())

    sum = 0
    for i, pair in enumerate(pairs):
        sum = sum + ((trees_distances[i]) - (print_mul[i]/lens[i]) )**2 #TODO: ADD a Length
    print(sitemultipliers)
    print(sum/len(pairs))
    #print(np.median(sitemultipliers))

    return sum/len(pairs)

def concatination_list_of_gene(genes_alignment_dir="/Users/omersabary/Desktop/GenePlacement/3000/alignments/ffn/", query_species_list=[]):
    dic_concat={}
    for spec in query_species_list:
        dic_concat[spec] = ""
    #print(dic_concat)
    #exit(0)
    current_spec=""
    species_list = []
    line_next =""
    len_line= 0
    for filename in os.listdir(genes_alignment_dir):
        if filename.endswith(".fasta"):
            file_alignment = open(genes_alignment_dir+filename, "r")
            for i, line in enumerate(file_alignment.readlines()):
                if line.startswith('>'):
                    if i==0:
                        current_spec=line[1:-1]
                    species_list.append(line[1:-1])
                    if i != 0:
                        if current_spec in dic_concat.keys():
                            dic_concat[current_spec]= dic_concat[current_spec]+line_next
                        else:
                            if current_spec in query_species_list:
                                dic_concat[current_spec] = line_next
                            #else:
                            #    print(current_spec)
                        if current_spec in dic_concat.keys():
                            if len(dic_concat[current_spec]) > len_line:
                                len_line = len(dic_concat[current_spec])
                        current_spec = line[1:-1]
                        line_next = ""
                else:
                    line_next = line_next + line.rstrip()
                    if (line_next.startswith('>')):
                        print(line_next)
                        exit(0)
            for spec in dic_concat.keys():
                if spec not in species_list or len(dic_concat[spec])<len_line:
                    gaps=""
                    for i in range(len_line-len(dic_concat[spec])):
                        gaps = gaps +"-"
                    #gaps = *len_line)
                    if spec in dic_concat.keys():
                        dic_concat[spec]=dic_concat[spec]+gaps
                    else:
                        dic_concat[spec] = gaps
            species_list = []
    prev_len=len(list(dic_concat.values())[0])
    for i, alignment in enumerate(dic_concat.values()):
        if len(alignment) != prev_len:
            print("ERO " + str(i))
            #print(i)
            #print(len(alignment))
            #print(str(prev_len)+"\n")
            #exit(0)
        prev_len=len(alignment)
    return dic_concat

#a, b, dic_genes = parseFasta()
#print(len(b))
cluster_num="202"

print("query starts")
query_species_file = "/Users/omersabary/Desktop/GenePlacement/clusters_omer/cluster"+cluster_num+".txt"
query_s = []
query_species_f = open(query_species_file, "r")
for line in query_species_f:
    query_s.append(line.rstrip()[:10])
print(query_s)
print(len(query_s))
print("query ends")
print("dic starts")

# a, b, dic_genes = parseFasta()
dic_genes = concatination_list_of_gene(genes_alignment_dir="/Users/omersabary/Desktop/GenePlacement/alignments/ffn/",
                                      query_species_list=query_s)
print("dic ends")

#Store data (serialize)
with open('dic_genes_cluster'+cluster_num+'.pickle', 'wb') as handle:
    pickle.dump(dic_genes, handle, protocol=pickle.HIGHEST_PROTOCOL)
# Load data (deserialize)
print("dic ends pickling")

#with open('dic_genes_cluster118.pickle', 'rb') as handle:
#    dic_genes = pickle.load(handle)
#exit(0)
print("lens array start")
lens = createLenArray(genome_dic=dic_genes, tree_dis_file="tree_distance_pairs_new_clusters_cluster"+cluster_num+".csv", lens_name="cluster"+cluster_num+"_lens.npy")
print(len(lens))
print("len ends")

np.save("cluster"+cluster_num+"_lens.npy", lens, allow_pickle=True)
print("lens array pickling ends")

with open('dic_genes_cluster'+cluster_num+'.pickle', 'rb') as handle:
    dic_genes = pickle.load(handle)
print("hm mat starts")
createHammingDisMatrix(genome_dic=dic_genes, matrix_name="cluster"+cluster_num+"_hamming.npy", tree_dis_file="tree_distance_pairs_new_clusters_cluster"+cluster_num+".csv")
print("hm mat ends")
exit(0)


