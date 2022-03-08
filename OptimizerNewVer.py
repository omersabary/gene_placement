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
    #sitemultipliers = sorted(sitemultipliers, reverse=True)
    #set_mult = sorted(list(set(multipliers))) # TODO: Taking this line out of the function
    print("len mul is :" +str(len(set_mult)))
    multipliers_new=np.zeros(len(multipliers))
    for i in range(len(set_mult)):
        multipliers_new[multipliers == set_mult[i] ] = sitemultipliers[i]
    #print(multipliers[30:40])
    return multipliers_new

def target_function_L1(sitemultipliers, *args):
    print(sitemultipliers)
    list_args=list(args[0].values())
    multipliers=list_args[1]
    trees_distances=list_args[2]
    hm_matrix=list_args[3]
    lens = list_args[4]
    pairs = list_args[5]

    multipliers_updated=set_multipliers(sitemultipliers, multipliers, len(sitemultipliers))
    print("max:")
    print(np.max(multipliers_updated))
    multipliers = np.transpose(multipliers_updated)
    print_mul = np.matmul(hm_matrix, multipliers_updated)
    #hm_dis = (print_mul.sum())

    sum = 0
    for i, pair in enumerate(pairs):
        sum = sum + abs((trees_distances[i]) - (print_mul[i]/lens[i]) ) #TODO: ADD a Length
    #print(sitemultipliers)
    #print(np.median(sitemultipliers))

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

    multipliers_updated=set_multipliers(sitemultipliers, multipliers, len(sitemultipliers), set_mult)
    #print("max:")
    #print(np.max(multipliers))
    multipliers_updated = np.transpose(multipliers_updated)
    print_mul = np.matmul(hm_matrix, multipliers_updated)
    hm_dis = (print_mul.sum())

    sum = 0
    for i, pair in enumerate(pairs):
        sum = sum + (jci(trees_distances[i]) - (print_mul[i]/lens[i]) )**2 #TODO: ADD a Length
    print(sitemultipliers)
    print(math.sqrt(sum))
    #print(np.median(sitemultipliers))

    return math.sqrt(sum)

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
clutser_num="118"
print("query starts")
query_species_file = "/Users/omersabary/Desktop/GenePlacement/clusters_omer/cluster"+clutser_num+".txt"
query_s = []
query_species_f = open(query_species_file, "r")
for line in query_species_f:
    query_s.append(line.rstrip()[:10])
print(query_s)
print(len(query_s))
print("query ends")

# a, b, dic_genes = parseFasta()
#dic_genes = concatination_list_of_gene(genes_alignment_dir="/Users/omersabary/Desktop/GenePlacement/alignments/ffn/",
                 #                      query_species_list=query_s)

## Store data (serialize)
#with open('dic_genes_cluster118.pickle', 'wb') as handle:
#    pickle.dump(dic_genes, handle, protocol=pickle.HIGHEST_PROTOCOL)
# Load data (deserialize)
with open('dic_genes_cluster'+clutser_num+'.pickle', 'rb') as handle:
    dic_genes = pickle.load(handle)
#exit(0)
#f = open("dic_genes_118.fasta", "w")
#f.write('\n'.join(">%s\n%s"%(k,v) for k,v in dic_genes.items()))
#exit(0)
print("lens array start")
#lens = createLenArray(genome_dic=dic_genes, tree_dis_file="tree_distance_pairs_new_clusters_cluster118.csv", lens_name="cluster118_lens.npy")
#print(len(lens))
#np.save("cluster118_lens.npy", lens, allow_pickle=True)

#print(lens)
lens = np.load("cluster"+clutser_num+"_lens.npy")
print("lens array ends")
print(lens.shape)
print("hm mat start")
rate_multipliers_files="site_rates/total_rates/sites_rates_of_total_7_levels_querires_cluster"+clutser_num+".csv"
df_multipliers = pd.read_csv(rate_multipliers_files)

#createHammingDisMatrix(genome_dic=dic_genes, matrix_name="cluster118_hamming.npy", tree_dis_file="tree_distance_pairs_new_clusters_cluster118.csv")
hm_matrix=np.load("cluster"+clutser_num+"_hamming.npy")
#print(hm_matrix)
print("hm mat ends")
print(hm_matrix.shape)




tree_dis_file = "tree_distance_pairs_new_clusters_cluster"+clutser_num+".csv"

multipliers_values =[0.1, 0.2, 1, 1.95, 1.99]

#multipliers_values = [1,1,1,1,1]

#multipliers_values = [0.5, 0.75, 1, 1.25, 1.5]
#df_lens = pd.read_csv("pairs_len_new_clusters.csv")
#lens = df_lens['total_len'].to_list()
df_trees = pd.read_csv(tree_dis_file)
pairs = df_trees['pair'].to_list()
tree_distances = df_trees['tree_dis'].to_list()
tree_distances = [jci(val) for val in tree_distances]
######multipliers_values=sorted(set(df_multipliers['rate'].to_list()))
multipliers_values=list(set(df_multipliers['rate'].to_list()))

print(multipliers_values)
#####multipliers_values= sorted([float( mul) for mul in multipliers_values], reverse=True)
print(multipliers_values)
#multipliers_values = [0.97, 0.98, 0.99, 1.0, 1.01, 1.02, 1.03]
#multipliers_values= [float( mul [2:-2] ) for mul in multipliers_values]
#multipliers_values = [1.22 , 1.22 , 1.22]
multipliers = df_multipliers['rate'].to_numpy()

print(multipliers_values)
sorted(list(set(multipliers)))
args_tuple={'dicgene': dic_genes,
      'mul': multipliers,
      'treedisfile': tree_distances,
            'hm_mat': hm_matrix,
            'lens': lens,
            'pairs': pairs,
            'set_mult': (list(set(multipliers)))}


cons = [{'type': 'ineq', 'fun': lambda x:  x },  #values are strictly positive
       # {'type': 'ineq', 'fun': lambda x: x[3] - 1},  # x[3] >=1
        #{'type': 'ineq', 'fun': lambda x: x[2] - 1 },  #x[2] >=1
        # {'type': 'ineq', 'fun': lambda x: x[4] - 1 },
        #{'type': 'ineq', 'fun': lambda x: 1-x[1]  },
        #{'type': 'ineq', 'fun': lambda x: 1-x[0]}, #x[0] <=1
        #{'type': 'eq', 'fun': lambda x: x[1] - 1 }, #x[1] == 1
       #{'type': 'ineq', 'fun': lambda x: x[0] - x[1] }, #x[0]>x[1]
       # {'type': 'ineq', 'fun': lambda x: x[1] - x[2]},
       # {'type': 'ineq', 'fun': lambda x: x[2] - x[3]},
       # {'type': 'ineq', 'fun': lambda x: x[3] - x[4]},
       # {'type': 'ineq', 'fun': lambda x: x[4] - x[5]},
       # {'type': 'ineq', 'fun': lambda x: x[5] - x[6]},

        #{'type': 'ineq', 'fun': lambda x: 1- np.median(x) } #median <=1
        # {'type': 'ineq', 'fun':lambda  x: x[4]-x[3]-x[2]-x[1]-x[0]}
        ]
b=(0,5) #sp.optimize.Bounds(lb=0,ub=3)

bnds=[b,b,b, b, b,b,b]#,b,b, b, b]
res = minimize(fun=target_function_L2, x0=multipliers_values, args=args_tuple, constraints=cons, method='trust-constr', bounds=bnds)

print(res)
print(res.x)
print(np.median(res.x))

