import os
import numpy as np
import scipy as sp
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import minimize
import math
from sklearn.cluster import KMeans
#from sklearn.cluster.tests.test_k_means import centers
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
                dic_genes[line[1:-1]] = line_next
                line_next=""
            #continue
        else:
            line_next=line_next+line.rstrip()
            if(line_next.startswith('>')):
                print(line_next)
                exit(0)
            #print(len(line_next))
            #print(line)
    return list_alignment, species_list, dic_genes

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

def gene_lists_alignment(alignment_file):
    list_alignment=[]
    species_list=[]
    file_alignment=open(alignment_file)
    line_next=""
    for i, line in enumerate(file_alignment):
        if line.startswith('>'):
            species_list.append(line)
            if i!=0:
                list_alignment.append(line_next)
                line_next=""
            #continue
        else:
            line_next=line_next+line.rstrip()
            if(line_next.startswith('>')):
                print(line_next)
                exit(0)
            #print(len(line_next))
            #print(line)
    return list_alignment, species_list

def computeEntropy(dic_i):
    entropy=0
    all_zeros=True
    for key in dic_i.keys():
        if dic_i[key]!=0:
            entropy=entropy+(dic_i[key] * np.log2(dic_i[key]))
            all_zeros=False
    if all_zeros:
        entropy=10000
    entropy = -1 *entropy
    return entropy

def computeEntropies(list_alignment, species_list):
    #print((list_alignment[1]))
    sites_entropies=[]
    for i in range(len(list_alignment[0])):
        dic_i={'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for j in range(len(list_alignment)):
            #if i >
            if list_alignment[j][i] in dic_i:
                dic_i[list_alignment[j][i]] = dic_i[list_alignment[j][i]] + 1
            else:
                if list_alignment[j][i] in {'-', 'N'}:
                    continue
                else:
                    print(list_alignment[j])
                    print(dic_i.get(list_alignment[j][i]))
                    print("wrong character "+ str(list_alignment[j][i]))
            #dic_i[list_alignment[j][i]]=dic_i[list_alignment[j][i]]+1
        dic_i_sum=np.sum(list(dic_i.values()))
        #dic_i_sum=dic_i_sum - dic_i['-']
        #dic_i_sum=dic_i_sum - dic_i['N']
        dic_i_new = {}
        DNA_alp = ['A', 'C', 'G', 'T']
        for key in DNA_alp:
            if dic_i_sum != 0:
                dic_i_new[key] = float(dic_i[key]) / dic_i_sum

        #if computeEntropy(dic_i) <1.01 and computeEntropy(dic_i)>0.99:
        #    print(dic_i)
        #    print(computeEntropy(dic_i))
        #if computeEntropy(dic_i) >2.01:
        #    print(dic_i)
        #    print(computeEntropy(dic_i))
        sites_entropies.append(computeEntropy(dic_i_new))
    return sites_entropies

def site_rate_calculator_m_levels(gene_alignment_folder='/Users/omersabary/Desktop/GenePlacement/3000/alignments/ffn', m=5):
    for filename in os.listdir(gene_alignment_folder):
        if filename.endswith(".fasta"):
            filepath = os.path.join(gene_alignment_folder, filename)
            list_alig, spec_list=gene_lists_alignment(filepath)
            sites_entropies=computeEntropies(list_alig, spec_list)
            sites_entropies_np = np.asarray(sites_entropies)
            number_of_parts=m
            percentiles=[]
            for i in range(1, m):
                p=np.percentile(sites_entropies_np, 100.*float(i)/number_of_parts, interpolation="nearest")
                percentiles.append(p)
            percentiles.append(max(sites_entropies_np))
            print(len(sites_entropies_np))
            data_splitted = []
            tocuhed = [False] * len(sites_entropies_np)
            for value_of_perc in percentiles:
                current_perc=[]
                for ind, val in enumerate(sites_entropies_np):
                    if val <=value_of_perc and not tocuhed[ind]:
                        current_perc.append(ind)
                        tocuhed[ind]=True
                data_splitted.append(current_perc)
            site_multipliers_options= np.arange(0.5, 1.55, 1.0/(m-1))
            site_multipliers = [0] * len(sites_entropies_np)
            for i in range(len(data_splitted)):
                for j in range(len(data_splitted[i])):
                    site_multipliers[data_splitted[i][j]]= site_multipliers_options[i]

            print(site_multipliers_options)
            print(site_multipliers)
            file_of_rates=open("site_rates/sites_rates_of_"+filename[:-6]+"_"+str(m)+"_levels.csv","w")
            file_of_rates.write("site,rate\n")
            for i,rate in enumerate(site_multipliers):
                file_of_rates.write(str(i)+","+str(rate)+"\n")
            file_of_rates.close()
            continue
        else:
            continue


def total_site_rate_calculator_m_queries(m, cluster_num):
    query_species_file = "/Users/omersabary/Desktop/GenePlacement/clusters_omer/cluster"+str(cluster_num)+".txt"
    query_s = []
    query_species_f = open(query_species_file, "r")
    for line in query_species_f:
        query_s.append(line.rstrip()[:10])
    print(query_s)
    print(len(query_s))
    #a, b, dic_genes = parseFasta()
    #dic_genes = concatination_list_of_gene(genes_alignment_dir="/Users/omersabary/Desktop/GenePlacement/alignments/ffn/",query_species_list=query_s)
    with open('dic_genes_cluster'+str(cluster_num)+'.pickle', 'rb') as handle:
        dic_genes = pickle.load(handle)
    #print(dic_genes)
    #exit(0)
    list_genes=[]
    for key in dic_genes.keys():
        if key in query_s:
            list_genes.append(dic_genes[key])
    site_enropies = computeEntropies(list_genes, 0)
    #print(site_enropies)
    sites_entropies_np = np.asarray(site_enropies)
    plt.hist(sites_entropies_np)
    plt.show()

    sites_entropies_np_reshaped = sites_entropies_np.reshape(-1,1)
    '''
    kmeans=KMeans(n_clusters=m, random_state=0).fit(sites_entropies_np_reshaped)

    centers =list(kmeans.cluster_centers_)
    print(centers)
    kmeans=kmeans.predict(sites_entropies_np_reshaped)
    print(len(kmeans))
    print(len(sites_entropies_np))
    #for val in kmeans:
    #    print(val)

    plt.hist(kmeans)
    plt.show()

    site_multi = []
    for ind in kmeans:
        site_multi.append(centers[ind])
    site_multi = np.asarray(site_multi)
    plt.hist(site_multi)
    site_multipliers = site_multi
    plt.show()

    #exit(0)
    '''
    number_of_parts = m
    percentiles = []
    sites_entropies_np_orig=sites_entropies_np
    sites_entropies_np = sites_entropies_np[sites_entropies_np > -100]
    plt.hist(sites_entropies_np)
    plt.show()
    for i in range(1, m):
        p = np.percentile(sites_entropies_np, 100. * 1.0 / (number_of_parts), interpolation="nearest")
        sites_entropies_np = sites_entropies_np[sites_entropies_np > p]
        percentiles.append(p)
    percentiles.append(max(sites_entropies_np))
    print(percentiles)
    #exit(0)
    sites_entropies_np=sites_entropies_np_orig
    print(len(sites_entropies_np))
    data_splitted = []
    tocuhed = [False] * len(sites_entropies_np)
    for value_of_perc in percentiles:
        current_perc = []
        for ind, val in enumerate(sites_entropies_np):
            if val <= value_of_perc and not tocuhed[ind]:
                current_perc.append(ind)
                tocuhed[ind] = True
        data_splitted.append(current_perc)
    site_multipliers_options = np.arange(0.5, 1.55, 1.0 / (m - 1))
    #print(set(site_multipliers_options))
    #print(site_multipliers_options)
    site_multipliers = [0] * len(sites_entropies_np)
    for i in range(len(data_splitted)):
        for j in range(len(data_splitted[i])):
            site_multipliers[data_splitted[i][j]] = site_multipliers_options[i]

    print(site_multipliers_options)
    #print(site_multipliers)
    print(set(site_multipliers))
    site_multipliers_np =  np.asarray(site_multipliers)
    plt.hist(site_multipliers_np)
    plt.show()

    file_of_rates = open("site_rates/total_rates/sites_rates_of_total_" + str(m) + "_levels_querires_cluster"+str(cluster_num)+".csv", "w")
    file_of_rates.write("site,rate\n")
    for i, rate in enumerate(site_multipliers):
        file_of_rates.write(str(i) + "," + str(rate) + "\n")
    file_of_rates.close()


def total_site_rate_calculator_m(m=5):
    a, b, dic_genes = parseFasta()
    site_enropies = computeEntropies(list(dic_genes.values()), 0)
    sites_entropies_np = np.asarray(site_enropies)
    number_of_parts = m
    percentiles = []
    for i in range(1, m):
        p = np.percentile(sites_entropies_np, 100. * float(i) / number_of_parts, interpolation="nearest")
        percentiles.append(p)
    percentiles.append(max(sites_entropies_np))
    print(len(sites_entropies_np))
    data_splitted = []
    tocuhed = [False] * len(sites_entropies_np)
    for value_of_perc in percentiles:
        current_perc = []
        for ind, val in enumerate(sites_entropies_np):
            if val <= value_of_perc and not tocuhed[ind]:
                current_perc.append(ind)
                tocuhed[ind] = True
        data_splitted.append(current_perc)
    site_multipliers_options = np.arange(0.5, 1.55, 1.0 / (m - 1))
    site_multipliers = [0] * len(sites_entropies_np)
    for i in range(len(data_splitted)):
        for j in range(len(data_splitted[i])):
            site_multipliers[data_splitted[i][j]] = site_multipliers_options[i]

    print(site_multipliers_options)
    print(site_multipliers)
    file_of_rates = open("site_rates/total_rates/sites_rates_of_total_" + str(m) + "_levels.csv", "w")
    file_of_rates.write("site,rate\n")
    for i, rate in enumerate(site_multipliers):
        file_of_rates.write(str(i) + "," + str(rate) + "\n")
    file_of_rates.close()



def hamming_distance_with_site_multipliers(g1, g2, site_multipliers, multipliers):
    values_list=set(sorted(multipliers))
    key_rates={}
    #print(multipliers[0])
    #print(multipliers[10])


    for key, val in zip(values_list, site_multipliers):
        key_rates[key] = val
    for i, val in enumerate(multipliers):
        if val in key_rates.keys():
            multipliers[i] = key_rates[val]
    mid = int((len(multipliers)+1)/2)
    multipliers[mid]=1
    #print(multipliers[0])
    #print(multipliers[10])
    #exit(0)
    if len(g1) != len(g2):
        print("wrong lengths")
        return 0
    g1=g1.upper()
    g2=g2.upper()
    dis=0.
    count=0.
    for i in range(len(g1)):
        if (g1[i])!='-' and (g2[i])!='-':
            count=count+1
            if g1[i] != g2[i]:
                dis = dis + 1
                #if upperCase(g1[i]) != upperCase(g2[i]):
    dis_with_gaps=0.
    count_with_gaps=0.
    for i in range(len(g1)):
        count_with_gaps=count_with_gaps+1
        if g1[i] != g2[i]:
            dis_with_gaps = dis_with_gaps + 1 * multipliers[i]

    if count>0 and count_with_gaps>0:
        return dis_with_gaps/count_with_gaps#, count_with_gaps
    else:
        return dis_with_gaps #, count_with_gaps
    return min(site_multipliers)

def target_function_L1(sitemultipliers):
    tree_dis_file = "trees_dis_pairs.csv"
    df_multipliers = pd.read_csv("site_rates/total_rates/sites_rates_of_total_5_levels.csv")
    multipliers = df_multipliers['rate'].to_list()
    #print(len(multipliers))
    #query_species_file = "/Users/omersabary/Desktop/GenePlacement/3000/queries.txt"
    #query_s = []
    #query_species_f = open(query_species_file, "r")
    #for line in query_species_f:
        #query_s.append(line.rstrip())
    #dic_genes = concatination_list_of_gene(query_species_list=query_s)
    a, b, dic_genes = parseFasta()
    #print(tree_dis_file)
    #print(sitemultipliers)
    df_trees=pd.read_csv(tree_dis_file)
    pairs=df_trees['pair'].to_list()
    sum = 0
    #print(dic_genes.keys())
    for i, pair in enumerate(pairs):
        #print(100*i/len(pairs))
        spec1=pair[2:12]
        spec2=pair[15:-2]
        #for c in dic_genes[spec1]:
        #    if c!='-':
        #        print(c)

        query=df_trees[df_trees['pair']==pair]
        trees_distance = (float(query['tree_dis']))
        #print("Omer")
        sum = sum + abs(jci(trees_distance) - hamming_distance_with_site_multipliers(dic_genes[spec1], dic_genes[spec2], sitemultipliers, multipliers))
    return sum


#target_function_L1(sitemultipliers=)
#target_function_L1()

x0 = [1, 2, 3, 4, 5]
cons = ({'type': 'ineq', 'fun': lambda x:  x[0] - 2 * x[1] + 2},
        {'type': 'ineq', 'fun': lambda x: -x[0] - 2 * x[1] + 6},
        {'type': 'ineq', 'fun': lambda x: -x[0] + 2 * x[1] + 2})

'''
query_species_file = "/Users/omersabary/Desktop/GenePlacement/3000/refs.txt"
query_s = []
query_species_f = open(query_species_file, "r")
for line in query_species_f:
    query_s.append(line.rstrip())
'''
total_site_rate_calculator_m_queries(m=7, cluster_num=118)
exit(0)


multipliers_values = [0.5, 0.75, 1, 1.25, 1.5]
#print(multipliers)
#exit(0)
#a, b, dic_genes = parseFasta()
#for al in a:
#    print(len(al))
#print(b)
#exit()

#dic_genes = concatination_list_of_gene(query_species_list=query_s)
#site_enropies=computeEntropies(list(dic_genes.values()),  0)
#print(b)
#exit(0)
res = minimize(fun=target_function_L1, x0=multipliers_values, method='Nelder-Mead')
print(res.x)
#site_rate_calculator_m_levels()

