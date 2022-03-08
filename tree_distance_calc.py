
import random
import matplotlib.pyplot as plt
import os
import numpy as np
from ratesCalculator import get_diameters
from ratesCalculator import get_species_num
from ratesCalculator import get_species_branches
from ratesCalculator import reCalculate
import json

import pandas

import math
import subprocess as sp



def jc(dis):
    return (-0.75)*math.log(1-dis*4.0/3.0)

def upperCase(c):
    dic={'a': 'A', 'A':'A', 'C':'C', 'c':'C', 'g':'G', 'G':'G',  't':'T', 'T':'T'}
    return dic[c]

def isValidChar(c):
    if c not in {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'}:
        return False
    return True


def hamming_distance(g1, g2):
    if len(g1) != len(g2):
        print("wrong lengths")
        return 0
    dis=0
    for i in range(len(g1)):
        if isValidChar(g1[i]) and isValidChar(g2[i]):
            if upperCase(g1[i]) != upperCase(g2[i]):
                dis=dis+1
    return dis



def hamming_distance_norm(g1, g2):
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

    if count>0:
        #print(dis)
        #print(count)
        return dis/count, count
    else:
        return dis, count

def hammingDistanceDict(dicSpecies):
    dic_res={}
    for key1 in dicSpecies.keys():
        for key2 in dicSpecies.keys():
            if key1==key2:
                continue
            else:
                dis=hamming_distance_norm(dicSpecies[key1], dicSpecies[key2])
                dic_res[str(key1) + " " +str(key2)]=(-.75)*math.log((1-(4.0/3.0)*(dis)))

    print(dic_res)
    return dic_res

def get_species_branches(directory = "/Users/omersabary/Desktop/GenePlacement/3000/trees/genes/rooted/"):
    branches_lengths=[]
    filenames=[]
    for filename in sorted(os.listdir(directory)):
        if filename.endswith(".nwk"):
            filepath = os.path.join(directory, filename)
            filenames.append(filename[:-4])
            f=(sp.getoutput('nw_distance -s l -m p ' + filepath)).split()
            f=[float(len) for len in f]
            #print(f)
            branches_lengths.append(f)
            #print(branches_lengths)
            #exit(0)
            continue
        else:
            continue
    return branches_lengths, filenames

#reCalculate()
#exit(0)
## X - Axis
cluster_num="29"
query_species_file="/Users/omersabary/Desktop/GenePlacement/3000/omer100.txt"
query_species_file="/Users/omersabary/Desktop/GenePlacement/clusters_omer/cluster"+cluster_num+".txt"
query_s=[]
query_species_f=open(query_species_file, "r")
for line in query_species_f:
    query_s.append(line.rstrip()[:10])
print(query_s)

pairs_species=[]
for i in range(len(query_s)):
    for j in range(i+1, len(query_s)):
        if i != j:
            pair = []
            pair.append(query_s[i][:10])
            pair.append(query_s[j][:10])
            pairs_species.append(pair)
    #r = random.sample(set(query_s), 2)
    #pairs_species.append(r)
print("pairs_species")
print(len(pairs_species))


'''
distances=open("/Users/omersabary/Desktop/GenePlacement/3000/trees/distances-omer.txt", "r")
dic_distances={}
for line in distances:
    line_splitted=line.split()
    #print(line_splitted)
    dic_distances[str(line_splitted[0])]=float(line_splitted[1])

pair_dis_tree= []
for pair in pairs_species:
    spec1, spec2 = pair
    dis = abs(dic_distances[spec1] - dic_distances[spec2])
    pair_dis_tree.append(dis)
'''
pair_dis_tree= []
pair_dic_dis={}
tree_path="/Users/omersabary/Desktop/GenePlacement/3000/trees/"
tree_path="/Users/omersabary/Desktop/GenePlacement/"

for pair in pairs_species:
    if len(pair)==2:
        filepath = os.path.join(tree_path, "astral_mp.nwk")
        filepath = os.path.join(tree_path, "backbone.nwk")
        filepath = os.path.join(tree_path, "astral.rand.lpp.nwk")
        #os.system("nw_distance -mm "+filepath+" "+pair[0]+" "+pair[1])
        m=sp.getoutput("nw_distance -sa -mm "+filepath+" "+pair[0]+" "+pair[1])
        dis=float(m[2:].split()[0])
        pair_dis_tree.append(dis)
        pair_dic_dis[str(pair)] = dis

with open('tree_distance_pairs_new_clusters_cluster'+cluster_num+'.csv', 'w') as outfile:
    outfile.write('pair,tree_dis\n')
    for key in pair_dic_dis.keys():
        outfile.write(key.replace(', ','_')+","+str(pair_dic_dis[key])+"\n")


exit(0)
## Y - Axis

## Gene - Alignment Distance - No Rates

def createDictionaryAlignment(file_path, cluster):
    dic={}
    file=open(file_path, "r")
    line = file.readline()
    while True:
        if line[1:-1] in cluster:
            spec=line[1:-1]
            align_gene=""
            line=file.readline()
            while not line.startswith(">"):
                align_gene=align_gene+line.replace('\n', '')
                line = file.readline()
            dic[spec]=align_gene
        else:
            line=file.readline()
        if line=="":
            break
    #print(len(dic.keys()))
    file.close()
    return dic

'''
branches_diameter, filenames_d = get_diameters(directory="/Users/omersabary/Desktop/GenePlacement/3000/trees/genes/rooted/")
diameters=[float(b)/np.max(branches_diameter) for b in branches_diameter]
dic_rates_diameter={}
for ind, name in enumerate(filenames_d):
    dic_rates_diameter[name]=diameters[ind]
print(dic_rates_diameter)
#exit(0)
'''

#plt.plot(pair_dis_tree, pair_dis_tree, marker='o')
#plt.xticks(pair_dis_tree, pairs_species, rotation='90', size=8)
alignment_dir ="/Users/omersabary/Desktop/GenePlacement/3000/alignments/ffn/"
dic_res={}


genes_hamming_dis=open("genes_hamming_distances_new_clusters.csv", "w")
gene_list=[]
for filename in os.listdir(alignment_dir):
    if filename.endswith(".fasta"):
        gene_list.append(filename[:-6])
line= "pair,"
for gene_name in gene_list:
    line = line + gene_name + "_dis," + gene_name+"_len,"
line = line[:-1] + "\n"
genes_hamming_dis.write(line)
#print(line)
#exit(0)



for i in range(0, int(len(pairs_species))):
    genes_distance_uncorrected=[]
    gene_distance_corrected=[]
    hamming_distances=[]
    lengths=[]
    y_sum=0.0
    len_sum=0
    #y2_sum=0.0
    #len2_sum=0
    hamming_dis_sum=0.0
    for filename in os.listdir(alignment_dir):
        if filename.endswith(".fasta"):
            filepath = os.path.join(alignment_dir, filename)
            open(filepath, "r")
            filenameshort=filename[:-6]
            dic_gene = createDictionaryAlignment(filepath, pairs_species[i])


            #exit(0)
            if int(len(dic_gene)) == 2:
                dic_res[filename[:-6]] = dic_gene
                y, gene_len = hamming_distance_norm(dic_gene[pairs_species[i][0]], dic_gene[pairs_species[i][1]])
                hamming_dis_sum=hamming_dis_sum+y*gene_len
                #print(str(y)+","+str(gene_len))
                hamming_distances.append(y)
                lengths.append(gene_len)
            else:
                hamming_distances.append(-1)
                lengths.append(-1)
                #exit(0)
                '''
                if y>=.75:
                    #print(y)
                    y=.749999
                hamming_dis = y # 1 - (4.0 / 3.0) * (y)
                #print(hamming_dis)

                if hamming_dis!=0:
                    y=jc(hamming_dis)
                    y_sum=y_sum+y*gene_len
                    len_sum=len_sum+gene_len
                    genes_distance_uncorrected.append(y)
                    y2= y*(1.0/(dic_rates[filenameshort]))
                    y2_sum=y2_sum+y2*gene_len
                    len2_sum=len2_sum+gene_len
                    gene_distance_corrected.append(y2)
                    y3=y*(1.0/(dic_rates_diameter[filenameshort]))
                else:
                    y = jc(0.000000001)
                    y_sum = y_sum + y * gene_len
                    len_sum = len_sum + gene_len
                    y2 = y*(1.0/(dic_rates[filenameshort]))
                    y2_sum = y2_sum + y2 * gene_len
                    len2_sum = len2_sum + gene_len
                    y3=y*(1.0/(dic_rates_diameter[filenameshort]))
                rgb = np.random.rand(3, )
                '''
                #plt.scatter(pair_dis_tree[i], y, marker='o', color=[rgb], label=filename[:-6])
                #plt.scatter(pair_dis_tree[i], y2, marker='*', color=[rgb], label=filename[:-6])
                #plt.scatter(pair_dis_tree[i], y3, marker='x', color=[rgb], label=filename[:-6])
                #plt.legend()
            #plt.boxplot(pair_dis_tree[i]-0.05, gene_distance_corrected)
            #plt.boxplot(pair_dis_tree[i]+0.05, genes_distance_uncorrected)
            continue
        else:
            continue
    line=str(pairs_species[i]).replace(", ",'_')+","
    for dis, length in zip(hamming_distances, lengths):
        line=line+(str(dis)+","+str(length)+",")
    line = line[:-1]
    line = line+("\n")
    print(line)
    genes_hamming_dis.write(line)
    #exit(0)
    #if len_sum != 0:
    #    plt.scatter(pair_dis_tree[i], y_sum / len_sum, marker='o', color=[rgb], label=filename[:-6])
    #    plt.scatter(pair_dis_tree[i], y2_sum / len2_sum, marker='*', color=[rgb], label=filename[:-6])
    #    plt.scatter(pair_dis_tree[i], jc(hamming_dis_sum/len2_sum), marker='x', color=[rgb], label=filename[:-6])


#print(dic_res)


#print(pair_dis_tree)

#exit(0)
'''
gene_tree_dir ="/Users/omersabary/Desktop/GenePlacement/3000/trees/genes/"
for x_ind, pair in enumerate(pairs_species):
    #x_ind=pairs_species.find(pair)
    x=pair_dis_tree[x_ind]
    for filename in os.listdir(gene_tree_dir)[0:10]:
        #file=open(filename, "r")
        m = sp.getoutput("nw_distance -sa -mm " + gene_tree_dir+ filename + " " + pair[0] + " " + pair[1])
        #print(m)
        #print(m[0])
        if m[0] != 'W':
            y=float(m[2:].split()[0])
            plt.scatter(x, y, marker='X', label=filename[:-3])
            y2 = (dic_rates[filename[:-4]])*y
            plt.scatter(x, y2, marker='+', label=filename[:-3])


'''
'''
plt.grid()
#plt.ylim(0,2)
plt.xlabel("Tree Distance")
plt.ylabel("Conc. Dist. - Naive $o$ vs Rate $*$")
plt.show()



'''