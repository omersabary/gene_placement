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
import random

def jci(dis):
    return .75 * (1 - math.exp((-4.0* dis)/ 3 ))

def generatelendf():
    #df = pd.read_csv("pairs_distances.csv")
    df = pd.read_csv("pairs_distances_new_clusters.csv")
    #df_with_gaps = pd.read_csv("pairs_distances_with_gaps_new.csv")
    df_with_gaps = pd.read_csv("pairs_distances_new_clusters.csv")

    df=df.mask(df==-1)
    df_with_gaps=df_with_gaps.mask(df==-1)
    len_column=[]
    len_column_with_gaps=[]
    for col in df.columns:
        if col.endswith('len'):
            len_column.append(col)
    for col in df_with_gaps.columns:
        if col.endswith('len'):
            len_column_with_gaps.append(col)
    df['total_len'] = df[len_column].sum(axis=1)
    df_with_gaps['total_len_with_gaps'] = df_with_gaps[len_column_with_gaps].sum(axis=1)
    coluns = ['pair', 'total_len']
    df =(df[coluns])
    coluns = ['pair', 'total_len_with_gaps']
    df_with_gaps = (df_with_gaps[coluns])
    df_with_gaps.to_csv("pairs_len_with_gaps.csv")
    df.to_csv("pairs_len_new_clusters.csv")
    exit(0)

def set_multipliers(sitemultipliers, multipliers, m):
    print(sitemultipliers)
    #sitemultipliers = sorted(sitemultipliers, reverse=True)
    #sitemultipliers = list(sitemultipliers)
    print(sitemultipliers)
    set_mult = sorted(list(set(multipliers)))
    print(set_mult)
    for i in range(len(set_mult)):
        multipliers[multipliers == set_mult[i] ] = sitemultipliers[i]
    return multipliers

def plot_tree_vs_hm():
    hm_matrix=np.load("test_np_file.npy")
    rate_multipliers_files="site_rates/total_rates/sites_rates_of_total_5_levels.csv"
    tree_dis_file = "trees_dis_pairs.csv"
    df_lens = pd.read_csv("pairs_len.csv")
    lens = df_lens['total_len'].to_list()
    df_trees = pd.read_csv(tree_dis_file)
    df_multipliers = pd.read_csv(rate_multipliers_files)
    multipliers = df_multipliers['rate'].to_numpy()

    multipliers = np.transpose(multipliers)
    #print(hm_matrix.shape)
    print_mul = np.matmul(hm_matrix, multipliers)
    ones = [1 for i in range(multipliers.shape[0])]
    non_normalized = np.matmul(hm_matrix, ones)
    site_multipliers= [1.00000000e+00, 7.01417630e-07, 1.00000003e+00, 1.00000000e+00, 5.69963860e-05] #inital [1, 1, 1, 1, 1]
    site_multipliers =[1.00000000e+00, 1.00000000e+00, 1.00000001e+00, 1.20045261e-05, 4.06999049e-05] #inital [.6, .9, 1, 1.3, 1.6]
    site_multipliers = [9.99999999e-01, 2.15341300e-05, 1.00000000e+00, 1.00000003e+00, 5.55985583e-05]
    site_multipliers= [0.99865852, 0.99964191, 1,         1.459614,  1.45983787] # [.5, .75, 1, 1.25, 1.5]
    site_multipliers = [0.99848304, 0.99947748, 1.        , 1.45975471, 1.45976027] #[.98, .99, 1, 1.01, 1.02]
    site_multipliers_L2 = [0.99999894, 0.99999972, 1.        , 1.28762778, 1.28762797] #[.98, .99, 1, 1.01, 1.02] L2
    #site_multipliers_L2 = [0.99999881, 0.99999969, 1.        , 1.28762785, 1.28762786]  #[.5, .75, 1, 1.25, 1.5] L2

    site_multipliers= set_multipliers(site_multipliers, df_multipliers['rate'].to_numpy(), 0)
    print_mul_opt = np.matmul(hm_matrix, site_multipliers)

    site_multipliers_L2= set_multipliers(site_multipliers_L2, df_multipliers['rate'].to_numpy(), 0)
    print_mul_opt_L2 = np.matmul(hm_matrix, site_multipliers)

    pairs = df_trees['pair'].to_list()
    selected_pairs=random.sample(range(len(pairs)), 10)
    selected_pairs=sorted(selected_pairs)
    tree_dist=df_trees['tree_dis'].to_list()
    #tree_dist = [(tree_dist[i]) for i in selected_pairs]
    #plt.plot(tree_dist, marker='o')
    #tree_dist = [jci(val) for val in tree_dist]
    #plt.plot(tree_dist, marker='x')
    tree_dist = [jci(tree_dist[i]) for i in selected_pairs]
    plt.plot(tree_dist, marker='o', label ='Tree Dis')
    print_mul = [print_mul[i]/lens[i] for i in selected_pairs]
    print_mul_opt = [print_mul_opt[i]/lens[i] for i in selected_pairs]
    print_mul_opt_L2 = [print_mul_opt_L2[i]/lens[i] for i in selected_pairs]

    plt.scatter(x=range(len(selected_pairs)), y=print_mul, marker='x', color='blue', label ='Naive Multipliers')
    plt.scatter(x=range(len(selected_pairs)), y=print_mul_opt, marker='*', color='orange', label ='Optimal Multipliers - L1')
    plt.scatter(x=range(len(selected_pairs)), y=print_mul_opt_L2, marker='s', color='violet', label ='Optimal Multipliers - L2')
    #plt.xticks(ticks=range(len(selected_pairs)), labels=[pairs[i] for i in selected_pairs], rotation=90)
    plt.xticks(ticks=range(len(selected_pairs)), labels= selected_pairs)

    non_normalized = [non_normalized[i]/lens[i] for i in selected_pairs]
    plt.scatter(x=range(len(selected_pairs)), y=non_normalized, marker='o', color='red', label ='Without Multipliers')
    plt.grid()
    plt.legend()
    plt.show()

def plot_tree_vs_hm7():
    hm_matrix=np.load("test_np_file.npy")
    rate_multipliers_files="site_rates/total_rates/sites_rates_of_total_7_levels.csv"
    tree_dis_file = "trees_dis_pairs.csv"
    df_lens = pd.read_csv("pairs_len.csv")
    lens = df_lens['total_len'].to_list()
    df_trees = pd.read_csv(tree_dis_file)
    df_multipliers = pd.read_csv(rate_multipliers_files)
    multipliers = df_multipliers['rate'].to_numpy()

    multipliers = np.transpose(multipliers)
    #print(hm_matrix.shape)
    print_mul = np.matmul(hm_matrix, multipliers)
    ones = [1 for i in range(multipliers.shape[0])]
    non_normalized = np.matmul(hm_matrix, ones)
    site_multipliers= [1.00000000e+00, 7.01417630e-07, 1.00000003e+00, 1.00000000e+00, 5.69963860e-05] #inital [1, 1, 1, 1, 1]
    site_multipliers =[1.00000000e+00, 1.00000000e+00, 1.00000001e+00, 1.20045261e-05, 4.06999049e-05] #inital [.6, .9, 1, 1.3, 1.6]
    site_multipliers = [9.99999999e-01, 2.15341300e-05, 1.00000000e+00, 1.00000003e+00, 5.55985583e-05]
    site_multipliers= [0.99865852, 0.99964191, 1,         1.459614,  1.45983787] # [.5, .75, 1, 1.25, 1.5]
    site_multipliers = [0.99315445, 0.99830454, 0.99950888, 1.        , 1.44353745, 1.44440051, 1.44486912] #[.5, .75, 1, 1.25, 1.5]
    site_multipliers_L2 = [0.99979432, 0.99995071, 0.99998412, 1.        , 1.27741436, 1.27744719, 1.27746513] #Naive L2
    #site_multipliers_L2 = [1.00000001, 1.00000001, 1.00000021, 1.5       , 1.50000002,1.50000011, 1.50000948] #Naive L2 - mean 1.5
    #site_multipliers_L2 = [0.99999881, 0.99999969, 1.        , 1.28762785, 1.28762786]  #[.5, .75, 1, 1.25, 1.5] L2

    site_multipliers= set_multipliers(site_multipliers, df_multipliers['rate'].to_numpy(), 0)
    print_mul_opt = np.matmul(hm_matrix, site_multipliers)

    site_multipliers_L2= set_multipliers(site_multipliers_L2, df_multipliers['rate'].to_numpy(), 0)
    print_mul_opt_L2 = np.matmul(hm_matrix, site_multipliers)

    pairs = df_trees['pair'].to_list()
    selected_pairs=random.sample(range(len(pairs)), 10)
    selected_pairs=sorted(selected_pairs)
    tree_dist=df_trees['tree_dis'].to_list()
    #tree_dist = [(tree_dist[i]) for i in selected_pairs]
    #plt.plot(tree_dist, marker='o')
    #tree_dist = [jci(val) for val in tree_dist]
    #plt.plot(tree_dist, marker='x')
    tree_dist = [jci(tree_dist[i]) for i in selected_pairs]
    plt.plot(tree_dist, marker='o', label ='Tree Dis')
    print_mul = [print_mul[i]/lens[i] for i in selected_pairs]
    print_mul_opt = [print_mul_opt[i]/lens[i] for i in selected_pairs]
    print_mul_opt_L2 = [print_mul_opt_L2[i]/lens[i] for i in selected_pairs]

    plt.scatter(x=range(len(selected_pairs)), y=print_mul, marker='x', color='blue', label ='Naive Multipliers')
    plt.scatter(x=range(len(selected_pairs)), y=print_mul_opt, marker='*', color='orange', label ='Optimal Multipliers - L1')
    plt.scatter(x=range(len(selected_pairs)), y=print_mul_opt_L2, marker='s', color='violet', label ='Optimal Multipliers - L2')
    #plt.xticks(ticks=range(len(selected_pairs)), labels=[pairs[i] for i in selected_pairs], rotation=90)
    plt.xticks(ticks=range(len(selected_pairs)), labels= selected_pairs)

    non_normalized = [non_normalized[i]/lens[i] for i in selected_pairs]
    plt.scatter(x=range(len(selected_pairs)), y=non_normalized, marker='o', color='red', label ='Without Multipliers')
    plt.grid()
    plt.legend()
    plt.show()

def plot_tree_vs_hm3():
    hm_matrix=np.load("test_np_file.npy")
    rate_multipliers_files="site_rates/total_rates/sites_rates_of_total_3_levels.csv"
    tree_dis_file = "trees_dis_pairs.csv"
    df_lens = pd.read_csv("pairs_len.csv")
    lens = df_lens['total_len'].to_list()
    df_trees = pd.read_csv(tree_dis_file)
    df_multipliers = pd.read_csv(rate_multipliers_files)
    multipliers = df_multipliers['rate'].to_numpy()

    multipliers = np.transpose(multipliers)
    #print(hm_matrix.shape)
    print_mul = np.matmul(hm_matrix, multipliers)
    ones = [1 for i in range(multipliers.shape[0])]
    non_normalized = np.matmul(hm_matrix, ones)
    site_multipliers= [1, 1, 1.5] #inital [1, 1, 1, 1, 1]
    site_multipliers_L2 = [1, 1, 1.32]
    site_multipliers= set_multipliers(site_multipliers, df_multipliers['rate'].to_numpy(), 0)
    print_mul_opt = np.matmul(hm_matrix, site_multipliers)

    site_multipliers_L2= set_multipliers(site_multipliers_L2, df_multipliers['rate'].to_numpy(), 0)
    print_mul_opt_L2 = np.matmul(hm_matrix, site_multipliers)

    pairs = df_trees['pair'].to_list()
    selected_pairs=random.sample(range(len(pairs)), 10)
    selected_pairs=sorted(selected_pairs)
    tree_dist=df_trees['tree_dis'].to_list()
    #tree_dist = [(tree_dist[i]) for i in selected_pairs]
    #plt.plot(tree_dist, marker='o')
    #tree_dist = [jci(val) for val in tree_dist]
    #plt.plot(tree_dist, marker='x')
    tree_dist = [jci(tree_dist[i]) for i in selected_pairs]
    plt.plot(tree_dist, marker='o', label ='Tree Dis')
    print_mul = [print_mul[i]/lens[i] for i in selected_pairs]
    print_mul_opt = [print_mul_opt[i]/lens[i] for i in selected_pairs]
    print_mul_opt_L2 = [print_mul_opt_L2[i]/lens[i] for i in selected_pairs]

    plt.scatter(x=range(len(selected_pairs)), y=print_mul, marker='x', color='blue', label ='Naive Multipliers')
    plt.scatter(x=range(len(selected_pairs)), y=print_mul_opt, marker='*', color='orange', label ='Optimal Multipliers - L1')
    plt.scatter(x=range(len(selected_pairs)), y=print_mul_opt_L2, marker='s', color='violet', label ='Optimal Multipliers - L2')
    #plt.xticks(ticks=range(len(selected_pairs)), labels=[pairs[i] for i in selected_pairs], rotation=90)
    plt.xticks(ticks=range(len(selected_pairs)), labels= selected_pairs)

    non_normalized = [non_normalized[i]/lens[i] for i in selected_pairs]
    plt.scatter(x=range(len(selected_pairs)), y=non_normalized, marker='o', color='red', label ='Without Multipliers')
    plt.grid()
    plt.legend()
    plt.show()


def plot_tree_vs_hm_new(cluster_num="", site_multipliers_L2=[]):
    hm_matrix = np.load("cluster"+cluster_num+"_hamming.npy")
    rate_multipliers_files = "site_rates/total_rates/sites_rates_of_total_7_levels_querires_cluster"+cluster_num+".csv"
    tree_dis_file = "tree_distance_pairs_new_clusters_cluster"+cluster_num+".csv"
    #df_lens = pd.read_csv("pairs_len.csv")
    lens = np.load("cluster"+cluster_num+"_lens.npy")
    df_trees = pd.read_csv(tree_dis_file)

    pairs = df_trees['pair'].to_list()
    selected_pairs = random.sample(range(len(pairs)), 50)
    selected_pairs = sorted(selected_pairs)
    hm_matrix = hm_matrix[selected_pairs, :]
    selected_pairs = range(len(selected_pairs))
    df_multipliers = pd.read_csv(rate_multipliers_files)


    multipliers_orig = df_multipliers['rate'].to_numpy()
    multipliers_orig = np.transpose(multipliers_orig)
    multipliers = np.copy(multipliers_orig)
    print(set(multipliers))
    print("set of multipliers is")
    multipliers = np.transpose(multipliers)
    # print(hm_matrix.shape)
    print_mul_no_reverse = np.matmul(hm_matrix, multipliers)
    site_multipliers_L2 = set_multipliers(site_multipliers_L2, np.transpose(df_multipliers['rate'].to_numpy()), 0)
    set_multi = sorted(list(set(multipliers)), reverse=True)

    multipliers = set_multipliers(set_multi, multipliers, 0 )
    print_mul = np.matmul(hm_matrix, multipliers)
    ones = [1 for i in range(multipliers.shape[0])]
    non_normalized = np.matmul(hm_matrix, ones)
    print_mul_opt_L2 = np.matmul(hm_matrix, site_multipliers_L2)


    tree_dist = df_trees['tree_dis'].to_list()
    # tree_dist = [(tree_dist[i]) for i in selected_pairs]
    # plt.plot(tree_dist, marker='o')
    # tree_dist = [jci(val) for val in tree_dist]
    # plt.plot(tree_dist, marker='x')
    tree_dist = [jci(tree_dist[i]) for i in selected_pairs]
    plt.plot(tree_dist, marker='o', label='Tree Dis')
    print_mul = [print_mul[i] / lens[i] for i in selected_pairs]
    print_mul_no_reverse = [print_mul_no_reverse[i]/lens[i] for i in selected_pairs]
    #print_mul_opt = [print_mul_opt[i] / lens[i] for i in selected_pairs]
    print_mul_opt_L2 = [print_mul_opt_L2[i] / lens[i] for i in selected_pairs]

    plt.scatter(x=range(len(selected_pairs)), y=print_mul, marker='x', color='blue', label='Naive Multipliers')
    #plt.scatter(x=range(len(selected_pairs)), y=print_mul_opt, marker='*', color='orange',
    #            label='Optimal Multipliers - L1')
    plt.scatter(x=range(len(selected_pairs)), y=print_mul_opt_L2, marker='s', color='violet',
                label='Optimal Multipliers - L2')
    plt.scatter(x=range(len(selected_pairs)), y=print_mul_no_reverse, marker='x', color='brown',
                label='Naive Multipliers - No Reverse')
    # plt.xticks(ticks=range(len(selected_pairs)), labels=[pairs[i] for i in selected_pairs], rotation=90)
    plt.xticks(ticks=range(len(selected_pairs)), labels=selected_pairs)

    non_normalized = [non_normalized[i] / lens[i] for i in selected_pairs]
    plt.scatter(x=range(len(selected_pairs)), y=non_normalized, marker='o', color='red', label='Without Multipliers')
    plt.grid()
    plt.legend()
    plt.show()


def average_L2_dis(cluster_num="", site_multipliers_L2=[]):
    hm_matrix = np.load("cluster"+cluster_num+"_hamming.npy")
    rate_multipliers_files = "site_rates/total_rates/sites_rates_of_total_7_levels_querires_cluster"+cluster_num+".csv"
    tree_dis_file = "tree_distance_pairs_new_clusters_cluster"+cluster_num+".csv"
    #df_lens = pd.read_csv("pairs_len.csv")
    lens = np.load("cluster"+cluster_num+"_lens.npy")
    df_trees = pd.read_csv(tree_dis_file)


    df_multipliers = pd.read_csv(rate_multipliers_files)
    multipliers_orig = df_multipliers['rate'].to_numpy()
    multipliers_orig = np.transpose(multipliers_orig)
    multipliers = np.copy(multipliers_orig)
    print(set(multipliers))
    print("set of multipliers is")
    multipliers = np.transpose(multipliers)
    # print(hm_matrix.shape)
    print_mul_no_reverse = np.matmul(hm_matrix, multipliers)
    site_multipliers_L2 = set_multipliers(site_multipliers_L2, np.transpose(df_multipliers['rate'].to_numpy()), 0)
    print_mul_opt_L2 = np.matmul(hm_matrix, site_multipliers_L2)
    tree_dist = df_trees['tree_dis'].to_list()
    set_multi = sorted(list(set(multipliers)), reverse=True)
    multipliers = set_multipliers(set_multi, multipliers, 0 )
    print_mul = np.matmul(hm_matrix, multipliers)
    diff_naive= [(jci(tree_dist[i])-print_mul[i]/lens[i])**2 for i in range(len(print_mul))]
    #print(diff_naive)
    diff_naive_np=np.copy(diff_naive)
    diff_naive_np_sum=diff_naive_np.sum()/len(diff_naive_np)
    print("Naive")
    print(len(diff_naive_np))
    print(diff_naive_np_sum)

    diff_naive_nr= [(jci(tree_dist[i])-print_mul_no_reverse[i]/lens[i])**2 for i in range(len(print_mul_no_reverse))]
    #print(diff_naive)
    diff_naive_nr_np=np.copy(diff_naive_nr)
    diff_naive_nr_np_sum=diff_naive_nr_np.sum()/len(diff_naive_nr_np)
    print("Naive No Reverse")
    print(len(diff_naive_nr_np))
    print(diff_naive_nr_np_sum)



    diff_naive_L2= [(jci(tree_dist[i])-print_mul_opt_L2[i]/lens[i])**2 for i in range(len(print_mul_opt_L2))]
    diff_naive_L2_np = np.copy(diff_naive_L2)
    diff_naive_L2_np_sum=diff_naive_L2_np.sum()/len(diff_naive_L2_np)
    print("Optimal")
    print((len(diff_naive_L2_np)))
    print(diff_naive_L2_np_sum)


    ones = [1 for i in range(multipliers.shape[0])]
    mult_ones = np.matmul(hm_matrix, ones)
    diff_naive= [(jci(tree_dist[i])-mult_ones[i]/lens[i])**2 for i in range(len(mult_ones))]
    #print(diff_naive)
    diff_naive_np=np.copy(diff_naive)
    diff_naive_np_sum=diff_naive_np.sum()/len(diff_naive_np)
    print("Ones")
    print(diff_naive_np_sum)




#plot_tree_vs_hm_new(cluster_num="118", site_multipliers_L2= [0.64056631, 0.99112782, 1.01509641, 1.00057839, 0.6179216 ,
#       0.97855129, 0.00750239])
#average_L2_dis(cluster_num="118", site_multipliers_L2= [2.087847, 1.54307874, 0.553876, 0.8157747, 0.77180169, 1.2321442, 1.01655312])
average_L2_dis(cluster_num="118", site_multipliers_L2= [0.48392136, 0.7829391, 1.10477099, 1.27281974, 0.94540131, 1.0422456,0.63708118])


#exit(0)
#plot_tree_vs_hm_new(cluster_num="118", site_multipliers_L2= [2.087847, 1.54307874, 0.553876, 0.8157747, 0.77180169, 1.2321442, 1.01655312])
#[0.52886575, 0.86250947, 0.88112116, 0.99474585, 0.99600333,
#       1.00024134, 1.01509641]

#plot_tree_vs_hm_new(cluster_num="118", site_multipliers_L2= [0.48392136, 0.7829391, 1.10477099, 1.27281974, 0.94540131, 1.0422456,0.63708118])
# 0.63708118])
#0.48392136 0.7829391  1.10477099 1.27281974 0.94540131 1.04224565
# 0.63708118
#plot_tree_vs_hm_new(cluster_num="118", site_multipliers_L2= [1.56410309, 1.32228623, 1.19225745, 1.08637933, 1.02606133, 1.05657028, 1.04465448])

#generatelendf()
#plot_tree_vs_hm3()

