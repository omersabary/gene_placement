# gene_placement

**Gene Placement Project - UCSD 
**This repository summerizes the Gene Placement project I did for Prof. Siavash Mirarab at UCSD. 

Generally speaking, the goal was to approximate the distance between two species in a given phylogenetic tree. The approximation is based on their weighted hamming distance, where the weights are given to the different site. 
The site weight, also called site multipliers are computed based on the site entropies. 
Formal definitions can be found in the attached document and presentation. 


The python files are organized as follows. 
The species in the tree are assumed to be clusteres (with TreeCluster) by the average disance from a clade. Then, I organized a few clusters in a file, and picked species from one of these clusters (see omer-clusters.txt). 

ClusterPickleMaker.py - recieves the number of cluster, and create a pickled dictionary with the concatinated sequences of the species in the cluster, their length, 
and pickled hamming distance binary matrix that is of M X N dimension, where $N$ is the concatinated sequence length and $M$ is the number of pairs of species in the cluster. That is, for cluster of size $m$, $M = m \choose 2$. 


OptimizerNewVer.py - The optimizer the finds the optimal site multipliers. The optimizer receives the files that were generate by the ClusterPickleMaker.py, together with an array of the inital values and returns optimal site multipliers. 

SiteEntropies.py - Calculates the site entropies and the naive multipliers (see more details in the documents). 

Tree_distance_calc.py - utiliy function to calculate the tree distance based on the nw_distance function of nw_utils. 

PlotNewVer.py - visualiztion of the results. Picks 10-20 random pairs and presents the tree distance, the regular hamming distance and the wieght hamming distance and compare them. 


References: 
1. nw_utils: https://github.com/nwutils/nw-utils
2. TreeClluster: https://github.com/niemasd/TreeCluster
