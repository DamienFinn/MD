# Metagenomic diversity (MD) calculator 
# Measures size and (dis)similarity of protein clusters from 
# MMSeqs2 cluster and alignment outputs

# Note: this code assumes you have installed MMseqs2 

# Usage without random subsampling: python MD.py -i <predicted_protein-encoding_genes_file.faa> -r F 
# With random subsampling: python MD.py -i <predicted_protein-encoding_genes_file.faa> -r T -N xxx
# where xxx is the total number of observations to randomly subsample
# Where users may wish to define their own MMSeqs2 E value and/or sequence identity cut-offs: python MD.py -i <predicted_protein-encoding_genes_file.faa> -r F -e xxx -I yyy
# where xxx and yyy are the user degined inputs for these cut-offs

# Written by Damien Finn, damien.finn@thuenen.de  

import os
import argparse
import numpy as np
from collections import defaultdict
import multiprocessing as mp
import random

# Define input data

parser = argparse.ArgumentParser(description= 'Derive a Metagenomic Diversity index from clustered amino acid sequences')
parser.add_argument('-i', '--input-file', required = True, nargs = 1, 
                    help = 'Amino acid sequences as .faa')
parser.add_argument('-r', '--random-sampling', required = False, nargs = 1, 
                    help = 'To randomly subsample observations or not, as T if randomly subsampling')
parser.add_argument('-N', '--seq-obs', required = False, nargs = 1, type = int, 
                    help = 'If randomly subsampling, N denotes the number of observations to subsample without replacement')
parser.add_argument('-e', '--Eval', required = False, nargs = 1, type = str, 
                    help = 'Redefine the MMSeqs2 minimum Evalue cutoff score')
parser.add_argument('-I', '--identity', required = False, nargs = 1, type = str, 
                    help = 'Redefine the MMSeqs2 minimum sequence identity cutoff score')
args = parser.parse_args()

input_file = os.path.join(os.getcwd(), str(args.input_file[0]))

# Run MMseqs2 
os.system('mmseqs createdb ' + input_file + ' DB')

if args.Eval is None and args.identity is None: 
    os.system('mmseqs cluster DB Clu tmp')
elif args.Eval is not None and args.identity is None:
    os.system('mmseqs cluster DB Clu tmp -e ' + args.Eval[0] )
elif args.Eval is None and args.identity is not None:
    os.system('mmseqs cluster DB Clu tmp --min-seq-id ' + args.identity[0] )
elif args.Eval is not None and args.identity is not None:
    os.system('mmseqs cluster DB Clu tmp -e ' + args.Eval[0] + ' --min-seq-id ' + args.identity[0] )

os.system('mmseqs createtsv DB DB Clu Clu.tsv')
os.system('mmseqs search DB DB res tmp')
os.system('mmseqs convertalis DB DB res res.m8')

input_clust = os.path.join(os.getcwd(), str('Clu.tsv'))
input_aln = os.path.join(os.getcwd(), str('res.m8'))

# Generate a dictionary of clusters 
print(' -- Making Protein dictionary -- ')
clustlist = [] 

with open(input_clust) as file:
    for l in file:
        clustlist.append(l.rstrip())

if str(args.random_sampling[0])  == 'T':
    rarvals = random.sample(clustlist, int(args.seq_obs[0]))
    clustlist = rarvals

clusters = defaultdict(list) 

pool = mp.Pool(mp.cpu_count())

def make_clust_dict(clustlist):
    for l in clustlist:
        col1 = l.split('\t')[0]
        col2 = l.split('\t')[1]  
        clusters[col1].append(col2) 

pool.apply_async(make_clust_dict(clustlist))

# Generate a dictionary of pairwise distances
print(' -- Collating pair-wise dissimilarities -- ')
alnlist = [] 

with open(input_aln) as file:
    for l in file:
        alnlist.append(l.rstrip())

pairsim = defaultdict(list)

def make_pairsim_dict(alnlist):
    for l in alnlist:
        col1 = l.split('\t')[0]
        col2 = l.split('\t')[1]
        col3 = l.split('\t')[2]
        tmp = [col2, col3]
        pairsim[col1].append(tmp)

pool.apply_async(make_pairsim_dict(alnlist))

# Now to bring them together
print(' -- Bringing things together -- ')
clustdist = defaultdict(list)

def sorter(clusters, pairsim):
    for c in clusters.keys(): # key name
        v = clusters[c] # values in c
        for q in pairsim.keys():
            if c == q: # if clusters key and pairsim key match
                valq = pairsim[q] # pairwise values in q
                for a in v:
                    for b in valq:
                        if a == b[0]: # matching cluster pair with pairsim value pair
                            clustdist[c].append(1 - float(b[1]))

pool.apply_async(sorter(clusters, pairsim))

pool.close()
pool.join()

# Derive indices from the cluster dictionary

print(' -- Deriving diversity indices -- ')

# Total Contigs
O = len(clustlist)
print('Total protein-encoding genes: ', O)

# Protein Richness
P = len(clustdist)
print('Protein Richness: ', P)

CDlens = [ len(k) for k in clustdist.values()] 
tmp = np.sum(CDlens)

# Shannon Diversity
H = [ (len(k)/tmp)*np.log(len(k)/tmp) for k in clustdist.values() ]
print('Shannon Diversity: ', np.sum(H)*-1) 

# Simpson Diversity
S = [ np.square(len(k)/tmp) for k in clustdist.values() ]
print('Simpson Evenness: ', 1 - np.sum(S))

# Metagenomic diversity
MD = [ (1 + (np.sum(k))/len(k)) for k in clustdist.values()]
print('Log10 Protein Dissimilarity:', np.log10(np.sum(MD)))
print('Metagenomic Diversity index:', (1/O)*np.sum(MD))

# A little Aufraumen 
os.system('rm -rf tmp')
os.system('rm Clu*')
os.system('rm res*')
os.system('rm DB*')
