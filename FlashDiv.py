# Taxonomic alpha-diversity from phyloFlash outputs

# usage without rarefying: python FlashDiv.py -i <NTUfull_abundance.csv> -r F
# with rarefying: python FlashDiv.py -i <NTUfull_abundance.csv> -r T -N xxx
# where xxx is the seq depth for rarefying 

# Written by Damien Finn, damien.finn@thuenen.de 

import os 
import argparse
import numpy as np
import random
from collections import defaultdict

# Define input data
parser = argparse.ArgumentParser(description = 'Derive taxonomic alpha-diversity indices from PhyloFlash 16S rRNA gene summaries')
parser.add_argument('-i', '--input-file', required = True, nargs = 1, 
                    help = 'List of taxa counts derived from PhyloFlash identification as NTUfull_abundance.csv file')
parser.add_argument('-r', '--rarefy', required = False, nargs = 1, 
                    help = 'To rarefy taxa observations or not, as T if rarefying')
parser.add_argument('-N', '--taxa-obs', required = False, nargs = 1, type = int, 
                    help = 'If rarefying, N denotes the number of observations')

args = parser.parse_args()

input_file = os.path.join(os.getcwd(), str(args.input_file[0]))

# Generate a list of observation counts
taxlist = []

with open(input_file) as file:
    for l in file:
        taxlist.append(l.rstrip())

obs = [] 

for l in taxlist:
    C = l.split(',')[1]
    obs.append(int(C)) 

# For rarefying, the observations must become singletons in a dictionary
tmplist = []  

countvec = 0

for l in obs:
    countvec += 1
    for i in range(l):
        tmplist.append('sp'+str(countvec))

# If rarefying
if str(args.rarefy[0]) == 'T':
    rarvals = random.sample(tmplist, int(args.taxa_obs[0]))  
    tmplist = rarvals

taxdict = defaultdict(list)

for i in tmplist:
    taxdict[i].append(tmplist.count(i)) 

# Deriving alpha-div indices from the taxa observations

# Total observations 
O = len(tmplist)
print('Total observations: ', int(O))

# Taxa Richness
R = len(taxdict)
print('Taxon Richness: ', R)

props = [ len(k) for k in taxdict.values()]
tmp = np.sum(props) 

# Shannon Diversity
H = [ (len(k)/tmp)*np.log(len(k)/tmp) for k in taxdict.values() ]
print('Shannon Diversity: ', np.sum(H)*-1)

# Simpson Diversity 
S = [ np.square(len(k)/tmp) for k in taxdict.values() ]
print('Simpson Evenness: ', 1 - np.sum(S)) 
