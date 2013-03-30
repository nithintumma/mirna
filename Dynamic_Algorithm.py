'''
Created on Jan 24, 2012
Aims to solve the maximum covering problem for miRNA's and Kegg Pathways using a dynamic progamming approach. Generally, the problem will be broken into sub-problems.
After finding the miRNA that targets the greatest number of genes in the pathway, the mirna is added to the list. The solution space is reduced to the genes that are not targetted by that mir, and the process is continued until the desired length is reached.
@author: nithintumma
'''
import random
import numpy
import csv
from Inhibitor_Modules import *

def getMirsPath (path):
    tmp_mirs = []
    tmp_mirs.extend(path_mirs[path])
    mir_list = []
    for i in tmp_mirs:
        if i in mirnas:
            mir_list.append(i)
    return mir_list


curr_path = "hsa04910"    # Just chose the Tryptophan metabolism, should be able to change this to any pathway.
# --- Open requisite files ---    
mirlist_file = open("/Users/nithintumma/Desktop/miRNAs.csv", "rU")
mirna_file = open("/Users/nithintumma/Desktop/miRNA_ordered.csv", "rU")
scores_file = open("/Users/nithintumma/Desktop/score_ordered.csv", "rU")
pathways_file = open("/Users/nithintumma/Desktop/Pathways_ordered.csv", "rU")
genes_file = open("/Users/nithintumma/Desktop/genes_ordered.csv", "rU")
paths_file = open("/Users/nithintumma/Desktop/pathways.csv")

pathways = list(csv.reader(paths_file))    # read unique pathways into list
pathways = [ x[0] for x in pathways ]

# --- READ miRNA MASTER LIST (only those in pathways) ---
mirnas = list(csv.reader(mirlist_file))    # read unique miRNAs into list
mirnas = [ x[0] for x in mirnas ]

# --- READ RELATIONS FILE INTO LIST ---
# "Pathway","Gene ID","Gene Name","Copies in Path","Transcript","miRNA","Score","Start","End"
#  0         1         2           3                4            5       6       7       8
relation = list(csv.reader(open("/Users/nithintumma/Desktop/relation.csv", "rU")))  # get relations
relation = relation[1:]                                 # discard header

# --- BUILD EMPTY DICTIONARY OF miRNAs and Pathways FROM RELATIONS (there may be some we're not interested in, i.e., not in pathways) ---
allRelMirs = set([])
allRelPaths = set([])
for r in relation:
    #print r
    allRelMirs.add(r[5])
    allRelPaths.add(r[0])
print 'Unique miRNAs in relation file = ',len(allRelMirs)
print 'Unique Pathways in relation file = ',len(allRelPaths)
mir_scores  = dict([ (x,[]) for x in list(allRelMirs) ])        # make dictionary from list of miRNAs (one list for each mir)
mir_targets = dict([ (x,[]) for x in list(allRelMirs) ])
path_genes = dict([ (x,set([])) for x in list(allRelPaths) ])   # make dictionary from list of Paths (one set for each path)
path_mirs  = dict([ (x,set([])) for x in list(allRelPaths) ])

# --- POPULATE DICTIONARIES ---
for r in relation:                                              # go through relations once and populate all dictionaries
    mir_scores[ r[5] ].append( [r[0],r[6]] )
    tmp_name = r[0] + " " + r[2]
    mir_targets[ r[5] ].append( [tmp_name, r[6]] )
    path_genes[ r[0] ].add( r[2] )
    path_mirs[ r[0] ].add( r[5] )
for key in path_genes.keys():                                   # convert sets to lists (just to be nice)
    path_genes[key] = list(path_genes[key])
    path_mirs[key] = list(path_mirs[key])

# --- mir_targets:= mirna ---> ("path gene", score), mir_scores := mirna ----> ("path", score)
# ---- path_mirs := path ----> mirna, path_genes:= path ---> gene
# pathways:= list of pathways (unique human), mirnas:= list of mirnas (unique human)
# -----Initialize output files ----- 