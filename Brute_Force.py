import random
import numpy
import csv
from Inhibitor_Modules import *
import os
import Graph as g
import math

def Brute():
	#------Create the list of all possible two-mirs -------
	mir_list = getMirsPath(curr_path)
	print len(mir_list)
	for i in range(len(mir_list)):
		for j in range(len(mir_list) - i):
			if i == i+j: continue
			else: population.append([mir_list[i+j], mir_list[i]])
	pop=population
	
	#----output the population at this point------
	#for i in pop:
	#	print i
	#print len(pop)
	
	#-----Ensure that the populations are binomial coefficients--------
	print len(pop)

	#-----Sort Each Population by Score to find Maximum Combinations ------- 
	sort_pop = Sort_Fitness(pop)


	#print sort_pop and the scores for the various sort pops.
	print sort_pop[0]
	print Score_Covering_No_Inh(sort_pop[0])
	
	return sort_pop
	
def Sort_Fitness (pop): # Don't worry about warning in Eclipse
    sort_pop = []
    pop_fit = ScoreFunction (pop, curr_path, mir_scores)
    # ---- Use Lambda function to sort list efficiently
    for mir_gene, score in sorted(zip(pop, pop_fit), key = lambda x:x[1], reverse = 1):
        sort_pop.append(mir_gene)
    
    return sort_pop
   
def ScoreFunction ( pop , path , mir_scores, inhibitors={} ):  # Takes in the population, pathway, and the mir_scores dictionary as inputs. Uses the Score function to build the fitness list. 
    inhibitors = {}
    curr_fitness = []
    for i in range(len(pop)):
        curr_fitness.append(Score_Covering_No_Inh(pop[i]))  #Need to Change to selectable type, for Score - need to add (path, mir_scores) as paramaters
	
    return curr_fitness

def Score_Covering_No_Inh (mir_gene):
	# --- Ensure that the no-length paths are given a score of zero
	if len(mir_gene) == 0: return 0.0 
	
	end_score = 0.0
	tmp_genes = []
	for mir in mir_gene:
		for gene_name, score in mir_targets[mir]:
			gene_name = gene_name.split()
			if gene_name[0] == curr_path:
				tmp_genes.append(gene_name[1])
	tmp_genes = list(set(tmp_genes))
	# ------- Normalize score to length of the current path ----
	score = float(len(tmp_genes))/float(len(path_genes[curr_path]))
	# --- Select against long combinations of miRNA's
	# ----- Normalizing to the eight root of the length of the combination appears to select for the 3-5 mirs----
	return score/float(math.sqrt(math.sqrt(math.sqrt((len(mir_gene))))))
	
def getMirsPath (path):
    tmp_mirs = []
    tmp_mirs.extend(path_mirs[path])
    mir_list = []
    for i in tmp_mirs:
        if i in mirnas:
            mir_list.append(i)
    return mir_list


# ------ Main Method Below, called as soon as the program begins -------- 

curr_path = "hsa04150"    # This is the test pathway

# --- Open requisite files ---    
os.chdir("/Users/nithintumma/Desktop")
mirlist_file = open("miRNAs.csv", "rU")
#genes_file = open("genes.csv", "rU")
paths_file = open("pathways.csv", "rU")

pathways = list(csv.reader(paths_file))    # read unique pathways into list
pathways = [ x[0] for x in pathways ]

# --- READ miRNA MASTER LIST (only those in pathways) ---
mirnas = list(csv.reader(mirlist_file, quotechar='"'))    # read unique miRNAs into list
mirnas = [ x[0] for x in mirnas ]

# --- READ RELATIONS FILE INTO LIST ---
# "Pathway","Gene ID","Gene Name","Copies in Path","Transcript","miRNA","Score","Start","End"
#  0         1         2           3                4            5       6       7       8
relation = list(csv.reader(open("relation.csv", "rU"), quotechar='"'))  # get relations
relation = relation[1:]                                 # discard header

# --- BUILD EMPTY DICTIONARY OF miRNAs and Pathways FROM RELATIONS (there may be some we're not interested in, i.e., not in pathways) ---
allRelMirs = set([])
allRelPaths = set([])
for r in relation:
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
    mir_scores[ r[5] ].append( [r[0], float(r[6])] )
    tmp_name = r[0] + " " + r[2]
    mir_targets[ r[5] ].append( [tmp_name, float(r[6])] )
    path_genes[ r[0] ].add( r[2] )
    path_mirs[ r[0] ].add( r[5] )
for key in path_genes.keys():                                   # convert sets to lists (just to be nice)
    path_genes[key] = list(path_genes[key])
    path_mirs[key] = list(path_mirs[key])

# --- mir_targets:= mirna ---> ("path gene", score), mir_scores := mirna ----> ("path", score)
# ---- path_mirs := path ----> mirna, path_genes:= path ---> gene
# pathways:= list of pathways (unique human), mirnas:= list of mirnas (unique human)
# -----Initialize output files ----- 
print("Initialization Complete")

#-----------Initialize population list
population = []
# --- Contains all of the human mirs that target the curr_path ----
mir_list = getMirsPath(curr_path)
lengths = [len(path) for path in path_genes]
print lengths
sort_len = []
for path, length in sorted(zip(lengths, pathways), key = lambda x:x[1], reverse = 0):
        sort_len.append(path)
print sort_len

#Brute()

