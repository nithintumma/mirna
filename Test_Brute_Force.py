import random
import numpy
import csv
from Inhibitor_Modules import *
import os
import Graph as g
import math

def Brute():
	#------Create the list of all possible two-mirs -------
	for i in range(len(mirnas)):
		for j in range(len(mirnas) - i):
			if i == i+j: continue
			else: population.append([mirnas[i+j], mirnas[i]])
	pop=population
	
	pop_3 = []
	#---Create the list of all possible three-mirs----
	for i in range(len(mirnas)):
		for j in mirnas[i::]:
			for k in mirnas[mirnas.index(j)::]:
				if len(set([mirnas[i], j, k])) < 3: continue
				else: pop_3.append([mirnas[i], j, k])
	
	pop_4 = []
	#---Create the list of all possible four-mirs----
	for i in mirnas:
		for j in mirnas[mirnas.index(i)::]:
			for k in mirnas[mirnas.index(j)::]:
				for l in mirnas[mirnas.index(k)::]:
					if len(set([i, j, k, l])) < 4: continue
					else: pop_4.append([i, j, k, l])
	pop_5 = []
	#---Create the list of all possible four-mirs----
	for i in mirnas:
		for j in mirnas[mirnas.index(i)::]:
			for k in mirnas[mirnas.index(j)::]:
				for l in mirnas[mirnas.index(k)::]:
					for m in mirnas[mirnas.index(l)::]:
						if len(set([i, j, k, l, m])) < 5: continue
						else: pop_5.append([i, j, k, l, m])
	
	#----output the population at this point------
	#for i in pop:
	#	print i
	#print len(pop)
	
	#-----Ensure that the populations are binomial coefficients--------
	print len(pop)
	print len(pop_3)
	print len(pop_4)
	print len(pop_5)
	
	#-----Sort Each Population by Score to find Maximum Combinations ------- 
	sort_pop = Sort_Fitness(pop)
	sort_pop_3 = Sort_Fitness(pop_3)
	sort_pop_4 = Sort_Fitness(pop_4)
	sort_pop_5 = Sort_Fitness(pop_5)

	#print sort_pop and the scores for the various sort pops.
	print sort_pop[0]
	print Score_Covering_No_Inh(sort_pop[0])
	
	print sort_pop_3[0]
	print Score_Covering_No_Inh(sort_pop_3[0])

	print sort_pop_4[0]
	print Score_Covering_No_Inh(sort_pop_4[0])
	
	print sort_pop_5[0]
	print Score_Covering_No_Inh(sort_pop_5[0])
	
	
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

curr_path = "hsa2000"    # This is the test pathway

# --- Open requisite files ---    
os.chdir("/home/nithin/Desktop/Data")
mirlist_file = open("miRNAs_test.csv", "rU")
paths_file = open("pathways_test.csv", "rU")

pathways = list(csv.reader(paths_file))    # read unique pathways into list
pathways = [ x[0] for x in pathways ]

# --- READ miRNA MASTER LIST (only those in pathways) ---
mirnas = list(csv.reader(mirlist_file, quotechar='"'))    # read unique miRNAs into list
mirnas = [ x[0] for x in mirnas ]

# --- READ RELATIONS FILE INTO LIST ---
# "Pathway","Gene ID","Gene Name","Copies in Path","Transcript","miRNA","Score","Start","End"
#  0         1         2           3                4            5       6       7       8
relation = list(csv.reader(open("relations_test.csv", "rU"), quotechar='"'))  # get relations
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
o_file = open("/home/nithin/Desktop/output.txt", "r+")   # The output file
fit_file = open("/home/nithin/Desktop/fitness.csv", "r+") # The fitness output file
maxfit_file = open("/home/nithin/Desktop/maxfit.csv", "r+")
print("Initialization Complete")

#-----------Initialize population list
population = []


# --- Contains all of the human mirs that target the curr_path ----
mir_list = getMirsPath(curr_path)
Brute()

